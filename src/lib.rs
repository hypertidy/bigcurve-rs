//! Adaptive densification of lines under a forward map.
//!
//! An edge between two longitude, latitude points is a great circle arc.
//! Move it into a map projection and its true image is a curve; joining
//! the moved ends draws the wrong line. bigcurve adds the vertices an edge
//! needs to keep its shape in the target plane: more where the map bends
//! it, none where it doesn't.
//!
//! The method is the adaptive resampling of d3-geo (`resampleLineTo`, Mike
//! Bostock): bisect an edge while the moved midpoint sits further than a
//! tolerance from the chord of the moved ends, or lands far from the
//! chord's middle. The tolerance is in target units.
//!
//! Three choices drive it:
//!
//! - an [`Edge`] model says what the middle of an edge is where it was
//!   drawn: [`GreatCircle`] (the default, longitude and latitude in
//!   degrees) or [`Planar`] (straight lines in the source coordinates);
//! - a forward map ([`Forward`], one point at a time, or [`ForwardBatch`],
//!   arrays of points per call) takes source coordinates to the target
//!   plane: a PROJ transformation, an analytic projection, or the map from
//!   positions to the index space of a grid. Points it cannot take never
//!   error; edges touching them are left alone, and a midpoint it cannot
//!   take stops that edge's refinement there;
//! - [`Options`]: the tolerance and a maximum depth.
//!
//! With [`GreatCircle`] this is a port of the engine of the R package
//! bigcurve (<https://github.com/hypertidy/bigcurve>, `src/bigcurve.h`):
//! it gives bitwise identical output for the same forward map, and the
//! tests replay fixtures recorded from the C++ engine with PROJ.
//!
//! Refinement runs level by level: every edge that needs a midpoint at
//! one level gets it from a single [`ForwardBatch`] call. Each edge's
//! decisions depend only on its own end points, so the vertices are
//! bitwise those of d3-geo's depth-first recursion
//! ([`Resampler::edge`]), in the same order.
//!
//! ```
//! use bigcurve::{densify_path, from_fn, Options};
//!
//! // a spherical Mercator, standing in for a projection library
//! let merc = from_fn(|lon: f64, lat: f64| {
//!     let r = 6378137.0;
//!     let y = (std::f64::consts::FRAC_PI_4 + lat.to_radians() / 2.0).tan().ln();
//!     Some([r * lon.to_radians(), r * y])
//! });
//! let line = [[0.0, 60.0], [100.0, 60.0]];
//! let d = densify_path(merc, &line, &Options::new(5e4));
//! assert!(d.len() > 2);
//! assert_eq!(d[0], line[0]);
//! assert_eq!(d[d.len() - 1], line[1]);
//! ```
//!
//! Planar edges into a curved space, here the index space of a polar
//! sector grid (radius by angle), with the tolerance a fraction of a
//! cell:
//!
//! ```
//! use bigcurve::{densify_path_with, from_fn, Options, Planar};
//!
//! let to_index = from_fn(|x: f64, y: f64| {
//!     Some([x.hypot(y) / 10.0, y.atan2(x).to_degrees()])
//! });
//! let d = densify_path_with(to_index, Planar::new(), &[[100.0, 0.0], [0.0, 100.0]], &Options::new(1e-3));
//! assert!(d.len() > 10);
//! ```
//!
//! An exact map is needed for refinement to trigger: a forward map that
//! itself interpolates linearly between samples (an approximate
//! transformer) puts every midpoint on its chord, and nothing is added.

#![forbid(unsafe_code)]

use std::collections::HashMap;
use std::fmt;

mod edge;
mod forward;

pub use edge::{Edge, GreatCircle, Planar, Span};
pub use forward::{from_fn, Forward, ForwardBatch, FromFn};

pub(crate) const DEG2RAD: f64 = 0.017453292519943295;
pub(crate) const RAD2DEG: f64 = 57.29577951308232;

/// The default maximum depth of bisection per edge, as in d3-geo.
pub const DEFAULT_MAX_DEPTH: u32 = 16;

/// Densification settings.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Options {
    /// How far, in target units, a chord may stray from the true curve.
    pub tolerance: f64,
    /// Maximum levels of bisection per input edge (at most
    /// `2^max_depth - 1` new vertices on one edge).
    pub max_depth: u32,
}

impl Options {
    /// Settings with the given tolerance and the default maximum depth.
    pub fn new(tolerance: f64) -> Self {
        Options {
            tolerance,
            max_depth: DEFAULT_MAX_DEPTH,
        }
    }

    /// The same settings with another maximum depth.
    pub fn max_depth(self, max_depth: u32) -> Self {
        Options { max_depth, ..self }
    }
}

/// A vertex in both spaces, so that no point is moved twice.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Vertex {
    /// Source x: longitude in degrees for [`GreatCircle`].
    pub lon: f64,
    /// Source y: latitude in degrees for [`GreatCircle`].
    pub lat: f64,
    /// Target plane x (meaningless when `ok` is false).
    pub x: f64,
    /// Target plane y (meaningless when `ok` is false).
    pub y: f64,
    /// Whether the forward map took this point.
    pub ok: bool,
}

impl Vertex {
    /// A vertex from source and target coordinates, for callers that
    /// moved the point themselves (a lattice transformed as scanlines).
    /// `ok` is whether both target coordinates are finite.
    pub fn new(source: [f64; 2], target: [f64; 2]) -> Self {
        Vertex {
            lon: source[0],
            lat: source[1],
            x: target[0],
            y: target[1],
            ok: target[0].is_finite() && target[1].is_finite(),
        }
    }

    fn src(&self) -> [f64; 2] {
        [self.lon, self.lat]
    }
}

/// The densification engine: a forward map, an edge model, a tolerance
/// and a depth.
///
/// Most uses want [`densify_path`], [`densify_mesh`] or their `_with`
/// forms; the resampler is for drivers of other shapes (rings held
/// elsewhere, streaming input, edges whose ends were moved by the
/// caller).
pub struct Resampler<T, E = GreatCircle> {
    forward: T,
    edge: E,
    tol2: f64,
    max_depth: u32,
    floor: f64,
    // scratch for the level-order engine
    work: Vec<Item>,
    next: Vec<Item>,
    nodes: Vec<Node>,
    bx: Vec<f64>,
    by: Vec<f64>,
    near: Vec<[f64; 2]>,
    ok: Vec<bool>,
}

impl<T, E> fmt::Debug for Resampler<T, E> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("Resampler")
            .field("tolerance", &self.tol2.sqrt())
            .field("max_depth", &self.max_depth)
            .finish_non_exhaustive()
    }
}

/// Where a refined midpoint hangs in its edge's bisection tree.
#[derive(Clone, Copy, Debug)]
enum Attach {
    Root(usize),
    Left(usize),
    Right(usize),
}

/// A pending bisection: the sub-edge `a -> b` with `depth` levels left.
#[derive(Clone, Copy, Debug)]
struct Item {
    a: Vertex,
    b: Vertex,
    depth: u32,
    at: Attach,
}

const NONE: usize = usize::MAX;

#[derive(Clone, Copy, Debug)]
struct Node {
    v: Vertex,
    left: usize,
    right: usize,
}

impl<T: ForwardBatch> Resampler<T, GreatCircle> {
    /// A resampler for great circle edges, the given forward map and
    /// settings.
    ///
    /// An arc spanning less ground than the tolerance can never need
    /// interior vertices, so refinement stops there whatever the planar
    /// test says (see [`GreatCircle`]).
    pub fn new(forward: T, options: &Options) -> Self {
        Self::with_edge(forward, GreatCircle, options)
    }
}

impl<T: ForwardBatch, E: Edge> Resampler<T, E> {
    /// A resampler for the given forward map, edge model and settings.
    pub fn with_edge(forward: T, edge: E, options: &Options) -> Self {
        let tol = options.tolerance;
        Resampler {
            floor: edge.floor(tol),
            forward,
            edge,
            tol2: tol * tol,
            max_depth: options.max_depth,
            work: Vec::new(),
            next: Vec::new(),
            nodes: Vec::new(),
            bx: Vec::new(),
            by: Vec::new(),
            near: Vec::new(),
            ok: Vec::new(),
        }
    }

    /// Give back the forward map.
    pub fn into_inner(self) -> T {
        self.forward
    }

    /// Move one point, keeping both coordinates.
    pub fn project(&mut self, lon: f64, lat: f64) -> Vertex {
        let (mut x, mut y, mut ok) = ([lon], [lat], [false]);
        self.forward
            .forward_batch(&mut x, &mut y, &[[f64::NAN, f64::NAN]], &mut ok);
        vertex(lon, lat, x[0], y[0], ok[0])
    }

    /// Move many points, keeping both coordinates, in one call.
    pub fn project_all(&mut self, coords: &[[f64; 2]]) -> Vec<Vertex> {
        let mut x: Vec<f64> = coords.iter().map(|c| c[0]).collect();
        let mut y: Vec<f64> = coords.iter().map(|c| c[1]).collect();
        let near = vec![[f64::NAN, f64::NAN]; coords.len()];
        let mut ok = vec![false; coords.len()];
        self.forward.forward_batch(&mut x, &mut y, &near, &mut ok);
        coords
            .iter()
            .enumerate()
            .map(|(i, c)| vertex(c[0], c[1], x[i], y[i], ok[i]))
            .collect()
    }

    /// Densify one edge `a -> b`, appending its interior vertices (in
    /// order, endpoints excluded) to `out`, depth first, one forward call
    /// per midpoint: d3-geo's recursion.
    ///
    /// Edges with an end the map cannot take are left alone, and a
    /// midpoint it cannot take stops that branch. Interior vertices are in
    /// the frame of the written start `a` (see [`Edge::align`]).
    ///
    /// The arithmetic always runs in one canonical direction, from the
    /// endpoint with the lexicographically smaller source coordinates,
    /// and the chain is reversed afterwards if the edge was written the
    /// other way. So `a -> b` and `b -> a` give bitwise identical vertices
    /// (up to an exact multiple of 360 in longitude when a great circle
    /// edge crosses the antimeridian), and an edge shared by two rings
    /// refines to the same vertices from both sides.
    pub fn edge(&mut self, a: &Vertex, b: &Vertex, out: &mut Vec<Vertex>) {
        let Some((p, q, rev)) = self.canonical(a, b) else {
            return;
        };
        let start = out.len();
        self.split(&p, &q, self.max_depth, out);
        if rev {
            self.unreverse(a, &q, &mut out[start..]);
        }
    }

    /// Densify many edges level by level: each level's midpoints are
    /// moved in one [`ForwardBatch`] call. Appends each edge's interior
    /// vertices to `out` and its end offset to `ends` (so edge `k`'s
    /// vertices are `out[ends[k - 1]..ends[k]]`, from where `out` stood
    /// on entry for `k = 0`).
    ///
    /// The vertices are bitwise those of [`Resampler::edge`] edge by
    /// edge; only the order of the forward calls differs.
    pub fn edges(
        &mut self,
        edges: &[(Vertex, Vertex)],
        out: &mut Vec<Vertex>,
        ends: &mut Vec<usize>,
    ) {
        let mut canon = Vec::with_capacity(edges.len());
        self.work.clear();
        self.nodes.clear();
        let mut roots = vec![NONE; edges.len()];
        for (k, (a, b)) in edges.iter().enumerate() {
            let c = self.canonical(a, b);
            if let Some((p, q, _)) = c {
                self.work.push(Item {
                    a: p,
                    b: q,
                    depth: self.max_depth,
                    at: Attach::Root(k),
                });
            }
            canon.push(c);
        }

        while !self.work.is_empty() {
            // pre-midpoint tests; survivors queue a midpoint
            self.bx.clear();
            self.by.clear();
            self.near.clear();
            let mut keep = 0;
            for i in 0..self.work.len() {
                let it = self.work[i];
                if let Some(m) = self.wants_midpoint(&it.a, &it.b, it.depth) {
                    self.bx.push(m[0]);
                    self.by.push(m[1]);
                    self.near.push([it.a.x, it.a.y]);
                    self.work[keep] = it;
                    keep += 1;
                }
            }
            self.work.truncate(keep);
            if self.work.is_empty() {
                break;
            }
            let n = self.work.len();
            let src: Vec<[f64; 2]> = (0..n).map(|i| [self.bx[i], self.by[i]]).collect();
            self.ok.clear();
            self.ok.resize(n, false);
            self.forward
                .forward_batch(&mut self.bx, &mut self.by, &self.near, &mut self.ok);

            self.next.clear();
            for (i, s) in src.iter().enumerate() {
                let it = self.work[i];
                let m = vertex(s[0], s[1], self.bx[i], self.by[i], self.ok[i]);
                if !m.ok || !self.go(&it.a, &it.b, &m) {
                    continue;
                }
                let id = self.nodes.len();
                self.nodes.push(Node {
                    v: m,
                    left: NONE,
                    right: NONE,
                });
                match it.at {
                    Attach::Root(k) => roots[k] = id,
                    Attach::Left(p) => self.nodes[p].left = id,
                    Attach::Right(p) => self.nodes[p].right = id,
                }
                self.next.push(Item {
                    a: it.a,
                    b: m,
                    depth: it.depth - 1,
                    at: Attach::Left(id),
                });
                self.next.push(Item {
                    a: m,
                    b: it.b,
                    depth: it.depth - 1,
                    at: Attach::Right(id),
                });
            }
            std::mem::swap(&mut self.work, &mut self.next);
        }

        // in-order walk of each edge's tree
        let mut stack = Vec::new();
        for (k, (a, _)) in edges.iter().enumerate() {
            let start = out.len();
            let mut cur = roots[k];
            loop {
                while cur != NONE {
                    stack.push(cur);
                    cur = self.nodes[cur].left;
                }
                let Some(top) = stack.pop() else { break };
                out.push(self.nodes[top].v);
                cur = self.nodes[top].right;
            }
            if let Some((_, q, true)) = canon[k] {
                self.unreverse(a, &q, &mut out[start..]);
            }
            ends.push(out.len());
        }
    }

    /// The edge in canonical direction, its far end aligned, and whether
    /// it was reversed; `None` when an end cannot be moved.
    fn canonical(&self, a: &Vertex, b: &Vertex) -> Option<(Vertex, Vertex, bool)> {
        if !a.ok || !b.ok {
            return None;
        }
        let rev = b.lon < a.lon || (b.lon == a.lon && b.lat < a.lat);
        let (p, mut q) = if rev { (*b, *a) } else { (*a, *b) };
        let [lon, lat] = self.edge.align(p.src(), q.src());
        q.lon = lon;
        q.lat = lat;
        Some((p, q, rev))
    }

    /// Reverse a chain refined from `b` towards `q` (the written start
    /// `a` aligned to `b`) and move it into `a`'s frame.
    fn unreverse(&self, a: &Vertex, q: &Vertex, chain: &mut [Vertex]) {
        chain.reverse();
        let [sx, sy] = self.edge.reframe(a.src(), q.src());
        if sx != 0.0 || sy != 0.0 {
            for v in chain {
                v.lon += sx;
                v.lat += sy;
            }
        }
    }

    /// The tests before a midpoint is moved: depth, the edge model's
    /// span, and a chord already within twice the tolerance. Returns the
    /// midpoint to move, if any.
    fn wants_midpoint(&self, a: &Vertex, b: &Vertex, depth: u32) -> Option<[f64; 2]> {
        if depth == 0 {
            return None;
        }
        let span = self.edge.span(a.src(), b.src(), self.floor);
        if span == Span::Negligible {
            return None;
        }
        let dx = b.x - a.x;
        let dy = b.y - a.y;
        let d2 = dx * dx + dy * dy;
        if span != Span::Long && d2 <= 4.0 * self.tol2 {
            return None;
        }
        Some(self.edge.mid(a.src(), b.src()))
    }

    /// d3-geo's test on the moved midpoint: off the chord by more than
    /// the tolerance, or outside its middle three fifths. Long edges
    /// always split.
    fn go(&self, a: &Vertex, b: &Vertex, m: &Vertex) -> bool {
        if self.edge.span(a.src(), b.src(), self.floor) == Span::Long {
            return true;
        }
        let dx = b.x - a.x;
        let dy = b.y - a.y;
        let d2 = dx * dx + dy * dy;
        if d2 <= 0.0 {
            return false;
        }
        let dx2 = m.x - a.x;
        let dy2 = m.y - a.y;
        let dz = dx * dy2 - dy * dx2; // twice the signed area
        let perp2 = (dz * dz) / d2; // squared perpendicular distance
        let tpar = (dx * dx2 + dy * dy2) / d2; // position along the chord
        perp2 > self.tol2 || !(0.2..=0.8).contains(&tpar)
    }

    fn split(&mut self, a: &Vertex, b: &Vertex, depth: u32, out: &mut Vec<Vertex>) {
        let Some([lonm, latm]) = self.wants_midpoint(a, b, depth) else {
            return;
        };
        let m = self.project(lonm, latm);
        if !m.ok {
            return;
        }
        if self.go(a, b, &m) {
            self.split(a, &m, depth - 1, out);
            out.push(m);
            self.split(&m, b, depth - 1, out);
        }
    }
}

fn vertex(lon: f64, lat: f64, x: f64, y: f64, ok: bool) -> Vertex {
    let ok = ok && x.is_finite() && y.is_finite();
    Vertex {
        lon,
        lat,
        x: if ok || x.is_finite() { x } else { f64::NAN },
        y: if ok || y.is_finite() { y } else { f64::NAN },
        ok,
    }
}

/// Densify one path of longitude, latitude vertices with great circle
/// edges. See [`densify_path_with`].
pub fn densify_path<T: ForwardBatch>(
    forward: T,
    coords: &[[f64; 2]],
    options: &Options,
) -> Vec<[f64; 2]> {
    densify_path_with(forward, GreatCircle, coords, options)
}

/// Densify one path under an edge model.
///
/// Returns the path with new vertices in place, still in source
/// coordinates; the input vertices are kept exactly. Each edge is refined
/// on its own (a path that revisits an edge refines it again, to the same
/// vertices).
pub fn densify_path_with<T: ForwardBatch, E: Edge>(
    forward: T,
    edge: E,
    coords: &[[f64; 2]],
    options: &Options,
) -> Vec<[f64; 2]> {
    let mut rs = Resampler::with_edge(forward, edge, options);
    let v = rs.project_all(coords);
    let pairs: Vec<(Vertex, Vertex)> = v.windows(2).map(|w| (w[0], w[1])).collect();
    let (mut interior, mut ends) = (Vec::new(), Vec::with_capacity(pairs.len()));
    rs.edges(&pairs, &mut interior, &mut ends);
    let mut out = Vec::with_capacity(coords.len() + interior.len());
    let mut from = 0;
    for (k, &end) in ends.iter().enumerate() {
        out.push(coords[k]);
        out.extend(interior[from..end].iter().map(|m| [m.lon, m.lat]));
        from = end;
    }
    if let Some(last) = coords.last() {
        out.push(*last);
    }
    out
}

/// A densified segment mesh.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct Mesh {
    /// The input vertices, unchanged, followed by the new ones.
    pub vertices: Vec<[f64; 2]>,
    /// Each input segment replaced by its chain of refined segments, in
    /// input order.
    pub segments: Vec<[usize; 2]>,
    /// For each output segment, the input segment it descends from; the
    /// hook for carrying per-segment attributes (feature, ring, arc)
    /// through refinement.
    pub parent: Vec<usize>,
}

/// An error from [`densify_mesh`].
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Error {
    /// A segment refers to a vertex that does not exist.
    IndexOutOfRange {
        segment: usize,
        index: usize,
        vertices: usize,
    },
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Error::IndexOutOfRange {
                segment,
                index,
                vertices,
            } => write!(
                f,
                "segment {segment} refers to vertex {index}, but there are {vertices} vertices"
            ),
        }
    }
}

impl std::error::Error for Error {}

/// Densify a segment mesh of longitude, latitude vertices with great
/// circle edges. See [`densify_mesh_with`].
pub fn densify_mesh<T: ForwardBatch>(
    forward: T,
    vertices: &[[f64; 2]],
    segments: &[[usize; 2]],
    options: &Options,
) -> Result<Mesh, Error> {
    densify_mesh_with(forward, GreatCircle, vertices, segments, options)
}

/// Densify a segment mesh under an edge model: vertices and index pairs
/// into them.
///
/// A segment that occurs more than once with the same two vertex indices
/// (in either direction, as a shared boundary does when neighbouring
/// rings share vertices) is refined once; later occurrences reuse the
/// vertices made the first time, reversed as needed, so neighbours stay
/// joined.
pub fn densify_mesh_with<T: ForwardBatch, E: Edge>(
    forward: T,
    edge: E,
    vertices: &[[f64; 2]],
    segments: &[[usize; 2]],
    options: &Options,
) -> Result<Mesh, Error> {
    let nv = vertices.len();
    for (k, s) in segments.iter().enumerate() {
        for &index in s {
            if index >= nv {
                return Err(Error::IndexOutOfRange {
                    segment: k,
                    index,
                    vertices: nv,
                });
            }
        }
    }

    let mut rs = Resampler::with_edge(forward, edge, options);
    let v = rs.project_all(vertices);

    // one refinement per undirected pair, in the direction first seen
    let mut first: HashMap<(usize, usize), usize> = HashMap::new();
    let mut unique: Vec<[usize; 2]> = Vec::new();
    for &[a, b] in segments {
        first.entry((a.min(b), a.max(b))).or_insert_with(|| {
            unique.push([a, b]);
            unique.len() - 1
        });
    }
    let pairs: Vec<(Vertex, Vertex)> = unique.iter().map(|&[a, b]| (v[a], v[b])).collect();
    let (mut interior, mut ends) = (Vec::new(), Vec::with_capacity(pairs.len()));
    rs.edges(&pairs, &mut interior, &mut ends);

    let mut out = Mesh {
        vertices: vertices.to_vec(),
        segments: Vec::with_capacity(segments.len() * 2),
        parent: Vec::with_capacity(segments.len() * 2),
    };
    // new vertices are numbered in order of first occurrence
    let mut placed: Vec<usize> = vec![NONE; unique.len()];
    let mut chain = Vec::new();
    for (k, &[a, b]) in segments.iter().enumerate() {
        let u = first[&(a.min(b), a.max(b))];
        let from = if u == 0 { 0 } else { ends[u - 1] };
        let n = ends[u] - from;
        if placed[u] == NONE {
            placed[u] = out.vertices.len();
            out.vertices
                .extend(interior[from..ends[u]].iter().map(|m| [m.lon, m.lat]));
        }
        chain.clear();
        chain.extend(placed[u]..placed[u] + n);
        if unique[u][0] != a {
            chain.reverse();
        }
        let mut prev = a;
        for &idx in &chain {
            out.segments.push([prev, idx]);
            out.parent.push(k);
            prev = idx;
        }
        out.segments.push([prev, b]);
        out.parent.push(k);
    }
    Ok(out)
}

/// A tolerance of one part in `pixels` of the larger side of the moved
/// extent of `coords`: about a pixel when the result is drawn `pixels`
/// wide (the R package uses 2048).
///
/// The larger side, so that a single parallel or meridian (zero extent in
/// one direction) still gets a tolerance. `None` when no point can be
/// moved or the extent is degenerate.
pub fn default_tolerance<T: ForwardBatch>(
    mut forward: T,
    coords: &[[f64; 2]],
    pixels: f64,
) -> Option<f64> {
    let mut x: Vec<f64> = coords.iter().map(|c| c[0]).collect();
    let mut y: Vec<f64> = coords.iter().map(|c| c[1]).collect();
    let near = vec![[f64::NAN, f64::NAN]; coords.len()];
    let mut ok = vec![false; coords.len()];
    forward.forward_batch(&mut x, &mut y, &near, &mut ok);
    let (mut xmin, mut xmax, mut ymin, mut ymax) = (
        f64::INFINITY,
        f64::NEG_INFINITY,
        f64::INFINITY,
        f64::NEG_INFINITY,
    );
    for i in 0..coords.len() {
        if ok[i] && x[i].is_finite() && y[i].is_finite() {
            xmin = xmin.min(x[i]);
            xmax = xmax.max(x[i]);
            ymin = ymin.min(y[i]);
            ymax = ymax.max(y[i]);
        }
    }
    let side = (xmax - xmin).max(ymax - ymin);
    if side.is_finite() && side > 0.0 {
        Some(side / pixels)
    } else {
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn mesh_rejects_bad_indices() {
        let id = from_fn(|lon, lat| Some([lon, lat]));
        let e = densify_mesh(id, &[[0.0, 0.0]], &[[0, 1]], &Options::new(1.0)).unwrap_err();
        assert_eq!(
            e,
            Error::IndexOutOfRange {
                segment: 0,
                index: 1,
                vertices: 1
            }
        );
    }
}
