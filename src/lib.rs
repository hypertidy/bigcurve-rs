//! Adaptive densification of lines under a forward map.
//!
//! An edge between two longitude, latitude points is a great circle arc.
//! Move it into a map projection and its true image is a curve; joining
//! the moved ends draws the wrong line. bigcurve adds the vertices an edge
//! needs to keep its shape in the target plane: more where the map bends
//! it, none where it doesn't.
//!
//! The method is the adaptive resampling of d3-geo (`resampleLineTo`, Mike
//! Bostock): bisect an edge while the moved great circle midpoint sits
//! further than a tolerance from the chord of the moved ends, or lands far
//! from the chord's middle. The tolerance is in target units.
//!
//! This is a port of the engine of the R package bigcurve
//! (<https://github.com/hypertidy/bigcurve>, `src/bigcurve.h`). It gives
//! bitwise identical output for the same forward map; the tests replay
//! fixtures recorded from the C++ engine with PROJ.
//!
//! The forward map is anything that implements [`Forward`]: a PROJ
//! transformation, an analytic projection, or the map from positions to
//! the index space of a grid. Points it cannot take never error; edges
//! touching them are left alone.
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

#![forbid(unsafe_code)]

use std::collections::HashMap;
use std::fmt;

const DEG2RAD: f64 = 0.017453292519943295;
const RAD2DEG: f64 = 57.29577951308232;

/// Mean Earth radius in metres, used to turn the tolerance into an arc
/// floor (see [`Resampler::new`]).
const R_MEAN: f64 = 6371008.8;

/// cos(30 degrees): arcs wider than this always split (d3-geo's
/// `cosMinDistance`).
const COS_MIN_DISTANCE: f64 = 0.8660254037844387;

/// The default maximum depth of bisection per edge, as in d3-geo.
pub const DEFAULT_MAX_DEPTH: u32 = 16;

/// A forward map from longitude, latitude (degrees) to the target plane.
///
/// Return `None` for a point that cannot be moved (beyond a horizon, out
/// of the map's domain). Non-finite output is treated the same way.
pub trait Forward {
    fn forward(&mut self, lon: f64, lat: f64) -> Option<[f64; 2]>;
}

impl<T: Forward + ?Sized> Forward for &mut T {
    fn forward(&mut self, lon: f64, lat: f64) -> Option<[f64; 2]> {
        (**self).forward(lon, lat)
    }
}

/// A [`Forward`] map made from a closure; see [`from_fn`].
#[derive(Clone, Copy, Debug)]
pub struct FromFn<F>(pub F);

impl<F: FnMut(f64, f64) -> Option<[f64; 2]>> Forward for FromFn<F> {
    fn forward(&mut self, lon: f64, lat: f64) -> Option<[f64; 2]> {
        (self.0)(lon, lat)
    }
}

/// Wrap a closure `(lon, lat) -> Option<[x, y]>` as a [`Forward`] map.
pub fn from_fn<F: FnMut(f64, f64) -> Option<[f64; 2]>>(f: F) -> FromFn<F> {
    FromFn(f)
}

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
    /// Longitude, degrees.
    pub lon: f64,
    /// Latitude, degrees.
    pub lat: f64,
    /// Target plane x (meaningless when `ok` is false).
    pub x: f64,
    /// Target plane y (meaningless when `ok` is false).
    pub y: f64,
    /// Whether the forward map took this point.
    pub ok: bool,
}

/// Normalise a longitude difference to [-180, 180], with the tie at
/// exactly +-180 keeping the sign of the input. An edge spanning exactly
/// 180 degrees is the one case where the two representatives are
/// equidistant, and the written direction is the only correct tiebreaker:
/// (0 -> 180) stays eastward, (0 -> -180) westward.
fn unwrap_delta(delta: f64) -> f64 {
    let mut d = delta - 360.0 * (delta / 360.0).round();
    if d.abs() == 180.0 {
        d = 180.0_f64.copysign(delta);
    }
    d
}

/// Spherical midpoint of the great circle arc between two points, degrees
/// in and out. The longitude representative is the one nearest the
/// parents' mean longitude (callers unwrap `lon2` relative to `lon1`
/// first), so each edge keeps its own convention: an edge along 180
/// refines at 180, along -180 at -180, and an edge written as crossing
/// gets a locally continuous chain, possibly outside [-180, 180). The
/// engine has no seam.
fn gc_midpoint(lon1: f64, lat1: f64, lon2: f64, lat2: f64) -> (f64, f64) {
    let frame_mid = (lon1 + lon2) * 0.5;
    let rlat1 = lat1 * DEG2RAD;
    let rlat2 = lat2 * DEG2RAD;
    let dlon = (lon2 - lon1) * DEG2RAD;
    let bx = rlat2.cos() * dlon.cos();
    let by = rlat2.cos() * dlon.sin();
    let cl1bx = rlat1.cos() + bx;
    let latm = (rlat1.sin() + rlat2.sin()).atan2((cl1bx * cl1bx + by * by).sqrt()) * RAD2DEG;
    // degenerate: same meridian, or an arc through or near a pole where the
    // azimuthal term is indeterminate; longitude is the frame mean (d3's
    // own fallback)
    if by.abs() < 1e-14 && cl1bx.abs() < 1e-14 {
        return (frame_mid, latm);
    }
    let lon = lon1 + by.atan2(cl1bx) * RAD2DEG;
    (frame_mid + unwrap_delta(lon - frame_mid), latm)
}

/// Cosine of the angular distance between two points (degrees).
fn cos_angular(lon1: f64, lat1: f64, lon2: f64, lat2: f64) -> f64 {
    let rlat1 = lat1 * DEG2RAD;
    let rlat2 = lat2 * DEG2RAD;
    rlat1.sin() * rlat2.sin() + rlat1.cos() * rlat2.cos() * ((lon2 - lon1) * DEG2RAD).cos()
}

/// The densification engine: a forward map, a tolerance and a depth.
///
/// Most uses want [`densify_path`] or [`densify_mesh`]; the resampler is
/// for drivers of other shapes (rings held elsewhere, streaming input).
pub struct Resampler<T> {
    forward: T,
    tol2: f64,
    max_depth: u32,
    cos_floor: f64,
}

impl<T: Forward> Resampler<T> {
    /// A resampler for the given forward map and settings.
    ///
    /// An arc spanning less ground than the tolerance can never need
    /// interior vertices, so refinement stops there whatever the planar
    /// test says. This bounds recursion on edges that cross a map's
    /// discontinuity, whose moved midpoints jump to the far side. The
    /// conversion to an angle assumes the tolerance is in metres; for other
    /// target units the floor is merely conservative.
    pub fn new(forward: T, options: &Options) -> Self {
        let tol = options.tolerance;
        let ang = if tol > 0.0 { tol / R_MEAN } else { 0.0 };
        Resampler {
            forward,
            tol2: tol * tol,
            max_depth: options.max_depth,
            cos_floor: ang.cos(),
        }
    }

    /// Give back the forward map.
    pub fn into_inner(self) -> T {
        self.forward
    }

    /// Move one point, keeping both coordinates.
    pub fn project(&mut self, lon: f64, lat: f64) -> Vertex {
        match self.forward.forward(lon, lat) {
            Some([x, y]) if x.is_finite() && y.is_finite() => Vertex {
                lon,
                lat,
                x,
                y,
                ok: true,
            },
            Some([x, y]) => Vertex {
                lon,
                lat,
                x,
                y,
                ok: false,
            },
            None => Vertex {
                lon,
                lat,
                x: f64::NAN,
                y: f64::NAN,
                ok: false,
            },
        }
    }

    /// Densify one edge `a -> b`, appending its interior vertices (in
    /// order, endpoints excluded) to `out`.
    ///
    /// Edges with an end the map cannot take are left alone, and a
    /// midpoint it cannot take stops that branch. Interior vertices are in
    /// the longitude frame of the written start `a`.
    ///
    /// The arithmetic always runs in one canonical direction, from the
    /// endpoint with the lexicographically smaller (lon, lat), and the
    /// chain is reversed afterwards if the edge was written the other way.
    /// So `a -> b` and `b -> a` give bitwise identical vertices (up to an
    /// exact multiple of 360 in longitude when the edge crosses the
    /// antimeridian), and an edge shared by two rings refines to the same
    /// vertices from both sides.
    pub fn edge(&mut self, a: &Vertex, b: &Vertex, out: &mut Vec<Vertex>) {
        if !a.ok || !b.ok {
            return;
        }
        let rev = b.lon < a.lon || (b.lon == a.lon && b.lat < a.lat);
        let (p, mut q) = if rev { (*b, *a) } else { (*a, *b) };
        q.lon = p.lon + unwrap_delta(q.lon - p.lon);
        let start = out.len();
        self.split(&p, &q, self.max_depth, out);
        if !rev {
            return;
        }
        out[start..].reverse();
        // p is b here and q is a unwrapped into b's frame; move the chain
        // into a's frame by the exact multiple of 360 between them
        let shift = 360.0 * ((a.lon - q.lon) / 360.0).round();
        if shift != 0.0 {
            for v in &mut out[start..] {
                v.lon += shift;
            }
        }
    }

    fn split(&mut self, a: &Vertex, b: &Vertex, depth: u32, out: &mut Vec<Vertex>) {
        if depth == 0 {
            return;
        }
        // arcs wider than 30 degrees always split: the midpoint test below
        // is blind to symmetric curves whose midpoint lands on the chord, so
        // it only adjudicates once arcs are short
        let cosang = cos_angular(a.lon, a.lat, b.lon, b.lat);
        // arc floor: shorter than the tolerance on the ground (also ends
        // coincident endpoints, cosang == 1)
        if cosang >= self.cos_floor {
            return;
        }
        let wide = cosang < COS_MIN_DISTANCE;

        let dx = b.x - a.x;
        let dy = b.y - a.y;
        let d2 = dx * dx + dy * dy;
        // chord shorter than twice the tolerance and arc not wide
        if !wide && d2 <= 4.0 * self.tol2 {
            return;
        }

        let (lonm, latm) = gc_midpoint(a.lon, a.lat, b.lon, b.lat);
        let m = self.project(lonm, latm);
        if !m.ok {
            return;
        }

        let mut go = wide;
        if !go && d2 > 0.0 {
            let dx2 = m.x - a.x;
            let dy2 = m.y - a.y;
            let dz = dx * dy2 - dy * dx2; // twice the signed area
            let perp2 = (dz * dz) / d2; // squared perpendicular distance
            let tpar = (dx * dx2 + dy * dy2) / d2; // position along the chord
            go = perp2 > self.tol2 || !(0.2..=0.8).contains(&tpar);
        }

        if go {
            self.split(a, &m, depth - 1, out);
            out.push(m);
            self.split(&m, b, depth - 1, out);
        }
    }
}

/// Densify one path of longitude, latitude vertices.
///
/// Returns the path with new vertices in place, still in longitude and
/// latitude; the input vertices are kept exactly. Each edge is refined on
/// its own (a path that revisits an edge refines it again, to the same
/// vertices).
pub fn densify_path<T: Forward>(
    forward: T,
    coords: &[[f64; 2]],
    options: &Options,
) -> Vec<[f64; 2]> {
    let mut rs = Resampler::new(forward, options);
    let v: Vec<Vertex> = coords.iter().map(|c| rs.project(c[0], c[1])).collect();
    let mut out = Vec::with_capacity(coords.len() * 2);
    let mut interior = Vec::new();
    for w in v.windows(2) {
        out.push([w[0].lon, w[0].lat]);
        interior.clear();
        rs.edge(&w[0], &w[1], &mut interior);
        out.extend(interior.iter().map(|m| [m.lon, m.lat]));
    }
    if let Some(last) = v.last() {
        out.push([last.lon, last.lat]);
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

/// Densify a segment mesh: vertices and index pairs into them.
///
/// A segment that occurs more than once with the same two vertex indices
/// (in either direction, as a shared boundary does when neighbouring
/// rings share vertices) is refined once; later occurrences reuse the
/// vertices made the first time, reversed as needed, so neighbours stay
/// joined.
pub fn densify_mesh<T: Forward>(
    forward: T,
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

    let mut rs = Resampler::new(forward, options);
    let v: Vec<Vertex> = vertices.iter().map(|c| rs.project(c[0], c[1])).collect();

    let mut out = Mesh {
        vertices: vertices.to_vec(),
        segments: Vec::with_capacity(segments.len() * 2),
        parent: Vec::with_capacity(segments.len() * 2),
    };
    // undirected pair -> (start vertex as first seen, first new index, count)
    let mut seen: HashMap<(usize, usize), (usize, usize, usize)> = HashMap::new();
    let mut interior = Vec::new();
    let mut chain = Vec::new();

    for (k, &[a, b]) in segments.iter().enumerate() {
        let key = (a.min(b), a.max(b));
        chain.clear();
        if let Some(&(from, first, n)) = seen.get(&key) {
            chain.extend(first..first + n);
            if from != a {
                chain.reverse();
            }
        } else {
            interior.clear();
            rs.edge(&v[a], &v[b], &mut interior);
            let first = out.vertices.len();
            out.vertices.extend(interior.iter().map(|m| [m.lon, m.lat]));
            chain.extend(first..first + interior.len());
            seen.insert(key, (a, first, interior.len()));
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
pub fn default_tolerance<T: Forward>(
    mut forward: T,
    coords: &[[f64; 2]],
    pixels: f64,
) -> Option<f64> {
    let (mut xmin, mut xmax, mut ymin, mut ymax) = (
        f64::INFINITY,
        f64::NEG_INFINITY,
        f64::INFINITY,
        f64::NEG_INFINITY,
    );
    for c in coords {
        if let Some([x, y]) = forward.forward(c[0], c[1]) {
            if x.is_finite() && y.is_finite() {
                xmin = xmin.min(x);
                xmax = xmax.max(x);
                ymin = ymin.min(y);
                ymax = ymax.max(y);
            }
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
    fn unwrap_delta_keeps_the_written_direction_at_180() {
        assert_eq!(unwrap_delta(180.0), 180.0);
        assert_eq!(unwrap_delta(-180.0), -180.0);
        assert_eq!(unwrap_delta(540.0), 180.0);
        assert_eq!(unwrap_delta(190.0), -170.0);
        assert_eq!(unwrap_delta(-190.0), 170.0);
    }

    #[test]
    fn midpoint_of_a_meridian_arc_stays_on_it() {
        let (lon, lat) = gc_midpoint(180.0, -60.0, 180.0, 60.0);
        assert_eq!(lon, 180.0);
        assert!(lat.abs() < 1e-12);
    }

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
