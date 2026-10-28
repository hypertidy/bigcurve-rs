# 01: Plan

Status: proposed. Each stage ends with something that runs and is
tested; later stages are sketches and will be revised as earlier ones
land.

## Summary

bigcurve is an adaptive densifier: it inserts vertices along edges only
where a forward map bends them, judged in the target plane. The R
package does this for longitude, latitude lines under PROJ map
projections, with edges read as great circles. Its engine is a short
header with no R or PROJ in it, written to be ported.

This crate ports that engine and generalises it in one direction: the
edge model becomes a parameter, and the forward map becomes any
function, not only a map projection. That turns bigcurve from a drawing
tool into a general piece of geometry plumbing: the step that keeps a
shape honest whenever it moves between spaces.

What it is not: a geometry library (no predicates, areas or overlay;
that is s2's and geo's work), a projection library (PROJ is plugged in,
not wrapped), or a clipper (points the map cannot take are left alone,
never errors).

## Why a crate of its own

Nothing in Rust combines an arbitrary forward map, a pluggable edge
model, a tolerance in the target plane, and identical refinement of
shared edges.

| crate | what it does | the gap |
|---|---|---|
| `geo` `Densify` (Euclidean, Haversine, Geodesic, Rhumb metric spaces) | uniform densify to a maximum segment length in the source metric | not adaptive, no target map; its `InterpolatePoint` metric spaces are the right shape for edge models and can be adapted |
| `d3_geo_rs` | Rust port of d3-geo, resampler included | the resampler lives inside d3's projection stream and its own projections; no arbitrary map, no shared-edge identity, no mesh input |
| `s2` (Rust port) | cells, points, regions | no edge tessellator (C++ s2's `S2EdgeTessellator` is the closest prior art, plate carree and Mercator only) |
| `geographiclib-rs` | ellipsoidal geodesics | an edge model, not an engine |
| PROJ `proj_trans_bounds` | densifies the sides of a box | uniform, boxes only |

bigcurve and s2 divide the work as in the R vignette: s2 answers what a
shape is on the sphere; bigcurve answers how many vertices it needs in
one target plane.

## The engine

Four parts, all traits or plain values:

1. **Edge model**: `fn mid(a, b) -> P` in source space, plus a way to
   say an edge is "long" (forcing a split, because the midpoint test can
   be fooled by symmetric curves) and "negligible" (an arc floor, so
   edges across a map's discontinuity stop recursing). Built in:
   - `GreatCircle`: the R package's model, the default for longitude,
     latitude; long means wider than 30 degrees, negligible means
     shorter than the tolerance on the ground.
   - `Planar`: straight lines in the source CRS, as GEOS, GDAL and most
     file formats read them; long and negligible restated in source
     units.
   - `Geodesic` (feature `geographiclib`): the sphere chooses where to
     split, each vertex is then moved to the same fraction of the WGS84
     geodesic, as `densify(ellipsoid = TRUE)` does in R.
2. **Forward map**: `fn fwd(&mut self, p) -> Option<(f64, f64)>`.
   `None` and non-finite output both mean "cannot be moved". Each point
   is mapped once; vertices carry source and target coordinates.
3. **Tolerance** in target units, and a maximum depth (16, as d3-geo).
4. **The loop**: d3-geo's test (perpendicular distance of the moved
   midpoint from the chord, and the midpoint landing between 0.2 and 0.8
   along it), run in one canonical direction per edge so `a -> b` and
   `b -> a` give bitwise identical vertices, with the longitude frame
   rules of the R engine (no seam, the input's representative kept).

Output: interior vertices per edge, and for mesh input, new segments
with the index of their parent segment.

## Uses this is shaped for

- **Drawing**: the R package's job. Lines and polygons in any PROJ
  projection, with a tolerance of about a pixel.
- **Index space of a non-regular grid.** To burn a polygon onto a grid
  that is not a regular lattice (curvilinear model grids, swaths, and any
  grid in a CRS other than the polygon's), the polygon is moved into the
  grid's index space (fractional column, row) and its edges must follow
  the curves they become there. With index space as the target plane:
  - the tolerance has a natural unit, a fraction of a cell (1e-3, say),
    the same for every grid and CRS;
  - each new vertex costs an inverse cell search on a curvilinear grid,
    so adding vertices only where the grid bends matters, and midpoints
    land next to their parents, so a walking search can start from the
    parent's cell;
  - shared zone boundaries refine identically, so zonal weights over a
    partition still sum to one.
  It is not the tool for rectilinear grids (cell walls are known, so
  cutting at them is exact) and it cannot remove the error that remains
  on curvilinear grids, which is about each cell's area, not the
  polygon's shape. Background:
  [grid front ends](https://github.com/hypertidy/controlledburn-rs)
  (design note on coverage for grids that are not a regular lattice).
- **Cell outlines in another CRS**: the exact world-space clip of a
  curvilinear cell traces its outline in the polygon's CRS; those
  outlines are curves there too.
- **Extents and windows**: the boundary of a grid moved into another
  CRS, to find the region a reprojection needs, without a fixed number
  of points per side.

## Stages

### 1. Faithful port (done)

Landed as `Forward` (plus `from_fn` for closures), `Resampler`,
`densify_path`, `densify_mesh` and `default_tolerance`; there is no edge
model trait yet, the great circle is built in. The oracle is the C++
header itself, run with PROJ through copies of the R package's path and
mesh drivers (`oracle/`); every forward call is recorded, so the tests
replay PROJ without linking it and fail if the Rust engine asks for a
point the C++ engine did not. All fixtures match bit for bit on Linux
with glibc; elsewhere the platform's trigonometry may differ in the last
bits, and the tests compare to 1e-9 degrees. `densify(ellipsoid = TRUE)`
is R-side in the package and moves to stage 2 with `Geodesic`.

Planned scope: `GreatCircle` only, the `Forward` trait, path and segment-mesh input,
parent indices, canonical direction. No dependencies. The d3-geo notice
ships with it.

Done when: outputs match the R package for the README and test cases
(vertex counts exactly; coordinates to the last bit on the same
platform and libm, otherwise to 1 ulp), including the shared-edge and
antimeridian cases. The oracle is a small R script that writes fixtures
with bigcurve and PROJ; it runs once on a machine with R, and the
fixtures are committed with the forward map's output so the Rust tests
need no PROJ.

### 2. Edge models and a batched forward map

The edge model trait, `Planar`, and `Geodesic` behind a feature. The
long and negligible rules restated per model. Tests: `Planar` under an
affine map adds nothing; `Planar` under a known nonlinear map converges
at the square-root rate the R README documents (ten times tighter
tolerance, about three times the vertices).

A batched forward map alongside `Forward`: slices of points in, per-point
success out, the shape of rwarp's `Transformer::transform`. The engine
then refines breadth first: all edges that fail the test at one level are
split together, their midpoints moved in one call. Results stay identical
to the one-point engine (same tests, same order of decisions per edge),
so the R fixtures still pin it. `Forward` implements the batched map
point by point, so existing callers are unchanged.

The first consumer driving this stage is rwarp's exact footprint
(wxylib/design m19, section 5.4): destination pixel edges, straight in
destination pixel space (`Planar`), moved to source pixel space by the
warp's transformer, with the tolerance a fraction of a source cell. rwarp
transforms its pixel lattice once as scanlines and hands bigcurve only
the edges that fail the test, a row band at a time, which is why batching
matters there. Its needs, beyond `Planar` and batching:

- the forward map is any transformer (geolocation arrays, GCPs), not a
  CRS: already so, through `from_fn`;
- refinement points must go through an exact map: a map that itself
  interpolates linearly between samples (rwarp's `ApproxTransformer`)
  makes every midpoint lie on its chord, so refinement never triggers.
  Worth a sentence in the docs, since the failure is silent;
- shared edges refined once and read by both neighbours: `densify_mesh`
  already does this per segment;
- a stated rule for an edge whose midpoint fails (near a map's horizon
  or a pole): stop refining that edge, or bisect towards the boundary.
  The R engine leaves such points alone; rwarp would skip the pixel, so
  "stop" is enough for it, but the choice should be explicit;
- pure Rust, builds for wasm.

### 3. Adapters and examples

- feature `proj`: a forward map from a PROJ transformation, normalised
  to longitude, latitude order as in R;
- feature `geo-types`: densify `LineString`, `Polygon` and their multis;
- example: an analytic curvilinear grid (a polar sector) as the target
  plane, polygons moved into index space, compared with uniform
  densification at equal vertex counts.

### 4. Consumers

Kept outside this crate, listed so the API serves them:

- grid front ends for controlledburn (curvilinear and reprojected
  polygons), which would take a `Forward` and an edge model;
- rwarp's exact footprint (adaptive densify of destination pixel edges,
  see stage 2) and, opt-in, its source-window planning, where an adaptive
  boundary must still guarantee the window contains every needed pixel;
- wxylib's coverage and warp-window code;
- the R package, if and when it makes sense for it to link this crate
  rather than its own header (the header stays the R package's source
  until then; the fixtures keep the two in step).

## Open questions

- Single precision or generic float: the R engine is `f64` only; generic
  `f32` support probably is not worth its cost.
- The R default tolerance (one part in 2048 of the projected extent)
  needs the whole input first. Keep it as a helper, not an engine rule.
- Whether `Planar` edges in longitude, latitude should warn: the model
  is what GEOS assumes, but on the sphere it is a choice, and silently
  densifying "straight" lines in degrees is the trap the R vignette
  explains.
- Kinks: a forward map that is only piecewise smooth (a bilinear
  curvilinear grid has kinks at cell walls) converges to tolerance but
  never hits the kink exactly. A hook for known breakpoints (cut here
  too) may belong in the edge model.
- Name of the edge model trait (`Edge`, `Interpolant`), and whether to
  implement it directly for geo's metric spaces.
