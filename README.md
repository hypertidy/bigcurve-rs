# bigcurve

A line between two points is only straight in the space it was drawn in.
Move it into another space (a map projection, another CRS, the index
space of a model grid) and its true image is a curve. Joining the moved
ends draws the wrong line. bigcurve adds the vertices a line needs to
keep its true shape in the target space: more where the map bends it,
none where it doesn't.

This is the Rust home of [bigcurve](https://github.com/hypertidy/bigcurve),
the R package. The method is the adaptive resampling Mike Bostock wrote
for [d3-geo](https://github.com/d3/d3-geo): bisect an edge while its
moved midpoint sits off the chord of its moved ends, with the error
measured where it matters, in the target plane.

Status: stage 2 of [the plan](design/01-plan.md). Stage 1 is a faithful
port of the R package's engine: for the same forward map it gives bitwise
identical output to the C++ engine (the tests replay fixtures recorded
from it with PROJ; see [oracle/](oracle/README.md)). Stage 2 adds the
edge model as a parameter (`GreatCircle`, the default, and `Planar`) and
a batched forward map (`ForwardBatch`): refinement runs level by level,
one forward call per level, with bitwise the same vertices as the
depth-first recursion. The geodesic edge model is still to come.

```rust
use bigcurve::{densify_path, from_fn, Options};

// any forward map: a PROJ transformation, an analytic projection,
// or positions to the index space of a grid
let merc = from_fn(|lon: f64, lat: f64| {
    let y = (std::f64::consts::FRAC_PI_4 + lat.to_radians() / 2.0).tan().ln();
    Some([6378137.0 * lon.to_radians(), 6378137.0 * y])
});
let line = [[0.0, 60.0], [100.0, 60.0]];
let d = densify_path(merc, &line, &Options::new(5e4)); // tolerance in metres
```

`densify_mesh` takes vertices and segment index pairs, refines a
segment that appears twice (a shared boundary) once, and reports each
new segment's parent; `default_tolerance` gives the R package's "one
part in 2048 of the projected extent".

Planar edges, a batched map, and a tolerance in the target's own units
(here a thousandth of a cell of a grid's index space):

```rust
use bigcurve::{densify_path_with, ForwardBatch, Options, Planar};

struct ToIndex; // stands in for a grid's world-to-index map
impl ForwardBatch for ToIndex {
    fn forward_batch(&mut self, x: &mut [f64], y: &mut [f64], _near: &[[f64; 2]], ok: &mut [bool]) {
        for i in 0..x.len() {
            let (r, a) = (x[i].hypot(y[i]) / 10.0, y[i].atan2(x[i]).to_degrees());
            (x[i], y[i], ok[i]) = (r, a, true);
        }
    }
}
let d = densify_path_with(ToIndex, Planar::new(), &[[100.0, 0.0], [0.0, 100.0]], &Options::new(1e-3));
```

`near` carries the target position of the edge being refined, so a map
whose inverse is a search (locating a point in a curvilinear grid) can
start next to the answer. Use an exact map: one that itself interpolates
between samples puts every midpoint on its chord, and nothing is added.

## The idea in one paragraph

Densifying is three choices and a loop. An **edge model** says what the
middle of an edge is where the line was drawn (a great circle on the
sphere, a straight line in a projected CRS, a geodesic on the
ellipsoid). A **forward map** takes a position to the target space, and
may fail. A **tolerance** says how far, in target units, a chord may
stray from the true curve. The loop bisects until every chord is within
tolerance, transforming each point once. Edges shared by neighbouring
polygons refine to identical vertices, so a set of shapes that tiled the
plane before still tiles it after.

## Shape

- a core with no dependencies: the edge model and forward map are traits
  you implement or pick from the built-in ones;
- inputs as coordinate slices and segment meshes (a vertex pool plus
  index pairs), with provenance from each new segment to its parent;
- optional adapters: PROJ for map projections, geodesics on the
  ellipsoid, `geo-types` geometry.

## Licence

MIT. The refinement test follows d3-geo's `resampleLineTo` (ISC,
Copyright 2010-2024 Mike Bostock); see [NOTICE](NOTICE).
