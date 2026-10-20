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

Status: stage 1 of [the plan](design/01-plan.md), a faithful port of
the R package's engine. For the same forward map it gives bitwise
identical output to the C++ engine (the tests replay fixtures recorded
from it with PROJ; see [oracle/](oracle/README.md)). Edges are great
circles; other edge models come next.

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
