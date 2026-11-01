# Changelog

## 0.3.0 (unreleased)

- This crate is now the canonical bigcurve engine; the R package's C++
  header is frozen (hypertidy/bigcurve#14).
- Feature `geodesic`: `densify_path_ellipsoid` and
  `densify_mesh_ellipsoid`, the R package's `densify(ellipsoid = TRUE)`.
  The sphere chooses how many vertices each edge gets and at what fraction
  of its arc; each new vertex moves to the same fraction of the WGS84
  geodesic (`geographiclib-rs`, pure Rust). Fractions come from the
  bisection tree, so they are exact. Input vertices never move, an edge
  gives the same vertices whichever way it is written, shared segments
  stay shared, and each vertex keeps its longitude frame.
- Fixtures: oracle/ellipsoid.py applies R/ellipsoid.R's replay to the C++
  engine's output with Python geographiclib (tests/fixtures/ellipsoid);
  the Rust drivers agree to 4e-14 degrees.

## 0.2.0 (2026-11-01)

- Edge models: the `Edge` trait, with `GreatCircle` (the stage 1 engine,
  still the default) and `Planar` (straight lines in the source
  coordinates; optional floor and forced-split lengths in source units).
  `densify_path_with`, `densify_mesh_with` and `Resampler::with_edge` take
  a model; `densify_path`, `densify_mesh` and `Resampler::new` keep great
  circles.
- A batched forward map, `ForwardBatch`: coordinate arrays moved in place
  with per-point success and a `near` hint per point (the target position
  of the start of the edge being refined). Every `Forward` is one, point
  by point.
- `Resampler::edges` refines many edges level by level, one forward call
  per level; the drivers use it. The vertices are bitwise those of the
  depth-first `Resampler::edge` (checked against the R fixtures and on 300
  random edges); only the order of forward calls changes.
- `Resampler::project_all`, and `Vertex::new` for edges whose ends the
  caller moved.

## 0.1.0

- Faithful port of the R package's engine (great circle edges), bitwise
  against fixtures recorded from the C++ header with PROJ.
