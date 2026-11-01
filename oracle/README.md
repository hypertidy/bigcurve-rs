# Fixture oracle

`oracle.cpp` runs the C++ engine of the R package bigcurve through the
same path and mesh drivers as the R package, with the same PROJ adapter,
and writes `tests/fixtures/*.txt`. The Rust tests replay the recorded PROJ
calls, so they need neither a C++ compiler nor PROJ.

- `bigcurve.h` is copied unchanged from hypertidy/bigcurve at commit
  402da30e16423f551b539f7db527d3ba4c558d0a (src/bigcurve.h, 0.6.0).
- The cases follow the R package's tests and README: the antimeridian
  and exactly-180-degree rules, unprojectable ends, the shared edge and
  repeated mesh edge, parent indices, the default tolerance, and a
  merged southern graticule.

Regenerate from the repo root (needs a C++ compiler and libproj-dev):

    ./oracle/run.sh

CI regenerates the fixtures and fails if they differ from the committed
ones, so a change to the header or to PROJ's answers shows up there.

## Fixture format

One record per line; doubles are IEEE 754 bit patterns in hex, so
comparisons are exact.

    kind path|mesh
    pixels <hex>             only when the tolerance is the default
    tolerance <hex>
    max_depth <n>
    in <lon> <lat>           input vertices
    seg <a> <b>              input segments (mesh)
    fwd <lon> <lat> ok|na <x> <y>   every forward call the engine made
    out <lon> <lat>          output vertices
    oseg <a> <b> <parent>    output segments (mesh)

## Ellipsoid fixtures

`ellipsoid.py` applies the R package's `densify(ellipsoid = TRUE)` rule
(`R/ellipsoid.R`, line for line) to the C++ engine's output in
`tests/fixtures`, using Python geographiclib for the geodesics, and
writes `tests/fixtures/ellipsoid/*.txt` (`eout <lon> <lat>` per output
vertex, doubles in hex). `run.sh` runs it after the C++ oracle; it needs
`python3 -m pip install geographiclib==2.1`. The Rust drivers are held
to 1e-9 degrees against these, since the two sides use different ports
of GeographicLib.
