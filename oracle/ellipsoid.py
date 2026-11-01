"""Ellipsoid fixtures: the R package's `densify(ellipsoid = TRUE)` replay.

Reads each sphere fixture in tests/fixtures (the C++ engine's output) and
applies the rule of R/ellipsoid.R in hypertidy/bigcurve, line for line:
each new vertex's fraction t along its input edge is the ratio of unit
sphere arc angles from the canonical start, snapped to the dyadic grid
2^min(max_depth, 50); the vertex moves to that fraction of the WGS84
geodesic between the edge's ends; its longitude keeps the sphere's frame.

GeographicLib's geodesic routines come from the Python package
geographiclib (Karney's own port; R's geographiclib wraps the C++ one).
Writes tests/fixtures/ellipsoid/<name>.txt with the moved vertices.

    python3 -m pip install geographiclib==2.1
    python3 oracle/ellipsoid.py
"""

import math
import os
import struct
import sys

from geographiclib.geodesic import Geodesic

G = Geodesic.WGS84


def unhex(s):
    return struct.unpack(">d", bytes.fromhex(s))[0]


def tohex(x):
    return struct.pack(">d", x).hex()


def read(path):
    f = {"input": [], "segments": [], "output": [], "oseg": []}
    with open(path) as fh:
        for line in fh:
            w = line.split()
            if not w or w[0] == "#":
                continue
            if w[0] == "kind":
                f["kind"] = w[1]
            elif w[0] == "max_depth":
                f["max_depth"] = int(w[1])
            elif w[0] == "in":
                f["input"].append((unhex(w[1]), unhex(w[2])))
            elif w[0] == "seg":
                f["segments"].append((int(w[1]), int(w[2])))
            elif w[0] == "out":
                f["output"].append((unhex(w[1]), unhex(w[2])))
            elif w[0] == "oseg":
                f["oseg"].append((int(w[1]), int(w[2]), int(w[3])))
    return f


def arc_angle(lon1, lat1, lon2, lat2):
    r = math.pi / 180
    x1 = math.cos(lat1 * r) * math.cos(lon1 * r)
    y1 = math.cos(lat1 * r) * math.sin(lon1 * r)
    z1 = math.sin(lat1 * r)
    x2 = math.cos(lat2 * r) * math.cos(lon2 * r)
    y2 = math.cos(lat2 * r) * math.sin(lon2 * r)
    z2 = math.sin(lat2 * r)
    cx = y1 * z2 - z1 * y2
    cy = z1 * x2 - x1 * z2
    cz = x1 * y2 - y1 * x2
    return math.atan2(math.sqrt(cx * cx + cy * cy + cz * cz), x1 * x2 + y1 * y2 + z1 * z2)


def move(a, b, v, max_depth):
    """R's ellipsoid_replay for one new vertex v of input edge a -> b."""
    ax, ay = a
    bx, by = b
    rev = bx < ax or (bx == ax and by < ay)
    px, py = (bx, by) if rev else (ax, ay)
    qx, qy = (ax, ay) if rev else (bx, by)
    vx, vy = v
    grid = 2.0 ** min(max_depth, 50)
    t = round(arc_angle(px, py, vx, vy) / arc_angle(px, py, qx, qy) * grid) / grid
    inv = G.Inverse(py, px, qy, qx)
    d = G.Direct(py, px, inv["azi1"], t * inv["s12"])
    lon = d["lon2"] + 360 * round((vx - d["lon2"]) / 360)
    return (lon, d["lat2"])


def mesh_vertices(f):
    """New vertices moved; each one starts an output segment whose parent
    is its input segment (R: k <- match(new_id, out$s0))."""
    verts = list(f["output"])
    n_old = len(f["input"])
    starts = {}
    for s0, _s1, parent in f["oseg"]:
        starts.setdefault(s0, parent)
    for i in range(n_old, len(verts)):
        a, b = f["segments"][starts[i]]
        verts[i] = move(verts[a], verts[b], verts[i], f["max_depth"])
    return verts


def path_vertices(f):
    """A path is a mesh of consecutive segments; walk the output, matching
    input vertices in order."""
    inp, out = f["input"], list(f["output"])
    j = 0
    k = -1
    for i, v in enumerate(out):
        if j < len(inp) and v == inp[j]:
            k = j
            j += 1
            continue
        out[i] = move(inp[k], inp[k + 1], v, f["max_depth"])
    assert j == len(inp), "every input vertex is in the output"
    return out


def main():
    root = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "tests", "fixtures")
    dest = os.path.join(root, "ellipsoid")
    os.makedirs(dest, exist_ok=True)
    for old in os.listdir(dest):
        if old.endswith(".txt"):
            os.remove(os.path.join(dest, old))
    n = 0
    for name in sorted(os.listdir(root)):
        if not name.endswith(".txt"):
            continue
        f = read(os.path.join(root, name))
        verts = mesh_vertices(f) if f["kind"] == "mesh" else path_vertices(f)
        if len(f["output"]) == len(f["input"]) and f["kind"] == "path":
            continue  # nothing was added, nothing moves
        with open(os.path.join(dest, name), "w") as fh:
            fh.write("# bigcurve ellipsoid fixture: %s\n" % name[:-4])
            fh.write("# R/ellipsoid.R replay with geographiclib (Python) %s\n" % getattr(
                sys.modules["geographiclib"], "__version__", "?"))
            for lon, lat in verts:
                fh.write("eout %s %s\n" % (tohex(lon), tohex(lat)))
        n += 1
    print("wrote %d ellipsoid fixtures to %s" % (n, os.path.relpath(dest)))


if __name__ == "__main__":
    import geographiclib  # noqa: F401  (for the version line)
    main()
