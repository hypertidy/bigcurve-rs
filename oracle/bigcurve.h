#ifndef BIGCURVE_H
#define BIGCURVE_H

// bigcurve core: adaptive densification of segments under a nonlinear
// forward transform, after Bostock's d3-geo resampling.
//
// This header has no dependency on R or PROJ. The transform is injected
// via the Transform interface; bindings supply an adapter (see densify.cpp
// for the PROJ one). Keep it that way: this file is the piece that gets
// reused in Python (or rewritten in Rust) later.

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

namespace bigcurve {

const double DEG2RAD = 0.017453292519943295;
const double RAD2DEG = 57.29577951308232;

// Forward transform lon/lat (degrees) -> planar x/y (projected units).
// Return false if the point cannot be transformed (out of domain etc.);
// implementations should also return false for non-finite output.
struct Transform {
  virtual bool fwd(double lon, double lat, double* x, double* y) = 0;
  virtual ~Transform() {}
};

// Normalize a longitude difference to [-180, 180], with the tie at
// exactly +-180 keeping the SIGN OF THE INPUT. An edge spanning exactly
// 180 degrees is the one case where the two representatives are
// equidistant, and the input's written direction is the only correct
// tiebreaker: (0 -> 180) must stay eastward, (0 -> -180) westward.
// (std::round alone rounds ties away from zero and silently flips
// such edges.)
inline double unwrap_delta(double delta) {
  double d = delta - 360.0 * std::round(delta / 360.0);
  if (std::fabs(d) == 180.0) d = std::copysign(180.0, delta);
  return d;
}

// Spherical midpoint of the great circle arc between two lon/lat points,
// degrees in and out. The longitude REPRESENTATIVE is chosen in the frame
// of the parents: nearest to the parents' mean longitude (callers unwrap
// lon2 relative to lon1 first). This preserves each input edge's own
// convention -- an edge along lon 180 refines at 180, an edge along -180
// refines at -180, and an edge written as crossing gets a locally
// continuous chain (possibly outside [-180, 180)). The engine itself has
// no seam; the representation the input chose is respected, never
// canonicalized.
inline void gc_midpoint(double lon1, double lat1, double lon2, double lat2,
                        double* lonm, double* latm) {
  double frame_mid = (lon1 + lon2) * 0.5;
  double rlat1 = lat1 * DEG2RAD;
  double rlat2 = lat2 * DEG2RAD;
  double dlon = (lon2 - lon1) * DEG2RAD;
  double bx = std::cos(rlat2) * std::cos(dlon);
  double by = std::cos(rlat2) * std::sin(dlon);
  double cl1bx = std::cos(rlat1) + bx;
  *latm = std::atan2(std::sin(rlat1) + std::sin(rlat2),
                     std::sqrt(cl1bx * cl1bx + by * by)) * RAD2DEG;
  // degenerate: same meridian, or arc through/near a pole where the
  // azimuthal term is indeterminate -- longitude is the frame mean
  // (d3's own fallback for these cases)
  if (std::fabs(by) < 1e-14 && std::fabs(cl1bx) < 1e-14) {
    *lonm = frame_mid;
    return;
  }
  double lon = lon1 + std::atan2(by, cl1bx) * RAD2DEG;
  // snap to the representative nearest the parents' frame (ties keep
  // the computed azimuthal side)
  *lonm = frame_mid + unwrap_delta(lon - frame_mid);
}

// A vertex carries both its geographic and projected coordinates so that
// no point is ever transformed twice.
struct Vertex {
  double lon, lat;  // geographic, degrees
  double x, y;      // projected
  bool ok;          // projectable and finite
};

// cosine of the angular distance between two lon/lat points (degrees)
inline double cos_angular(double lon1, double lat1, double lon2, double lat2) {
  double rlat1 = lat1 * DEG2RAD;
  double rlat2 = lat2 * DEG2RAD;
  return std::sin(rlat1) * std::sin(rlat2) +
         std::cos(rlat1) * std::cos(rlat2) *
         std::cos((lon2 - lon1) * DEG2RAD);
}

class Resampler {
public:
  // tol: tolerance in projected (target) units. A segment is split when the
  // projected great-circle midpoint sits further than tol from the planar
  // chord, or lands parametrically far from the chord's middle (shear).
  Resampler(Transform& t, double tol, int max_depth)
    : t_(t), tol2_(tol * tol), max_depth_(max_depth) {
    // arc floor: an arc spanning less ground than tol can never need
    // interior vertices, so refinement stops there regardless of the
    // planar metric. This bounds recursion on edges that cross the
    // target projection's discontinuity (whose projected midpoints jump
    // to the far side of the map and would otherwise recurse to
    // max_depth). Assumes tol is in metres when converting to an angle;
    // for non-metric target units the floor is merely conservative.
    const double R_MEAN = 6371008.8;
    double ang = tol > 0.0 ? tol / R_MEAN : 0.0;
    cos_floor_ = std::cos(ang);
  }

  Vertex project(double lon, double lat) {
    Vertex v;
    v.lon = lon;
    v.lat = lat;
    v.ok = t_.fwd(lon, lat, &v.x, &v.y) &&
           std::isfinite(v.x) && std::isfinite(v.y);
    return v;
  }

  // Densify one segment a-b, appending interior vertices (in order, not
  // including the endpoints) to out. Segments with an unprojectable
  // endpoint are left alone; unprojectable midpoints stop refinement of
  // that branch rather than erroring. All interior vertices are emitted
  // in the longitude frame of the written start vertex a (see
  // gc_midpoint).
  //
  // The arithmetic always runs in one canonical direction, from the
  // endpoint with the lexicographically smaller (lon, lat) to the other,
  // and the chain is reversed afterwards if a-b was written the other
  // way. So a-b and b-a yield bitwise identical interior vertices (up to
  // an exact multiple of 360 in longitude when the edge crosses the
  // antimeridian), and an edge shared by two neighbouring rings refines
  // to the same vertices from both sides.
  void edge(const Vertex& a, const Vertex& b, std::vector<Vertex>& out) {
    if (!a.ok || !b.ok) return;
    bool rev = (b.lon < a.lon) || (b.lon == a.lon && b.lat < a.lat);
    const Vertex& p = rev ? b : a;
    Vertex q = rev ? a : b;
    q.lon = p.lon + unwrap_delta(q.lon - p.lon);
    std::size_t start = out.size();
    split(p, q, max_depth_, out);
    if (!rev) return;
    std::reverse(out.begin() + start, out.end());
    // p is b here and q is a unwrapped into b's frame; move the chain
    // into a's frame by the exact multiple of 360 that separates them
    double shift = 360.0 * std::round((a.lon - q.lon) / 360.0);
    if (shift != 0.0) {
      for (std::size_t k = start; k < out.size(); k++) out[k].lon += shift;
    }
  }

private:
  void split(const Vertex& a, const Vertex& b, int depth,
             std::vector<Vertex>& out) {
    if (depth <= 0) return;

    // arcs wider than ~30 degrees always split: the midpoint metric below
    // is blind to symmetric curves whose midpoint lands on the chord
    // (e.g. an oblique great circle through its own inflection), so it is
    // only allowed to adjudicate once arcs are short (d3's cosMinDistance)
    const double COS_MIN_DISTANCE = 0.8660254037844387;  // cos(30 deg)
    double cosang = cos_angular(a.lon, a.lat, b.lon, b.lat);
    // arc floor: shorter than tolerance on the ground, nothing to add
    // (also terminates on coincident/duplicate endpoints, cosang == 1)
    if (cosang >= cos_floor_) return;
    bool wide = cosang < COS_MIN_DISTANCE;

    double dx = b.x - a.x;
    double dy = b.y - a.y;
    double d2 = dx * dx + dy * dy;
    // chord shorter than 2 * tol and arc not wide: nothing left to resolve
    if (!wide && d2 <= 4.0 * tol2_) return;

    double lonm, latm;
    gc_midpoint(a.lon, a.lat, b.lon, b.lat, &lonm, &latm);
    Vertex m = project(lonm, latm);
    if (!m.ok) return;

    bool go = wide;
    if (!go && d2 > 0.0) {
      double dx2 = m.x - a.x;
      double dy2 = m.y - a.y;
      double dz = dx * dy2 - dy * dx2;           // 2 * signed area
      double perp2 = (dz * dz) / d2;             // squared perpendicular dist
      double tpar = (dx * dx2 + dy * dy2) / d2;  // parametric position
      go = perp2 > tol2_ || tpar < 0.2 || tpar > 0.8;
    }

    if (go) {
      split(a, m, depth - 1, out);
      out.push_back(m);
      split(m, b, depth - 1, out);
    }
  }

  Transform& t_;
  double tol2_;
  int max_depth_;
  double cos_floor_;
};

}  // namespace bigcurve

#endif
