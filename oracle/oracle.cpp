// Fixture oracle for the Rust port.
//
// Runs the C++ engine of the R package (bigcurve.h, copied unchanged from
// hypertidy/bigcurve) through the same path and mesh drivers as the R
// package's src/densify.cpp, with the same PROJ adapter, and writes one
// fixture file per case. Every forward transform call is recorded, so the
// Rust tests replay PROJ's answers without linking PROJ: if the Rust
// engine ever asks for a point the C++ engine did not, the replay fails.
//
// Build and run with oracle/run.sh; see oracle/README.md.

#include <algorithm>
#include <cinttypes>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <map>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include <proj.h>

#include "bigcurve.h"

namespace {

std::uint64_t bits(double v) {
  std::uint64_t u;
  std::memcpy(&u, &v, sizeof u);
  return u;
}

struct Call {
  bool ok;
  double x, y;
};

// The R package's ProjTransform, plus a record of every call.
class RecordingProj : public bigcurve::Transform {
public:
  RecordingProj(const std::string& from, const std::string& to) {
    PJ_CONTEXT* ctx = PJ_DEFAULT_CTX;
    PJ* raw = proj_create_crs_to_crs(ctx, from.c_str(), to.c_str(), NULL);
    if (raw == NULL) throw std::runtime_error("cannot create " + to);
    PJ* norm = proj_normalize_for_visualization(ctx, raw);
    if (norm != NULL) {
      proj_destroy(raw);
      trans_ = norm;
    } else {
      trans_ = raw;
    }
  }
  ~RecordingProj() { proj_destroy(trans_); }

  bool fwd(double lon, double lat, double* x, double* y) override {
    PJ_COORD in;
    in.xyzt.x = lon;
    in.xyzt.y = lat;
    in.xyzt.z = 0.0;
    in.xyzt.t = HUGE_VAL;
    PJ_COORD out = proj_trans(trans_, PJ_FWD, in);
    *x = out.xy.x;
    *y = out.xy.y;
    bool ok = std::isfinite(*x) && std::isfinite(*y);
    calls[{bits(lon), bits(lat)}] = Call{ok, *x, *y};
    return ok;
  }

  std::map<std::pair<std::uint64_t, std::uint64_t>, Call> calls;

private:
  PJ* trans_;
};

struct Out {
  std::vector<double> x, y;
  std::vector<int> s0, s1, parent;
};

// densify_path_cpp from the R package, without R
Out path(bigcurve::Transform& t, const std::vector<double>& lon,
         const std::vector<double>& lat, double tol, int max_depth) {
  bigcurve::Resampler rs(t, tol, max_depth);
  std::size_t n = lon.size();
  std::vector<bigcurve::Vertex> v(n);
  for (std::size_t i = 0; i < n; i++) v[i] = rs.project(lon[i], lat[i]);
  Out o;
  std::vector<bigcurve::Vertex> interior;
  for (std::size_t i = 0; i + 1 < n; i++) {
    o.x.push_back(v[i].lon);
    o.y.push_back(v[i].lat);
    interior.clear();
    rs.edge(v[i], v[i + 1], interior);
    for (auto& w : interior) {
      o.x.push_back(w.lon);
      o.y.push_back(w.lat);
    }
  }
  if (n > 0) {
    o.x.push_back(v[n - 1].lon);
    o.y.push_back(v[n - 1].lat);
  }
  return o;
}

// densify_mesh_cpp from the R package, without R
Out mesh(bigcurve::Transform& t, const std::vector<double>& x,
         const std::vector<double>& y, const std::vector<int>& s0,
         const std::vector<int>& s1, double tol, int max_depth) {
  bigcurve::Resampler rs(t, tol, max_depth);
  int nv = (int)x.size();
  std::vector<bigcurve::Vertex> v((std::size_t)nv);
  for (int i = 0; i < nv; i++) v[(std::size_t)i] = rs.project(x[(std::size_t)i], y[(std::size_t)i]);
  Out o;
  o.x = x;
  o.y = y;
  std::vector<bigcurve::Vertex> interior;
  struct Seen { int from; int first; int n; };
  std::unordered_map<long long, Seen> seen;
  std::vector<int> chain;
  for (std::size_t k = 0; k < s0.size(); k++) {
    int a = s0[k];
    int b = s1[k];
    long long key = (long long)std::min(a, b) * (long long)nv +
                    (long long)std::max(a, b);
    chain.clear();
    auto hit = seen.find(key);
    if (hit != seen.end()) {
      const Seen& h = hit->second;
      for (int j = 0; j < h.n; j++) chain.push_back(h.first + j);
      if (h.from != a) std::reverse(chain.begin(), chain.end());
    } else {
      interior.clear();
      rs.edge(v[(std::size_t)a], v[(std::size_t)b], interior);
      int first = (int)o.x.size();
      for (std::size_t j = 0; j < interior.size(); j++) {
        o.x.push_back(interior[j].lon);
        o.y.push_back(interior[j].lat);
        chain.push_back(first + (int)j);
      }
      seen.emplace(key, Seen{a, first, (int)interior.size()});
    }
    int prev = a;
    for (int idx : chain) {
      o.s0.push_back(prev);
      o.s1.push_back(idx);
      o.parent.push_back((int)k);
      prev = idx;
    }
    o.s0.push_back(prev);
    o.s1.push_back(b);
    o.parent.push_back((int)k);
  }
  return o;
}

// default_tolerance from the R package's densify.R: one part in 'pixels'
// of the larger side of the projected extent of the input vertices
double default_tolerance(bigcurve::Transform& t, const std::vector<double>& lon,
                         const std::vector<double>& lat, double pixels) {
  double xmin = INFINITY, xmax = -INFINITY, ymin = INFINITY, ymax = -INFINITY;
  for (std::size_t i = 0; i < lon.size(); i++) {
    double x, y;
    if (!t.fwd(lon[i], lat[i], &x, &y)) continue;
    xmin = std::min(xmin, x);
    xmax = std::max(xmax, x);
    ymin = std::min(ymin, y);
    ymax = std::max(ymax, y);
  }
  return std::max(xmax - xmin, ymax - ymin) / pixels;
}

struct Case {
  std::string name;
  std::string target;
  std::vector<double> lon, lat;
  std::vector<int> s0, s1;  // empty for a path
  double tol;               // NAN: use the default with 'pixels'
  int max_depth;
  double pixels;
};

void write(const std::string& dir, const Case& c) {
  RecordingProj t("EPSG:4326", c.target);
  double tol = c.tol;
  bool dflt = std::isnan(tol);
  if (dflt) tol = default_tolerance(t, c.lon, c.lat, c.pixels);
  bool is_mesh = !c.s0.empty();
  Out o = is_mesh ? mesh(t, c.lon, c.lat, c.s0, c.s1, tol, c.max_depth)
                  : path(t, c.lon, c.lat, tol, c.max_depth);

  std::string file = dir + "/" + c.name + ".txt";
  FILE* f = std::fopen(file.c_str(), "w");
  if (!f) throw std::runtime_error("cannot write " + file);
  std::fprintf(f, "# bigcurve fixture: %s\n", c.name.c_str());
  std::fprintf(f, "# target %s\n", c.target.c_str());
  std::fprintf(f, "# doubles are IEEE 754 bit patterns in hex\n");
  std::fprintf(f, "kind %s\n", is_mesh ? "mesh" : "path");
  if (dflt) {
    std::fprintf(f, "pixels %016" PRIx64 "\n", bits(c.pixels));
  }
  std::fprintf(f, "tolerance %016" PRIx64 "\n", bits(tol));
  std::fprintf(f, "max_depth %d\n", c.max_depth);
  for (std::size_t i = 0; i < c.lon.size(); i++) {
    std::fprintf(f, "in %016" PRIx64 " %016" PRIx64 "\n", bits(c.lon[i]), bits(c.lat[i]));
  }
  for (std::size_t i = 0; i < c.s0.size(); i++) {
    std::fprintf(f, "seg %d %d\n", c.s0[i], c.s1[i]);
  }
  for (auto& kv : t.calls) {
    std::fprintf(f, "fwd %016" PRIx64 " %016" PRIx64 " %s %016" PRIx64 " %016" PRIx64 "\n",
                 kv.first.first, kv.first.second, kv.second.ok ? "ok" : "na",
                 bits(kv.second.x), bits(kv.second.y));
  }
  for (std::size_t i = 0; i < o.x.size(); i++) {
    std::fprintf(f, "out %016" PRIx64 " %016" PRIx64 "\n", bits(o.x[i]), bits(o.y[i]));
  }
  for (std::size_t i = 0; i < o.s0.size(); i++) {
    std::fprintf(f, "oseg %d %d %d\n", o.s0[i], o.s1[i], o.parent[i]);
  }
  std::fclose(f);
  std::printf("%-28s %-40s tol %-12.6g out %zu calls %zu\n", c.name.c_str(),
              c.target.c_str(), tol, o.x.size(), t.calls.size());
}

Case seg(const std::string& name, const std::string& target, double lon0,
         double lon1, double lat0, double lat1, double tol, int depth = 16) {
  return Case{name, target, {lon0, lon1}, {lat0, lat1}, {}, {}, tol, depth, 2048.0};
}

}  // namespace

int main(int argc, char** argv) {
  std::string dir = argc > 1 ? argv[1] : "tests/fixtures";
  const double D = NAN;  // default tolerance
  const std::string tmerc = "+proj=tmerc +lon_0=147 +lat_0=-42";
  const std::string laea147 = "+proj=laea +lon_0=147 +lat_0=-42";
  const std::string laea50 = "+proj=laea +lon_0=50 +lat_0=30";
  const std::string merc = "+proj=merc";

  std::vector<Case> cases;

  // the R package's test-densify.R, case by case
  cases.push_back(seg("crossing_tmerc_default", tmerc, 140, -164, -89, 80, D));
  cases.push_back(seg("equator_laea_straight", "+proj=laea +lon_0=0 +lat_0=0", -10, 10, 0, 0, 1));
  cases.push_back(seg("crossing_laea_1e5", laea147, 140, -164, -89, 80, 1e5));
  cases.push_back(seg("crossing_laea_1e3", laea147, 140, -164, -89, 80, 1e3));
  cases.push_back(seg("crossing_laea_depth2", laea147, 140, -164, -89, 80, 1, 2));
  cases.push_back(seg("ortho_far_side", "+proj=ortho +lon_0=-170 +lat_0=0", 0, 10, 0, 10, 1000));
  cases.push_back(seg("diagonal_forward", laea50, 0, 100, 0, 50, 2e4));
  cases.push_back(seg("diagonal_reverse", laea50, 100, 0, 50, 0, 2e4));
  cases.push_back(Case{"mesh_repeated_edge", laea50, {0, 100}, {0, 50}, {0, 1}, {1, 0}, 2e4, 16, 2048});
  cases.push_back(Case{"mesh_parent", laea147, {140, -164, -89}, {-89, 80, 80}, {0, 1}, {1, 2}, 1e4, 16, 2048});
  cases.push_back(seg("merc_along_180", merc, 180, 180, -60, 60, 5e4));
  cases.push_back(seg("merc_along_m180", merc, -180, -180, -60, 60, 5e4));
  cases.push_back(seg("merc_written_crossing", merc, 170, -170, 20, 20, 5e4));
  cases.push_back(seg("merc_short_seam", merc, 179.9, -179.9, 0, 0, 5e4));
  cases.push_back(seg("merc_equator_e1", merc, 0, 180, 0, 0, 5e4));
  cases.push_back(seg("merc_equator_e2", merc, 180, 0, 0, 0, 5e4));
  cases.push_back(seg("merc_equator_w1", merc, -180, 0, 0, 0, 5e4));
  cases.push_back(seg("merc_equator_w2", merc, 0, -180, 0, 0, 5e4));
  cases.push_back(seg("parallel_default", merc, 0, 100, 60, 60, D));
  cases.push_back(seg("meridian_default", "+proj=laea +lon_0=150", 180, 180, -60, 60, D));
  cases.push_back(Case{"neighbours_shared_diagonal", laea50,
                       {0, 100, 0, 100}, {0, 50, 60, 0},
                       {0, 1, 2, 0, 3, 1}, {1, 2, 0, 3, 1, 0}, 2e4, 16, 2048});

  // the README: one route in several projections, and its tolerance ladder
  cases.push_back(seg("hobart_london_eqearth", "+proj=eqearth", 147.3, -0.1, -42.9, 51.5, D));
  cases.push_back(seg("hobart_london_laea", "+proj=laea +lon_0=70 +lat_0=10", 147.3, -0.1, -42.9, 51.5, D));
  // in a gnomonic projection every great circle is straight: only the
  // forced splits of a wide arc are added
  cases.push_back(seg("hobart_kolkata_gnom", "+proj=gnom +lon_0=100 +lat_0=0", 147.3, 88.4, -42.9, 22.6, D));
  cases.push_back(seg("hobart_kolkata_ortho", "+proj=ortho +lon_0=100 +lat_0=0", 147.3, 88.4, -42.9, 22.6, D));
  // one end beyond the horizon: the segment is left alone
  cases.push_back(seg("hobart_london_ortho_horizon", "+proj=ortho +lon_0=100 +lat_0=0", 147.3, -0.1, -42.9, 51.5, 1000));
  const double tols[] = {1e5, 1e4, 1e3, 1e2};
  const char* tname[] = {"1e5", "1e4", "1e3", "1e2"};
  for (int i = 0; i < 4; i++) {
    cases.push_back(seg(std::string("hobart_london_eqearth_") + tname[i],
                        "+proj=eqearth", 147.3, -0.1, -42.9, 51.5, tols[i]));
  }
  // a multi-vertex path, the README hero's routes joined end to end
  cases.push_back(Case{"routes_path_eqearth", "+proj=eqearth",
                       {147.3, -0.1, -70.7, 28.0, -74.0, 37.6, 115.9},
                       {-42.9, 51.5, -33.4, -26.2, 40.7, 55.8, -31.9},
                       {}, {}, D, 16, 2048});

  // a southern graticule (30 x 15 degree cells), merged: shared edges
  // appear twice in opposite directions and must be refined once
  {
    Case g{"graticule_merged", "+proj=laea +lat_0=-90", {}, {}, {}, {}, D, 16, 2048};
    const int nx = 12, ny = 4;
    auto id = [&](int i, int j) { return j * (nx + 1) + i; };
    for (int j = 0; j <= ny; j++) {
      for (int i = 0; i <= nx; i++) {
        g.lon.push_back(-180.0 + 30.0 * i);
        g.lat.push_back(-75.0 + 15.0 * j);
      }
    }
    for (int j = 0; j < ny; j++) {
      for (int i = 0; i < nx; i++) {
        int ring[5] = {id(i, j), id(i + 1, j), id(i + 1, j + 1), id(i, j + 1), id(i, j)};
        for (int k = 0; k < 4; k++) {
          g.s0.push_back(ring[k]);
          g.s1.push_back(ring[k + 1]);
        }
      }
    }
    cases.push_back(g);
  }

  for (auto& c : cases) write(dir, c);
  return 0;
}
