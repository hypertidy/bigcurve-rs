//! Ellipsoidal placement of new vertices (feature `geodesic`).
//!
//! The sphere decides how many vertices each edge gets and at what
//! fraction of its great circle arc; each new vertex is then moved to the
//! same fraction of the WGS84 geodesic between the same two ends. Input
//! vertices never move. This is the rule of the R package's
//! `densify(ellipsoid = TRUE)` (`R/ellipsoid.R`), with the fractions taken
//! from the bisection tree, where they are exact, instead of recovered
//! from arc angles and snapped to the dyadic grid.

use geographiclib_rs::{DirectGeodesic, Geodesic, InverseGeodesic};

pub(crate) struct Replay {
    g: Geodesic,
}

impl Replay {
    pub(crate) fn wgs84() -> Self {
        Replay {
            g: Geodesic::wgs84(),
        }
    }

    /// Move the new vertices `chain` of the written edge `a -> b` onto the
    /// geodesic. `t` holds each vertex's fraction from the canonical start
    /// (the end with the smaller longitude, then latitude), so an edge
    /// moves to the same vertices whichever way it is written. Each moved
    /// vertex keeps the longitude frame the sphere gave it.
    pub(crate) fn chain(&self, a: [f64; 2], b: [f64; 2], chain: &mut [[f64; 2]], t: &[f64]) {
        debug_assert_eq!(chain.len(), t.len());
        if chain.is_empty() {
            return;
        }
        let rev = b[0] < a[0] || (b[0] == a[0] && b[1] < a[1]);
        let (p, q) = if rev { (b, a) } else { (a, b) };
        let (s12, azi1, _azi2, _a12): (f64, f64, f64, f64) = self.g.inverse(p[1], p[0], q[1], q[0]);
        for (v, &f) in chain.iter_mut().zip(t) {
            let (lat2, lon2): (f64, f64) = self.g.direct(p[1], p[0], azi1, f * s12);
            v[0] = lon2 + 360.0 * ((v[0] - lon2) / 360.0).round();
            v[1] = lat2;
        }
    }
}
