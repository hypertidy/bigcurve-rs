//! Edge models: what the middle of an edge is in the space it was drawn in.

use crate::{DEG2RAD, RAD2DEG};

/// Mean Earth radius in metres, used to turn the tolerance into an arc
/// floor for [`GreatCircle`].
const R_MEAN: f64 = 6371008.8;

/// cos(30 degrees): arcs wider than this always split (d3-geo's
/// `cosMinDistance`).
const COS_MIN_DISTANCE: f64 = 0.8660254037844387;

/// How an edge model classes an edge before the midpoint test.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Span {
    /// Too short to ever need interior vertices: refinement stops,
    /// whatever the target plane says. This bounds recursion on edges
    /// across a map's discontinuity, whose moved midpoints jump away.
    Negligible,
    /// Decided by the midpoint test in the target plane.
    Normal,
    /// Long enough that the edge always splits: the midpoint test is
    /// blind to symmetric curves whose moved midpoint lands on the chord,
    /// so it only adjudicates once edges are short.
    Long,
}

/// An edge model: the shape of an edge between two points in the source
/// space, as far as the engine needs to know it.
///
/// Every method sees an edge in its canonical direction (the end with the
/// lexicographically smaller coordinates first) with the far end already
/// [`align`](Edge::align)ed to the near end, so a model that is
/// deterministic gives bitwise identical vertices for `a -> b` and
/// `b -> a`.
pub trait Edge {
    /// The midpoint of the edge `a -> b`, in source coordinates.
    fn mid(&self, a: [f64; 2], b: [f64; 2]) -> [f64; 2];

    /// Class the edge `a -> b` before its midpoint is moved. `floor` is
    /// what [`floor`](Edge::floor) returned for the tolerance in use.
    fn span(&self, a: [f64; 2], b: [f64; 2], floor: f64) -> Span;

    /// A per-tolerance constant for [`span`](Edge::span), computed once
    /// per engine. The default is 0.
    fn floor(&self, tolerance: f64) -> f64 {
        let _ = tolerance;
        0.0
    }

    /// `to` expressed in the frame of `from`, for spaces with more than
    /// one representative of a point (longitudes differing by 360). The
    /// default is `to` unchanged.
    fn align(&self, from: [f64; 2], to: [f64; 2]) -> [f64; 2] {
        let _ = from;
        to
    }

    /// The offset that moves a vertex refined in the frame of `framed`
    /// into the frame of `written`, two representatives of one point. The
    /// default is no offset.
    fn reframe(&self, written: [f64; 2], framed: [f64; 2]) -> [f64; 2] {
        let _ = (written, framed);
        [0.0, 0.0]
    }
}

impl<E: Edge + ?Sized> Edge for &E {
    fn mid(&self, a: [f64; 2], b: [f64; 2]) -> [f64; 2] {
        (**self).mid(a, b)
    }
    fn span(&self, a: [f64; 2], b: [f64; 2], floor: f64) -> Span {
        (**self).span(a, b, floor)
    }
    fn floor(&self, tolerance: f64) -> f64 {
        (**self).floor(tolerance)
    }
    fn align(&self, from: [f64; 2], to: [f64; 2]) -> [f64; 2] {
        (**self).align(from, to)
    }
    fn reframe(&self, written: [f64; 2], framed: [f64; 2]) -> [f64; 2] {
        (**self).reframe(written, framed)
    }
}

/// Great circle arcs between longitude, latitude points in degrees: the R
/// package's model and the default.
///
/// - The midpoint is the spherical midpoint, in the longitude frame of
///   the edge's start (no seam: an edge along 180 refines at 180, along
///   -180 at -180).
/// - Arcs wider than 30 degrees always split ([`Span::Long`]).
/// - Arcs spanning less ground than the tolerance never split
///   ([`Span::Negligible`]). The conversion to an angle assumes the
///   tolerance is in metres on a sphere of mean Earth radius; for other
///   target units the floor is merely conservative.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct GreatCircle;

impl Edge for GreatCircle {
    fn mid(&self, a: [f64; 2], b: [f64; 2]) -> [f64; 2] {
        let (lon, lat) = gc_midpoint(a[0], a[1], b[0], b[1]);
        [lon, lat]
    }

    fn span(&self, a: [f64; 2], b: [f64; 2], cos_floor: f64) -> Span {
        let cosang = cos_angular(a[0], a[1], b[0], b[1]);
        // also ends coincident endpoints, cosang == 1
        if cosang >= cos_floor {
            Span::Negligible
        } else if cosang < COS_MIN_DISTANCE {
            Span::Long
        } else {
            Span::Normal
        }
    }

    fn floor(&self, tolerance: f64) -> f64 {
        let ang = if tolerance > 0.0 {
            tolerance / R_MEAN
        } else {
            0.0
        };
        ang.cos()
    }

    fn align(&self, from: [f64; 2], to: [f64; 2]) -> [f64; 2] {
        [from[0] + unwrap_delta(to[0] - from[0]), to[1]]
    }

    fn reframe(&self, written: [f64; 2], framed: [f64; 2]) -> [f64; 2] {
        [360.0 * ((written[0] - framed[0]) / 360.0).round(), 0.0]
    }
}

/// Straight lines in the source coordinates, as GEOS, GDAL and most file
/// formats read edges, and as a destination pixel's sides are in pixel
/// space.
///
/// - The midpoint is the arithmetic mean of the ends.
/// - An edge at least [`long`](Planar::long) long (source units) always
///   splits; the default is never.
/// - An edge no longer than [`floor`](Planar::floor) (source units) never
///   splits; the default 0 stops only at coincident ends, so refinement
///   across a discontinuity is bounded by the maximum depth alone.
///
/// Planar edges in longitude and latitude are a choice, not a default:
/// "straight in degrees" is a rhumb-like line that no great circle
/// follows.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Planar {
    /// Edges no longer than this (source units) are negligible.
    pub floor: f64,
    /// Edges at least this long (source units) always split.
    pub long: f64,
}

impl Default for Planar {
    fn default() -> Self {
        Planar {
            floor: 0.0,
            long: f64::INFINITY,
        }
    }
}

impl Planar {
    /// Planar edges with no floor and no forced splits.
    pub fn new() -> Self {
        Self::default()
    }

    /// The same model with another floor (source units).
    pub fn with_floor(self, floor: f64) -> Self {
        Planar { floor, ..self }
    }

    /// The same model with another length above which edges always split
    /// (source units).
    pub fn with_long(self, long: f64) -> Self {
        Planar { long, ..self }
    }
}

impl Edge for Planar {
    fn mid(&self, a: [f64; 2], b: [f64; 2]) -> [f64; 2] {
        [(a[0] + b[0]) * 0.5, (a[1] + b[1]) * 0.5]
    }

    fn span(&self, a: [f64; 2], b: [f64; 2], _floor: f64) -> Span {
        let (dx, dy) = (b[0] - a[0], b[1] - a[1]);
        let d2 = dx * dx + dy * dy;
        if d2 <= self.floor * self.floor {
            Span::Negligible
        } else if d2 >= self.long * self.long {
            Span::Long
        } else {
            Span::Normal
        }
    }
}

/// Normalise a longitude difference to [-180, 180], with the tie at
/// exactly +-180 keeping the sign of the input. An edge spanning exactly
/// 180 degrees is the one case where the two representatives are
/// equidistant, and the written direction is the only correct tiebreaker:
/// (0 -> 180) stays eastward, (0 -> -180) westward.
pub(crate) fn unwrap_delta(delta: f64) -> f64 {
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
pub(crate) fn gc_midpoint(lon1: f64, lat1: f64, lon2: f64, lat2: f64) -> (f64, f64) {
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
    fn planar_spans() {
        let p = Planar::new().with_floor(1.0).with_long(10.0);
        assert_eq!(p.span([0.0, 0.0], [0.5, 0.5], 0.0), Span::Negligible);
        assert_eq!(p.span([0.0, 0.0], [5.0, 0.0], 0.0), Span::Normal);
        assert_eq!(p.span([0.0, 0.0], [10.0, 0.0], 0.0), Span::Long);
        assert_eq!(
            Planar::new().span([1.0, 1.0], [1.0, 1.0], 0.0),
            Span::Negligible
        );
        assert_eq!(
            Planar::new().span([0.0, 0.0], [1e100, 0.0], 0.0),
            Span::Normal
        );
    }
}
