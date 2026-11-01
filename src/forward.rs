//! Forward maps: one point at a time, or in batches.

/// A forward map from source coordinates to the target plane.
///
/// Source coordinates are whatever the [`Edge`](crate::Edge) model reads:
/// longitude, latitude in degrees for [`GreatCircle`](crate::GreatCircle),
/// any plane for [`Planar`](crate::Planar).
///
/// Return `None` for a point that cannot be moved (beyond a horizon, out
/// of the map's domain). Non-finite output is treated the same way.
pub trait Forward {
    fn forward(&mut self, lon: f64, lat: f64) -> Option<[f64; 2]>;
}

impl<T: Forward + ?Sized> Forward for &mut T {
    fn forward(&mut self, lon: f64, lat: f64) -> Option<[f64; 2]> {
        (**self).forward(lon, lat)
    }
}

/// A [`Forward`] map made from a closure; see [`from_fn`].
#[derive(Clone, Copy, Debug)]
pub struct FromFn<F>(pub F);

impl<F: FnMut(f64, f64) -> Option<[f64; 2]>> Forward for FromFn<F> {
    fn forward(&mut self, lon: f64, lat: f64) -> Option<[f64; 2]> {
        (self.0)(lon, lat)
    }
}

/// Wrap a closure `(x, y) -> Option<[x, y]>` as a [`Forward`] map.
pub fn from_fn<F: FnMut(f64, f64) -> Option<[f64; 2]>>(f: F) -> FromFn<F> {
    FromFn(f)
}

/// A forward map that moves many points per call: the shape of a
/// transformer that works on coordinate arrays (PROJ's `proj_trans_generic`,
/// rwarp's `Transformer::transform`).
///
/// The engine calls it once per level of refinement with every midpoint
/// that level needs, so the number of calls is the depth reached, not the
/// number of points.
///
/// Every [`Forward`] map is a `ForwardBatch` (point by point), so the
/// drivers take either.
pub trait ForwardBatch {
    /// Move `x[i], y[i]` in place and set `ok[i]` to whether the point
    /// could be moved (non-finite output is treated as failure by the
    /// engine whatever `ok` says).
    ///
    /// `near[i]` is the target plane position of the start of the edge
    /// the point refines: a neighbour of the answer, for maps whose
    /// inverse is a search (a curvilinear grid's cell location) to start
    /// from. Maps that do not search ignore it. All slices have the same
    /// length.
    fn forward_batch(&mut self, x: &mut [f64], y: &mut [f64], near: &[[f64; 2]], ok: &mut [bool]);
}

impl<T: Forward> ForwardBatch for T {
    fn forward_batch(&mut self, x: &mut [f64], y: &mut [f64], _near: &[[f64; 2]], ok: &mut [bool]) {
        for i in 0..x.len() {
            match self.forward(x[i], y[i]) {
                Some([px, py]) => {
                    x[i] = px;
                    y[i] = py;
                    ok[i] = true;
                }
                None => {
                    x[i] = f64::NAN;
                    y[i] = f64::NAN;
                    ok[i] = false;
                }
            }
        }
    }
}
