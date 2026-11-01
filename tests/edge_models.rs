//! Edge models and the level-order (batched) engine.

use bigcurve::{
    densify_mesh_with, densify_path, densify_path_with, from_fn, Forward, ForwardBatch,
    GreatCircle, Options, Planar, Resampler, Vertex,
};

/// Spherical Lambert azimuthal equal area centred on (lon0, lat0), in
/// metres; `None` at the antipode.
fn laea(lon0: f64, lat0: f64) -> impl Forward + Copy {
    let r = 6371008.8;
    let (sp0, cp0) = lat0.to_radians().sin_cos();
    from_fn(move |lon: f64, lat: f64| {
        let (sp, cp) = lat.to_radians().sin_cos();
        let (sl, cl) = (lon - lon0).to_radians().sin_cos();
        let d = 1.0 + sp0 * sp + cp0 * cp * cl;
        if d <= 1e-12 {
            return None;
        }
        let k = (2.0 / d).sqrt();
        Some([r * k * cp * sl, r * k * (cp0 * sp - sp0 * cp * cl)])
    })
}

/// A polar sector grid's index space: column from the radius, row from
/// the angle (degrees). Planar edges in x, y become curves here.
fn polar_index(x: f64, y: f64) -> Option<[f64; 2]> {
    Some([x.hypot(y) / 10.0, y.atan2(x).to_degrees() / 2.0])
}

fn bits(v: &[Vertex]) -> Vec<[u64; 4]> {
    v.iter()
        .map(|p| {
            [
                p.lon.to_bits(),
                p.lat.to_bits(),
                p.x.to_bits(),
                p.y.to_bits(),
            ]
        })
        .collect()
}

/// A tiny deterministic generator, so the test needs no dependencies.
struct Lcg(u64);
impl Lcg {
    fn next(&mut self) -> f64 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        (self.0 >> 11) as f64 / (1u64 << 53) as f64
    }
}

#[test]
fn level_order_matches_depth_first_bitwise() {
    let mut g = Lcg(42);
    let mut edges = Vec::new();
    for _ in 0..300 {
        let a = [g.next() * 360.0 - 180.0, g.next() * 170.0 - 85.0];
        let b = [g.next() * 360.0 - 180.0, g.next() * 170.0 - 85.0];
        edges.push((a, b));
    }
    let opts = Options::new(1e4);
    let mut rs = Resampler::new(laea(147.0, -42.0), &opts);
    let pairs: Vec<(Vertex, Vertex)> = edges
        .iter()
        .map(|(a, b)| (rs.project(a[0], a[1]), rs.project(b[0], b[1])))
        .collect();

    let (mut all, mut ends) = (Vec::new(), Vec::new());
    rs.edges(&pairs, &mut all, &mut ends);
    assert_eq!(ends.len(), pairs.len());
    let mut from = 0;
    let mut total = 0;
    for (k, (a, b)) in pairs.iter().enumerate() {
        let mut one = Vec::new();
        rs.edge(a, b, &mut one);
        assert_eq!(bits(&all[from..ends[k]]), bits(&one), "edge {k}");
        total += one.len();
        from = ends[k];
    }
    assert!(total > 1000, "the case should refine: {total}");
}

/// Counts calls, implements only the batched trait (so this also checks a
/// downstream type can implement it beside the blanket impl), and checks
/// the `near` hints.
struct Counting {
    calls: usize,
    points: usize,
}

impl ForwardBatch for Counting {
    fn forward_batch(&mut self, x: &mut [f64], y: &mut [f64], near: &[[f64; 2]], ok: &mut [bool]) {
        self.calls += 1;
        self.points += x.len();
        for i in 0..x.len() {
            let p = polar_index(x[i], y[i]).unwrap();
            // the hint is the start of the refined edge: a neighbour, within
            // the edge's own span of index space
            if near[i][0].is_finite() {
                assert!((p[0] - near[i][0]).abs() < 20.0 && (p[1] - near[i][1]).abs() < 50.0);
            }
            x[i] = p[0];
            y[i] = p[1];
            ok[i] = true;
        }
    }
}

#[test]
fn one_forward_call_per_level() {
    // clear of the angle's branch cut along the negative x axis
    let ring = [[100.0, 0.0], [0.0, 100.0], [-100.0, 10.0]];
    let opts = Options::new(1e-3);
    let mut f = Counting {
        calls: 0,
        points: 0,
    };
    let d = densify_path_with(&mut f, Planar::new(), &ring, &opts);
    // one call for the input vertices, then one per level reached
    assert!(
        f.calls >= 2 && f.calls <= 1 + opts.max_depth as usize,
        "{}",
        f.calls
    );
    assert!(f.points >= d.len());
    // and the same vertices as the point by point map
    let e = densify_path_with(from_fn(polar_index), Planar::new(), &ring, &opts);
    assert_eq!(d, e);
}

impl ForwardBatch for &mut Counting {
    fn forward_batch(&mut self, x: &mut [f64], y: &mut [f64], near: &[[f64; 2]], ok: &mut [bool]) {
        (**self).forward_batch(x, y, near, ok)
    }
}

#[test]
fn planar_under_an_affine_map_adds_nothing() {
    let affine =
        from_fn(|x: f64, y: f64| Some([2.0 * x + 0.5 * y + 3.0, -0.25 * x + 1.5 * y - 7.0]));
    let ring = [
        [0.0, 0.0],
        [100.0, 0.0],
        [100.0, 50.0],
        [0.0, 50.0],
        [0.0, 0.0],
    ];
    let d = densify_path_with(affine, Planar::new(), &ring, &Options::new(1e-9));
    assert_eq!(d, ring.to_vec());
}

#[test]
fn planar_converges_at_the_square_root_rate() {
    // a quarter circle of radius 100 into polar index space, and back: the
    // vertex count grows about sqrt(10) times per tenfold tolerance
    let line = [[100.0, 0.0], [0.0, 100.0]];
    let n: Vec<usize> = [1e-2, 1e-3, 1e-4, 1e-5]
        .iter()
        .map(|&t| {
            densify_path_with(
                from_fn(polar_index),
                Planar::new(),
                &line,
                &Options::new(t).max_depth(30),
            )
            .len()
        })
        .collect();
    // the symmetric curve refines whole levels at a time (counts are
    // 2^k + 1), so judge the rate over the three decades
    let rate = (n[3] as f64 / n[0] as f64).powf(1.0 / 3.0);
    assert!((2.6..3.8).contains(&rate), "{n:?}: {rate}");
}

#[test]
fn planar_long_forces_splits() {
    // under the identity nothing bends, but edges longer than 10 split
    let id = from_fn(|x, y| Some([x, y]));
    let d = densify_path_with(
        id,
        Planar::new().with_long(10.0),
        &[[0.0, 0.0], [100.0, 0.0]],
        &Options::new(1.0),
    );
    assert_eq!(d.len(), 17);
    assert!(d.windows(2).all(|w| (w[1][0] - w[0][0]) < 10.0));
}

#[test]
fn planar_shared_edges_refine_once() {
    // two squares sharing the edge 1 -> 2, written in opposite directions
    let v = [
        [100.0, 1.0],
        [150.0, 1.0],
        [150.0, 60.0],
        [100.0, 60.0],
        [200.0, 1.0],
        [200.0, 60.0],
    ];
    let segs = [
        [0, 1],
        [1, 2],
        [2, 3],
        [3, 0],
        [1, 4],
        [4, 5],
        [5, 2],
        [2, 1],
    ];
    let m = densify_mesh_with(
        from_fn(polar_index),
        Planar::new(),
        &v,
        &segs,
        &Options::new(1e-3),
    )
    .unwrap();
    let chain = |k: usize| -> Vec<usize> {
        let mut c: Vec<usize> = m
            .segments
            .iter()
            .zip(&m.parent)
            .filter(|(_, &p)| p == k)
            .map(|(s, _)| s[0])
            .collect();
        c.push(
            m.segments
                .iter()
                .zip(&m.parent)
                .rfind(|(_, &p)| p == k)
                .unwrap()
                .0[1],
        );
        c
    };
    let mut back = chain(7);
    back.reverse();
    assert!(chain(1).len() > 2);
    assert_eq!(chain(1), back);
}

#[test]
fn the_great_circle_default_is_unchanged() {
    let line = [[140.0, -60.0], [-150.0, -50.0], [-30.0, -70.0]];
    let a = densify_path(laea(147.0, -42.0), &line, &Options::new(2e4));
    let b = densify_path_with(laea(147.0, -42.0), GreatCircle, &line, &Options::new(2e4));
    assert_eq!(a, b);
    assert!(a.len() > line.len());
}
