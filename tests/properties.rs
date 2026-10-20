//! Properties of the engine that hold for any forward map, checked with
//! an analytic spherical Lambert azimuthal equal area projection.

use bigcurve::{densify_mesh, densify_path, from_fn, Forward, Options};

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

fn bits(v: &[[f64; 2]]) -> Vec<[u64; 2]> {
    v.iter().map(|p| [p[0].to_bits(), p[1].to_bits()]).collect()
}

#[test]
fn an_edge_refines_identically_in_either_direction() {
    let opts = Options::new(2e4);
    let f = densify_path(laea(50.0, 30.0), &[[0.0, 0.0], [100.0, 50.0]], &opts);
    let mut r = densify_path(laea(50.0, 30.0), &[[100.0, 50.0], [0.0, 0.0]], &opts);
    assert!(f.len() > 2);
    r.reverse();
    assert_eq!(bits(&f), bits(&r));
}

#[test]
fn max_depth_caps_refinement() {
    let line = [[140.0, -89.0], [-164.0, 80.0]];
    for depth in 0..5u32 {
        let d = densify_path(
            laea(147.0, -42.0),
            &line,
            &Options::new(1.0).max_depth(depth),
        );
        assert!(
            d.len() < 2 + (1usize << depth),
            "depth {depth}: {}",
            d.len()
        );
    }
}

#[test]
fn unmovable_points_leave_edges_alone() {
    // the antipode of the projection centre cannot be moved
    let line = [[-130.0, -30.0], [0.0, 0.0], [10.0, 10.0]];
    let d = densify_path(laea(50.0, 30.0), &line, &Options::new(1e3));
    assert_eq!(d[0], line[0]);
    assert!(d.len() > 3, "the second edge still refines");
    let none = from_fn(|_, _| None);
    assert_eq!(densify_path(none, &line, &Options::new(1e3)), line.to_vec());
}

/// Merging shared vertices before densifying (each shared edge refined
/// once) gives the same new vertices as densifying every cell on its own:
/// the twins agree to the last bit.
#[test]
fn shared_edges_get_identical_vertices() {
    let (nx, ny) = (6usize, 3usize);
    let corner = |i: usize, j: usize| [-180.0 + 60.0 * i as f64, -75.0 + 15.0 * j as f64];
    let opts = Options::new(5e3);
    let prj = laea(0.0, -90.0);

    // merged: one vertex per grid corner
    let merged_v: Vec<[f64; 2]> = (0..=ny)
        .flat_map(|j| (0..=nx).map(move |i| corner(i, j)))
        .collect();
    let id = |i: usize, j: usize| j * (nx + 1) + i;
    let mut merged_s = Vec::new();
    // unmerged: four vertices per cell
    let mut cells_v = Vec::new();
    let mut cells_s = Vec::new();
    for j in 0..ny {
        for i in 0..nx {
            let ring = [(i, j), (i + 1, j), (i + 1, j + 1), (i, j + 1)];
            let base = cells_v.len();
            for k in 0..4 {
                let (a, b) = (ring[k], ring[(k + 1) % 4]);
                merged_s.push([id(a.0, a.1), id(b.0, b.1)]);
                cells_v.push(corner(a.0, a.1));
                cells_s.push([base + k, base + (k + 1) % 4]);
            }
        }
    }
    let merged = densify_mesh(prj, &merged_v, &merged_s, &opts).unwrap();
    let cells = densify_mesh(prj, &cells_v, &cells_s, &opts).unwrap();

    let mut a = bits(&merged.vertices[merged_v.len()..]);
    let mut b = bits(&cells.vertices[cells_v.len()..]);
    assert!(!a.is_empty());
    // every interior edge appears twice among the cells
    a.sort();
    b.sort();
    b.dedup();
    assert_eq!(a, b);
}
