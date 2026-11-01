//! `densify_*_ellipsoid` against the R package's `densify(ellipsoid =
//! TRUE)` rule.
//!
//! tests/fixtures/ellipsoid/*.txt were written by oracle/ellipsoid.py: the
//! C++ engine's sphere output from tests/fixtures, with each new vertex
//! moved as R/ellipsoid.R moves it (Python geographiclib for the
//! geodesics). Here the same input runs through the Rust drivers, with the
//! recorded forward calls as the map. Vertex counts must agree exactly;
//! coordinates to 1e-9 degrees, since the two sides use different ports of
//! GeographicLib and the R rule recovers fractions from arc angles where
//! the Rust engine knows them exactly.
#![cfg(feature = "geodesic")]

use std::collections::HashMap;
use std::fs;
use std::path::Path;

use bigcurve::{
    densify_mesh_ellipsoid, densify_path, densify_path_ellipsoid, from_fn, Forward, Options,
};

const CLOSE: f64 = 1e-9;

fn hex(s: &str) -> f64 {
    f64::from_bits(u64::from_str_radix(s, 16).expect("hex double"))
}

#[derive(Default)]
struct Fixture {
    kind: String,
    tolerance: f64,
    max_depth: u32,
    input: Vec<[f64; 2]>,
    segments: Vec<[usize; 2]>,
    calls: Vec<([f64; 2], Option<[f64; 2]>)>,
    exact: HashMap<(u64, u64), Option<[f64; 2]>>,
}

fn read(path: &Path) -> Fixture {
    let mut f = Fixture::default();
    for line in fs::read_to_string(path).unwrap().lines() {
        let w: Vec<&str> = line.split_whitespace().collect();
        match w.first().copied() {
            Some("kind") => f.kind = w[1].to_string(),
            Some("tolerance") => f.tolerance = hex(w[1]),
            Some("max_depth") => f.max_depth = w[1].parse().unwrap(),
            Some("in") => f.input.push([hex(w[1]), hex(w[2])]),
            Some("seg") => f
                .segments
                .push([w[1].parse().unwrap(), w[2].parse().unwrap()]),
            Some("fwd") => {
                let p = [hex(w[1]), hex(w[2])];
                let v = (w[3] == "ok").then(|| [hex(w[4]), hex(w[5])]);
                f.exact.insert((p[0].to_bits(), p[1].to_bits()), v);
                f.calls.push((p, v));
            }
            _ => {}
        }
    }
    f
}

fn read_moved(path: &Path) -> Vec<[f64; 2]> {
    fs::read_to_string(path)
        .unwrap()
        .lines()
        .filter_map(|l| {
            let w: Vec<&str> = l.split_whitespace().collect();
            (w.first() == Some(&"eout")).then(|| [hex(w[1]), hex(w[2])])
        })
        .collect()
}

struct Replay<'a>(&'a Fixture);

impl Forward for Replay<'_> {
    fn forward(&mut self, lon: f64, lat: f64) -> Option<[f64; 2]> {
        if let Some(v) = self.0.exact.get(&(lon.to_bits(), lat.to_bits())) {
            return *v;
        }
        // off glibc the sphere midpoints differ in the last bits
        self.0
            .calls
            .iter()
            .find(|(p, _)| (p[0] - lon).abs() < CLOSE && (p[1] - lat).abs() < CLOSE)
            .map(|(_, v)| *v)
            .unwrap_or_else(|| panic!("the C++ engine never moved ({lon:?}, {lat:?})"))
    }
}

#[test]
fn matches_the_r_ellipsoid_rule() {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures");
    let mut files: Vec<_> = fs::read_dir(root.join("ellipsoid"))
        .unwrap()
        .map(|e| e.unwrap().path())
        .filter(|p| p.extension().is_some_and(|e| e == "txt"))
        .collect();
    files.sort();
    assert!(files.len() >= 20, "expected oracle/ellipsoid.py's fixtures");
    let mut worst = 0.0f64;
    for p in &files {
        let name = p.file_stem().unwrap().to_str().unwrap();
        let f = read(&root.join(p.file_name().unwrap()));
        let want = read_moved(p);
        let opts = Options::new(f.tolerance).max_depth(f.max_depth);
        let got = match f.kind.as_str() {
            "path" => densify_path_ellipsoid(Replay(&f), &f.input, &opts),
            _ => {
                densify_mesh_ellipsoid(Replay(&f), &f.input, &f.segments, &opts)
                    .unwrap()
                    .vertices
            }
        };
        assert_eq!(got.len(), want.len(), "{name}: vertex count");
        for (i, (a, b)) in got.iter().zip(&want).enumerate() {
            let d = (a[0] - b[0]).abs().max((a[1] - b[1]).abs());
            worst = worst.max(d);
            assert!(d < CLOSE, "{name}: vertex {i}: {a:?} != {b:?}");
        }
        // input vertices never move
        if f.kind == "path" {
            assert_eq!(got[0], f.input[0], "{name}: first vertex");
            assert_eq!(
                got[got.len() - 1],
                f.input[f.input.len() - 1],
                "{name}: last vertex"
            );
        } else {
            assert_eq!(
                &got[..f.input.len()],
                &f.input[..],
                "{name}: input vertices"
            );
        }
    }
    eprintln!("ellipsoid fixtures: worst difference {worst:e} degrees");
}

/// An equal-earth-like analytic map (no PROJ), enough to need vertices.
fn eqc() -> impl Forward {
    from_fn(|lon: f64, lat: f64| {
        let r = 6378137.0;
        let p = lat.to_radians();
        Some([r * lon.to_radians() * p.cos(), r * p])
    })
}

#[test]
fn same_vertices_whichever_way_an_edge_is_written() {
    let s = [[147.3, -42.9], [-0.1, 51.5]];
    let opts = Options::new(1e4);
    let fwd = densify_path_ellipsoid(eqc(), &s, &opts);
    let mut rev = densify_path_ellipsoid(eqc(), &[s[1], s[0]], &opts);
    rev.reverse();
    assert!(fwd.len() > 10);
    assert_eq!(fwd, rev);
}

#[test]
fn same_count_as_the_sphere_and_off_the_great_circle() {
    let s = [[147.3, -42.9], [-0.1, 51.5]];
    let opts = Options::new(1e4);
    let sph = densify_path(eqc(), &s, &opts);
    let ell = densify_path_ellipsoid(eqc(), &s, &opts);
    assert_eq!(sph.len(), ell.len());
    let moved = sph
        .iter()
        .zip(&ell)
        .map(|(a, b)| (a[0] - b[0]).abs().max((a[1] - b[1]).abs()))
        .fold(0.0, f64::max);
    // the geodesic and the great circle part by tenths of a degree here
    assert!(moved > 0.05, "moved {moved}");
}

#[test]
fn keeps_the_longitude_frame_across_180() {
    let x = [[170.0, -40.0], [-170.0, -45.0]];
    let merc180 = from_fn(|lon: f64, lat: f64| {
        let l = (lon - 180.0 + 540.0).rem_euclid(360.0) - 180.0;
        let y = (std::f64::consts::FRAC_PI_4 + lat.to_radians() / 2.0)
            .tan()
            .ln();
        Some([6378137.0 * l.to_radians(), 6378137.0 * y])
    });
    let d = densify_path_ellipsoid(merc180, &x, &Options::new(1000.0));
    assert!(d.len() > 2);
    for v in &d[1..d.len() - 1] {
        assert!(v[0] > 170.0 && v[0] < 190.0, "{v:?}");
    }
}
