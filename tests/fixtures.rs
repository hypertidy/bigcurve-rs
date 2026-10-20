//! Bitwise agreement with the C++ engine of the R package.
//!
//! Each fixture in tests/fixtures was written by oracle/oracle.cpp: the
//! input, the settings, every forward call the C++ engine made through
//! PROJ (point in, point out), and the output. The test replays those
//! calls as the forward map. A point the C++ engine never asked for is a
//! failure, and so is any difference in the output bits.
//!
//! The fixtures were recorded on Linux with glibc, and the engine's
//! trigonometry comes from the platform's maths library, which is not
//! correctly rounded. So the comparison is bitwise on Linux with glibc,
//! and elsewhere matches each requested point to the nearest recorded one
//! and compares vertex counts exactly and coordinates to 1e-9 degrees.

use std::collections::HashMap;
use std::fs;
use std::path::Path;

use bigcurve::{default_tolerance, densify_mesh, densify_path, Forward, Options};

#[derive(Default)]
struct Fixture {
    kind: String,
    pixels: Option<f64>,
    tolerance: f64,
    max_depth: u32,
    input: Vec<[f64; 2]>,
    segments: Vec<[usize; 2]>,
    calls: HashMap<(u64, u64), Option<[f64; 2]>>,
    output: Vec<[f64; 2]>,
    out_segments: Vec<([usize; 2], usize)>,
}

fn hex(s: &str) -> f64 {
    f64::from_bits(u64::from_str_radix(s, 16).expect("hex double"))
}

fn read(path: &Path) -> Fixture {
    let mut f = Fixture::default();
    for line in fs::read_to_string(path).unwrap().lines() {
        let w: Vec<&str> = line.split_whitespace().collect();
        match w.first().copied() {
            None | Some("#") => {}
            Some("kind") => f.kind = w[1].to_string(),
            Some("pixels") => f.pixels = Some(hex(w[1])),
            Some("tolerance") => f.tolerance = hex(w[1]),
            Some("max_depth") => f.max_depth = w[1].parse().unwrap(),
            Some("in") => f.input.push([hex(w[1]), hex(w[2])]),
            Some("seg") => f
                .segments
                .push([w[1].parse().unwrap(), w[2].parse().unwrap()]),
            Some("fwd") => {
                let key = (
                    u64::from_str_radix(w[1], 16).unwrap(),
                    u64::from_str_radix(w[2], 16).unwrap(),
                );
                let val = if w[3] == "ok" {
                    Some([hex(w[4]), hex(w[5])])
                } else {
                    None
                };
                f.calls.insert(key, val);
            }
            Some("out") => f.output.push([hex(w[1]), hex(w[2])]),
            Some("oseg") => f.out_segments.push((
                [w[1].parse().unwrap(), w[2].parse().unwrap()],
                w[3].parse().unwrap(),
            )),
            Some(other) => panic!("{}: unknown record {other}", path.display()),
        }
    }
    f
}

/// Bitwise where the fixtures were recorded.
const STRICT: bool = cfg!(all(target_os = "linux", target_env = "gnu"));

/// Agreement in degrees when not bitwise.
const CLOSE: f64 = 1e-9;

/// Replays the recorded forward calls.
struct Replay<'a> {
    name: &'a str,
    calls: &'a HashMap<(u64, u64), Option<[f64; 2]>>,
}

impl Forward for Replay<'_> {
    fn forward(&mut self, lon: f64, lat: f64) -> Option<[f64; 2]> {
        if let Some(v) = self.calls.get(&(lon.to_bits(), lat.to_bits())) {
            return *v;
        }
        if !STRICT {
            let near = self.calls.iter().find(|((x, y), _)| {
                (f64::from_bits(*x) - lon).abs() < CLOSE && (f64::from_bits(*y) - lat).abs() < CLOSE
            });
            if let Some((_, v)) = near {
                return *v;
            }
        }
        panic!(
            "{}: the C++ engine never moved ({lon:?}, {lat:?})",
            self.name
        )
    }
}

fn same(name: &str, what: &str, a: &[[f64; 2]], b: &[[f64; 2]]) {
    assert_eq!(a.len(), b.len(), "{name}: {what} count");
    for (i, (p, q)) in a.iter().zip(b).enumerate() {
        let ok = if STRICT {
            p[0].to_bits() == q[0].to_bits() && p[1].to_bits() == q[1].to_bits()
        } else {
            (p[0] - q[0]).abs() < CLOSE && (p[1] - q[1]).abs() < CLOSE
        };
        assert!(ok, "{name}: {what} {i}: {p:?} != {q:?}");
    }
}

fn check(path: &Path) {
    let name = path.file_stem().unwrap().to_str().unwrap();
    let f = read(path);
    let replay = || Replay {
        name,
        calls: &f.calls,
    };

    if let Some(pixels) = f.pixels {
        let tol = default_tolerance(replay(), &f.input, pixels).expect("a default tolerance");
        if STRICT {
            assert_eq!(
                tol.to_bits(),
                f.tolerance.to_bits(),
                "{name}: default tolerance"
            );
        } else {
            assert!(
                (tol / f.tolerance - 1.0).abs() < 1e-12,
                "{name}: default tolerance"
            );
        }
    }
    let opts = Options::new(f.tolerance).max_depth(f.max_depth);

    match f.kind.as_str() {
        "path" => {
            let out = densify_path(replay(), &f.input, &opts);
            same(name, "vertex", &out, &f.output);
        }
        "mesh" => {
            let out = densify_mesh(replay(), &f.input, &f.segments, &opts).unwrap();
            same(name, "vertex", &out.vertices, &f.output);
            let segs: Vec<([usize; 2], usize)> = out
                .segments
                .iter()
                .copied()
                .zip(out.parent.iter().copied())
                .collect();
            assert_eq!(segs, f.out_segments, "{name}: segments and parents");
        }
        k => panic!("{name}: unknown kind {k}"),
    }
}

#[test]
fn matches_the_cpp_engine() {
    let dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures");
    let mut files: Vec<_> = fs::read_dir(&dir)
        .unwrap()
        .map(|e| e.unwrap().path())
        .filter(|p| p.extension().is_some_and(|e| e == "txt"))
        .collect();
    files.sort();
    assert!(
        files.len() >= 30,
        "expected the oracle's fixtures in {}",
        dir.display()
    );
    for p in &files {
        check(p);
    }
}
