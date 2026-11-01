//! Canon tables: the cases this crate pins and the C++ header follows.
//!
//! This crate is the canonical bigcurve engine. `tests/canon/canon-tables.txt`
//! holds a fixed set of cases (inputs, settings, every forward call the
//! engine made, and the output), written by this file and copied into the
//! R package (hypertidy/bigcurve, `canon/`), whose C++ header must
//! reproduce every great circle case. The tables are the contract between
//! the two engines: change the Rust, regenerate, copy, port the C++ until
//! its check is green (or mark the case Rust-only by its edge or place).
//!
//! The forward maps here are analytic, on the sphere, so no PROJ is
//! needed to write the tables; the forward calls are recorded, so neither
//! side needs the maps to check them. The inputs are built with exact
//! arithmetic (no trigonometry), so they are the same everywhere.
//!
//! Regenerate (on Linux with glibc, where the tables are bitwise):
//!
//!     CANON_WRITE=1 cargo test --all-features --test canon
//!
//! Without `CANON_WRITE` the test replays the recorded calls and requires
//! the same output: bitwise on Linux with glibc, elsewhere to 1e-9 in
//! the source coordinates (the engine's own trigonometry comes from the
//! platform's maths library).
//!
//! Format, one record per line, doubles as IEEE 754 bit patterns in hex:
//!
//!     case <name>
//!     kind path|mesh
//!     edge greatcircle|planar
//!     place sphere|ellipsoid
//!     tolerance <hex>
//!     max_depth <n>
//!     in <x> <y>                 input vertices
//!     seg <a> <b>                input segments (mesh)
//!     fwd <x> <y> ok|na <X> <Y>  every forward call, sorted by input bits
//!     out <x> <y>                output vertices
//!     oseg <a> <b> <parent>      output segments (mesh)
//!     end
#![cfg(feature = "geodesic")]

use std::collections::BTreeMap;
use std::fmt::Write as _;
use std::fs;
use std::path::PathBuf;

use bigcurve::{
    densify_mesh_ellipsoid, densify_mesh_with, densify_path_ellipsoid, densify_path_with, Forward,
    GreatCircle, Mesh, Options, Planar,
};

const R: f64 = 6378137.0;

/// Bitwise where the tables are written.
const STRICT: bool = cfg!(all(target_os = "linux", target_env = "gnu"));
const CLOSE: f64 = 1e-9;

#[derive(Clone, Copy, Debug, PartialEq)]
enum Map {
    /// spherical Mercator, longitude wrapped to [-180, 180) as PROJ does
    Merc,
    /// orthographic on the sphere; the far side cannot be moved
    Ortho(f64, f64),
    /// Lambert azimuthal equal area on the sphere
    Laea(f64, f64),
    /// sinusoidal on the sphere, wrapped
    Sinu,
    /// planar input: radius and angle (degrees) over a cell of 10 units,
    /// the index space of a polar sector grid
    Sector,
}

impl Map {
    fn name(self) -> String {
        match self {
            Map::Merc => "merc".into(),
            Map::Ortho(a, b) => format!("ortho{a}_{b}"),
            Map::Laea(a, b) => format!("laea{a}_{b}"),
            Map::Sinu => "sinu".into(),
            Map::Sector => "sector".into(),
        }
    }

    fn apply(self, x: f64, y: f64) -> Option<[f64; 2]> {
        let wrap = |l: f64| (l + 180.0).rem_euclid(360.0) - 180.0;
        match self {
            Map::Merc => {
                if y.abs() > 85.0 {
                    return None;
                }
                let t = (std::f64::consts::FRAC_PI_4 + y.to_radians() / 2.0)
                    .tan()
                    .ln();
                Some([R * wrap(x).to_radians(), R * t])
            }
            Map::Sinu => {
                let p = y.to_radians();
                Some([R * wrap(x).to_radians() * p.cos(), R * p])
            }
            Map::Ortho(l0, p0) | Map::Laea(l0, p0) => {
                let (l, p) = ((x - l0).to_radians(), y.to_radians());
                let p0 = p0.to_radians();
                let cosc = p0.sin() * p.sin() + p0.cos() * p.cos() * l.cos();
                let (u, v) = (
                    p.cos() * l.sin(),
                    p0.cos() * p.sin() - p0.sin() * p.cos() * l.cos(),
                );
                if let Map::Ortho(..) = self {
                    (cosc >= 0.0).then_some([R * u, R * v])
                } else {
                    if cosc <= -1.0 + 1e-12 {
                        return None;
                    }
                    let k = (2.0 / (1.0 + cosc)).sqrt();
                    Some([R * k * u, R * k * v])
                }
            }
            Map::Sector => Some([x.hypot(y) / 10.0, y.atan2(x).to_degrees()]),
        }
    }
}

struct Rec {
    map: Map,
    calls: BTreeMap<(u64, u64), Option<[f64; 2]>>,
}

impl Forward for Rec {
    fn forward(&mut self, x: f64, y: f64) -> Option<[f64; 2]> {
        let v = self
            .map
            .apply(x, y)
            .filter(|v| v[0].is_finite() && v[1].is_finite());
        self.calls.insert((x.to_bits(), y.to_bits()), v);
        v
    }
}

#[derive(Clone, Debug, Default)]
struct Case {
    name: String,
    mesh: bool,
    planar: bool,
    ellipsoid: bool,
    tolerance: f64,
    max_depth: u32,
    input: Vec<[f64; 2]>,
    segments: Vec<[usize; 2]>,
    calls: BTreeMap<(u64, u64), Option<[f64; 2]>>,
    output: Vec<[f64; 2]>,
    out_segments: Vec<[usize; 3]>,
}

fn run<F: Forward>(c: &Case, f: F) -> (Vec<[f64; 2]>, Vec<[usize; 3]>) {
    let o = Options::new(c.tolerance).max_depth(c.max_depth);
    let flat = |m: Mesh| {
        let segs = m
            .segments
            .iter()
            .zip(&m.parent)
            .map(|(s, &p)| [s[0], s[1], p])
            .collect();
        (m.vertices, segs)
    };
    match (c.mesh, c.planar, c.ellipsoid) {
        (false, false, false) => (densify_path_with(f, GreatCircle, &c.input, &o), vec![]),
        (false, true, _) => (densify_path_with(f, Planar::new(), &c.input, &o), vec![]),
        (false, false, true) => (densify_path_ellipsoid(f, &c.input, &o), vec![]),
        (true, false, false) => {
            flat(densify_mesh_with(f, GreatCircle, &c.input, &c.segments, &o).unwrap())
        }
        (true, true, _) => {
            flat(densify_mesh_with(f, Planar::new(), &c.input, &c.segments, &o).unwrap())
        }
        (true, false, true) => flat(densify_mesh_ellipsoid(f, &c.input, &c.segments, &o).unwrap()),
    }
}

/// A path through exact points: whole and half degrees only.
fn path(name: &str, pts: &[[f64; 2]]) -> Case {
    Case {
        name: name.into(),
        input: pts.to_vec(),
        max_depth: 16,
        ..Default::default()
    }
}

fn inputs() -> Vec<Case> {
    let mut v = vec![
        path("hobart_london", &[[147.5, -43.0], [-0.0, 51.5]]),
        path(
            "routes",
            &[
                [147.5, -43.0],
                [-0.5, 51.5],
                [-70.5, -33.5],
                [28.0, -26.0],
                [-74.0, 40.5],
                [37.5, 56.0],
                [116.0, -32.0],
            ],
        ),
        path("crossing", &[[140.0, -89.0], [-164.0, 80.0]]),
        path("seam_written", &[[170.0, 20.0], [-170.0, 20.0]]),
        path("seam_short", &[[179.5, 0.0], [-179.5, 0.0]]),
        path("along_180", &[[180.0, -60.0], [180.0, 60.0]]),
        path("along_m180", &[[-180.0, -60.0], [-180.0, 60.0]]),
        path("equator_e", &[[0.0, 0.0], [180.0, 0.0]]),
        path("equator_w", &[[0.0, 0.0], [-180.0, 0.0]]),
        path("parallel", &[[0.0, 60.0], [100.0, 60.0]]),
        path("polar_hop", &[[-45.0, -70.0], [135.0, -70.0]]),
        path("tiny", &[[10.0, 10.0], [10.5, 10.0]]),
        path("repeat_vertex", &[[10.0, 10.0], [10.0, 10.0], [60.0, 30.0]]),
    ];
    // a scatter of edges from exact integer arithmetic
    let mut pts = Vec::new();
    for i in 0..12u32 {
        let lon = -180.0 + f64::from((i * 37) % 360) + 0.5 * f64::from(i % 3);
        let lat = -80.0 + f64::from((i * 23) % 161);
        pts.push([lon, lat]);
    }
    v.push(path("scatter", &pts));
    // a merged southern graticule: shared edges both ways
    let (nx, ny) = (8usize, 3usize);
    let mut g = Case {
        name: "graticule".into(),
        mesh: true,
        max_depth: 16,
        ..Default::default()
    };
    for j in 0..=ny {
        for i in 0..=nx {
            g.input
                .push([-180.0 + 45.0 * i as f64, -80.0 + 20.0 * j as f64]);
        }
    }
    let id = |i: usize, j: usize| j * (nx + 1) + i;
    for j in 0..ny {
        for i in 0..nx {
            let r = [
                id(i, j),
                id(i + 1, j),
                id(i + 1, j + 1),
                id(i, j + 1),
                id(i, j),
            ];
            for k in 0..4 {
                g.segments.push([r[k], r[k + 1]]);
            }
        }
    }
    v.push(g);
    v.push(Case {
        name: "neighbours".into(),
        mesh: true,
        max_depth: 16,
        input: vec![[0.0, 0.0], [100.0, 50.0], [0.0, 60.0], [100.0, 0.0]],
        segments: vec![[0, 1], [1, 2], [2, 0], [0, 3], [3, 1], [1, 0]],
        ..Default::default()
    });
    v
}

fn planar_inputs() -> Vec<Case> {
    vec![
        path("quarter", &[[100.0, 0.0], [0.0, 100.0]]),
        path("chord", &[[100.0, 10.0], [-50.0, 80.0], [-90.0, -40.0]]),
        path("radial", &[[10.0, 10.0], [200.0, 200.0]]),
        Case {
            name: "square".into(),
            mesh: true,
            max_depth: 16,
            input: vec![[20.0, 20.0], [120.0, 20.0], [120.0, 120.0], [20.0, 120.0]],
            segments: vec![[0, 1], [1, 2], [2, 3], [3, 0], [0, 2], [2, 0]],
            ..Default::default()
        },
    ]
}

fn cases() -> Vec<(Case, Map)> {
    let maps = [
        Map::Merc,
        Map::Ortho(100.0, 0.0),
        Map::Laea(147.0, -42.0),
        Map::Laea(0.0, -90.0),
        Map::Sinu,
    ];
    let mut out = Vec::new();
    for base in inputs() {
        for &m in &maps {
            for (tn, tol) in [("1e5", 1e5), ("1e4", 1e4)] {
                let mut c = base.clone();
                c.name = format!("{}_{}_{tn}", base.name, m.name());
                c.tolerance = tol;
                out.push((c, m));
            }
        }
    }
    // depth limits, and the ellipsoid placement (Rust only for now)
    for base in inputs().into_iter().take(4) {
        let mut c = base.clone();
        c.name = format!("{}_laea147_-42_depth3", base.name);
        c.tolerance = 1.0;
        c.max_depth = 3;
        out.push((c, Map::Laea(147.0, -42.0)));
        let mut c = base.clone();
        c.name = format!("{}_sinu_ellipsoid", base.name);
        c.tolerance = 1e4;
        c.ellipsoid = true;
        out.push((c, Map::Sinu));
    }
    for base in planar_inputs() {
        for (tn, tol) in [("1e-1", 1e-1), ("1e-3", 1e-3)] {
            let mut c = base.clone();
            c.name = format!("{}_sector_{tn}", base.name);
            c.tolerance = tol;
            c.planar = true;
            out.push((c, Map::Sector));
        }
    }
    out
}

fn h(x: f64) -> String {
    format!("{:016x}", x.to_bits())
}

fn write(path: &PathBuf) {
    let mut s = String::new();
    s.push_str("# canonical tables from bigcurve-rs tests/canon.rs; do not edit\n");
    s.push_str("# doubles are IEEE 754 bit patterns in hex\n");
    let mut n = 0;
    for (mut c, map) in cases() {
        let mut rec = Rec {
            map,
            calls: BTreeMap::new(),
        };
        let (o, os) = run(&c, &mut rec);
        c.calls = rec.calls;
        c.output = o;
        c.out_segments = os;
        let _ = writeln!(s, "case {}", c.name);
        let _ = writeln!(s, "kind {}", if c.mesh { "mesh" } else { "path" });
        let _ = writeln!(
            s,
            "edge {}",
            if c.planar { "planar" } else { "greatcircle" }
        );
        let _ = writeln!(
            s,
            "place {}",
            if c.ellipsoid { "ellipsoid" } else { "sphere" }
        );
        let _ = writeln!(s, "tolerance {}", h(c.tolerance));
        let _ = writeln!(s, "max_depth {}", c.max_depth);
        for p in &c.input {
            let _ = writeln!(s, "in {} {}", h(p[0]), h(p[1]));
        }
        for g in &c.segments {
            let _ = writeln!(s, "seg {} {}", g[0], g[1]);
        }
        for (&(x, y), v) in &c.calls {
            match v {
                Some(p) => {
                    let _ = writeln!(s, "fwd {x:016x} {y:016x} ok {} {}", h(p[0]), h(p[1]));
                }
                None => {
                    let _ = writeln!(
                        s,
                        "fwd {x:016x} {y:016x} na {} {}",
                        h(f64::NAN),
                        h(f64::NAN)
                    );
                }
            }
        }
        for p in &c.output {
            let _ = writeln!(s, "out {} {}", h(p[0]), h(p[1]));
        }
        for g in &c.out_segments {
            let _ = writeln!(s, "oseg {} {} {}", g[0], g[1], g[2]);
        }
        s.push_str("end\n");
        n += 1;
    }
    fs::write(path, s).unwrap();
    eprintln!("wrote {n} cases to {}", path.display());
}

fn hex(s: &str) -> f64 {
    f64::from_bits(u64::from_str_radix(s, 16).expect("hex double"))
}

fn read(path: &PathBuf) -> Vec<Case> {
    let mut out = Vec::new();
    let mut c = Case::default();
    for line in fs::read_to_string(path).unwrap().lines() {
        let w: Vec<&str> = line.split_whitespace().collect();
        match w.first().copied() {
            None | Some("#") => {}
            Some("case") => {
                c = Case {
                    name: w[1].into(),
                    ..Default::default()
                }
            }
            Some("kind") => c.mesh = w[1] == "mesh",
            Some("edge") => c.planar = w[1] == "planar",
            Some("place") => c.ellipsoid = w[1] == "ellipsoid",
            Some("tolerance") => c.tolerance = hex(w[1]),
            Some("max_depth") => c.max_depth = w[1].parse().unwrap(),
            Some("in") => c.input.push([hex(w[1]), hex(w[2])]),
            Some("seg") => c
                .segments
                .push([w[1].parse().unwrap(), w[2].parse().unwrap()]),
            Some("fwd") => {
                let k = (
                    u64::from_str_radix(w[1], 16).unwrap(),
                    u64::from_str_radix(w[2], 16).unwrap(),
                );
                let v = (w[3] == "ok").then(|| [hex(w[4]), hex(w[5])]);
                c.calls.insert(k, v);
            }
            Some("out") => c.output.push([hex(w[1]), hex(w[2])]),
            Some("oseg") => c.out_segments.push([
                w[1].parse().unwrap(),
                w[2].parse().unwrap(),
                w[3].parse().unwrap(),
            ]),
            Some("end") => out.push(std::mem::take(&mut c)),
            Some(other) => panic!("unknown record {other}"),
        }
    }
    out
}

struct Replay<'a>(&'a Case);

impl Forward for Replay<'_> {
    fn forward(&mut self, x: f64, y: f64) -> Option<[f64; 2]> {
        if let Some(v) = self.0.calls.get(&(x.to_bits(), y.to_bits())) {
            return *v;
        }
        if !STRICT {
            let near = self.0.calls.iter().find(|((a, b), _)| {
                (f64::from_bits(*a) - x).abs() < CLOSE && (f64::from_bits(*b) - y).abs() < CLOSE
            });
            if let Some((_, v)) = near {
                return *v;
            }
        }
        panic!("{}: no recorded call for ({x:?}, {y:?})", self.0.name)
    }
}

#[test]
fn canon_tables() {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/canon/canon-tables.txt");
    if std::env::var_os("CANON_WRITE").is_some() {
        write(&path);
    }
    let cases = read(&path);
    assert_eq!(
        cases.len(),
        self::cases().len(),
        "case count; regenerate with CANON_WRITE=1"
    );
    for c in &cases {
        let (o, os) = run(c, Replay(c));
        assert_eq!(o.len(), c.output.len(), "{}: vertex count", c.name);
        for (i, (p, q)) in o.iter().zip(&c.output).enumerate() {
            let ok = if STRICT {
                p[0].to_bits() == q[0].to_bits() && p[1].to_bits() == q[1].to_bits()
            } else {
                (p[0] - q[0]).abs() < CLOSE && (p[1] - q[1]).abs() < CLOSE
            };
            assert!(ok, "{}: vertex {i}: {p:?} != {q:?}", c.name);
        }
        assert_eq!(os, c.out_segments, "{}: segments", c.name);
    }
}
