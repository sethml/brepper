//! Generates test models that cannot be produced via ccad, using OCCT
//! directly. The STL/STEP pairs are written to tests/occt/ and tests/nurbs/
//! and committed to the repo, so this tool only needs to be re-run when a
//! model changes.
//!
//! Usage: gen_test_models [<occt_dir>] [<nurbs_dir>]
//!        (defaults: tests/occt tests/nurbs)
//!
//! tests/occt/ models (analytic, exercise the existing pipeline):
//! - hollow_cube: 20 mm cube with a fully enclosed 10 mm cubic void at its
//!   center. Exercises multi-shell solid construction: the mesh has two
//!   closed shells and the STEP output must be ONE solid with an internal
//!   cavity (volume 7000 mm³), not two nested solids.
//! - full_torus: complete donut (R=10, r=3). Exercises the doubly-periodic
//!   seam handling in stage 3.4: a complete torus face needs both a U seam
//!   (around the major circle) and a V seam (around the tube).
//!
//! tests/nurbs/ models (freeform, targets for the DEVELOPMENT_PLAN §2.7
//! NURBS fallback — expected to FAIL --compare until that stage exists):
//!
//! Every reference STEP is built from EXACT geometry — explicit-pole
//! Geom_BSplineSurface lattices, conic/B-spline extrusions, and booleans of
//! those. No surface is produced by fitting sampled points: the recognizer
//! will use OCCT's fitting APIs, so a fitted reference would be circular
//! (shared biases cancel) and merely approximate.
//!
//! STLs come in two styles. The doubly-curved patch models are scripted
//! regular grids with vertices sampled exactly on the reference surface —
//! the cleanest possible first inputs for a fitter. The ruled/extruded
//! models (elliptical_arch, elliptical_cylinder, bspline_wave) arrive as
//! STRIP meshes — full-length quads with NO interior vertices along the
//! ruling, which is how real CAD exporters tessellate extrusions — kept
//! deliberately so the NURBS recognizer is built against realistic strip
//! meshes from day one. The two conic models get their strips from genuine
//! OCCT tessellation; the wave's strips are scripted because BRepMesh
//! internally grid-refines B-spline-based surfaces (interior nodes appear
//! at any deflection/angle setting). hill_with_hole is a plain OCCT
//! tessellation: one realistically irregular CAD-export-style mesh. Every
//! written STL is parsed back and validated (closed oriented manifold,
//! volume against the B-rep).
//!
//! In roughly increasing order of recognizer difficulty:
//! - elliptical_arch: half-ellipse profile extruded — the simplest freeform
//!   target: one non-periodic singly-curved patch with sharp boundaries
//!   (arriving as a strip mesh: 2 vertex rings along the extrusion).
//! - bspline_hill: doubly-curved dome patch on a slab (explicit pole
//!   lattice), the canonical height-field-style B-spline target.
//! - bspline_saddle: negative Gaussian curvature everywhere (hyperbolic-
//!   paraboloid-shaped pole net).
//! - bspline_wave: B-spline profile curve extruded along Y (strip mesh).
//!   Singly curved with alternating curvature sign; crest/trough radii
//!   ~3 mm sit in fillet territory and bait the cylinder detector into
//!   over-segmenting.
//! - bspline_terrain: bump + dip + ripple in one patch; wide curvature
//!   variation in sign and magnitude.
//! - bspline_blister: raised center on an OTHERWISE EXACTLY PLANAR top
//!   (outer pole rings all flat; B-spline local support makes the margin
//!   exactly planar). Tangent freeform→plane junction: the plane detector
//!   will absorb the flat margin and the boundary is ambiguous.
//! - hill_with_hole: bspline_hill minus an off-center vertical cylinder.
//!   Freeform × analytic interaction: inner loop on the B-spline face,
//!   cylinder∩B-spline edge curves, tessellated (irregular) mesh.
//! - elliptical_cylinder: full elliptical prism — a closed periodic
//!   freeform wall (stretch case: needs a periodic B-spline face + seam).

use brepper::read_step_file;
use opencascade_sys::{
    b_rep_algo_api, b_rep_check, b_rep_fillet_api, b_rep_g_prop, b_rep_mesh, b_rep_prim_api,
    b_rep_builder_api, b_rep_tools, ch_fi3d, g_prop, geom, gp, message, stl_api,
    n_collection, top_abs, top_exp, topo_ds, OwnedPtr,
};
use std::collections::HashMap;
use std::f64::consts::PI;

// ---------------------------------------------------------------------------
// Shared helpers
// ---------------------------------------------------------------------------

fn write_step(shape: &topo_ds::Shape, dir: &str, name: &str) {
    let step_path = format!("{dir}/{name}.step");
    brepper::write_step_file(shape, &step_path).expect("failed to write STEP");
    println!("wrote {step_path}");
}

/// Tessellate a shape and write its STL (used for occt/ models and for
/// nurbs models whose mesh should look like a real CAD export).
fn write_tessellated_stl(shape: &topo_ds::Shape, dir: &str, name: &str, deflection: f64, angle_deg: f64) {
    let _mesh = b_rep_mesh::IncrementalMesh::new_shape_real_bool_real_bool(
        shape, deflection, false, angle_deg.to_radians(), false,
    );
    let stl_path = format!("{dir}/{name}.stl");
    let mut writer = stl_api::Writer::new();
    *writer.ascii_mode() = false;
    let progress = message::ProgressRange::new();
    if !writer.write_shape_charptr_progressrange(shape, &stl_path, &progress) {
        panic!("failed to write {stl_path}");
    }
    println!("wrote {stl_path}");
}

fn write_pair(shape: &topo_ds::Shape, dir: &str, name: &str) {
    write_tessellated_stl(shape, dir, name, 0.1, 28.6);
    write_step(shape, dir, name);
}

/// Parse a just-written STL back and validate it: closed oriented manifold,
/// and enclosed volume within `rel_tol` of the closest of `brep_vols`
/// (tessellation loses volume to chordal sag, and the two BRepGProp
/// quadratures disagree on freeform faces — see shape_volume /
/// shape_volume_adaptive). Returns the mesh volume.
fn validate_written_stl(dir: &str, name: &str, brep_vols: &[f64], rel_tol: f64) -> f64 {
    let mesh = MeshBuf::from_stl_file(&format!("{dir}/{name}.stl"));
    mesh.assert_closed_manifold(name);
    let vol = mesh.volume();
    let best = brep_vols
        .iter()
        .copied()
        .min_by(|a, b| (a - vol).abs().total_cmp(&(b - vol).abs()))
        .expect("brep_vols must be non-empty");
    let rel = (vol - best).abs() / best;
    assert!(
        rel < rel_tol,
        "{name}: STL volume {vol:.3} vs B-rep volumes {brep_vols:?} (rel {rel:.2e})"
    );
    vol
}

/// Volume via fixed-order Gauss quadrature: exact for planar/conic faces
/// (used with tight asserts on the analytic prism models), but
/// under-integrates B-spline faces with localized bumps by ~0.5%.
fn shape_volume(shape: &topo_ds::Shape) -> f64 {
    let mut gprops = g_prop::GProps::new();
    b_rep_g_prop::volume_properties_shape_gprops_bool3(shape, &mut gprops, true, false, false);
    gprops.mass()
}

/// Volume via the adaptive-accuracy variant, for solids with B-spline faces.
/// (On conic faces this is oddly *less* accurate than the fixed-order path —
/// ~2e-5 relative on the elliptical prism — so it is only used where needed.)
fn shape_volume_adaptive(shape: &topo_ds::Shape) -> f64 {
    let mut gprops = g_prop::GProps::new();
    b_rep_g_prop::volume_properties_shape_gprops_real_bool2(shape, &mut gprops, 1e-7, true, false);
    gprops.mass()
}

/// Assert a generated reference shape passes BRepCheck.
fn assert_shape_valid(shape: &topo_ds::Shape, name: &str) {
    let analyzer = b_rep_check::Analyzer::new_shape_bool(shape, true);
    assert!(analyzer.is_valid(), "{name}: reference shape failed BRepCheck");
}

// ---------------------------------------------------------------------------
// Scripted STL meshes
//
// The nurbs meshes are built by sampling the exact surfaces, so vertices lie
// exactly on the reference geometry. MeshBuf validates that the result is a
// closed, consistently oriented manifold (every directed edge is matched by
// its reverse — this also catches T-junctions and flipped faces) and
// computes the enclosed volume for cross-checking against the B-rep.
// ---------------------------------------------------------------------------

struct MeshBuf {
    tris: Vec<[[f64; 3]; 3]>,
}

fn vkey(p: [f64; 3]) -> [u64; 3] {
    [p[0].to_bits(), p[1].to_bits(), p[2].to_bits()]
}

impl MeshBuf {
    fn new() -> Self {
        MeshBuf { tris: Vec::new() }
    }

    /// Parse a binary STL back into a MeshBuf (all triangles kept verbatim,
    /// bypassing `tri`'s zero-area filter, so closure checks see exactly
    /// what was written).
    fn from_stl_file(path: &str) -> Self {
        let data = std::fs::read(path).expect("failed to read back STL");
        assert!(data.len() >= 84, "{path}: truncated binary STL");
        let n = u32::from_le_bytes(data[80..84].try_into().unwrap()) as usize;
        assert_eq!(data.len(), 84 + 50 * n, "{path}: bad binary STL length");
        let mut tris = Vec::with_capacity(n);
        for t in 0..n {
            let base = 84 + 50 * t + 12; // skip the stored normal
            let mut tri = [[0.0f64; 3]; 3];
            for (v, corner) in tri.iter_mut().enumerate() {
                for (c, coord) in corner.iter_mut().enumerate() {
                    let off = base + 12 * v + 4 * c;
                    *coord = f32::from_le_bytes(data[off..off + 4].try_into().unwrap()) as f64;
                }
            }
            tris.push(tri);
        }
        MeshBuf { tris }
    }

    fn tri(&mut self, a: [f64; 3], b: [f64; 3], c: [f64; 3]) {
        let u = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
        let v = [c[0] - a[0], c[1] - a[1], c[2] - a[2]];
        let n = cross(u, v);
        // Chain-merged triangulations legitimately emit zero-area triangles
        // where two chains share a vertex; drop them.
        if norm(n) < 1e-12 {
            return;
        }
        self.tris.push([a, b, c]);
    }

    /// Quad split into two triangles; `alt` flips the diagonal so grid
    /// meshes alternate diagonals checkerboard-style.
    fn quad(&mut self, a: [f64; 3], b: [f64; 3], c: [f64; 3], d: [f64; 3], alt: bool) {
        if alt {
            self.tri(b, c, d);
            self.tri(b, d, a);
        } else {
            self.tri(a, b, c);
            self.tri(a, c, d);
        }
    }

    /// Signed volume via the divergence theorem. Positive for a closed mesh
    /// with outward-facing normals.
    fn volume(&self) -> f64 {
        let mut v6 = 0.0;
        for t in &self.tris {
            v6 += dot(t[0], cross(t[1], t[2]));
        }
        v6 / 6.0
    }

    /// Every directed edge must be matched by its reverse exactly once.
    /// Catches open boundaries, T-junctions, duplicated and flipped faces.
    fn assert_closed_manifold(&self, name: &str) {
        let mut edges: HashMap<([u64; 3], [u64; 3]), i64> = HashMap::new();
        for t in &self.tris {
            for (p, q) in [(t[0], t[1]), (t[1], t[2]), (t[2], t[0])] {
                let (kp, kq) = (vkey(p), vkey(q));
                // A zero-length edge (two corners rounding to the same f32
                // triple in a parsed STL) pairs with nothing; skip it — the
                // sliver's remaining edges still participate in closure.
                if kp == kq {
                    continue;
                }
                // +1 in lexicographically increasing direction, -1 reverse.
                if kp < kq {
                    *edges.entry((kp, kq)).or_insert(0) += 1;
                } else {
                    *edges.entry((kq, kp)).or_insert(0) -= 1;
                }
            }
        }
        let bad = edges.values().filter(|&&c| c != 0).count();
        assert_eq!(bad, 0, "{name}: mesh is not a closed oriented manifold ({bad} unmatched edges)");
    }

    fn write_stl(&self, dir: &str, name: &str) {
        let path = format!("{dir}/{name}.stl");
        let mut buf: Vec<u8> = Vec::with_capacity(84 + 50 * self.tris.len());
        let mut header = [0u8; 80];
        let tag = b"brepper gen-test-models scripted mesh";
        header[..tag.len()].copy_from_slice(tag);
        buf.extend_from_slice(&header);
        buf.extend_from_slice(&(self.tris.len() as u32).to_le_bytes());
        for t in &self.tris {
            let u = [t[1][0] - t[0][0], t[1][1] - t[0][1], t[1][2] - t[0][2]];
            let v = [t[2][0] - t[0][0], t[2][1] - t[0][1], t[2][2] - t[0][2]];
            let mut n = cross(u, v);
            let len = norm(n);
            n = [n[0] / len, n[1] / len, n[2] / len];
            for x in n {
                buf.extend_from_slice(&(x as f32).to_le_bytes());
            }
            for p in t {
                for x in p {
                    buf.extend_from_slice(&(*x as f32).to_le_bytes());
                }
            }
            buf.extend_from_slice(&0u16.to_le_bytes());
        }
        std::fs::write(&path, buf).expect("failed to write STL");
        println!("wrote {path} ({} triangles)", self.tris.len());
    }
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn norm(a: [f64; 3]) -> f64 {
    dot(a, a).sqrt()
}

// ---------------------------------------------------------------------------
// Exact freeform surface construction
// ---------------------------------------------------------------------------

/// Clamped-uniform cubic knot vector for `n` control points: distinct knots
/// 0..n-3 with end multiplicity 4.
fn clamped_cubic_knots(n: usize) -> (OwnedPtr<n_collection::Array1_double>, OwnedPtr<n_collection::Array1_int>) {
    let nk = n - 2; // distinct knots for degree 3
    let mut knots = n_collection::Array1_double::new_int2(1, nk as i32);
    let mut mults = n_collection::Array1_int::new_int2(1, nk as i32);
    for k in 0..nk {
        knots.set_value(k as i32 + 1, &(k as f64));
        let m = if k == 0 || k == nk - 1 { 4 } else { 1 };
        mults.set_value(k as i32 + 1, &m);
    }
    (knots, mults)
}

/// Exact bicubic B-spline surface from an n×n pole lattice on a regular XY
/// grid (spacing `d`, lower corner `x0`,`y0`), pole heights from `zf`.
/// Poles on a regular lattice make x depend only on u and y only on v.
fn lattice_bspline_surface(
    n: usize,
    x0: f64,
    y0: f64,
    d: f64,
    zf: &dyn Fn(f64, f64) -> f64,
) -> OwnedPtr<geom::Handle_Geom_BSplineSurface> {
    let mut poles = n_collection::Array2_gp_Pnt::new_int4(1, n as i32, 1, n as i32);
    for i in 0..n {
        for j in 0..n {
            let x = x0 + i as f64 * d;
            let y = y0 + j as f64 * d;
            let p = gp::Pnt::new_real3(x, y, zf(x, y));
            poles.set_value(i as i32 + 1, j as i32 + 1, &p);
        }
    }
    let (uknots, umults) = clamped_cubic_knots(n);
    let (vknots, vmults) = clamped_cubic_knots(n);
    let surf = geom::BSplineSurface::new_arr2_pnt_arr1_double2_arr1_int2_int2(
        &poles, &uknots, &vknots, &umults, &vmults, 3, 3,
    );
    geom::BSplineSurface::to_handle(surf)
}

/// Exact cubic B-spline curve in the XZ plane (y=0) from poles at
/// x = x0 + i*d, z from `zs`.
fn profile_bspline_curve(x0: f64, d: f64, zs: &[f64]) -> OwnedPtr<geom::Handle_Geom_BSplineCurve> {
    let n = zs.len();
    let mut poles = n_collection::Array1_gp_Pnt::new_int2(1, n as i32);
    for (i, z) in zs.iter().enumerate() {
        let p = gp::Pnt::new_real3(x0 + i as f64 * d, 0.0, *z);
        poles.set_value(i as i32 + 1, &p);
    }
    let (knots, mults) = clamped_cubic_knots(n);
    let curve = geom::BSplineCurve::new_arr1_pnt_arr1_double_arr1_int_int(
        &poles, &knots, &mults, 3,
    );
    geom::BSplineCurve::to_handle(curve)
}

/// Find the parameter where the monotone function `f` equals `target`.
fn bisect_param(f: &dyn Fn(f64) -> f64, mut lo: f64, mut hi: f64, target: f64) -> f64 {
    let (flo, fhi) = (f(lo), f(hi));
    assert!(
        (flo <= target && target <= fhi) || (fhi <= target && target <= flo),
        "bisect_param: target {target} not bracketed by [{flo}, {fhi}]"
    );
    let increasing = fhi >= flo;
    for _ in 0..100 {
        let mid = 0.5 * (lo + hi);
        if (f(mid) < target) == increasing {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    0.5 * (lo + hi)
}

// ---------------------------------------------------------------------------
// Slab-under-surface models (bspline_hill / saddle / terrain / blister)
//
// The exact surface overhangs the box footprint (pole lattice padded beyond
// it). The reference solid is Cut(box, Prism(face(surface), +Z)): the box
// material above the surface is removed, leaving a slab whose top is the
// exact freeform surface and whose walls/bottom are box planes. The STL
// samples the same surface on a regular parameter grid.
// ---------------------------------------------------------------------------

struct SlabSpec<'a> {
    name: &'a str,
    /// Exact surface evaluator (u,v) → point on the reference surface.
    eval: &'a dyn Fn(f64, f64) -> [f64; 3],
    /// Parameter windows corresponding to the box footprint [0,bx]×[0,by].
    u_range: (f64, f64),
    v_range: (f64, f64),
    /// Box dimensions.
    bx: f64,
    by: f64,
    bz: f64,
    /// STL sampling density (vertex counts along u and v).
    nu: usize,
    nv: usize,
}

/// Build the scripted STL for a slab model and return (mesh, top z range).
fn slab_mesh(spec: &SlabSpec) -> (MeshBuf, f64, f64) {
    let SlabSpec { eval, nu, nv, bx, by, .. } = *spec;
    let (u0, u1) = spec.u_range;
    let (v0, v1) = spec.v_range;

    // Top grid: vertices exactly on the surface, with boundary coordinates
    // snapped onto the box walls (the parameter windows solve to the walls
    // within ~1e-12; snapping makes the walls exactly planar).
    let mut top = vec![vec![[0.0f64; 3]; nv + 1]; nu + 1];
    let (mut zmin, mut zmax) = (f64::INFINITY, f64::NEG_INFINITY);
    for (i, row) in top.iter_mut().enumerate() {
        let u = u0 + (u1 - u0) * i as f64 / nu as f64;
        for (j, slot) in row.iter_mut().enumerate() {
            let v = v0 + (v1 - v0) * j as f64 / nv as f64;
            let mut p = eval(u, v);
            for (cond, coord, want) in [
                (i == 0, 0, 0.0),
                (i == nu, 0, bx),
                (j == 0, 1, 0.0),
                (j == nv, 1, by),
            ] {
                if cond {
                    assert!(
                        (p[coord] - want).abs() < 1e-6,
                        "{}: boundary sample off wall by {:.2e}",
                        spec.name,
                        (p[coord] - want).abs()
                    );
                    p[coord] = want;
                }
            }
            zmin = zmin.min(p[2]);
            zmax = zmax.max(p[2]);
            *slot = p;
        }
    }
    assert!(zmin > 0.5, "{}: surface dips too close to the box bottom (zmin {zmin})", spec.name);
    assert!(zmax < spec.bz - 0.5, "{}: surface too close to box top (zmax {zmax})", spec.name);

    let bot = |i: usize, j: usize| [top[i][j][0], top[i][j][1], 0.0];

    let mut m = MeshBuf::new();
    // Top (normals up).
    for i in 0..nu {
        for j in 0..nv {
            m.quad(top[i][j], top[i + 1][j], top[i + 1][j + 1], top[i][j + 1], (i + j) % 2 == 1);
        }
    }
    // Bottom (normals down).
    for i in 0..nu {
        for j in 0..nv {
            m.quad(bot(i, j), bot(i, j + 1), bot(i + 1, j + 1), bot(i + 1, j), (i + j) % 2 == 1);
        }
    }
    // Walls x=0 and x=bx (top edge follows the surface).
    for j in 0..nv {
        m.quad(top[0][j], top[0][j + 1], bot(0, j + 1), bot(0, j), false);
        m.quad(top[nu][j + 1], top[nu][j], bot(nu, j), bot(nu, j + 1), false);
    }
    // Walls y=0 and y=by.
    for i in 0..nu {
        m.quad(top[i + 1][0], top[i][0], bot(i, 0), bot(i + 1, 0), false);
        m.quad(top[i][nv], top[i + 1][nv], bot(i + 1, nv), bot(i, nv), false);
    }
    (m, zmin, zmax)
}

/// Prism the (overhanging) exact-surface face upward and cut it from the
/// box: the box material above the surface is removed, leaving a slab whose
/// top is the untouched exact surface and whose walls/bottom are box planes.
fn slab_solid(
    name: &str,
    face_shape: &topo_ds::Shape,
    dims: (f64, f64, f64),
) -> OwnedPtr<b_rep_algo_api::Cut> {
    let (bx, by, bz) = dims;
    let prism_vec = gp::Vec::new_real3(0.0, 0.0, bz + 30.0);
    let mut prism = b_rep_prim_api::MakePrism::new_shape_vec(face_shape, &prism_vec);
    let mut boxm = b_rep_prim_api::MakeBox::new_real3(bx, by, bz);
    let progress = message::ProgressRange::new();
    let mut cut =
        b_rep_algo_api::Cut::new_shape2_progressrange(boxm.shape(), prism.shape(), &progress);
    assert!(cut.is_done(), "{name}: box-minus-prism boolean failed");
    assert_shape_valid(cut.shape(), name);
    cut
}

/// Build the reference solid, cross-validate, and write the STL/STEP pair.
fn emit_slab_model(dir: &str, spec: &SlabSpec, face_shape: &topo_ds::Shape) {
    let (mesh, zmin, zmax) = slab_mesh(spec);
    mesh.assert_closed_manifold(spec.name);

    let mut cut = slab_solid(spec.name, face_shape, (spec.bx, spec.by, spec.bz));
    let solid = cut.shape();

    // Cross-check the B-rep against the mesh (whose vertices lie exactly on
    // the surface, with only chordal-sag error). Neither BRepGProp variant
    // is reliable on every freeform face — fixed-order Gauss is ~0.5% low on
    // the hill's localized bump, while the adaptive variant is ~4% HIGH on
    // the extrusion-of-B-spline solid — so accept the B-rep if either
    // independent quadrature agrees with the mesh.
    let stl_vol = mesh.volume();
    let fixed_vol = shape_volume(solid);
    let adaptive_vol = shape_volume_adaptive(solid);
    let brep_vol = if (fixed_vol - stl_vol).abs() < (adaptive_vol - stl_vol).abs() {
        fixed_vol
    } else {
        adaptive_vol
    };
    let rel = (brep_vol - stl_vol).abs() / brep_vol;
    assert!(
        rel < 0.01,
        "{}: STL volume {stl_vol:.3} vs B-rep volumes (fixed {fixed_vol:.3}, adaptive {adaptive_vol:.3})",
        spec.name
    );

    mesh.write_stl(dir, spec.name);
    write_step(solid, dir, spec.name);
    println!(
        "  {}: volume {brep_vol:.2} mm³ (mesh {stl_vol:.2}, rel diff {rel:.1e}), top z ∈ [{zmin:.2}, {zmax:.2}]",
        spec.name
    );
}

/// Common driver for the pole-lattice patch models.
fn lattice_slab_model(
    dir: &str,
    name: &str,
    n: usize,
    lattice_origin: f64,
    lattice_spacing: f64,
    zf: &dyn Fn(f64, f64) -> f64,
    boxdims: (f64, f64, f64),
) -> OwnedPtr<geom::Handle_Geom_BSplineSurface> {
    let (bx, by, bz) = boxdims;
    let handle = lattice_bspline_surface(n, lattice_origin, lattice_origin, lattice_spacing, zf);

    // Regular lattice ⇒ x depends only on u, y only on v. Solve the
    // parameter windows that map to the box footprint.
    let kmax = (n - 3) as f64;
    let vmid = kmax / 2.0;
    let x_of_u = |u: f64| handle.get().value(u, vmid).x();
    let y_of_v = |v: f64| handle.get().value(vmid, v).y();
    let u_range = (bisect_param(&x_of_u, 0.0, kmax, 0.0), bisect_param(&x_of_u, 0.0, kmax, bx));
    let v_range = (bisect_param(&y_of_v, 0.0, kmax, 0.0), bisect_param(&y_of_v, 0.0, kmax, by));

    let eval = |u: f64, v: f64| {
        let p = handle.get().value(u, v);
        [p.x(), p.y(), p.z()]
    };
    let spec = SlabSpec {
        name,
        eval: &eval,
        u_range,
        v_range,
        bx,
        by,
        bz,
        nu: 40,
        nv: 40,
    };
    let mut face =
        b_rep_builder_api::MakeFace::new_handlegeomsurface_real(&handle.to_handle_surface(), 1e-7);
    emit_slab_model(dir, &spec, face.shape());
    handle
}

// ---------------------------------------------------------------------------
// Re-mesh mode (Phase 0.2)
// ---------------------------------------------------------------------------
//
// `gen_test_models remesh [<out_dir>] [<step> ...]` loads existing reference
// STEPs and re-tessellates them with OCCT's BRepMesh at several densities,
// writing STLs next to a copy of the reference each `--compare` run needs.
//
// Why: the corpus is Parasolid-heavy (Onshape exports), and Parasolid
// tessellates by sweeping STRIPS — full-length quads with no interior
// vertices along a ruling. BRepMesh is constrained Delaunay instead, which
// breaks two properties the pipeline has quietly been able to rely on: a
// facet no longer spans two rulings of a cylinder (so a facet normal is no
// longer exactly perpendicular to the axis — the Phase 1.1 initializer
// caveat), and facets are no longer grid-aligned (so meridian normals on a
// torus are only approximate). Re-meshing the SAME references keeps the
// ground truth exact while changing only the provenance.
//
// Off-origin variants exist to make the float32 quantum testable. An STL
// stores coordinates as f32, so vertex precision is a *relative* quantity:
// max|coordinate| * 1.19e-7. At the origin-centred corpus sizes that is
// ~1e-6 mm; translated a metre away it is ~1e-4 mm, which is larger than the
// 5e-5 mm default --vertex-tolerance. Both the STL and a translated copy of
// the reference STEP are written so --compare still applies.

/// Tessellation settings for the re-mesh variants. These are inputs to the
/// mesh GENERATOR, not pipeline tolerances: they describe what a CAD exporter
/// was asked for, exactly as Onshape's coarse/medium/fine presets do.
/// Deflection is chordal deviation in mm; the angle bounds the turn between
/// adjacent facets.
const REMESH_DENSITIES: &[(&str, f64, f64)] = &[
    ("coarse", 0.5, 30.0),
    ("medium", 0.1, 15.0),
    ("fine", 0.02, 6.0),
];

/// Translation for the off-origin variants, ~1 m in length. Deliberately not
/// uniform across the axes and deliberately leaving z near the origin: the
/// f32 quantum of a *distance* is direction-weighted, so a part that is far
/// away in x and y but not in z is the informative case — and it is also what
/// a part positioned inside a real assembly looks like.
const REMESH_OFFSET: [f64; 3] = [600.0, 800.0, 0.0];

/// Reference STEPs re-meshed by default. Chosen for surface-type coverage
/// (cylinder, cone, sphere, torus, fillet networks, graded plane/hole pairs)
/// among models whose current STLs are all Parasolid strips.
const REMESH_MODELS: &[(&str, &str)] = &[
    ("tests/onshape/cylinder_10x30.step", "cylinder_10x30"),
    ("tests/onshape/sphere_25.step", "sphere_25"),
    ("tests/onshape/pill.step", "pill"),
    ("tests/onshape/plate_with_hole_100x50.step", "plate_with_hole_100x50"),
    ("tests/onshape/rounded_cube_10_r2.step", "rounded_cube_10_r2"),
    ("tests/onshape/cone_tangent_planes.step", "cone_tangent_planes"),
    ("tests/onshape/pipe_elbow_10.step", "pipe_elbow_10"),
    ("tests/onshape/stepped_block_toroidal_fillets.step", "stepped_block_toroidal"),
];

/// Tessellate `shape` at (deflection, angle) and write `<dir>/<name>.stl`,
/// then read it back and report closure, volume error against the B-rep, and
/// the f32 quantum implied by its coordinates.
fn write_remesh_variant(
    shape: &topo_ds::Shape,
    dir: &str,
    name: &str,
    deflection: f64,
    angle_deg: f64,
    brep_vol: f64,
) {
    // BRepMesh is *incremental*: it reuses any triangulation already cached on
    // the shape. Clean first or the second density silently returns the first.
    b_rep_tools::clean(shape, true);
    write_tessellated_stl(shape, dir, name, deflection, angle_deg);

    let mesh = MeshBuf::from_stl_file(&format!("{dir}/{name}.stl"));
    mesh.assert_closed_manifold(name);
    let vol = mesh.volume();
    let rel = (vol - brep_vol).abs() / brep_vol;
    let max_coord = mesh
        .tris
        .iter()
        .flatten()
        .flatten()
        .fold(0.0f64, |m, c| m.max(c.abs()));
    println!(
        "  {name}: {} tris, volume {vol:.4} vs B-rep {brep_vol:.4} (rel {rel:.2e}), \
         max|coord| {max_coord:.1} mm -> f32 quantum {:.2e} mm",
        mesh.tris.len(),
        max_coord * 1.192_092_9e-7,
    );
}

/// Re-mesh one reference STEP at every density, plus an off-origin variant.
fn remesh_model(out_dir: &str, step_path: &str, name: &str) {
    let shape = read_step_file(step_path)
        .unwrap_or_else(|e| panic!("failed to read {step_path}: {e}"));
    let brep_vol = shape_volume(&shape).abs();
    println!("remeshing {step_path} (B-rep volume {brep_vol:.4} mm³)");

    for (label, deflection, angle) in REMESH_DENSITIES {
        write_remesh_variant(
            &shape,
            out_dir,
            &format!("{name}_occt_{label}"),
            *deflection,
            *angle,
            brep_vol,
        );
    }

    // Off-origin: a translated reference STEP plus one medium tessellation of
    // it. Copy the geometry (rather than only setting a location) so the
    // written STEP carries genuinely translated surfaces.
    let mut trsf = gp::Trsf::new();
    let offset = gp::Vec::new_real3(REMESH_OFFSET[0], REMESH_OFFSET[1], REMESH_OFFSET[2]);
    trsf.set_translation_part(&offset);
    let mut xform = b_rep_builder_api::Transform::new_shape_trsf_bool(&shape, &trsf, true);
    let moved_name = format!("{name}_offset");
    {
        let moved = xform.as_modify_shape_mut().shape();
        write_step(moved, out_dir, &moved_name);
    }
    let moved = xform.as_modify_shape_mut().shape();
    let (_, deflection, angle) = REMESH_DENSITIES[1];
    write_remesh_variant(
        moved,
        out_dir,
        &format!("{moved_name}_occt_medium"),
        deflection,
        angle,
        brep_vol,
    );
}

fn remesh_corpus(out_dir: &str, models: &[(String, String)]) {
    std::fs::create_dir_all(out_dir).expect("failed to create remesh output dir");
    for (step_path, name) in models {
        remesh_model(out_dir, step_path, name);
    }

    // Record which reference each model came from. A re-mesh's reference lives
    // in another directory, and guessing it by searching tests/ for a matching
    // stem is ambiguous (plate_with_hole_100x50.step exists under both
    // fusion/ and onshape/) and filesystem-order dependent. tools/sweep.sh
    // reads this file, so the pairing is recorded rather than inferred.
    let manifest: String = models
        .iter()
        .map(|(step_path, name)| format!("{name}\t{step_path}\n"))
        .collect();
    let path = format!("{out_dir}/sources.tsv");
    std::fs::write(&path, format!("# model\treference STEP it was tessellated from\n{manifest}"))
        .expect("failed to write remesh sources manifest");
    println!("wrote {path}");
}

// ---------------------------------------------------------------------------
// Randomized variants (Phase 0.3)
// ---------------------------------------------------------------------------
//
// `gen_test_models random [<out_dir>] [<per_family>] [<seed>]` writes a seeded
// population of parameterized models: STL + STEP pair each, at a tessellation
// density that is also drawn.
//
// Why random: fourteen hand-watched models can be satisfied by machinery that
// happens to suit their particular radii and fillet sizes. A seeded population
// cannot be memorized — the parameters land where no constant was ever tuned.
// The families deliberately sweep the quantities the tuned constants care
// about: fillet radius against face size, torus r/R across the whole range up
// to and past `MAX_TORUS_RADIUS_RATIO`, cone half-angles including the
// near-cylindrical end, and mesh density independently of all of it.
//
// The seed and every drawn parameter go into `<out_dir>/params.tsv`, so a
// surprising row can be reproduced and read without re-deriving it.

/// An owned copy of a shape, so an intermediate boolean result can outlive the
/// builder that produced it. `Oriented` with the shape's own orientation is
/// the identity: TopoDS_Shape is a refcounted handle, so this shares the
/// underlying geometry rather than duplicating it.
fn own_shape(shape: &topo_ds::Shape) -> OwnedPtr<topo_ds::Shape> {
    let orientation = shape.orientation();
    shape.oriented(orientation)
}

/// SplitMix64. A deterministic generator lives here rather than in a
/// dependency so the corpus is reproducible from the seed alone, forever.
struct Rng(u64);

impl Rng {
    fn next_u64(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    /// Uniform in [0, 1).
    fn unit(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64
    }

    /// Uniform in [lo, hi), rounded to 0.01 mm so the recorded parameters read
    /// like something a person would have typed.
    fn range(&mut self, lo: f64, hi: f64) -> f64 {
        ((lo + (hi - lo) * self.unit()) * 100.0).round() / 100.0
    }

    fn int(&mut self, lo: i32, hi: i32) -> i32 {
        lo + (self.next_u64() % ((hi - lo + 1) as u64)) as i32
    }
}

/// Round a derived dimension to 0.01 mm, for the same reason `Rng::range`
/// does: the recorded parameters should read like typed input.
fn round2(x: f64) -> f64 {
    (x * 100.0).round() / 100.0
}

/// One randomized model: write the pair and return its params.tsv line.
fn emit_random(
    dir: &str,
    name: &str,
    shape: &topo_ds::Shape,
    deflection: f64,
    angle: f64,
    params: &str,
) -> String {
    assert_shape_valid(shape, name);
    let brep_vol = shape_volume(shape).abs();
    write_tessellated_stl(shape, dir, name, deflection, angle);
    write_step(shape, dir, name);
    // A coarse tessellation of a small-radius feature legitimately loses a few
    // percent of the volume, so this bound only catches a broken mesh.
    let mesh_vol = validate_written_stl(dir, name, &[brep_vol], 0.1);
    format!(
        "{name}\t{params}\tdeflection={deflection:.3}\tangle={angle:.1}\
         \tbrep_vol={brep_vol:.4}\tmesh_vol={mesh_vol:.4}\n"
    )
}

/// A box with every edge filleted: fillet radius against face size is the
/// ratio the flat-region seeding constants are sensitive to.
fn random_filleted_box(dir: &str, rng: &mut Rng, i: usize) -> String {
    let (dx, dy, dz) = (rng.range(8.0, 40.0), rng.range(8.0, 40.0), rng.range(8.0, 40.0));
    let min_dim = dx.min(dy).min(dz);
    // Up to a quarter of the smallest dimension: any larger and opposite
    // fillets on the same face meet, which OCCT rejects.
    let r = rng.range(0.4, min_dim * 0.24);
    let mut boxm = b_rep_prim_api::MakeBox::new_real3(dx, dy, dz);
    let mut fillet = b_rep_fillet_api::MakeFillet::new_shape_filletshape(
        boxm.shape(),
        ch_fi3d::FilletShape::Rational,
    );
    // A shape MAP, not TopExp_Explorer: the explorer visits an edge once per
    // adjoining face, so a solid's every edge would be added twice and the
    // fillet would fail as over-specified.
    let mut edges = n_collection::IndexedMap_TopoDS_Shape_TopTools_ShapeMapHasher::new();
    top_exp::map_shapes_shape_shapeenum_imap_shape(
        boxm.shape(),
        top_abs::ShapeEnum::Edge,
        &mut edges,
    );
    for k in 1..=edges.size() {
        fillet.add_real_edge(r, topo_ds::edge(edges.find_key(k)));
    }
    let name = format!("rand_filleted_box_{i}");
    // IsDone is only meaningful after the build, which Shape() triggers.
    let built = own_shape(fillet.shape());
    assert!(fillet.is_done(), "{name}: fillet failed (r={r})");
    let (deflection, angle) = (rng.range(0.02, 0.4), rng.range(8.0, 28.0));
    emit_random(
        dir,
        &name,
        &built,
        deflection,
        angle,
        &format!("box={dx}x{dy}x{dz}\tfillet_r={r}"),
    )
}

/// A stack of truncated cones sharing radii at the joins. Consecutive cones
/// with different half-angles make a sharp ring; equal half-angles make the
/// join tangent and invisible, which is the case `cone_tangent_planes` covers
/// by hand at exactly one set of angles.
fn random_cone_chain(dir: &str, rng: &mut Rng, i: usize) -> String {
    let n = rng.int(2, 4) as usize;
    let mut radii = vec![rng.range(4.0, 18.0)];
    let mut heights = Vec::new();
    for _ in 0..n {
        // 0.35..1.0 of the previous radius: keeps every frustum a real cone
        // (never a cylinder by accident) and the stack monotonically tapering.
        radii.push(round2((radii[radii.len() - 1] * rng.range(0.35, 1.0)).max(1.0)));
        heights.push(rng.range(4.0, 20.0));
    }
    let mut acc: Option<OwnedPtr<topo_ds::Shape>> = None;
    let mut z = 0.0;
    for k in 0..n {
        let axis = gp::Ax2::new_pnt_dir(
            &gp::Pnt::new_real3(0.0, 0.0, z),
            &gp::Dir::new_real3(0.0, 0.0, 1.0),
        );
        let mut cone =
            b_rep_prim_api::MakeCone::new_ax2_real3(&axis, radii[k], radii[k + 1], heights[k]);
        z += heights[k];
        acc = Some(match acc {
            None => own_shape(cone.shape()),
            Some(prev) => {
                let progress = message::ProgressRange::new();
                let mut fuse =
                    b_rep_algo_api::Fuse::new_shape2_progressrange(&prev, cone.shape(), &progress);
                assert!(fuse.is_done(), "rand_cone_chain_{i}: fuse failed");
                own_shape(fuse.shape())
            }
        });
    }
    let name = format!("rand_cone_chain_{i}");
    let (deflection, angle) = (rng.range(0.02, 0.4), rng.range(8.0, 28.0));
    let radii_s = radii.iter().map(|r| format!("{r}")).collect::<Vec<_>>().join(",");
    let heights_s = heights.iter().map(|h| format!("{h}")).collect::<Vec<_>>().join(",");
    emit_random(
        dir,
        &name,
        acc.as_ref().unwrap(),
        deflection,
        angle,
        &format!("cones={n}\tradii={radii_s}\theights={heights_s}"),
    )
}

/// A full torus with r/R drawn across the whole range, deliberately including
/// values above `MAX_TORUS_RADIUS_RATIO = 0.7` — a legitimate fat torus must
/// not be rejected by a guard aimed at spurious fits.
fn random_torus(dir: &str, rng: &mut Rng, i: usize) -> String {
    let major = rng.range(6.0, 30.0);
    let ratio = rng.range(0.08, 0.85);
    let minor = round2((major * ratio).max(0.5));
    let mut torus = b_rep_prim_api::MakeTorus::new_real2(major, minor);
    let name = format!("rand_torus_{i}");
    let (deflection, angle) = (rng.range(0.02, 0.4), rng.range(8.0, 28.0));
    emit_random(
        dir,
        &name,
        torus.shape(),
        deflection,
        angle,
        &format!("R={major}\tr={minor}\tratio={:.3}", minor / major),
    )
}

/// A tube: two coaxial cylinders and two annular end planes. The wall
/// thickness sweeps down to values where the two cylinders' supports are
/// narrow strips.
fn random_tube(dir: &str, rng: &mut Rng, i: usize) -> String {
    let outer = rng.range(5.0, 25.0);
    let inner = round2((outer * rng.range(0.25, 0.9)).max(1.0));
    let height = rng.range(5.0, 40.0);
    let mut outer_cyl = b_rep_prim_api::MakeCylinder::new_real2(outer, height);
    let mut inner_cyl = b_rep_prim_api::MakeCylinder::new_real2(inner, height);
    let progress = message::ProgressRange::new();
    let mut cut = b_rep_algo_api::Cut::new_shape2_progressrange(
        outer_cyl.shape(),
        inner_cyl.shape(),
        &progress,
    );
    let name = format!("rand_tube_{i}");
    assert!(cut.is_done(), "{name}: cut failed");
    let (deflection, angle) = (rng.range(0.02, 0.4), rng.range(8.0, 28.0));
    emit_random(
        dir,
        &name,
        cut.shape(),
        deflection,
        angle,
        &format!("outer_r={outer}\tinner_r={inner}\th={height}"),
    )
}

/// A sphere sitting on a cylinder, fused. Where the sphere radius approaches
/// the cylinder radius the junction becomes near-tangent, which is the
/// ambiguity the earn loop's branching is meant to survive.
fn random_ball_on_cylinder(dir: &str, rng: &mut Rng, i: usize) -> String {
    let cyl_r = rng.range(3.0, 15.0);
    let height = rng.range(6.0, 30.0);
    // 0.6..1.6 of the cylinder radius spans "ball narrower than the post"
    // through "ball nearly flush with it" to "ball overhanging".
    let mut ball_r = round2((cyl_r * rng.range(0.6, 1.6)).max(1.0));
    // Sink the centre so the sphere always intersects the cylinder's top;
    // centre_z + ball_r > height, so it also always protrudes.
    let sink = rng.range(0.1, 0.7);
    let mut centre_z = round2(height - ball_r * sink);
    // A big ball on a short post can CONTAIN the cylinder, and then the fused
    // solid is a plain sphere — a wasted harness row (measured: two of the
    // first three drawn were degenerate this way). The cylinder's bottom rim
    // is its farthest point from the centre, so keeping the ball clear of that
    // rim keeps the post visible. Shrink until it is.
    while centre_z > 0.0 && ball_r > 1.0 && ball_r >= 0.9 * (cyl_r * cyl_r + centre_z * centre_z).sqrt() {
        ball_r = round2(ball_r * 0.8);
        centre_z = round2(height - ball_r * sink);
    }
    let (ball_r, centre_z) = (ball_r, centre_z);
    let mut cyl = b_rep_prim_api::MakeCylinder::new_real2(cyl_r, height);
    let centre = gp::Pnt::new_real3(0.0, 0.0, centre_z);
    let mut ball = b_rep_prim_api::MakeSphere::new_pnt_real(&centre, ball_r);
    let progress = message::ProgressRange::new();
    let mut fuse =
        b_rep_algo_api::Fuse::new_shape2_progressrange(cyl.shape(), ball.shape(), &progress);
    let name = format!("rand_ball_on_cylinder_{i}");
    assert!(fuse.is_done(), "{name}: fuse failed");
    let (deflection, angle) = (rng.range(0.02, 0.4), rng.range(8.0, 28.0));
    emit_random(
        dir,
        &name,
        fuse.shape(),
        deflection,
        angle,
        &format!("cyl_r={cyl_r}\th={height}\tball_r={ball_r}\tball_z={centre_z:.2}"),
    )
}

fn random_corpus(out_dir: &str, per_family: usize, seed: u64) {
    std::fs::create_dir_all(out_dir).expect("failed to create random output dir");
    let mut rng = Rng(seed);
    let mut manifest = format!("# seeded by `gen-test-models random <dir> {per_family} {seed}`\n");
    for i in 0..per_family {
        manifest.push_str(&random_filleted_box(out_dir, &mut rng, i));
        manifest.push_str(&random_cone_chain(out_dir, &mut rng, i));
        manifest.push_str(&random_torus(out_dir, &mut rng, i));
        manifest.push_str(&random_tube(out_dir, &mut rng, i));
        manifest.push_str(&random_ball_on_cylinder(out_dir, &mut rng, i));
    }
    let path = format!("{out_dir}/params.tsv");
    std::fs::write(&path, manifest).expect("failed to write random params manifest");
    println!("wrote {path}");
}

// ---------------------------------------------------------------------------
// Models
// ---------------------------------------------------------------------------

fn hollow_cube(dir: &str) {
    let mut outer = b_rep_prim_api::MakeBox::new_real3(20.0, 20.0, 20.0);
    let inner_corner = gp::Pnt::new_real3(5.0, 5.0, 5.0);
    let mut inner = b_rep_prim_api::MakeBox::new_pnt_real3(&inner_corner, 10.0, 10.0, 10.0);

    let progress = message::ProgressRange::new();
    let mut cut = b_rep_algo_api::Cut::new_shape2_progressrange(
        outer.shape(),
        inner.shape(),
        &progress,
    );
    assert!(cut.is_done(), "hollow_cube boolean cut failed");
    write_pair(cut.shape(), dir, "hollow_cube");
}

fn full_torus(dir: &str) {
    let mut torus = b_rep_prim_api::MakeTorus::new_real2(10.0, 3.0);
    write_pair(torus.shape(), dir, "full_torus");
}

/// Half-elliptical vault (a=15, b=8) extruded 20 mm along +Y. Faces: one
/// freeform vault (linear extrusion of an elliptical arc), two planar
/// half-ellipse caps, one planar bottom.
fn elliptical_arch(dir: &str) {
    let (a, b, len) = (15.0, 8.0, 20.0);
    let center = gp::Pnt::new_real3(20.0, 0.0, 0.0);
    // XDir = -x̂ so the parameter-space upper half (sin t > 0) points to +z:
    // P(t) = (20 - a·cos t, 0, b·sin t), t ∈ [0, π] runs (5,0,0)→(35,0,0).
    let normal = gp::Dir::new_real3(0.0, 1.0, 0.0);
    let xdir = gp::Dir::new_real3(-1.0, 0.0, 0.0);
    let ax2 = gp::Ax2::new_pnt_dir2(&center, &normal, &xdir);
    let elips = gp::Elips::new_ax2_real2(&ax2, a, b);

    let mut arc = b_rep_builder_api::MakeEdge::new_elips_real2(&elips, 0.0, PI);
    let p_right = gp::Pnt::new_real3(35.0, 0.0, 0.0);
    let p_left = gp::Pnt::new_real3(5.0, 0.0, 0.0);
    let mut chord = b_rep_builder_api::MakeEdge::new_pnt2(&p_right, &p_left);
    let mut wire = b_rep_builder_api::MakeWire::new_edge2(arc.edge(), chord.edge());
    let mut face = b_rep_builder_api::MakeFace::new_wire_bool(wire.wire(), true);
    let sweep = gp::Vec::new_real3(0.0, len, 0.0);
    let mut prism = b_rep_prim_api::MakePrism::new_shape_vec(face.shape(), &sweep);
    let solid = prism.shape();
    assert_shape_valid(solid, "elliptical_arch");

    let vol = shape_volume(solid);
    let expected = 0.5 * PI * a * b * len;
    assert!(
        (vol - expected).abs() / expected < 1e-9,
        "elliptical_arch: volume {vol} vs analytic {expected}"
    );

    // Real CAD exporters subdivide an extrusion only along its profile, so
    // the vault arrives as full-length strips (2 vertex rings along Y); the
    // NURBS recognizer must cope with ruling-direction sampling this sparse.
    write_tessellated_stl(solid, dir, "elliptical_arch", 0.05, 10.0);
    let stl_vol = validate_written_stl(dir, "elliptical_arch", &[vol], 0.02);
    write_step(solid, dir, "elliptical_arch");
    println!("  elliptical_arch: volume {vol:.2} mm³ (mesh {stl_vol:.2})");
}

/// Full elliptical prism (a=12, b=7, h=18): closed periodic freeform wall.
fn elliptical_cylinder(dir: &str) {
    let (a, b, h) = (12.0, 7.0, 18.0);
    let (cx, cy) = (15.0, 10.0);
    let center = gp::Pnt::new_real3(cx, cy, 0.0);
    let normal = gp::Dir::new_real3(0.0, 0.0, 1.0);
    let xdir = gp::Dir::new_real3(1.0, 0.0, 0.0);
    let ax2 = gp::Ax2::new_pnt_dir2(&center, &normal, &xdir);
    let elips = gp::Elips::new_ax2_real2(&ax2, a, b);

    let mut edge = b_rep_builder_api::MakeEdge::new_elips(&elips);
    let mut wire = b_rep_builder_api::MakeWire::new_edge(edge.edge());
    let mut face = b_rep_builder_api::MakeFace::new_wire_bool(wire.wire(), true);
    let sweep = gp::Vec::new_real3(0.0, 0.0, h);
    let mut prism = b_rep_prim_api::MakePrism::new_shape_vec(face.shape(), &sweep);
    let solid = prism.shape();
    assert_shape_valid(solid, "elliptical_cylinder");

    let vol = shape_volume(solid);
    let expected = PI * a * b * h;
    assert!(
        (vol - expected).abs() / expected < 1e-9,
        "elliptical_cylinder: volume {vol} vs analytic {expected}"
    );

    // Strip mesh: the wall tessellates into full-height strips (2 vertex
    // rings), exactly how CAD exporters emit extrusions.
    write_tessellated_stl(solid, dir, "elliptical_cylinder", 0.05, 10.0);
    let stl_vol = validate_written_stl(dir, "elliptical_cylinder", &[vol], 0.02);
    write_step(solid, dir, "elliptical_cylinder");
    println!("  elliptical_cylinder: volume {vol:.2} mm³ (mesh {stl_vol:.2})");
}

fn bspline_hill(dir: &str) -> OwnedPtr<geom::Handle_Geom_BSplineSurface> {
    let zf = |x: f64, y: f64| {
        let r2 = (x - 20.0).powi(2) + (y - 20.0).powi(2);
        5.0 + 9.0 * (-r2 / (2.0 * 64.0)).exp()
    };
    lattice_slab_model(dir, "bspline_hill", 15, -8.0, 4.0, &zf, (40.0, 40.0, 16.0))
}

fn bspline_saddle(dir: &str) {
    let zf = |x: f64, y: f64| 12.0 + ((x - 20.0).powi(2) - (y - 20.0).powi(2)) / 55.0;
    lattice_slab_model(dir, "bspline_saddle", 15, -8.0, 4.0, &zf, (40.0, 40.0, 22.0));
}

fn bspline_terrain(dir: &str) {
    let zf = |x: f64, y: f64| {
        let bump = 6.0 * (-((x - 14.0).powi(2) + (y - 15.0).powi(2)) / 72.0).exp();
        let dip = -5.0 * (-((x - 27.0).powi(2) + (y - 26.0).powi(2)) / 98.0).exp();
        let ripple = 1.3 * (x / 3.2).sin() * (y / 4.1).cos();
        10.0 + bump + dip + ripple
    };
    lattice_slab_model(dir, "bspline_terrain", 15, -8.0, 4.0, &zf, (40.0, 40.0, 20.0));
}

fn bspline_blister(dir: &str) {
    // Raised poles only within r<9 of the center; every other pole is at
    // z=8 exactly. A raised pole's cubic basis reaches 6 mm past it in each
    // axis, so outside the union of those support boxes — Chebyshev distance
    // > 12 from the center — the surface is EXACTLY planar z=8: a tangent
    // freeform→plane junction with an ambiguous segmentation boundary.
    let zf = |x: f64, y: f64| {
        let r2 = (x - 20.0).powi(2) + (y - 20.0).powi(2);
        if r2 < 81.0 {
            let s = 1.0 - r2 / 81.0;
            8.0 + 6.0 * s * s
        } else {
            8.0
        }
    };
    let handle = lattice_slab_model(dir, "bspline_blister", 19, -7.0, 3.0, &zf, (40.0, 40.0, 16.0));

    // Verify the exact-flat annulus and a real bump.
    let kmax = 16.0;
    let vmid = kmax / 2.0;
    let x_of_u = |u: f64| handle.get().value(u, vmid).x();
    let y_of_v = |v: f64| handle.get().value(vmid, v).y();
    let mut max_flat_dev = 0.0f64;
    let mut x: f64 = 0.5;
    while x < 39.5 {
        let mut y: f64 = 0.5;
        while y < 39.5 {
            if (x - 20.0).abs().max((y - 20.0).abs()) > 12.5 {
                let u = bisect_param(&x_of_u, 0.0, kmax, x);
                let v = bisect_param(&y_of_v, 0.0, kmax, y);
                max_flat_dev = max_flat_dev.max((handle.get().value(u, v).z() - 8.0).abs());
            }
            y += 1.0;
        }
        x += 1.0;
    }
    assert!(
        max_flat_dev < 1e-9,
        "bspline_blister: flat margin deviates by {max_flat_dev:.2e}"
    );
    let u_c = bisect_param(&x_of_u, 0.0, kmax, 20.0);
    let v_c = bisect_param(&y_of_v, 0.0, kmax, 20.0);
    let peak = handle.get().value(u_c, v_c).z();
    assert!(peak > 10.5, "bspline_blister: bump too shallow (peak z {peak})");
    println!("  bspline_blister: flat-margin max deviation {max_flat_dev:.1e}, peak z {peak:.2}");
}

fn bspline_wave(dir: &str) {
    // Exact B-spline profile curve extruded along +Y: crest/trough radii in
    // fillet territory to bait cylinder over-segmentation. Built like the
    // arch — closed profile wire swept with MakePrism — so the reference is
    // a plain prism solid. The STL is a SCRIPTED strip mesh (top quads span
    // the full 30 mm of Y between 2 vertex rings): mainstream CAD exporters
    // emit exactly such strips for extrusions, but OCCT's own mesher cannot
    // be used here — BRepMesh internally grid-refines B-spline-based
    // surfaces (interior nodes appear even at 0.1 mm / 20° settings),
    // unlike the conic extrusions of elliptical_arch/elliptical_cylinder
    // which it tessellates into true strips.
    let (bx, by) = (40.0, 30.0);
    let pattern = [0.0, 1.0, 0.0, -1.0];
    let zs: Vec<f64> = (0..11).map(|i| 8.0 + 3.2 * pattern[i % 4]).collect();
    let curve = profile_bspline_curve(0.0, bx / 10.0, &zs);
    let kmax = (zs.len() - 3) as f64;

    // z-range sanity (z depends only on the curve parameter), and the exact
    // profile area by dense trapezoid integration of z·dx — an analytic-
    // quality volume reference independent of BRepGProp.
    let samples = 200_000usize;
    let (mut zmin, mut zmax, mut area) = (f64::INFINITY, f64::NEG_INFINITY, 0.0);
    let mut prev: Option<(f64, f64)> = None;
    for i in 0..=samples {
        let p = curve.get().value(kmax * i as f64 / samples as f64);
        let (x, z) = (p.x(), p.z());
        zmin = zmin.min(z);
        zmax = zmax.max(z);
        if let Some((px, pz)) = prev {
            area += 0.5 * (z + pz) * (x - px);
        }
        prev = Some((x, z));
    }
    assert!(zmin > 0.5, "bspline_wave: profile z dips to {zmin:.2}");

    // Closed profile wire in the XZ plane: wavy top from (0,0,8) to
    // (40,0,8), then down, across, and up three straight edges.
    let p_ne = gp::Pnt::new_real3(bx, 0.0, 8.0);
    let p_se = gp::Pnt::new_real3(bx, 0.0, 0.0);
    let p_sw = gp::Pnt::new_real3(0.0, 0.0, 0.0);
    let p_nw = gp::Pnt::new_real3(0.0, 0.0, 8.0);
    let mut top = b_rep_builder_api::MakeEdge::new_handlegeomcurve(&curve.to_handle_curve());
    let mut down = b_rep_builder_api::MakeEdge::new_pnt2(&p_ne, &p_se);
    let mut across = b_rep_builder_api::MakeEdge::new_pnt2(&p_se, &p_sw);
    let mut up = b_rep_builder_api::MakeEdge::new_pnt2(&p_sw, &p_nw);
    let mut wire = b_rep_builder_api::MakeWire::new_edge4(
        top.edge(),
        down.edge(),
        across.edge(),
        up.edge(),
    );
    let mut face = b_rep_builder_api::MakeFace::new_wire_bool(wire.wire(), true);
    let sweep = gp::Vec::new_real3(0.0, by, 0.0);
    let mut prism = b_rep_prim_api::MakePrism::new_shape_vec(face.shape(), &sweep);
    let solid = prism.shape();
    assert_shape_valid(solid, "bspline_wave");

    let expected = area * by;
    let vol = shape_volume(solid);
    // The profile integral is machine-exact; fixed-order Gauss carries
    // ~4e-5 relative error on the B-spline-extrusion face (quadrature
    // cells don't align with knot spans), so the bound is 1e-4, not 1e-9.
    assert!(
        (vol - expected).abs() / expected < 1e-4,
        "bspline_wave: volume {vol} vs profile integral {expected}"
    );

    // Scripted strip mesh: profile sampled at nu chords (~0.5 mm at the
    // crests, comparable to a fine CAD export), every quad spanning the
    // full extrusion length. Caps are triangulated as vertical slices on
    // the same x samples (a corner fan would poke above the troughs), the
    // bottom reuses those x samples, and each side wall is a single quad —
    // so every shared boundary vertex is bit-identical and the mesh closes.
    let nu = 80usize;
    let profile: Vec<[f64; 3]> = (0..=nu)
        .map(|k| {
            let p = curve.get().value(kmax * k as f64 / nu as f64);
            let mut pt = [p.x(), 0.0, p.z()];
            if k == 0 || k == nu {
                let want_x = if k == 0 { 0.0 } else { bx };
                assert!(
                    (pt[0] - want_x).abs() < 1e-9 && (pt[2] - 8.0).abs() < 1e-9,
                    "bspline_wave: profile endpoint off by {:.2e}",
                    (pt[0] - want_x).abs().max((pt[2] - 8.0).abs())
                );
                pt = [want_x, 0.0, 8.0];
            }
            pt
        })
        .collect();
    let mut m = MeshBuf::new();
    for k in 0..nu {
        let ([x0, _, z0], [x1, _, z1]) = (profile[k], profile[k + 1]);
        // Top strip (normal up, tilting with the wave).
        m.quad([x0, 0.0, z0], [x1, 0.0, z1], [x1, by, z1], [x0, by, z0], false);
        // Caps at y=0 (normal -y) and y=by (normal +y), as vertical slices.
        m.quad([x0, 0.0, z0], [x0, 0.0, 0.0], [x1, 0.0, 0.0], [x1, 0.0, z1], false);
        m.quad([x1, by, z1], [x1, by, 0.0], [x0, by, 0.0], [x0, by, z0], false);
        // Bottom strip (normal down).
        m.quad([x0, 0.0, 0.0], [x0, by, 0.0], [x1, by, 0.0], [x1, 0.0, 0.0], false);
    }
    // Side walls x=0 (normal -x) and x=bx (normal +x).
    m.quad([0.0, 0.0, 0.0], [0.0, 0.0, 8.0], [0.0, by, 8.0], [0.0, by, 0.0], false);
    m.quad([bx, 0.0, 0.0], [bx, by, 0.0], [bx, by, 8.0], [bx, 0.0, 8.0], false);
    m.assert_closed_manifold("bspline_wave");
    m.write_stl(dir, "bspline_wave");
    let stl_vol = validate_written_stl(dir, "bspline_wave", &[expected], 0.02);
    write_step(solid, dir, "bspline_wave");
    println!(
        "  bspline_wave: volume {vol:.2} mm³ (adaptive {:.2}, mesh {stl_vol:.2}), top z ∈ [{zmin:.2}, {zmax:.2}]",
        shape_volume_adaptive(solid)
    );
}

/// bspline_hill minus an off-center vertical cylinder: inner loop on the
/// B-spline face, cylinder∩B-spline edge curves, and a realistically
/// irregular tessellated mesh.
fn hill_with_hole(dir: &str, hill: &geom::Handle_Geom_BSplineSurface) {
    let (bx, by, bz) = (40.0, 40.0, 16.0);
    let mut face =
        b_rep_builder_api::MakeFace::new_handlegeomsurface_real(&hill.to_handle_surface(), 1e-7);
    let mut slab = slab_solid("hill_with_hole slab", face.shape(), (bx, by, bz));
    let slab_vol = shape_volume_adaptive(slab.shape());

    let hole_center = gp::Pnt::new_real3(28.0, 14.0, -2.0);
    let zdir = gp::Dir::new_real3(0.0, 0.0, 1.0);
    let ax2 = gp::Ax2::new_pnt_dir(&hole_center, &zdir);
    let mut cyl = b_rep_prim_api::MakeCylinder::new_ax2_real2(&ax2, 5.0, bz + 8.0);
    let progress2 = message::ProgressRange::new();
    let mut cut =
        b_rep_algo_api::Cut::new_shape2_progressrange(slab.shape(), cyl.shape(), &progress2);
    assert!(cut.is_done(), "hill_with_hole: hole boolean failed");
    let solid = cut.shape();
    assert_shape_valid(solid, "hill_with_hole");

    // The removed material is the surface height integrated over the disk.
    let kmax = 12.0;
    let vmid = kmax / 2.0;
    let x_of_u = |u: f64| hill.get().value(u, vmid).x();
    let y_of_v = |v: f64| hill.get().value(vmid, v).y();
    let step: f64 = 0.1;
    let mut removed = 0.0;
    let mut gx: f64 = 28.0 - 5.0 + step / 2.0;
    while gx < 33.0 {
        let mut gy = 14.0 - 5.0 + step / 2.0;
        while gy < 19.0 {
            if (gx - 28.0).powi(2) + (gy - 14.0).powi(2) < 25.0 {
                let u = bisect_param(&x_of_u, 0.0, kmax, gx);
                let v = bisect_param(&y_of_v, 0.0, kmax, gy);
                removed += hill.get().value(u, v).z() * step * step;
            }
            gy += step;
        }
        gx += step;
    }
    let vol = shape_volume_adaptive(solid);
    let rel = (slab_vol - vol - removed).abs() / removed;
    assert!(
        rel < 0.01,
        "hill_with_hole: removed volume {:.2} vs integral {removed:.2} (rel {rel:.2e})",
        slab_vol - vol
    );

    write_tessellated_stl(solid, dir, "hill_with_hole", 0.05, 10.0);
    validate_written_stl(dir, "hill_with_hole", &[shape_volume(solid), vol], 0.02);
    write_step(solid, dir, "hill_with_hole");
    println!("  hill_with_hole: volume {vol:.2} mm³ (hill {slab_vol:.2}, hole {removed:.2})");
}

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    if args.first().map(String::as_str) == Some("remesh") {
        let out_dir = args.get(1).cloned().unwrap_or_else(|| "tests/remesh".to_string());
        let models: Vec<(String, String)> = if args.len() > 2 {
            args[2..]
                .iter()
                .map(|p| {
                    let stem = std::path::Path::new(p)
                        .file_stem()
                        .expect("STEP path has no file name")
                        .to_string_lossy()
                        .to_string();
                    (p.clone(), stem)
                })
                .collect()
        } else {
            REMESH_MODELS
                .iter()
                .map(|(p, n)| (p.to_string(), n.to_string()))
                .collect()
        };
        remesh_corpus(&out_dir, &models);
        return;
    }
    if args.first().map(String::as_str) == Some("random") {
        let out_dir = args.get(1).cloned().unwrap_or_else(|| "tests/random".to_string());
        let per_family = args.get(2).map_or(3, |s| s.parse().expect("per_family must be a number"));
        let seed = args.get(3).map_or(20_260_729, |s| s.parse().expect("seed must be a number"));
        random_corpus(&out_dir, per_family, seed);
        return;
    }

    let occt_dir = args.first().cloned().unwrap_or_else(|| "tests/occt".to_string());
    let nurbs_dir = args.get(1).cloned().unwrap_or_else(|| "tests/nurbs".to_string());
    std::fs::create_dir_all(&occt_dir).expect("failed to create occt output dir");
    std::fs::create_dir_all(&nurbs_dir).expect("failed to create nurbs output dir");

    hollow_cube(&occt_dir);
    full_torus(&occt_dir);

    elliptical_arch(&nurbs_dir);
    elliptical_cylinder(&nurbs_dir);
    let hill = bspline_hill(&nurbs_dir);
    bspline_saddle(&nurbs_dir);
    bspline_terrain(&nurbs_dir);
    bspline_blister(&nurbs_dir);
    bspline_wave(&nurbs_dir);
    hill_with_hole(&nurbs_dir, &hill);
}
