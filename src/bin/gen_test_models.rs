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
//! (shared biases cancel) and merely approximate. STLs are derived from the
//! same exact geometry by direct sampling (vertices lie exactly on the
//! reference surface), except hill_with_hole which uses OCCT tessellation
//! to provide one realistically irregular CAD-export-style mesh.
//!
//! In roughly increasing order of recognizer difficulty:
//! - elliptical_arch: half-ellipse profile extruded — the simplest freeform
//!   target: one non-periodic singly-curved patch with sharp boundaries.
//! - bspline_hill: doubly-curved dome patch on a slab (explicit pole
//!   lattice), the canonical height-field-style B-spline target.
//! - bspline_saddle: negative Gaussian curvature everywhere (hyperbolic-
//!   paraboloid-shaped pole net).
//! - bspline_wave: B-spline profile curve extruded along Y. Singly curved
//!   with alternating curvature sign; crest/trough radii ~3 mm sit in
//!   fillet territory and bait the cylinder detector into over-segmenting.
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

use opencascade_sys::{
    b_rep_algo_api, b_rep_check, b_rep_g_prop, b_rep_mesh, b_rep_prim_api, b_rep_builder_api,
    g_prop, geom, gp, message, stl_api, t_col_std, t_colgp, topo_ds, OwnedPtr,
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
    if !writer.write(shape, &stl_path, &progress) {
        panic!("failed to write {stl_path}");
    }
    println!("wrote {stl_path}");
}

fn write_pair(shape: &topo_ds::Shape, dir: &str, name: &str) {
    write_tessellated_stl(shape, dir, name, 0.1, 28.6);
    write_step(shape, dir, name);
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

    fn tri_oriented(&mut self, a: [f64; 3], b: [f64; 3], c: [f64; 3], flip: bool) {
        if flip {
            self.tri(a, c, b);
        } else {
            self.tri(a, b, c);
        }
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
                assert_ne!(kp, kq, "{name}: degenerate edge in mesh");
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

/// Triangulate the region between two x-monotone chains that share their
/// first and last points (e.g. an elliptical arc and its chord). Emits
/// triangles wound (upper, upper_next, lower) / (upper, lower_next, lower);
/// pass `flip` to reverse.
fn strip_between_chains(m: &mut MeshBuf, upper: &[[f64; 3]], lower: &[[f64; 3]], flip: bool) {
    let (mut i, mut j) = (0usize, 0usize);
    while i + 1 < upper.len() || j + 1 < lower.len() {
        let advance_upper = if j + 1 >= lower.len() {
            true
        } else if i + 1 >= upper.len() {
            false
        } else {
            upper[i + 1][0] <= lower[j + 1][0]
        };
        if advance_upper {
            m.tri_oriented(upper[i], upper[i + 1], lower[j], flip);
            i += 1;
        } else {
            m.tri_oriented(upper[i], lower[j + 1], lower[j], flip);
            j += 1;
        }
    }
}

// ---------------------------------------------------------------------------
// Exact freeform surface construction
// ---------------------------------------------------------------------------

/// Clamped-uniform cubic knot vector for `n` control points: distinct knots
/// 0..n-3 with end multiplicity 4.
fn clamped_cubic_knots(n: usize) -> (OwnedPtr<t_col_std::Array1OfReal>, OwnedPtr<t_col_std::Array1OfInteger>) {
    let nk = n - 2; // distinct knots for degree 3
    let mut knots = t_col_std::Array1OfReal::new_int2(1, nk as i32);
    let mut mults = t_col_std::Array1OfInteger::new_int2(1, nk as i32);
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
) -> OwnedPtr<geom::HandleGeomBSplineSurface> {
    let mut poles = t_colgp::Array2OfPnt::new_int4(1, n as i32, 1, n as i32);
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
    let surf = geom::BSplineSurface::new_array2ofpnt_array1ofreal2_array1ofinteger2_int2(
        &poles, &uknots, &vknots, &umults, &vmults, 3, 3,
    );
    geom::BSplineSurface::to_handle(surf)
}

/// Exact cubic B-spline curve in the XZ plane (y=0) from poles at
/// x = x0 + i*d, z from `zs`.
fn profile_bspline_curve(x0: f64, d: f64, zs: &[f64]) -> OwnedPtr<geom::HandleGeomBSplineCurve> {
    let n = zs.len();
    let mut poles = t_colgp::Array1OfPnt::new_int2(1, n as i32);
    for (i, z) in zs.iter().enumerate() {
        let p = gp::Pnt::new_real3(x0 + i as f64 * d, 0.0, *z);
        poles.set_value(i as i32 + 1, &p);
    }
    let (knots, mults) = clamped_cubic_knots(n);
    let curve = geom::BSplineCurve::new_array1ofpnt_array1ofreal_array1ofinteger_int(
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
// Slab-under-surface models (bspline_hill / saddle / terrain / blister / wave)
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

/// Trapezoid-rule volume under the sampled top grid (for cross-checking the
/// B-rep volume; grid x/y spacing is slightly non-uniform in x,y).
fn mesh_volume_estimate(spec: &SlabSpec, mesh: &MeshBuf) -> f64 {
    let _ = spec;
    mesh.volume()
}

/// Build the reference solid, cross-validate, and write the STL/STEP pair.
fn emit_slab_model(dir: &str, spec: &SlabSpec, face_shape: &topo_ds::Shape) {
    let (mesh, zmin, zmax) = slab_mesh(spec);
    mesh.assert_closed_manifold(spec.name);

    let prism_vec = gp::Vec::new_real3(0.0, 0.0, spec.bz + 30.0);
    let mut prism = b_rep_prim_api::MakePrism::new_shape_vec(face_shape, &prism_vec);
    let mut boxm = b_rep_prim_api::MakeBox::new_real3(spec.bx, spec.by, spec.bz);
    let progress = message::ProgressRange::new();
    let mut cut =
        b_rep_algo_api::Cut::new_shape2_progressrange(boxm.shape(), prism.shape(), &progress);
    assert!(cut.is_done(), "{}: box-minus-prism boolean failed", spec.name);
    let solid = cut.shape();
    assert_shape_valid(solid, spec.name);

    // Cross-check the B-rep against the mesh (whose vertices lie exactly on
    // the surface, with only chordal-sag error). Neither BRepGProp variant
    // is reliable on every freeform face — fixed-order Gauss is ~0.5% low on
    // the hill's localized bump, while the adaptive variant is ~4% HIGH on
    // the extrusion-of-B-spline solid — so accept the B-rep if either
    // independent quadrature agrees with the mesh.
    let stl_vol = mesh_volume_estimate(spec, &mesh);
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
) -> OwnedPtr<geom::HandleGeomBSplineSurface> {
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

    // Scripted mesh: vault rings sampled on the exact ellipse.
    let (nt, ny, ncx) = (48usize, 20usize, 30usize);
    let ring = |k: usize| -> [f64; 2] {
        // (x, z) of ring k; endpoints snapped exactly onto the chord plane.
        if k == 0 {
            return [5.0, 0.0];
        }
        if k == nt {
            return [35.0, 0.0];
        }
        let t = PI * k as f64 / nt as f64;
        [20.0 - a * t.cos(), b * t.sin()]
    };
    let ys: Vec<f64> = (0..=ny).map(|j| len * j as f64 / ny as f64).collect();
    let mut m = MeshBuf::new();
    // Vault: outward normal ≈ P'(t) × ŷ.
    for k in 0..nt {
        let ([x0, z0], [x1, z1]) = (ring(k), ring(k + 1));
        for j in 0..ny {
            m.quad(
                [x0, ys[j], z0],
                [x1, ys[j], z1],
                [x1, ys[j + 1], z1],
                [x0, ys[j + 1], z0],
                false,
            );
        }
    }
    // Caps: strip between the arc chain and the chord chain.
    let xs_chord: Vec<f64> = (0..=ncx).map(|i| 5.0 + 30.0 * i as f64 / ncx as f64).collect();
    for (y, flip) in [(0.0, true), (len, false)] {
        let arc_chain: Vec<[f64; 3]> = (0..=nt)
            .map(|k| {
                let [x, z] = ring(k);
                [x, y, z]
            })
            .collect();
        let chord_chain: Vec<[f64; 3]> = xs_chord.iter().map(|&x| [x, y, 0.0]).collect();
        strip_between_chains(&mut m, &arc_chain, &chord_chain, flip);
    }
    // Bottom: grid at z=0, normals down. Columns reuse chord xs, rows the
    // vault ys, so boundaries match the caps and vault exactly.
    for i in 0..ncx {
        for j in 0..ny {
            m.quad(
                [xs_chord[i], ys[j], 0.0],
                [xs_chord[i], ys[j + 1], 0.0],
                [xs_chord[i + 1], ys[j + 1], 0.0],
                [xs_chord[i + 1], ys[j], 0.0],
                false,
            );
        }
    }
    m.assert_closed_manifold("elliptical_arch");
    let stl_vol = m.volume();
    assert!(
        (stl_vol - vol).abs() / vol < 0.01,
        "elliptical_arch: STL volume {stl_vol} vs B-rep {vol}"
    );
    m.write_stl(dir, "elliptical_arch");
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

    let (nt, nz) = (76usize, 10usize);
    let ring = |k: usize| -> [f64; 2] {
        let t = 2.0 * PI * (k % nt) as f64 / nt as f64;
        [cx + a * t.cos(), cy + b * t.sin()]
    };
    let zs: Vec<f64> = (0..=nz).map(|j| h * j as f64 / nz as f64).collect();
    let mut m = MeshBuf::new();
    for k in 0..nt {
        let ([x0, y0], [x1, y1]) = (ring(k), ring(k + 1));
        for j in 0..nz {
            m.quad(
                [x0, y0, zs[j]],
                [x1, y1, zs[j]],
                [x1, y1, zs[j + 1]],
                [x0, y0, zs[j + 1]],
                false,
            );
        }
    }
    for (z, flip) in [(h, false), (0.0, true)] {
        let c = [cx, cy, z];
        for k in 0..nt {
            let ([x0, y0], [x1, y1]) = (ring(k), ring(k + 1));
            m.tri_oriented(c, [x0, y0, z], [x1, y1, z], flip);
        }
    }
    m.assert_closed_manifold("elliptical_cylinder");
    let stl_vol = m.volume();
    assert!(
        (stl_vol - vol).abs() / vol < 0.01,
        "elliptical_cylinder: STL volume {stl_vol} vs B-rep {vol}"
    );
    m.write_stl(dir, "elliptical_cylinder");
    write_step(solid, dir, "elliptical_cylinder");
    println!("  elliptical_cylinder: volume {vol:.2} mm³ (mesh {stl_vol:.2})");
}

fn bspline_hill(dir: &str) -> OwnedPtr<geom::HandleGeomBSplineSurface> {
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
    // fillet territory to bait cylinder over-segmentation.
    let pattern = [0.0, 1.0, 0.0, -1.0];
    let zs: Vec<f64> = (0..15).map(|i| 8.0 + 3.2 * pattern[i % 4]).collect();
    let curve = profile_bspline_curve(-8.0, 4.0, &zs);
    let ydir = gp::Dir::new_real3(0.0, 1.0, 0.0);
    let surf = geom::SurfaceOfLinearExtrusion::new_handlegeomcurve_dir(
        &curve.to_handle_curve(),
        &ydir,
    );
    let handle = geom::SurfaceOfLinearExtrusion::to_handle(surf);

    let (bx, by, bz) = (40.0, 30.0, 12.0);
    let kmax = 12.0;
    let x_of_u = |u: f64| curve.get().value(u).x();
    let u_range = (bisect_param(&x_of_u, 0.0, kmax, 0.0), bisect_param(&x_of_u, 0.0, kmax, bx));
    // Extrusion parameterization: (u, v) → curve(u) + v·ŷ, so v is y.
    let eval = |u: f64, v: f64| {
        let p = curve.get().value(u);
        [p.x(), v, p.z()]
    };
    let spec = SlabSpec {
        name: "bspline_wave",
        eval: &eval,
        u_range,
        v_range: (0.0, by),
        bx,
        by,
        bz,
        nu: 80,
        nv: 15,
    };
    let mut face = b_rep_builder_api::MakeFace::new_handlegeomsurface_real5(
        &handle.to_handle_surface(),
        0.0,
        kmax,
        -5.0,
        by + 5.0,
        1e-7,
    );
    emit_slab_model(dir, &spec, face.shape());
}

/// bspline_hill minus an off-center vertical cylinder: inner loop on the
/// B-spline face, cylinder∩B-spline edge curves, and a realistically
/// irregular tessellated mesh.
fn hill_with_hole(dir: &str, hill: &geom::HandleGeomBSplineSurface) {
    let (bx, by, bz) = (40.0, 40.0, 16.0);
    let mut face =
        b_rep_builder_api::MakeFace::new_handlegeomsurface_real(&hill.to_handle_surface(), 1e-7);
    let prism_vec = gp::Vec::new_real3(0.0, 0.0, bz + 30.0);
    let mut prism = b_rep_prim_api::MakePrism::new_shape_vec(face.shape(), &prism_vec);
    let mut boxm = b_rep_prim_api::MakeBox::new_real3(bx, by, bz);
    let progress = message::ProgressRange::new();
    let mut slab =
        b_rep_algo_api::Cut::new_shape2_progressrange(boxm.shape(), prism.shape(), &progress);
    assert!(slab.is_done(), "hill_with_hole: slab boolean failed");
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
    write_step(solid, dir, "hill_with_hole");
    println!("  hill_with_hole: volume {vol:.2} mm³ (hill {slab_vol:.2}, hole {removed:.2})");
}

fn main() {
    let occt_dir = std::env::args().nth(1).unwrap_or_else(|| "tests/occt".to_string());
    let nurbs_dir = std::env::args().nth(2).unwrap_or_else(|| "tests/nurbs".to_string());
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
