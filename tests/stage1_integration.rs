//! Integration tests for Stage 1: mesh input and validation.
//!
//! Tests that all STL/STEP test file pairs pass stage1 with consistency checks
//! and --compare validation, and that intentionally bad files fail correctly.

use brepper::config;
use brepper::stage1;

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

fn manifest_dir() -> String {
    env!("CARGO_MANIFEST_DIR").to_string()
}

fn config_for_stl(stl_path: &str) -> config::Config {
    config::parse_config_from(["brepper", stl_path, "--stage=1", "-q"]).unwrap()
}

fn config_for_compare(stl_path: &str, step_path: &str) -> config::Config {
    let mut config =
        config::parse_config_from(["brepper", stl_path, "--compare", step_path, "--stage=1", "-q"])
            .unwrap();
    config.load_compare_step().unwrap();
    config
}

// ---------------------------------------------------------------------------
// Macros for generating test cases
// ---------------------------------------------------------------------------

/// Generate a test that runs stage1 (read + validate) on a given STL file.
macro_rules! test_stl_stage1 {
    ($name:ident, $rel_path:literal) => {
        #[test]
        fn $name() {
            let stl = format!("{}/{}", manifest_dir(), $rel_path);
            let config = config_for_stl(&stl);
            stage1::stage1(&config).expect("stage1 should pass");
        }
    };
}

/// Generate a test that runs stage1 with --compare against a STEP file.
macro_rules! test_stl_step_compare {
    ($name:ident, $stl_path:literal, $step_path:literal) => {
        #[test]
        fn $name() {
            let stl = format!("{}/{}", manifest_dir(), $stl_path);
            let step = format!("{}/{}", manifest_dir(), $step_path);
            let config = config_for_compare(&stl, &step);
            stage1::stage1(&config).expect("stage1 with --compare should pass");
        }
    };
}

// ===========================================================================
// Good cases: all test STL files should pass stage1 validation
// ===========================================================================

// --- tests/manual/ ---
test_stl_stage1!(manual_cube_stage1, "tests/manual/cube.stl");

// --- tests/ccad/generated/ ---
test_stl_stage1!(ccad_cube_stage1, "tests/ccad/generated/cube.stl");
test_stl_stage1!(ccad_wedge_stage1, "tests/ccad/generated/wedge.stl");
test_stl_stage1!(ccad_t_shape_stage1, "tests/ccad/generated/t_shape.stl");
test_stl_stage1!(ccad_staircase_stage1, "tests/ccad/generated/staircase.stl");
test_stl_stage1!(ccad_simple_cylinder_stage1, "tests/ccad/generated/simple_cylinder.stl");
test_stl_stage1!(ccad_block_with_hole_stage1, "tests/ccad/generated/block_with_hole.stl");
test_stl_stage1!(ccad_pipe_stage1, "tests/ccad/generated/pipe.stl");
test_stl_stage1!(ccad_stepped_cylinder_stage1, "tests/ccad/generated/stepped_cylinder.stl");
test_stl_stage1!(ccad_two_holes_stage1, "tests/ccad/generated/two_holes.stl");
test_stl_stage1!(ccad_simple_sphere_stage1, "tests/ccad/generated/simple_sphere.stl");
test_stl_stage1!(ccad_hemisphere_stage1, "tests/ccad/generated/hemisphere.stl");
test_stl_stage1!(ccad_spherical_pocket_stage1, "tests/ccad/generated/spherical_pocket.stl");
test_stl_stage1!(ccad_ball_on_cylinder_stage1, "tests/ccad/generated/ball_on_cylinder.stl");

// --- tests/onshape/ ---
test_stl_stage1!(
    onshape_chamfered_cube_stage1,
    "tests/onshape/chamfered_cube_10_c1_medium.stl"
);
test_stl_stage1!(onshape_cone_stage1, "tests/onshape/cone_15x20_medium.stl");
test_stl_stage1!(
    onshape_cylinder_stage1,
    "tests/onshape/cylinder_10x30_medium.stl"
);
test_stl_stage1!(
    onshape_dome_hemisphere_stage1,
    "tests/onshape/dome_hemisphere_20_fine.stl"
);
test_stl_stage1!(
    onshape_l_bracket_stage1,
    "tests/onshape/l_bracket_simple_medium.stl"
);
test_stl_stage1!(
    onshape_pipe_elbow_stage1,
    "tests/onshape/pipe_elbow_10_fine.stl"
);
test_stl_stage1!(
    onshape_plate_with_hole_stage1,
    "tests/onshape/plate_with_hole_100x50_coarse.stl"
);
test_stl_stage1!(
    onshape_rounded_cube_fine_stage1,
    "tests/onshape/rounded_cube_10_r2_fine.stl"
);
test_stl_stage1!(
    onshape_rounded_cube_medium_stage1,
    "tests/onshape/rounded_cube_10_r2_medium.stl"
);
test_stl_stage1!(
    onshape_rounded_cube_coarse_stage1,
    "tests/onshape/rounded_cube_10_r2_coarse.stl"
);
test_stl_stage1!(
    onshape_sphere_stage1,
    "tests/onshape/sphere_25_fine.stl"
);
test_stl_stage1!(
    onshape_stepped_block_stage1,
    "tests/onshape/stepped_block_coarse.stl"
);

// --- tests/fusion/ ---
test_stl_stage1!(
    fusion_plate_high_stage1,
    "tests/fusion/plate_with_hole_100x50_high.stl"
);
test_stl_stage1!(
    fusion_plate_low_stage1,
    "tests/fusion/plate_with_hole_100x50_low.stl"
);
test_stl_stage1!(
    fusion_plate_medium_stage1,
    "tests/fusion/plate_with_hole_100x50_medium.stl"
);

// --- tests/nurbs/ (freeform models; scripted grid and tessellated strip
// meshes must both be valid solids) ---
test_stl_stage1!(nurbs_elliptical_arch_stage1, "tests/nurbs/elliptical_arch.stl");
test_stl_stage1!(
    nurbs_elliptical_cylinder_stage1,
    "tests/nurbs/elliptical_cylinder.stl"
);
test_stl_stage1!(nurbs_bspline_hill_stage1, "tests/nurbs/bspline_hill.stl");
test_stl_stage1!(nurbs_bspline_saddle_stage1, "tests/nurbs/bspline_saddle.stl");
test_stl_stage1!(nurbs_bspline_terrain_stage1, "tests/nurbs/bspline_terrain.stl");
test_stl_stage1!(nurbs_bspline_blister_stage1, "tests/nurbs/bspline_blister.stl");
test_stl_stage1!(nurbs_bspline_wave_stage1, "tests/nurbs/bspline_wave.stl");
test_stl_stage1!(nurbs_hill_with_hole_stage1, "tests/nurbs/hill_with_hole.stl");

// ===========================================================================
// Good cases: all test STL/STEP pairs should pass --compare
// ===========================================================================

// --- tests/manual/ ---
test_stl_step_compare!(
    manual_cube_compare,
    "tests/manual/cube.stl",
    "tests/manual/cube.step"
);

// --- tests/ccad/generated/ ---
test_stl_step_compare!(
    ccad_cube_compare,
    "tests/ccad/generated/cube.stl",
    "tests/ccad/generated/cube.step"
);
test_stl_step_compare!(
    ccad_wedge_compare,
    "tests/ccad/generated/wedge.stl",
    "tests/ccad/generated/wedge.step"
);
test_stl_step_compare!(
    ccad_t_shape_compare,
    "tests/ccad/generated/t_shape.stl",
    "tests/ccad/generated/t_shape.step"
);
test_stl_step_compare!(
    ccad_staircase_compare,
    "tests/ccad/generated/staircase.stl",
    "tests/ccad/generated/staircase.step"
);
test_stl_step_compare!(
    ccad_simple_cylinder_compare,
    "tests/ccad/generated/simple_cylinder.stl",
    "tests/ccad/generated/simple_cylinder.step"
);
test_stl_step_compare!(
    ccad_block_with_hole_compare,
    "tests/ccad/generated/block_with_hole.stl",
    "tests/ccad/generated/block_with_hole.step"
);
test_stl_step_compare!(
    ccad_pipe_compare,
    "tests/ccad/generated/pipe.stl",
    "tests/ccad/generated/pipe.step"
);
test_stl_step_compare!(
    ccad_stepped_cylinder_compare,
    "tests/ccad/generated/stepped_cylinder.stl",
    "tests/ccad/generated/stepped_cylinder.step"
);
test_stl_step_compare!(
    ccad_two_holes_compare,
    "tests/ccad/generated/two_holes.stl",
    "tests/ccad/generated/two_holes.step"
);
test_stl_step_compare!(
    ccad_simple_sphere_compare,
    "tests/ccad/generated/simple_sphere.stl",
    "tests/ccad/generated/simple_sphere.step"
);
test_stl_step_compare!(
    ccad_hemisphere_compare,
    "tests/ccad/generated/hemisphere.stl",
    "tests/ccad/generated/hemisphere.step"
);
test_stl_step_compare!(
    ccad_spherical_pocket_compare,
    "tests/ccad/generated/spherical_pocket.stl",
    "tests/ccad/generated/spherical_pocket.step"
);
test_stl_step_compare!(
    ccad_ball_on_cylinder_compare,
    "tests/ccad/generated/ball_on_cylinder.stl",
    "tests/ccad/generated/ball_on_cylinder.step"
);

// --- tests/onshape/ ---
test_stl_step_compare!(
    onshape_chamfered_cube_compare,
    "tests/onshape/chamfered_cube_10_c1_medium.stl",
    "tests/onshape/chamfered_cube_10_c1.step"
);
// NOTE (Phase 0.1 audit, 2026-07-29): cone_15x20.step is a CORPUS DEFECT — the
// reference is itself faceted, 59 planar faces rather than a cone plus a base
// plane, and its volume (1175.794 mm³) is the inscribed polyhedron's, 1.96e-3
// below the true πr²h/3 = 1178.097 mm³. brepper fits the cone exactly (2
// surfaces, 1178.097251 mm³) and therefore "loses" against it by that 1.96e-3.
// Never tune the pipeline toward this reference; use it only as a stage-1
// vertices-on-surface check, which is valid because the STL was tessellated
// from the same faceted body.
test_stl_step_compare!(
    onshape_cone_compare,
    "tests/onshape/cone_15x20_medium.stl",
    "tests/onshape/cone_15x20.step"
);
test_stl_step_compare!(
    onshape_cylinder_compare,
    "tests/onshape/cylinder_10x30_medium.stl",
    "tests/onshape/cylinder_10x30.step"
);
test_stl_step_compare!(
    onshape_dome_hemisphere_compare,
    "tests/onshape/dome_hemisphere_20_fine.stl",
    "tests/onshape/dome_hemisphere_20.step"
);
test_stl_step_compare!(
    onshape_l_bracket_compare,
    "tests/onshape/l_bracket_simple_medium.stl",
    "tests/onshape/l_bracket_simple.step"
);
test_stl_step_compare!(
    onshape_pipe_elbow_compare,
    "tests/onshape/pipe_elbow_10_fine.stl",
    "tests/onshape/pipe_elbow_10.step"
);
test_stl_step_compare!(
    onshape_plate_with_hole_compare,
    "tests/onshape/plate_with_hole_100x50_coarse.stl",
    "tests/onshape/plate_with_hole_100x50.step"
);
test_stl_step_compare!(
    onshape_rounded_cube_fine_compare,
    "tests/onshape/rounded_cube_10_r2_fine.stl",
    "tests/onshape/rounded_cube_10_r2.step"
);
test_stl_step_compare!(
    onshape_rounded_cube_medium_compare,
    "tests/onshape/rounded_cube_10_r2_medium.stl",
    "tests/onshape/rounded_cube_10_r2.step"
);
test_stl_step_compare!(
    onshape_rounded_cube_coarse_compare,
    "tests/onshape/rounded_cube_10_r2_coarse.stl",
    "tests/onshape/rounded_cube_10_r2.step"
);
test_stl_step_compare!(
    onshape_sphere_compare,
    "tests/onshape/sphere_25_fine.stl",
    "tests/onshape/sphere_25.step"
);
test_stl_step_compare!(
    onshape_stepped_block_compare,
    "tests/onshape/stepped_block_coarse.stl",
    "tests/onshape/stepped_block.step"
);

// --- tests/occt/ ---
test_stl_step_compare!(
    occt_hollow_cube_compare,
    "tests/occt/hollow_cube.stl",
    "tests/occt/hollow_cube.step"
);

#[test]
fn occt_hollow_cube_void_stats() {
    // The hollow cube has 2 closed shells; the inner one (negative signed
    // volume, contained in the outer shell's bounding box) is a legitimate
    // internal void, not an inverted-winding error.
    let stl = format!("{}/tests/occt/hollow_cube.stl", manifest_dir());
    let config = config_for_stl(&stl);
    let mesh = stage1::stage1(&config).expect("hollow cube should validate");
    assert_eq!(mesh.stats.connected_shells, 2, "should have 2 shells");
    assert_eq!(mesh.stats.voids_within_solids, 1, "should have 1 void");
}

// --- tests/fusion/ ---
test_stl_step_compare!(
    fusion_plate_high_compare,
    "tests/fusion/plate_with_hole_100x50_high.stl",
    "tests/fusion/plate_with_hole_100x50.step"
);
test_stl_step_compare!(
    fusion_plate_low_compare,
    "tests/fusion/plate_with_hole_100x50_low.stl",
    "tests/fusion/plate_with_hole_100x50.step"
);
test_stl_step_compare!(
    fusion_plate_medium_compare,
    "tests/fusion/plate_with_hole_100x50_medium.stl",
    "tests/fusion/plate_with_hole_100x50.step"
);

// --- tests/nurbs/ (mesh vertices lie on the exact reference surfaces) ---
test_stl_step_compare!(
    nurbs_elliptical_arch_compare,
    "tests/nurbs/elliptical_arch.stl",
    "tests/nurbs/elliptical_arch.step"
);
test_stl_step_compare!(
    nurbs_elliptical_cylinder_compare,
    "tests/nurbs/elliptical_cylinder.stl",
    "tests/nurbs/elliptical_cylinder.step"
);
test_stl_step_compare!(
    nurbs_bspline_hill_compare,
    "tests/nurbs/bspline_hill.stl",
    "tests/nurbs/bspline_hill.step"
);
test_stl_step_compare!(
    nurbs_bspline_saddle_compare,
    "tests/nurbs/bspline_saddle.stl",
    "tests/nurbs/bspline_saddle.step"
);
test_stl_step_compare!(
    nurbs_bspline_terrain_compare,
    "tests/nurbs/bspline_terrain.stl",
    "tests/nurbs/bspline_terrain.step"
);
test_stl_step_compare!(
    nurbs_bspline_blister_compare,
    "tests/nurbs/bspline_blister.stl",
    "tests/nurbs/bspline_blister.step"
);
test_stl_step_compare!(
    nurbs_bspline_wave_compare,
    "tests/nurbs/bspline_wave.stl",
    "tests/nurbs/bspline_wave.step"
);
test_stl_step_compare!(
    nurbs_hill_with_hole_compare,
    "tests/nurbs/hill_with_hole.stl",
    "tests/nurbs/hill_with_hole.step"
);

// ===========================================================================
// Bad cases: mesh validation failures
// ===========================================================================

#[test]
fn bad_degenerate_face_strict() {
    let stl = format!("{}/tests/bad/degenerate_face.stl", manifest_dir());
    let config =
        config::parse_config_from(["brepper", &stl, "--stage=1", "-q", "--strict-mesh"]).unwrap();
    let err = stage1::stage1(&config).expect_err("should fail with degenerate face");
    let msg = err.to_string();
    assert!(
        msg.contains("degenerate"),
        "error should mention degenerate faces, got: {msg}"
    );
}

#[test]
fn bad_degenerate_face_repaired_by_default() {
    let stl = format!("{}/tests/bad/degenerate_face.stl", manifest_dir());
    let config = config_for_stl(&stl);
    let mesh = stage1::stage1(&config).expect("degenerate face should be repaired by default");
    assert_eq!(mesh.stats.mesh_faces_degenerate, 0);
    // The sliver is removed and its neighbor split; the fixture's genuinely
    // missing bottom triangle remains as an open hole.
    assert_eq!(mesh.stats.mesh_edges_open, 4);
}

#[test]
fn bad_non_manifold_edge() {
    let stl = format!("{}/tests/bad/non_manifold_edge.stl", manifest_dir());
    let config = config_for_stl(&stl);
    let err = stage1::stage1(&config).expect_err("should fail with non-manifold edge");
    let msg = err.to_string();
    assert!(
        msg.contains("non-manifold"),
        "error should mention non-manifold edges, got: {msg}"
    );
}

#[test]
fn bad_inconsistent_winding_repaired_by_default() {
    let stl = format!("{}/tests/bad/inconsistent_winding.stl", manifest_dir());
    let config = config_for_stl(&stl);
    let mesh = stage1::stage1(&config).expect("flipped face should be repaired by default");
    assert_eq!(mesh.stats.mesh_edges_inconsistent_orientation, 0);
    assert_eq!(mesh.stats.solids, 1);
}

#[test]
fn bad_inconsistent_winding_strict() {
    let stl = format!("{}/tests/bad/inconsistent_winding.stl", manifest_dir());
    let config =
        config::parse_config_from(["brepper", &stl, "--stage=1", "-q", "--strict-mesh"]).unwrap();
    let err = stage1::stage1(&config).expect_err("should fail with inconsistent orientation");
    let msg = err.to_string();
    assert!(
        msg.contains("flipped orientation"),
        "error should mention flipped orientation, got: {msg}"
    );
}

#[test]
fn float32_position_bound_scales_with_distance_from_origin() {
    // The same body, tessellated by the same tessellator at the same settings,
    // once at the origin and once translated ~1 m (tests/remesh). An STL stores
    // f32, so the far copy arrives with coarser vertex precision — the point of
    // the off-origin variants, and the reason the vertex tolerance needs a
    // floor (Phase 0.6).
    let dir = manifest_dir();
    let near = stage1::stage1(&config_for_stl(&format!(
        "{dir}/tests/remesh/pipe_elbow_10_occt_medium.stl"
    )))
    .expect("origin-centred stage1 should pass");
    let far = stage1::stage1(&config_for_stl(&format!(
        "{dir}/tests/remesh/pipe_elbow_10_offset_occt_medium.stl"
    )))
    .expect("off-origin stage1 should pass");

    let near_bound = stage1::float32_position_bound_mm(&near);
    let far_bound = stage1::float32_position_bound_mm(&far);
    // Measured 2026-07-29: 4.5e-6 mm at the origin, 6.0e-5 mm at (600, 800, 0).
    assert!(
        near_bound > 1.0e-6 && near_bound < 1.0e-5,
        "origin-centred position bound {near_bound:.3e} mm outside the expected 1e-6..1e-5"
    );
    assert!(
        far_bound > 10.0 * near_bound,
        "translating 1 m should coarsen precision by more than 10x: {near_bound:.3e} -> {far_bound:.3e}"
    );

    // The clamp is INERT ON THE WHOLE CORPUS, including at 1 m: the bound
    // there is 4.3e-5 mm, still under the 5e-5 mm default vertex tolerance.
    // That is the measurement, not a target — it says the default has ~14%
    // margin at a metre out, and none beyond about 1.02 m, where the f32
    // exponent steps and the bound roughly doubles.
    let cfg = config_for_stl(&format!("{dir}/tests/remesh/pipe_elbow_10_offset_occt_medium.stl"));
    for (label, mesh) in [("origin-centred", &near), ("1 m off-origin", &far)] {
        assert_eq!(
            stage1::effective_vertex_tolerance_mm(&cfg, mesh),
            cfg.vertex_tolerance_mm,
            "clamp should still be inert on the {label} corpus mesh"
        );
    }

    // Where it does bite: a part positioned several metres out, as in a large
    // assembly. Synthetic vertices rather than a corpus file, because the point
    // is the arithmetic of f32 storage, not any particular model.
    let far_5m = brepper::stage1::ConnectedMesh {
        vertices: vec![
            brepper::stage1::MeshVertex::from_xyz(5000.0, 4000.0, 3000.0),
            brepper::stage1::MeshVertex::from_xyz(5010.0, 4010.0, 3010.0),
        ],
        ..Default::default()
    };
    let far_bound_5m = stage1::float32_position_bound_mm(&far_5m);
    assert!(
        far_bound_5m > cfg.vertex_tolerance_mm,
        "at 5 m the f32 position bound {far_bound_5m:.3e} mm should exceed the {:.3e} mm default",
        cfg.vertex_tolerance_mm
    );
    assert_eq!(
        stage1::effective_vertex_tolerance_mm(&cfg, &far_5m),
        far_bound_5m,
        "the clamp must raise the tolerance to the bound once the bound is the larger"
    );
}

#[test]
fn chamfered_cube_ascii_matches_binary() {
    // tests/onshape/chamfered_cube_10_c1_medium.ascii.stl is the ASCII export of
    // the same body as the binary variant, and the Phase 0.1 coverage audit
    // found it referenced by no test at all. Its purpose is the ASCII reader, so
    // pin exactly that: identical topology and the same geometry, rather than a
    // duplicate run of the whole pipeline.
    //
    // The bar is 1e-4 mm because that is the ASCII file's own resolution — it
    // carries 6 significant digits ("48.0012"), coarser than the binary file's
    // float32 (~6.7e-6 mm at these coordinates). It is a file-format precision,
    // not a tolerance to tune.
    let dir = manifest_dir();
    let ascii = stage1::stage1(&config_for_stl(&format!(
        "{dir}/tests/onshape/chamfered_cube_10_c1_medium.ascii.stl"
    )))
    .expect("ASCII stage1 should pass");
    let binary = stage1::stage1(&config_for_stl(&format!(
        "{dir}/tests/onshape/chamfered_cube_10_c1_medium.stl"
    )))
    .expect("binary stage1 should pass");

    assert_eq!(ascii.vertices.len(), binary.vertices.len(), "vertex count");
    assert_eq!(ascii.faces.len(), binary.faces.len(), "face count");
    assert_eq!(ascii.stats.mesh_quads, binary.stats.mesh_quads, "quad count");
    assert_eq!(ascii.stats.solids, binary.stats.solids, "solid count");
    assert_eq!(ascii.stats.mesh_edges_open, binary.stats.mesh_edges_open, "open edges");

    // Order-independent: every ASCII vertex has a binary counterpart nearby.
    for (i, a) in ascii.vertices.iter().enumerate() {
        let best = binary
            .vertices
            .iter()
            .map(|b| ((a.x - b.x).powi(2) + (a.y - b.y).powi(2) + (a.z - b.z).powi(2)).sqrt())
            .fold(f64::INFINITY, f64::min);
        assert!(
            best < 1.0e-4,
            "ASCII vertex {i} at ({}, {}, {}) is {best:.3e} mm from the nearest binary vertex",
            a.x,
            a.y,
            a.z
        );
    }
}

// ===========================================================================
// Bad cases: --compare failures
// ===========================================================================

#[test]
fn bad_compare_cube_shifted() {
    let stl = format!("{}/tests/bad/cube_shifted.stl", manifest_dir());
    let step = format!("{}/tests/bad/cube_shifted.step", manifest_dir());
    let config = config_for_compare(&stl, &step);
    let err = stage1::stage1(&config).expect_err("should fail with shifted vertex");
    let msg = err.to_string();
    assert!(
        msg.contains("comparison failed"),
        "error should mention comparison failure, got: {msg}"
    );
}

#[test]
fn bad_compare_cube_on_plane() {
    // Vertex shifted along an infinite surface plane (y=1, z=1) but far off
    // the bounded face. Verifies that --compare uses bounded face distance.
    let stl = format!("{}/tests/bad/cube_on_plane.stl", manifest_dir());
    let step = format!("{}/tests/bad/cube_on_plane.step", manifest_dir());
    let config = config_for_compare(&stl, &step);
    let err = stage1::stage1(&config).expect_err("should fail: vertex on infinite plane but off bounded face");
    let msg = err.to_string();
    assert!(
        msg.contains("comparison failed"),
        "error should mention comparison failure, got: {msg}"
    );
}
