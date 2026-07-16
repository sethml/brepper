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
fn bad_inconsistent_winding() {
    let stl = format!("{}/tests/bad/inconsistent_winding.stl", manifest_dir());
    let config = config_for_stl(&stl);
    let err = stage1::stage1(&config).expect_err("should fail with inconsistent orientation");
    let msg = err.to_string();
    assert!(
        msg.contains("flipped orientation"),
        "error should mention flipped orientation, got: {msg}"
    );
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
