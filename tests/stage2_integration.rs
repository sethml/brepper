//! Integration tests for Stage 2: surface fitting.
//!
//! Tests that planar hypothesis deduction produces correct results for
//! known planar models, and that --compare validation passes.

use brepper::config;
use brepper::{stage1, stage2};

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

fn manifest_dir() -> String {
    env!("CARGO_MANIFEST_DIR").to_string()
}

fn config_for_stl(stl_path: &str) -> config::Config {
    config::parse_config_from(["brepper", stl_path, "--stage=2.1", "-q"]).unwrap()
}

fn config_for_compare(stl_path: &str, step_path: &str) -> config::Config {
    let mut config =
        config::parse_config_from(["brepper", stl_path, "--compare", step_path, "--stage=2.1", "-q"])
            .unwrap();
    config.load_compare_step().unwrap();
    config
}

/// Run stages 1 and 2.1, returning the Stage2Output.
fn run_stage2(config: &config::Config) -> stage2::Stage2Output {
    let mesh = stage1::stage1(config).expect("stage1 should pass");
    stage2::stage2(config, mesh).expect("stage2 should pass")
}

// ---------------------------------------------------------------------------
// Macros for generating test cases
// ---------------------------------------------------------------------------

/// Generate a test that runs stage 2.1 with --compare against a STEP file.
macro_rules! test_stl_step_stage2_compare {
    ($name:ident, $stl_path:literal, $step_path:literal) => {
        #[test]
        fn $name() {
            let stl = format!("{}/{}", manifest_dir(), $stl_path);
            let step = format!("{}/{}", manifest_dir(), $step_path);
            let config = config_for_compare(&stl, &step);
            run_stage2(&config);
        }
    };
}

// ===========================================================================
// Planar hypothesis count tests
// ===========================================================================

#[test]
fn cube_produces_six_planar_hypotheses() {
    let stl = format!("{}/tests/manual/cube.stl", manifest_dir());
    let config = config_for_stl(&stl);
    let output = run_stage2(&config);

    assert_eq!(
        output.planar_hypotheses.len(),
        6,
        "cube should have 6 planar hypotheses"
    );
    for h in &output.planar_hypotheses {
        assert_eq!(h.faces.len(), 2, "each cube face hypothesis should have 2 triangles");
    }
}

#[test]
fn wedge_produces_six_planar_hypotheses() {
    let stl = format!("{}/tests/ccad/generated/wedge.stl", manifest_dir());
    let config = config_for_stl(&stl);
    let output = run_stage2(&config);

    assert_eq!(
        output.planar_hypotheses.len(),
        6,
        "wedge should have 6 planar hypotheses"
    );
}

#[test]
fn chamfered_cube_produces_26_planar_hypotheses() {
    let stl = format!(
        "{}/tests/onshape/chamfered_cube_10_c1_medium.stl",
        manifest_dir()
    );
    let config = config_for_stl(&stl);
    let output = run_stage2(&config);

    // 6 main faces + 12 edge chamfers + 8 corner chamfers = 26
    assert_eq!(
        output.planar_hypotheses.len(),
        26,
        "chamfered cube should have 26 planar hypotheses"
    );
}

#[test]
fn all_faces_assigned_to_planar_hypotheses() {
    // For a fully planar model, every face should belong to some hypothesis
    let stl = format!("{}/tests/manual/cube.stl", manifest_dir());
    let config = config_for_stl(&stl);
    let output = run_stage2(&config);

    let total_faces: usize = output.planar_hypotheses.iter().map(|h| h.faces.len()).sum();
    assert_eq!(
        total_faces,
        output.mesh.faces.len(),
        "all faces should be assigned"
    );
}

// ===========================================================================
// Compare tests: all test STL/STEP pairs should pass --compare at stage 2.1
// ===========================================================================

// --- tests/manual/ ---
test_stl_step_stage2_compare!(
    manual_cube_stage2_compare,
    "tests/manual/cube.stl",
    "tests/manual/cube.step"
);

// --- tests/ccad/generated/ ---
test_stl_step_stage2_compare!(
    ccad_cube_stage2_compare,
    "tests/ccad/generated/cube.stl",
    "tests/ccad/generated/cube.step"
);
test_stl_step_stage2_compare!(
    ccad_wedge_stage2_compare,
    "tests/ccad/generated/wedge.stl",
    "tests/ccad/generated/wedge.step"
);
test_stl_step_stage2_compare!(
    ccad_t_shape_stage2_compare,
    "tests/ccad/generated/t_shape.stl",
    "tests/ccad/generated/t_shape.step"
);
test_stl_step_stage2_compare!(
    ccad_staircase_stage2_compare,
    "tests/ccad/generated/staircase.stl",
    "tests/ccad/generated/staircase.step"
);

// --- tests/onshape/ (planar models) ---
test_stl_step_stage2_compare!(
    onshape_chamfered_cube_stage2_compare,
    "tests/onshape/chamfered_cube_10_c1_medium.stl",
    "tests/onshape/chamfered_cube_10_c1_medium.step"
);
test_stl_step_stage2_compare!(
    onshape_stepped_block_stage2_compare,
    "tests/onshape/stepped_block_coarse.stl",
    "tests/onshape/stepped_block_coarse.step"
);
test_stl_step_stage2_compare!(
    onshape_l_bracket_stage2_compare,
    "tests/onshape/l_bracket_simple_medium.stl",
    "tests/onshape/l_bracket_simple_medium.step"
);

// --- tests/onshape/ (mixed planar/curved models — should still pass) ---
test_stl_step_stage2_compare!(
    onshape_cone_stage2_compare,
    "tests/onshape/cone_15x20_medium.stl",
    "tests/onshape/cone_15x20_medium.step"
);
test_stl_step_stage2_compare!(
    onshape_cylinder_stage2_compare,
    "tests/onshape/cylinder_10x30_medium.stl",
    "tests/onshape/cylinder_10x30_medium.step"
);
test_stl_step_stage2_compare!(
    onshape_sphere_stage2_compare,
    "tests/onshape/sphere_25_fine.stl",
    "tests/onshape/sphere_25_fine.step"
);
test_stl_step_stage2_compare!(
    onshape_dome_hemisphere_stage2_compare,
    "tests/onshape/dome_hemisphere_20_fine.stl",
    "tests/onshape/dome_hemisphere_20_fine.step"
);
test_stl_step_stage2_compare!(
    onshape_pipe_elbow_stage2_compare,
    "tests/onshape/pipe_elbow_10_fine.stl",
    "tests/onshape/pipe_elbow_10_fine.step"
);
test_stl_step_stage2_compare!(
    onshape_plate_with_hole_stage2_compare,
    "tests/onshape/plate_with_hole_100x50_coarse.stl",
    "tests/onshape/plate_with_hole_100x50_coarse.step"
);
test_stl_step_stage2_compare!(
    onshape_rounded_cube_stage2_compare,
    "tests/onshape/rounded_cube_10_r2_fine.stl",
    "tests/onshape/rounded_cube_10_r2_fine.step"
);

// --- tests/fusion/ ---
test_stl_step_stage2_compare!(
    fusion_plate_high_stage2_compare,
    "tests/fusion/plate_with_hole_100x50_high.stl",
    "tests/fusion/plate_with_hole_100x50_high.step"
);
test_stl_step_stage2_compare!(
    fusion_plate_low_stage2_compare,
    "tests/fusion/plate_with_hole_100x50_low.stl",
    "tests/fusion/plate_with_hole_100x50_low.step"
);
test_stl_step_stage2_compare!(
    fusion_plate_medium_stage2_compare,
    "tests/fusion/plate_with_hole_100x50_medium.stl",
    "tests/fusion/plate_with_hole_100x50_medium.step"
);
