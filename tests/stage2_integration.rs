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
fn config_for_stl_stage22(stl_path: &str) -> config::Config {
    config::parse_config_from(["brepper", stl_path, "--stage=2.2", "-q"]).unwrap()
}

fn config_for_compare_stage22(stl_path: &str, step_path: &str) -> config::Config {
    let mut config =
        config::parse_config_from(["brepper", stl_path, "--compare", step_path, "--stage=2.2", "-q"])
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

/// Generate a test that runs stage 2.2 with --compare against a STEP file.
macro_rules! test_stl_step_stage22_compare {
    ($name:ident, $stl_path:literal, $step_path:literal) => {
        #[test]
        fn $name() {
            let stl = format!("{}/{}", manifest_dir(), $stl_path);
            let step = format!("{}/{}", manifest_dir(), $step_path);
            let config = config_for_compare_stage22(&stl, &step);
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
        assert_eq!(h.faces.len(), 1, "each cube face hypothesis should have 1 quad");
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
test_stl_step_stage2_compare!(
    ccad_simple_cylinder_stage2_compare,
    "tests/ccad/generated/simple_cylinder.stl",
    "tests/ccad/generated/simple_cylinder.step"
);
test_stl_step_stage2_compare!(
    ccad_block_with_hole_stage2_compare,
    "tests/ccad/generated/block_with_hole.stl",
    "tests/ccad/generated/block_with_hole.step"
);
test_stl_step_stage2_compare!(
    ccad_pipe_stage2_compare,
    "tests/ccad/generated/pipe.stl",
    "tests/ccad/generated/pipe.step"
);
test_stl_step_stage2_compare!(
    ccad_stepped_cylinder_stage2_compare,
    "tests/ccad/generated/stepped_cylinder.stl",
    "tests/ccad/generated/stepped_cylinder.step"
);
test_stl_step_stage2_compare!(
    ccad_two_holes_stage2_compare,
    "tests/ccad/generated/two_holes.stl",
    "tests/ccad/generated/two_holes.step"
);
test_stl_step_stage2_compare!(
    ccad_simple_sphere_stage2_compare,
    "tests/ccad/generated/simple_sphere.stl",
    "tests/ccad/generated/simple_sphere.step"
);
test_stl_step_stage2_compare!(
    ccad_hemisphere_stage2_compare,
    "tests/ccad/generated/hemisphere.stl",
    "tests/ccad/generated/hemisphere.step"
);
test_stl_step_stage2_compare!(
    ccad_spherical_pocket_stage2_compare,
    "tests/ccad/generated/spherical_pocket.stl",
    "tests/ccad/generated/spherical_pocket.step"
);
test_stl_step_stage2_compare!(
    ccad_ball_on_cylinder_stage2_compare,
    "tests/ccad/generated/ball_on_cylinder.stl",
    "tests/ccad/generated/ball_on_cylinder.step"
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


// ===========================================================================
// Stage 2.2: Cylindrical hypothesis count tests
// ===========================================================================

#[test]
fn simple_cylinder_produces_one_cylindrical_hypothesis() {
    let stl = format!("{}/tests/ccad/generated/simple_cylinder.stl", manifest_dir());
    let config = config_for_stl_stage22(&stl);
    let output = run_stage2(&config);

    assert_eq!(
        output.cylindrical_hypotheses.len(),
        1,
        "simple_cylinder should have 1 cylindrical hypothesis"
    );
    let h = &output.cylindrical_hypotheses[0];
    assert!(h.convex, "simple_cylinder should be convex");
    assert!((h.radius - 10.0).abs() < 0.01, "radius should be ~10.0, got {}", h.radius);
}

#[test]
fn block_with_hole_produces_one_concave_cylinder() {
    let stl = format!("{}/tests/ccad/generated/block_with_hole.stl", manifest_dir());
    let config = config_for_stl_stage22(&stl);
    let output = run_stage2(&config);

    assert_eq!(
        output.cylindrical_hypotheses.len(),
        1,
        "block_with_hole should have 1 cylindrical hypothesis"
    );
    let h = &output.cylindrical_hypotheses[0];
    assert!(!h.convex, "block_with_hole cylinder should be concave");
    assert!((h.radius - 6.0).abs() < 0.01, "radius should be ~6.0, got {}", h.radius);
}

#[test]
fn pipe_produces_two_cylindrical_hypotheses() {
    let stl = format!("{}/tests/ccad/generated/pipe.stl", manifest_dir());
    let config = config_for_stl_stage22(&stl);
    let output = run_stage2(&config);

    assert_eq!(
        output.cylindrical_hypotheses.len(),
        2,
        "pipe should have 2 cylindrical hypotheses"
    );
    // One convex (outer) and one concave (inner)
    let convex_count = output.cylindrical_hypotheses.iter().filter(|h| h.convex).count();
    let concave_count = output.cylindrical_hypotheses.iter().filter(|h| !h.convex).count();
    assert_eq!(convex_count, 1, "pipe should have 1 convex cylinder");
    assert_eq!(concave_count, 1, "pipe should have 1 concave cylinder");
}

#[test]
fn two_holes_produces_two_concave_cylinders() {
    let stl = format!("{}/tests/ccad/generated/two_holes.stl", manifest_dir());
    let config = config_for_stl_stage22(&stl);
    let output = run_stage2(&config);

    assert_eq!(
        output.cylindrical_hypotheses.len(),
        2,
        "two_holes should have 2 cylindrical hypotheses"
    );
    let concave_count = output.cylindrical_hypotheses.iter().filter(|h| !h.convex).count();
    assert_eq!(concave_count, 2, "both cylinders should be concave");
    let mut radii: Vec<f64> = output.cylindrical_hypotheses.iter().map(|h| h.radius).collect();
    radii.sort_by(|a, b| a.partial_cmp(b).unwrap());
    assert!((radii[0] - 3.0).abs() < 0.01, "smaller radius should be ~3.0, got {}", radii[0]);
    assert!((radii[1] - 5.0).abs() < 0.01, "larger radius should be ~5.0, got {}", radii[1]);
}

#[test]
fn chamfered_cube_produces_no_cylindrical_hypotheses() {
    let stl = format!(
        "{}/tests/onshape/chamfered_cube_10_c1_medium.stl",
        manifest_dir()
    );
    let config = config_for_stl_stage22(&stl);
    let output = run_stage2(&config);

    assert_eq!(
        output.cylindrical_hypotheses.len(),
        0,
        "chamfered cube should have no cylindrical hypotheses"
    );
}

// ===========================================================================
// Stage 2.2: Compare tests (all models should pass --compare at stage 2.2)
// ===========================================================================

// --- tests/ccad/generated/ ---
test_stl_step_stage22_compare!(
    ccad_cube_stage22_compare,
    "tests/ccad/generated/cube.stl",
    "tests/ccad/generated/cube.step"
);
test_stl_step_stage22_compare!(
    ccad_simple_cylinder_stage22_compare,
    "tests/ccad/generated/simple_cylinder.stl",
    "tests/ccad/generated/simple_cylinder.step"
);
test_stl_step_stage22_compare!(
    ccad_block_with_hole_stage22_compare,
    "tests/ccad/generated/block_with_hole.stl",
    "tests/ccad/generated/block_with_hole.step"
);
test_stl_step_stage22_compare!(
    ccad_pipe_stage22_compare,
    "tests/ccad/generated/pipe.stl",
    "tests/ccad/generated/pipe.step"
);
test_stl_step_stage22_compare!(
    ccad_stepped_cylinder_stage22_compare,
    "tests/ccad/generated/stepped_cylinder.stl",
    "tests/ccad/generated/stepped_cylinder.step"
);
test_stl_step_stage22_compare!(
    ccad_two_holes_stage22_compare,
    "tests/ccad/generated/two_holes.stl",
    "tests/ccad/generated/two_holes.step"
);
test_stl_step_stage22_compare!(
    ccad_simple_sphere_stage22_compare,
    "tests/ccad/generated/simple_sphere.stl",
    "tests/ccad/generated/simple_sphere.step"
);
test_stl_step_stage22_compare!(
    ccad_hemisphere_stage22_compare,
    "tests/ccad/generated/hemisphere.stl",
    "tests/ccad/generated/hemisphere.step"
);
test_stl_step_stage22_compare!(
    ccad_spherical_pocket_stage22_compare,
    "tests/ccad/generated/spherical_pocket.stl",
    "tests/ccad/generated/spherical_pocket.step"
);
test_stl_step_stage22_compare!(
    ccad_ball_on_cylinder_stage22_compare,
    "tests/ccad/generated/ball_on_cylinder.stl",
    "tests/ccad/generated/ball_on_cylinder.step"
);

// --- tests/onshape/ ---
test_stl_step_stage22_compare!(
    onshape_cylinder_stage22_compare,
    "tests/onshape/cylinder_10x30_medium.stl",
    "tests/onshape/cylinder_10x30_medium.step"
);
test_stl_step_stage22_compare!(
    onshape_l_bracket_stage22_compare,
    "tests/onshape/l_bracket_simple_medium.stl",
    "tests/onshape/l_bracket_simple_medium.step"
);
test_stl_step_stage22_compare!(
    onshape_plate_with_hole_stage22_compare,
    "tests/onshape/plate_with_hole_100x50_coarse.stl",
    "tests/onshape/plate_with_hole_100x50_coarse.step"
);
test_stl_step_stage22_compare!(
    onshape_pipe_elbow_stage22_compare,
    "tests/onshape/pipe_elbow_10_fine.stl",
    "tests/onshape/pipe_elbow_10_fine.step"
);

// --- tests/fusion/ ---
test_stl_step_stage22_compare!(
    fusion_plate_low_stage22_compare,
    "tests/fusion/plate_with_hole_100x50_low.stl",
    "tests/fusion/plate_with_hole_100x50_low.step"
);
test_stl_step_stage22_compare!(
    fusion_plate_medium_stage22_compare,
    "tests/fusion/plate_with_hole_100x50_medium.stl",
    "tests/fusion/plate_with_hole_100x50_medium.step"
);
test_stl_step_stage22_compare!(
    fusion_plate_high_stage22_compare,
    "tests/fusion/plate_with_hole_100x50_high.stl",
    "tests/fusion/plate_with_hole_100x50_high.step"
);

// ===========================================================================
// Stage 2.3: Spherical hypothesis count tests
// ===========================================================================

fn config_for_stl_stage23(stl_path: &str) -> config::Config {
    config::parse_config_from(["brepper", stl_path, "--stage=2.3", "-q"]).unwrap()
}

fn config_for_compare_stage23(stl_path: &str, step_path: &str) -> config::Config {
    let mut config =
        config::parse_config_from(["brepper", stl_path, "--compare", step_path, "--stage=2.3", "-q"])
            .unwrap();
    config.load_compare_step().unwrap();
    config
}

/// Generate a test that runs stage 2.3 with --compare against a STEP file.
macro_rules! test_stl_step_stage23_compare {
    ($name:ident, $stl_path:literal, $step_path:literal) => {
        #[test]
        fn $name() {
            let stl = format!("{}/{}", manifest_dir(), $stl_path);
            let step = format!("{}/{}", manifest_dir(), $step_path);
            let config = config_for_compare_stage23(&stl, &step);
            run_stage2(&config);
        }
    };
}

#[test]
fn simple_sphere_produces_one_spherical_hypothesis() {
    let stl = format!("{}/tests/ccad/generated/simple_sphere.stl", manifest_dir());
    let config = config_for_stl_stage23(&stl);
    let output = run_stage2(&config);

    assert_eq!(
        output.spherical_hypotheses.len(),
        1,
        "simple_sphere should have 1 spherical hypothesis"
    );
    let h = &output.spherical_hypotheses[0];
    assert!(h.convex, "simple_sphere should be convex");
    assert!((h.radius - 10.0).abs() < 0.01, "radius should be ~10.0, got {}", h.radius);
}

#[test]
fn hemisphere_produces_one_spherical_hypothesis() {
    let stl = format!("{}/tests/ccad/generated/hemisphere.stl", manifest_dir());
    let config = config_for_stl_stage23(&stl);
    let output = run_stage2(&config);

    assert_eq!(
        output.spherical_hypotheses.len(),
        1,
        "hemisphere should have 1 spherical hypothesis"
    );
    let h = &output.spherical_hypotheses[0];
    assert!(h.convex, "hemisphere should be convex");
    assert!((h.radius - 10.0).abs() < 0.01, "radius should be ~10.0, got {}", h.radius);
}

#[test]
fn spherical_pocket_produces_one_concave_sphere() {
    let stl = format!("{}/tests/ccad/generated/spherical_pocket.stl", manifest_dir());
    let config = config_for_stl_stage23(&stl);
    let output = run_stage2(&config);

    assert_eq!(
        output.spherical_hypotheses.len(),
        1,
        "spherical_pocket should have 1 spherical hypothesis"
    );
    let h = &output.spherical_hypotheses[0];
    assert!(!h.convex, "spherical_pocket sphere should be concave");
    assert!((h.radius - 8.0).abs() < 0.01, "radius should be ~8.0, got {}", h.radius);
}

#[test]
fn ball_on_cylinder_produces_one_sphere_one_cylinder() {
    let stl = format!("{}/tests/ccad/generated/ball_on_cylinder.stl", manifest_dir());
    let config = config_for_stl_stage23(&stl);
    let output = run_stage2(&config);

    assert_eq!(
        output.cylindrical_hypotheses.len(),
        1,
        "ball_on_cylinder should have 1 cylindrical hypothesis"
    );
    assert_eq!(
        output.spherical_hypotheses.len(),
        1,
        "ball_on_cylinder should have 1 spherical hypothesis"
    );
    let cyl = &output.cylindrical_hypotheses[0];
    assert!(cyl.convex, "cylinder should be convex");
    assert!((cyl.radius - 5.0).abs() < 0.01, "cylinder radius should be ~5.0, got {}", cyl.radius);
    let sph = &output.spherical_hypotheses[0];
    assert!(sph.convex, "sphere should be convex");
    assert!((sph.radius - 8.0).abs() < 0.01, "sphere radius should be ~8.0, got {}", sph.radius);
}

#[test]
fn chamfered_cube_produces_no_spherical_hypotheses() {
    let stl = format!(
        "{}/tests/onshape/chamfered_cube_10_c1_medium.stl",
        manifest_dir()
    );
    let config = config_for_stl_stage23(&stl);
    let output = run_stage2(&config);

    assert_eq!(
        output.spherical_hypotheses.len(),
        0,
        "chamfered cube should have no spherical hypotheses"
    );
}

#[test]
fn onshape_sphere_produces_one_spherical_hypothesis() {
    let stl = format!("{}/tests/onshape/sphere_25_fine.stl", manifest_dir());
    let config = config_for_stl_stage23(&stl);
    let output = run_stage2(&config);

    assert_eq!(
        output.spherical_hypotheses.len(),
        1,
        "onshape sphere should have 1 spherical hypothesis"
    );
    let h = &output.spherical_hypotheses[0];
    assert!(h.convex, "onshape sphere should be convex");
    assert!((h.radius - 12.5).abs() < 0.1, "radius should be ~12.5, got {}", h.radius);
}

#[test]
fn onshape_dome_hemisphere_produces_one_spherical_hypothesis() {
    let stl = format!("{}/tests/onshape/dome_hemisphere_20_fine.stl", manifest_dir());
    let config = config_for_stl_stage23(&stl);
    let output = run_stage2(&config);

    assert_eq!(
        output.spherical_hypotheses.len(),
        1,
        "onshape dome should have 1 spherical hypothesis"
    );
    let h = &output.spherical_hypotheses[0];
    assert!(h.convex, "onshape dome should be convex");
    assert!((h.radius - 10.0).abs() < 0.1, "radius should be ~10.0, got {}", h.radius);
}


// ===========================================================================
// Stage 2.3: Compare tests (all models should pass --compare at stage 2.3)
// ===========================================================================

// --- tests/ccad/generated/ ---
test_stl_step_stage23_compare!(
    ccad_cube_stage23_compare,
    "tests/ccad/generated/cube.stl",
    "tests/ccad/generated/cube.step"
);
test_stl_step_stage23_compare!(
    ccad_simple_cylinder_stage23_compare,
    "tests/ccad/generated/simple_cylinder.stl",
    "tests/ccad/generated/simple_cylinder.step"
);
test_stl_step_stage23_compare!(
    ccad_block_with_hole_stage23_compare,
    "tests/ccad/generated/block_with_hole.stl",
    "tests/ccad/generated/block_with_hole.step"
);
test_stl_step_stage23_compare!(
    ccad_pipe_stage23_compare,
    "tests/ccad/generated/pipe.stl",
    "tests/ccad/generated/pipe.step"
);
test_stl_step_stage23_compare!(
    ccad_stepped_cylinder_stage23_compare,
    "tests/ccad/generated/stepped_cylinder.stl",
    "tests/ccad/generated/stepped_cylinder.step"
);
test_stl_step_stage23_compare!(
    ccad_two_holes_stage23_compare,
    "tests/ccad/generated/two_holes.stl",
    "tests/ccad/generated/two_holes.step"
);
test_stl_step_stage23_compare!(
    ccad_simple_sphere_stage23_compare,
    "tests/ccad/generated/simple_sphere.stl",
    "tests/ccad/generated/simple_sphere.step"
);
test_stl_step_stage23_compare!(
    ccad_hemisphere_stage23_compare,
    "tests/ccad/generated/hemisphere.stl",
    "tests/ccad/generated/hemisphere.step"
);
test_stl_step_stage23_compare!(
    ccad_spherical_pocket_stage23_compare,
    "tests/ccad/generated/spherical_pocket.stl",
    "tests/ccad/generated/spherical_pocket.step"
);
test_stl_step_stage23_compare!(
    ccad_ball_on_cylinder_stage23_compare,
    "tests/ccad/generated/ball_on_cylinder.stl",
    "tests/ccad/generated/ball_on_cylinder.step"
);

// --- tests/onshape/ ---
test_stl_step_stage23_compare!(
    onshape_sphere_stage23_compare,
    "tests/onshape/sphere_25_fine.stl",
    "tests/onshape/sphere_25_fine.step"
);
test_stl_step_stage23_compare!(
    onshape_dome_hemisphere_stage23_compare,
    "tests/onshape/dome_hemisphere_20_fine.stl",
    "tests/onshape/dome_hemisphere_20_fine.step"
);
test_stl_step_stage23_compare!(
    onshape_cylinder_stage23_compare,
    "tests/onshape/cylinder_10x30_medium.stl",
    "tests/onshape/cylinder_10x30_medium.step"
);
test_stl_step_stage23_compare!(
    onshape_l_bracket_stage23_compare,
    "tests/onshape/l_bracket_simple_medium.stl",
    "tests/onshape/l_bracket_simple_medium.step"
);
test_stl_step_stage23_compare!(
    onshape_plate_with_hole_stage23_compare,
    "tests/onshape/plate_with_hole_100x50_coarse.stl",
    "tests/onshape/plate_with_hole_100x50_coarse.step"
);
test_stl_step_stage23_compare!(
    onshape_pipe_elbow_stage23_compare,
    "tests/onshape/pipe_elbow_10_fine.stl",
    "tests/onshape/pipe_elbow_10_fine.step"
);

// --- tests/fusion/ ---
test_stl_step_stage23_compare!(
    fusion_plate_low_stage23_compare,
    "tests/fusion/plate_with_hole_100x50_low.stl",
    "tests/fusion/plate_with_hole_100x50_low.step"
);
test_stl_step_stage23_compare!(
    fusion_plate_medium_stage23_compare,
    "tests/fusion/plate_with_hole_100x50_medium.stl",
    "tests/fusion/plate_with_hole_100x50_medium.step"
);
test_stl_step_stage23_compare!(
    fusion_plate_high_stage23_compare,
    "tests/fusion/plate_with_hole_100x50_high.stl",
    "tests/fusion/plate_with_hole_100x50_high.step"
);

// ===========================================================================
// Stage 2.6: Surface selection tests
// ===========================================================================

fn config_for_stl_stage26(stl_path: &str) -> config::Config {
    config::parse_config_from(["brepper", stl_path, "--stage=2.6", "-q"]).unwrap()
}

fn config_for_compare_stage26(stl_path: &str, step_path: &str) -> config::Config {
    let mut config =
        config::parse_config_from(["brepper", stl_path, "--compare", step_path, "--stage=2.6", "-q"])
            .unwrap();
    config.load_compare_step().unwrap();
    config
}

/// Generate a test that runs stage 2.6 with --compare against a STEP file.
macro_rules! test_stl_step_stage26_compare {
    ($name:ident, $stl_path:literal, $step_path:literal) => {
        #[test]
        fn $name() {
            let stl = format!("{}/{}", manifest_dir(), $stl_path);
            let step = format!("{}/{}", manifest_dir(), $step_path);
            let config = config_for_compare_stage26(&stl, &step);
            run_stage2(&config);
        }
    };
}

// ---------------------------------------------------------------------------
// Surface selection count tests
// ---------------------------------------------------------------------------

#[test]
fn cube_surface_selection_all_planar() {
    let stl = format!("{}/tests/ccad/generated/cube.stl", manifest_dir());
    let config = config_for_stl_stage26(&stl);
    let output = run_stage2(&config);

    assert!(!output.selected_surfaces.is_empty(), "should have selected surfaces");
    // Cube: 6 faces, all single-face planar hypotheses, no other hypotheses
    let planar_count = output.selected_surfaces.iter()
        .filter(|s| matches!(s, stage2::SelectedSurface::Planar(_)))
        .count();
    assert_eq!(planar_count, output.selected_surfaces.len(), "cube should select only planar surfaces");
}

#[test]
fn block_with_hole_surface_selection() {
    let stl = format!("{}/tests/ccad/generated/block_with_hole.stl", manifest_dir());
    let config = config_for_stl_stage26(&stl);
    let output = run_stage2(&config);

    // block_with_hole: 6 planar faces + 1 cylindrical hole = 7 selected surfaces
    let planar_count = output.selected_surfaces.iter()
        .filter(|s| matches!(s, stage2::SelectedSurface::Planar(_)))
        .count();
    let cyl_count = output.selected_surfaces.iter()
        .filter(|s| matches!(s, stage2::SelectedSurface::Cylindrical(_)))
        .count();
    assert_eq!(planar_count, 6, "block_with_hole should have 6 planar surfaces");
    assert_eq!(cyl_count, 1, "block_with_hole should have 1 cylindrical surface");
    assert_eq!(output.selected_surfaces.len(), 7, "block_with_hole should have 7 total selected surfaces");
}

#[test]
fn pipe_surface_selection() {
    let stl = format!("{}/tests/ccad/generated/pipe.stl", manifest_dir());
    let config = config_for_stl_stage26(&stl);
    let output = run_stage2(&config);

    // pipe: 2 planar annular ends + 2 cylindrical surfaces (inner + outer) = 4
    let planar_count = output.selected_surfaces.iter()
        .filter(|s| matches!(s, stage2::SelectedSurface::Planar(_)))
        .count();
    let cyl_count = output.selected_surfaces.iter()
        .filter(|s| matches!(s, stage2::SelectedSurface::Cylindrical(_)))
        .count();
    assert_eq!(planar_count, 2, "pipe should have 2 planar surfaces");
    assert_eq!(cyl_count, 2, "pipe should have 2 cylindrical surfaces");
}

#[test]
fn ball_on_cylinder_surface_selection() {
    let stl = format!("{}/tests/ccad/generated/ball_on_cylinder.stl", manifest_dir());
    let config = config_for_stl_stage26(&stl);
    let output = run_stage2(&config);

    // ball_on_cylinder: 1 planar base + 1 cylindrical stalk + 1 spherical ball = 3
    let planar_count = output.selected_surfaces.iter()
        .filter(|s| matches!(s, stage2::SelectedSurface::Planar(_)))
        .count();
    let cyl_count = output.selected_surfaces.iter()
        .filter(|s| matches!(s, stage2::SelectedSurface::Cylindrical(_)))
        .count();
    let sph_count = output.selected_surfaces.iter()
        .filter(|s| matches!(s, stage2::SelectedSurface::Spherical(_)))
        .count();
    assert_eq!(planar_count, 1, "ball_on_cylinder should have 1 planar surface");
    assert_eq!(cyl_count, 1, "ball_on_cylinder should have 1 cylindrical surface");
    assert_eq!(sph_count, 1, "ball_on_cylinder should have 1 spherical surface");
}

#[test]
fn hemisphere_surface_selection() {
    let stl = format!("{}/tests/ccad/generated/hemisphere.stl", manifest_dir());
    let config = config_for_stl_stage26(&stl);
    let output = run_stage2(&config);

    // hemisphere: 1 planar base + 1 spherical dome = 2
    let planar_count = output.selected_surfaces.iter()
        .filter(|s| matches!(s, stage2::SelectedSurface::Planar(_)))
        .count();
    let sph_count = output.selected_surfaces.iter()
        .filter(|s| matches!(s, stage2::SelectedSurface::Spherical(_)))
        .count();
    assert_eq!(planar_count, 1, "hemisphere should have 1 planar surface");
    assert_eq!(sph_count, 1, "hemisphere should have 1 spherical surface");
    assert_eq!(output.selected_surfaces.len(), 2, "hemisphere should have 2 total");
}

#[test]
fn simple_sphere_surface_selection() {
    let stl = format!("{}/tests/ccad/generated/simple_sphere.stl", manifest_dir());
    let config = config_for_stl_stage26(&stl);
    let output = run_stage2(&config);

    // simple_sphere: 1 spherical surface
    let sph_count = output.selected_surfaces.iter()
        .filter(|s| matches!(s, stage2::SelectedSurface::Spherical(_)))
        .count();
    assert_eq!(sph_count, 1, "simple_sphere should have 1 spherical surface");
    assert_eq!(output.selected_surfaces.len(), 1, "simple_sphere should have 1 total");
}

#[test]
fn all_faces_covered_by_selection() {
    use std::collections::HashSet;
    // For every model, every face must be covered by exactly one selected surface
    let models = [
        "tests/ccad/generated/cube.stl",
        "tests/ccad/generated/block_with_hole.stl",
        "tests/ccad/generated/pipe.stl",
        "tests/ccad/generated/ball_on_cylinder.stl",
        "tests/ccad/generated/simple_sphere.stl",
        "tests/ccad/generated/hemisphere.stl",
    ];
    for model in &models {
        let stl = format!("{}/{}", manifest_dir(), model);
        let config = config_for_stl_stage26(&stl);
        let output = run_stage2(&config);

        // Collect all unique face indices from selected surfaces
        let mut covered_faces = HashSet::new();
        for s in &output.selected_surfaces {
            let faces = match s {
                stage2::SelectedSurface::Planar(i) => &output.planar_hypotheses[*i].faces,
                stage2::SelectedSurface::Cylindrical(i) => &output.cylindrical_hypotheses[*i].faces,
                stage2::SelectedSurface::Spherical(i) => &output.spherical_hypotheses[*i].faces,
            };
            for &fi in faces {
                covered_faces.insert(fi);
            }
        }

        assert_eq!(
            covered_faces.len(),
            output.mesh.faces.len(),
            "{}: expected {} unique faces covered, got {}",
            model, output.mesh.faces.len(), covered_faces.len()
        );
    }
}

// ---------------------------------------------------------------------------
// Compare tests: all test STL/STEP pairs should pass --compare at stage 2.6
// ---------------------------------------------------------------------------

// --- tests/manual/ ---
test_stl_step_stage26_compare!(
    manual_cube_stage26_compare,
    "tests/manual/cube.stl",
    "tests/manual/cube.step"
);

// --- tests/ccad/generated/ ---
test_stl_step_stage26_compare!(
    ccad_cube_stage26_compare,
    "tests/ccad/generated/cube.stl",
    "tests/ccad/generated/cube.step"
);
test_stl_step_stage26_compare!(
    ccad_wedge_stage26_compare,
    "tests/ccad/generated/wedge.stl",
    "tests/ccad/generated/wedge.step"
);
test_stl_step_stage26_compare!(
    ccad_t_shape_stage26_compare,
    "tests/ccad/generated/t_shape.stl",
    "tests/ccad/generated/t_shape.step"
);
test_stl_step_stage26_compare!(
    ccad_staircase_stage26_compare,
    "tests/ccad/generated/staircase.stl",
    "tests/ccad/generated/staircase.step"
);
test_stl_step_stage26_compare!(
    ccad_simple_cylinder_stage26_compare,
    "tests/ccad/generated/simple_cylinder.stl",
    "tests/ccad/generated/simple_cylinder.step"
);
test_stl_step_stage26_compare!(
    ccad_block_with_hole_stage26_compare,
    "tests/ccad/generated/block_with_hole.stl",
    "tests/ccad/generated/block_with_hole.step"
);
test_stl_step_stage26_compare!(
    ccad_pipe_stage26_compare,
    "tests/ccad/generated/pipe.stl",
    "tests/ccad/generated/pipe.step"
);
test_stl_step_stage26_compare!(
    ccad_stepped_cylinder_stage26_compare,
    "tests/ccad/generated/stepped_cylinder.stl",
    "tests/ccad/generated/stepped_cylinder.step"
);
test_stl_step_stage26_compare!(
    ccad_two_holes_stage26_compare,
    "tests/ccad/generated/two_holes.stl",
    "tests/ccad/generated/two_holes.step"
);
test_stl_step_stage26_compare!(
    ccad_simple_sphere_stage26_compare,
    "tests/ccad/generated/simple_sphere.stl",
    "tests/ccad/generated/simple_sphere.step"
);
test_stl_step_stage26_compare!(
    ccad_hemisphere_stage26_compare,
    "tests/ccad/generated/hemisphere.stl",
    "tests/ccad/generated/hemisphere.step"
);
test_stl_step_stage26_compare!(
    ccad_spherical_pocket_stage26_compare,
    "tests/ccad/generated/spherical_pocket.stl",
    "tests/ccad/generated/spherical_pocket.step"
);
test_stl_step_stage26_compare!(
    ccad_ball_on_cylinder_stage26_compare,
    "tests/ccad/generated/ball_on_cylinder.stl",
    "tests/ccad/generated/ball_on_cylinder.step"
);

// --- tests/onshape/ ---
test_stl_step_stage26_compare!(
    onshape_chamfered_cube_stage26_compare,
    "tests/onshape/chamfered_cube_10_c1_medium.stl",
    "tests/onshape/chamfered_cube_10_c1_medium.step"
);
test_stl_step_stage26_compare!(
    onshape_stepped_block_stage26_compare,
    "tests/onshape/stepped_block_coarse.stl",
    "tests/onshape/stepped_block_coarse.step"
);
test_stl_step_stage26_compare!(
    onshape_l_bracket_stage26_compare,
    "tests/onshape/l_bracket_simple_medium.stl",
    "tests/onshape/l_bracket_simple_medium.step"
);
test_stl_step_stage26_compare!(
    onshape_cylinder_stage26_compare,
    "tests/onshape/cylinder_10x30_medium.stl",
    "tests/onshape/cylinder_10x30_medium.step"
);
test_stl_step_stage26_compare!(
    onshape_sphere_stage26_compare,
    "tests/onshape/sphere_25_fine.stl",
    "tests/onshape/sphere_25_fine.step"
);
test_stl_step_stage26_compare!(
    onshape_dome_hemisphere_stage26_compare,
    "tests/onshape/dome_hemisphere_20_fine.stl",
    "tests/onshape/dome_hemisphere_20_fine.step"
);
test_stl_step_stage26_compare!(
    onshape_pipe_elbow_stage26_compare,
    "tests/onshape/pipe_elbow_10_fine.stl",
    "tests/onshape/pipe_elbow_10_fine.step"
);
test_stl_step_stage26_compare!(
    onshape_plate_with_hole_stage26_compare,
    "tests/onshape/plate_with_hole_100x50_coarse.stl",
    "tests/onshape/plate_with_hole_100x50_coarse.step"
);
test_stl_step_stage26_compare!(
    onshape_rounded_cube_stage26_compare,
    "tests/onshape/rounded_cube_10_r2_fine.stl",
    "tests/onshape/rounded_cube_10_r2_fine.step"
);
test_stl_step_stage26_compare!(
    onshape_cone_stage26_compare,
    "tests/onshape/cone_15x20_medium.stl",
    "tests/onshape/cone_15x20_medium.step"
);

// --- tests/fusion/ ---
test_stl_step_stage26_compare!(
    fusion_plate_low_stage26_compare,
    "tests/fusion/plate_with_hole_100x50_low.stl",
    "tests/fusion/plate_with_hole_100x50_low.step"
);
test_stl_step_stage26_compare!(
    fusion_plate_medium_stage26_compare,
    "tests/fusion/plate_with_hole_100x50_medium.stl",
    "tests/fusion/plate_with_hole_100x50_medium.step"
);
test_stl_step_stage26_compare!(
    fusion_plate_high_stage26_compare,
    "tests/fusion/plate_with_hole_100x50_high.stl",
    "tests/fusion/plate_with_hole_100x50_high.step"
);
