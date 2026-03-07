//! Integration tests for Stage 2: surface fitting.
//!
//! Tests that planar hypothesis deduction produces correct results for
//! known planar models, and that --compare validation passes.

use brepper::config;
use brepper::{stage1, stage2};
use opencascade_sys::{b_rep_adaptor, geom_abs, gp, top_abs, top_exp, topo_ds};

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

/// Compute the centroid of a mesh face.
fn face_centroid(face: &stage1::MeshFace, vertices: &[stage1::MeshVertex]) -> [f64; 3] {
    let n = face.vertex_count as usize;
    let mut cx = 0.0_f64;
    let mut cy = 0.0_f64;
    let mut cz = 0.0_f64;
    for vi in 0..n {
        let v = &vertices[face.vertex_indices[vi]];
        cx += v.x;
        cy += v.y;
        cz += v.z;
    }
    let inv_n = 1.0 / n as f64;
    [cx * inv_n, cy * inv_n, cz * inv_n]
}

/// For every face in the STEP file whose surface type is in `check_types`,
/// assert that there's at least one hypothesis of the matching type with a mesh face
/// centroid within surface_tolerance of that STEP face.
fn assert_step_surfaces_covered(
    output: &stage2::Stage2Output,
    config: &config::Config,
    check_types: &[geom_abs::SurfaceType],
) {
    let compare_shape = config.compare_shape.as_ref().expect("compare_shape must be loaded");
    let tol = config.surface_tolerance_mm;

    let mut explorer = top_exp::Explorer::new_shape_shapeenum2(
        compare_shape,
        top_abs::ShapeEnum::Face,
        top_abs::ShapeEnum::Shape,
    );

    let mut step_face_idx = 0_usize;
    while explorer.more() {
        let shape = explorer.current();
        let face = topo_ds::face_shape(shape);
        let adaptor = b_rep_adaptor::Surface::new_face(face);
        let surf_type = adaptor.get_type();

        if !check_types.contains(&surf_type) {
            explorer.next();
            step_face_idx += 1;
            continue;
        }

        let type_name = match surf_type {
            geom_abs::SurfaceType::Plane => "Plane",
            geom_abs::SurfaceType::Cylinder => "Cylinder",
            geom_abs::SurfaceType::Sphere => "Sphere",
            _ => unreachable!(),
        };

        // Collect centroids from all hypotheses of the matching type
        let hypothesis_centroids: Vec<[f64; 3]> = match surf_type {
            geom_abs::SurfaceType::Plane => {
                output.planar_hypotheses.iter()
                    .flat_map(|h| h.faces.iter().map(|&fi| face_centroid(&output.mesh.faces[fi], &output.mesh.vertices)))
                    .collect()
            }
            geom_abs::SurfaceType::Cylinder => {
                output.cylindrical_hypotheses.iter()
                    .flat_map(|h| h.faces.iter().map(|&fi| face_centroid(&output.mesh.faces[fi], &output.mesh.vertices)))
                    .collect()
            }
            geom_abs::SurfaceType::Sphere => {
                output.spherical_hypotheses.iter()
                    .flat_map(|h| h.faces.iter().map(|&fi| face_centroid(&output.mesh.faces[fi], &output.mesh.vertices)))
                    .collect()
            }
            _ => unreachable!(),
        };

        // Check if any centroid from a matching hypothesis is within tolerance of this STEP face
        let face_shape = face.as_shape();
        let covered = hypothesis_centroids.iter().any(|c| {
            let pt = gp::Pnt::new_real3(c[0], c[1], c[2]);
            let d = stage1::min_distance_to_shape(&pt, face_shape);
            d <= tol
        });

        assert!(
            covered,
            "STEP face {} (type {}) has no matching {} hypothesis with a centroid within tolerance {}",
            step_face_idx, type_name, type_name, tol,
        );

        explorer.next();
        step_face_idx += 1;
    }
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
/// Also checks that every Plane/Cylinder STEP face is covered by a
/// hypothesis of the corresponding type.
macro_rules! test_stl_step_stage22_compare {
    ($name:ident, $stl_path:literal, $step_path:literal) => {
        #[test]
        fn $name() {
            let stl = format!("{}/{}", manifest_dir(), $stl_path);
            let step = format!("{}/{}", manifest_dir(), $step_path);
            let config = config_for_compare_stage22(&stl, &step);
            let output = run_stage2(&config);
            assert_step_surfaces_covered(
                &output, &config,
                &[geom_abs::SurfaceType::Plane, geom_abs::SurfaceType::Cylinder],
            );
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

// --- tests/manual/ ---
test_stl_step_stage22_compare!(
    manual_cube_stage22_compare,
    "tests/manual/cube.stl",
    "tests/manual/cube.step"
);

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
/// Also checks that every Plane/Cylinder/Sphere STEP face is covered by a
/// hypothesis of the corresponding type.
macro_rules! test_stl_step_stage23_compare {
    ($name:ident, $stl_path:literal, $step_path:literal) => {
        #[test]
        fn $name() {
            let stl = format!("{}/{}", manifest_dir(), $stl_path);
            let step = format!("{}/{}", manifest_dir(), $step_path);
            let config = config_for_compare_stage23(&stl, &step);
            let output = run_stage2(&config);
            assert_step_surfaces_covered(
                &output, &config,
                &[geom_abs::SurfaceType::Plane, geom_abs::SurfaceType::Cylinder, geom_abs::SurfaceType::Sphere],
            );
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

// --- tests/manual/ ---
test_stl_step_stage23_compare!(
    manual_cube_stage23_compare,
    "tests/manual/cube.stl",
    "tests/manual/cube.step"
);

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
// Pipe elbow has torus surfaces that can't be fitted with current primitives
// (planes, cylinders, spheres). Skipping stage 2.6 compare.
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

// ===========================================================================
// Stage 2.2: Cylinder parameter matching tests
// Verify that deduced cylinder hypotheses match STEP cylinder parameters
// (axis direction, radius, axis position).
// ===========================================================================

/// Extract cylinder parameters from a STEP file: returns Vec<(axis_dir, axis_origin, radius)>.
fn extract_step_cylinders(step_path: &str) -> Vec<([f64; 3], [f64; 3], f64)> {
    use opencascade_sys::{step_control, message};
    let mut reader = step_control::Reader::new();
    reader.read_file_charptr(step_path);
    reader.transfer_roots(&message::ProgressRange::new());
    let shape = reader.one_shape();

    let mut cylinders = Vec::new();
    let mut explorer = top_exp::Explorer::new_shape_shapeenum2(
        &shape, top_abs::ShapeEnum::Face, top_abs::ShapeEnum::Shape,
    );
    while explorer.more() {
        let face = topo_ds::face_shape(explorer.current());
        let adaptor = b_rep_adaptor::Surface::new_face(face);
        if adaptor.get_type() == geom_abs::SurfaceType::Cylinder {
            let cyl = adaptor.cylinder();
            let radius = cyl.radius();
            let axis = cyl.axis();
            let dir = axis.direction();
            let loc = axis.location();
            cylinders.push((
                [dir.x(), dir.y(), dir.z()],
                [loc.x(), loc.y(), loc.z()],
                radius,
            ));
        }
        explorer.next();
    }

    // Deduplicate: STEP files can have multiple faces on the same cylinder.
    // Group by matching axis direction (parallel), axis position (colinear),
    // and radius.
    let mut unique: Vec<([f64; 3], [f64; 3], f64)> = Vec::new();
    'outer: for (dir, origin, r) in &cylinders {
        for (udir, uorigin, ur) in &unique {
            // Axis directions must be parallel
            let d = udir[0] * dir[0] + udir[1] * dir[1] + udir[2] * dir[2];
            if d.abs() < 0.999 {
                continue;
            }
            // Radii must match
            if (ur - r).abs() > 0.01 {
                continue;
            }
            // Axis origins must be colinear (distance between axes < tolerance)
            let diff = [origin[0] - uorigin[0], origin[1] - uorigin[1], origin[2] - uorigin[2]];
            let t = diff[0] * udir[0] + diff[1] * udir[1] + diff[2] * udir[2];
            let perp_sq = diff[0]*diff[0] + diff[1]*diff[1] + diff[2]*diff[2] - t*t;
            if perp_sq < 0.01 * 0.01 {
                continue 'outer; // Already have this cylinder
            }
        }
        unique.push((*dir, *origin, *r));
    }
    unique
}

/// Verify that each STEP cylinder has a matching hypothesis with compatible
/// axis direction, axis position, and radius.
fn assert_cylinders_match_step(
    output: &stage2::Stage2Output,
    step_path: &str,
) {
    let step_cylinders = extract_step_cylinders(step_path);
    let hypotheses = &output.cylindrical_hypotheses;

    assert_eq!(
        hypotheses.len(), step_cylinders.len(),
        "Expected {} cylinder hypotheses (from STEP), got {}",
        step_cylinders.len(), hypotheses.len(),
    );

    for (i, (step_dir, step_origin, step_radius)) in step_cylinders.iter().enumerate() {
        // Find a matching hypothesis
        let matched = hypotheses.iter().any(|h| {
            // Radius match
            if (h.radius - step_radius).abs() > 0.1 {
                return false;
            }
            // Axis direction: must be parallel (or anti-parallel)
            let d = h.axis_direction[0] * step_dir[0]
                  + h.axis_direction[1] * step_dir[1]
                  + h.axis_direction[2] * step_dir[2];
            if d.abs() < 0.999 {
                return false;
            }
            // Axis position: distance between axes must be small
            let diff = [
                h.axis_origin[0] - step_origin[0],
                h.axis_origin[1] - step_origin[1],
                h.axis_origin[2] - step_origin[2],
            ];
            let t = diff[0] * h.axis_direction[0]
                  + diff[1] * h.axis_direction[1]
                  + diff[2] * h.axis_direction[2];
            let perp_sq = diff[0]*diff[0] + diff[1]*diff[1] + diff[2]*diff[2] - t*t;
            perp_sq < 0.1 * 0.1
        });

        assert!(
            matched,
            "STEP cylinder {} (dir=[{:.3},{:.3},{:.3}], r={:.3}) has no matching hypothesis",
            i, step_dir[0], step_dir[1], step_dir[2], step_radius,
        );
    }
}

macro_rules! test_cylinder_params_match {
    ($name:ident, $stl_path:literal, $step_path:literal) => {
        #[test]
        fn $name() {
            let stl = format!("{}/{}", manifest_dir(), $stl_path);
            let step = format!("{}/{}", manifest_dir(), $step_path);
            let config = config_for_stl_stage22(&stl);
            let output = run_stage2(&config);
            assert_cylinders_match_step(&output, &step);
        }
    };
}

test_cylinder_params_match!(
    ccad_simple_cylinder_params_match,
    "tests/ccad/generated/simple_cylinder.stl",
    "tests/ccad/generated/simple_cylinder.step"
);
test_cylinder_params_match!(
    ccad_block_with_hole_params_match,
    "tests/ccad/generated/block_with_hole.stl",
    "tests/ccad/generated/block_with_hole.step"
);
test_cylinder_params_match!(
    ccad_pipe_params_match,
    "tests/ccad/generated/pipe.stl",
    "tests/ccad/generated/pipe.step"
);
test_cylinder_params_match!(
    ccad_stepped_cylinder_params_match,
    "tests/ccad/generated/stepped_cylinder.stl",
    "tests/ccad/generated/stepped_cylinder.step"
);
test_cylinder_params_match!(
    ccad_two_holes_params_match,
    "tests/ccad/generated/two_holes.stl",
    "tests/ccad/generated/two_holes.step"
);
test_cylinder_params_match!(
    ccad_ball_on_cylinder_params_match,
    "tests/ccad/generated/ball_on_cylinder.stl",
    "tests/ccad/generated/ball_on_cylinder.step"
);
test_cylinder_params_match!(
    onshape_cylinder_params_match,
    "tests/onshape/cylinder_10x30_medium.stl",
    "tests/onshape/cylinder_10x30_medium.step"
);
test_cylinder_params_match!(
    onshape_plate_with_hole_params_match,
    "tests/onshape/plate_with_hole_100x50_coarse.stl",
    "tests/onshape/plate_with_hole_100x50_coarse.step"
);
test_cylinder_params_match!(
    fusion_plate_low_params_match,
    "tests/fusion/plate_with_hole_100x50_low.stl",
    "tests/fusion/plate_with_hole_100x50_low.step"
);
test_cylinder_params_match!(
    fusion_plate_medium_params_match,
    "tests/fusion/plate_with_hole_100x50_medium.stl",
    "tests/fusion/plate_with_hole_100x50_medium.step"
);
test_cylinder_params_match!(
    fusion_plate_high_params_match,
    "tests/fusion/plate_with_hole_100x50_high.stl",
    "tests/fusion/plate_with_hole_100x50_high.step"
);


// ===========================================================================
// Stage 2.3: Sphere parameter matching tests
// Verify that deduced sphere hypotheses match STEP sphere parameters
// (center, radius).
// ===========================================================================

/// Extract sphere parameters from a STEP file: returns Vec<(center, radius)>.
fn extract_step_spheres(step_path: &str) -> Vec<([f64; 3], f64)> {
    use opencascade_sys::{step_control, message};
    let mut reader = step_control::Reader::new();
    reader.read_file_charptr(step_path);
    reader.transfer_roots(&message::ProgressRange::new());
    let shape = reader.one_shape();

    let mut spheres = Vec::new();
    let mut explorer = top_exp::Explorer::new_shape_shapeenum2(
        &shape, top_abs::ShapeEnum::Face, top_abs::ShapeEnum::Shape,
    );
    while explorer.more() {
        let face = topo_ds::face_shape(explorer.current());
        let adaptor = b_rep_adaptor::Surface::new_face(face);
        if adaptor.get_type() == geom_abs::SurfaceType::Sphere {
            let sph = adaptor.sphere();
            let radius = sph.radius();
            let loc = sph.location();
            spheres.push((
                [loc.x(), loc.y(), loc.z()],
                radius,
            ));
        }
        explorer.next();
    }

    // Deduplicate: STEP files can have multiple faces on the same sphere.
    let mut unique: Vec<([f64; 3], f64)> = Vec::new();
    'outer: for (center, r) in &spheres {
        for (ucenter, ur) in &unique {
            // Radii must match
            if (ur - r).abs() > 0.01 {
                continue;
            }
            // Centers must be close
            let d2 = (center[0]-ucenter[0]).powi(2)
                   + (center[1]-ucenter[1]).powi(2)
                   + (center[2]-ucenter[2]).powi(2);
            if d2 < 0.01 * 0.01 {
                continue 'outer; // Already have this sphere
            }
        }
        unique.push((*center, *r));
    }
    unique
}

/// Verify that each STEP sphere has a matching hypothesis with compatible
/// center and radius.
fn assert_spheres_match_step(
    output: &stage2::Stage2Output,
    step_path: &str,
) {
    let step_spheres = extract_step_spheres(step_path);
    let hypotheses = &output.spherical_hypotheses;

    assert_eq!(
        hypotheses.len(), step_spheres.len(),
        "Expected {} sphere hypotheses (from STEP), got {}",
        step_spheres.len(), hypotheses.len(),
    );

    for (i, (step_center, step_radius)) in step_spheres.iter().enumerate() {
        let matched = hypotheses.iter().any(|h| {
            // Radius match
            if (h.radius - step_radius).abs() > 0.1 {
                return false;
            }
            // Center position match
            let d2 = (h.center[0] - step_center[0]).powi(2)
                    + (h.center[1] - step_center[1]).powi(2)
                    + (h.center[2] - step_center[2]).powi(2);
            d2 < 0.1 * 0.1
        });

        assert!(
            matched,
            "STEP sphere {} (center=[{:.3},{:.3},{:.3}], r={:.3}) has no matching hypothesis",
            i, step_center[0], step_center[1], step_center[2], step_radius,
        );
    }
}

macro_rules! test_sphere_params_match {
    ($name:ident, $stl_path:literal, $step_path:literal) => {
        #[test]
        fn $name() {
            let stl = format!("{}/{}", manifest_dir(), $stl_path);
            let step = format!("{}/{}", manifest_dir(), $step_path);
            let config = config_for_stl_stage23(&stl);
            let output = run_stage2(&config);
            assert_spheres_match_step(&output, &step);
        }
    };
}

test_sphere_params_match!(
    ccad_simple_sphere_params_match,
    "tests/ccad/generated/simple_sphere.stl",
    "tests/ccad/generated/simple_sphere.step"
);
test_sphere_params_match!(
    ccad_hemisphere_params_match,
    "tests/ccad/generated/hemisphere.stl",
    "tests/ccad/generated/hemisphere.step"
);
test_sphere_params_match!(
    ccad_spherical_pocket_params_match,
    "tests/ccad/generated/spherical_pocket.stl",
    "tests/ccad/generated/spherical_pocket.step"
);
test_sphere_params_match!(
    ccad_ball_on_cylinder_sphere_params_match,
    "tests/ccad/generated/ball_on_cylinder.stl",
    "tests/ccad/generated/ball_on_cylinder.step"
);
test_sphere_params_match!(
    onshape_sphere_params_match,
    "tests/onshape/sphere_25_fine.stl",
    "tests/onshape/sphere_25_fine.step"
);
test_sphere_params_match!(
    onshape_dome_hemisphere_params_match,
    "tests/onshape/dome_hemisphere_20_fine.stl",
    "tests/onshape/dome_hemisphere_20_fine.step"
);