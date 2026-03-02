//! Integration tests for Stage 3: surface reconstruction.
//!
//! Tests that stage 3.1 (OCCT surface creation and adjacency graph construction)
//! produces correct topology for known models.

use brepper::config;
use brepper::{stage1, stage2, stage3};

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

fn manifest_dir() -> String {
    env!("CARGO_MANIFEST_DIR").to_string()
}

fn config_for_stl(stl_path: &str) -> config::Config {
    config::parse_config_from(["brepper", stl_path, "--stage=3.1", "-q"]).unwrap()
}

fn config_for_compare(stl_path: &str, step_path: &str) -> config::Config {
    let mut config =
        config::parse_config_from(["brepper", stl_path, "--compare", step_path, "--stage=3.1", "-q"])
            .unwrap();
    config.load_compare_step().unwrap();
    config
}

/// Run stages 1, 2, and 3.1, returning the Stage3Output.
fn run_stage3(config: &config::Config) -> stage3::Stage3Output {
    let mesh = stage1::stage1(config).expect("stage1 should pass");
    let surfaces = stage2::stage2(config, mesh).expect("stage2 should pass");
    stage3::stage3(config, surfaces).expect("stage3 should pass")
}

// ---------------------------------------------------------------------------
// Topology invariant checks
// ---------------------------------------------------------------------------

/// Verify Euler's formula: V - E + F = 2 * (shells - holes).
/// For simple closed solids with no holes: V - E + F = 2.
fn check_euler_formula(output: &stage3::Stage3Output) {
    let v = output.vertices.len() as i64;
    let e = output.edges.len() as i64;
    let f = output.face_descriptors.len() as i64;
    let euler = v - e + f;
    // For all current test models, the Euler characteristic should be 2
    // (single connected closed solid, no through-holes in the B-Rep topology)
    assert_eq!(euler, 2, "Euler formula V({v})-E({e})+F({f}) = {euler}, expected 2");
}

/// Verify that every edge references valid face and vertex indices.
fn check_edge_validity(output: &stage3::Stage3Output) {
    let nf = output.face_descriptors.len();
    let nv = output.vertices.len();
    for (ei, edge) in output.edges.iter().enumerate() {
        for &fi in &edge.face_indices {
            assert!(fi < nf, "edge {ei} references invalid face {fi}");
        }
        for &vi in &edge.vertex_indices {
            if vi != usize::MAX {
                assert!(vi < nv, "edge {ei} references invalid vertex {vi}");
            }
        }
        assert!(
            edge.mesh_boundary_vertices.len() >= 2,
            "edge {ei} has < 2 mesh boundary vertices"
        );
    }
}

/// Verify face adjacency is symmetric: if face A lists B as adjacent, then B lists A.
fn check_adjacency_symmetry(output: &stage3::Stage3Output) {
    for (fi, fd) in output.face_descriptors.iter().enumerate() {
        for &adj_fi in &fd.adjacent_faces {
            let adj_fd = &output.face_descriptors[adj_fi];
            assert!(
                adj_fd.adjacent_faces.contains(&fi),
                "face {fi} lists face {adj_fi} as adjacent, but {adj_fi} does not list {fi}"
            );
        }
    }
}

/// Verify that edge face_indices are consistent with face adjacency.
fn check_edge_face_consistency(output: &stage3::Stage3Output) {
    for (ei, edge) in output.edges.iter().enumerate() {
        let [f0, f1] = edge.face_indices;
        let fd0 = &output.face_descriptors[f0];
        let fd1 = &output.face_descriptors[f1];
        assert!(
            fd0.edge_indices.contains(&ei),
            "edge {ei} lists face {f0}, but face {f0} doesn't reference edge {ei}"
        );
        assert!(
            fd1.edge_indices.contains(&ei),
            "edge {ei} lists face {f1}, but face {f1} doesn't reference edge {ei}"
        );
    }
}

/// Run all standard topology checks on a Stage3Output.
fn check_topology(output: &stage3::Stage3Output) {
    check_edge_validity(output);
    check_adjacency_symmetry(output);
    check_edge_face_consistency(output);
}

// ---------------------------------------------------------------------------
// Cube tests
// ---------------------------------------------------------------------------

#[test]
fn cube_topology() {
    let stl = format!("{}/tests/ccad/generated/cube.stl", manifest_dir());
    let config = config_for_stl(&stl);
    let output = run_stage3(&config);

    assert_eq!(output.face_descriptors.len(), 6, "cube should have 6 faces");
    assert_eq!(output.edges.len(), 12, "cube should have 12 edges");
    assert_eq!(output.vertices.len(), 8, "cube should have 8 vertices");

    // Each face should have exactly 4 adjacent faces, 4 edges, 4 vertices
    for (fi, fd) in output.face_descriptors.iter().enumerate() {
        assert_eq!(fd.adjacent_faces.len(), 4, "cube face {fi} should have 4 adjacent faces");
        assert_eq!(fd.edge_indices.len(), 4, "cube face {fi} should have 4 edges");
        assert_eq!(fd.vertex_indices.len(), 4, "cube face {fi} should have 4 vertices");
    }

    // Each vertex should have exactly 3 adjacent faces and 3 edges
    for (vi, v) in output.vertices.iter().enumerate() {
        assert_eq!(v.adjacent_faces.len(), 3, "cube vertex {vi} should have 3 adjacent faces");
        assert_eq!(v.adjacent_edges.len(), 3, "cube vertex {vi} should have 3 adjacent edges");
    }

    check_euler_formula(&output);
    check_topology(&output);
}

#[test]
fn wedge_topology() {
    let stl = format!("{}/tests/ccad/generated/wedge.stl", manifest_dir());
    let config = config_for_stl(&stl);
    let output = run_stage3(&config);

    // Wedge (truncated box with narrowed top): 6 faces, 12 edges, 8 vertices
    assert_eq!(output.face_descriptors.len(), 6, "wedge should have 6 faces");
    assert_eq!(output.edges.len(), 12, "wedge should have 12 edges");
    assert_eq!(output.vertices.len(), 8, "wedge should have 8 vertices");

    check_euler_formula(&output);
    check_topology(&output);
}

// ---------------------------------------------------------------------------
// Cylinder tests
// ---------------------------------------------------------------------------

#[test]
fn simple_cylinder_topology() {
    let stl = format!("{}/tests/ccad/generated/simple_cylinder.stl", manifest_dir());
    let config = config_for_stl(&stl);
    let output = run_stage3(&config);

    // Cylinder: 3 surfaces (cylinder + 2 end caps), 2 edges (circles), 0 vertices
    assert_eq!(output.face_descriptors.len(), 3, "cylinder should have 3 faces");
    assert_eq!(output.edges.len(), 2, "cylinder should have 2 edges");
    assert_eq!(output.vertices.len(), 0, "cylinder should have 0 vertices");

    check_topology(&output);
}

#[test]
fn block_with_hole_topology() {
    let stl = format!("{}/tests/ccad/generated/block_with_hole.stl", manifest_dir());
    let config = config_for_stl(&stl);
    let output = run_stage3(&config);

    // Block with hole: 7 surfaces (6 planar + 1 cylindrical)
    assert_eq!(output.face_descriptors.len(), 7, "block_with_hole should have 7 faces");

    // The cylindrical face should exist
    let cyl_faces: Vec<_> = output.face_descriptors.iter().enumerate()
        .filter(|(_, fd)| {
            matches!(
                output.stage2.selected_surfaces[fd.selected_surface_idx],
                stage2::SelectedSurface::Cylindrical(_)
            )
        })
        .collect();
    assert_eq!(cyl_faces.len(), 1, "should have 1 cylindrical face");

    check_topology(&output);
}

#[test]
fn pipe_topology() {
    let stl = format!("{}/tests/ccad/generated/pipe.stl", manifest_dir());
    let config = config_for_stl(&stl);
    let output = run_stage3(&config);

    // Pipe: 4 surfaces (2 cylindrical + 2 annular planar)
    assert_eq!(output.face_descriptors.len(), 4, "pipe should have 4 faces");

    check_topology(&output);
}

// ---------------------------------------------------------------------------
// Sphere tests
// ---------------------------------------------------------------------------

#[test]
fn simple_sphere_topology() {
    let stl = format!("{}/tests/ccad/generated/simple_sphere.stl", manifest_dir());
    let config = config_for_stl(&stl);
    let output = run_stage3(&config);

    // Sphere: 1 surface, 0 edges, 0 vertices
    assert_eq!(output.face_descriptors.len(), 1, "sphere should have 1 face");
    assert_eq!(output.edges.len(), 0, "sphere should have 0 edges");
    assert_eq!(output.vertices.len(), 0, "sphere should have 0 vertices");

    check_topology(&output);
}

#[test]
fn hemisphere_topology() {
    let stl = format!("{}/tests/ccad/generated/hemisphere.stl", manifest_dir());
    let config = config_for_stl(&stl);
    let output = run_stage3(&config);

    // Hemisphere: 2 surfaces (sphere + plane), 1 edge, 0 vertices
    assert_eq!(output.face_descriptors.len(), 2, "hemisphere should have 2 faces");
    assert_eq!(output.edges.len(), 1, "hemisphere should have 1 edge");
    assert_eq!(output.vertices.len(), 0, "hemisphere should have 0 vertices");

    check_topology(&output);
}

#[test]
fn ball_on_cylinder_topology() {
    let stl = format!("{}/tests/ccad/generated/ball_on_cylinder.stl", manifest_dir());
    let config = config_for_stl(&stl);
    let output = run_stage3(&config);

    // Ball on cylinder: 3 surfaces (sphere + cylinder + plane), 2 edges, 0 vertices
    assert_eq!(output.face_descriptors.len(), 3, "ball_on_cylinder should have 3 faces");
    assert_eq!(output.edges.len(), 2, "ball_on_cylinder should have 2 edges");
    assert_eq!(output.vertices.len(), 0, "ball_on_cylinder should have 0 vertices");

    check_topology(&output);
}

#[test]
fn spherical_pocket_topology() {
    let stl = format!("{}/tests/ccad/generated/spherical_pocket.stl", manifest_dir());
    let config = config_for_stl(&stl);
    let output = run_stage3(&config);

    // Spherical pocket: 7 surfaces (6 planar + 1 spherical)
    assert_eq!(output.face_descriptors.len(), 7, "spherical_pocket should have 7 faces");

    let sphere_count = output.face_descriptors.iter()
        .filter(|fd| {
            matches!(
                output.stage2.selected_surfaces[fd.selected_surface_idx],
                stage2::SelectedSurface::Spherical(_)
            )
        })
        .count();
    assert_eq!(sphere_count, 1, "should have 1 spherical face");

    check_topology(&output);
}

// ---------------------------------------------------------------------------
// Chamfered cube (complex planar test)
// ---------------------------------------------------------------------------

#[test]
fn chamfered_cube_topology() {
    let stl = format!("{}/tests/onshape/chamfered_cube_10_c1_medium.stl", manifest_dir());
    let config = config_for_stl(&stl);
    let output = run_stage3(&config);

    // Chamfered cube: 26 faces, 48 edges, 24 vertices
    assert_eq!(output.face_descriptors.len(), 26, "chamfered cube should have 26 faces");
    assert_eq!(output.edges.len(), 48, "chamfered cube should have 48 edges");
    assert_eq!(output.vertices.len(), 24, "chamfered cube should have 24 vertices");

    check_euler_formula(&output);
    check_topology(&output);
}

// ---------------------------------------------------------------------------
// Compare tests: ensure stage 3.1 passes --compare validation
// (These verify that the OCCT surfaces correctly represent the geometry)
// ---------------------------------------------------------------------------

#[test]
fn ccad_cube_stage31_compare() {
    let dir = manifest_dir();
    let config = config_for_compare(
        &format!("{dir}/tests/ccad/generated/cube.stl"),
        &format!("{dir}/tests/ccad/generated/cube.step"),
    );
    run_stage3(&config);
}

#[test]
fn ccad_simple_cylinder_stage31_compare() {
    let dir = manifest_dir();
    let config = config_for_compare(
        &format!("{dir}/tests/ccad/generated/simple_cylinder.stl"),
        &format!("{dir}/tests/ccad/generated/simple_cylinder.step"),
    );
    run_stage3(&config);
}

#[test]
fn ccad_hemisphere_stage31_compare() {
    let dir = manifest_dir();
    let config = config_for_compare(
        &format!("{dir}/tests/ccad/generated/hemisphere.stl"),
        &format!("{dir}/tests/ccad/generated/hemisphere.step"),
    );
    run_stage3(&config);
}

#[test]
fn ccad_ball_on_cylinder_stage31_compare() {
    let dir = manifest_dir();
    let config = config_for_compare(
        &format!("{dir}/tests/ccad/generated/ball_on_cylinder.stl"),
        &format!("{dir}/tests/ccad/generated/ball_on_cylinder.step"),
    );
    run_stage3(&config);
}
