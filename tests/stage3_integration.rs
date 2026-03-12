//! Integration tests for Stage 3: surface reconstruction.
//!
//! Tests that stage 3.1 (OCCT surface creation and adjacency graph construction)
//! produces correct topology for known models.

use brepper::config;
use brepper::{stage1, stage2, stage3};
use opencascade_sys::b_rep_check;

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

fn manifest_dir() -> String {
    env!("CARGO_MANIFEST_DIR").to_string()
}

fn config_for_stl(stl_path: &str) -> config::Config {
    config::parse_config_from(["brepper", stl_path, "--stage=3.3", "-q"]).unwrap()
}

fn config_for_compare(stl_path: &str, step_path: &str) -> config::Config {
    let mut config =
        config::parse_config_from(["brepper", stl_path, "--compare", step_path, "--stage=3.3", "-q"])
            .unwrap();
    config.load_compare_step().unwrap();
    config
}

/// Run stages 1, 2, and 3 (through current max substage), returning the Stage3Output.
fn run_stage3(config: &config::Config) -> stage3::Stage3Output {
    let mesh = stage1::stage1(config).expect("stage1 should pass");
    let surfaces = stage2::stage2(config, mesh, None).expect("stage2 should pass");
    stage3::stage3(config, surfaces, None).expect("stage3 should pass")
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

#[test]
fn ccad_wedge_stage31_compare() {
    let dir = manifest_dir();
    let config = config_for_compare(
        &format!("{dir}/tests/ccad/generated/wedge.stl"),
        &format!("{dir}/tests/ccad/generated/wedge.step"),
    );
    run_stage3(&config);
}

#[test]
fn ccad_block_with_hole_stage31_compare() {
    let dir = manifest_dir();
    let config = config_for_compare(
        &format!("{dir}/tests/ccad/generated/block_with_hole.stl"),
        &format!("{dir}/tests/ccad/generated/block_with_hole.step"),
    );
    run_stage3(&config);
}

#[test]
fn ccad_pipe_stage31_compare() {
    let dir = manifest_dir();
    let config = config_for_compare(
        &format!("{dir}/tests/ccad/generated/pipe.stl"),
        &format!("{dir}/tests/ccad/generated/pipe.step"),
    );
    run_stage3(&config);
}

#[test]
fn ccad_simple_sphere_stage31_compare() {
    let dir = manifest_dir();
    let config = config_for_compare(
        &format!("{dir}/tests/ccad/generated/simple_sphere.stl"),
        &format!("{dir}/tests/ccad/generated/simple_sphere.step"),
    );
    run_stage3(&config);
}

#[test]
fn ccad_spherical_pocket_stage31_compare() {
    let dir = manifest_dir();
    let config = config_for_compare(
        &format!("{dir}/tests/ccad/generated/spherical_pocket.stl"),
        &format!("{dir}/tests/ccad/generated/spherical_pocket.step"),
    );
    run_stage3(&config);
}

#[test]
fn onshape_chamfered_cube_stage31_compare() {
    let dir = manifest_dir();
    let config = config_for_compare(
        &format!("{dir}/tests/onshape/chamfered_cube_10_c1_medium.stl"),
        &format!("{dir}/tests/onshape/chamfered_cube_10_c1.step"),
    );
    run_stage3(&config);
}

#[test]
fn manual_cube_stage31_compare() {
    let dir = manifest_dir();
    let config = config_for_compare(
        &format!("{dir}/tests/manual/cube.stl"),
        &format!("{dir}/tests/manual/cube.step"),
    );
    run_stage3(&config);
}

#[test]
fn onshape_part_rounded_cube_stage31_compare() {
    let dir = manifest_dir();
    let config = config_for_compare(
        &format!("{dir}/tests/onshape/part_rounded_cube_10_r2_coarse.stl"),
        &format!("{dir}/tests/onshape/part_rounded_cube_10_r2.step"),
    );
    run_stage3(&config);
}

// ---------------------------------------------------------------------------
// Tangency detection tests (stage 3.2)
// ---------------------------------------------------------------------------

/// For models composed only of planar surfaces meeting at angles, no edges should be tangent.
#[test]
fn cube_no_tangent_edges() {
    let stl = format!("{}/tests/ccad/generated/cube.stl", manifest_dir());
    let config = config_for_stl(&stl);
    let output = run_stage3(&config);
    for (ei, edge) in output.edges.iter().enumerate() {
        assert!(!edge.tangent, "cube edge {ei} should not be tangent");
    }
}

/// Chamfered cube: all surfaces are planar at distinct angles, no tangency.
#[test]
fn chamfered_cube_no_tangent_edges() {
    let stl = format!("{}/tests/onshape/chamfered_cube_10_c1_medium.stl", manifest_dir());
    let config = config_for_stl(&stl);
    let output = run_stage3(&config);
    for (ei, edge) in output.edges.iter().enumerate() {
        assert!(!edge.tangent, "chamfered cube edge {ei} should not be tangent");
    }
}

/// Cylinder with planar end caps: plane-cylinder edges are not tangent.
#[test]
fn cylinder_no_tangent_edges() {
    let stl = format!("{}/tests/ccad/generated/simple_cylinder.stl", manifest_dir());
    let config = config_for_stl(&stl);
    let output = run_stage3(&config);
    for (ei, edge) in output.edges.iter().enumerate() {
        assert!(!edge.tangent, "cylinder edge {ei} should not be tangent");
    }
}

/// Block with hole: plane-cylinder edges are not tangent.
#[test]
fn block_with_hole_no_tangent_edges() {
    let stl = format!("{}/tests/ccad/generated/block_with_hole.stl", manifest_dir());
    let config = config_for_stl(&stl);
    let output = run_stage3(&config);
    for (ei, edge) in output.edges.iter().enumerate() {
        assert!(!edge.tangent, "block_with_hole edge {ei} should not be tangent");
    }
}

/// Hemisphere: sphere-plane edge is not tangent (they meet at 90° around the equator).
#[test]
fn hemisphere_no_tangent_edges() {
    let stl = format!("{}/tests/ccad/generated/hemisphere.stl", manifest_dir());
    let config = config_for_stl(&stl);
    let output = run_stage3(&config);
    for (ei, edge) in output.edges.iter().enumerate() {
        assert!(!edge.tangent, "hemisphere edge {ei} should not be tangent");
    }
}

/// Part rounded cube: 8 tangent edges where cylinder fillets meet planar faces.
#[test]
fn part_rounded_cube_tangent_edges() {
    let stl = format!("{}/tests/onshape/part_rounded_cube_10_r2_coarse.stl", manifest_dir());
    let config = config_for_stl(&stl);
    let output = run_stage3(&config);
    let tangent_count = output.edges.iter().filter(|e| e.tangent).count();
    assert_eq!(tangent_count, 8, "part_rounded_cube should have 8 tangent edges (got {tangent_count})");
    // All edges should have curves computed (including tangent ones)
    check_all_edges_have_curves(&output);
}
/// Full rounded cube: 48 tangent edges on medium mesh (plane-cylinder + sphere-cylinder).
#[test]
fn rounded_cube_tangent_edges() {
    let stl = format!("{}/tests/onshape/rounded_cube_10_r2_medium.stl", manifest_dir());
    let config = config_for_stl(&stl);
    let output = run_stage3(&config);
    let tangent_count = output.edges.iter().filter(|e| e.tangent).count();
    assert_eq!(tangent_count, 48, "rounded_cube should have 48 tangent edges (got {tangent_count})");
    // All edges should have curves computed (including sphere-cylinder tangent ones)
    check_all_edges_have_curves(&output);
}


// ---------------------------------------------------------------------------
// Edge curve tests (stage 3.3)
// ---------------------------------------------------------------------------

/// Verify that all edges have computed 3D curves.
fn check_all_edges_have_curves(output: &stage3::Stage3Output) {
    for (ei, edge) in output.edges.iter().enumerate() {
        assert!(
            edge.curve_3d.is_some(),
            "edge {ei} should have a 3D curve after stage 3.3"
        );
    }
}

/// Cube: all 12 plane-plane intersection edges should produce line curves.
#[test]
fn cube_edge_curves() {
    let stl = format!("{}/tests/ccad/generated/cube.stl", manifest_dir());
    let config = config_for_stl(&stl);
    let output = run_stage3(&config);
    check_all_edges_have_curves(&output);
    assert_eq!(output.edges.len(), 12);
}

/// Cylinder: 2 plane-cylinder intersection edges should produce circle curves.
#[test]
fn cylinder_edge_curves() {
    let stl = format!("{}/tests/ccad/generated/simple_cylinder.stl", manifest_dir());
    let config = config_for_stl(&stl);
    let output = run_stage3(&config);
    check_all_edges_have_curves(&output);
    assert_eq!(output.edges.len(), 2);
}

/// Block with hole: 15 edges (12 planar + 2 plane-cylinder circles + 1 extra).
#[test]
fn block_with_hole_edge_curves() {
    let stl = format!("{}/tests/ccad/generated/block_with_hole.stl", manifest_dir());
    let config = config_for_stl(&stl);
    let output = run_stage3(&config);
    check_all_edges_have_curves(&output);
}

/// Hemisphere: sphere-plane intersection should produce a circle curve.
#[test]
fn hemisphere_edge_curves() {
    let stl = format!("{}/tests/ccad/generated/hemisphere.stl", manifest_dir());
    let config = config_for_stl(&stl);
    let output = run_stage3(&config);
    check_all_edges_have_curves(&output);
    assert_eq!(output.edges.len(), 1);
}

/// Chamfered cube: all 48 plane-plane edges should have curves.
#[test]
fn chamfered_cube_edge_curves() {
    let stl = format!("{}/tests/onshape/chamfered_cube_10_c1_medium.stl", manifest_dir());
    let config = config_for_stl(&stl);
    let output = run_stage3(&config);
    check_all_edges_have_curves(&output);
    assert_eq!(output.edges.len(), 48);
}

/// Spherical pocket: 13 edges mixing plane-plane and plane-sphere intersections.
#[test]
fn spherical_pocket_edge_curves() {
    let stl = format!("{}/tests/ccad/generated/spherical_pocket.stl", manifest_dir());
    let config = config_for_stl(&stl);
    let output = run_stage3(&config);
    check_all_edges_have_curves(&output);
}

/// Ball on cylinder: sphere-cylinder and cylinder-plane intersections.
#[test]
fn ball_on_cylinder_edge_curves() {
    let stl = format!("{}/tests/ccad/generated/ball_on_cylinder.stl", manifest_dir());
    let config = config_for_stl(&stl);
    let output = run_stage3(&config);
    check_all_edges_have_curves(&output);
    assert_eq!(output.edges.len(), 2);
}

// ---------------------------------------------------------------------------
// Face creation tests (stage 3.4)
// ---------------------------------------------------------------------------

fn config_for_stl_34(stl_path: &str) -> config::Config {
    config::parse_config_from(["brepper", stl_path, "--stage=3.4", "-q"]).unwrap()
}

fn config_for_compare_34(stl_path: &str, step_path: &str) -> config::Config {
    let mut config =
        config::parse_config_from(["brepper", stl_path, "--compare", step_path, "--stage=3.4", "-q"])
            .unwrap();
    config.load_compare_step().unwrap();
    config
}

fn run_stage3_34(config: &config::Config) -> stage3::Stage3Output {
    let mesh = stage1::stage1(config).expect("stage1 should pass");
    let surfaces = stage2::stage2(config, mesh, None).expect("stage2 should pass");
    stage3::stage3(config, surfaces, None).expect("stage3 should pass")
}

/// Verify that every face descriptor has a corresponding OCCT face.
fn check_all_faces_created(output: &stage3::Stage3Output) {
    assert_eq!(
        output.make_faces.len(),
        output.face_descriptors.len(),
        "make_faces count should match face_descriptors count"
    );
}

// --- Cube ---

#[test]
fn cube_face_creation() {
    let stl = format!("{}/tests/ccad/generated/cube.stl", manifest_dir());
    let config = config_for_stl_34(&stl);
    let output = run_stage3_34(&config);
    check_all_faces_created(&output);
    assert_eq!(output.make_faces.len(), 6);
}

// --- Wedge ---

#[test]
fn wedge_face_creation() {
    let stl = format!("{}/tests/ccad/generated/wedge.stl", manifest_dir());
    let config = config_for_stl_34(&stl);
    let output = run_stage3_34(&config);
    check_all_faces_created(&output);
    assert_eq!(output.make_faces.len(), 6);
}

// --- Cylinder ---

#[test]
fn simple_cylinder_face_creation() {
    let stl = format!("{}/tests/ccad/generated/simple_cylinder.stl", manifest_dir());
    let config = config_for_stl_34(&stl);
    let output = run_stage3_34(&config);
    check_all_faces_created(&output);
    assert_eq!(output.make_faces.len(), 3);
}

// --- Block with hole ---

#[test]
fn block_with_hole_face_creation() {
    let stl = format!("{}/tests/ccad/generated/block_with_hole.stl", manifest_dir());
    let config = config_for_stl_34(&stl);
    let output = run_stage3_34(&config);
    check_all_faces_created(&output);
    assert_eq!(output.make_faces.len(), 7);
}

// --- Pipe ---

#[test]
fn pipe_face_creation() {
    let stl = format!("{}/tests/ccad/generated/pipe.stl", manifest_dir());
    let config = config_for_stl_34(&stl);
    let output = run_stage3_34(&config);
    check_all_faces_created(&output);
    assert_eq!(output.make_faces.len(), 4);
}

// --- Sphere ---

#[test]
fn simple_sphere_face_creation() {
    let stl = format!("{}/tests/ccad/generated/simple_sphere.stl", manifest_dir());
    let config = config_for_stl_34(&stl);
    let output = run_stage3_34(&config);
    check_all_faces_created(&output);
    assert_eq!(output.make_faces.len(), 1);
}

// --- Hemisphere ---

#[test]
fn hemisphere_face_creation() {
    let stl = format!("{}/tests/ccad/generated/hemisphere.stl", manifest_dir());
    let config = config_for_stl_34(&stl);
    let output = run_stage3_34(&config);
    check_all_faces_created(&output);
    assert_eq!(output.make_faces.len(), 2);
}

// --- Ball on cylinder ---

#[test]
fn ball_on_cylinder_face_creation() {
    let stl = format!("{}/tests/ccad/generated/ball_on_cylinder.stl", manifest_dir());
    let config = config_for_stl_34(&stl);
    let output = run_stage3_34(&config);
    check_all_faces_created(&output);
    assert_eq!(output.make_faces.len(), 3);
}

// --- Spherical pocket ---

#[test]
fn spherical_pocket_face_creation() {
    let stl = format!("{}/tests/ccad/generated/spherical_pocket.stl", manifest_dir());
    let config = config_for_stl_34(&stl);
    let output = run_stage3_34(&config);
    check_all_faces_created(&output);
    assert_eq!(output.make_faces.len(), 7);
}

// --- Chamfered cube ---

#[test]
fn chamfered_cube_face_creation() {
    let stl = format!("{}/tests/onshape/chamfered_cube_10_c1_medium.stl", manifest_dir());
    let config = config_for_stl_34(&stl);
    let output = run_stage3_34(&config);
    check_all_faces_created(&output);
    assert_eq!(output.make_faces.len(), 26);
}

// ---------------------------------------------------------------------------
// Stage 3.4 compare tests
// ---------------------------------------------------------------------------

#[test]
fn ccad_cube_stage34_compare() {
    let dir = manifest_dir();
    let config = config_for_compare_34(
        &format!("{dir}/tests/ccad/generated/cube.stl"),
        &format!("{dir}/tests/ccad/generated/cube.step"),
    );
    run_stage3_34(&config);
}

#[test]
fn ccad_wedge_stage34_compare() {
    let dir = manifest_dir();
    let config = config_for_compare_34(
        &format!("{dir}/tests/ccad/generated/wedge.stl"),
        &format!("{dir}/tests/ccad/generated/wedge.step"),
    );
    run_stage3_34(&config);
}

#[test]
fn ccad_simple_cylinder_stage34_compare() {
    let dir = manifest_dir();
    let config = config_for_compare_34(
        &format!("{dir}/tests/ccad/generated/simple_cylinder.stl"),
        &format!("{dir}/tests/ccad/generated/simple_cylinder.step"),
    );
    run_stage3_34(&config);
}

#[test]
fn ccad_block_with_hole_stage34_compare() {
    let dir = manifest_dir();
    let config = config_for_compare_34(
        &format!("{dir}/tests/ccad/generated/block_with_hole.stl"),
        &format!("{dir}/tests/ccad/generated/block_with_hole.step"),
    );
    run_stage3_34(&config);
}

#[test]
fn ccad_pipe_stage34_compare() {
    let dir = manifest_dir();
    let config = config_for_compare_34(
        &format!("{dir}/tests/ccad/generated/pipe.stl"),
        &format!("{dir}/tests/ccad/generated/pipe.step"),
    );
    run_stage3_34(&config);
}

#[test]
fn ccad_simple_sphere_stage34_compare() {
    let dir = manifest_dir();
    let config = config_for_compare_34(
        &format!("{dir}/tests/ccad/generated/simple_sphere.stl"),
        &format!("{dir}/tests/ccad/generated/simple_sphere.step"),
    );
    run_stage3_34(&config);
}

#[test]
fn ccad_hemisphere_stage34_compare() {
    let dir = manifest_dir();
    let config = config_for_compare_34(
        &format!("{dir}/tests/ccad/generated/hemisphere.stl"),
        &format!("{dir}/tests/ccad/generated/hemisphere.step"),
    );
    run_stage3_34(&config);
}

#[test]
fn ccad_ball_on_cylinder_stage34_compare() {
    let dir = manifest_dir();
    let config = config_for_compare_34(
        &format!("{dir}/tests/ccad/generated/ball_on_cylinder.stl"),
        &format!("{dir}/tests/ccad/generated/ball_on_cylinder.step"),
    );
    run_stage3_34(&config);
}

#[test]
fn ccad_spherical_pocket_stage34_compare() {
    let dir = manifest_dir();
    let config = config_for_compare_34(
        &format!("{dir}/tests/ccad/generated/spherical_pocket.stl"),
        &format!("{dir}/tests/ccad/generated/spherical_pocket.step"),
    );
    run_stage3_34(&config);
}

#[test]
fn onshape_chamfered_cube_stage34_compare() {
    let dir = manifest_dir();
    let config = config_for_compare_34(
        &format!("{dir}/tests/onshape/chamfered_cube_10_c1_medium.stl"),
        &format!("{dir}/tests/onshape/chamfered_cube_10_c1.step"),
    );
    run_stage3_34(&config);
}

#[test]
fn manual_cube_stage34_compare() {
    let dir = manifest_dir();
    let config = config_for_compare_34(
        &format!("{dir}/tests/manual/cube.stl"),
        &format!("{dir}/tests/manual/cube.step"),
    );
    run_stage3_34(&config);
}

#[test]
fn onshape_part_rounded_cube_stage34_compare() {
    let dir = manifest_dir();
    let config = config_for_compare_34(
        &format!("{dir}/tests/onshape/part_rounded_cube_10_r2_coarse.stl"),
        &format!("{dir}/tests/onshape/part_rounded_cube_10_r2.step"),
    );
    run_stage3_34(&config);
}

// ===========================================================================
// Stage 3.5 tests — shell construction via BRepBuilderAPI_Sewing
// ===========================================================================

fn config_for_stl_35(stl_path: &str) -> config::Config {
    config::parse_config_from(["brepper", stl_path, "--stage=3.5", "-q"]).unwrap()
}

fn config_for_compare_35(stl_path: &str, step_path: &str) -> config::Config {
    let mut config =
        config::parse_config_from(["brepper", stl_path, "--compare", step_path, "--stage=3.5", "-q"])
            .unwrap();
    config.load_compare_step().unwrap();
    config
}

fn run_stage3_35(config: &config::Config) -> stage3::Stage3Output {
    let mesh = stage1::stage1(config).expect("stage1 should pass");
    let surfaces = stage2::stage2(config, mesh, None).expect("stage2 should pass");
    stage3::stage3(config, surfaces, None).expect("stage3 should pass")
}

/// Verify that at least one shell was constructed.
fn check_shells_constructed(output: &stage3::Stage3Output) {
    assert!(
        !output.shells.is_empty(),
        "at least one shell should be constructed"
    );
}

// --- Shell construction tests ---

#[test]
fn cube_shell_construction() {
    let stl = format!("{}/tests/ccad/generated/cube.stl", manifest_dir());
    let config = config_for_stl_35(&stl);
    let output = run_stage3_35(&config);
    check_shells_constructed(&output);
    assert_eq!(output.shells.len(), 1);
}

#[test]
fn simple_cylinder_shell_construction() {
    let stl = format!("{}/tests/ccad/generated/simple_cylinder.stl", manifest_dir());
    let config = config_for_stl_35(&stl);
    let output = run_stage3_35(&config);
    check_shells_constructed(&output);
    assert_eq!(output.shells.len(), 1);
}

#[test]
fn simple_sphere_shell_construction() {
    let stl = format!("{}/tests/ccad/generated/simple_sphere.stl", manifest_dir());
    let config = config_for_stl_35(&stl);
    let output = run_stage3_35(&config);
    check_shells_constructed(&output);
    assert_eq!(output.shells.len(), 1);
}

// --- Stage 3.5 compare tests ---

#[test]
fn ccad_cube_stage35_compare() {
    let dir = manifest_dir();
    let config = config_for_compare_35(
        &format!("{dir}/tests/ccad/generated/cube.stl"),
        &format!("{dir}/tests/ccad/generated/cube.step"),
    );
    run_stage3_35(&config);
}

#[test]
fn ccad_wedge_stage35_compare() {
    let dir = manifest_dir();
    let config = config_for_compare_35(
        &format!("{dir}/tests/ccad/generated/wedge.stl"),
        &format!("{dir}/tests/ccad/generated/wedge.step"),
    );
    run_stage3_35(&config);
}

#[test]
fn ccad_simple_cylinder_stage35_compare() {
    let dir = manifest_dir();
    let config = config_for_compare_35(
        &format!("{dir}/tests/ccad/generated/simple_cylinder.stl"),
        &format!("{dir}/tests/ccad/generated/simple_cylinder.step"),
    );
    run_stage3_35(&config);
}

#[test]
fn ccad_block_with_hole_stage35_compare() {
    let dir = manifest_dir();
    let config = config_for_compare_35(
        &format!("{dir}/tests/ccad/generated/block_with_hole.stl"),
        &format!("{dir}/tests/ccad/generated/block_with_hole.step"),
    );
    run_stage3_35(&config);
}

#[test]
fn ccad_pipe_stage35_compare() {
    let dir = manifest_dir();
    let config = config_for_compare_35(
        &format!("{dir}/tests/ccad/generated/pipe.stl"),
        &format!("{dir}/tests/ccad/generated/pipe.step"),
    );
    run_stage3_35(&config);
}

#[test]
fn ccad_simple_sphere_stage35_compare() {
    let dir = manifest_dir();
    let config = config_for_compare_35(
        &format!("{dir}/tests/ccad/generated/simple_sphere.stl"),
        &format!("{dir}/tests/ccad/generated/simple_sphere.step"),
    );
    run_stage3_35(&config);
}

#[test]
fn ccad_hemisphere_stage35_compare() {
    let dir = manifest_dir();
    let config = config_for_compare_35(
        &format!("{dir}/tests/ccad/generated/hemisphere.stl"),
        &format!("{dir}/tests/ccad/generated/hemisphere.step"),
    );
    run_stage3_35(&config);
}

#[test]
fn ccad_ball_on_cylinder_stage35_compare() {
    let dir = manifest_dir();
    let config = config_for_compare_35(
        &format!("{dir}/tests/ccad/generated/ball_on_cylinder.stl"),
        &format!("{dir}/tests/ccad/generated/ball_on_cylinder.step"),
    );
    run_stage3_35(&config);
}

#[test]
fn ccad_spherical_pocket_stage35_compare() {
    let dir = manifest_dir();
    let config = config_for_compare_35(
        &format!("{dir}/tests/ccad/generated/spherical_pocket.stl"),
        &format!("{dir}/tests/ccad/generated/spherical_pocket.step"),
    );
    run_stage3_35(&config);
}

#[test]
fn onshape_chamfered_cube_stage35_compare() {
    let dir = manifest_dir();
    let config = config_for_compare_35(
        &format!("{dir}/tests/onshape/chamfered_cube_10_c1_medium.stl"),
        &format!("{dir}/tests/onshape/chamfered_cube_10_c1.step"),
    );
    run_stage3_35(&config);
}

#[test]
fn manual_cube_stage35_compare() {
    let dir = manifest_dir();
    let config = config_for_compare_35(
        &format!("{dir}/tests/manual/cube.stl"),
        &format!("{dir}/tests/manual/cube.step"),
    );
    run_stage3_35(&config);
}

#[test]
fn onshape_part_rounded_cube_stage35_compare() {
    let dir = manifest_dir();
    let config = config_for_compare_35(
        &format!("{dir}/tests/onshape/part_rounded_cube_10_r2_coarse.stl"),
        &format!("{dir}/tests/onshape/part_rounded_cube_10_r2.step"),
    );
    run_stage3_35(&config);
}

// ---------------------------------------------------------------------------
// Stage 3.6 helpers
// ---------------------------------------------------------------------------

fn config_for_stl_36(stl_path: &str) -> config::Config {
    config::parse_config_from(["brepper", stl_path, "--stage=3.6", "-q"]).unwrap()
}

fn config_for_compare_36(stl_path: &str, step_path: &str) -> config::Config {
    let mut config =
        config::parse_config_from(["brepper", stl_path, "--compare", step_path, "--stage=3.6", "-q"])
            .unwrap();
    config.load_compare_step().unwrap();
    config
}

fn run_stage3_36(config: &config::Config) -> stage3::Stage3Output {
    let mesh = stage1::stage1(config).expect("stage1 should pass");
    let surfaces = stage2::stage2(config, mesh, None).expect("stage2 should pass");
    stage3::stage3(config, surfaces, None).expect("stage3 should pass")
}

/// Verify that at least one solid was constructed.
fn check_solids_constructed(output: &stage3::Stage3Output) {
    assert!(
        !output.solids.is_empty(),
        "at least one solid should be constructed"
    );
}

// --- Solid construction tests ---

#[test]
fn cube_solid_construction() {
    let stl = format!("{}/tests/ccad/generated/cube.stl", manifest_dir());
    let config = config_for_stl_36(&stl);
    let output = run_stage3_36(&config);
    check_solids_constructed(&output);
    assert_eq!(output.solids.len(), 1);
}

#[test]
fn simple_cylinder_solid_construction() {
    let stl = format!("{}/tests/ccad/generated/simple_cylinder.stl", manifest_dir());
    let config = config_for_stl_36(&stl);
    let output = run_stage3_36(&config);
    check_solids_constructed(&output);
    assert_eq!(output.solids.len(), 1);
}

#[test]
fn simple_sphere_solid_construction() {
    let stl = format!("{}/tests/ccad/generated/simple_sphere.stl", manifest_dir());
    let config = config_for_stl_36(&stl);
    let output = run_stage3_36(&config);
    check_solids_constructed(&output);
    assert_eq!(output.solids.len(), 1);
}

// --- Stage 3.6 compare tests ---

#[test]
fn ccad_cube_stage36_compare() {
    let dir = manifest_dir();
    let config = config_for_compare_36(
        &format!("{dir}/tests/ccad/generated/cube.stl"),
        &format!("{dir}/tests/ccad/generated/cube.step"),
    );
    run_stage3_36(&config);
}

#[test]
fn ccad_wedge_stage36_compare() {
    let dir = manifest_dir();
    let config = config_for_compare_36(
        &format!("{dir}/tests/ccad/generated/wedge.stl"),
        &format!("{dir}/tests/ccad/generated/wedge.step"),
    );
    run_stage3_36(&config);
}

#[test]
fn ccad_simple_cylinder_stage36_compare() {
    let dir = manifest_dir();
    let config = config_for_compare_36(
        &format!("{dir}/tests/ccad/generated/simple_cylinder.stl"),
        &format!("{dir}/tests/ccad/generated/simple_cylinder.step"),
    );
    run_stage3_36(&config);
}

#[test]
fn ccad_block_with_hole_stage36_compare() {
    let dir = manifest_dir();
    let config = config_for_compare_36(
        &format!("{dir}/tests/ccad/generated/block_with_hole.stl"),
        &format!("{dir}/tests/ccad/generated/block_with_hole.step"),
    );
    run_stage3_36(&config);
}

#[test]
fn ccad_pipe_stage36_compare() {
    let dir = manifest_dir();
    let config = config_for_compare_36(
        &format!("{dir}/tests/ccad/generated/pipe.stl"),
        &format!("{dir}/tests/ccad/generated/pipe.step"),
    );
    run_stage3_36(&config);
}

#[test]
fn ccad_simple_sphere_stage36_compare() {
    let dir = manifest_dir();
    let config = config_for_compare_36(
        &format!("{dir}/tests/ccad/generated/simple_sphere.stl"),
        &format!("{dir}/tests/ccad/generated/simple_sphere.step"),
    );
    run_stage3_36(&config);
}

#[test]
fn ccad_hemisphere_stage36_compare() {
    let dir = manifest_dir();
    let config = config_for_compare_36(
        &format!("{dir}/tests/ccad/generated/hemisphere.stl"),
        &format!("{dir}/tests/ccad/generated/hemisphere.step"),
    );
    run_stage3_36(&config);
}

#[test]
fn ccad_ball_on_cylinder_stage36_compare() {
    let dir = manifest_dir();
    let config = config_for_compare_36(
        &format!("{dir}/tests/ccad/generated/ball_on_cylinder.stl"),
        &format!("{dir}/tests/ccad/generated/ball_on_cylinder.step"),
    );
    run_stage3_36(&config);
}

#[test]
fn ccad_spherical_pocket_stage36_compare() {
    let dir = manifest_dir();
    let config = config_for_compare_36(
        &format!("{dir}/tests/ccad/generated/spherical_pocket.stl"),
        &format!("{dir}/tests/ccad/generated/spherical_pocket.step"),
    );
    run_stage3_36(&config);
}

#[test]
fn onshape_chamfered_cube_stage36_compare() {
    let dir = manifest_dir();
    let config = config_for_compare_36(
        &format!("{dir}/tests/onshape/chamfered_cube_10_c1_medium.stl"),
        &format!("{dir}/tests/onshape/chamfered_cube_10_c1.step"),
    );
    run_stage3_36(&config);
}

#[test]
fn manual_cube_stage36_compare() {
    let dir = manifest_dir();
    let config = config_for_compare_36(
        &format!("{dir}/tests/manual/cube.stl"),
        &format!("{dir}/tests/manual/cube.step"),
    );
    run_stage3_36(&config);
}

#[test]
fn onshape_part_rounded_cube_stage36_compare() {
    let dir = manifest_dir();
    let config = config_for_compare_36(
        &format!("{dir}/tests/onshape/part_rounded_cube_10_r2_coarse.stl"),
        &format!("{dir}/tests/onshape/part_rounded_cube_10_r2.step"),
    );
    run_stage3_36(&config);
}

// ---------------------------------------------------------------------------
// BRepCheck validation tests
// ---------------------------------------------------------------------------

/// Verify that all solids pass BRepCheck validation (no self-intersecting wires, etc.).
fn check_solids_brep_valid(output: &stage3::Stage3Output) {
    for (si, solid) in output.solids.iter().enumerate() {
        let analyzer = b_rep_check::Analyzer::new_shape_bool(
            solid.as_shape(),
            true,
        );
        assert!(
            analyzer.is_valid(),
            "solid {si} failed BRepCheck validation"
        );
    }
}

/// Rounded cube (coarse) with default parameters must produce a BRepCheck-valid solid.
/// This catches sewing-induced wire corruption (e.g., sphere face wires getting
/// cylinder edges inserted by BRepBuilderAPI_Sewing).
#[test]
fn rounded_cube_coarse_brep_check() {
    let stl = format!("{}/tests/onshape/rounded_cube_10_r2_coarse.stl", manifest_dir());
    let config = config_for_stl_36(&stl);
    let output = run_stage3_36(&config);
    check_solids_constructed(&output);
    check_solids_brep_valid(&output);
}

// ---------------------------------------------------------------------------
// Cone topology and pipeline tests
// ---------------------------------------------------------------------------

#[test]
fn simple_cone_topology() {
    let stl = format!("{}/tests/ccad/generated/simple_cone.stl", manifest_dir());
    let config = config_for_stl(&stl);
    let output = run_stage3(&config);

    // Truncated cone: 3 faces (cone + 2 planar caps), 2 edges (circles), 0 vertices
    assert_eq!(output.face_descriptors.len(), 3, "simple_cone should have 3 faces");
    assert_eq!(output.edges.len(), 2, "simple_cone should have 2 edges");
    assert_eq!(output.vertices.len(), 0, "simple_cone should have 0 vertices");

    let cone_faces: Vec<_> = output.face_descriptors.iter().enumerate()
        .filter(|(_, fd)| {
            matches!(
                output.stage2.selected_surfaces[fd.selected_surface_idx],
                stage2::SelectedSurface::Conical(_)
            )
        })
        .collect();
    assert_eq!(cone_faces.len(), 1, "should have 1 conical face");

    check_topology(&output);
}

#[test]
fn simple_cone_stage31_compare() {
    let dir = manifest_dir();
    let config = config_for_compare(
        &format!("{dir}/tests/ccad/generated/simple_cone.stl"),
        &format!("{dir}/tests/ccad/generated/simple_cone.step"),
    );
    run_stage3(&config);
}

#[test]
fn simple_cone_face_creation() {
    let stl = format!("{}/tests/ccad/generated/simple_cone.stl", manifest_dir());
    let config = config_for_stl_34(&stl);
    let output = run_stage3_34(&config);
    check_all_faces_created(&output);
    assert_eq!(output.make_faces.len(), 3);
}

#[test]
fn simple_cone_stage34_compare() {
    let dir = manifest_dir();
    let config = config_for_compare_34(
        &format!("{dir}/tests/ccad/generated/simple_cone.stl"),
        &format!("{dir}/tests/ccad/generated/simple_cone.step"),
    );
    run_stage3_34(&config);
}

#[test]
fn simple_cone_stage35_compare() {
    let dir = manifest_dir();
    let config = config_for_compare_35(
        &format!("{dir}/tests/ccad/generated/simple_cone.stl"),
        &format!("{dir}/tests/ccad/generated/simple_cone.step"),
    );
    run_stage3_35(&config);
}

#[test]
fn simple_cone_stage36_compare() {
    let dir = manifest_dir();
    let config = config_for_compare_36(
        &format!("{dir}/tests/ccad/generated/simple_cone.stl"),
        &format!("{dir}/tests/ccad/generated/simple_cone.step"),
    );
    run_stage3_36(&config);
}

#[test]
fn cone_cylinder_stage31_compare() {
    let dir = manifest_dir();
    let config = config_for_compare(
        &format!("{dir}/tests/ccad/generated/cone_cylinder.stl"),
        &format!("{dir}/tests/ccad/generated/cone_cylinder.step"),
    );
    run_stage3(&config);
}

#[test]
fn cone_cylinder_stage36_compare() {
    let dir = manifest_dir();
    let config = config_for_compare_36(
        &format!("{dir}/tests/ccad/generated/cone_cylinder.stl"),
        &format!("{dir}/tests/ccad/generated/cone_cylinder.step"),
    );
    run_stage3_36(&config);
}

#[test]
fn nosecone_stage31_compare() {
    let dir = manifest_dir();
    let config = config_for_compare(
        &format!("{dir}/tests/ccad/generated/nosecone.stl"),
        &format!("{dir}/tests/ccad/generated/nosecone.step"),
    );
    run_stage3(&config);
}

#[test]
fn nosecone_stage36_compare() {
    let dir = manifest_dir();
    let config = config_for_compare_36(
        &format!("{dir}/tests/ccad/generated/nosecone.stl"),
        &format!("{dir}/tests/ccad/generated/nosecone.step"),
    );
    run_stage3_36(&config);
}