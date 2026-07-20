//! Integration tests for the stage-2 rewrite frontier: edge classification
//! (2.1) and region extraction (2.2) — the fit-free partition.
//!
//! Expectations are verified against the reference STEP structure of known
//! models (cube = 6 faces, block_with_hole = 7, etc.). Models that need
//! threshold calibration (dense tangent-fillet models, coarse fillets) are
//! intentionally absent: they are the calibration corpus of Rewrite Roadmap
//! step 1, tracked in DEVELOPMENT_PLAN.md.

use brepper::config;
use brepper::stage2::{EdgeClass, Stage2Output};
use brepper::{stage1, stage2};

fn manifest_dir() -> String {
    env!("CARGO_MANIFEST_DIR").to_string()
}

/// Run stages 1 and 2 (through the 2.2 frontier) on a test STL.
fn run_segmentation(stl_rel: &str) -> Stage2Output {
    let stl = format!("{}/{stl_rel}", manifest_dir());
    let config = config::parse_config_from(["brepper", &stl, "-q"]).unwrap();
    let mesh = stage1::stage1(&config).expect("stage1 should pass");
    stage2::stage2(&config, mesh, None).expect("stage2 should pass")
}

/// Assert the post-2.2 invariants (DEVELOPMENT_PLAN §2.4): total coverage
/// (every face in exactly one region) and per-region mesh connectivity.
fn assert_partition_invariants(output: &Stage2Output) {
    let seg = &output.segmentation;
    let n_faces = output.mesh.faces.len();

    // Total coverage: every face has a region, and region membership lists
    // cover every face exactly once.
    assert_eq!(seg.face_region.len(), n_faces);
    let mut membership = vec![0usize; n_faces];
    for (rid, region) in seg.regions.iter().enumerate() {
        for &fi in &region.faces {
            assert_eq!(seg.face_region[fi], rid, "face_region/region mismatch");
            membership[fi] += 1;
        }
    }
    assert!(
        membership.iter().all(|&m| m == 1),
        "every face must belong to exactly one region"
    );

    // Connectivity: each region's faces form one component over mesh edges.
    for (rid, region) in seg.regions.iter().enumerate() {
        let mut seen = vec![false; n_faces];
        let mut stack = vec![region.faces[0]];
        seen[region.faces[0]] = true;
        while let Some(fi) = stack.pop() {
            let vc = output.mesh.faces[fi].vertex_count as usize;
            for &ni in &output.mesh.faces[fi].neighbors[..vc] {
                if ni >= 0 && !seen[ni as usize] && seg.face_region[ni as usize] == rid {
                    seen[ni as usize] = true;
                    stack.push(ni as usize);
                }
            }
        }
        assert!(
            region.faces.iter().all(|&fi| seen[fi]),
            "region {rid} is not mesh-connected ({} faces)",
            region.faces.len()
        );
    }
}

#[test]
fn segmentation_cube() {
    let out = run_segmentation("tests/ccad/generated/cube.stl");
    let seg = &out.segmentation;
    assert_eq!(seg.regions.len(), 6, "cube = 6 planar regions");
    assert_eq!(seg.edge_classes.sharp_count, 12);
    assert_eq!(seg.edge_classes.jump_count, 0);
    assert_eq!(seg.edge_classes.interior_count, 0);
    assert_partition_invariants(&out);
}

#[test]
fn segmentation_block_with_hole() {
    let out = run_segmentation("tests/ccad/generated/block_with_hole.stl");
    let seg = &out.segmentation;
    assert_eq!(seg.regions.len(), 7, "block_with_hole = 6 planar + 1 cylindrical");
    assert_eq!(seg.edge_classes.sharp_count, 70);
    assert_eq!(seg.edge_classes.jump_count, 0);
    assert_partition_invariants(&out);
}

#[test]
fn segmentation_simple_cylinder() {
    let out = run_segmentation("tests/ccad/generated/simple_cylinder.stl");
    let seg = &out.segmentation;
    assert_eq!(seg.regions.len(), 3, "cylinder = wall + 2 caps");
    assert_eq!(seg.edge_classes.jump_count, 0);
    assert_partition_invariants(&out);
}

#[test]
fn segmentation_simple_sphere() {
    let out = run_segmentation("tests/ccad/generated/simple_sphere.stl");
    let seg = &out.segmentation;
    assert_eq!(seg.regions.len(), 1, "sphere = 1 region, no boundaries");
    assert_eq!(seg.edge_classes.sharp_count, 0);
    assert_eq!(seg.edge_classes.jump_count, 0);
    assert_partition_invariants(&out);
}

#[test]
fn segmentation_hemisphere() {
    let out = run_segmentation("tests/ccad/generated/hemisphere.stl");
    assert_eq!(out.segmentation.regions.len(), 2, "hemisphere = dome + base");
    assert_partition_invariants(&out);
}

#[test]
fn segmentation_ball_on_cylinder() {
    let out = run_segmentation("tests/ccad/generated/ball_on_cylinder.stl");
    assert_eq!(
        out.segmentation.regions.len(),
        3,
        "ball_on_cylinder = sphere + cylinder + base"
    );
    assert_partition_invariants(&out);
}

#[test]
fn segmentation_pill_fine_tangent_junction() {
    // The capsule's hemisphere↔cylinder junction is G1-smooth (no sharp
    // edges) and the sphere and cylinder share σ_max = 1/r — it is caught
    // purely by the two-sided curvature descriptor (σ_min: 1/r vs 0).
    let out = run_segmentation("tests/onshape/pill_fine.stl");
    let seg = &out.segmentation;
    assert_eq!(seg.regions.len(), 3, "pill = 2 hemisphere caps + cylinder wall");
    assert_eq!(seg.edge_classes.sharp_count, 0);
    assert!(
        seg.edge_classes.jump_count > 0,
        "tangent junction must produce curvature-jump edges"
    );
    // Every boundary edge is jump-classified (nothing sharp on a capsule).
    for (fi, classes) in seg.edge_classes.class.iter().enumerate() {
        let vc = out.mesh.faces[fi].vertex_count as usize;
        for &c in &classes[..vc] {
            assert_ne!(c, EdgeClass::Sharp);
        }
    }
    assert_partition_invariants(&out);
}

#[test]
fn segmentation_pill_coarse_tangent_junction() {
    // Same capsule at coarse tessellation: the 3 true surfaces (2
    // hemisphere caps + cylinder wall) must be clean regions; any residual
    // regions are tiny junction-circle singletons (faces on the tangent
    // seam, whose seam-circle neighborhood is planar — repair-loop input).
    let out = run_segmentation("tests/onshape/pill_coarse.stl");
    let seg = &out.segmentation;
    let big = seg.regions.iter().filter(|r| r.faces.len() >= 100).count();
    assert_eq!(big, 3, "pill = 3 substantial regions, got {big}");
    let max_small = seg
        .regions
        .iter()
        .filter(|r| r.faces.len() < 100)
        .map(|r| r.faces.len())
        .max()
        .unwrap_or(0);
    assert!(max_small <= 2, "junction confetti should be <=2 faces, got {max_small}");
    assert_eq!(seg.edge_classes.sharp_count, 0);
    assert_partition_invariants(&out);
}

#[test]
fn segmentation_pipe_elbow_torus_bend() {
    // The pipe elbow's torus bend is ONE region: smooth interior edges and
    // a continuous curvature descriptor along the tube, so no boundary
    // fires inside the bend. (The abandoned BFS pipeline could only
    // decompose it into ~290 cylindrical strips.) The cylinder walls
    // currently merge with the bend through the similar-radii outer
    // junction (2.4 split-repair territory); the 2 end caps are separate.
    let out = run_segmentation("tests/onshape/pipe_elbow_10_medium.stl");
    let seg = &out.segmentation;
    let huge = seg.regions.iter().filter(|r| r.faces.len() >= 1000).count();
    assert_eq!(huge, 1, "the bend (+ merged walls) is one region, got {huge}");
    let caps = seg
        .regions
        .iter()
        .filter(|r| (20..1000).contains(&r.faces.len()))
        .count();
    assert_eq!(caps, 2, "two end-cap regions, got {caps}");
    let max_small = seg
        .regions
        .iter()
        .filter(|r| r.faces.len() < 20)
        .map(|r| r.faces.len())
        .max()
        .unwrap_or(0);
    assert!(max_small <= 2, "junction confetti should be <=2 faces, got {max_small}");
    assert_partition_invariants(&out);
}

#[test]
fn segmentation_chamfered_cube() {
    // 26 planar faces (6 faces + 12 edge chamfers + 8 corner chamfers),
    // all sharp-bounded.
    let out = run_segmentation("tests/onshape/chamfered_cube_10_c1_medium.stl");
    let seg = &out.segmentation;
    assert_eq!(seg.regions.len(), 26, "chamfered cube = 26 planar regions");
    assert_eq!(seg.edge_classes.jump_count, 0);
    assert_partition_invariants(&out);
}

#[test]
fn segmentation_pipe_and_two_holes() {
    let out = run_segmentation("tests/ccad/generated/pipe.stl");
    assert_eq!(
        out.segmentation.regions.len(),
        4,
        "pipe = 2 planar annulus caps + 2 cylinders"
    );
    assert_partition_invariants(&out);

    let out = run_segmentation("tests/ccad/generated/two_holes.stl");
    assert_eq!(out.segmentation.regions.len(), 8, "two_holes = 6 planar + 2 cylindrical");
    assert_partition_invariants(&out);
}

#[test]
fn segmentation_spherical_pocket() {
    let out = run_segmentation("tests/ccad/generated/spherical_pocket.stl");
    assert_eq!(
        out.segmentation.regions.len(),
        7,
        "spherical_pocket = 6 planar + 1 spherical"
    );
    assert_partition_invariants(&out);
}

#[test]
fn segmentation_rounded_cube_fine_true_regions_present() {
    // The 20 curved true surfaces must appear as substantial regions: 8
    // corner spheres (~1073 faces) + 12 edge fillets (~160). The 6 planes
    // and junction-zone slivers at the triple points are small (plane
    // consolidation and confetti cleanup are the 2.4 repair loop's job).
    let out = run_segmentation("tests/onshape/rounded_cube_10_r2_fine.stl");
    let seg = &out.segmentation;
    assert_eq!(seg.edge_classes.sharp_count, 0, "rounded cube is fully tangent");
    let big = seg.regions.iter().filter(|r| r.faces.len() >= 1000).count();
    let medium = seg
        .regions
        .iter()
        .filter(|r| (100..1000).contains(&r.faces.len()))
        .count();
    assert_eq!(big, 8, "expected 8 corner-sphere regions, got {big}");
    assert_eq!(medium, 12, "expected 12 edge-fillet regions, got {medium}");
    // Planes and confetti: every remaining region is tiny.
    let max_small = seg
        .regions
        .iter()
        .filter(|r| r.faces.len() < 100)
        .map(|r| r.faces.len())
        .max()
        .unwrap_or(0);
    assert!(max_small <= 15, "small regions should be <=15 faces, got {max_small}");
    assert_partition_invariants(&out);
}

#[test]
fn segmentation_deterministic() {
    let a = run_segmentation("tests/ccad/generated/block_with_hole.stl");
    let b = run_segmentation("tests/ccad/generated/block_with_hole.stl");
    assert_eq!(a.segmentation.face_region, b.segmentation.face_region);
    assert_eq!(a.segmentation.regions.len(), b.segmentation.regions.len());
}
