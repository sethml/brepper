//! Integration tests for the stage-2 rewrite frontier: sharp-edge
//! classification (2.1(a)) and smooth-region extraction (2.2) — the
//! sharp-only partition.
//!
//! Region counts here are **smooth-region** counts, not final surface
//! counts: a smooth region is a maximal set of faces reachable without
//! crossing a sharp edge, so tangent chains merge into one region (the pill
//! capsule is ONE region, the rounded-cube shell is ONE region, the pipe
//! elbow's tube is ONE region). Decomposing those chains into their true
//! surfaces is stage 2.3's fit competition, tested elsewhere. Sharp-edged
//! analytic models already segment exactly here (cube = 6, block_with_hole
//! = 7, …).

use brepper::config;
use brepper::stage2::Stage2Output;
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
    assert_eq!(seg.edge_classes.interior_count, 0);
    assert_partition_invariants(&out);
}

#[test]
fn segmentation_block_with_hole() {
    let out = run_segmentation("tests/ccad/generated/block_with_hole.stl");
    let seg = &out.segmentation;
    assert_eq!(seg.regions.len(), 7, "block_with_hole = 6 planar + 1 cylindrical");
    assert_eq!(seg.edge_classes.sharp_count, 70);
    assert_partition_invariants(&out);
}

#[test]
fn segmentation_simple_cylinder() {
    let out = run_segmentation("tests/ccad/generated/simple_cylinder.stl");
    let seg = &out.segmentation;
    assert_eq!(seg.regions.len(), 3, "cylinder = wall + 2 caps");
    assert_partition_invariants(&out);
}

#[test]
fn segmentation_simple_sphere() {
    let out = run_segmentation("tests/ccad/generated/simple_sphere.stl");
    let seg = &out.segmentation;
    assert_eq!(seg.regions.len(), 1, "sphere = 1 region, no boundaries");
    assert_eq!(seg.edge_classes.sharp_count, 0);
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
fn segmentation_pill_fine_one_smooth_region() {
    // The capsule's hemisphere↔cylinder junctions are G1-smooth (no sharp
    // edges), so the entire outer surface is ONE smooth region. Splitting it
    // into the 2 hemisphere caps + cylinder wall is stage 2.3's job.
    let out = run_segmentation("tests/onshape/pill_fine.stl");
    let seg = &out.segmentation;
    assert_eq!(seg.regions.len(), 1, "capsule surface = one smooth region");
    assert_eq!(seg.edge_classes.sharp_count, 0, "a capsule has no sharp edges");
    assert!(seg.edge_classes.interior_count > 0);
    assert_partition_invariants(&out);
}

#[test]
fn segmentation_pill_coarse_one_smooth_region() {
    // Same capsule at coarse tessellation: still one tangent chain, one
    // smooth region.
    let out = run_segmentation("tests/onshape/pill_coarse.stl");
    let seg = &out.segmentation;
    assert_eq!(seg.regions.len(), 1, "capsule surface = one smooth region");
    assert_eq!(seg.edge_classes.sharp_count, 0);
    assert_partition_invariants(&out);
}

#[test]
fn segmentation_pipe_elbow_tube_one_region() {
    // The pipe elbow's tube surface — the torus bend tangent to both
    // straight cylinder walls — is ONE smooth region (no sharp edge along
    // the whole sweep). The 2 flat end caps are sharp-bounded and separate.
    // Decomposing the tube into torus + 2 cylinders is stage 2.3's job.
    let out = run_segmentation("tests/onshape/pipe_elbow_10_medium.stl");
    let seg = &out.segmentation;
    assert_eq!(seg.regions.len(), 3, "tube + 2 end caps");
    let huge = seg.regions.iter().filter(|r| r.faces.len() >= 1000).count();
    assert_eq!(huge, 1, "the whole tube (bend + walls) is one region, got {huge}");
    let caps = seg
        .regions
        .iter()
        .filter(|r| (20..1000).contains(&r.faces.len()))
        .count();
    assert_eq!(caps, 2, "two end-cap regions, got {caps}");
    assert_partition_invariants(&out);
}

#[test]
fn segmentation_chamfered_cube() {
    // 26 planar faces (6 faces + 12 edge chamfers + 8 corner chamfers),
    // all sharp-bounded.
    let out = run_segmentation("tests/onshape/chamfered_cube_10_c1_medium.stl");
    let seg = &out.segmentation;
    assert_eq!(seg.regions.len(), 26, "chamfered cube = 26 planar regions");
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
fn segmentation_rounded_cube_fine_one_smooth_region() {
    // The rounded cube's entire outer shell is one tangent chain: 6 planes,
    // 12 edge-fillet cylinders, and 8 corner spheres, every junction
    // G1-smooth. So it is ONE smooth region. Decomposing it into 6+12+8 is
    // stage 2.3's fit competition.
    let out = run_segmentation("tests/onshape/rounded_cube_10_r2_fine.stl");
    let seg = &out.segmentation;
    assert_eq!(seg.edge_classes.sharp_count, 0, "rounded cube is fully tangent");
    assert_eq!(seg.regions.len(), 1, "the whole outer shell is one smooth region");
    assert_partition_invariants(&out);
}

#[test]
fn segmentation_deterministic() {
    let a = run_segmentation("tests/ccad/generated/block_with_hole.stl");
    let b = run_segmentation("tests/ccad/generated/block_with_hole.stl");
    assert_eq!(a.segmentation.face_region, b.segmentation.face_region);
    assert_eq!(a.segmentation.regions.len(), b.segmentation.regions.len());
}
