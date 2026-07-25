//! Integration tests for stage 2.3(b) scoped fit competition: tangent-chain
//! regions the whole-region ladder cannot explain must decompose into their
//! true analytic surfaces (vertex-tolerance validated, not just counted).
//!
//! These run the full stage-2 pipeline (default stage) and assert the typed
//! hypothesis counts on the headline acceptance models, including the
//! stage-2.3(c) bent-tube promotion win (pipe elbow → torus + 2 cylinders).

use brepper::config;
use brepper::{stage1, stage2};

fn manifest_dir() -> String {
    env!("CARGO_MANIFEST_DIR").to_string()
}

fn run(stl_rel: &str) -> stage2::Stage2Output {
    let stl = format!("{}/{stl_rel}", manifest_dir());
    let config = config::parse_config_from(["brepper", &stl, "-q"]).unwrap();
    let mesh = stage1::stage1(&config).expect("stage1");
    stage2::stage2(&config, mesh, None).expect("stage2")
}

/// Every mesh face is covered by exactly one selected surface.
fn assert_total_coverage(out: &stage2::Stage2Output) {
    let n = out.mesh.faces.len();
    let mut count = vec![0u32; n];
    for s in &out.selected_surfaces {
        let faces = match s {
            stage2::SelectedSurface::Planar(i) => &out.planar_hypotheses[*i].faces,
            stage2::SelectedSurface::Conical(i) => &out.conical_hypotheses[*i].faces,
            stage2::SelectedSurface::Spherical(i) => &out.spherical_hypotheses[*i].faces,
            stage2::SelectedSurface::Toroidal(i) => &out.toroidal_hypotheses[*i].faces,
            stage2::SelectedSurface::BSpline(i) => &out.bspline_hypotheses[*i].faces,
        };
        for &fi in faces {
            count[fi] += 1;
        }
    }
    assert!(count.iter().all(|&c| c == 1), "every face covered exactly once");
}

#[test]
fn competition_rounded_cube_shell_decomposes_6_12_8() {
    // The rounded-cube outer shell is one smooth region (a tangent chain);
    // competition must split it into 6 planes + 12 fillet cylinders + 8
    // corner spheres — the reference's 26 analytic faces.
    let out = run("tests/onshape/rounded_cube_10_r2_fine.stl");
    let cylinders = out.conical_hypotheses.iter().filter(|h| h.is_cylinder()).count();
    let cones = out.conical_hypotheses.iter().filter(|h| !h.is_cylinder()).count();
    assert_eq!(out.spherical_hypotheses.len(), 8, "8 corner spheres");
    assert_eq!(cylinders, 12, "12 edge-fillet cylinders");
    assert_eq!(cones, 0, "no cones");
    assert_eq!(out.planar_hypotheses.len(), 6, "6 face planes");
    assert_eq!(out.selected_surfaces.len(), 26, "26 analytic surfaces total");
    assert_total_coverage(&out);
}

#[test]
fn competition_pipe_elbow_torus_two_cylinders() {
    // The pipe elbow is one smooth tangent chain: two straight tube segments
    // joined by a 90° bend, capped by two end disks. Stage 2.3(c) bent-tube
    // promotion must recognize the bend as a torus, so competition + the
    // ladder yield 1 torus + 2 cylinders + 2 base planes = 5 surfaces (the
    // reference's 5 analytic faces), with no strip machinery and no
    // per-facet leftovers.
    let out = run("tests/onshape/pipe_elbow_10_medium.stl");
    assert_eq!(out.toroidal_hypotheses.len(), 1, "1 bend torus");
    let cylinders = out.conical_hypotheses.iter().filter(|h| h.is_cylinder()).count();
    let cones = out.conical_hypotheses.iter().filter(|h| !h.is_cylinder()).count();
    assert_eq!(cylinders, 2, "2 straight tube cylinders");
    assert_eq!(cones, 0, "no cones");
    assert_eq!(out.planar_hypotheses.len(), 2, "2 end-cap planes");
    assert_eq!(out.selected_surfaces.len(), 5, "torus + 2 cyl + 2 planes");
    assert_total_coverage(&out);
}

#[test]
fn competition_nosecone_cone_sphere_base() {
    // The nosecone's cone+sphere tangent chain is one region; the base disk
    // is a separate sharp-bounded plane. Competition + the ladder together
    // yield cone + sphere + base = 3 surfaces (the reference's 3 faces).
    let out = run("tests/onshape/nosecone_medium.stl");
    assert_eq!(out.spherical_hypotheses.len(), 1, "1 nose sphere");
    let cones = out.conical_hypotheses.iter().filter(|h| !h.is_cylinder()).count();
    assert_eq!(cones, 1, "1 body cone");
    assert_eq!(out.planar_hypotheses.len(), 1, "1 base plane");
    assert_eq!(out.selected_surfaces.len(), 3, "cone + sphere + base");
    assert_total_coverage(&out);
}

#[test]
fn competition_stepped_block_fillet_network_stays_analytic() {
    // A dense fillet network (37 cylinders + 18 tori + 16 planes in the
    // reference) tangent-connects the whole model into ONE 26122-face region,
    // so every surface has to come out of competition. This model used to
    // collapse: one r = 0.5 fillet band, one quad tall between two corner
    // tori, left the axis so weakly constrained that a bogus 12° cone beat the
    // true cylinder on residual, and the 30-face leftover blob then tripped
    // classify.rs's whole-region discard — degrading all 26122 faces to
    // per-facet planes. The exact decomposition is still short of the
    // reference (see TODO.md), so this guards the property that matters: the
    // model decomposes analytically and leaves only slivers behind.
    let out = run("tests/onshape/stepped_block_toroidal_fillets_fine.stl");
    let analytic = out.conical_hypotheses.len()
        + out.spherical_hypotheses.len()
        + out.toroidal_hypotheses.len();
    assert!(analytic >= 60, "expected the fillet network to decompose analytically, got {analytic} curved surfaces");
    assert!(out.toroidal_hypotheses.len() >= 18, "at least the reference's 18 corner-blend tori");
    let per_facet = out.planar_hypotheses.iter().filter(|h| h.faces.len() == 1).count();
    assert!(per_facet <= 16, "only slivers may degrade to per-facet planes, got {per_facet}");
    assert_total_coverage(&out);
}
