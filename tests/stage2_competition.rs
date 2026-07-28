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
    run_with_config(stl_rel).0
}

fn run_with_config(stl_rel: &str) -> (stage2::Stage2Output, config::Config) {
    let stl = format!("{}/{stl_rel}", manifest_dir());
    let config = config::parse_config_from(["brepper", &stl, "-q"]).unwrap();
    let mesh = stage1::stage1(&config).expect("stage1");
    let out = stage2::stage2(&config, mesh, None).expect("stage2");
    (out, config)
}

/// Count selected surfaces by type: (planar, cylinder, cone, spherical,
/// toroidal, bspline). Per project policy, assertions are on stage 2.6
/// *selected* surfaces — the hypothesis vectors also carry entries that repair
/// emptied out.
fn selected_counts(out: &stage2::Stage2Output) -> [usize; 6] {
    let mut n = [0usize; 6];
    for s in &out.selected_surfaces {
        match s {
            stage2::SelectedSurface::Planar(_) => n[0] += 1,
            stage2::SelectedSurface::Conical(i) => {
                if out.conical_hypotheses[*i].is_cylinder() {
                    n[1] += 1
                } else {
                    n[2] += 1
                }
            }
            stage2::SelectedSurface::Spherical(_) => n[3] += 1,
            stage2::SelectedSurface::Toroidal(_) => n[4] += 1,
            stage2::SelectedSurface::BSpline(_) => n[5] += 1,
        }
    }
    n
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
    let [planar, cylinders, cones, spheres, tori, _] = selected_counts(&out);
    assert_eq!(spheres, 8, "8 corner spheres");
    assert_eq!(cylinders, 12, "12 edge-fillet cylinders");
    assert_eq!((cones, tori), (0, 0), "no cones or tori");
    assert_eq!(planar, 6, "6 face planes");
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
    let [planar, cylinders, cones, _, tori, _] = selected_counts(&out);
    assert_eq!(tori, 1, "1 bend torus");
    assert_eq!(cylinders, 2, "2 straight tube cylinders");
    assert_eq!(cones, 0, "no cones");
    assert_eq!(planar, 2, "2 end-cap planes");
    assert_eq!(out.selected_surfaces.len(), 5, "torus + 2 cyl + 2 planes");
    assert_total_coverage(&out);
}

#[test]
fn competition_nosecone_cone_sphere_base() {
    // The nosecone's cone+sphere tangent chain is one region; the base disk
    // is a separate sharp-bounded plane. Competition + the ladder together
    // yield cone + sphere + base = 3 surfaces (the reference's 3 faces).
    let out = run("tests/onshape/nosecone_medium.stl");
    let [planar, _, cones, spheres, _, _] = selected_counts(&out);
    assert_eq!(spheres, 1, "1 nose sphere");
    assert_eq!(cones, 1, "1 body cone");
    assert_eq!(planar, 1, "1 base plane");
    assert_eq!(out.selected_surfaces.len(), 3, "cone + sphere + base");
    assert_total_coverage(&out);
}

#[test]
fn competition_stepped_block_fillet_network_stays_analytic() {
    // A dense fillet network — 37 cylinders + 18 tori + 16 planes in the
    // reference — tangent-connects the whole model into ONE 26122-face region,
    // so every surface has to come out of competition.
    //
    // This model used to collapse entirely to per-facet planes: one r = 0.5
    // fillet band, a single quad tall between two corner tori, left the axis so
    // weakly constrained that a bogus 12° cone beat the true cylinder on
    // residual, and the leftover blob that left behind tripped classify.rs's
    // whole-region discard. Then it decomposed but over-fitted: near-cylindrical
    // tori ate perimeter fillets and near-planar cylinders ate wall strips,
    // each winning on area as a super-model of the surface it displaced
    // (see `MIN_AXIAL_SWEEP_DEG`).
    //
    // The decomposition now matches the reference exactly, surface for surface,
    // and `onshape_stepped_block_toroidal_fine_stage41_compare` checks the
    // geometry that comes out of it.
    let out = run("tests/onshape/stepped_block_toroidal_fillets_fine.stl");
    let [planar, cylinders, cones, spheres, tori, bsplines] = selected_counts(&out);
    assert_eq!(tori, 18, "the reference's 18 corner-blend tori");
    assert_eq!(cylinders, 37, "the reference's 37 fillet cylinders");
    assert_eq!(planar, 16, "the reference's 16 planes");
    assert_eq!((cones, spheres, bsplines), (0, 0, 0), "no cones, spheres or freeform");
    assert_eq!(out.selected_surfaces.len(), 71, "the reference's 71 surfaces");
    assert_total_coverage(&out);
}

/// Every tangent junction must end up on the mesh-edge chain the two fitted
/// surfaces' junction geometry picks out (stage 2.4 boundary snap,
/// `locus.rs`). This is the property snap exists to establish, and it is
/// stronger than a surface count: a boundary one tessellation row off leaves
/// the counts intact and shows up only here (and, much later, as a stage-3
/// `--compare` edge distance).
///
/// The second number is how many adjacent pairs have a junction geometry
/// `locus.rs` does not implement — cone-plane ruling tangency, non-coaxial
/// pairs other than the cylinder-to-torus tube junction, point tangencies.
/// Every corpus model here reaches zero of both; a rise in the second is a new
/// junction type the module should learn, and stage 2.4 reports it by name.
#[test]
fn boundaries_sit_on_their_junctions() {
    for model in [
        "tests/onshape/stepped_block_toroidal_fillets_fine.stl",
        "tests/onshape/stepped_block_toroidal_fillets_medium.stl",
        "tests/onshape/rounded_cube_10_r2_fine.stl",
        "tests/onshape/pipe_elbow_10_medium.stl",
        "tests/onshape/pill_fine.stl",
        "tests/ccad/generated/filleted_cylinder.stl",
        "tests/ccad/generated/filleted_pipe.stl",
    ] {
        let (out, config) = run_with_config(model);
        let (unhandled, wrong_side) = stage2::locus_summary(&out, &config);
        assert_eq!(
            wrong_side, 0,
            "{model}: {wrong_side} face(s) sit on the wrong side of a tangent junction",
        );
        assert_eq!(
            unhandled, 0,
            "{model}: {unhandled} pair(s) with a tangent junction locus.rs cannot build",
        );
    }
}
