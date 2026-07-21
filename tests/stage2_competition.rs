//! Integration tests for stage 2.3(b) scoped fit competition: tangent-chain
//! regions the whole-region ladder cannot explain must decompose into their
//! true analytic surfaces (vertex-tolerance validated, not just counted).
//!
//! These run the full stage-2 pipeline (default stage) and assert the typed
//! hypothesis counts on the two headline acceptance models. Torus-requiring
//! chains (pipe elbow, toroidal fillets) are covered once stage 2.3(c)
//! promotion lands.

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
