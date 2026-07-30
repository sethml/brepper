use brepper::{config, stage1, stage2, stage3, stage4};

fn manifest_dir() -> String {
    env!("CARGO_MANIFEST_DIR").to_string()
}

fn config_for_output(stl_path: &str, output_path: &str) -> config::Config {
    config::parse_config_from(["brepper", stl_path, "-o", output_path, "-q"]).unwrap()
}

fn config_for_compare(stl_path: &str, output_path: &str, step_path: &str) -> config::Config {
    let mut config = config::parse_config_from([
        "brepper", stl_path, "-o", output_path, "--compare", step_path, "-q",
    ])
    .unwrap();
    config.load_compare_step().unwrap();
    config
}

fn run_stage4(config: &config::Config) {
    let mesh = stage1::stage1(config).expect("stage1 should pass");
    let surfaces = stage2::stage2(config, mesh, None).expect("stage2 should pass");
    let brep = stage3::stage3(config, surfaces, None).expect("stage3 should pass");
    stage4::stage4(config, brep).expect("stage4 should pass");
}

// --- Basic STEP output tests ---

#[test]
fn cube_step_output() {
    let dir = manifest_dir();
    let stl = format!("{dir}/tests/ccad/generated/cube.stl");
    let out = format!("{dir}/tmp/test_cube_output.step");
    let config = config_for_output(&stl, &out);
    run_stage4(&config);
    assert!(std::path::Path::new(&out).exists());
}

#[test]
fn simple_cylinder_step_output() {
    let dir = manifest_dir();
    let stl = format!("{dir}/tests/ccad/generated/simple_cylinder.stl");
    let out = format!("{dir}/tmp/test_cylinder_output.step");
    let config = config_for_output(&stl, &out);
    run_stage4(&config);
    assert!(std::path::Path::new(&out).exists());
}

#[test]
fn simple_sphere_step_output() {
    let dir = manifest_dir();
    let stl = format!("{dir}/tests/ccad/generated/simple_sphere.stl");
    let out = format!("{dir}/tmp/test_sphere_output.step");
    let config = config_for_output(&stl, &out);
    run_stage4(&config);
    assert!(std::path::Path::new(&out).exists());
}
#[test]
fn rounded_cube_step_output() {
    // --angular-tolerance=20: the coarse tessellation exceeds the 17.5°
    // default dihedral limit (the same flag tools/sweep.sh applies to this model).
    let dir = manifest_dir();
    let stl = format!("{dir}/tests/onshape/rounded_cube_10_r2_coarse.stl");
    let out = format!("{dir}/tmp/test_rounded_cube_output.step");
    let config = config::parse_config_from(
        ["brepper", &stl, "-o", &out, "--angular-tolerance=20", "-q"],
    ).unwrap();
    run_stage4(&config);
    assert!(std::path::Path::new(&out).exists());
}


// --- Stage 4.1 compare tests ---

#[test]
fn ccad_cube_stage41_compare() {
    let dir = manifest_dir();
    let config = config_for_compare(
        &format!("{dir}/tests/ccad/generated/cube.stl"),
        &format!("{dir}/tmp/test_cube_cmp.step"),
        &format!("{dir}/tests/ccad/generated/cube.step"),
    );
    run_stage4(&config);
}

#[test]
fn ccad_wedge_stage41_compare() {
    let dir = manifest_dir();
    let config = config_for_compare(
        &format!("{dir}/tests/ccad/generated/wedge.stl"),
        &format!("{dir}/tmp/test_wedge_cmp.step"),
        &format!("{dir}/tests/ccad/generated/wedge.step"),
    );
    run_stage4(&config);
}

#[test]
fn ccad_simple_cylinder_stage41_compare() {
    let dir = manifest_dir();
    let config = config_for_compare(
        &format!("{dir}/tests/ccad/generated/simple_cylinder.stl"),
        &format!("{dir}/tmp/test_simple_cylinder_cmp.step"),
        &format!("{dir}/tests/ccad/generated/simple_cylinder.step"),
    );
    run_stage4(&config);
}

#[test]
fn ccad_block_with_hole_stage41_compare() {
    let dir = manifest_dir();
    let config = config_for_compare(
        &format!("{dir}/tests/ccad/generated/block_with_hole.stl"),
        &format!("{dir}/tmp/test_block_with_hole_cmp.step"),
        &format!("{dir}/tests/ccad/generated/block_with_hole.step"),
    );
    run_stage4(&config);
}

#[test]
fn ccad_pipe_stage41_compare() {
    let dir = manifest_dir();
    let config = config_for_compare(
        &format!("{dir}/tests/ccad/generated/pipe.stl"),
        &format!("{dir}/tmp/test_pipe_cmp.step"),
        &format!("{dir}/tests/ccad/generated/pipe.step"),
    );
    run_stage4(&config);
}

#[test]
fn ccad_simple_sphere_stage41_compare() {
    let dir = manifest_dir();
    let config = config_for_compare(
        &format!("{dir}/tests/ccad/generated/simple_sphere.stl"),
        &format!("{dir}/tmp/test_simple_sphere_cmp.step"),
        &format!("{dir}/tests/ccad/generated/simple_sphere.step"),
    );
    run_stage4(&config);
}

#[test]
fn ccad_hemisphere_stage41_compare() {
    let dir = manifest_dir();
    let config = config_for_compare(
        &format!("{dir}/tests/ccad/generated/hemisphere.stl"),
        &format!("{dir}/tmp/test_hemisphere_cmp.step"),
        &format!("{dir}/tests/ccad/generated/hemisphere.step"),
    );
    run_stage4(&config);
}

#[test]
fn ccad_ball_on_cylinder_stage41_compare() {
    let dir = manifest_dir();
    let config = config_for_compare(
        &format!("{dir}/tests/ccad/generated/ball_on_cylinder.stl"),
        &format!("{dir}/tmp/test_ball_on_cylinder_cmp.step"),
        &format!("{dir}/tests/ccad/generated/ball_on_cylinder.step"),
    );
    run_stage4(&config);
}

#[test]
fn ccad_spherical_pocket_stage41_compare() {
    let dir = manifest_dir();
    let config = config_for_compare(
        &format!("{dir}/tests/ccad/generated/spherical_pocket.stl"),
        &format!("{dir}/tmp/test_spherical_pocket_cmp.step"),
        &format!("{dir}/tests/ccad/generated/spherical_pocket.step"),
    );
    run_stage4(&config);
}

#[test]
fn onshape_chamfered_cube_stage41_compare() {
    let dir = manifest_dir();
    let config = config_for_compare(
        &format!("{dir}/tests/onshape/chamfered_cube_10_c1_medium.stl"),
        &format!("{dir}/tmp/test_chamfered_cube_cmp.step"),
        &format!("{dir}/tests/onshape/chamfered_cube_10_c1.step"),
    );
    run_stage4(&config);
}

#[test]
fn manual_cube_stage41_compare() {
    let dir = manifest_dir();
    let config = config_for_compare(
        &format!("{dir}/tests/manual/cube.stl"),
        &format!("{dir}/tmp/test_manual_cube_cmp.step"),
        &format!("{dir}/tests/manual/cube.step"),
    );
    run_stage4(&config);
}

#[test]
fn onshape_part_rounded_cube_stage41_compare() {
    let dir = manifest_dir();
    let config = config_for_compare(
        &format!("{dir}/tests/onshape/part_rounded_cube_10_r2_coarse.stl"),
        &format!("{dir}/tmp/test_part_rounded_cube_cmp.step"),
        &format!("{dir}/tests/onshape/part_rounded_cube_10_r2.step"),
    );
    run_stage4(&config);
}

#[test]
fn onshape_pill_coarse_stage41_compare() {
    let dir = manifest_dir();
    let config = config_for_compare(
        &format!("{dir}/tests/onshape/pill_coarse.stl"),
        &format!("{dir}/tmp/test_pill_coarse_cmp.step"),
        &format!("{dir}/tests/onshape/pill.step"),
    );
    run_stage4(&config);
}

#[test]
fn onshape_pill_fine_stage41_compare() {
    let dir = manifest_dir();
    let config = config_for_compare(
        &format!("{dir}/tests/onshape/pill_fine.stl"),
        &format!("{dir}/tmp/test_pill_fine_cmp.step"),
        &format!("{dir}/tests/onshape/pill.step"),
    );
    run_stage4(&config);
}

#[test]
fn onshape_dome_hemisphere_stage41_compare() {
    let dir = manifest_dir();
    let config = config_for_compare(
        &format!("{dir}/tests/onshape/dome_hemisphere_20_fine.stl"),
        &format!("{dir}/tmp/test_dome_hemisphere_cmp.step"),
        &format!("{dir}/tests/onshape/dome_hemisphere_20.step"),
    );
    run_stage4(&config);
}

// --- Cone pipeline tests ---

#[test]
fn simple_cone_step_output() {
    let dir = manifest_dir();
    let stl = format!("{dir}/tests/ccad/generated/simple_cone.stl");
    let out = format!("{dir}/tmp/test_simple_cone_output.step");
    let config = config_for_output(&stl, &out);
    run_stage4(&config);
    assert!(std::path::Path::new(&out).exists());
}

#[test]
fn ccad_simple_cone_stage41_compare() {
    let dir = manifest_dir();
    let config = config_for_compare(
        &format!("{dir}/tests/ccad/generated/simple_cone.stl"),
        &format!("{dir}/tmp/test_simple_cone_cmp.step"),
        &format!("{dir}/tests/ccad/generated/simple_cone.step"),
    );
    run_stage4(&config);
}

#[test]
fn ccad_cone_cylinder_stage41_compare() {
    let dir = manifest_dir();
    let config = config_for_compare(
        &format!("{dir}/tests/ccad/generated/cone_cylinder.stl"),
        &format!("{dir}/tmp/test_cone_cylinder_cmp.step"),
        &format!("{dir}/tests/ccad/generated/cone_cylinder.step"),
    );
    run_stage4(&config);
}

#[test]
fn ccad_nosecone_stage41_compare() {
    let dir = manifest_dir();
    let config = config_for_compare(
        &format!("{dir}/tests/ccad/generated/nosecone.stl"),
        &format!("{dir}/tmp/test_nosecone_cmp.step"),
        &format!("{dir}/tests/ccad/generated/nosecone.step"),
    );
    run_stage4(&config);
}

#[test]
fn ccad_block_with_conical_hole_stage41_compare() {
    let dir = manifest_dir();
    let config = config_for_compare(
        &format!("{dir}/tests/ccad/generated/block_with_conical_hole.stl"),
        &format!("{dir}/tmp/test_block_with_conical_hole_cmp.step"),
        &format!("{dir}/tests/ccad/generated/block_with_conical_hole.step"),
    );
    run_stage4(&config);
}

// --- Torus fillet models: full pipeline with compare ---

#[test]
fn ccad_filleted_cylinder_stage41_compare() {
    let dir = manifest_dir();
    let config = config_for_compare(
        &format!("{dir}/tests/ccad/generated/filleted_cylinder.stl"),
        &format!("{dir}/tmp/test_filleted_cylinder_cmp.step"),
        &format!("{dir}/tests/ccad/generated/filleted_cylinder.step"),
    );
    run_stage4(&config);
}

#[test]
fn ccad_filleted_pipe_stage41_compare() {
    let dir = manifest_dir();
    let config = config_for_compare(
        &format!("{dir}/tests/ccad/generated/filleted_pipe.stl"),
        &format!("{dir}/tmp/test_filleted_pipe_cmp.step"),
        &format!("{dir}/tests/ccad/generated/filleted_pipe.step"),
    );
    run_stage4(&config);
}

#[test]
fn ccad_filleted_hole_block_stage41_compare() {
    let dir = manifest_dir();
    let config = config_for_compare(
        &format!("{dir}/tests/ccad/generated/filleted_hole_block.stl"),
        &format!("{dir}/tmp/test_filleted_hole_block_cmp.step"),
        &format!("{dir}/tests/ccad/generated/filleted_hole_block.step"),
    );
    run_stage4(&config);
}

#[test]
#[ignore = "stepped_block_toroidal fine regressed 71 -> 138 surfaces when the self-calibrating acceptance tolerance was removed (2026-07-28). The tolerance was NOT doing its stated job -- it was masking a seed-conditioning defect: 323 seeds produce NO GROWTH, including FLAT seeds failing to grow a PLANE, which cannot be a tolerance problem since a wall's vertices sit at ~1e-15. The cause is condition_patch giving up at MIN_CONDITION_FACES without checking whether the seed validates, on a model whose flat regions are narrow strips between fillets so a 32-face BFS patch always straddles. THE CONSTANT CANNOT BE RETUNED: measured, the conical conditioning floor at 2 gives stepped_block 75 surfaces (cylinders and tori exactly right at 37/18) but hands filleted_pipe a spurious 18-face cone; at 6 filleted_pipe is correct and stepped_block is 136; at 4 both are wrong. It is a lever, not a bound -- same pathology as the tolerance it was hiding behind. Needs the seed-patch mechanism reworked so a patch cannot straddle, not a new value. See DEVELOPMENT_PLAN."]
fn onshape_stepped_block_toroidal_fine_stage41_compare() {
    // The heaviest model (51K triangles, 72 surfaces incl. 18 tori in a
    // fillet network). Exercises curvature-guided seeding, boundary
    // refinement, and sharp-corner vertex refinement end to end.
    let dir = manifest_dir();
    let config = config_for_compare(
        &format!("{dir}/tests/onshape/stepped_block_toroidal_fillets_fine.stl"),
        &format!("{dir}/tmp/test_stepped_block_fine_cmp.step"),
        &format!("{dir}/tests/onshape/stepped_block_toroidal_fillets.step"),
    );
    run_stage4(&config);
}

#[test]
#[ignore = "coarse sampling: the fine mesh passes; at 2615 faces this is 82 surfaces against the reference 71 (~4 fillet cylinders missing, ~14 arc facets still per-facet) and a 1.6e-1 mm reference-edge compare. Blocked on LEFTOVER_BLOB_MAX's whole-region discard, which makes the coarse path untunable — see TODO.md."]
fn onshape_stepped_block_toroidal_medium_stage41_compare() {
    // Same model at MEDIUM resolution (2.3K triangles): small-radius fillets
    // are tessellated with only ~4 rows across each 90° tube arc, so the
    // ~24° row-to-row steps require the coarse-tessellation coverage
    // allowance (angular_coverage_valid max_gap) and the raised default
    // angular tolerance to detect the corner-blend tori. Fast (~3 s), unlike
    // the fine variant, so it guards the medium-resolution path cheaply.
    let dir = manifest_dir();
    let config = config_for_compare(
        &format!("{dir}/tests/onshape/stepped_block_toroidal_fillets_medium.stl"),
        &format!("{dir}/tmp/test_stepped_block_medium_cmp.step"),
        &format!("{dir}/tests/onshape/stepped_block_toroidal_fillets.step"),
    );
    run_stage4(&config);
}

#[test]
fn occt_full_torus_stage41_compare() {
    // Complete donut: exercises both U and V seams on the doubly-periodic
    // torus face (natural-bounds construction, no boundary edges).
    let dir = manifest_dir();
    let config = config_for_compare(
        &format!("{dir}/tests/occt/full_torus.stl"),
        &format!("{dir}/tmp/test_full_torus_cmp.step"),
        &format!("{dir}/tests/occt/full_torus.step"),
    );
    run_stage4(&config);
}

#[test]
fn occt_hollow_cube_stage41_compare() {
    // Multi-shell solid: one solid with an internal void (volume 7000 mm^3).
    let dir = manifest_dir();
    let config = config_for_compare(
        &format!("{dir}/tests/occt/hollow_cube.stl"),
        &format!("{dir}/tmp/test_hollow_cube_cmp.step"),
        &format!("{dir}/tests/occt/hollow_cube.step"),
    );
    run_stage4(&config);
}

// --- Onshape cone models: full pipeline with compare ---

#[test]
fn onshape_nosecone_medium_stage41_compare() {
    let dir = manifest_dir();
    let config = config_for_compare(
        &format!("{dir}/tests/onshape/nosecone_medium.stl"),
        &format!("{dir}/tmp/test_nosecone_medium_cmp.step"),
        &format!("{dir}/tests/onshape/nosecone.step"),
    );
    run_stage4(&config);
}

#[test]
fn onshape_pipe_elbow_medium_stage41_compare() {
    // The elbow torus is tessellated into narrow strips that stage 2.2 fits
    // as ~83 cylinders; strip-merge torus seeding (stage 2.5) must recover
    // the single quarter-torus for the full pipeline to pass.
    let dir = manifest_dir();
    let config = config_for_compare(
        &format!("{dir}/tests/onshape/pipe_elbow_10_medium.stl"),
        &format!("{dir}/tmp/test_pipe_elbow_medium_cmp.step"),
        &format!("{dir}/tests/onshape/pipe_elbow_10.step"),
    );
    run_stage4(&config);
}

#[test]
fn onshape_pipe_elbow_fine_stage41_compare() {
    // Fine tessellation of the same quarter-torus elbow. Under OCCT 7.9.3 the
    // in-memory solid was correct but the WRITTEN STEP read back 14% low
    // (partial-torus pcurve bug); OCCT 8.0.0_p1 fixes the round-trip, so the
    // --compare STEP round-trip self-check (stage4.rs) now passes here.
    let dir = manifest_dir();
    let config = config_for_compare(
        &format!("{dir}/tests/onshape/pipe_elbow_10_fine.stl"),
        &format!("{dir}/tmp/test_pipe_elbow_fine_cmp.step"),
        &format!("{dir}/tests/onshape/pipe_elbow_10.step"),
    );
    run_stage4(&config);
}

#[test]
fn onshape_cone_truncated_stage41_compare() {
    let dir = manifest_dir();
    let config = config_for_compare(
        &format!("{dir}/tests/onshape/cone_truncated_medium.stl"),
        &format!("{dir}/tmp/test_cone_truncated_cmp.step"),
        &format!("{dir}/tests/onshape/cone_truncated.step"),
    );
    run_stage4(&config);
}

#[test]
fn onshape_cone_stage41_compare() {
    let dir = manifest_dir();
    let config = config_for_compare(
        &format!("{dir}/tests/onshape/cone_15x20_medium.stl"),
        &format!("{dir}/tmp/test_onshape_cone_cmp.step"),
        &format!("{dir}/tests/onshape/cone_15x20.step"),
    );
    run_stage4(&config);
}

// --- Missing output path test ---

#[test]
fn missing_output_path_error() {
    let dir = manifest_dir();
    let stl = format!("{dir}/tests/ccad/generated/cube.stl");
    let config = config::parse_config_from(["brepper", &stl, "-q"]).unwrap();
    let mesh = stage1::stage1(&config).expect("stage1 should pass");
    let surfaces = stage2::stage2(&config, mesh, None).expect("stage2 should pass");
    let brep = stage3::stage3(&config, surfaces, None).expect("stage3 should pass");
    let result = stage4::stage4(&config, brep);
    assert!(result.is_err());
    let err = result.unwrap_err();
    assert!(
        err.to_string().contains("no output STEP file path"),
        "expected MissingOutputPath error, got: {err}"
    );
}

// ---------------------------------------------------------------------------
// NURBS corpus: end-to-end --compare
// ---------------------------------------------------------------------------
//
// Added 2026-07-28. Before this the NURBS models had stage-2 *selection* tests
// and nothing end to end, which is how three of the eight came to be silently
// broken: `bspline_blister` and `bspline_saddle` had collapsed to all-per-facet
// polyhedra and `hill_with_hole` was being shredded into 47 surfaces against a
// reference 7, all unnoticed because only stage 1 was checked. Every model has
// a reference STEP sitting next to it, so there was never a reason not to.

macro_rules! nurbs_compare {
    ($name:ident, $model:literal) => {
        #[test]
        fn $name() {
            let dir = manifest_dir();
            let config = config_for_compare(
                &format!("{dir}/tests/nurbs/{}.stl", $model),
                &format!("{dir}/tmp/test_nurbs_{}.step", $model),
                &format!("{dir}/tests/nurbs/{}.step", $model),
            );
            run_stage4(&config);
        }
    };
    ($name:ident, $model:literal, ignore = $why:literal) => {
        #[test]
        #[ignore = $why]
        fn $name() {
            let dir = manifest_dir();
            let config = config_for_compare(
                &format!("{dir}/tests/nurbs/{}.stl", $model),
                &format!("{dir}/tmp/test_nurbs_{}.step", $model),
                &format!("{dir}/tests/nurbs/{}.step", $model),
            );
            run_stage4(&config);
        }
    };
}

nurbs_compare!(nurbs_bspline_hill_stage41_compare, "bspline_hill");
nurbs_compare!(nurbs_bspline_terrain_stage41_compare, "bspline_terrain");
nurbs_compare!(nurbs_bspline_wave_stage41_compare, "bspline_wave");
nurbs_compare!(nurbs_elliptical_arch_stage41_compare, "elliptical_arch");
nurbs_compare!(nurbs_elliptical_cylinder_stage41_compare, "elliptical_cylinder");

nurbs_compare!(
    nurbs_hill_with_hole_stage41_compare,
    "hill_with_hole",
    ignore = "8 surfaces against a reference 7 (5 planes + 1 cylinder + 1 B-spline); stage 3.1 edge \
compare 7.67 mm. Was 47 surfaces before 2026-07-28 — removing the self-calibrating acceptance \
tolerance stopped competition carving the freeform hill into spurious analytic patches, and it \
briefly reached exactly 7 surfaces (failing then at stage 4.1 on volume, 6.96e-2 relative) before \
the 2.3 rework steps took it to 8. One surplus surface to find."
);
nurbs_compare!(
    nurbs_bspline_blister_stage41_compare,
    "bspline_blister",
    ignore = "7 surfaces against a reference 6 (1 B-spline + 5 planes); stage 3.1 edge compare \
9.04 mm. Was 968 surfaces (959 of them per-facet) before 2026-07-28. Also builds 4 invalid faces, \
all single-mesh-face per-facet planes — the micro-sliver family, tracked in stage3_integration."
);
nurbs_compare!(
    nurbs_bspline_saddle_stage41_compare,
    "bspline_saddle",
    ignore = "still collapsed to an all-per-facet polyhedron (1605 planar surfaces): three of its \
freeform components fail the §2.7 pole-bound check ('control points 246.7 / 19378.7 / 573.3 mm \
outside the data bounding box'), so nothing explains them. Fix the freeform fit first."
);

// ---------------------------------------------------------------------------
// Coverage lock (Phase 0.1, 2026-07-29)
// ---------------------------------------------------------------------------
//
// Every .stl under tests/ (excluding tests/bad/) was audited for coverage past
// stage 1. Nineteen models had a reference STEP sitting next to them and no
// end-to-end test — the same hole that let three NURBS models rot silently.
// All nineteen pass today; measured relative volume errors are recorded per
// case below, so a future change that breaks one is visible rather than silent.
//
// Not covered here, deliberately:
//   * tests/downloaded/ (3DBenchy, OrcaCube, xyz-cube-test) — real-world
//     meshes with no verified provenance; Phase 0.8 baselines them as harness
//     rows, not as suite gates.
//   * tests/onshape/chamfered_cube_10_c1_medium.ascii.stl — the same geometry
//     as the binary variant, which is covered through stage 4. Its purpose is
//     the ASCII reader; pinned by `chamfered_cube_ascii_matches_binary` in
//     stage1_integration.rs rather than by a duplicate pipeline run.

macro_rules! stage41_compare {
    ($name:ident, $stl:literal, $step:literal) => {
        #[test]
        fn $name() {
            let dir = manifest_dir();
            let config = config_for_compare(
                &format!("{dir}/tests/{}.stl", $stl),
                &format!("{dir}/tmp/test_cover_{}.step", stringify!($name)),
                &format!("{dir}/tests/{}.step", $step),
            );
            run_stage4(&config);
        }
    };
    ($name:ident, $stl:literal, $step:literal, ignore = $why:literal) => {
        #[test]
        #[ignore = $why]
        fn $name() {
            let dir = manifest_dir();
            let config = config_for_compare(
                &format!("{dir}/tests/{}.stl", $stl),
                &format!("{dir}/tmp/test_cover_{}.step", stringify!($name)),
                &format!("{dir}/tests/{}.step", $step),
            );
            run_stage4(&config);
        }
    };
}

// --- tests/ccad/generated/ (previously stage 1 + stage 2 only) ---
// staircase 12 planes, volume rel diff 1.8e-16; t_shape 10 planes, exact.
stage41_compare!(
    ccad_staircase_stage41_compare,
    "ccad/generated/staircase",
    "ccad/generated/staircase"
);
stage41_compare!(
    ccad_t_shape_stage41_compare,
    "ccad/generated/t_shape",
    "ccad/generated/t_shape"
);
stage41_compare!(
    ccad_stepped_cylinder_stage41_compare,
    "ccad/generated/stepped_cylinder",
    "ccad/generated/stepped_cylinder"
);
stage41_compare!(
    ccad_two_holes_stage41_compare,
    "ccad/generated/two_holes",
    "ccad/generated/two_holes"
);
// block_corner_hole had NO test at any stage before this.
stage41_compare!(
    ccad_block_corner_hole_stage41_compare,
    "ccad/generated/block_corner_hole",
    "ccad/generated/block_corner_hole"
);

// --- plate_with_hole across two exporters and four densities ---
// The only model in the corpus with both Fusion and Onshape tessellations of
// one reference; all four decompose to 6 planes + 1 cylinder (rel diff
// 2.8e-10 low / 2.9e-10 medium / 6.3e-10 high / 1.6e-9 Onshape coarse).
stage41_compare!(
    fusion_plate_low_stage41_compare,
    "fusion/plate_with_hole_100x50_low",
    "fusion/plate_with_hole_100x50"
);
stage41_compare!(
    fusion_plate_medium_stage41_compare,
    "fusion/plate_with_hole_100x50_medium",
    "fusion/plate_with_hole_100x50"
);
stage41_compare!(
    fusion_plate_high_stage41_compare,
    "fusion/plate_with_hole_100x50_high",
    "fusion/plate_with_hole_100x50"
);
stage41_compare!(
    onshape_plate_coarse_stage41_compare,
    "onshape/plate_with_hole_100x50_coarse",
    "onshape/plate_with_hole_100x50"
);

// --- tests/onshape/ analytic primitives (previously stage 1 + stage 2 only) ---
stage41_compare!(
    onshape_cylinder_10x30_stage41_compare,
    "onshape/cylinder_10x30_medium",
    "onshape/cylinder_10x30"
);
stage41_compare!(
    onshape_cylinder_10_r2_stage41_compare,
    "onshape/cylinder_10_r2_fine",
    "onshape/cylinder_10_r2"
);
stage41_compare!(
    onshape_sphere_25_stage41_compare,
    "onshape/sphere_25_fine",
    "onshape/sphere_25"
);
stage41_compare!(
    onshape_l_bracket_stage41_compare,
    "onshape/l_bracket_simple_medium",
    "onshape/l_bracket_simple"
);
stage41_compare!(
    onshape_stepped_block_coarse_stage41_compare,
    "onshape/stepped_block_coarse",
    "onshape/stepped_block"
);

// --- rounded_cube: the 6-plane/12-cylinder/8-sphere fillet network ---
// Both densities reach exactly 26 surfaces (rel diff 9.4e-8 fine, 1.2e-7
// medium). The fine variant costs ~27 s, the heaviest non-ignored test here;
// it is worth it — this is the model the sphere/cylinder-fillet machinery was
// built against, and it had no end-to-end gate.
stage41_compare!(
    onshape_rounded_cube_fine_stage41_compare,
    "onshape/rounded_cube_10_r2_fine",
    "onshape/rounded_cube_10_r2"
);
stage41_compare!(
    onshape_rounded_cube_medium_stage41_compare,
    "onshape/rounded_cube_10_r2_medium",
    "onshape/rounded_cube_10_r2"
);
stage41_compare!(
    onshape_part_rounded_cube_fine_stage41_compare,
    "onshape/part_rounded_cube_10_r2_fine",
    "onshape/part_rounded_cube_10_r2"
);

// nosecone_truncated: an older measurement had this missing 1.4e-2 of the
// reference volume. Re-measured 2026-07-29 at 1.08e-7 (2 planes + 1 cone +
// 1 sphere) — the old number is stale, whatever fixed it went unrecorded.
stage41_compare!(
    onshape_nosecone_truncated_stage41_compare,
    "onshape/nosecone_truncated_medium",
    "onshape/nosecone_truncated"
);

// cone_tangent_planes: passes at rel diff 5.2e-3 against a 1% gate — three
// decimals of margin where every other model here has six or more. Its
// freeform rung fails as ill-conditioned, so the analytic decomposition wins
// by accident (DEVELOPMENT_PLAN, Phase 5c); this test is the canary for that.
stage41_compare!(
    onshape_cone_tangent_planes_stage41_compare,
    "onshape/cone_tangent_planes_medium",
    "onshape/cone_tangent_planes"
);

// ---------------------------------------------------------------------------
// OCCT re-meshes (Phase 0.2, 2026-07-29)
// ---------------------------------------------------------------------------
//
// tests/remesh/ re-tessellates references that already exist elsewhere in
// tests/, so the ground truth is identical and only the mesh PROVENANCE
// changes: constrained Delaunay instead of Parasolid strips, and (for the
// *_offset* variants) coordinates ~1 m from the origin, which raises the STL's
// f32 quantum from ~3e-6 mm to ~1e-4 mm. See tests/remesh/README.md.
//
// 25 of 32 pass. The 7 that fail are findings, not corpus defects — the strip
// assumption becoming measurable — and are ignored with their numbers below.

stage41_compare!(
    remesh_cylinder_coarse,
    "remesh/cylinder_10x30_occt_coarse",
    "onshape/cylinder_10x30"
);
stage41_compare!(
    remesh_cylinder_medium,
    "remesh/cylinder_10x30_occt_medium",
    "onshape/cylinder_10x30"
);
stage41_compare!(
    remesh_cylinder_fine,
    "remesh/cylinder_10x30_occt_fine",
    "onshape/cylinder_10x30"
);
stage41_compare!(
    remesh_cylinder_offset,
    "remesh/cylinder_10x30_offset_occt_medium",
    "remesh/cylinder_10x30_offset"
);

stage41_compare!(
    remesh_sphere_medium,
    "remesh/sphere_25_occt_medium",
    "onshape/sphere_25"
);
stage41_compare!(remesh_sphere_fine, "remesh/sphere_25_occt_fine", "onshape/sphere_25");
stage41_compare!(
    remesh_sphere_offset,
    "remesh/sphere_25_offset_occt_medium",
    "remesh/sphere_25_offset"
);
stage41_compare!(
    remesh_sphere_coarse,
    "remesh/sphere_25_occt_coarse",
    "onshape/sphere_25",
    ignore = "304 Delaunay triangles over a whole r=12.5 sphere (deflection 0.5 mm, 4.65e-2 \
tessellation volume loss): decomposes to 6 planes + 2 spheres instead of 1 sphere, failing stage \
3.1 with an 8.40 mm edge distance. The medium (1256 tris) and fine (7302 tris) variants are exact. \
A sphere this coarse has no facet pair whose normals separate enough to be spherical rather than \
polyhedral, which is an honest identifiability limit, not obviously a bug — measure before fixing."
);

stage41_compare!(remesh_pill_coarse, "remesh/pill_occt_coarse", "onshape/pill");
stage41_compare!(remesh_pill_medium, "remesh/pill_occt_medium", "onshape/pill");
stage41_compare!(remesh_pill_fine, "remesh/pill_occt_fine", "onshape/pill");
stage41_compare!(
    remesh_pill_offset,
    "remesh/pill_offset_occt_medium",
    "remesh/pill_offset"
);

stage41_compare!(
    remesh_plate_coarse,
    "remesh/plate_with_hole_100x50_occt_coarse",
    "onshape/plate_with_hole_100x50"
);
stage41_compare!(
    remesh_plate_medium,
    "remesh/plate_with_hole_100x50_occt_medium",
    "onshape/plate_with_hole_100x50"
);
stage41_compare!(
    remesh_plate_fine,
    "remesh/plate_with_hole_100x50_occt_fine",
    "onshape/plate_with_hole_100x50"
);
stage41_compare!(
    remesh_plate_offset,
    "remesh/plate_with_hole_100x50_offset_occt_medium",
    "remesh/plate_with_hole_100x50_offset"
);

stage41_compare!(
    remesh_rounded_cube_medium,
    "remesh/rounded_cube_10_r2_occt_medium",
    "onshape/rounded_cube_10_r2"
);
stage41_compare!(
    remesh_rounded_cube_fine,
    "remesh/rounded_cube_10_r2_occt_fine",
    "onshape/rounded_cube_10_r2"
);
stage41_compare!(
    remesh_rounded_cube_offset,
    "remesh/rounded_cube_10_r2_offset_occt_medium",
    "remesh/rounded_cube_10_r2_offset"
);
stage41_compare!(
    remesh_rounded_cube_coarse,
    "remesh/rounded_cube_10_r2_occt_coarse",
    "onshape/rounded_cube_10_r2",
    ignore = "34 surfaces against the reference 26 (8 planes / 8 cylinders / 14 spheres / 4 tori — \
note the 4 spurious tori and 6 surplus spheres where the reference has 12 cylinders and 8 \
spheres), failing stage 3.1 with a 0.504 mm edge distance at a sphere|torus interface. 580 \
triangles for the whole fillet network: each r=2 fillet quadrant gets ~2 rows. Medium (1668) and \
fine (8676) are exact at 26. The Parasolid coarse mesh of the same body passes, so this is \
provenance-specific."
);

stage41_compare!(
    remesh_cone_tangent_medium,
    "remesh/cone_tangent_planes_occt_medium",
    "onshape/cone_tangent_planes"
);
stage41_compare!(
    remesh_cone_tangent_fine,
    "remesh/cone_tangent_planes_occt_fine",
    "onshape/cone_tangent_planes"
);
stage41_compare!(
    remesh_cone_tangent_offset,
    "remesh/cone_tangent_planes_offset_occt_medium",
    "remesh/cone_tangent_planes_offset"
);
stage41_compare!(
    remesh_cone_tangent_coarse,
    "remesh/cone_tangent_planes_occt_coarse",
    "onshape/cone_tangent_planes",
    ignore = "10 surfaces against a reference 7, failing stage 3.1 with a 3.71 mm edge distance at \
a cone|plane interface. 202 triangles across three tangent cone-plane junctions. Stage 2.4 also \
reports 'cone-plane ruling tangency is not implemented' twice (NO LOCUS), so the snap has nothing \
to work with here. Medium (544 tris) and fine (2392) pass — note they resolve as 4 planes + 1 cone \
+ 2 B-SPLINES, a different decomposition from the Parasolid mesh's 4 planes + 2 cones + 1 sphere, \
both inside the 1% volume gate. Phase 5c."
);

stage41_compare!(
    remesh_pipe_elbow_coarse,
    "remesh/pipe_elbow_10_occt_coarse",
    "onshape/pipe_elbow_10"
);
stage41_compare!(
    remesh_pipe_elbow_medium,
    "remesh/pipe_elbow_10_occt_medium",
    "onshape/pipe_elbow_10"
);
stage41_compare!(
    remesh_pipe_elbow_fine,
    "remesh/pipe_elbow_10_occt_fine",
    "onshape/pipe_elbow_10"
);
stage41_compare!(
    remesh_pipe_elbow_offset,
    "remesh/pipe_elbow_10_offset_occt_medium",
    "remesh/pipe_elbow_10_offset"
);

// stepped_block_toroidal fails at EVERY OCCT density while the Parasolid fine
// mesh of the same body passes at exactly 71 surfaces. The 18 tori are
// recovered every time; the planes and cylinders around them fragment. This is
// the clearest measurement in the re-mesh set of what the strip assumption is
// worth, and it is a Phase 1 acceptance target.
stage41_compare!(
    remesh_stepped_block_coarse,
    "remesh/stepped_block_toroidal_occt_coarse",
    "onshape/stepped_block_toroidal_fillets",
    ignore = "183 surfaces against the reference 71 (149 planes / 16 cylinders / 18 tori), stage \
3.1 vertex distance 0.783 mm. 1948 triangles. All 18 tori are found; the 149 planes are the \
fragmentation. Stage 2.4 reports many 'NO LOCUS: surfaces do not share an axis' snaps."
);
stage41_compare!(
    remesh_stepped_block_medium,
    "remesh/stepped_block_toroidal_occt_medium",
    "onshape/stepped_block_toroidal_fillets",
    ignore = "88 surfaces against the reference 71 (34 planes / 33 cylinders / 18 tori / 3 \
B-splines), stage 3.1 vertex distance 0.131 mm at vertex 92. 6464 triangles — comparable to the \
Parasolid MEDIUM mesh (2615 faces) which reaches 82, and far from the Parasolid FINE mesh which \
reaches exactly 71. The tori are all correct; cylinders are 33 against 37 and the surplus is in \
planes plus 3 spurious B-splines."
);
stage41_compare!(
    remesh_stepped_block_fine,
    "remesh/stepped_block_toroidal_occt_fine",
    "onshape/stepped_block_toroidal_fillets",
    ignore = "97 surfaces against the reference 71 (41 planes / 36 cylinders / 18 tori / 2 \
B-splines) at 35596 triangles — MORE surfaces than the same body's medium re-mesh, so this is not \
a resolution problem. Stage 3.1 vertex distance 0.618 mm. The Parasolid fine mesh of this exact \
reference passes at 71. The single most valuable failing row in tests/remesh/: whatever the strip \
tessellation was providing, more Delaunay triangles do not supply it."
);
stage41_compare!(
    remesh_stepped_block_offset,
    "remesh/stepped_block_toroidal_offset_occt_medium",
    "remesh/stepped_block_toroidal_offset",
    ignore = "96 surfaces against the reference 71, stage 3.1 vertex distance 0.131 mm — the same \
failure as the origin-centred medium re-mesh (88 surfaces, 0.131 mm), plus 8 more surfaces. The \
extra 8 are the only part attributable to the 1 m translation; fix the medium case first, then \
re-measure this one to isolate the quantum's contribution."
);

// ---------------------------------------------------------------------------
// Fusion 360 re-exports (Phase 0.5a, 2026-07-30)
// ---------------------------------------------------------------------------
//
// A THIRD tessellator over references that already live in tests/ — Parasolid
// strips (Onshape) and constrained Delaunay (OCCT, tests/remesh/) being the
// other two. Generated by tools/fusion/BrepperReexport.py at the SAME requested
// tolerances tests/remesh asks of OCCT (surface deviation 0.5/0.1/0.02 mm at
// 30/15/6°), so Fusion-vs-OCCT is a controlled comparison rather than a
// comparison of two vendors' preset names.
//
// All 31 pass the Phase 0.5c ground-truth gate: mesh vertices lie on the
// original reference surfaces to within 1.8e-7..7.3e-6 mm.
//
// FOUR of these fail at STAGE 1, not in the decomposition: a mesh exported at
// 0.5 mm chordal deviation has facet-centroid sag of 0.51..0.83 mm, over the
// 0.4 mm default --surface-tolerance. That is the tolerance regime rejecting an
// ordinary coarse export, and it is a Phase 3 datum, not a bug in stage 2/3.

stage41_compare!(
    fusion_chamfered_cube_10_c1_coarse,
    "fusion_reexport/chamfered_cube_10_c1_fusion_coarse",
    "onshape/chamfered_cube_10_c1"
);
stage41_compare!(
    fusion_chamfered_cube_10_c1_fine,
    "fusion_reexport/chamfered_cube_10_c1_fusion_fine",
    "onshape/chamfered_cube_10_c1"
);
stage41_compare!(
    fusion_chamfered_cube_10_c1_medium,
    "fusion_reexport/chamfered_cube_10_c1_fusion_medium",
    "onshape/chamfered_cube_10_c1"
);
stage41_compare!(
    fusion_cone_tangent_planes_coarse,
    "fusion_reexport/cone_tangent_planes_fusion_coarse",
    "onshape/cone_tangent_planes",
    ignore = "14 surfaces against a reference 7 (12/0/1/1/0/0). stage 3.1 compare: vertex 1 is 5.018723e-5 mm \
from nearest STEP vertex (tolerance: 5.000000e-5 mm)"
);
stage41_compare!(
    fusion_cone_tangent_planes_fine,
    "fusion_reexport/cone_tangent_planes_fusion_fine",
    "onshape/cone_tangent_planes"
);
stage41_compare!(
    fusion_cone_tangent_planes_medium,
    "fusion_reexport/cone_tangent_planes_fusion_medium",
    "onshape/cone_tangent_planes",
    ignore = "5 surfaces against a reference 7 (2/0/1/1/0/1). stage 3.1 adjacency error: face 2: winding loops \
have incompatible seam anchors (0.0000 vs 2.1581 rad); edge splitting not implemented"
);
stage41_compare!(
    fusion_cylinder_10x30_coarse,
    "fusion_reexport/cylinder_10x30_fusion_coarse",
    "onshape/cylinder_10x30",
    ignore = "14 surfaces against a reference 3 (14/0/0/0/0/0). stage 3.1 compare: vertex 0 is 2.588190e0 mm \
from nearest STEP vertex (tolerance: 5.000000e-5 mm)"
);
stage41_compare!(
    fusion_l_bracket_simple_coarse,
    "fusion_reexport/l_bracket_simple_fusion_coarse",
    "onshape/l_bracket_simple",
    ignore = "32 surfaces against a reference 10 (32/0/0/0/0/0). stage 3.1 compare: vertex 5 is 6.000000e0 mm \
from nearest STEP vertex (tolerance: 5.000000e-5 mm)"
);
stage41_compare!(
    fusion_l_bracket_simple_fine,
    "fusion_reexport/l_bracket_simple_fusion_fine",
    "onshape/l_bracket_simple"
);
stage41_compare!(
    fusion_l_bracket_simple_medium,
    "fusion_reexport/l_bracket_simple_fusion_medium",
    "onshape/l_bracket_simple"
);
stage41_compare!(
    fusion_nosecone_coarse,
    "fusion_reexport/nosecone_fusion_coarse",
    "onshape/nosecone",
    ignore = "rejected before stage 2. stage 1.2 (compare): comparison failed: face 2 centroid is 8.301656e-1 \
mm from nearest STEP surface (tolerance: 4.000000e-1 mm)"
);
stage41_compare!(
    fusion_nosecone_fine,
    "fusion_reexport/nosecone_fusion_fine",
    "onshape/nosecone",
    ignore = "4 surfaces against a reference 3 (2/0/1/1/0/0). stage 3.1 adjacency error: face 1: 3 winding \
boundary loops (at most 2 supported)"
);
stage41_compare!(
    fusion_nosecone_medium,
    "fusion_reexport/nosecone_fusion_medium",
    "onshape/nosecone"
);
stage41_compare!(
    fusion_part_rounded_cube_10_r2_coarse,
    "fusion_reexport/part_rounded_cube_10_r2_fusion_coarse",
    "onshape/part_rounded_cube_10_r2",
    ignore = "10 surfaces against a reference 10 (6/4/0/0/0/0). stage 3.1 compare: vertex 0 is 1.035277e0 mm \
from nearest STEP vertex (tolerance: 5.000000e-5 mm)"
);
stage41_compare!(
    fusion_part_rounded_cube_10_r2_fine,
    "fusion_reexport/part_rounded_cube_10_r2_fusion_fine",
    "onshape/part_rounded_cube_10_r2"
);
stage41_compare!(
    fusion_part_rounded_cube_10_r2_medium,
    "fusion_reexport/part_rounded_cube_10_r2_fusion_medium",
    "onshape/part_rounded_cube_10_r2"
);
stage41_compare!(
    fusion_pill_coarse,
    "fusion_reexport/pill_fusion_coarse",
    "onshape/pill",
    ignore = "14 surfaces against a reference 3 (12/0/0/2/0/0). stage 3.1 compare: vertex 0 is 1.035275e0 mm \
from nearest STEP vertex (tolerance: 5.000000e-5 mm)"
);
stage41_compare!(
    fusion_pipe_elbow_10_coarse,
    "fusion_reexport/pipe_elbow_10_fusion_coarse",
    "onshape/pipe_elbow_10",
    ignore = "rejected before stage 2. stage 1.2 (compare): comparison failed: face 21 centroid is 6.511208e-1 \
mm from nearest STEP surface (tolerance: 4.000000e-1 mm)"
);
stage41_compare!(
    fusion_plate_with_hole_100x50_coarse,
    "fusion_reexport/plate_with_hole_100x50_fusion_coarse",
    "fusion/plate_with_hole_100x50",
    ignore = "18 surfaces against a reference 7 (18/0/0/0/0/0). stage 3.1 compare: vertex 8 is 2.588192e0 mm \
from nearest STEP vertex (tolerance: 5.000000e-5 mm)"
);
stage41_compare!(
    fusion_rounded_cube_10_r2_coarse,
    "fusion_reexport/rounded_cube_10_r2_fusion_coarse",
    "onshape/rounded_cube_10_r2",
    ignore = "70 surfaces against a reference 26 (61/3/0/4/0/2). stage 3.1 compare: edge 2 is 8.881986e-1 mm \
from nearest STEP edge (tolerance: 5.000000e-5 mm)"
);
stage41_compare!(
    fusion_rounded_cube_10_r2_fine,
    "fusion_reexport/rounded_cube_10_r2_fusion_fine",
    "onshape/rounded_cube_10_r2"
);
stage41_compare!(
    fusion_rounded_cube_10_r2_medium,
    "fusion_reexport/rounded_cube_10_r2_fusion_medium",
    "onshape/rounded_cube_10_r2",
    ignore = "9 surfaces against a reference 26 (0/0/0/8/0/1). stage 3.6 compare: solid 0 volume differs from \
STEP by 1.94e-2 relative (tolerance 1e-2)"
);
stage41_compare!(
    fusion_sphere_25_coarse,
    "fusion_reexport/sphere_25_fusion_coarse",
    "onshape/sphere_25",
    ignore = "rejected before stage 2. stage 1.2 (compare): comparison failed: face 37 centroid is 8.093293e-1 \
mm from nearest STEP surface (tolerance: 4.000000e-1 mm)"
);
stage41_compare!(
    fusion_sphere_25_fine,
    "fusion_reexport/sphere_25_fusion_fine",
    "onshape/sphere_25"
);
stage41_compare!(
    fusion_sphere_25_medium,
    "fusion_reexport/sphere_25_fusion_medium",
    "onshape/sphere_25"
);
stage41_compare!(
    fusion_stepped_block_coarse,
    "fusion_reexport/stepped_block_fusion_coarse",
    "onshape/stepped_block"
);
stage41_compare!(
    fusion_stepped_block_fine,
    "fusion_reexport/stepped_block_fusion_fine",
    "onshape/stepped_block"
);
stage41_compare!(
    fusion_stepped_block_medium,
    "fusion_reexport/stepped_block_fusion_medium",
    "onshape/stepped_block"
);
stage41_compare!(
    fusion_stepped_block_toroidal_fillets_coarse,
    "fusion_reexport/stepped_block_toroidal_fillets_fusion_coarse",
    "onshape/stepped_block_toroidal_fillets",
    ignore = "rejected before stage 2. stage 1.2 (compare): comparison failed: face 141 centroid is \
5.111130e-1 mm from nearest STEP surface (tolerance: 4.000000e-1 mm)"
);
stage41_compare!(
    fusion_stepped_block_toroidal_fillets_fine,
    "fusion_reexport/stepped_block_toroidal_fillets_fusion_fine",
    "onshape/stepped_block_toroidal_fillets",
    ignore = "82 surfaces against a reference 71 (27/36/0/0/18/1). stage 3.1 compare: vertex 92 is 1.046676e-1 \
mm from nearest STEP vertex (tolerance: 5.000000e-5 mm)"
);
stage41_compare!(
    fusion_stepped_block_toroidal_fillets_medium,
    "fusion_reexport/stepped_block_toroidal_fillets_fusion_medium",
    "onshape/stepped_block_toroidal_fillets",
    ignore = "129 surfaces against a reference 71 (82/29/0/0/18/0). stage 3.1 compare: edge 282 is 2.610530e-1 \
mm from nearest STEP edge (tolerance: 5.000000e-5 mm)"
);
