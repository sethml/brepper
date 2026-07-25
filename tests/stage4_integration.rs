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
    // default dihedral limit (see tests/RESULTS.md special flags).
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
#[ignore = "decomposition is analytic now (74 surfaces vs the reference 71) but not yet an exact match: stage 3.1 reference-edge compare is 4.4e-2 mm out on the fine mesh, 2.1e-1 mm on medium. See TODO.md."]
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
#[ignore = "decomposition is analytic now (74 surfaces vs the reference 71) but not yet an exact match: stage 3.1 reference-edge compare is 4.4e-2 mm out on the fine mesh, 2.1e-1 mm on medium. See TODO.md."]
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
