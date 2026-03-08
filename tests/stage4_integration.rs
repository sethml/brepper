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
    let surfaces = stage2::stage2(config, mesh).expect("stage2 should pass");
    let brep = stage3::stage3(config, surfaces).expect("stage3 should pass");
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
        &format!("{dir}/tests/onshape/chamfered_cube_10_c1_medium.step"),
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

// --- Missing output path test ---

#[test]
fn missing_output_path_error() {
    let dir = manifest_dir();
    let stl = format!("{dir}/tests/ccad/generated/cube.stl");
    let config = config::parse_config_from(["brepper", &stl, "-q"]).unwrap();
    let mesh = stage1::stage1(&config).expect("stage1 should pass");
    let surfaces = stage2::stage2(&config, mesh).expect("stage2 should pass");
    let brep = stage3::stage3(&config, surfaces).expect("stage3 should pass");
    let result = stage4::stage4(&config, brep);
    assert!(result.is_err());
    let err = result.unwrap_err();
    assert!(
        err.to_string().contains("no output STEP file path"),
        "expected MissingOutputPath error, got: {err}"
    );
}
