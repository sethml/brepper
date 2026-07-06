//! Generates test models that cannot be produced via ccad, using OCCT
//! directly. The STL/STEP pairs are written to tests/occt/ and committed to
//! the repo, so this tool only needs to be re-run when a model changes.
//!
//! Usage: gen_test_models [<output_dir>]   (default: tests/occt)
//!
//! Models:
//! - hollow_cube: 20 mm cube with a fully enclosed 10 mm cubic void at its
//!   center. Exercises multi-shell solid construction: the mesh has two
//!   closed shells and the STEP output must be ONE solid with an internal
//!   cavity (volume 7000 mm³), not two nested solids.
//! - full_torus: complete donut (R=10, r=3). Exercises the doubly-periodic
//!   seam handling in stage 3.4: a complete torus face needs both a U seam
//!   (around the major circle) and a V seam (around the tube).

use opencascade_sys::{b_rep_algo_api, b_rep_mesh, b_rep_prim_api, gp, message, stl_api, topo_ds};

fn write_pair(shape: &topo_ds::Shape, dir: &str, name: &str) {
    // Tessellate for STL (also used by STEP viewers, harmless for STEP).
    let _mesh = b_rep_mesh::IncrementalMesh::new_shape_real_bool_real_bool(
        shape, 0.1, false, 28.6_f64.to_radians(), false,
    );

    let stl_path = format!("{dir}/{name}.stl");
    let mut writer = stl_api::Writer::new();
    *writer.ascii_mode() = false;
    let progress = message::ProgressRange::new();
    if !writer.write(shape, &stl_path, &progress) {
        panic!("failed to write {stl_path}");
    }
    println!("wrote {stl_path}");

    let step_path = format!("{dir}/{name}.step");
    brepper::write_step_file(shape, &step_path).expect("failed to write STEP");
    println!("wrote {step_path}");
}

fn hollow_cube(dir: &str) {
    let mut outer = b_rep_prim_api::MakeBox::new_real3(20.0, 20.0, 20.0);
    let inner_corner = gp::Pnt::new_real3(5.0, 5.0, 5.0);
    let mut inner = b_rep_prim_api::MakeBox::new_pnt_real3(&inner_corner, 10.0, 10.0, 10.0);

    let progress = message::ProgressRange::new();
    let mut cut = b_rep_algo_api::Cut::new_shape2_progressrange(
        outer.shape(),
        inner.shape(),
        &progress,
    );
    assert!(cut.is_done(), "hollow_cube boolean cut failed");
    write_pair(cut.shape(), dir, "hollow_cube");
}

fn full_torus(dir: &str) {
    let mut torus = b_rep_prim_api::MakeTorus::new_real2(10.0, 3.0);
    write_pair(torus.shape(), dir, "full_torus");
}

fn main() {
    let dir = std::env::args().nth(1).unwrap_or_else(|| "tests/occt".to_string());
    std::fs::create_dir_all(&dir).expect("failed to create output dir");
    hollow_cube(&dir);
    full_torus(&dir);
}
