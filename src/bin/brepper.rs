//! brepper — Convert STL mesh to STEP with fitted analytic and freeform surfaces.

use brepper::config::{self, Config};
use brepper::{stage1, stage2, stage3, stage4};
use brepper::viz::{self, VizSender};
use std::process;

fn run(config: &Config, viz: Option<&VizSender>) -> Result<(), Box<dyn std::error::Error>> {
    // Stage 1: Mesh Input & Preprocessing
    let mesh = stage1::stage1(config)?;

    if !config.stage.at_least(2, 1) {
        return Ok(());
    }

    // Stage 2: Surface Fitting
    let surfaces = stage2::stage2(config, mesh, viz)?;

    if !config.stage.at_least(3, 1) {
        return Ok(());
    }

    // Stage 3: Surface Reconstruction
    let brep = stage3::stage3(config, surfaces, viz)?;

    if !config.stage.at_least(4, 1) {
        return Ok(());
    }

    // Stage 4: Output
    stage4::stage4(config, brep)?;

    Ok(())
}

fn main() {
    let mut config = config::parse_config();

    if let Err(e) = config.load_compare_step() {
        eprintln!("Error loading comparison STEP file: {e}");
        process::exit(1);
    }

    if config.viz_stages.is_empty() {
        // No visualization: run directly
        if let Err(e) = run(&config, None) {
            eprintln!("Error: {e}");
            process::exit(1);
        }
    } else {
        // Visualization mode: build mesh data on main thread, spawn pipeline
        let mesh = match stage1::stage1(&config) {
            Ok(m) => m,
            Err(e) => {
                eprintln!("Error in stage 1: {e}");
                process::exit(1);
            }
        };

        let (mesh_data, face_indices, face_vc) = viz::meshdata_from_connected_mesh(&mesh);

        // Load compare STEP mesh for visualization (not the same as config.compare_shape)
        let compare_mesh = config.compare_step.as_ref().map(|p| {
            viz::load_step_meshdata(std::path::Path::new(p), 0.5, 0.5)
        });

        let setup = viz::VizSetup {
            base_mesh: mesh_data,
            compare_mesh,
            face_indices,
            face_vertex_counts: face_vc,
        };

        // Move config into the pipeline thread (re-run stage1 there)
        viz::run_viz_window(setup, move |viz_sender| {
            let mesh = match stage1::stage1(&config) {
                Ok(m) => m,
                Err(e) => {
                    eprintln!("Error in stage 1 (viz thread): {e}");
                    return;
                }
            };

            if !config.stage.at_least(2, 1) {
                return;
            }

            let surfaces = match stage2::stage2(&config, mesh, Some(&viz_sender)) {
                Ok(s) => s,
                Err(e) => {
                    eprintln!("Error: {e}");
                    return;
                }
            };

            if !config.stage.at_least(3, 1) {
                return;
            }

            let brep = match stage3::stage3(&config, surfaces, Some(&viz_sender)) {
                Ok(b) => b,
                Err(e) => {
                    eprintln!("Error: {e}");
                    return;
                }
            };

            if !config.stage.at_least(4, 1) {
                return;
            }

            if let Err(e) = stage4::stage4(&config, brep) {
                eprintln!("Error: {e}");
            }
        });
    }
}