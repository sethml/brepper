//! brepper — Convert STL mesh to STEP with fitted analytic and freeform surfaces.

use brepper::config::{self, Config};
use brepper::{stage1, stage2, stage3, stage4};
use std::env;
use std::process;

fn run(config: &Config) -> Result<(), Box<dyn std::error::Error>> {
    // Stage 1: Mesh Input & Preprocessing
    let mesh = stage1::stage1(config)?;

    if !config.stage.at_least(2, 1) {
        return Ok(());
    }

    // Stage 2: Surface Fitting
    let surfaces = stage2::stage2(config, mesh)?;

    if !config.stage.at_least(3, 1) {
        return Ok(());
    }

    // Stage 3: Surface Reconstruction
    let brep = stage3::stage3(config, surfaces)?;

    if !config.stage.at_least(4, 1) {
        return Ok(());
    }

    // Stage 4: Output
    stage4::stage4(config, brep)?;

    Ok(())
}

fn main() {
    let args: Vec<String> = env::args().collect();
    let program = args[0].clone();

    let config = match config::parse_args(args.into_iter().skip(1)) {
        Ok(c) => c,
        Err(e) => {
            eprintln!("Error: {e}");
            eprintln!();
            config::print_usage(&program);
            process::exit(1);
        }
    };

    if let Err(e) = run(&config) {
        eprintln!("Error: {e}");
        process::exit(1);
    }
}
