//! Stage 4: Output
//!
//! 4.1: Write constructed objects to a STEP file.

use crate::config::Config;
use crate::stage3::Stage3Output;
use std::error::Error;
use std::fmt::{Display, Formatter};

// ---------------------------------------------------------------------------
// Errors
// ---------------------------------------------------------------------------

#[derive(Debug)]
pub enum Stage4Error {
    // TODO: specific error variants as stage is implemented
    NotImplemented,
    MissingOutputPath,
}

impl Display for Stage4Error {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        match self {
            Stage4Error::NotImplemented => write!(f, "stage 4 not yet implemented"),
            Stage4Error::MissingOutputPath => write!(f, "no output STEP file path specified (-o)"),
        }
    }
}

impl Error for Stage4Error {}

// ---------------------------------------------------------------------------
// Stage 4 entry point
// ---------------------------------------------------------------------------

/// Run stage 4: write the reconstructed B-Rep to a STEP file.
pub fn stage4(config: &Config, _input: Stage3Output) -> Result<(), Stage4Error> {
    let _output_path = config
        .output_step
        .as_deref()
        .ok_or(Stage4Error::MissingOutputPath)?;

    // Stage 4.1: Write STEP file
    // TODO: Use STEPControl_Writer to export the solid(s)
    // TODO: Set STEP header metadata (author, organization)
    // TODO: Consider also supporting BREP format for debugging
    if !config.quiet {
        eprintln!("Stage 4.1: Write STEP file (not yet implemented)");
    }

    Err(Stage4Error::NotImplemented)
}
