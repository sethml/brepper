//! Stage 4: Output
//!
//! 4.1: Write constructed objects to a STEP file.

use crate::config::Config;
use crate::stage3::Stage3Output;
use opencascade_sys::{
    b_rep_extrema, b_rep_g_prop, extrema, g_prop, message,
    topo_ds, OwnedPtr,
};
use std::error::Error;
use std::fmt::{Display, Formatter};
use std::time::Instant;

// ---------------------------------------------------------------------------
// Errors
// ---------------------------------------------------------------------------

#[derive(Debug)]
pub enum Stage4Error {
    MissingOutputPath,
    TransferFailed(String),
    WriteFailed(String),
    Compare(Stage4CompareError),
}

impl Display for Stage4Error {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        match self {
            Stage4Error::MissingOutputPath => write!(f, "no output STEP file path specified (-o)"),
            Stage4Error::TransferFailed(msg) => write!(f, "stage 4.1 STEP transfer failed: {msg}"),
            Stage4Error::WriteFailed(msg) => write!(f, "stage 4.1 STEP write failed: {msg}"),
            Stage4Error::Compare(e) => write!(f, "stage 4.1 compare: {e}"),
        }
    }
}

impl Error for Stage4Error {}

impl From<Stage4CompareError> for Stage4Error {
    fn from(e: Stage4CompareError) -> Self {
        Stage4Error::Compare(e)
    }
}

#[derive(Debug)]
pub enum Stage4CompareError {
    /// The written STEP file could not be read back for validation.
    RoundTrip { path: String, error: String },
    /// The output's volume differs from the reference beyond 1% relative.
    Volume { rel_diff: f64 },
    /// The output shape is farther from the reference than --surface-tolerance.
    Distance { max_distance: f64, tolerance: f64 },
}

impl Display for Stage4CompareError {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        match self {
            Stage4CompareError::RoundTrip { path, error } => {
                write!(f, "failed to re-read written STEP file {path}: {error}")
            }
            Stage4CompareError::Volume { rel_diff } => {
                write!(
                    f,
                    "output volume differs from STEP by {rel_diff:.2e} relative (tolerance 1e-2)"
                )
            }
            Stage4CompareError::Distance { max_distance, tolerance } => {
                write!(
                    f,
                    "output STEP distance {max_distance:.6e} exceeds tolerance {tolerance:.6e}"
                )
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Stage 4 entry point
// ---------------------------------------------------------------------------

/// Run stage 4: write the reconstructed B-Rep to a STEP file.
pub fn stage4(config: &Config, input: Stage3Output) -> Result<(), Stage4Error> {
    let t = Instant::now();
    let output_path = config
        .output_step
        .as_deref()
        .ok_or(Stage4Error::MissingOutputPath)?;

    // Build a compound of all solids
    let compound = build_solid_compound(&input.solids);

    // Write STEP file (serialized via StepWriter's built-in mutex)
    crate::write_step_file(compound.as_shape(), output_path)
        .map_err(|e| {
            if e.contains("Transfer") {
                Stage4Error::TransferFailed(e)
            } else {
                Stage4Error::WriteFailed(e)
            }
        })?;

    if !config.quiet {
        eprintln!("Stage 4.1 ({:.3}s): Wrote {} solid(s) to {output_path}", t.elapsed().as_secs_f64(), input.solids.len());
    }

    // Compare against reference STEP if --compare
    if config.compare_shape.is_some() {
        compare_output_to_step(output_path, config)?;
    }

    Ok(())
}

/// Build a compound containing all solids.
fn build_solid_compound(solids: &[OwnedPtr<topo_ds::Solid>]) -> OwnedPtr<topo_ds::Compound> {
    let builder = topo_ds::Builder::new();
    let mut compound = topo_ds::Compound::new();
    builder.make_compound(&mut compound);
    for solid in solids {
        builder.add(compound.as_shape_mut(), solid.as_shape());
    }
    compound
}

/// Compare the written STEP file against reference STEP shape.
fn compare_output_to_step(
    output_path: &str,
    config: &Config,
) -> Result<(), Stage4CompareError> {
    let compare_shape = config.compare_shape.as_ref().unwrap();

    // Re-read the written STEP file to validate it round-trips correctly.
    // An unreadable output (including -o /dev/null) means the round-trip
    // cannot be validated, which fails --compare.
    let output_shape = match crate::read_step_file(output_path) {
        Ok(shape) => shape,
        Err(e) => {
            eprintln!("  Compare 4.1: ERROR failed to re-read written STEP file {output_path}: {e}");
            return Err(Stage4CompareError::RoundTrip {
                path: output_path.to_string(),
                error: e,
            });
        }
    };

    // Compare volumes
    let output_vol = compute_volume(&output_shape);
    let ref_vol = compute_volume(compare_shape);
    let rel_diff = if ref_vol > 0.0 {
        (output_vol - ref_vol).abs() / ref_vol
    } else {
        0.0
    };

    if !config.quiet {
        eprintln!(
            "  Compare 4.1: output volume {output_vol:.6} mm\u{00b3} vs STEP {ref_vol:.6} mm\u{00b3} (rel diff {rel_diff:.2e})",
        );
    }

    if rel_diff > 0.01 {
        eprintln!(
            "  Compare 4.1: ERROR volume differs by {rel_diff:.2e} (>1%)",
        );
        return Err(Stage4CompareError::Volume { rel_diff });
    }

    // Compute BRepExtrema distance between output and reference
    let progress = message::ProgressRange::new();
    let dist_calc = b_rep_extrema::DistShapeShape::new_shape2_extflag_extalgo_progressrange(
        &output_shape,
        compare_shape,
        0,
        extrema::ExtAlgo::Grad,
        &progress,
    );
    if dist_calc.is_done() && dist_calc.nb_solution() > 0 {
        let dist = dist_calc.value();
        if !config.quiet {
            eprintln!("  Compare 4.1: output distance to STEP = {dist:.6e} mm");
        }
        if dist > config.surface_tolerance_mm {
            return Err(Stage4CompareError::Distance {
                max_distance: dist,
                tolerance: config.surface_tolerance_mm,
            });
        }
    }

    Ok(())
}

fn compute_volume(shape: &topo_ds::Shape) -> f64 {
    let mut gprops = g_prop::GProps::new();
    b_rep_g_prop::volume_properties_shape_gprops_bool3(
        shape,
        &mut gprops,
        true, false, false,
    );
    gprops.mass().abs()
}
