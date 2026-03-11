//! Stage 4: Output
//!
//! 4.1: Write constructed objects to a STEP file.

use crate::config::Config;
use crate::stage3::Stage3Output;
use opencascade_sys::{
    b_rep_extrema, b_rep_g_prop, extrema, g_prop, if_select, interface, message,
    step_control, topo_ds, OwnedPtr,
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
pub struct Stage4CompareError {
    pub check_type: &'static str,
    pub max_distance: f64,
    pub tolerance: f64,
}

impl Display for Stage4CompareError {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{} distance {:.6e} exceeds tolerance {:.6e}",
            self.check_type, self.max_distance, self.tolerance
        )
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

    // Set STEP metadata
    interface::Static::set_c_val("write.step.schema", "AP214");
    interface::Static::set_c_val("write.step.product.name", "brepper");

    // Create writer and transfer the compound
    let mut writer = step_control::Writer::new();
    let progress = message::ProgressRange::new();
    let status = writer.transfer_shape_stepmodeltype_bool_progressrange(
        compound.as_shape(),
        step_control::StepModelType::Asis,
        true,
        &progress,
    );
    if status != if_select::ReturnStatus::Retdone {
        return Err(Stage4Error::TransferFailed(format!(
            "STEPControl_Writer::Transfer returned {:?}",
            status
        )));
    }

    // Write to file
    let status = writer.write(output_path);
    if status != if_select::ReturnStatus::Retdone {
        return Err(Stage4Error::WriteFailed(format!(
            "STEPControl_Writer::Write returned {:?}",
            status
        )));
    }

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

    // Re-read the written STEP file to validate it round-trips correctly
    let mut reader = step_control::Reader::new();
    let read_status = reader.read_file_charptr(output_path);
    if read_status != if_select::ReturnStatus::Retdone {
        eprintln!("  Compare 4.1: WARNING failed to re-read written STEP file {output_path}: {read_status:?}");
        return Ok(());
    }
    let progress = message::ProgressRange::new();
    reader.transfer_roots(&progress);
    let output_shape = reader.one_shape();

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
            "  Compare 4.1: WARNING volume differs by {rel_diff:.2e} (>1%)",
        );
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
            return Err(Stage4CompareError {
                check_type: "output STEP",
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
