// brepper - B-Rep from Mesh

pub mod config;
pub mod stage1;
pub mod stage2;
pub mod stage3;
pub mod stage4;
pub mod viz;

use opencascade_sys::{if_select, interface, message, step_control, topo_ds, OwnedPtr};
use std::sync::Mutex;

/// Global mutex to serialize OCCT's STEP Reader/Writer operations.
///
/// OCCT's STEP data exchange has process-global mutable state that is not
/// thread-safe: `STEPControl_Controller::Init()` (called by every Reader/Writer
/// constructor) writes to `static XSAlgo_AlgoContainer` and calls many
/// `Interface_Static::Init()` / `SetIVal()` on global settings. The STEP file
/// parser also writes to the global `Message_Messenger`.
///
/// Other OCCT APIs (geometry construction, BRep operations, ShapeFix, Sewing,
/// BRepExtrema, STL reading) operate only on owned/local data and are safe to
/// call concurrently.
static OCCT_STEP_MUTEX: Mutex<()> = Mutex::new(());

/// Read a STEP file and return the loaded shape.
///
/// Serializes access via `OCCT_STEP_MUTEX` to avoid data races in OCCT's
/// STEP reader global state.
pub fn read_step_file(path: &str) -> Result<OwnedPtr<topo_ds::Shape>, String> {
    let _lock = OCCT_STEP_MUTEX.lock().unwrap();
    let mut reader = step_control::Reader::new();
    let status = reader.read_file_charptr(path);
    if status != if_select::ReturnStatus::Retdone {
        return Err(format!("Failed to read STEP file {path}: {status:?}"));
    }
    let progress = message::ProgressRange::new();
    reader.transfer_roots(&progress);
    Ok(reader.one_shape())
}

/// Write a shape to a STEP file.
///
/// Serializes access via `OCCT_STEP_MUTEX` to avoid data races in OCCT's
/// STEP writer global state (including `Interface_Static` settings).
pub fn write_step_file(shape: &topo_ds::Shape, path: &str) -> Result<(), String> {
    let _lock = OCCT_STEP_MUTEX.lock().unwrap();
    interface::Static::set_c_val("write.step.schema", "AP214");
    interface::Static::set_c_val("write.step.product.name", "brepper");
    let mut writer = step_control::Writer::new();
    let progress = message::ProgressRange::new();
    let status = writer.transfer_shape_stepmodeltype_bool_progressrange(
        shape,
        step_control::StepModelType::Asis,
        true,
        &progress,
    );
    if status != if_select::ReturnStatus::Retdone {
        return Err(format!("STEPControl_Writer::Transfer returned {status:?}"));
    }
    let status = writer.write(path);
    if status != if_select::ReturnStatus::Retdone {
        return Err(format!("STEPControl_Writer::Write returned {status:?}"));
    }
    Ok(())
}

/// Suppress OCCT's default console output (transfer statistics, etc.).
/// Replaces the default stdout printer with one that only prints Fail-level messages.
pub fn suppress_occt_messages() {
    let mut messenger_handle = message::default_messenger();
    let messenger = messenger_handle.get_mut();
    // Remove the default PrinterOStream (which prints to stdout)
    let type_desc = message::PrinterOStream::get_type_descriptor();
    messenger.remove_printers(type_desc);
    // Add a printer that only shows failures
    let printer = message::PrinterOStream::new_gravity(message::Gravity::Fail);
    let printer_handle = message::PrinterOStream::to_handle(printer);
    let printer_base = printer_handle.to_handle_printer();
    messenger.add_printer(&printer_base);
}
