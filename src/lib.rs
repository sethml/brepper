// brepper - B-Rep from Mesh

pub mod config;
pub mod stage1;
pub mod stage2;
pub mod stage3;
pub mod stage4;
pub mod viz;

use opencascade_sys::message;

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
