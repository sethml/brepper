#include "step_writer.hpp"
#include "common/logging.hpp"

// OpenCASCADE STEP export headers
#include <STEPControl_Writer.hxx>
#include <IFSelect_ReturnStatus.hxx>
#include <Message.hxx>
#include <Message_PrinterOStream.hxx>

namespace brepper {

STEPWriter::STEPWriter(const Config& config) : config_(config) {}

bool STEPWriter::write(const TopoDS_Shape& shape, const std::string& filename) {
    LOG_DEBUG("Exporting STEP file: ", filename);
    
    try {
        // Suppress OCCT verbose "Statistics on Transfer" messages
        Message::DefaultMessenger()->RemovePrinters(STANDARD_TYPE(Message_PrinterOStream));
        
        STEPControl_Writer writer;
        
        // Transfer shape
        IFSelect_ReturnStatus status = writer.Transfer(shape, STEPControl_AsIs);
        if (status != IFSelect_RetDone) {
            LOG_ERROR("Failed to transfer shape to STEP writer");
            return false;
        }
        
        // Write file
        status = writer.Write(filename.c_str());
        if (status != IFSelect_RetDone) {
            LOG_ERROR("Failed to write STEP file");
            return false;
        }
        
        LOG_INFO("Successfully exported STEP file: ", filename);
        return true;
        
    } catch (const std::exception& e) {
        LOG_ERROR("Exception during STEP export: ", e.what());
        return false;
    }
}

} // namespace brepper