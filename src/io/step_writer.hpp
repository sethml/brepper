#pragma once

#include "common/config.hpp"
#include <TopoDS_Shape.hxx>

namespace brepper {

class STEPWriter {
public:
    explicit STEPWriter(const Config& config);
    
    // Export OCCT shape to STEP file
    bool write(const TopoDS_Shape& shape, const std::string& filename);
    
private:
    const Config& config_;
};

} // namespace brepper