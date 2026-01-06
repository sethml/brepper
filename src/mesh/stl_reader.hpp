#pragma once

#include "common/types.hpp"
#include "common/config.hpp"

namespace brepper {

class STLReader {
public:
    explicit STLReader(const Config& config);
    
    // Load STL file into PCL mesh
    bool load(const std::string& filename, pcl::PolygonMesh& mesh);
    
private:
    const Config& config_;
};

} // namespace brepper