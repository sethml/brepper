#pragma once

#include "common/types.hpp"
#include "common/config.hpp"
#include <pcl/PolygonMesh.h>

namespace brepper {

class STLReader {
public:
    explicit STLReader(const Config& config);
    
    // Load STL file into PCL mesh
    bool load(const std::string& filename, pcl::PolygonMesh& mesh);
    
    // Get mesh statistics after loading
    size_t getNumVertices() const { return num_vertices_; }
    size_t getNumTriangles() const { return num_triangles_; }
    
private:
    const Config& config_;
    size_t num_vertices_ = 0;
    size_t num_triangles_ = 0;
};

} // namespace brepper