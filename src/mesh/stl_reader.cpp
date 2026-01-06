#include "stl_reader.hpp"
#include "common/logging.hpp"
#include <pcl/io/vtk_lib_io.h>
#include <filesystem>

namespace brepper {

STLReader::STLReader(const Config& config) : config_(config) {}

bool STLReader::load(const std::string& filename, pcl::PolygonMesh& mesh) {
    LOG_DEBUG("Loading STL file: ", filename);
    
    // Check file exists
    if (!std::filesystem::exists(filename)) {
        LOG_ERROR("STL file not found: ", filename);
        return false;
    }
    
    // Load using PCL's STL reader
    int result = pcl::io::loadPolygonFileSTL(filename, mesh);
    
    if (result < 0) {
        LOG_ERROR("Failed to load STL file: ", filename);
        return false;
    }
    
    // Extract statistics
    num_vertices_ = result;  // loadPolygonFileSTL returns number of vertices
    num_triangles_ = mesh.polygons.size();
    
    // Validate mesh
    if (num_triangles_ == 0) {
        LOG_ERROR("STL file contains no triangles");
        return false;
    }
    
    // Check all polygons are triangles
    for (const auto& polygon : mesh.polygons) {
        if (polygon.vertices.size() != 3) {
            LOG_ERROR("STL file contains non-triangular faces");
            return false;
        }
    }
    
    LOG_INFO("Loaded STL: ", num_vertices_, " vertices, ", num_triangles_, " triangles");
    
    return true;
}

} // namespace brepper