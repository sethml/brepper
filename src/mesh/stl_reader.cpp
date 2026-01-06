#include "stl_reader.hpp"
#include "common/logging.hpp"
#include <pcl/io/vtk_lib_io.h>
#include <pcl/conversions.h>
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
    
    // Convert coordinates to mm (OpenCASCADE works internally in mm)
    if (config_.stl_units != Units::Millimeters) {
        double scale = units_to_mm(config_.stl_units);
        LOG_DEBUG("Converting from ", 
                  (config_.stl_units == Units::Meters ? "meters" :
                   config_.stl_units == Units::Centimeters ? "centimeters" : "inches"),
                  " to mm (scale: ", scale, ")");
        
        pcl::PointCloud<pcl::PointXYZ>::Ptr vertices(new pcl::PointCloud<pcl::PointXYZ>);
        pcl::fromPCLPointCloud2(mesh.cloud, *vertices);
        
        for (auto& pt : *vertices) {
            pt.x *= scale;
            pt.y *= scale;
            pt.z *= scale;
        }
        
        pcl::toPCLPointCloud2(*vertices, mesh.cloud);
    }
    
    LOG_INFO("Loaded STL: ", num_vertices_, " vertices, ", num_triangles_, " triangles");
    
    return true;
}

} // namespace brepper