#pragma once

#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
#include <pcl/PolygonMesh.h>
#include <Eigen/Dense>
#include <vector>
#include <memory>

namespace brepper {

// Common point types
using Point = pcl::PointXYZ;
using PointNormal = pcl::PointNormal;
using PointCloud = pcl::PointCloud<Point>;
using PointCloudNormal = pcl::PointCloud<PointNormal>;
using PointCloudPtr = PointCloud::Ptr;
using PointCloudNormalPtr = PointCloudNormal::Ptr;

// Surface types that can be fitted
enum class SurfaceType {
    UNKNOWN = 0,
    PLANE,
    CYLINDER, 
    SPHERE,
    CONE,
    TORUS,
    NURBS
};

// Fitted surface representation
struct FittedSurface {
    SurfaceType type;
    std::vector<double> coefficients;  // Model-specific parameters
    PointCloudNormalPtr points;        // Points belonging to this surface
    std::vector<int> triangle_ids;     // Original mesh triangles
    double fitting_error;             // RMS fitting error
    int surface_id;                    // Unique ID
};

// Triangle assignment result
struct TriangleAssignment {
    int triangle_id;
    int surface_id;
    double distance;
    double normal_deviation_degrees;
};

// Boundary curve between two surfaces
struct BoundaryCurve {
    std::vector<Eigen::Vector3d> points;
    int surface_id_left;
    int surface_id_right;
    std::vector<int> edge_ids;  // Mesh edges forming this curve
};

// Processing results at each stage
struct ProcessingResults {
    // Stage 1: Mesh input
    pcl::PolygonMesh input_mesh;
    PointCloudNormalPtr sampled_cloud;
    
    // Stage 2: Surface segmentation  
    std::vector<FittedSurface> fitted_surfaces;
    PointCloudNormalPtr remaining_points;
    
    // Stage 3: Triangle assignment
    std::vector<TriangleAssignment> assignments;
    
    // Stage 4: Boundary detection
    std::vector<BoundaryCurve> boundary_curves;
};

} // namespace brepper