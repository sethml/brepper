#pragma once

#include "common/types.hpp"
#include "common/config.hpp"
#include <pcl/PolygonMesh.h>
#include <random>

namespace brepper {

class MeshSampler {
public:
    explicit MeshSampler(const Config& config);
    
    // Sample points uniformly from mesh triangles and compute normals
    // Returns a point cloud with XYZ + Normal data
    bool sample(const pcl::PolygonMesh& mesh, PointCloudNormalPtr& cloud);
    
    // Get statistics after sampling
    size_t getNumSampledPoints() const { return num_sampled_points_; }
    double getAverageTriangleArea() const { return avg_triangle_area_; }
    
private:
    const Config& config_;
    std::mt19937 rng_;
    
    size_t num_sampled_points_ = 0;
    double avg_triangle_area_ = 0.0;
    
    // Compute face normal from 3 vertices
    Eigen::Vector3f computeFaceNormal(
        const Eigen::Vector3f& v0,
        const Eigen::Vector3f& v1,
        const Eigen::Vector3f& v2
    ) const;
    
    // Compute triangle area
    double computeTriangleArea(
        const Eigen::Vector3f& v0,
        const Eigen::Vector3f& v1,
        const Eigen::Vector3f& v2
    ) const;
    
    // Sample points within a single triangle using barycentric coordinates
    void sampleTriangle(
        const Eigen::Vector3f& v0,
        const Eigen::Vector3f& v1, 
        const Eigen::Vector3f& v2,
        const Eigen::Vector3f& normal,
        int num_samples,
        PointCloudNormal& cloud
    );
    
    // Generate random point in triangle using barycentric coords
    Eigen::Vector3f randomPointInTriangle(
        const Eigen::Vector3f& v0,
        const Eigen::Vector3f& v1,
        const Eigen::Vector3f& v2
    );
};

} // namespace brepper