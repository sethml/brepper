#pragma once

#include "common/types.hpp"
#include "common/config.hpp"

namespace brepper {

class MeshSampler {
public:
    explicit MeshSampler(const Config& config);
    
    // Sample points uniformly from mesh triangles
    bool sample(const pcl::PolygonMesh& mesh, PointCloudNormalPtr& cloud);
    
private:
    const Config& config_;
    
    // Sample points within a single triangle
    void sample_triangle(
        const Eigen::Vector3f& v0,
        const Eigen::Vector3f& v1, 
        const Eigen::Vector3f& v2,
        const Eigen::Vector3f& normal,
        int num_samples,
        PointCloudNormalPtr& cloud
    );
};

} // namespace brepper