#pragma once

#include "common/types.hpp"
#include "common/config.hpp"

namespace brepper {

class TriangleAssigner {
public:
    explicit TriangleAssigner(const Config& config);
    
    // Assign each mesh triangle to a fitted surface
    bool assign(
        const pcl::PolygonMesh& mesh,
        const std::vector<FittedSurface>& surfaces,
        std::vector<TriangleAssignment>& assignments
    );
    
private:
    const Config& config_;
    
    // Compute distance from triangle to surface
    double compute_distance(
        const Eigen::Vector3f& centroid,
        const Eigen::Vector3f& normal,
        const FittedSurface& surface
    );
};

} // namespace brepper