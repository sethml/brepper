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
    
    // Surface-specific distance computations
    double compute_plane_distance(const Eigen::Vector3f& point, const FittedSurface& surface);
    double compute_cylinder_distance(const Eigen::Vector3f& point, const FittedSurface& surface);
    double compute_sphere_distance(const Eigen::Vector3f& point, const FittedSurface& surface);
    double compute_cone_distance(const Eigen::Vector3f& point, const FittedSurface& surface);
    
    // Compute expected surface normal at a given point
    Eigen::Vector3f compute_surface_normal(const Eigen::Vector3f& point, const FittedSurface& surface);
};

} // namespace brepper