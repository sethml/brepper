#include "triangle_assignment.hpp"
#include "common/logging.hpp"

namespace brepper {

TriangleAssigner::TriangleAssigner(const Config& config) : config_(config) {}

bool TriangleAssigner::assign(
    const pcl::PolygonMesh& mesh,
    const std::vector<FittedSurface>& surfaces,
    std::vector<TriangleAssignment>& assignments
) {
    LOG_DEBUG("Assigning triangles to surfaces");
    
    // TODO: Implement triangle assignment
    LOG_WARN("Triangle assignment not implemented yet");
    return true;
}

double TriangleAssigner::compute_distance(
    const Eigen::Vector3f& centroid,
    const Eigen::Vector3f& normal,
    const FittedSurface& surface
) {
    // TODO: Implement distance computation for different surface types
    return 0.0;
}

} // namespace brepper