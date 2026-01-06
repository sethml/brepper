#include "mesh_sampling.hpp"
#include "common/logging.hpp"

namespace brepper {

MeshSampler::MeshSampler(const Config& config) : config_(config) {}

bool MeshSampler::sample(const pcl::PolygonMesh& mesh, PointCloudNormalPtr& cloud) {
    LOG_DEBUG("Sampling mesh triangles");
    
    // TODO: Implement mesh sampling
    LOG_WARN("Mesh sampling not implemented yet");
    return true;
}

void MeshSampler::sample_triangle(
    const Eigen::Vector3f& v0,
    const Eigen::Vector3f& v1, 
    const Eigen::Vector3f& v2,
    const Eigen::Vector3f& normal,
    int num_samples,
    PointCloudNormalPtr& cloud
) {
    // TODO: Implement triangle sampling using barycentric coordinates
}

} // namespace brepper