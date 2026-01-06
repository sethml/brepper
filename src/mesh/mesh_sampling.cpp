#include "mesh_sampling.hpp"
#include "common/logging.hpp"
#include <pcl/conversions.h>
#include <cmath>

namespace brepper {

MeshSampler::MeshSampler(const Config& config) 
    : config_(config)
    , rng_(std::random_device{}()) 
{}

Eigen::Vector3f MeshSampler::computeFaceNormal(
    const Eigen::Vector3f& v0,
    const Eigen::Vector3f& v1,
    const Eigen::Vector3f& v2
) const {
    Eigen::Vector3f edge1 = v1 - v0;
    Eigen::Vector3f edge2 = v2 - v0;
    Eigen::Vector3f normal = edge1.cross(edge2);
    
    float len = normal.norm();
    if (len > 1e-10f) {
        normal /= len;
    } else {
        normal = Eigen::Vector3f(0, 0, 1);  // Degenerate triangle fallback
    }
    return normal;
}

double MeshSampler::computeTriangleArea(
    const Eigen::Vector3f& v0,
    const Eigen::Vector3f& v1,
    const Eigen::Vector3f& v2
) const {
    Eigen::Vector3f edge1 = v1 - v0;
    Eigen::Vector3f edge2 = v2 - v0;
    return 0.5 * edge1.cross(edge2).norm();
}

Eigen::Vector3f MeshSampler::randomPointInTriangle(
    const Eigen::Vector3f& v0,
    const Eigen::Vector3f& v1,
    const Eigen::Vector3f& v2
) {
    // Generate random barycentric coordinates
    std::uniform_real_distribution<float> dist(0.0f, 1.0f);
    float r1 = dist(rng_);
    float r2 = dist(rng_);
    
    // Ensure point is inside triangle (fold if outside)
    if (r1 + r2 > 1.0f) {
        r1 = 1.0f - r1;
        r2 = 1.0f - r2;
    }
    
    float r3 = 1.0f - r1 - r2;
    
    return r1 * v0 + r2 * v1 + r3 * v2;
}

void MeshSampler::sampleTriangle(
    const Eigen::Vector3f& v0,
    const Eigen::Vector3f& v1, 
    const Eigen::Vector3f& v2,
    const Eigen::Vector3f& normal,
    int num_samples,
    PointCloudNormal& cloud
) {
    for (int i = 0; i < num_samples; ++i) {
        Eigen::Vector3f pt = randomPointInTriangle(v0, v1, v2);
        
        pcl::PointNormal pn;
        pn.x = pt.x();
        pn.y = pt.y();
        pn.z = pt.z();
        pn.normal_x = normal.x();
        pn.normal_y = normal.y();
        pn.normal_z = normal.z();
        pn.curvature = 0.0f;  // Will be computed later if needed
        
        cloud.push_back(pn);
    }
}

bool MeshSampler::sample(const pcl::PolygonMesh& mesh, PointCloudNormalPtr& cloud) {
    LOG_DEBUG("Sampling mesh triangles");
    
    // Convert mesh cloud to PointCloud<PointXYZ>
    pcl::PointCloud<pcl::PointXYZ>::Ptr vertices(new pcl::PointCloud<pcl::PointXYZ>);
    pcl::fromPCLPointCloud2(mesh.cloud, *vertices);
    
    if (vertices->empty()) {
        LOG_ERROR("Mesh has no vertices");
        return false;
    }
    
    // Initialize output cloud
    cloud.reset(new PointCloudNormal);
    cloud->reserve(mesh.polygons.size() * (config_.min_samples_per_triangle + 1));
    
    // First pass: compute total area and average triangle area
    double total_area = 0.0;
    std::vector<double> triangle_areas;
    triangle_areas.reserve(mesh.polygons.size());
    
    for (const auto& polygon : mesh.polygons) {
        if (polygon.vertices.size() != 3) continue;
        
        const auto& pt0 = (*vertices)[polygon.vertices[0]];
        const auto& pt1 = (*vertices)[polygon.vertices[1]];
        const auto& pt2 = (*vertices)[polygon.vertices[2]];
        
        Eigen::Vector3f v0(pt0.x, pt0.y, pt0.z);
        Eigen::Vector3f v1(pt1.x, pt1.y, pt1.z);
        Eigen::Vector3f v2(pt2.x, pt2.y, pt2.z);
        
        double area = computeTriangleArea(v0, v1, v2);
        triangle_areas.push_back(area);
        total_area += area;
    }
    
    avg_triangle_area_ = total_area / std::max(size_t(1), mesh.polygons.size());
    
    // Determine sampling density
    double sample_density = config_.sample_density;
    if (sample_density <= 0.0) {
        // Auto-compute: aim for roughly uniform point density
        // Target approximately 10 points per average triangle
        sample_density = 10.0 / std::max(avg_triangle_area_, 1e-10);
        LOG_DEBUG("Auto-computed sample density: ", sample_density, " points/unit²");
    }
    
    // Second pass: sample each triangle
    size_t tri_idx = 0;
    for (const auto& polygon : mesh.polygons) {
        if (polygon.vertices.size() != 3) continue;
        
        const auto& pt0 = (*vertices)[polygon.vertices[0]];
        const auto& pt1 = (*vertices)[polygon.vertices[1]];
        const auto& pt2 = (*vertices)[polygon.vertices[2]];
        
        Eigen::Vector3f v0(pt0.x, pt0.y, pt0.z);
        Eigen::Vector3f v1(pt1.x, pt1.y, pt1.z);
        Eigen::Vector3f v2(pt2.x, pt2.y, pt2.z);
        
        // Compute face normal
        Eigen::Vector3f normal = computeFaceNormal(v0, v1, v2);
        
        // Always add triangle centroid with face normal
        Eigen::Vector3f centroid = (v0 + v1 + v2) / 3.0f;
        pcl::PointNormal pn;
        pn.x = centroid.x();
        pn.y = centroid.y();
        pn.z = centroid.z();
        pn.normal_x = normal.x();
        pn.normal_y = normal.y();
        pn.normal_z = normal.z();
        pn.curvature = 0.0f;
        cloud->push_back(pn);
        
        // Add vertices with face normal
        for (int i = 0; i < 3; ++i) {
            const auto& pt = (*vertices)[polygon.vertices[i]];
            pcl::PointNormal vpn;
            vpn.x = pt.x;
            vpn.y = pt.y;
            vpn.z = pt.z;
            vpn.normal_x = normal.x();
            vpn.normal_y = normal.y();
            vpn.normal_z = normal.z();
            vpn.curvature = 0.0f;
            cloud->push_back(vpn);
        }
        
        // Compute additional samples based on area
        double area = triangle_areas[tri_idx];
        int num_additional = static_cast<int>(area * sample_density);
        num_additional = std::max(num_additional, config_.min_samples_per_triangle - 4);
        
        if (num_additional > 0) {
            sampleTriangle(v0, v1, v2, normal, num_additional, *cloud);
        }
        
        ++tri_idx;
    }
    
    // Remove duplicate points (optional, can be expensive)
    // For now we skip this to preserve all samples
    
    num_sampled_points_ = cloud->size();
    cloud->width = num_sampled_points_;
    cloud->height = 1;
    cloud->is_dense = true;
    
    LOG_INFO("Sampled ", num_sampled_points_, " points from ", mesh.polygons.size(), " triangles");
    LOG_DEBUG("Average triangle area: ", avg_triangle_area_);
    
    return true;
}

} // namespace brepper