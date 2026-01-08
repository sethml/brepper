#include "mesh_sampling.hpp"
#include "common/logging.hpp"
#include <pcl/conversions.h>
#include <omp.h>
#include <cmath>
#include <vector>

namespace brepper {

MeshSampler::MeshSampler(const Config& config) 
    : config_(config)
    , rng_(config.random_seed >= 0 ? static_cast<unsigned int>(config.random_seed) : std::random_device{}()) 
{
    if (config.random_seed >= 0) {
        LOG_DEBUG("Using fixed random seed: ", config.random_seed);
    }
}

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
    
    const size_t num_triangles = mesh.polygons.size();
    
    // First pass: compute triangle areas (sequential, fast)
    std::vector<double> triangle_areas(num_triangles);
    double total_area = 0.0;
    
    for (size_t tri_idx = 0; tri_idx < num_triangles; ++tri_idx) {
        const auto& polygon = mesh.polygons[tri_idx];
        if (polygon.vertices.size() != 3) {
            triangle_areas[tri_idx] = 0.0;
            continue;
        }
        
        const auto& pt0 = (*vertices)[polygon.vertices[0]];
        const auto& pt1 = (*vertices)[polygon.vertices[1]];
        const auto& pt2 = (*vertices)[polygon.vertices[2]];
        
        Eigen::Vector3f v0(pt0.x, pt0.y, pt0.z);
        Eigen::Vector3f v1(pt1.x, pt1.y, pt1.z);
        Eigen::Vector3f v2(pt2.x, pt2.y, pt2.z);
        
        double area = computeTriangleArea(v0, v1, v2);
        triangle_areas[tri_idx] = area;
        total_area += area;
    }
    
    avg_triangle_area_ = total_area / std::max(size_t(1), num_triangles);
    
    double max_dist = config_.max_point_distance_mm;
    LOG_DEBUG("Max point distance: ", max_dist);
    
    // Pre-compute sample counts per triangle
    std::vector<int> samples_per_triangle(num_triangles);
    size_t total_estimated_points = 0;
    
    for (size_t tri_idx = 0; tri_idx < num_triangles; ++tri_idx) {
        double area = triangle_areas[tri_idx];
        double char_length = std::sqrt(area);
        int samples_per_length = static_cast<int>(std::ceil(char_length / max_dist));
        int num_samples = samples_per_length * samples_per_length;
        num_samples = std::max(num_samples, config_.min_samples_per_triangle);
        samples_per_triangle[tri_idx] = num_samples;
        total_estimated_points += num_samples;
    }
    
    // Check if we would exceed the maximum allowed samples
    if (total_estimated_points > config_.max_total_samples) {
        LOG_WARN("Estimated ", total_estimated_points, " samples exceeds max_total_samples (", 
                 config_.max_total_samples, "). Reducing sampling density.");
        
        // Scale down sample counts proportionally
        double scale_factor = static_cast<double>(config_.max_total_samples) / total_estimated_points;
        total_estimated_points = 0;
        for (size_t tri_idx = 0; tri_idx < num_triangles; ++tri_idx) {
            int scaled = static_cast<int>(samples_per_triangle[tri_idx] * scale_factor);
            samples_per_triangle[tri_idx] = std::max(scaled, config_.min_samples_per_triangle);
            total_estimated_points += samples_per_triangle[tri_idx];
        }
        LOG_DEBUG("Adjusted to ", total_estimated_points, " samples");
    }
    
    // Parallel sampling using per-thread clouds
    int num_threads = config_.num_threads > 0 ? config_.num_threads : omp_get_max_threads();
    omp_set_num_threads(num_threads);
    LOG_DEBUG("Using ", num_threads, " OpenMP threads for sampling");
    
    std::vector<PointCloudNormal> thread_clouds(num_threads);
    for (auto& tc : thread_clouds) {
        tc.reserve(total_estimated_points / num_threads + 1000);
    }
    
    #pragma omp parallel
    {
        int thread_id = omp_get_thread_num();
        PointCloudNormal& local_cloud = thread_clouds[thread_id];
        
        // Base seed for deterministic results
        unsigned int base_seed = (config_.random_seed >= 0) 
            ? static_cast<unsigned int>(config_.random_seed)
            : (std::random_device{}() + thread_id);  // Non-deterministic only when no seed set
        
        // Use static scheduling for deterministic results - each thread always
        // processes the same triangles, ensuring reproducible point generation
        #pragma omp for schedule(static)
        for (size_t tri_idx = 0; tri_idx < num_triangles; ++tri_idx) {
            const auto& polygon = mesh.polygons[tri_idx];
            if (polygon.vertices.size() != 3) continue;
            
            const auto& pt0 = (*vertices)[polygon.vertices[0]];
            const auto& pt1 = (*vertices)[polygon.vertices[1]];
            const auto& pt2 = (*vertices)[polygon.vertices[2]];
            
            Eigen::Vector3f v0(pt0.x, pt0.y, pt0.z);
            Eigen::Vector3f v1(pt1.x, pt1.y, pt1.z);
            Eigen::Vector3f v2(pt2.x, pt2.y, pt2.z);
            
            Eigen::Vector3f normal = computeFaceNormal(v0, v1, v2);
            
            // Add centroid
            Eigen::Vector3f centroid = (v0 + v1 + v2) / 3.0f;
            pcl::PointNormal pn;
            pn.x = centroid.x(); pn.y = centroid.y(); pn.z = centroid.z();
            pn.normal_x = normal.x(); pn.normal_y = normal.y(); pn.normal_z = normal.z();
            pn.curvature = 0.0f;
            local_cloud.push_back(pn);
            
            // Add vertices
            for (int i = 0; i < 3; ++i) {
                const auto& pt = (*vertices)[polygon.vertices[i]];
                pcl::PointNormal vpn;
                vpn.x = pt.x; vpn.y = pt.y; vpn.z = pt.z;
                vpn.normal_x = normal.x(); vpn.normal_y = normal.y(); vpn.normal_z = normal.z();
                vpn.curvature = 0.0f;
                local_cloud.push_back(vpn);
            }
            
            // Random samples - seed based on triangle index for thread-count independence
            int num_additional = samples_per_triangle[tri_idx] - 4;
            if (num_additional > 0) {
                // Each triangle gets its own deterministic seed based on triangle index
                std::mt19937 tri_rng(base_seed + static_cast<unsigned int>(tri_idx));
                std::uniform_real_distribution<float> dist(0.0f, 1.0f);
                
                for (int i = 0; i < num_additional; ++i) {
                    float r1 = dist(tri_rng);
                    float r2 = dist(tri_rng);
                    if (r1 + r2 > 1.0f) { r1 = 1.0f - r1; r2 = 1.0f - r2; }
                    float r3 = 1.0f - r1 - r2;
                    
                    Eigen::Vector3f pt = r1 * v0 + r2 * v1 + r3 * v2;
                    pcl::PointNormal spn;
                    spn.x = pt.x(); spn.y = pt.y(); spn.z = pt.z();
                    spn.normal_x = normal.x(); spn.normal_y = normal.y(); spn.normal_z = normal.z();
                    spn.curvature = 0.0f;
                    local_cloud.push_back(spn);
                }
            }
        }
    }
    
    // Merge thread-local clouds
    cloud.reset(new PointCloudNormal);
    size_t total_points = 0;
    for (const auto& tc : thread_clouds) {
        total_points += tc.size();
    }
    cloud->reserve(total_points);
    
    for (const auto& tc : thread_clouds) {
        cloud->insert(cloud->end(), tc.begin(), tc.end());
    }
    
    num_sampled_points_ = cloud->size();
    cloud->width = num_sampled_points_;
    cloud->height = 1;
    cloud->is_dense = true;
    
    LOG_INFO("Sampled ", num_sampled_points_, " points from ", num_triangles, " triangles");
    LOG_DEBUG("Average triangle area: ", avg_triangle_area_);
    
    return true;
}

} // namespace brepper