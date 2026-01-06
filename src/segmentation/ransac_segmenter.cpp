#include "ransac_segmenter.hpp"
#include "common/logging.hpp"
#include <pcl/segmentation/sac_segmentation.h>
#include <pcl/filters/extract_indices.h>
#include <pcl/sample_consensus/model_types.h>
#include <pcl/sample_consensus/method_types.h>
#include <pcl/console/print.h>

namespace brepper {

RANSACSegmenter::RANSACSegmenter(const Config& config) : config_(config) {}

bool RANSACSegmenter::segment(PointCloudNormalPtr cloud, std::vector<FittedSurface>& surfaces) {
    LOG_DEBUG("Starting RANSAC segmentation with ", cloud->size(), " points");
    
    // Suppress all PCL console output during RANSAC
    // These warnings are normal when fitting surface types that don't match the data
    pcl::console::setVerbosityLevel(pcl::console::L_ALWAYS);
    
    surfaces.clear();
    int surface_id = 0;
    
    // Work with a copy so we can remove points
    PointCloudNormalPtr remaining(new PointCloudNormal(*cloud));
    
    // Minimum points to continue segmentation
    const size_t min_remaining = static_cast<size_t>(
        std::max(config_.min_inliers, 
                 static_cast<int>(cloud->size() * config_.min_inlier_ratio)));
    
    LOG_DEBUG("Minimum remaining points threshold: ", min_remaining);
    
    while (remaining->size() > min_remaining) {
        LOG_DEBUG("Iteration ", surfaces.size() + 1, ": ", remaining->size(), " points remaining");
        
        // Try to fit each enabled surface type and keep the best one
        FittedSurface best_surface;
        best_surface.type = SurfaceType::UNKNOWN;
        best_surface.points.reset(new PointCloudNormal);
        pcl::PointIndices::Ptr best_inliers(new pcl::PointIndices);
        
        // Try planes first (most common)
        if (config_.fit_planes) {
            FittedSurface plane_result;
            pcl::PointIndices::Ptr plane_inliers(new pcl::PointIndices);
            if (fit_plane(remaining, plane_result, plane_inliers)) {
                if (plane_result.points->size() > best_surface.points->size()) {
                    best_surface = std::move(plane_result);
                    best_inliers = plane_inliers;
                }
            }
        }
        
        // Try cylinders
        if (config_.fit_cylinders) {
            FittedSurface cylinder_result;
            pcl::PointIndices::Ptr cylinder_inliers(new pcl::PointIndices);
            if (fit_cylinder(remaining, cylinder_result, cylinder_inliers)) {
                if (cylinder_result.points->size() > best_surface.points->size()) {
                    best_surface = std::move(cylinder_result);
                    best_inliers = cylinder_inliers;
                }
            }
        }
        
        // Try spheres
        if (config_.fit_spheres) {
            FittedSurface sphere_result;
            pcl::PointIndices::Ptr sphere_inliers(new pcl::PointIndices);
            if (fit_sphere(remaining, sphere_result, sphere_inliers)) {
                if (sphere_result.points->size() > best_surface.points->size()) {
                    best_surface = std::move(sphere_result);
                    best_inliers = sphere_inliers;
                }
            }
        }
        
        // Try cones
        if (config_.fit_cones) {
            FittedSurface cone_result;
            pcl::PointIndices::Ptr cone_inliers(new pcl::PointIndices);
            if (fit_cone(remaining, cone_result, cone_inliers)) {
                if (cone_result.points->size() > best_surface.points->size()) {
                    best_surface = std::move(cone_result);
                    best_inliers = cone_inliers;
                }
            }
        }
        
        // Check if we found a valid surface
        if (best_surface.type == SurfaceType::UNKNOWN || 
            best_surface.points->size() < static_cast<size_t>(config_.min_inliers)) {
            LOG_DEBUG("No more valid surfaces found");
            break;
        }
        
        // Assign surface ID and add to results
        best_surface.surface_id = surface_id++;
        
        LOG_INFO("Found ", surfaceTypeToString(best_surface.type), 
                 " surface #", best_surface.surface_id, 
                 " with ", best_surface.points->size(), " points",
                 ", fitting error: ", best_surface.fitting_error);
        
        surfaces.push_back(std::move(best_surface));
        
        // Remove inliers from remaining cloud
        pcl::ExtractIndices<PointNormal> extract;
        extract.setInputCloud(remaining);
        extract.setIndices(best_inliers);
        extract.setNegative(true);  // Keep non-inliers
        
        PointCloudNormalPtr filtered(new PointCloudNormal);
        extract.filter(*filtered);
        remaining = filtered;
    }
    
    LOG_INFO("RANSAC segmentation complete: found ", surfaces.size(), " surfaces, ",
             remaining->size(), " points remaining");
    
    // Restore default PCL verbosity
    pcl::console::setVerbosityLevel(pcl::console::L_INFO);
    
    return true;
}

bool RANSACSegmenter::fit_surface_type(
    PointCloudNormalPtr cloud,
    SurfaceType type,
    FittedSurface& surface
) {
    pcl::PointIndices::Ptr inliers(new pcl::PointIndices);
    switch (type) {
        case SurfaceType::PLANE:
            return fit_plane(cloud, surface, inliers);
        case SurfaceType::CYLINDER:
            return fit_cylinder(cloud, surface, inliers);
        case SurfaceType::SPHERE:
            return fit_sphere(cloud, surface, inliers);
        case SurfaceType::CONE:
            return fit_cone(cloud, surface, inliers);
        default:
            return false;
    }
}

bool RANSACSegmenter::fit_plane(PointCloudNormalPtr cloud, FittedSurface& surface,
                                 pcl::PointIndices::Ptr& out_inliers) {
    if (cloud->size() < 3) return false;
    
    // Setup SAC segmentation for planes with normals
    pcl::SACSegmentationFromNormals<PointNormal, PointNormal> seg;
    seg.setOptimizeCoefficients(true);
    seg.setModelType(pcl::SACMODEL_NORMAL_PLANE);
    seg.setMethodType(pcl::SAC_RANSAC);
    seg.setMaxIterations(config_.ransac_max_iterations);
    seg.setDistanceThreshold(config_.ransac_distance_threshold);
    seg.setNormalDistanceWeight(config_.normal_distance_weight);
    
    seg.setInputCloud(cloud);
    seg.setInputNormals(cloud);  // PointNormal has both XYZ and normals
    
    pcl::ModelCoefficients::Ptr coefficients(new pcl::ModelCoefficients);
    
    seg.segment(*out_inliers, *coefficients);
    
    if (out_inliers->indices.empty()) {
        return false;
    }
    
    // Extract inlier points
    surface.type = SurfaceType::PLANE;
    surface.points.reset(new PointCloudNormal);
    
    pcl::ExtractIndices<PointNormal> extract;
    extract.setInputCloud(cloud);
    extract.setIndices(out_inliers);
    extract.setNegative(false);
    extract.filter(*surface.points);
    
    // Store coefficients: [nx, ny, nz, d] for plane ax + by + cz + d = 0
    surface.coefficients.clear();
    for (const auto& c : coefficients->values) {
        surface.coefficients.push_back(c);
    }
    
    // Compute fitting error (RMS distance to plane)
    double sum_sq_dist = 0.0;
    for (const auto& pt : *surface.points) {
        double dist = std::abs(coefficients->values[0] * pt.x +
                               coefficients->values[1] * pt.y +
                               coefficients->values[2] * pt.z +
                               coefficients->values[3]);
        sum_sq_dist += dist * dist;
    }
    surface.fitting_error = std::sqrt(sum_sq_dist / surface.points->size());
    
    return true;
}

bool RANSACSegmenter::fit_cylinder(PointCloudNormalPtr cloud, FittedSurface& surface,
                                    pcl::PointIndices::Ptr& out_inliers) {
    if (cloud->size() < 10) return false;
    
    pcl::SACSegmentationFromNormals<PointNormal, PointNormal> seg;
    seg.setOptimizeCoefficients(true);
    seg.setModelType(pcl::SACMODEL_CYLINDER);
    seg.setMethodType(pcl::SAC_RANSAC);
    seg.setMaxIterations(config_.ransac_max_iterations);
    seg.setDistanceThreshold(config_.ransac_distance_threshold);
    seg.setNormalDistanceWeight(config_.normal_distance_weight);
    seg.setRadiusLimits(0.0, 1000.0);  // Allow any radius
    
    seg.setInputCloud(cloud);
    seg.setInputNormals(cloud);
    
    pcl::ModelCoefficients::Ptr coefficients(new pcl::ModelCoefficients);
    
    seg.segment(*out_inliers, *coefficients);
    
    if (out_inliers->indices.empty()) {
        return false;
    }
    
    surface.type = SurfaceType::CYLINDER;
    surface.points.reset(new PointCloudNormal);
    
    pcl::ExtractIndices<PointNormal> extract;
    extract.setInputCloud(cloud);
    extract.setIndices(out_inliers);
    extract.setNegative(false);
    extract.filter(*surface.points);
    
    // Cylinder coefficients: [point_on_axis.x, .y, .z, axis.x, .y, .z, radius]
    surface.coefficients.clear();
    for (const auto& c : coefficients->values) {
        surface.coefficients.push_back(c);
    }
    
    // Compute fitting error (RMS distance to cylinder surface)
    double sum_sq_dist = 0.0;
    Eigen::Vector3d axis_point(coefficients->values[0], coefficients->values[1], coefficients->values[2]);
    Eigen::Vector3d axis_dir(coefficients->values[3], coefficients->values[4], coefficients->values[5]);
    axis_dir.normalize();
    double radius = coefficients->values[6];
    
    for (const auto& pt : *surface.points) {
        Eigen::Vector3d p(pt.x, pt.y, pt.z);
        Eigen::Vector3d v = p - axis_point;
        double dist_along_axis = v.dot(axis_dir);
        Eigen::Vector3d closest_on_axis = axis_point + dist_along_axis * axis_dir;
        double dist_to_axis = (p - closest_on_axis).norm();
        double dist_to_surface = std::abs(dist_to_axis - radius);
        sum_sq_dist += dist_to_surface * dist_to_surface;
    }
    surface.fitting_error = std::sqrt(sum_sq_dist / surface.points->size());
    
    return true;
}

bool RANSACSegmenter::fit_sphere(PointCloudNormalPtr cloud, FittedSurface& surface,
                                  pcl::PointIndices::Ptr& out_inliers) {
    if (cloud->size() < 10) return false;
    
    pcl::SACSegmentationFromNormals<PointNormal, PointNormal> seg;
    seg.setOptimizeCoefficients(true);
    seg.setModelType(pcl::SACMODEL_NORMAL_SPHERE);
    seg.setMethodType(pcl::SAC_RANSAC);
    seg.setMaxIterations(config_.ransac_max_iterations);
    seg.setDistanceThreshold(config_.ransac_distance_threshold);
    seg.setNormalDistanceWeight(config_.normal_distance_weight);
    seg.setRadiusLimits(0.0, 1000.0);
    
    seg.setInputCloud(cloud);
    seg.setInputNormals(cloud);
    
    pcl::ModelCoefficients::Ptr coefficients(new pcl::ModelCoefficients);
    
    seg.segment(*out_inliers, *coefficients);
    
    if (out_inliers->indices.empty()) {
        return false;
    }
    
    surface.type = SurfaceType::SPHERE;
    surface.points.reset(new PointCloudNormal);
    
    pcl::ExtractIndices<PointNormal> extract;
    extract.setInputCloud(cloud);
    extract.setIndices(out_inliers);
    extract.setNegative(false);
    extract.filter(*surface.points);
    
    // Sphere coefficients: [center.x, .y, .z, radius]
    surface.coefficients.clear();
    for (const auto& c : coefficients->values) {
        surface.coefficients.push_back(c);
    }
    
    // Compute fitting error
    double sum_sq_dist = 0.0;
    Eigen::Vector3d center(coefficients->values[0], coefficients->values[1], coefficients->values[2]);
    double radius = coefficients->values[3];
    
    for (const auto& pt : *surface.points) {
        Eigen::Vector3d p(pt.x, pt.y, pt.z);
        double dist_to_center = (p - center).norm();
        double dist_to_surface = std::abs(dist_to_center - radius);
        sum_sq_dist += dist_to_surface * dist_to_surface;
    }
    surface.fitting_error = std::sqrt(sum_sq_dist / surface.points->size());
    
    return true;
}

bool RANSACSegmenter::fit_cone(PointCloudNormalPtr cloud, FittedSurface& surface,
                                pcl::PointIndices::Ptr& out_inliers) {
    if (cloud->size() < 10) return false;
    
    pcl::SACSegmentationFromNormals<PointNormal, PointNormal> seg;
    seg.setOptimizeCoefficients(true);
    seg.setModelType(pcl::SACMODEL_CONE);
    seg.setMethodType(pcl::SAC_RANSAC);
    seg.setMaxIterations(config_.ransac_max_iterations);
    seg.setDistanceThreshold(config_.ransac_distance_threshold);
    seg.setNormalDistanceWeight(config_.normal_distance_weight);
    seg.setMinMaxOpeningAngle(0.0, M_PI / 2.0);  // 0 to 90 degrees half-angle
    
    seg.setInputCloud(cloud);
    seg.setInputNormals(cloud);
    
    pcl::ModelCoefficients::Ptr coefficients(new pcl::ModelCoefficients);
    
    seg.segment(*out_inliers, *coefficients);
    
    if (out_inliers->indices.empty()) {
        return false;
    }
    
    surface.type = SurfaceType::CONE;
    surface.points.reset(new PointCloudNormal);
    
    pcl::ExtractIndices<PointNormal> extract;
    extract.setInputCloud(cloud);
    extract.setIndices(out_inliers);
    extract.setNegative(false);
    extract.filter(*surface.points);
    
    // Cone coefficients: [apex.x, .y, .z, axis.x, .y, .z, opening_angle]
    surface.coefficients.clear();
    for (const auto& c : coefficients->values) {
        surface.coefficients.push_back(c);
    }
    
    // Compute fitting error (approximate - distance to cone surface)
    double sum_sq_dist = 0.0;
    Eigen::Vector3d apex(coefficients->values[0], coefficients->values[1], coefficients->values[2]);
    Eigen::Vector3d axis_dir(coefficients->values[3], coefficients->values[4], coefficients->values[5]);
    axis_dir.normalize();
    double opening_angle = coefficients->values[6];
    
    for (const auto& pt : *surface.points) {
        Eigen::Vector3d p(pt.x, pt.y, pt.z);
        Eigen::Vector3d v = p - apex;
        double dist_along_axis = v.dot(axis_dir);
        double expected_radius = std::abs(dist_along_axis) * std::tan(opening_angle);
        Eigen::Vector3d closest_on_axis = apex + dist_along_axis * axis_dir;
        double actual_radius = (p - closest_on_axis).norm();
        double dist_to_surface = std::abs(actual_radius - expected_radius);
        sum_sq_dist += dist_to_surface * dist_to_surface;
    }
    surface.fitting_error = std::sqrt(sum_sq_dist / surface.points->size());
    
    return true;
}

std::string RANSACSegmenter::surfaceTypeToString(SurfaceType type) {
    switch (type) {
        case SurfaceType::PLANE: return "plane";
        case SurfaceType::CYLINDER: return "cylinder";
        case SurfaceType::SPHERE: return "sphere";
        case SurfaceType::CONE: return "cone";
        case SurfaceType::TORUS: return "torus";
        case SurfaceType::NURBS: return "NURBS";
        default: return "unknown";
    }
}

} // namespace brepper