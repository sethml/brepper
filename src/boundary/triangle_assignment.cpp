#include "triangle_assignment.hpp"
#include "common/logging.hpp"
#include <pcl/conversions.h>
#include <limits>
#include <cmath>

namespace brepper {

TriangleAssigner::TriangleAssigner(const Config& config) : config_(config) {}

bool TriangleAssigner::assign(
    const pcl::PolygonMesh& mesh,
    const std::vector<FittedSurface>& surfaces,
    std::vector<TriangleAssignment>& assignments
) {
    LOG_DEBUG("Assigning triangles to surfaces");
    
    if (surfaces.empty()) {
        LOG_WARN("No surfaces to assign triangles to");
        return true;
    }
    
    // Extract vertices from mesh cloud
    pcl::PointCloud<pcl::PointXYZ> vertices;
    pcl::fromPCLPointCloud2(mesh.cloud, vertices);
    
    assignments.clear();
    assignments.reserve(mesh.polygons.size());
    
    int assigned_count = 0;
    int unassigned_count = 0;
    
    for (size_t tri_idx = 0; tri_idx < mesh.polygons.size(); ++tri_idx) {
        const auto& polygon = mesh.polygons[tri_idx];
        
        if (polygon.vertices.size() < 3) {
            continue;
        }
        
        // Get triangle vertices
        const auto& v0 = vertices[polygon.vertices[0]];
        const auto& v1 = vertices[polygon.vertices[1]];
        const auto& v2 = vertices[polygon.vertices[2]];
        
        // Compute centroid
        Eigen::Vector3f centroid(
            (v0.x + v1.x + v2.x) / 3.0f,
            (v0.y + v1.y + v2.y) / 3.0f,
            (v0.z + v1.z + v2.z) / 3.0f
        );
        
        // Compute normal (cross product of edges)
        Eigen::Vector3f e1(v1.x - v0.x, v1.y - v0.y, v1.z - v0.z);
        Eigen::Vector3f e2(v2.x - v0.x, v2.y - v0.y, v2.z - v0.z);
        Eigen::Vector3f normal = e1.cross(e2);
        float normal_len = normal.norm();
        if (normal_len > 1e-6f) {
            normal /= normal_len;
        } else {
            // Degenerate triangle
            normal = Eigen::Vector3f(0, 0, 1);
        }
        
        // Find best matching surface
        int best_surface_id = -1;
        double best_distance = std::numeric_limits<double>::max();
        double best_angle_deviation = 180.0;
        
        for (const auto& surface : surfaces) {
            double distance = compute_distance(centroid, normal, surface);
            
            // Compute normal deviation
            Eigen::Vector3f surface_normal = compute_surface_normal(centroid, surface);
            double dot = std::abs(normal.dot(surface_normal));
            dot = std::min(1.0f, std::max(-1.0f, (float)dot));  // Clamp to [-1, 1]
            double angle_deviation = std::acos(dot) * 180.0 / M_PI;
            
            // Check if this surface is better
            if (distance < config_.assignment_distance_threshold &&
                angle_deviation < config_.assignment_angle_threshold_degrees &&
                distance < best_distance) {
                best_distance = distance;
                best_surface_id = surface.surface_id;
                best_angle_deviation = angle_deviation;
            }
        }
        
        TriangleAssignment assignment;
        assignment.triangle_id = static_cast<int>(tri_idx);
        assignment.surface_id = best_surface_id;
        assignment.distance = best_distance;
        assignment.normal_deviation_degrees = best_angle_deviation;
        
        assignments.push_back(assignment);
        
        if (best_surface_id >= 0) {
            ++assigned_count;
        } else {
            ++unassigned_count;
        }
    }
    
    LOG_INFO("Triangle assignment: ", assigned_count, " assigned, ", 
             unassigned_count, " unassigned out of ", mesh.polygons.size(), " triangles");
    
    return true;
}

double TriangleAssigner::compute_distance(
    const Eigen::Vector3f& centroid,
    [[maybe_unused]] const Eigen::Vector3f& normal,
    const FittedSurface& surface
) {
    switch (surface.type) {
        case SurfaceType::PLANE:
            return compute_plane_distance(centroid, surface);
        case SurfaceType::CYLINDER:
            return compute_cylinder_distance(centroid, surface);
        case SurfaceType::SPHERE:
            return compute_sphere_distance(centroid, surface);
        case SurfaceType::CONE:
            return compute_cone_distance(centroid, surface);
        default:
            return std::numeric_limits<double>::max();
    }
}

double TriangleAssigner::compute_plane_distance(
    const Eigen::Vector3f& point,
    const FittedSurface& surface
) {
    // Plane coefficients: [nx, ny, nz, d] for ax + by + cz + d = 0
    if (surface.coefficients.size() < 4) {
        return std::numeric_limits<double>::max();
    }
    
    double a = surface.coefficients[0];
    double b = surface.coefficients[1];
    double c = surface.coefficients[2];
    double d = surface.coefficients[3];
    
    // Point-to-plane distance
    return std::abs(a * point.x() + b * point.y() + c * point.z() + d);
}

double TriangleAssigner::compute_cylinder_distance(
    const Eigen::Vector3f& point,
    const FittedSurface& surface
) {
    // Cylinder coefficients: [point_on_axis.x, .y, .z, axis.x, .y, .z, radius]
    if (surface.coefficients.size() < 7) {
        return std::numeric_limits<double>::max();
    }
    
    Eigen::Vector3d axis_point(surface.coefficients[0], surface.coefficients[1], surface.coefficients[2]);
    Eigen::Vector3d axis_dir(surface.coefficients[3], surface.coefficients[4], surface.coefficients[5]);
    axis_dir.normalize();
    double radius = surface.coefficients[6];
    
    Eigen::Vector3d p(point.x(), point.y(), point.z());
    Eigen::Vector3d v = p - axis_point;
    double dist_along_axis = v.dot(axis_dir);
    Eigen::Vector3d closest_on_axis = axis_point + dist_along_axis * axis_dir;
    double dist_to_axis = (p - closest_on_axis).norm();
    
    return std::abs(dist_to_axis - radius);
}

double TriangleAssigner::compute_sphere_distance(
    const Eigen::Vector3f& point,
    const FittedSurface& surface
) {
    // Sphere coefficients: [center.x, .y, .z, radius]
    if (surface.coefficients.size() < 4) {
        return std::numeric_limits<double>::max();
    }
    
    Eigen::Vector3d center(surface.coefficients[0], surface.coefficients[1], surface.coefficients[2]);
    double radius = surface.coefficients[3];
    
    Eigen::Vector3d p(point.x(), point.y(), point.z());
    double dist_to_center = (p - center).norm();
    
    return std::abs(dist_to_center - radius);
}

double TriangleAssigner::compute_cone_distance(
    const Eigen::Vector3f& point,
    const FittedSurface& surface
) {
    // Cone coefficients: [apex.x, .y, .z, axis.x, .y, .z, opening_angle]
    if (surface.coefficients.size() < 7) {
        return std::numeric_limits<double>::max();
    }
    
    Eigen::Vector3d apex(surface.coefficients[0], surface.coefficients[1], surface.coefficients[2]);
    Eigen::Vector3d axis_dir(surface.coefficients[3], surface.coefficients[4], surface.coefficients[5]);
    axis_dir.normalize();
    double opening_angle = surface.coefficients[6];
    
    Eigen::Vector3d p(point.x(), point.y(), point.z());
    Eigen::Vector3d v = p - apex;
    double dist_along_axis = v.dot(axis_dir);
    double expected_radius = std::abs(dist_along_axis) * std::tan(opening_angle);
    Eigen::Vector3d closest_on_axis = apex + dist_along_axis * axis_dir;
    double actual_radius = (p - closest_on_axis).norm();
    
    return std::abs(actual_radius - expected_radius);
}

Eigen::Vector3f TriangleAssigner::compute_surface_normal(
    const Eigen::Vector3f& point,
    const FittedSurface& surface
) {
    switch (surface.type) {
        case SurfaceType::PLANE: {
            // Plane normal is constant
            if (surface.coefficients.size() >= 3) {
                return Eigen::Vector3f(
                    static_cast<float>(surface.coefficients[0]),
                    static_cast<float>(surface.coefficients[1]),
                    static_cast<float>(surface.coefficients[2])
                ).normalized();
            }
            break;
        }
        
        case SurfaceType::CYLINDER: {
            // Cylinder normal points radially outward from axis
            if (surface.coefficients.size() >= 7) {
                Eigen::Vector3d axis_point(surface.coefficients[0], surface.coefficients[1], surface.coefficients[2]);
                Eigen::Vector3d axis_dir(surface.coefficients[3], surface.coefficients[4], surface.coefficients[5]);
                axis_dir.normalize();
                
                Eigen::Vector3d p(point.x(), point.y(), point.z());
                Eigen::Vector3d v = p - axis_point;
                double dist_along_axis = v.dot(axis_dir);
                Eigen::Vector3d closest_on_axis = axis_point + dist_along_axis * axis_dir;
                Eigen::Vector3d radial = p - closest_on_axis;
                if (radial.norm() > 1e-6) {
                    radial.normalize();
                    return radial.cast<float>();
                }
            }
            break;
        }
        
        case SurfaceType::SPHERE: {
            // Sphere normal points radially outward from center
            if (surface.coefficients.size() >= 3) {
                Eigen::Vector3d center(surface.coefficients[0], surface.coefficients[1], surface.coefficients[2]);
                Eigen::Vector3d p(point.x(), point.y(), point.z());
                Eigen::Vector3d radial = p - center;
                if (radial.norm() > 1e-6) {
                    radial.normalize();
                    return radial.cast<float>();
                }
            }
            break;
        }
        
        case SurfaceType::CONE: {
            // Cone normal is perpendicular to surface, pointing outward
            if (surface.coefficients.size() >= 7) {
                Eigen::Vector3d apex(surface.coefficients[0], surface.coefficients[1], surface.coefficients[2]);
                Eigen::Vector3d axis_dir(surface.coefficients[3], surface.coefficients[4], surface.coefficients[5]);
                axis_dir.normalize();
                double opening_angle = surface.coefficients[6];
                
                Eigen::Vector3d p(point.x(), point.y(), point.z());
                Eigen::Vector3d v = p - apex;
                double dist_along_axis = v.dot(axis_dir);
                Eigen::Vector3d closest_on_axis = apex + dist_along_axis * axis_dir;
                Eigen::Vector3d radial = p - closest_on_axis;
                if (radial.norm() > 1e-6) {
                    radial.normalize();
                    // Cone normal is tilted by (90 - opening_angle) from radial direction
                    // towards the axis direction (away from apex)
                    double cos_tilt = std::sin(opening_angle);
                    double sin_tilt = std::cos(opening_angle);
                    Eigen::Vector3d normal = radial * cos_tilt;
                    if (dist_along_axis >= 0) {
                        normal -= axis_dir * sin_tilt;
                    } else {
                        normal += axis_dir * sin_tilt;
                    }
                    normal.normalize();
                    return normal.cast<float>();
                }
            }
            break;
        }
        
        default:
            break;
    }
    
    // Fallback: return up vector
    return Eigen::Vector3f(0, 0, 1);
}

} // namespace brepper