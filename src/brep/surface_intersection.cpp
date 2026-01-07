#include "surface_intersection.hpp"
#include "common/logging.hpp"

#include <gp_Pln.hxx>
#include <gp_Ax3.hxx>
#include <gp_Dir.hxx>
#include <gp_Pnt.hxx>
#include <gp_Cylinder.hxx>
#include <gp_Sphere.hxx>
#include <gp_Cone.hxx>
#include <Geom_Plane.hxx>
#include <Geom_CylindricalSurface.hxx>
#include <Geom_SphericalSurface.hxx>
#include <Geom_ConicalSurface.hxx>
#include <GeomAPI_IntSS.hxx>
#include <GeomAPI_ProjectPointOnSurf.hxx>
#include <GeomAdaptor_Surface.hxx>
#include <Extrema_ExtSS.hxx>
#include <Extrema_POnSurf.hxx>
#include <Adaptor3d_Surface.hxx>
#include <Standard_Failure.hxx>
#include <cmath>

namespace brepper {

SurfaceIntersectionAnalyzer::SurfaceIntersectionAnalyzer(const Config& config) 
    : config_(config) {}

void SurfaceIntersectionAnalyzer::build_adjacency_graph(
    const std::vector<FittedSurface>& surfaces,
    const std::vector<BoundaryCurve>& boundaries)
{
    adjacent_pairs_.clear();
    expected_boundaries_.clear();
    fitted_surfaces_.clear();
    
    // Store fitted surfaces by ID
    for (const auto& surface : surfaces) {
        fitted_surfaces_[surface.surface_id] = surface;
    }
    
    LOG_DEBUG("Processing ", boundaries.size(), " boundary curves");
    
    // Build adjacency from boundary curves
    for (const auto& curve : boundaries) {
        LOG_DEBUG("  Boundary curve: surface_id_left=", curve.surface_id_left, 
                 " surface_id_right=", curve.surface_id_right,
                 " points=", curve.points.size());
        
        int id1 = std::min(curve.surface_id_left, curve.surface_id_right);
        int id2 = std::max(curve.surface_id_left, curve.surface_id_right);
        
        // Skip mesh boundaries (where one side is -1)
        if (id1 < 0 || id2 < 0) {
            LOG_DEBUG("    -> Skipping (mesh boundary, id<0)");
            continue;
        }
        
        auto pair = std::make_pair(id1, id2);
        adjacent_pairs_.insert(pair);
        
        // Store expected boundary points
        for (const auto& pt : curve.points) {
            expected_boundaries_[pair].emplace_back(pt.x(), pt.y(), pt.z());
        }
    }
    
    LOG_INFO("Built adjacency graph: ", adjacent_pairs_.size(), " adjacent surface pairs");
}

void SurfaceIntersectionAnalyzer::create_surfaces(const std::vector<FittedSurface>& surfaces)
{
    geom_surfaces_.clear();
    
    for (const auto& surface : surfaces) {
        Handle(Geom_Surface) geom_surface = create_surface(surface);
        if (!geom_surface.IsNull()) {
            geom_surfaces_[surface.surface_id] = geom_surface;
        } else {
            LOG_WARN("Failed to create OCCT surface for surface #", surface.surface_id);
        }
    }
    
    LOG_INFO("Created ", geom_surfaces_.size(), " OCCT surfaces");
}

Handle(Geom_Surface) SurfaceIntersectionAnalyzer::create_surface(const FittedSurface& surface)
{
    switch (surface.type) {
        case SurfaceType::PLANE:
            return create_plane(surface);
        case SurfaceType::CYLINDER:
            return create_cylinder(surface);
        case SurfaceType::SPHERE:
            return create_sphere(surface);
        case SurfaceType::CONE:
            return create_cone(surface);
        default:
            return nullptr;
    }
}

Handle(Geom_Surface) SurfaceIntersectionAnalyzer::create_plane(const FittedSurface& surface)
{
    // Plane coefficients from PCL: [a, b, c, d] for ax + by + cz + d = 0
    if (surface.coefficients.size() < 4) return nullptr;
    
    double a = surface.coefficients[0];
    double b = surface.coefficients[1];
    double c = surface.coefficients[2];
    double d = surface.coefficients[3];
    
    try {
        return new Geom_Plane(a, b, c, d);
    } catch (...) {
        return nullptr;
    }
}

Handle(Geom_Surface) SurfaceIntersectionAnalyzer::create_cylinder(const FittedSurface& surface)
{
    if (surface.coefficients.size() < 7) return nullptr;
    
    gp_Pnt axis_point(surface.coefficients[0], surface.coefficients[1], surface.coefficients[2]);
    gp_Dir axis_dir(surface.coefficients[3], surface.coefficients[4], surface.coefficients[5]);
    double radius = surface.coefficients[6];
    
    if (radius <= 0) return nullptr;
    
    // Create reference direction perpendicular to axis
    gp_Dir ref_dir;
    if (std::abs(axis_dir.X()) < 0.9) {
        ref_dir = gp_Dir(1, 0, 0).Crossed(axis_dir);
    } else {
        ref_dir = gp_Dir(0, 1, 0).Crossed(axis_dir);
    }
    
    gp_Ax3 ax3(axis_point, axis_dir, ref_dir);
    gp_Cylinder cylinder(ax3, radius);
    
    return new Geom_CylindricalSurface(cylinder);
}

Handle(Geom_Surface) SurfaceIntersectionAnalyzer::create_sphere(const FittedSurface& surface)
{
    if (surface.coefficients.size() < 4) return nullptr;
    
    gp_Pnt center(surface.coefficients[0], surface.coefficients[1], surface.coefficients[2]);
    double radius = surface.coefficients[3];
    
    if (radius <= 0) return nullptr;
    
    gp_Ax3 ax3(center, gp_Dir(0, 0, 1), gp_Dir(1, 0, 0));
    gp_Sphere sphere(ax3, radius);
    
    return new Geom_SphericalSurface(sphere);
}

Handle(Geom_Surface) SurfaceIntersectionAnalyzer::create_cone(const FittedSurface& surface)
{
    if (surface.coefficients.size() < 7) return nullptr;
    
    gp_Pnt apex(surface.coefficients[0], surface.coefficients[1], surface.coefficients[2]);
    gp_Dir axis_dir(surface.coefficients[3], surface.coefficients[4], surface.coefficients[5]);
    double half_angle = surface.coefficients[6];
    
    if (half_angle <= 0 || half_angle >= M_PI/2) return nullptr;
    
    // Create coordinate system with origin at apex and Z along the axis
    gp_Ax3 ax3(apex, axis_dir);
    
    // Geom_ConicalSurface needs a reference radius at the coordinate system origin.
    // Since PCL gives us the apex, we need a non-zero reference radius.
    // Use a small but reasonable reference radius.
    double ref_radius = 1.0;  // 1mm reference radius
    
    // If we have points, compute a better reference radius from the data
    if (surface.points && !surface.points->empty()) {
        double avg_dist = 0.0;
        for (const auto& pt : *surface.points) {
            Eigen::Vector3d p(pt.x, pt.y, pt.z);
            Eigen::Vector3d a(apex.X(), apex.Y(), apex.Z());
            Eigen::Vector3d d(axis_dir.X(), axis_dir.Y(), axis_dir.Z());
            avg_dist += std::abs((p - a).dot(d));
        }
        avg_dist /= surface.points->size();
        ref_radius = avg_dist * std::tan(half_angle);
        if (ref_radius < 1e-6) ref_radius = 1.0;
    }
    
    try {
        return new Geom_ConicalSurface(ax3, half_angle, ref_radius);
    } catch (...) {
        return nullptr;
    }
}

void SurfaceIntersectionAnalyzer::analyze_intersections()
{
    results_.clear();
    result_map_.clear();
    
    LOG_INFO("Analyzing ", adjacent_pairs_.size(), " surface-surface intersections...");
    
    int success_count = 0;
    int failure_count = 0;
    
    for (const auto& pair : adjacent_pairs_) {
        IntersectionResult result = try_intersection(pair.first, pair.second);
        results_.push_back(result);
        result_map_[pair] = result;
        
        if (result.success) {
            ++success_count;
        } else {
            ++failure_count;
        }
    }
    
    LOG_INFO("Intersection analysis complete: ", success_count, " succeeded, ", 
             failure_count, " failed");
}

IntersectionResult SurfaceIntersectionAnalyzer::try_intersection(int id1, int id2)
{
    IntersectionResult result;
    result.surface_id_1 = id1;
    result.surface_id_2 = id2;
    result.success = false;
    result.min_distance = std::numeric_limits<double>::max();
    
    // Get OCCT surfaces
    auto it1 = geom_surfaces_.find(id1);
    auto it2 = geom_surfaces_.find(id2);
    
    if (it1 == geom_surfaces_.end() || it2 == geom_surfaces_.end()) {
        LOG_WARN("Missing surface for intersection test: ", id1, " or ", id2);
        return result;
    }
    
    const Handle(Geom_Surface)& s1 = it1->second;
    const Handle(Geom_Surface)& s2 = it2->second;
    
    // Get expected boundary points
    auto pair = std::make_pair(std::min(id1, id2), std::max(id1, id2));
    auto boundary_it = expected_boundaries_.find(pair);
    if (boundary_it != expected_boundaries_.end()) {
        result.expected_points = boundary_it->second;
    }
    
    // Try intersection with tolerance
    double tol = config_.sewing_tolerance;  // Use sewing tolerance as intersection tolerance
    
    try {
        GeomAPI_IntSS intersector(s1, s2, tol);
        
        if (intersector.IsDone() && intersector.NbLines() > 0) {
            result.success = true;
            result.curves.reserve(intersector.NbLines());
            
            for (int i = 1; i <= intersector.NbLines(); ++i) {
                result.curves.push_back(intersector.Line(i));
            }
            
            LOG_DEBUG("Intersection success: surfaces ", id1, " and ", id2, 
                     " -> ", result.curves.size(), " curves");
        } else {
            // Intersection failed - compute minimum distance
            LOG_DEBUG("Intersection failed: surfaces ", id1, " and ", id2, 
                     " - computing extrema...");
            
            compute_extrema(s1, s2, result.expected_points, 
                          result.min_distance, result.closest_point_1, result.closest_point_2);
            
            LOG_DEBUG("  Min distance: ", result.min_distance, " mm");
        }
    } catch (const Standard_Failure& e) {
        LOG_WARN("Intersection exception for surfaces ", id1, " and ", id2, 
                ": ", e.GetMessageString());
        
        // Try to compute extrema anyway
        try {
            compute_extrema(s1, s2, result.expected_points,
                          result.min_distance, result.closest_point_1, result.closest_point_2);
        } catch (...) {
            // Extrema also failed
        }
    }
    
    return result;
}

void SurfaceIntersectionAnalyzer::compute_extrema(
    const Handle(Geom_Surface)& s1,
    const Handle(Geom_Surface)& s2,
    const std::vector<gp_Pnt>& hint_points,
    double& min_dist,
    gp_Pnt& closest_p1,
    gp_Pnt& closest_p2)
{
    min_dist = std::numeric_limits<double>::max();
    
    // First, try using hint points to find approximate closest points
    // Project hint points onto both surfaces and measure distances
    for (const auto& hint : hint_points) {
        try {
            GeomAPI_ProjectPointOnSurf proj1(hint, s1);
            GeomAPI_ProjectPointOnSurf proj2(hint, s2);
            
            if (proj1.IsDone() && proj1.NbPoints() > 0 &&
                proj2.IsDone() && proj2.NbPoints() > 0) {
                
                gp_Pnt p1 = proj1.NearestPoint();
                gp_Pnt p2 = proj2.NearestPoint();
                double dist = p1.Distance(p2);
                
                if (dist < min_dist) {
                    min_dist = dist;
                    closest_p1 = p1;
                    closest_p2 = p2;
                }
            }
        } catch (...) {
            // Continue with next hint point
        }
    }
    
    // Also try Extrema_ExtSS for a more thorough search
    // This requires Adaptor3d_Surface, which GeomAdaptor_Surface provides
    try {
        GeomAdaptor_Surface as1(s1);
        GeomAdaptor_Surface as2(s2);
        
        // Get parameter bounds
        double u1min, u1max, v1min, v1max;
        double u2min, u2max, v2min, v2max;
        
        s1->Bounds(u1min, u1max, v1min, v1max);
        s2->Bounds(u2min, u2max, v2min, v2max);
        
        // Clamp infinite bounds for numerical computation
        const double MAX_PARAM = 1000.0;
        u1min = std::max(u1min, -MAX_PARAM);
        u1max = std::min(u1max, MAX_PARAM);
        v1min = std::max(v1min, -MAX_PARAM);
        v1max = std::min(v1max, MAX_PARAM);
        u2min = std::max(u2min, -MAX_PARAM);
        u2max = std::min(u2max, MAX_PARAM);
        v2min = std::max(v2min, -MAX_PARAM);
        v2max = std::min(v2max, MAX_PARAM);
        
        // Use hint points to narrow down the search region
        if (!hint_points.empty()) {
            // Find approximate parameter ranges from hint points
            double u1_center = 0, v1_center = 0;
            double u2_center = 0, v2_center = 0;
            int count = 0;
            
            for (const auto& hint : hint_points) {
                try {
                    GeomAPI_ProjectPointOnSurf proj1(hint, s1);
                    GeomAPI_ProjectPointOnSurf proj2(hint, s2);
                    
                    if (proj1.IsDone() && proj1.NbPoints() > 0 &&
                        proj2.IsDone() && proj2.NbPoints() > 0) {
                        double u, v;
                        proj1.LowerDistanceParameters(u, v);
                        u1_center += u;
                        v1_center += v;
                        proj2.LowerDistanceParameters(u, v);
                        u2_center += u;
                        v2_center += v;
                        ++count;
                    }
                } catch (...) {}
            }
            
            if (count > 0) {
                u1_center /= count;
                v1_center /= count;
                u2_center /= count;
                v2_center /= count;
                
                // Search within a reasonable range around the hint area
                const double SEARCH_RANGE = 50.0;  // Adjust based on expected model size
                u1min = std::max(u1min, u1_center - SEARCH_RANGE);
                u1max = std::min(u1max, u1_center + SEARCH_RANGE);
                v1min = std::max(v1min, v1_center - SEARCH_RANGE);
                v1max = std::min(v1max, v1_center + SEARCH_RANGE);
                u2min = std::max(u2min, u2_center - SEARCH_RANGE);
                u2max = std::min(u2max, u2_center + SEARCH_RANGE);
                v2min = std::max(v2min, v2_center - SEARCH_RANGE);
                v2max = std::min(v2max, v2_center + SEARCH_RANGE);
            }
        }
        
        double tol = config_.sewing_tolerance;
        Extrema_ExtSS extrema(as1, as2, u1min, u1max, v1min, v1max,
                              u2min, u2max, v2min, v2max, tol, tol);
        
        if (extrema.IsDone() && extrema.NbExt() > 0) {
            // Find minimum
            for (int i = 1; i <= extrema.NbExt(); ++i) {
                double sq_dist = extrema.SquareDistance(i);
                double dist = std::sqrt(sq_dist);
                
                if (dist < min_dist) {
                    min_dist = dist;
                    Extrema_POnSurf p1, p2;
                    extrema.Points(i, p1, p2);
                    closest_p1 = p1.Value();
                    closest_p2 = p2.Value();
                }
            }
        }
    } catch (const Standard_Failure& e) {
        LOG_DEBUG("Extrema computation failed: ", e.GetMessageString());
    } catch (...) {
        LOG_DEBUG("Extrema computation failed with unknown exception");
    }
}

std::optional<SurfaceAdjustment> SurfaceIntersectionAnalyzer::compute_adjustment(
    const IntersectionResult& failed_result) const
{
    if (failed_result.success) {
        return std::nullopt;  // No adjustment needed
    }
    
    // Get surface types
    auto it1 = fitted_surfaces_.find(failed_result.surface_id_1);
    auto it2 = fitted_surfaces_.find(failed_result.surface_id_2);
    
    if (it1 == fitted_surfaces_.end() || it2 == fitted_surfaces_.end()) {
        return std::nullopt;
    }
    
    const FittedSurface& surf1 = it1->second;
    const FittedSurface& surf2 = it2->second;
    
    SurfaceAdjustment adj;
    adj.translation = Eigen::Vector3d::Zero();
    adj.rotation_axis = Eigen::Vector3d::UnitZ();
    adj.rotation_angle = 0;
    adj.size_delta = 0;
    adj.fitting_cost = 0;
    
    // Strategy: Prefer adjusting the more complex surface, or the smaller one
    // Planes are "simpler" than cylinders, which are simpler than cones/spheres
    
    auto surface_complexity = [](SurfaceType t) -> int {
        switch (t) {
            case SurfaceType::PLANE: return 1;
            case SurfaceType::CYLINDER: return 2;
            case SurfaceType::SPHERE: return 3;
            case SurfaceType::CONE: return 3;
            default: return 4;
        }
    };
    
    int c1 = surface_complexity(surf1.type);
    int c2 = surface_complexity(surf2.type);
    
    // Adjust the more complex surface, or surface 2 if equal
    const FittedSurface& to_adjust = (c1 > c2) ? surf1 : surf2;
    adj.surface_id = to_adjust.surface_id;
    adj.surface_type = to_adjust.type;
    
    // Compute adjustment direction from closest points
    gp_Vec gap_vec(failed_result.closest_point_1, failed_result.closest_point_2);
    double gap = gap_vec.Magnitude();
    
    if (gap < 1e-9) {
        // Surfaces are essentially touching
        adj.type = SurfaceAdjustment::AdjustmentType::NONE;
        return adj;
    }
    
    // For now, suggest a simple translation to close the gap
    // More sophisticated logic would analyze the geometry to decide between
    // translate, resize, rotate, or combined adjustments
    
    gp_Dir gap_dir = gap_vec;
    
    switch (to_adjust.type) {
        case SurfaceType::PLANE: {
            // For planes, translate along the gap direction
            adj.type = SurfaceAdjustment::AdjustmentType::TRANSLATE;
            // Move half the gap distance
            adj.translation = Eigen::Vector3d(
                gap_dir.X() * gap / 2,
                gap_dir.Y() * gap / 2,
                gap_dir.Z() * gap / 2
            );
            break;
        }
        
        case SurfaceType::CYLINDER: {
            // For cylinders, could be axis translation or radius adjustment
            // Check if gap direction is perpendicular to axis (suggests radius change)
            // or parallel to axis (suggests axis translation)
            gp_Dir axis_dir(to_adjust.coefficients[3], 
                           to_adjust.coefficients[4], 
                           to_adjust.coefficients[5]);
            
            double parallel_component = std::abs(gap_dir.Dot(axis_dir));
            
            if (parallel_component < 0.3) {
                // Gap is mostly perpendicular to axis -> radius adjustment likely
                // But first check if translation would work
                // TODO: More sophisticated analysis
                adj.type = SurfaceAdjustment::AdjustmentType::RESIZE;
                adj.size_delta = gap / 2;  // Increase radius by half the gap
                LOG_DEBUG("  Suggest cylinder radius adjustment: +", adj.size_delta);
            } else {
                // Gap has significant parallel component -> translate axis
                adj.type = SurfaceAdjustment::AdjustmentType::TRANSLATE;
                adj.translation = Eigen::Vector3d(
                    gap_dir.X() * gap / 2,
                    gap_dir.Y() * gap / 2,
                    gap_dir.Z() * gap / 2
                );
            }
            break;
        }
        
        case SurfaceType::SPHERE: {
            // For spheres, could be center translation or radius adjustment
            // If the other surface is a plane tangent to the sphere,
            // and the sphere center is correctly placed, need radius adjustment
            adj.type = SurfaceAdjustment::AdjustmentType::RESIZE;
            adj.size_delta = gap / 2;  // Adjust radius
            LOG_DEBUG("  Suggest sphere radius adjustment: +", adj.size_delta);
            break;
        }
        
        case SurfaceType::CONE: {
            // For cones, could be apex position, axis direction, or half-angle
            // Half-angle adjustment is like radius for cylinders
            adj.type = SurfaceAdjustment::AdjustmentType::RESIZE;
            // Convert gap to angle change is complex - simplified for now
            double current_half_angle = to_adjust.coefficients[6];
            // Rough approximation: small angle change
            adj.size_delta = std::atan(gap / 10.0);  // Very rough estimate
            LOG_DEBUG("  Suggest cone half-angle adjustment: +", adj.size_delta, " rad");
            break;
        }
        
        default:
            adj.type = SurfaceAdjustment::AdjustmentType::TRANSLATE;
            adj.translation = Eigen::Vector3d(
                gap_dir.X() * gap / 2,
                gap_dir.Y() * gap / 2,
                gap_dir.Z() * gap / 2
            );
            break;
    }
    
    return adj;
}

void SurfaceIntersectionAnalyzer::print_summary() const
{
    LOG_INFO("=== Surface Intersection Analysis Summary ===");
    LOG_INFO("Total adjacent pairs: ", adjacent_pairs_.size());
    
    int success = 0, failure = 0;
    double max_gap = 0;
    
    for (const auto& result : results_) {
        if (result.success) {
            ++success;
        } else {
            ++failure;
            max_gap = std::max(max_gap, result.min_distance);
        }
    }
    
    LOG_INFO("Intersections found: ", success);
    LOG_INFO("Intersections failed: ", failure);
    if (failure > 0) {
        LOG_INFO("Maximum gap: ", max_gap, " mm");
    }
    
    // Print details for failures
    for (const auto& result : results_) {
        if (!result.success && result.min_distance < std::numeric_limits<double>::max() / 2) {
            auto it1 = fitted_surfaces_.find(result.surface_id_1);
            auto it2 = fitted_surfaces_.find(result.surface_id_2);
            
            std::string type1 = "?", type2 = "?";
            if (it1 != fitted_surfaces_.end()) {
                switch (it1->second.type) {
                    case SurfaceType::PLANE: type1 = "PLANE"; break;
                    case SurfaceType::CYLINDER: type1 = "CYLINDER"; break;
                    case SurfaceType::SPHERE: type1 = "SPHERE"; break;
                    case SurfaceType::CONE: type1 = "CONE"; break;
                    default: type1 = "OTHER"; break;
                }
            }
            if (it2 != fitted_surfaces_.end()) {
                switch (it2->second.type) {
                    case SurfaceType::PLANE: type2 = "PLANE"; break;
                    case SurfaceType::CYLINDER: type2 = "CYLINDER"; break;
                    case SurfaceType::SPHERE: type2 = "SPHERE"; break;
                    case SurfaceType::CONE: type2 = "CONE"; break;
                    default: type2 = "OTHER"; break;
                }
            }
            
            LOG_INFO("  Gap between ", type1, " #", result.surface_id_1, 
                    " and ", type2, " #", result.surface_id_2, 
                    ": ", result.min_distance, " mm");
            
            // Compute and show suggested adjustment
            auto adj = compute_adjustment(result);
            if (adj.has_value()) {
                const auto& a = adj.value();
                switch (a.type) {
                    case SurfaceAdjustment::AdjustmentType::TRANSLATE:
                        LOG_INFO("    -> Suggest translating surface #", a.surface_id, 
                                " by (", a.translation.x(), ", ", a.translation.y(), 
                                ", ", a.translation.z(), ")");
                        break;
                    case SurfaceAdjustment::AdjustmentType::RESIZE:
                        LOG_INFO("    -> Suggest resizing surface #", a.surface_id, 
                                " by ", a.size_delta);
                        break;
                    default:
                        break;
                }
            }
        }
    }
    
    LOG_INFO("==============================================");
}

} // namespace brepper
