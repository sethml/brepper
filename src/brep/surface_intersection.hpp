#pragma once

#include "common/types.hpp"
#include "common/config.hpp"

#include <Geom_Surface.hxx>
#include <Geom_Curve.hxx>
#include <gp_Pnt.hxx>
#include <vector>
#include <map>
#include <set>
#include <optional>

namespace brepper {

// Result of attempting to intersect two surfaces
struct IntersectionResult {
    int surface_id_1;
    int surface_id_2;
    bool success;                           // Did GeomAPI_IntSS find an intersection?
    std::vector<Handle(Geom_Curve)> curves; // Intersection curves if successful
    
    // If intersection failed, diagnostic info:
    double min_distance;                    // Minimum distance between surfaces
    gp_Pnt closest_point_1;                 // Closest point on surface 1
    gp_Pnt closest_point_2;                 // Closest point on surface 2
    
    // Expected intersection location (from mesh boundary)
    std::vector<gp_Pnt> expected_points;
};

// Suggested adjustment to make surfaces intersect
struct SurfaceAdjustment {
    int surface_id;
    SurfaceType surface_type;
    
    // What kind of adjustment is needed
    enum class AdjustmentType {
        NONE,
        TRANSLATE,      // Shift position (axis for cylinder/cone, center for sphere)
        ROTATE,         // Rotate axis direction (cylinder/cone)
        RESIZE,         // Adjust radius (cylinder/sphere) or half-angle (cone)
        COMBINED        // Multiple adjustments needed
    };
    AdjustmentType type;
    
    // Adjustment parameters (interpretation depends on type and surface_type)
    Eigen::Vector3d translation;      // For TRANSLATE
    Eigen::Vector3d rotation_axis;    // For ROTATE
    double rotation_angle;            // For ROTATE (radians)
    double size_delta;                // For RESIZE (radius change or angle change)
    
    // How much does this adjustment deviate from the original fit?
    double fitting_cost;              // RMS change in surface position at sample points
};

// Analyzes surface-surface intersections and computes necessary adjustments
class SurfaceIntersectionAnalyzer {
public:
    explicit SurfaceIntersectionAnalyzer(const Config& config);
    
    // Build adjacency graph from boundary curves
    void build_adjacency_graph(
        const std::vector<FittedSurface>& surfaces,
        const std::vector<BoundaryCurve>& boundaries);
    
    // Create OCCT surfaces from fitted surfaces
    void create_surfaces(const std::vector<FittedSurface>& surfaces);
    
    // Attempt intersection for all adjacent surface pairs
    void analyze_intersections();
    
    // Get results
    const std::vector<IntersectionResult>& get_results() const { return results_; }
    const std::map<std::pair<int,int>, IntersectionResult>& get_result_map() const { return result_map_; }
    
    // Get adjacency information
    const std::set<std::pair<int,int>>& get_adjacent_pairs() const { return adjacent_pairs_; }
    
    // Compute suggested adjustment for a failed intersection
    std::optional<SurfaceAdjustment> compute_adjustment(
        const IntersectionResult& failed_result) const;
    
    // Print diagnostic summary
    void print_summary() const;
    
private:
    const Config& config_;
    
    // Surface data
    std::map<int, FittedSurface> fitted_surfaces_;
    std::map<int, Handle(Geom_Surface)> geom_surfaces_;
    
    // Adjacency graph: pairs of surface IDs that share a boundary
    std::set<std::pair<int,int>> adjacent_pairs_;
    
    // Expected boundary points between surface pairs (from mesh)
    std::map<std::pair<int,int>, std::vector<gp_Pnt>> expected_boundaries_;
    
    // Results
    std::vector<IntersectionResult> results_;
    std::map<std::pair<int,int>, IntersectionResult> result_map_;
    
    // Helper methods
    Handle(Geom_Surface) create_surface(const FittedSurface& surface);
    Handle(Geom_Surface) create_plane(const FittedSurface& surface);
    Handle(Geom_Surface) create_cylinder(const FittedSurface& surface);
    Handle(Geom_Surface) create_sphere(const FittedSurface& surface);
    Handle(Geom_Surface) create_cone(const FittedSurface& surface);
    
    // Attempt intersection between two surfaces
    IntersectionResult try_intersection(int id1, int id2);
    
    // Compute minimum distance when intersection fails
    void compute_extrema(
        const Handle(Geom_Surface)& s1,
        const Handle(Geom_Surface)& s2,
        const std::vector<gp_Pnt>& hint_points,
        double& min_dist,
        gp_Pnt& closest_p1,
        gp_Pnt& closest_p2);
};

} // namespace brepper
