#pragma once

#include "common/types.hpp"
#include "common/config.hpp"

#include <Geom_Curve.hxx>
#include <Geom_Surface.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Wire.hxx>
#include <TopoDS_Face.hxx>
#include <gp_Pnt.hxx>
#include <vector>
#include <map>
#include <set>
#include <optional>

namespace brepper {

/// Result of intersecting two surfaces
struct SurfaceIntersectionResult {
    int surface_id_a;
    int surface_id_b;
    std::vector<Handle(Geom_Curve)> curves;     // Intersection curves (may be multiple)
    std::vector<TopoDS_Edge> edges;              // Trimmed edges
    bool success = false;
    std::string error_message;
};

/// An edge shared between two faces, with proper pcurves
struct SharedEdge {
    TopoDS_Edge edge;
    int surface_id_a;
    int surface_id_b;
    double param_start;
    double param_end;
};

/// Computes surface-surface intersections and creates shared edges
class SurfaceIntersector {
public:
    explicit SurfaceIntersector(const Config& config);
    
    /// Compute intersection between two surfaces
    /// Returns curves trimmed to the region defined by boundary points
    SurfaceIntersectionResult intersect_surfaces(
        int surface_id_a,
        int surface_id_b,
        const Handle(Geom_Surface)& surf_a,
        const Handle(Geom_Surface)& surf_b,
        const std::vector<Eigen::Vector3d>& boundary_points
    );
    
    /// Process all boundaries and create shared edges
    /// Returns a map from (surface_id_a, surface_id_b) -> edges
    std::map<std::pair<int,int>, std::vector<SharedEdge>> create_all_shared_edges(
        const std::vector<BoundaryCurve>& boundaries,
        const std::vector<FittedSurface>& surfaces,
        const std::map<int, Handle(Geom_Surface)>& geom_surfaces
    );
    
    /// Get all edges for a specific surface (for building its wire)
    std::vector<TopoDS_Edge> get_edges_for_surface(
        int surface_id,
        const std::map<std::pair<int,int>, std::vector<SharedEdge>>& all_edges
    );
    
    /// Build a wire from edges for a surface
    std::optional<TopoDS_Wire> build_wire_for_surface(
        int surface_id,
        const std::map<std::pair<int,int>, std::vector<SharedEdge>>& all_edges
    );
    
    /// Build multiple wires for a surface (for multi-boundary faces like cylinders)
    std::vector<TopoDS_Wire> build_wires_for_surface(
        int surface_id,
        const std::map<std::pair<int,int>, std::vector<SharedEdge>>& all_edges
    );
    
    /// Create a trimmed face using the surface and its boundary wire
    std::optional<TopoDS_Face> create_trimmed_face(
        const Handle(Geom_Surface)& surface,
        const TopoDS_Wire& boundary_wire
    );
    
    /// Create a face with multiple boundary wires
    std::optional<TopoDS_Face> create_face_with_wires(
        const Handle(Geom_Surface)& surface,
        const std::vector<TopoDS_Wire>& wires
    );
    
private:
    const Config& config_;
    double tolerance_ = 1e-4;  // Intersection tolerance
    
    /// Find the curve parameter range that covers the boundary points
    std::pair<double, double> compute_trim_parameters(
        const Handle(Geom_Curve)& curve,
        const std::vector<Eigen::Vector3d>& boundary_points
    );
    
    /// Create an edge with pcurves on both surfaces
    TopoDS_Edge create_shared_edge(
        const Handle(Geom_Curve)& curve,
        const Handle(Geom_Surface)& surf_a,
        const Handle(Geom_Surface)& surf_b,
        double param_start,
        double param_end
    );
    
    /// Select the intersection curve closest to the boundary points
    /// (handles cases where intersection returns multiple branches)
    int select_best_curve(
        const std::vector<Handle(Geom_Curve)>& curves,
        const std::vector<Eigen::Vector3d>& boundary_points
    );
    
    /// Check if a curve segment is near the boundary points
    bool curve_near_boundary(
        const Handle(Geom_Curve)& curve,
        double t_start,
        double t_end,
        const std::vector<Eigen::Vector3d>& boundary_points,
        double max_distance
    );
};

} // namespace brepper
