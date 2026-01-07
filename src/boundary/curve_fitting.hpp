#pragma once

#include "common/types.hpp"
#include "common/config.hpp"

#include <Geom_Curve.hxx>
#include <Geom_Surface.hxx>
#include <TopoDS_Edge.hxx>
#include <gp_Pnt.hxx>
#include <vector>
#include <optional>
#include <map>

namespace brepper {

/// Result of fitting a curve to boundary points
struct CurveFitResult {
    Handle(Geom_Curve) curve;        // The fitted curve (may be null if fitting failed)
    TopoDS_Edge edge;                // The edge created from the curve
    double fitting_error;            // RMS error of fit
    std::string curve_type;          // "line", "circle", "ellipse", "bspline"
    bool success;
};

/// Types of analytic curves
enum class CurveType {
    LINE,
    CIRCLE, 
    ELLIPSE,
    BSPLINE
};

/// Fits analytic curves to boundary chains based on adjacent surface types
class BoundaryCurveFitter {
public:
    explicit BoundaryCurveFitter(const Config& config);
    
    /// Fit a curve to a boundary between two surfaces
    /// Uses surface type information to constrain the fit
    CurveFitResult fit_boundary_curve(
        const BoundaryCurve& boundary,
        const FittedSurface& surface_left,
        const FittedSurface& surface_right,
        const Handle(Geom_Surface)& geom_left,
        const Handle(Geom_Surface)& geom_right
    );
    
    /// Fit all boundaries and return edges keyed by boundary index
    std::map<int, CurveFitResult> fit_all_boundaries(
        const std::vector<BoundaryCurve>& boundaries,
        const std::vector<FittedSurface>& surfaces,
        const std::map<int, Handle(Geom_Surface)>& geom_surfaces
    );
    
private:
    const Config& config_;
    
    /// Determine expected curve type from two surface types
    std::vector<CurveType> get_candidate_curve_types(
        SurfaceType type_a, 
        SurfaceType type_b
    );
    
    /// Fit a line to points
    CurveFitResult fit_line(const std::vector<Eigen::Vector3d>& points);
    
    /// Fit a circle to points (3D)
    CurveFitResult fit_circle(const std::vector<Eigen::Vector3d>& points);
    
    /// Fit a circle constrained to a plane
    CurveFitResult fit_circle_on_plane(
        const std::vector<Eigen::Vector3d>& points,
        const FittedSurface& plane_surface
    );
    
    /// Fit a circle at plane-cylinder intersection
    CurveFitResult fit_plane_cylinder_intersection(
        const std::vector<Eigen::Vector3d>& points,
        const FittedSurface& plane,
        const FittedSurface& cylinder
    );
    
    /// Fit a circle at plane-sphere intersection
    CurveFitResult fit_plane_sphere_intersection(
        const std::vector<Eigen::Vector3d>& points,
        const FittedSurface& plane,
        const FittedSurface& sphere
    );
    
    /// Fit a circle at plane-cone intersection
    CurveFitResult fit_plane_cone_intersection(
        const std::vector<Eigen::Vector3d>& points,
        const FittedSurface& plane,
        const FittedSurface& cone
    );
    
    /// Fit an ellipse to points
    CurveFitResult fit_ellipse(const std::vector<Eigen::Vector3d>& points);
    
    /// Fit a B-spline curve (fallback)
    CurveFitResult fit_bspline(const std::vector<Eigen::Vector3d>& points);
    
    /// Compute RMS error of a curve fit
    double compute_fit_error(
        const Handle(Geom_Curve)& curve,
        const std::vector<Eigen::Vector3d>& points
    );
    
    /// Create an edge from a curve, trimmed to the point range
    TopoDS_Edge create_edge_from_curve(
        const Handle(Geom_Curve)& curve,
        const std::vector<Eigen::Vector3d>& points
    );
};

} // namespace brepper
