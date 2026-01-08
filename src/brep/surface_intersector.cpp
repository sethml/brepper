#include "surface_intersector.hpp"
#include "common/logging.hpp"

#include <GeomAPI_IntSS.hxx>
#include <GeomAPI_ProjectPointOnCurve.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <ShapeFix_Wire.hxx>
#include <ShapeFix_Edge.hxx>
#include <ShapeFix_Face.hxx>
#include <ShapeExtend_Status.hxx>
#include <BRep_Tool.hxx>
#include <BRep_Builder.hxx>
#include <TopoDS.hxx>
#include <Geom_Line.hxx>
#include <Geom_Circle.hxx>
#include <Geom_TrimmedCurve.hxx>
#include <Geom2d_Curve.hxx>
#include <Geom2d_Line.hxx>
#include <Geom2d_OffsetCurve.hxx>
#include <ShapeConstruct_ProjectCurveOnSurface.hxx>
#include <ShapeAnalysis_Surface.hxx>
#include <gp_Pnt.hxx>
#include <gp_Pnt2d.hxx>
#include <gp_Vec2d.hxx>
#include <gp_Lin.hxx>
#include <gp_Circ.hxx>
#include <Precision.hxx>
#include <Standard_Failure.hxx>
#include <TopLoc_Location.hxx>

#include <algorithm>
#include <cmath>
#include <limits>

namespace brepper {

// Helper to shift pcurve U values into the canonical range [0, 2π) for U-periodic surfaces
// We normalize based on the START point to ensure consistency across edges
static Handle(Geom2d_Curve) normalize_pcurve_for_periodic_surface(
    const Handle(Geom2d_Curve)& pcurve,
    const Handle(Geom_Surface)& surface,
    double param_start,
    double /* param_end */
) {
    if (pcurve.IsNull() || surface.IsNull()) {
        return pcurve;
    }
    
    // Check if surface is U-periodic (cylinders, cones, spheres, etc.)
    if (!surface->IsUPeriodic()) {
        return pcurve;
    }
    
    double u_period = surface->UPeriod();
    
    // Use the start point's U value to determine the shift
    // This ensures all pcurves on the same surface are shifted consistently
    gp_Pnt2d p_start = pcurve->Value(param_start);
    double u_start = p_start.X();
    
    // Calculate shift needed to put u_start in [0, u_period)
    double shift = -u_period * std::floor(u_start / u_period);
    
    // Only apply shift if it's significant
    if (std::abs(shift) < 1e-10) {
        return pcurve;
    }
    
    // Create a translated copy of the pcurve
    Handle(Geom2d_Geometry) geom = pcurve->Translated(gp_Vec2d(shift, 0.0));
    Handle(Geom2d_Curve) shifted = Handle(Geom2d_Curve)::DownCast(geom);
    
    return shifted.IsNull() ? pcurve : shifted;
}

SurfaceIntersector::SurfaceIntersector(const Config& config) 
    : config_(config) 
{
    tolerance_ = 1e-4;  // 0.1 micron default
}

SurfaceIntersectionResult SurfaceIntersector::intersect_surfaces(
    int surface_id_a,
    int surface_id_b,
    const Handle(Geom_Surface)& surf_a,
    const Handle(Geom_Surface)& surf_b,
    const std::vector<Eigen::Vector3d>& boundary_points
) {
    SurfaceIntersectionResult result;
    result.surface_id_a = surface_id_a;
    result.surface_id_b = surface_id_b;
    
    if (surf_a.IsNull() || surf_b.IsNull()) {
        result.error_message = "Null surface(s)";
        return result;
    }
    
    if (boundary_points.empty()) {
        result.error_message = "No boundary points";
        return result;
    }
    
    try {
        // Compute surface-surface intersection
        GeomAPI_IntSS intersector(surf_a, surf_b, tolerance_);
        
        if (!intersector.IsDone()) {
            result.error_message = "Intersection computation failed";
            return result;
        }
        
        int num_curves = intersector.NbLines();
        if (num_curves == 0) {
            result.error_message = "No intersection curves found";
            return result;
        }
        
        LOG_DEBUG("Surface intersection found ", num_curves, " curve(s) between surfaces #",
                  surface_id_a, " and #", surface_id_b);
        
        // Collect all intersection curves
        for (int i = 1; i <= num_curves; ++i) {
            Handle(Geom_Curve) curve = intersector.Line(i);
            if (!curve.IsNull()) {
                result.curves.push_back(curve);
            }
        }
        
        if (result.curves.empty()) {
            result.error_message = "All intersection curves were null";
            return result;
        }
        
        // Select the best curve (closest to boundary points)
        int best_idx = select_best_curve(result.curves, boundary_points);
        if (best_idx < 0) {
            result.error_message = "No intersection curve near boundary points";
            return result;
        }
        
        Handle(Geom_Curve) best_curve = result.curves[best_idx];
        
        // Compute trim parameters from boundary points
        auto [t_start, t_end] = compute_trim_parameters(best_curve, boundary_points);
        
        if (std::abs(t_end - t_start) < Precision::Confusion()) {
            result.error_message = "Degenerate trim range";
            return result;
        }
        
        LOG_DEBUG("Trim parameters: t=[", t_start, ", ", t_end, "]");
        
        // Create edge with pcurves on both surfaces
        TopoDS_Edge edge = create_shared_edge(best_curve, surf_a, surf_b, t_start, t_end);
        
        if (edge.IsNull()) {
            result.error_message = "Failed to create shared edge";
            return result;
        }
        
        result.edges.push_back(edge);
        result.success = true;
        
        // Log curve type for debugging
        if (best_curve->IsKind(STANDARD_TYPE(Geom_Line))) {
            LOG_DEBUG("Created LINE edge between surfaces #", surface_id_a, " and #", surface_id_b);
        } else if (best_curve->IsKind(STANDARD_TYPE(Geom_Circle))) {
            LOG_DEBUG("Created CIRCLE edge between surfaces #", surface_id_a, " and #", surface_id_b);
        } else {
            LOG_DEBUG("Created curve edge between surfaces #", surface_id_a, " and #", surface_id_b);
        }
        
    } catch (const Standard_Failure& e) {
        result.error_message = std::string("OCCT exception: ") + e.GetMessageString();
        LOG_ERROR("Surface intersection failed: ", result.error_message);
    } catch (const std::exception& e) {
        result.error_message = std::string("Exception: ") + e.what();
        LOG_ERROR("Surface intersection failed: ", result.error_message);
    }
    
    return result;
}

std::pair<double, double> SurfaceIntersector::compute_trim_parameters(
    const Handle(Geom_Curve)& curve,
    const std::vector<Eigen::Vector3d>& boundary_points
) {
    double t_min = std::numeric_limits<double>::max();
    double t_max = std::numeric_limits<double>::lowest();
    
    for (const auto& pt : boundary_points) {
        gp_Pnt point(pt.x(), pt.y(), pt.z());
        
        try {
            GeomAPI_ProjectPointOnCurve projector(point, curve);
            
            if (projector.NbPoints() > 0) {
                double t = projector.LowerDistanceParameter();
                t_min = std::min(t_min, t);
                t_max = std::max(t_max, t);
            }
        } catch (...) {
            // Skip points that can't be projected
            continue;
        }
    }
    
    // Handle infinite curves (lines) - extend slightly beyond points
    if (curve->IsKind(STANDARD_TYPE(Geom_Line))) {
        // For lines, ensure we have valid bounds
        if (t_min > t_max) {
            t_min = 0;
            t_max = 1;
        }
    }
    
    // Handle periodic curves (full circles)
    if (curve->IsPeriodic()) {
        double period = curve->Period();
        // If the parameter range spans nearly the full period, use full period
        if ((t_max - t_min) > 0.9 * period) {
            t_min = curve->FirstParameter();
            t_max = curve->LastParameter();
        }
    }
    
    return {t_min, t_max};
}

TopoDS_Edge SurfaceIntersector::create_shared_edge(
    const Handle(Geom_Curve)& curve,
    const Handle(Geom_Surface)& surf_a,
    const Handle(Geom_Surface)& surf_b,
    double param_start,
    double param_end
) {
    try {
        TopoDS_Edge edge;
        
        // Get the underlying curve if this is a trimmed curve
        Handle(Geom_Curve) basis_curve = curve;
        Handle(Geom_TrimmedCurve) trimmed = Handle(Geom_TrimmedCurve)::DownCast(curve);
        if (!trimmed.IsNull()) {
            basis_curve = trimmed->BasisCurve();
        }
        
        // For periodic curves (circles, ellipses) spanning nearly the full period,
        // create a closed edge without explicit parameter bounds
        if (basis_curve->IsPeriodic()) {
            double period = basis_curve->Period();
            double range = param_end - param_start;
            
            // If we span more than 90% of the period, use the full curve
            if (range > 0.9 * period) {
                
                // For a closed edge, use the basis curve with its natural bounds
                BRepBuilderAPI_MakeEdge edgeBuilder(basis_curve);
                
                if (!edgeBuilder.IsDone()) {
                    LOG_ERROR("Closed edge creation from periodic curve failed");
                    return TopoDS_Edge();
                }
                
                edge = edgeBuilder.Edge();
                
                // Update param_start/param_end to match the basis curve's natural bounds
                // This is important for pcurve computation below
                param_start = basis_curve->FirstParameter();
                param_end = basis_curve->LastParameter();
            }
        }
        
        // For non-periodic curves or partial periodic curves, use explicit bounds
        if (edge.IsNull()) {
            BRepBuilderAPI_MakeEdge edgeBuilder(curve, param_start, param_end);
            
            if (!edgeBuilder.IsDone()) {
                LOG_ERROR("Edge creation from curve failed");
                return TopoDS_Edge();
            }
            
            edge = edgeBuilder.Edge();
        }
        
        // Now add pcurves (2D parametric curves) for both surfaces
        // This is critical for proper face construction and sewing
        BRep_Builder builder;
        TopLoc_Location identity;  // No transformation
        
        // For pcurve projection, use the basis curve if available
        Handle(Geom_Curve) curve_for_pcurve = basis_curve.IsNull() ? curve : basis_curve;
        
        // Project 3D curve onto surface A to get pcurve A
        Handle(Geom2d_Curve) pcurve_a;
        {
            ShapeConstruct_ProjectCurveOnSurface projector;
            projector.Init(surf_a, tolerance_);
            
            Handle(Geom2d_Curve) c2d;
            if (projector.Perform(curve_for_pcurve, param_start, param_end, c2d)) {
                // Normalize pcurve to canonical parameter range for periodic surfaces
                pcurve_a = normalize_pcurve_for_periodic_surface(c2d, surf_a, param_start, param_end);
                LOG_DEBUG("Computed pcurve A for edge");
            } else {
                LOG_WARN("Failed to project curve onto surface A");
            }
        }
        
        // Project 3D curve onto surface B to get pcurve B
        Handle(Geom2d_Curve) pcurve_b;
        {
            ShapeConstruct_ProjectCurveOnSurface projector;
            projector.Init(surf_b, tolerance_);
            
            Handle(Geom2d_Curve) c2d;
            if (projector.Perform(curve_for_pcurve, param_start, param_end, c2d)) {
                // Normalize pcurve to canonical parameter range for periodic surfaces
                pcurve_b = normalize_pcurve_for_periodic_surface(c2d, surf_b, param_start, param_end);
                LOG_DEBUG("Computed pcurve B for edge");
            } else {
                LOG_WARN("Failed to project curve onto surface B");
            }
        }
        
        // Add pcurves to the edge using BRep_Builder::UpdateEdge
        if (!pcurve_a.IsNull()) {
            builder.UpdateEdge(edge, pcurve_a, surf_a, identity, tolerance_);
        }
        if (!pcurve_b.IsNull()) {
            builder.UpdateEdge(edge, pcurve_b, surf_b, identity, tolerance_);
        }
        
        return edge;
        
    } catch (const Standard_Failure& e) {
        LOG_ERROR("Exception creating shared edge: ", e.GetMessageString());
        return TopoDS_Edge();
    }
}

int SurfaceIntersector::select_best_curve(
    const std::vector<Handle(Geom_Curve)>& curves,
    const std::vector<Eigen::Vector3d>& boundary_points
) {
    if (curves.empty()) return -1;
    if (curves.size() == 1) return 0;
    
    // Find the curve closest to the boundary points
    double best_distance = std::numeric_limits<double>::max();
    int best_idx = -1;
    
    for (size_t i = 0; i < curves.size(); ++i) {
        const Handle(Geom_Curve)& curve = curves[i];
        if (curve.IsNull()) continue;
        
        double total_dist = 0;
        int count = 0;
        
        for (const auto& pt : boundary_points) {
            gp_Pnt point(pt.x(), pt.y(), pt.z());
            try {
                GeomAPI_ProjectPointOnCurve projector(point, curve);
                if (projector.NbPoints() > 0) {
                    total_dist += projector.LowerDistance();
                    count++;
                }
            } catch (...) {
                continue;
            }
        }
        
        if (count > 0) {
            double avg_dist = total_dist / count;
            if (avg_dist < best_distance) {
                best_distance = avg_dist;
                best_idx = static_cast<int>(i);
            }
        }
    }
    
    if (best_idx >= 0) {
        LOG_DEBUG("Selected curve ", best_idx, " with average distance ", best_distance);
    }
    
    return best_idx;
}

bool SurfaceIntersector::curve_near_boundary(
    const Handle(Geom_Curve)& curve,
    double t_start,
    double t_end,
    const std::vector<Eigen::Vector3d>& boundary_points,
    double max_distance
) {
    // Sample curve and check distances
    int num_samples = 10;
    double dt = (t_end - t_start) / (num_samples - 1);
    
    for (int i = 0; i < num_samples; ++i) {
        double t = t_start + i * dt;
        gp_Pnt curve_pt = curve->Value(t);
        
        // Find minimum distance to any boundary point
        double min_dist = std::numeric_limits<double>::max();
        for (const auto& bp : boundary_points) {
            double dx = curve_pt.X() - bp.x();
            double dy = curve_pt.Y() - bp.y();
            double dz = curve_pt.Z() - bp.z();
            double dist = std::sqrt(dx*dx + dy*dy + dz*dz);
            min_dist = std::min(min_dist, dist);
        }
        
        if (min_dist > max_distance) {
            return false;
        }
    }
    
    return true;
}

std::map<std::pair<int,int>, std::vector<SharedEdge>> SurfaceIntersector::create_all_shared_edges(
    const std::vector<BoundaryCurve>& boundaries,
    const std::vector<FittedSurface>& surfaces,
    const std::map<int, Handle(Geom_Surface)>& geom_surfaces
) {
    std::map<std::pair<int,int>, std::vector<SharedEdge>> result;
    
    // Create a lookup for surfaces by ID
    std::map<int, const FittedSurface*> surface_lookup;
    for (const auto& surf : surfaces) {
        surface_lookup[surf.surface_id] = &surf;
    }
    
    for (size_t i = 0; i < boundaries.size(); ++i) {
        const BoundaryCurve& boundary = boundaries[i];
        
        int id_a = boundary.surface_id_left;
        int id_b = boundary.surface_id_right;
        
        // Normalize the pair (smaller ID first)
        if (id_a > id_b) std::swap(id_a, id_b);
        auto key = std::make_pair(id_a, id_b);
        
        // Skip if we don't have both surfaces
        auto it_a = geom_surfaces.find(id_a);
        auto it_b = geom_surfaces.find(id_b);
        if (it_a == geom_surfaces.end() || it_b == geom_surfaces.end()) {
            LOG_WARN("Missing surface for boundary between ", id_a, " and ", id_b);
            continue;
        }
        
        // Intersect the surfaces
        auto intersection = intersect_surfaces(
            id_a, id_b,
            it_a->second, it_b->second,
            boundary.points
        );
        
        if (intersection.success) {
            for (const auto& edge : intersection.edges) {
                SharedEdge shared;
                shared.edge = edge;
                shared.surface_id_a = id_a;
                shared.surface_id_b = id_b;
                result[key].push_back(shared);
            }
        } else {
            LOG_WARN("Failed to create edge for boundary ", i, " between surfaces ", 
                     id_a, " and ", id_b, ": ", intersection.error_message);
        }
    }
    
    LOG_INFO("Created shared edges for ", result.size(), " surface pairs from ", 
             boundaries.size(), " boundaries");
    
    return result;
}

std::vector<TopoDS_Edge> SurfaceIntersector::get_edges_for_surface(
    int surface_id,
    const std::map<std::pair<int,int>, std::vector<SharedEdge>>& all_edges
) {
    std::vector<TopoDS_Edge> edges;
    
    for (const auto& [key, shared_edges] : all_edges) {
        if (key.first == surface_id || key.second == surface_id) {
            for (const auto& se : shared_edges) {
                edges.push_back(se.edge);
            }
        }
    }
    
    return edges;
}

std::optional<TopoDS_Wire> SurfaceIntersector::build_wire_for_surface(
    int surface_id,
    const std::map<std::pair<int,int>, std::vector<SharedEdge>>& all_edges
) {
    auto edges = get_edges_for_surface(surface_id, all_edges);
    
    if (edges.empty()) {
        LOG_DEBUG("No edges found for surface #", surface_id);
        return std::nullopt;
    }
    
    // If we have a single closed edge (like a circle), it forms its own wire
    if (edges.size() == 1) {
        try {
            BRepBuilderAPI_MakeWire wireBuilder;
            wireBuilder.Add(edges[0]);
            
            if (wireBuilder.IsDone()) {
                TopoDS_Wire wire = wireBuilder.Wire();
                if (!wire.Closed()) {
                    LOG_DEBUG("Wire not closed for surface #", surface_id, ", has 1 edge");
                }
                return wire;
            }
        } catch (const Standard_Failure& e) {
            LOG_ERROR("Exception building single-edge wire: ", e.GetMessageString());
        }
        return std::nullopt;
    }
    
    // For multiple edges, try to build connected wires
    // First, try adding all edges to one wire
    try {
        BRepBuilderAPI_MakeWire wireBuilder;
        
        for (const auto& edge : edges) {
            wireBuilder.Add(edge);
            if (wireBuilder.Error() != BRepBuilderAPI_WireDone) {
                // Edges don't connect - this is expected for cylinders with two circular edges
                LOG_DEBUG("Edges don't form connected wire for surface #", surface_id,
                          " (", edges.size(), " edges) - this may be multi-boundary surface");
                return std::nullopt;
            }
        }
        
        if (!wireBuilder.IsDone()) {
            LOG_DEBUG("Wire construction not done for surface #", surface_id);
            return std::nullopt;
        }
        
        TopoDS_Wire wire = wireBuilder.Wire();
        
        // Check if wire is closed
        if (!wire.Closed()) {
            LOG_DEBUG("Wire not closed for surface #", surface_id, ", has ", edges.size(), " edges");
        }
        
        return wire;
        
    } catch (const Standard_Failure& e) {
        LOG_ERROR("Exception building wire: ", e.GetMessageString());
        return std::nullopt;
    }
}

std::optional<TopoDS_Face> SurfaceIntersector::create_trimmed_face(
    const Handle(Geom_Surface)& surface,
    const TopoDS_Wire& boundary_wire
) {
    if (surface.IsNull()) {
        LOG_ERROR("Null surface for trimmed face");
        return std::nullopt;
    }
    
    try {
        BRepBuilderAPI_MakeFace faceBuilder(surface, boundary_wire, Standard_True);
        
        if (!faceBuilder.IsDone()) {
            LOG_ERROR("Face creation failed with error: ", faceBuilder.Error());
            return std::nullopt;
        }
        
        TopoDS_Face face = faceBuilder.Face();
        return face;
        
    } catch (const Standard_Failure& e) {
        LOG_ERROR("Exception creating trimmed face: ", e.GetMessageString());
        return std::nullopt;
    }
}

std::vector<TopoDS_Wire> SurfaceIntersector::build_wires_for_surface(
    int surface_id,
    const std::map<std::pair<int,int>, std::vector<SharedEdge>>& all_edges
) {
    std::vector<TopoDS_Wire> wires;
    auto edges = get_edges_for_surface(surface_id, all_edges);
    
    if (edges.empty()) {
        LOG_DEBUG("No edges found for surface #", surface_id);
        return wires;
    }
    
    // Each edge becomes its own wire (for now - assumes edges are closed curves like circles)
    for (const auto& edge : edges) {
        try {
            BRepBuilderAPI_MakeWire wireBuilder;
            wireBuilder.Add(edge);
            
            if (wireBuilder.IsDone()) {
                TopoDS_Wire wire = wireBuilder.Wire();
                wires.push_back(wire);
                LOG_DEBUG("Created wire from single edge for surface #", surface_id);
            }
        } catch (const Standard_Failure& e) {
            LOG_ERROR("Exception building wire from edge: ", e.GetMessageString());
        }
    }
    
    LOG_DEBUG("Built ", wires.size(), " wire(s) for surface #", surface_id);
    return wires;
}

std::optional<TopoDS_Face> SurfaceIntersector::create_face_with_wires(
    const Handle(Geom_Surface)& surface,
    const std::vector<TopoDS_Wire>& wires
) {
    if (surface.IsNull()) {
        LOG_ERROR("Null surface for face creation");
        return std::nullopt;
    }
    
    if (wires.empty()) {
        LOG_ERROR("No wires for face creation");
        return std::nullopt;
    }
    
    try {
        // Create face with first wire
        BRepBuilderAPI_MakeFace faceBuilder(surface, wires[0], Standard_True);
        
        if (!faceBuilder.IsDone()) {
            LOG_ERROR("Initial face creation failed with error: ", faceBuilder.Error());
            return std::nullopt;
        }
        
        // Add additional wires (inner boundaries)
        for (size_t i = 1; i < wires.size(); ++i) {
            faceBuilder.Add(wires[i]);
            if (!faceBuilder.IsDone()) {
                LOG_ERROR("Failed to add wire ", i, " to face, error: ", faceBuilder.Error());
            }
        }
        
        TopoDS_Face face = faceBuilder.Face();
        
        // For U-periodic surfaces, use ShapeFix_Face to fix the pcurve parameterization
        // This should handle cases where different edges have pcurves in different periods
        if (surface->IsUPeriodic()) {
            ShapeFix_Face fixer(face);
            fixer.SetPrecision(tolerance_);
            
            // Fix all issues
            fixer.FixOrientation();
            fixer.FixAddNaturalBound();
            fixer.FixPeriodicDegenerated();
            fixer.Perform();
            
            if (fixer.Status(ShapeExtend_DONE)) {
                face = fixer.Face();
                LOG_DEBUG("Applied ShapeFix_Face to periodic surface face");
            }
        }
        
        return face;
        
    } catch (const Standard_Failure& e) {
        LOG_ERROR("Exception creating face with wires: ", e.GetMessageString());
        return std::nullopt;
    }
}

} // namespace brepper
