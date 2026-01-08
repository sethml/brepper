#include "brep_builder.hpp"
#include "surface_intersector.hpp"
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
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepBuilderAPI_Sewing.hxx>
#include <BRepCheck_Shell.hxx>
#include <ShapeFix_Shape.hxx>
#include <ShapeFix_Shell.hxx>
#include <ShapeFix_Solid.hxx>
#include <ShapeFix_Wire.hxx>
#include <ShapeFix_Face.hxx>
#include <ShapeFix_Edge.hxx>
#include <ShapeExtend_Status.hxx>
#include <ShapeExtend_WireData.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Shell.hxx>
#include <TopoDS_Solid.hxx>
#include <BRep_Builder.hxx>
#include <GeomAPI_ProjectPointOnSurf.hxx>
#include <Standard_Failure.hxx>
#include <cmath>
#include <limits>
#include <set>

namespace brepper {

BRepBuilder::BRepBuilder(const Config& config) : config_(config) {}

bool BRepBuilder::build(
    const std::vector<FittedSurface>& surfaces,
    const std::vector<BoundaryCurve>& boundaries,
    TopoDS_Shape& result
) {
    if (surfaces.empty()) {
        LOG_ERROR("No surfaces provided for B-Rep construction");
        return false;
    }
    
    LOG_INFO("Building B-Rep from ", surfaces.size(), " surfaces and ", 
             boundaries.size(), " boundary curves");
    
    // Store created OCCT surfaces for reuse
    std::map<int, Handle(Geom_Surface)> surface_map;
    for (const auto& surface : surfaces) {
        Handle(Geom_Surface) geom_surface = create_surface(surface);
        if (!geom_surface.IsNull()) {
            surface_map[surface.surface_id] = geom_surface;
        }
    }
    
    std::vector<TopoDS_Face> faces;
    faces.reserve(surfaces.size());
    
    for (const auto& surface : surfaces) {
        auto it = surface_map.find(surface.surface_id);
        if (it == surface_map.end()) {
            LOG_WARN("Failed to create OCCT surface for surface #", surface.surface_id, 
                     " (type: ", static_cast<int>(surface.type), ")");
            continue;
        }
        
        const Handle(Geom_Surface)& geom_surface = it->second;
        
        // Use bbox-bounded faces for all surfaces
        // The sewing step will handle trimming where surfaces intersect
        TopoDS_Face face = create_face_with_bounds(geom_surface, surface);
        if (!face.IsNull()) {
            LOG_DEBUG("Created bbox-bounded face for surface #", surface.surface_id);
        }
        
        if (face.IsNull()) {
            LOG_WARN("Failed to create face for surface #", surface.surface_id);
            continue;
        }
        
        faces.push_back(face);
    }
    
    if (faces.empty()) {
        LOG_ERROR("No valid faces created");
        return false;
    }
    
    LOG_INFO("Created ", faces.size(), " faces from ", surfaces.size(), " surfaces");
    
    if (faces.size() == 1) {
        result = faces[0];
    } else {
        // Try to sew faces together
        if (!sew_faces(faces, result)) {
            LOG_WARN("Face sewing failed, creating compound of faces");
            BRep_Builder builder;
            TopoDS_Compound compound;
            builder.MakeCompound(compound);
            for (const auto& face : faces) {
                builder.Add(compound, face);
            }
            result = compound;
        }
    }
    
    // Apply shape healing
    if (!heal_shape(result)) {
        LOG_WARN("Shape healing had issues, but continuing");
    }
    
    LOG_INFO("B-Rep construction complete");
    return true;
}

bool BRepBuilder::build_with_intersections(
    const std::vector<FittedSurface>& surfaces,
    const std::vector<BoundaryCurve>& boundaries,
    TopoDS_Shape& result
) {
    if (surfaces.empty()) {
        LOG_ERROR("No surfaces provided for B-Rep construction");
        return false;
    }
    
    LOG_INFO("Building B-Rep using surface-surface intersections from ", 
             surfaces.size(), " surfaces and ", boundaries.size(), " boundary curves");
    
    // Step 1: Create OCCT surfaces for all fitted surfaces
    std::map<int, Handle(Geom_Surface)> surface_map;
    for (const auto& surface : surfaces) {
        Handle(Geom_Surface) geom_surface = create_surface(surface);
        if (!geom_surface.IsNull()) {
            surface_map[surface.surface_id] = geom_surface;
            LOG_DEBUG("Created OCCT surface #", surface.surface_id, 
                      " type=", static_cast<int>(surface.type));
        }
    }
    
    // Step 2: Compute surface-surface intersections for all boundaries
    // This gives us exact edges with pcurves on both surfaces
    SurfaceIntersector intersector(config_);
    auto shared_edges = intersector.create_all_shared_edges(
        boundaries, surfaces, surface_map
    );
    
    LOG_INFO("Created shared edges for ", shared_edges.size(), " surface pairs");
    
    // Step 3: Build faces - either trimmed faces if we have edges, or bbox-bounded
    std::vector<TopoDS_Face> faces;
    faces.reserve(surfaces.size());
    
    for (const auto& surface : surfaces) {
        auto it = surface_map.find(surface.surface_id);
        if (it == surface_map.end()) {
            LOG_WARN("No OCCT surface for surface #", surface.surface_id);
            continue;
        }
        
        const Handle(Geom_Surface)& geom_surface = it->second;
        
        // Get all edges for this surface
        auto edges = intersector.get_edges_for_surface(surface.surface_id, shared_edges);
        
        TopoDS_Face face;
        
        if (!edges.empty()) {
            // Try single-wire approach first (for surfaces with single boundary loop)
            auto wire_opt = intersector.build_wire_for_surface(surface.surface_id, shared_edges);
            
            if (wire_opt.has_value() && !wire_opt->IsNull()) {
                // Single wire - simple case
                auto face_opt = intersector.create_trimmed_face(geom_surface, *wire_opt);
                if (face_opt.has_value()) {
                    face = *face_opt;
                    LOG_DEBUG("Created single-wire trimmed face for surface #", surface.surface_id);
                }
            }
            
            // If single wire failed, try multi-wire approach (for cylinders, etc.)
            if (face.IsNull() && edges.size() > 1) {
                auto wires = intersector.build_wires_for_surface(surface.surface_id, shared_edges);
                if (!wires.empty()) {
                    auto face_opt = intersector.create_face_with_wires(geom_surface, wires);
                    if (face_opt.has_value()) {
                        face = *face_opt;
                        LOG_DEBUG("Created multi-wire face for surface #", surface.surface_id,
                                  " with ", wires.size(), " wires");
                    }
                }
            }
        }
        
        // Fallback to bbox-bounded face if trimmed face failed
        if (face.IsNull()) {
            face = create_face_with_bounds(geom_surface, surface);
            if (!face.IsNull()) {
                LOG_DEBUG("Created bbox-bounded face for surface #", surface.surface_id);
            }
        }
        
        if (!face.IsNull()) {
            faces.push_back(face);
        } else {
            LOG_WARN("Failed to create any face for surface #", surface.surface_id);
        }
    }
    
    if (faces.empty()) {
        LOG_ERROR("No valid faces created");
        return false;
    }
    
    LOG_INFO("Created ", faces.size(), " faces");
    
    // Step 4: Sew faces together (or create solid from single closed face)
    if (faces.size() == 1) {
        // Single face - try to make a solid if it's a closed surface (e.g., sphere)
        result = try_make_solid_from_face(faces[0]);
    } else {
        if (!sew_faces(faces, result)) {
            LOG_WARN("Face sewing failed, creating compound of faces");
            BRep_Builder builder;
            TopoDS_Compound compound;
            builder.MakeCompound(compound);
            for (const auto& face : faces) {
                builder.Add(compound, face);
            }
            result = compound;
        }
    }
    
    // Step 5: Apply shape healing
    if (!heal_shape(result)) {
        LOG_WARN("Shape healing had issues, but continuing");
    }
    
    LOG_INFO("B-Rep construction with intersections complete");
    return true;
}

Handle(Geom_Surface) BRepBuilder::create_surface(const FittedSurface& surface) {
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
            LOG_WARN("Unsupported surface type: ", static_cast<int>(surface.type));
            return Handle(Geom_Surface)();
    }
}

Handle(Geom_Surface) BRepBuilder::create_plane(const FittedSurface& surface) {
    // Plane coefficients from PCL: [nx, ny, nz, d] for ax + by + cz + d = 0
    if (surface.coefficients.size() < 4) {
        LOG_ERROR("Plane requires 4 coefficients, got ", surface.coefficients.size());
        return Handle(Geom_Surface)();
    }
    
    double a = surface.coefficients[0];
    double b = surface.coefficients[1];
    double c = surface.coefficients[2];
    double d = surface.coefficients[3];
    
    // Create OCCT plane from equation coefficients
    // Geom_Plane constructor: Geom_Plane(A, B, C, D) creates plane Ax + By + Cz + D = 0
    try {
        Handle(Geom_Plane) plane = new Geom_Plane(a, b, c, d);
        LOG_DEBUG("Created plane with normal (", a, ", ", b, ", ", c, "), d=", d);
        return plane;
    } catch (...) {
        LOG_ERROR("Failed to create Geom_Plane");
        return Handle(Geom_Surface)();
    }
}

Handle(Geom_Surface) BRepBuilder::create_cylinder(const FittedSurface& surface) {
    // Cylinder coefficients from PCL: [point_on_axis.x, .y, .z, axis.x, .y, .z, radius]
    if (surface.coefficients.size() < 7) {
        LOG_ERROR("Cylinder requires 7 coefficients, got ", surface.coefficients.size());
        return Handle(Geom_Surface)();
    }
    
    gp_Pnt axis_point(
        surface.coefficients[0],
        surface.coefficients[1],
        surface.coefficients[2]
    );
    
    gp_Dir axis_dir(
        surface.coefficients[3],
        surface.coefficients[4],
        surface.coefficients[5]
    );
    
    double radius = std::abs(surface.coefficients[6]);
    
    if (radius < 1e-10) {
        LOG_ERROR("Cylinder radius too small: ", radius);
        return Handle(Geom_Surface)();
    }
    
    // Create coordinate system with Z along the cylinder axis
    gp_Ax3 ax3(axis_point, axis_dir);
    
    try {
        Handle(Geom_CylindricalSurface) cylinder = new Geom_CylindricalSurface(ax3, radius);
        LOG_DEBUG("Created cylinder with radius ", radius, 
                  " at (", axis_point.X(), ", ", axis_point.Y(), ", ", axis_point.Z(), ")",
                  " axis (", axis_dir.X(), ", ", axis_dir.Y(), ", ", axis_dir.Z(), ")");
        return cylinder;
    } catch (...) {
        LOG_ERROR("Failed to create Geom_CylindricalSurface");
        return Handle(Geom_Surface)();
    }
}

Handle(Geom_Surface) BRepBuilder::create_sphere(const FittedSurface& surface) {
    // Sphere coefficients from PCL: [center.x, .y, .z, radius]
    if (surface.coefficients.size() < 4) {
        LOG_ERROR("Sphere requires 4 coefficients, got ", surface.coefficients.size());
        return Handle(Geom_Surface)();
    }
    
    gp_Pnt center(
        surface.coefficients[0],
        surface.coefficients[1],
        surface.coefficients[2]
    );
    
    double radius = std::abs(surface.coefficients[3]);
    
    if (radius < 1e-10) {
        LOG_ERROR("Sphere radius too small: ", radius);
        return Handle(Geom_Surface)();
    }
    
    // Create coordinate system centered at the sphere center
    gp_Ax3 ax3(center, gp_Dir(0, 0, 1));
    
    try {
        Handle(Geom_SphericalSurface) sphere = new Geom_SphericalSurface(ax3, radius);
        LOG_DEBUG("Created sphere with radius ", radius,
                  " at (", center.X(), ", ", center.Y(), ", ", center.Z(), ")");
        return sphere;
    } catch (...) {
        LOG_ERROR("Failed to create Geom_SphericalSurface");
        return Handle(Geom_Surface)();
    }
}

Handle(Geom_Surface) BRepBuilder::create_cone(const FittedSurface& surface) {
    // Cone coefficients from PCL: [apex.x, .y, .z, axis.x, .y, .z, opening_angle]
    if (surface.coefficients.size() < 7) {
        LOG_ERROR("Cone requires 7 coefficients, got ", surface.coefficients.size());
        return Handle(Geom_Surface)();
    }
    
    gp_Pnt apex(
        surface.coefficients[0],
        surface.coefficients[1],
        surface.coefficients[2]
    );
    
    gp_Dir axis_dir(
        surface.coefficients[3],
        surface.coefficients[4],
        surface.coefficients[5]
    );
    
    double opening_angle = surface.coefficients[6];  // Half-angle at apex in radians
    
    // OCCT Geom_ConicalSurface uses semi-angle (same as PCL's opening_angle)
    // But the cone surface has its apex at the coordinate system origin,
    // and the axis is the Z direction
    
    if (opening_angle <= 0 || opening_angle >= M_PI / 2) {
        LOG_ERROR("Cone opening angle out of range: ", opening_angle);
        return Handle(Geom_Surface)();
    }
    
    // Create coordinate system with origin at apex and Z along the axis
    gp_Ax3 ax3(apex, axis_dir);
    
    // Geom_ConicalSurface needs a reference radius at a reference location.
    // For OCCT cones, Radius is at the origin of the coordinate system.
    // Since PCL gives us the apex, we use radius=0 at the apex, but OCCT
    // doesn't support that. We need to offset the origin along the axis.
    // Let's use a small reference radius and compute the appropriate offset.
    
    // For a cone, at distance h from apex, radius = h * tan(opening_angle)
    // So if we want radius R at the coordinate origin, we need h = R / tan(opening_angle)
    // Let's pick a reference radius based on the point cloud extent
    
    double ref_radius = 1.0;  // 1mm reference radius
    if (surface.points && !surface.points->empty()) {
        // Find the average distance from apex along axis to determine reference
        double avg_dist = 0.0;
        for (const auto& pt : *surface.points) {
            Eigen::Vector3d p(pt.x, pt.y, pt.z);
            Eigen::Vector3d a(apex.X(), apex.Y(), apex.Z());
            Eigen::Vector3d d(axis_dir.X(), axis_dir.Y(), axis_dir.Z());
            avg_dist += std::abs((p - a).dot(d));
        }
        avg_dist /= surface.points->size();
        ref_radius = avg_dist * std::tan(opening_angle);
    }
    
    if (ref_radius < 1e-6) {
        ref_radius = 1e-6;  // Minimum radius
    }
    
    try {
        Handle(Geom_ConicalSurface) cone = new Geom_ConicalSurface(ax3, opening_angle, ref_radius);
        LOG_DEBUG("Created cone with opening angle ", opening_angle * 180.0 / M_PI, " degrees",
                  " at apex (", apex.X(), ", ", apex.Y(), ", ", apex.Z(), ")",
                  " axis (", axis_dir.X(), ", ", axis_dir.Y(), ", ", axis_dir.Z(), ")");
        return cone;
    } catch (...) {
        LOG_ERROR("Failed to create Geom_ConicalSurface");
        return Handle(Geom_Surface)();
    }
}

TopoDS_Face BRepBuilder::create_face_with_bounds(const Handle(Geom_Surface)& surface, const FittedSurface& fitted) {
    if (surface.IsNull()) {
        return TopoDS_Face();
    }
    
    // Compute bounds from the point cloud by projecting points onto the surface
    // and finding the UV parameter extents
    double umin = std::numeric_limits<double>::max();
    double umax = std::numeric_limits<double>::lowest();
    double vmin = std::numeric_limits<double>::max();
    double vmax = std::numeric_limits<double>::lowest();
    
    if (fitted.points && !fitted.points->empty()) {
        // Project points onto surface to get UV parameters
        for (const auto& pt : *fitted.points) {
            gp_Pnt point(pt.x, pt.y, pt.z);
            
            try {
                GeomAPI_ProjectPointOnSurf projector(point, surface);
                if (projector.NbPoints() > 0) {
                    double u, v;
                    projector.LowerDistanceParameters(u, v);
                    umin = std::min(umin, u);
                    umax = std::max(umax, u);
                    vmin = std::min(vmin, v);
                    vmax = std::max(vmax, v);
                }
            } catch (...) {
                // Skip points that fail to project
            }
        }
        
        // If projection failed for all points, fall back to surface-type defaults
        if (umin > umax || vmin > vmax) {
            LOG_DEBUG("UV projection failed, using surface-type defaults");
            
            double xmin = std::numeric_limits<double>::max();
            double xmax = std::numeric_limits<double>::lowest();
            double ymin = std::numeric_limits<double>::max();
            double ymax = std::numeric_limits<double>::lowest();
            double zmin = std::numeric_limits<double>::max();
            double zmax = std::numeric_limits<double>::lowest();
            
            for (const auto& pt : *fitted.points) {
                xmin = std::min(xmin, static_cast<double>(pt.x));
                xmax = std::max(xmax, static_cast<double>(pt.x));
                ymin = std::min(ymin, static_cast<double>(pt.y));
                ymax = std::max(ymax, static_cast<double>(pt.y));
                zmin = std::min(zmin, static_cast<double>(pt.z));
                zmax = std::max(zmax, static_cast<double>(pt.z));
            }
            
            double size = std::max({xmax - xmin, ymax - ymin, zmax - zmin});
            
            switch (fitted.type) {
                case SurfaceType::PLANE:
                    umin = -size;
                    umax = size;
                    vmin = -size;
                    vmax = size;
                    break;
                    
                case SurfaceType::CYLINDER:
                    umin = 0;
                    umax = 2.0 * M_PI;
                    vmin = -size;
                    vmax = size;
                    break;
                    
                case SurfaceType::SPHERE:
                    umin = 0;
                    umax = 2.0 * M_PI;
                    vmin = -M_PI / 2.0;
                    vmax = M_PI / 2.0;
                    break;
                    
                case SurfaceType::CONE:
                    umin = 0;
                    umax = 2.0 * M_PI;
                    vmin = 0;
                    vmax = size * 2.0;
                    break;
                    
                default:
                    umin = -size;
                    umax = size;
                    vmin = -size;
                    vmax = size;
                    break;
            }
        }
        
        LOG_DEBUG("Surface #", fitted.surface_id, " UV bounds from projection: ",
                  "U=[", umin, ", ", umax, "], V=[", vmin, ", ", vmax, "]");
    } else {
        // Default bounds
        umin = -100.0;
        umax = 100.0;
        vmin = -100.0;
        vmax = 100.0;
    }
    
    // Add a small margin to ensure we don't clip points at the boundary
    double u_range = umax - umin;
    double v_range = vmax - vmin;
    
    // Use per-dimension margin to respect partial fits (e.g. thin strips)
    // Ensure a minimum margin to handle degenerate/noisy cases
    double u_margin = std::max(0.01 * u_range, 1e-4);
    double v_margin = std::max(0.01 * v_range, 1e-4);
    
    umin -= u_margin;
    umax += u_margin;
    vmin -= v_margin;
    vmax += v_margin;
    
    // Clamp bounds to valid ranges for specific surface types
    if (fitted.type == SurfaceType::SPHERE) {
        // V must be in [-PI/2, PI/2]
        vmin = std::max(vmin, -M_PI / 2.0);
        vmax = std::min(vmax, M_PI / 2.0);
        
        // If U covers full circle (approx), snap to full range
        if (umax - umin >= 2.0 * M_PI - 0.1) {
            umin = 0.0;
            umax = 2.0 * M_PI;
        }
    } else if (fitted.type == SurfaceType::CYLINDER || fitted.type == SurfaceType::CONE) {
        // U is periodic 0..2PI
        if (umax - umin >= 2.0 * M_PI - 0.1) {
            umin = 0.0;
            umax = 2.0 * M_PI;
        }
    }
    
    try {
        // Create face with bounded parameters
        BRepBuilderAPI_MakeFace make_face(surface, umin, umax, vmin, vmax, 1e-6);
        
        if (!make_face.IsDone()) {
            LOG_ERROR("BRepBuilderAPI_MakeFace failed for surface #", fitted.surface_id);
            return TopoDS_Face();
        }
        
        return make_face.Face();
    } catch (Standard_Failure const& e) {
        LOG_ERROR("OCCT exception while creating face: ", e.GetMessageString());
        return TopoDS_Face();
    } catch (...) {
        LOG_ERROR("Unknown exception while creating face");
        return TopoDS_Face();
    }
}

std::vector<TopoDS_Wire> BRepBuilder::build_boundary_wires(
    int surface_id,
    const Handle(Geom_Surface)& surface,
    const std::vector<BoundaryCurve>& boundaries)
{
    std::vector<TopoDS_Wire> wires;
    
    // Find all boundary curves that touch this surface
    for (const auto& curve : boundaries) {
        if (curve.surface_id_left == surface_id || curve.surface_id_right == surface_id) {
            TopoDS_Wire wire = create_wire_from_curve(curve, surface);
            if (!wire.IsNull()) {
                wires.push_back(wire);
            }
        }
    }
    
    LOG_DEBUG("Built ", wires.size(), " boundary wires for surface #", surface_id);
    return wires;
}

TopoDS_Wire BRepBuilder::create_wire_from_curve(
    const BoundaryCurve& curve,
    const Handle(Geom_Surface)& surface)
{
    (void)surface;  // Will be used later for pcurve computation
    
    if (curve.points.size() < 2) {
        return TopoDS_Wire();
    }
    
    try {
        BRepBuilderAPI_MakeWire wire_maker;
        
        // Create edges from consecutive points
        for (size_t i = 0; i < curve.points.size() - 1; ++i) {
            const auto& p1 = curve.points[i];
            const auto& p2 = curve.points[i + 1];
            
            gp_Pnt gp_p1(p1.x(), p1.y(), p1.z());
            gp_Pnt gp_p2(p2.x(), p2.y(), p2.z());
            
            // Skip degenerate edges
            if (gp_p1.Distance(gp_p2) < 1e-9) {
                continue;
            }
            
            TopoDS_Edge edge = create_edge(gp_p1, gp_p2);
            if (!edge.IsNull()) {
                wire_maker.Add(edge);
            }
        }
        
        // Check if the curve is closed (first and last points close together)
        const auto& first = curve.points.front();
        const auto& last = curve.points.back();
        gp_Pnt gp_first(first.x(), first.y(), first.z());
        gp_Pnt gp_last(last.x(), last.y(), last.z());
        
        if (gp_first.Distance(gp_last) < config_.sewing_tolerance * 10) {
            // Close the wire by adding an edge from last to first
            TopoDS_Edge closing_edge = create_edge(gp_last, gp_first);
            if (!closing_edge.IsNull()) {
                wire_maker.Add(closing_edge);
            }
        }
        
        if (!wire_maker.IsDone()) {
            LOG_DEBUG("Wire construction incomplete for boundary curve");
            return TopoDS_Wire();
        }
        
        TopoDS_Wire wire = wire_maker.Wire();
        
        // Basic wire fixing (without face context - pcurves will be added later)
        // We only do basic connectivity fixes here
        ShapeFix_Wire wire_fixer;
        wire_fixer.Load(wire);
        wire_fixer.SetPrecision(config_.sewing_tolerance);
        wire_fixer.FixConnected(config_.sewing_tolerance);
        wire = wire_fixer.Wire();
        
        return wire;
    } catch (...) {
        LOG_DEBUG("Exception creating wire from boundary curve");
        return TopoDS_Wire();
    }
}

TopoDS_Edge BRepBuilder::create_edge(const gp_Pnt& p1, const gp_Pnt& p2) {
    try {
        BRepBuilderAPI_MakeEdge edge_maker(p1, p2);
        if (edge_maker.IsDone()) {
            return edge_maker.Edge();
        }
    } catch (...) {
        // Ignore exceptions from degenerate edges
    }
    return TopoDS_Edge();
}

TopoDS_Face BRepBuilder::create_trimmed_face(
    const Handle(Geom_Surface)& surface,
    const FittedSurface& fitted,
    const std::vector<TopoDS_Wire>& boundary_wires)
{
    if (surface.IsNull() || boundary_wires.empty()) {
        return TopoDS_Face();
    }
    
    try {
        // Create a base face from the surface
        // For non-planar surfaces, we create an unbounded/naturally-bounded face first
        BRepBuilderAPI_MakeFace face_maker(surface, config_.sewing_tolerance);
        
        if (!face_maker.IsDone()) {
            LOG_DEBUG("Failed to create base face for surface #", fitted.surface_id);
            return TopoDS_Face();
        }
        
        TopoDS_Face base_face = face_maker.Face();
        
        // For each boundary wire, add pcurves to its edges using ShapeFix_Wire
        // This projects the 3D edges onto the surface to create 2D pcurves
        std::vector<TopoDS_Wire> fixed_wires;
        
        for (const auto& wire : boundary_wires) {
            if (wire.IsNull()) continue;
            
            // Use ShapeFix_Wire to add pcurves to the wire edges
            ShapeFix_Wire wire_fixer(wire, base_face, config_.sewing_tolerance);
            
            // Enable pcurve fixing
            wire_fixer.FixAddPCurveMode() = 1;  // Add missing pcurves
            wire_fixer.FixEdgeCurvesMode() = 1; // Fix edge curves including pcurves
            wire_fixer.FixShiftedMode() = 1;    // Fix pcurves shifted on closed surfaces (cylinders)
            
            // Run fixes
            wire_fixer.FixEdgeCurves();
            
            // Get the fixed wire
            TopoDS_Wire fixed_wire = wire_fixer.Wire();
            if (!fixed_wire.IsNull()) {
                fixed_wires.push_back(fixed_wire);
                LOG_DEBUG("Fixed wire for surface #", fitted.surface_id, 
                         " - pcurves added");
            }
        }
        
        if (fixed_wires.empty()) {
            LOG_DEBUG("No wires could be fixed for surface #", fitted.surface_id);
            return TopoDS_Face();
        }
        
        // Now try to create a face with the fixed boundary wires
        // Create a new face with the wires as bounds
        BRepBuilderAPI_MakeFace bounded_face_maker(surface, config_.sewing_tolerance);
        
        if (!bounded_face_maker.IsDone()) {
            LOG_DEBUG("Failed to create bounded face for surface #", fitted.surface_id);
            return TopoDS_Face();
        }
        
        for (const auto& wire : fixed_wires) {
            bounded_face_maker.Add(wire);
            
            if (bounded_face_maker.Error() != BRepBuilderAPI_FaceDone) {
                LOG_DEBUG("Failed to add fixed wire to face, error: ", 
                         static_cast<int>(bounded_face_maker.Error()));
                // Continue trying other wires
            }
        }
        
        if (!bounded_face_maker.IsDone()) {
            LOG_DEBUG("Face construction incomplete after adding fixed wires");
            return TopoDS_Face();
        }
        
        TopoDS_Face face = bounded_face_maker.Face();
        
        // Apply face fixing
        ShapeFix_Face face_fixer(face);
        face_fixer.SetPrecision(config_.sewing_tolerance);
        face_fixer.FixAddNaturalBoundMode() = 1;  // Add natural bounds for periodic surfaces
        face_fixer.Perform();
        
        if (face_fixer.Status(ShapeExtend_DONE)) {
            face = face_fixer.Face();
        }
        
        return face;
    } catch (Standard_Failure const& e) {
        LOG_DEBUG("OCCT exception creating trimmed face: ", e.GetMessageString());
        return TopoDS_Face();
    } catch (...) {
        LOG_DEBUG("Unknown exception creating trimmed face");
        return TopoDS_Face();
    }
}

TopoDS_Shape BRepBuilder::try_make_solid_from_face(const TopoDS_Face& face) {
    if (face.IsNull()) {
        return face;
    }
    
    try {
        // Build a shell from the single face
        BRep_Builder builder;
        TopoDS_Shell shell;
        builder.MakeShell(shell);
        builder.Add(shell, face);
        
        // Fix the shell
        ShapeFix_Shell shell_fixer(shell);
        shell_fixer.SetPrecision(config_.sewing_tolerance);
        shell_fixer.Perform();
        
        if (shell_fixer.Status(ShapeExtend_DONE)) {
            shell = shell_fixer.Shell();
        }
        
        // Check if the shell is closed (e.g., a complete sphere)
        BRepCheck_Shell shell_checker(shell);
        if (shell_checker.Closed() == BRepCheck_NoError) {
            LOG_DEBUG("Single face forms closed shell, creating solid");
            
            // Create a solid from the closed shell
            ShapeFix_Solid solid_fixer;
            solid_fixer.SetPrecision(config_.sewing_tolerance);
            TopoDS_Solid solid = solid_fixer.SolidFromShell(shell);
            
            if (!solid.IsNull()) {
                LOG_DEBUG("Created solid from single closed face");
                return solid;
            }
        }
        
        // Not a closed surface, return the original face
        LOG_DEBUG("Single face is not closed, keeping as face");
        return face;
    } catch (Standard_Failure const& e) {
        LOG_DEBUG("Exception trying to make solid from face: ", e.GetMessageString());
        return face;
    } catch (...) {
        LOG_DEBUG("Unknown exception trying to make solid from face");
        return face;
    }
}

bool BRepBuilder::sew_faces(const std::vector<TopoDS_Face>& faces, TopoDS_Shape& result) {
    if (faces.empty()) {
        return false;
    }
    
    if (faces.size() == 1) {
        result = faces[0];
        return true;
    }
    
    try {
        // Use a more generous sewing tolerance to connect faces whose boundary
        // edges may not match exactly (e.g., polyline approximation of circles)
        double sewing_tol = config_.sewing_tolerance;
        
        // Adaptive tolerance based on mesh characteristics
        // For curved surfaces meeting planes, the boundary curves are approximated
        // as polylines which won't exactly match the analytic curves
        sewing_tol = std::max(sewing_tol, 0.1);  // At least 0.1mm tolerance
        
        BRepBuilderAPI_Sewing sewing(sewing_tol);
        sewing.SetNonManifoldMode(false);  // We want a manifold solid
        sewing.SetFloatingEdgesMode(true);  // Allow floating edges to be sewn
        
        for (const auto& face : faces) {
            sewing.Add(face);
        }
        
        sewing.Perform();
        
        result = sewing.SewedShape();
        
        int nb_free_edges = sewing.NbFreeEdges();
        int nb_multiple_edges = sewing.NbMultipleEdges();
        int nb_degenerated = sewing.NbDegeneratedShapes();
        
        LOG_DEBUG("Sewing complete: ", nb_free_edges, " free edges, ",
                 nb_multiple_edges, " multiple edges, ",
                 nb_degenerated, " degenerated shapes, tolerance=", sewing_tol);
        
        // If sewing produced a shell, try to make it into a solid
        if (!result.IsNull()) {
            TopoDS_Shell shell;
            bool found_shell = false;
            
            // Look for a shell in the result
            for (TopExp_Explorer exp(result, TopAbs_SHELL); exp.More(); exp.Next()) {
                shell = TopoDS::Shell(exp.Current());
                found_shell = true;
                break;
            }
            
            if (found_shell) {
                // Fix the shell orientation and try to make a solid
                ShapeFix_Shell shell_fixer(shell);
                shell_fixer.SetPrecision(config_.sewing_tolerance);
                shell_fixer.Perform();
                
                if (shell_fixer.Status(ShapeExtend_DONE)) {
                    shell = shell_fixer.Shell();
                    LOG_DEBUG("Shell fixing applied");
                }
                
                // Check if shell is closed
                BRepCheck_Shell shell_checker(shell);
                if (shell_checker.Closed() == BRepCheck_NoError) {
                    LOG_DEBUG("Shell is closed, creating solid");
                    
                    // Try to create a solid from the closed shell
                    ShapeFix_Solid solid_fixer;
                    solid_fixer.SetPrecision(config_.sewing_tolerance);
                    TopoDS_Solid solid = solid_fixer.SolidFromShell(shell);
                    
                    if (!solid.IsNull()) {
                        result = solid;
                        LOG_DEBUG("Created solid from shell");
                    }
                } else {
                    LOG_DEBUG("Shell is not closed, keeping as shell");
                    result = shell;
                }
            }
        }
        
        return !result.IsNull();
    } catch (...) {
        LOG_ERROR("Exception during face sewing");
        return false;
    }
}

bool BRepBuilder::heal_shape(TopoDS_Shape& shape) {
    if (shape.IsNull()) {
        return false;
    }
    
    try {
        TopAbs_ShapeEnum original_type = shape.ShapeType();
        LOG_DEBUG("Shape type before healing: ", static_cast<int>(original_type));

        ShapeFix_Shape fixer(shape);
        fixer.SetPrecision(config_.sewing_tolerance);
        fixer.Perform();

        TopoDS_Shape fixed = fixer.Shape();
        if (!fixed.IsNull()) {
            TopAbs_ShapeEnum fixed_type = fixed.ShapeType();
            LOG_DEBUG("Shape type after healing: ", static_cast<int>(fixed_type));

            // Warn if healing made any change to the shape
            if (!shape.IsSame(fixed)) {
                LOG_WARN("Shape healing modified the geometry. This indicates a problem in previous modeling steps. Please investigate and fix the root cause.");
            }

            // Don't allow healing to demote a solid to a shell
            // This can happen if ShapeFix finds issues but we'd rather keep the solid
            if (original_type == TopAbs_SOLID && fixed_type != TopAbs_SOLID) {
                LOG_WARN("Shape healing demoted solid to ", static_cast<int>(fixed_type), 
                         ", keeping original solid");
                return true;  // Keep original shape
            }

            shape = fixed;
            LOG_DEBUG("Shape healing applied successfully");
            return true;
        }

        return false;
    } catch (...) {
        LOG_ERROR("Exception during shape healing");
        return false;
    }
}

} // namespace brepper
