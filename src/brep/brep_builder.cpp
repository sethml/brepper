#include "brep_builder.hpp"
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
#include <ShapeFix_Shape.hxx>
#include <ShapeFix_Shell.hxx>
#include <ShapeFix_Solid.hxx>
#include <ShapeFix_Wire.hxx>
#include <ShapeFix_Face.hxx>
#include <ShapeExtend_Status.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Shell.hxx>
#include <TopoDS_Solid.hxx>
#include <BRep_Builder.hxx>
#include <GeomAPI_ProjectPointOnSurf.hxx>
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
        
        // Build boundary wires for this surface
        std::vector<TopoDS_Wire> boundary_wires = build_boundary_wires(
            surface.surface_id, geom_surface, boundaries);
        
        TopoDS_Face face;
        if (!boundary_wires.empty()) {
            // Try to create trimmed face with boundary wires
            face = create_trimmed_face(geom_surface, surface, boundary_wires);
            if (!face.IsNull()) {
                LOG_DEBUG("Created trimmed face for surface #", surface.surface_id,
                         " with ", boundary_wires.size(), " boundary wires");
            }
        }
        
        // Fall back to bounding box face if trimmed face creation failed
        if (face.IsNull()) {
            face = create_face_with_bounds(geom_surface, surface);
            if (!face.IsNull()) {
                LOG_DEBUG("Created bbox-bounded face for surface #", surface.surface_id);
            }
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
    
    // Compute bounds from the point cloud
    double umin = std::numeric_limits<double>::max();
    double umax = std::numeric_limits<double>::lowest();
    double vmin = std::numeric_limits<double>::max();
    double vmax = std::numeric_limits<double>::lowest();
    
    // If we have points, compute parametric bounds
    if (fitted.points && !fitted.points->empty()) {
        // For now, use a simpler approach: compute bounding box in 3D
        // and use that to estimate parametric bounds
        
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
        
        // Use bounding box size for parametric bounds (rough approximation)
        double size = std::max({xmax - xmin, ymax - ymin, zmax - zmin});
        
        switch (fitted.type) {
            case SurfaceType::PLANE:
                // Planes: parametric u,v correspond to distances along XDir/YDir
                umin = -size;
                umax = size;
                vmin = -size;
                vmax = size;
                break;
                
            case SurfaceType::CYLINDER:
                // Cylinder: u is angle [0, 2*pi], v is distance along axis
                umin = 0;
                umax = 2.0 * M_PI;
                vmin = -size;
                vmax = size;
                break;
                
            case SurfaceType::SPHERE:
                // Sphere: u is longitude [0, 2*pi], v is latitude [-pi/2, pi/2]
                umin = 0;
                umax = 2.0 * M_PI;
                vmin = -M_PI / 2.0;
                vmax = M_PI / 2.0;
                break;
                
            case SurfaceType::CONE:
                // Cone: u is angle [0, 2*pi], v is distance from apex
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
    } else {
        // Default bounds
        umin = -100.0;
        umax = 100.0;
        vmin = -100.0;
        vmax = 100.0;
    }
    
    // Add some margin
    double margin = 0.1 * std::max(umax - umin, vmax - vmin);
    umin -= margin;
    umax += margin;
    vmin -= margin;
    vmax += margin;
    
    try {
        // Create face with bounded parameters
        BRepBuilderAPI_MakeFace make_face(surface, umin, umax, vmin, vmax, 1e-6);
        
        if (!make_face.IsDone()) {
            LOG_ERROR("BRepBuilderAPI_MakeFace failed");
            return TopoDS_Face();
        }
        
        return make_face.Face();
    } catch (...) {
        LOG_ERROR("Exception while creating face");
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
    const Handle(Geom_Surface)& /*surface*/)
{
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
        
        // Apply wire fixing
        ShapeFix_Wire wire_fixer(wire, TopoDS_Face(), config_.sewing_tolerance);
        wire_fixer.Perform();
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
    const FittedSurface& /*fitted*/,
    const std::vector<TopoDS_Wire>& boundary_wires)
{
    if (surface.IsNull() || boundary_wires.empty()) {
        return TopoDS_Face();
    }
    
    try {
        // First create an unbounded face on the surface
        BRepBuilderAPI_MakeFace face_maker(surface, config_.sewing_tolerance);
        
        if (!face_maker.IsDone()) {
            LOG_DEBUG("Failed to create base face for trimming");
            return TopoDS_Face();
        }
        
        // Try to add boundary wires as face boundaries
        for (const auto& wire : boundary_wires) {
            if (wire.IsNull()) continue;
            
            // Add the wire to the face
            face_maker.Add(wire);
            
            if (face_maker.Error() != BRepBuilderAPI_FaceDone) {
                LOG_DEBUG("Failed to add boundary wire to face, error: ", 
                         static_cast<int>(face_maker.Error()));
                // Continue trying other wires
            }
        }
        
        if (!face_maker.IsDone()) {
            LOG_DEBUG("Face construction incomplete after adding wires");
            return TopoDS_Face();
        }
        
        TopoDS_Face face = face_maker.Face();
        
        // Apply face fixing
        ShapeFix_Face face_fixer(face);
        face_fixer.SetPrecision(config_.sewing_tolerance);
        face_fixer.Perform();
        
        if (face_fixer.Status(ShapeExtend_DONE)) {
            face = face_fixer.Face();
        }
        
        return face;
    } catch (...) {
        LOG_DEBUG("Exception creating trimmed face");
        return TopoDS_Face();
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
        BRepBuilderAPI_Sewing sewing(config_.sewing_tolerance);
        
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
                  nb_degenerated, " degenerated shapes");
        
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
        ShapeFix_Shape fixer(shape);
        fixer.SetPrecision(config_.sewing_tolerance);
        fixer.Perform();
        
        TopoDS_Shape fixed = fixer.Shape();
        if (!fixed.IsNull()) {
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
