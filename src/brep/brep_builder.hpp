#pragma once

#include "common/types.hpp"
#include "common/config.hpp"
#include "surface_intersector.hpp"

#include <TopoDS_Shape.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Wire.hxx>
#include <TopoDS_Edge.hxx>
#include <Geom_Surface.hxx>
#include <Geom_Curve.hxx>
#include <vector>
#include <map>

namespace brepper {

class BRepBuilder {
public:
    explicit BRepBuilder(const Config& config);
    
    // Build a B-Rep solid from fitted surfaces and boundary curves
    bool build(
        const std::vector<FittedSurface>& surfaces,
        const std::vector<BoundaryCurve>& boundaries,
        TopoDS_Shape& result
    );
    
    // Build using surface-surface intersection for exact edges
    bool build_with_intersections(
        const std::vector<FittedSurface>& surfaces,
        const std::vector<BoundaryCurve>& boundaries,
        TopoDS_Shape& result
    );
    
    // Create OCCT surface from fitted surface (public for reuse)
    Handle(Geom_Surface) create_surface(const FittedSurface& surface);
    
private:
    const Config& config_;
    
    // Create surface for each type
    Handle(Geom_Surface) create_plane(const FittedSurface& surface);
    Handle(Geom_Surface) create_cylinder(const FittedSurface& surface);
    Handle(Geom_Surface) create_sphere(const FittedSurface& surface);
    Handle(Geom_Surface) create_cone(const FittedSurface& surface);
    
    // Create a face from a surface with bounding box bounds (fallback)
    TopoDS_Face create_face_with_bounds(const Handle(Geom_Surface)& surface, 
                                         const FittedSurface& fitted);
    
    // Create a trimmed face from a surface and boundary wires
    TopoDS_Face create_trimmed_face(const Handle(Geom_Surface)& surface,
                                     const FittedSurface& fitted,
                                     const std::vector<TopoDS_Wire>& boundary_wires);
    
    // Build boundary wires for a surface from boundary curves
    std::vector<TopoDS_Wire> build_boundary_wires(
        int surface_id,
        const Handle(Geom_Surface)& surface,
        const std::vector<BoundaryCurve>& boundaries);
    
    // Create a wire from a boundary curve's points
    TopoDS_Wire create_wire_from_curve(const BoundaryCurve& curve,
                                        const Handle(Geom_Surface)& surface);
    
    // Create an edge from two 3D points (as a line segment)
    TopoDS_Edge create_edge(const gp_Pnt& p1, const gp_Pnt& p2);
    
    // Sew faces together into a shell/solid
    bool sew_faces(const std::vector<TopoDS_Face>& faces, TopoDS_Shape& result);
    
    // Try to make a solid from a single face (for closed surfaces like spheres)
    TopoDS_Shape try_make_solid_from_face(const TopoDS_Face& face);
    
    // Apply shape healing/fixing
    bool heal_shape(TopoDS_Shape& shape);
};

} // namespace brepper
