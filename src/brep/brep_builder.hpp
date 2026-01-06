#pragma once

#include "common/types.hpp"
#include "common/config.hpp"

#include <TopoDS_Shape.hxx>
#include <TopoDS_Face.hxx>
#include <Geom_Surface.hxx>
#include <Geom_Curve.hxx>
#include <vector>

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
    
private:
    const Config& config_;
    
    // Convert a FittedSurface to an OCCT Geom_Surface
    Handle(Geom_Surface) create_surface(const FittedSurface& surface);
    
    // Create surface for each type
    Handle(Geom_Surface) create_plane(const FittedSurface& surface);
    Handle(Geom_Surface) create_cylinder(const FittedSurface& surface);
    Handle(Geom_Surface) create_sphere(const FittedSurface& surface);
    Handle(Geom_Surface) create_cone(const FittedSurface& surface);
    
    // Create a face from a surface (untrimmed, infinite bounds clamped)
    TopoDS_Face create_face(const Handle(Geom_Surface)& surface, const FittedSurface& fitted);
    
    // Sew faces together into a shell/solid
    bool sew_faces(const std::vector<TopoDS_Face>& faces, TopoDS_Shape& result);
    
    // Apply shape healing/fixing
    bool heal_shape(TopoDS_Shape& shape);
};

} // namespace brepper
