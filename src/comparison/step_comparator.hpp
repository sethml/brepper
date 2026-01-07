#pragma once

#include <string>
#include <vector>
#include <optional>

#include <TopoDS_Shape.hxx>
#include <gp_Pnt.hxx>

namespace brepper {

/// Result of comparing two B-Rep shapes
struct ShapeComparisonResult {
    bool passed = false;
    
    // Topology comparison
    int ref_vertices = 0;
    int gen_vertices = 0;
    int ref_edges = 0;
    int gen_edges = 0;
    int ref_faces = 0;
    int gen_faces = 0;
    int ref_shells = 0;
    int gen_shells = 0;
    int ref_solids = 0;
    int gen_solids = 0;
    
    // Geometric comparison
    double ref_volume = 0.0;
    double gen_volume = 0.0;
    double volume_error_percent = 0.0;
    
    double ref_surface_area = 0.0;
    double gen_surface_area = 0.0;
    double area_error_percent = 0.0;
    
    // Bounding box comparison
    double ref_bbox_diagonal = 0.0;
    double gen_bbox_diagonal = 0.0;
    double bbox_error_percent = 0.0;
    
    // Center of mass comparison
    gp_Pnt ref_centroid{0, 0, 0};
    gp_Pnt gen_centroid{0, 0, 0};
    double centroid_distance = 0.0;
    
    // Tolerance used
    double tolerance = 0.0;
    
    // Error messages
    std::vector<std::string> errors;
    std::vector<std::string> warnings;
    
    /// Get a human-readable summary
    std::string summary() const;
};

/// Utility class for comparing OCCT B-Rep shapes
class STEPComparator {
public:
    /// Set comparison tolerance (relative to characteristic length)
    void set_tolerance(double tol) { tolerance_ = tol; }
    
    /// Read a STEP file and return the shape
    std::optional<TopoDS_Shape> read_step(const std::string& filename);
    
    /// Compare two shapes
    ShapeComparisonResult compare(const TopoDS_Shape& reference, 
                                  const TopoDS_Shape& generated);
    
    /// Full comparison: read both files and compare
    ShapeComparisonResult compare_files(const std::string& reference_step,
                                        const std::string& generated_step);
    
    /// Count topology entities in a shape
    static void count_topology(const TopoDS_Shape& shape,
                               int& vertices, int& edges, int& faces,
                               int& shells, int& solids);
    
    /// Compute volume of a shape (returns 0 for non-solid shapes)
    static double compute_volume(const TopoDS_Shape& shape);
    
    /// Compute surface area of a shape
    static double compute_surface_area(const TopoDS_Shape& shape);
    
    /// Compute bounding box diagonal
    static double compute_bbox_diagonal(const TopoDS_Shape& shape);
    
    /// Compute centroid (center of mass) of a shape
    static gp_Pnt compute_centroid(const TopoDS_Shape& shape);

private:
    double tolerance_ = 0.01;  // 1% relative tolerance by default
};

} // namespace brepper
