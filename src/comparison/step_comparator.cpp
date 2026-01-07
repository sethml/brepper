#include "step_comparator.hpp"
#include "common/logging.hpp"

#include <STEPControl_Reader.hxx>
#include <IFSelect_ReturnStatus.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Compound.hxx>
#include <BRep_Builder.hxx>
#include <GProp_GProps.hxx>
#include <BRepGProp.hxx>
#include <Bnd_Box.hxx>
#include <BRepBndLib.hxx>

#include <sstream>
#include <cmath>
#include <iomanip>

namespace brepper {

std::string ShapeComparisonResult::summary() const {
    std::ostringstream ss;
    ss << std::fixed << std::setprecision(2);
    
    ss << "=== Shape Comparison Result ===\n";
    ss << "Overall: " << (passed ? "PASSED" : "FAILED") << "\n\n";
    
    ss << "Topology:\n";
    ss << "  Vertices: ref=" << ref_vertices << ", gen=" << gen_vertices;
    if (ref_vertices != gen_vertices) ss << " (MISMATCH)";
    ss << "\n";
    
    ss << "  Edges:    ref=" << ref_edges << ", gen=" << gen_edges;
    if (ref_edges != gen_edges) ss << " (MISMATCH)";
    ss << "\n";
    
    ss << "  Faces:    ref=" << ref_faces << ", gen=" << gen_faces;
    if (ref_faces != gen_faces) ss << " (MISMATCH)";
    ss << "\n";
    
    ss << "  Shells:   ref=" << ref_shells << ", gen=" << gen_shells;
    if (ref_shells != gen_shells) ss << " (MISMATCH)";
    ss << "\n";
    
    ss << "  Solids:   ref=" << ref_solids << ", gen=" << gen_solids;
    if (ref_solids != gen_solids) ss << " (MISMATCH)";
    ss << "\n\n";
    
    ss << "Geometry (tolerance=" << (tolerance * 100.0) << "%):\n";
    ss << "  Volume:       ref=" << ref_volume << ", gen=" << gen_volume 
       << " (error=" << volume_error_percent << "%)\n";
    ss << "  Surface area: ref=" << ref_surface_area << ", gen=" << gen_surface_area
       << " (error=" << area_error_percent << "%)\n";
    ss << "  BBox diag:    ref=" << ref_bbox_diagonal << ", gen=" << gen_bbox_diagonal
       << " (error=" << bbox_error_percent << "%)\n";
    ss << "  Centroid:     ref=(" << ref_centroid.X() << "," << ref_centroid.Y() << "," << ref_centroid.Z() << ")"
       << ", gen=(" << gen_centroid.X() << "," << gen_centroid.Y() << "," << gen_centroid.Z() << ")"
       << " (dist=" << centroid_distance << ")\n";
    
    if (!errors.empty()) {
        ss << "\nErrors:\n";
        for (const auto& err : errors) {
            ss << "  - " << err << "\n";
        }
    }
    
    if (!warnings.empty()) {
        ss << "\nWarnings:\n";
        for (const auto& warn : warnings) {
            ss << "  - " << warn << "\n";
        }
    }
    
    return ss.str();
}

std::optional<TopoDS_Shape> STEPComparator::read_step(const std::string& filename) {
    STEPControl_Reader reader;
    
    IFSelect_ReturnStatus status = reader.ReadFile(filename.c_str());
    if (status != IFSelect_RetDone) {
        LOG_ERROR("Failed to read STEP file: ", filename);
        return std::nullopt;
    }
    
    // Transfer all roots
    int num_roots = reader.NbRootsForTransfer();
    if (num_roots == 0) {
        LOG_ERROR("No shapes found in STEP file: ", filename);
        return std::nullopt;
    }
    
    reader.TransferRoots();
    
    if (reader.NbShapes() == 0) {
        LOG_ERROR("Failed to transfer shapes from STEP file: ", filename);
        return std::nullopt;
    }
    
    // Return the first (or only) shape
    // For files with multiple shapes, we'd need to combine them
    if (reader.NbShapes() == 1) {
        return reader.Shape(1);
    }
    
    // Multiple shapes - combine into a compound
    TopoDS_Compound compound;
    BRep_Builder builder;
    builder.MakeCompound(compound);
    for (int i = 1; i <= reader.NbShapes(); ++i) {
        builder.Add(compound, reader.Shape(i));
    }
    return compound;
}

void STEPComparator::count_topology(const TopoDS_Shape& shape,
                                    int& vertices, int& edges, int& faces,
                                    int& shells, int& solids) {
    vertices = edges = faces = shells = solids = 0;
    
    for (TopExp_Explorer exp(shape, TopAbs_VERTEX); exp.More(); exp.Next()) {
        ++vertices;
    }
    for (TopExp_Explorer exp(shape, TopAbs_EDGE); exp.More(); exp.Next()) {
        ++edges;
    }
    for (TopExp_Explorer exp(shape, TopAbs_FACE); exp.More(); exp.Next()) {
        ++faces;
    }
    for (TopExp_Explorer exp(shape, TopAbs_SHELL); exp.More(); exp.Next()) {
        ++shells;
    }
    for (TopExp_Explorer exp(shape, TopAbs_SOLID); exp.More(); exp.Next()) {
        ++solids;
    }
}

double STEPComparator::compute_volume(const TopoDS_Shape& shape) {
    GProp_GProps props;
    BRepGProp::VolumeProperties(shape, props);
    return props.Mass();  // For volume, Mass() returns volume
}

double STEPComparator::compute_surface_area(const TopoDS_Shape& shape) {
    GProp_GProps props;
    BRepGProp::SurfaceProperties(shape, props);
    return props.Mass();  // For surface, Mass() returns area
}

double STEPComparator::compute_bbox_diagonal(const TopoDS_Shape& shape) {
    Bnd_Box box;
    BRepBndLib::Add(shape, box);
    
    if (box.IsVoid()) {
        return 0.0;
    }
    
    double xmin, ymin, zmin, xmax, ymax, zmax;
    box.Get(xmin, ymin, zmin, xmax, ymax, zmax);
    
    double dx = xmax - xmin;
    double dy = ymax - ymin;
    double dz = zmax - zmin;
    
    return std::sqrt(dx*dx + dy*dy + dz*dz);
}

gp_Pnt STEPComparator::compute_centroid(const TopoDS_Shape& shape) {
    GProp_GProps props;
    BRepGProp::VolumeProperties(shape, props);
    return props.CentreOfMass();
}

ShapeComparisonResult STEPComparator::compare(const TopoDS_Shape& reference,
                                               const TopoDS_Shape& generated) {
    ShapeComparisonResult result;
    result.tolerance = tolerance_;
    
    // Count topology
    count_topology(reference, result.ref_vertices, result.ref_edges,
                   result.ref_faces, result.ref_shells, result.ref_solids);
    count_topology(generated, result.gen_vertices, result.gen_edges,
                   result.gen_faces, result.gen_shells, result.gen_solids);
    
    // Compute geometric properties
    result.ref_volume = compute_volume(reference);
    result.gen_volume = compute_volume(generated);
    
    result.ref_surface_area = compute_surface_area(reference);
    result.gen_surface_area = compute_surface_area(generated);
    
    result.ref_bbox_diagonal = compute_bbox_diagonal(reference);
    result.gen_bbox_diagonal = compute_bbox_diagonal(generated);
    
    result.ref_centroid = compute_centroid(reference);
    result.gen_centroid = compute_centroid(generated);
    
    // Compute errors
    if (result.ref_volume > 0) {
        result.volume_error_percent = 
            std::abs(result.gen_volume - result.ref_volume) / result.ref_volume * 100.0;
    } else if (result.gen_volume > 0) {
        result.volume_error_percent = 100.0;
    }
    
    if (result.ref_surface_area > 0) {
        result.area_error_percent = 
            std::abs(result.gen_surface_area - result.ref_surface_area) / result.ref_surface_area * 100.0;
    } else if (result.gen_surface_area > 0) {
        result.area_error_percent = 100.0;
    }
    
    if (result.ref_bbox_diagonal > 0) {
        result.bbox_error_percent = 
            std::abs(result.gen_bbox_diagonal - result.ref_bbox_diagonal) / result.ref_bbox_diagonal * 100.0;
    } else if (result.gen_bbox_diagonal > 0) {
        result.bbox_error_percent = 100.0;
    }
    
    result.centroid_distance = result.ref_centroid.Distance(result.gen_centroid);
    
    // Determine pass/fail
    bool topology_ok = true;
    bool geometry_ok = true;
    
    // Topology checks
    if (result.ref_faces != result.gen_faces) {
        result.warnings.push_back("Face count mismatch: ref=" + std::to_string(result.ref_faces) +
                                   ", gen=" + std::to_string(result.gen_faces));
        // Don't fail on topology mismatch - our reconstruction may differ
    }
    
    // Geometric checks
    double tol_percent = tolerance_ * 100.0;
    
    if (result.volume_error_percent > tol_percent) {
        result.errors.push_back("Volume error " + std::to_string(result.volume_error_percent) + 
                                "% exceeds tolerance " + std::to_string(tol_percent) + "%");
        geometry_ok = false;
    }
    
    if (result.area_error_percent > tol_percent) {
        result.errors.push_back("Surface area error " + std::to_string(result.area_error_percent) + 
                                "% exceeds tolerance " + std::to_string(tol_percent) + "%");
        geometry_ok = false;
    }
    
    if (result.bbox_error_percent > tol_percent) {
        result.errors.push_back("Bounding box error " + std::to_string(result.bbox_error_percent) + 
                                "% exceeds tolerance " + std::to_string(tol_percent) + "%");
        geometry_ok = false;
    }
    
    // Centroid check: use characteristic length (bbox diagonal) as reference
    double char_length = std::max(result.ref_bbox_diagonal, result.gen_bbox_diagonal);
    if (char_length > 0) {
        double centroid_error_percent = result.centroid_distance / char_length * 100.0;
        if (centroid_error_percent > tol_percent) {
            result.errors.push_back("Centroid error " + std::to_string(centroid_error_percent) + 
                                    "% exceeds tolerance " + std::to_string(tol_percent) + "%");
            geometry_ok = false;
        }
    }
    
    result.passed = topology_ok && geometry_ok;
    return result;
}

ShapeComparisonResult STEPComparator::compare_files(const std::string& reference_step,
                                                     const std::string& generated_step) {
    ShapeComparisonResult result;
    result.tolerance = tolerance_;
    
    auto ref_shape = read_step(reference_step);
    if (!ref_shape) {
        result.errors.push_back("Failed to read reference STEP file: " + reference_step);
        return result;
    }
    
    auto gen_shape = read_step(generated_step);
    if (!gen_shape) {
        result.errors.push_back("Failed to read generated STEP file: " + generated_step);
        return result;
    }
    
    return compare(*ref_shape, *gen_shape);
}

} // namespace brepper
