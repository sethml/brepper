#pragma once

#include <string>

namespace brepper {

// Pipeline stages (for --stage option)
enum class PipelineStage {
    Load = 1,        // Stage 1: Load mesh and sample points
    Segment = 2,     // Stage 2: RANSAC surface segmentation
    Assign = 3,      // Stage 3: Assign triangles to surfaces
    Boundary = 4,    // Stage 4: Detect boundaries
    BRep = 5,        // Stage 5: Build B-Rep
    Export = 6       // Stage 6: Export STEP file
};

// Units for input STL files
enum class Units {
    Millimeters,  // mm (default)
    Centimeters,  // cm
    Meters,       // m
    Inches        // in
};

inline double units_to_mm(Units units) {
    switch (units) {
        case Units::Millimeters: return 1.0;
        case Units::Centimeters: return 10.0;
        case Units::Meters: return 1000.0;
        case Units::Inches: return 25.4;
    }
    return 1.0;
}

// Configuration parameters for all processing stages
struct Config {
    // Input/Output
    std::string input_file;
    std::string output_file;
    bool verbose = false;
    bool quiet = false;
    bool debug = false;
    int num_threads = 0;  // 0 = auto-detect
    Units stl_units = Units::Millimeters;  // Units of STL file; coordinates converted to mm internally

    // Mesh preprocessing 
    double max_point_distance_mm = 0.2;  // Maximum distance between sampled points (in mm)
    int min_samples_per_triangle = 1;
    size_t max_total_samples = 10000000;  // Cap to prevent runaway memory usage (10M points)
    
    // RANSAC segmentation
    double ransac_distance_threshold = 0.01;
    int ransac_max_iterations = 1000;
    double normal_distance_weight = 0.1;
    int min_inliers = 100;
    double min_inlier_ratio = 0.01;
    
    // Surface type toggles
    bool fit_planes = true;
    bool fit_cylinders = true;
    bool fit_spheres = true;
    bool fit_cones = true;
    bool fit_tori = false;
    
    // Surface merging thresholds
    double plane_merge_angle_degrees = 5.0;
    double plane_merge_distance = 0.01;
    double cylinder_radius_threshold = 0.01;
    
    // Clustering  
    double cluster_tolerance = 0.02;
    int min_cluster_size = 50;
    int max_cluster_size = 1000000;
    
    // NURBS fitting
    int nurbs_degree = 3;
    int nurbs_control_points = 10;
    double nurbs_fitting_tolerance = 0.001;
    
    // Triangle assignment
    double assignment_distance_threshold = 0.02;
    double assignment_angle_threshold_degrees = 15.0;
    
    // Curve fitting
    double curve_fitting_tolerance = 0.001;
    bool prefer_analytic_curves = true;
    
    // B-Rep construction
    double sewing_tolerance = 0.001;
    double healing_tolerance = 0.001;
    std::string step_schema = "AP214";
    
    // Debug outputs
    std::string save_point_cloud;
    std::string save_segmentation;  
    std::string save_boundaries;
    
    // Analysis options
    bool print_dimensions = false;  // Print mesh bounding box dimensions
    PipelineStage stop_after_stage = PipelineStage::Export;  // Run pipeline up to this stage
};

} // namespace brepper