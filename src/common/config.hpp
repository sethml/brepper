#pragma once

#include <string>

namespace brepper {

// Configuration parameters for all processing stages
struct Config {
    // Input/Output
    std::string input_file;
    std::string output_file;
    bool verbose = false;
    bool quiet = false;
    bool debug = false;
    int num_threads = 0;  // 0 = auto-detect

    // Mesh preprocessing 
    double sample_density = 0.0;  // 0 = auto-compute
    int min_samples_per_triangle = 1;
    
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
};

} // namespace brepper