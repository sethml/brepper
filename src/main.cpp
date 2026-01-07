#include <CLI/CLI.hpp>
#include <iostream>
#include <filesystem>

#include "brepper.hpp"
#include "common/config.hpp"
#include "common/logging.hpp"

using namespace brepper;

bool validate_files(const Config& config) {
    // Check input file exists
    if (!std::filesystem::exists(config.input_file)) {
        LOG_ERROR("Input file does not exist: ", config.input_file);
        return false;
    }
    
    // Check input file has .stl extension
    std::filesystem::path input_path(config.input_file);
    if (input_path.extension() != ".stl") {
        LOG_ERROR("Input file must have .stl extension");
        return false;
    }
    
    // Only require output file if running to export stage
    if (config.stop_after_stage != PipelineStage::Export) {
        return true;
    }
    
    // Require output file for export stage
    if (config.output_file.empty()) {
        LOG_ERROR("Output file required when running to export stage (use -o or --stage to stop earlier)");
        return false;
    }
    
    // Check output file has .step extension
    std::filesystem::path output_path(config.output_file);
    if (output_path.extension() != ".step" && output_path.extension() != ".stp") {
        LOG_ERROR("Output file must have .step or .stp extension");
        return false;
    }
    
    // Create output directory if needed
    std::filesystem::path output_dir = output_path.parent_path();
    if (!output_dir.empty() && !std::filesystem::exists(output_dir)) {
        try {
            std::filesystem::create_directories(output_dir);
        } catch (const std::exception& e) {
            LOG_ERROR("Failed to create output directory: ", e.what());
            return false;
        }
    }
    
    return true;
}

int main(int argc, char* argv[]) {
    CLI::App app{"brepper - Convert STL mesh to STEP with fitted surfaces", "brepper"};
    
    Config config;
    
    // Required arguments
    app.add_option("input", config.input_file, "Input STL file (binary or ASCII)")
        ->required()
        ->check(CLI::ExistingFile);
        
    app.add_option("-o,--output", config.output_file, "Output STEP file");
    
    // General options
    app.add_flag("-v,--verbose", config.verbose, "Enable verbose output");
    app.add_flag("-q,--quiet", config.quiet, "Suppress non-error output");
    app.add_flag("--debug", config.debug, "Enable debug output and intermediate files");
    app.add_option("--threads", config.num_threads, "Number of threads (default: auto)");
    app.add_flag("--dimensions", config.print_dimensions, "Print mesh bounding box dimensions");
    
    // Pipeline stage control
    int stage_num = 6;
    app.add_option("--stage", stage_num,
                  "Stop after stage: 1=load, 2=segment, 3=assign, 4=boundary, 5=brep, 6=export (default: 6)")
        ->check(CLI::Range(1, 6));
    
    // Units option
    std::string units_str = "mm";
    app.add_option("--units", units_str, 
                  "Units of STL file coordinates (converted to mm internally): mm|cm|m|in (default: mm)")
        ->check(CLI::IsMember({"mm", "cm", "m", "in"}));
    
    // Mesh preprocessing
    auto* mesh_group = app.add_option_group("Mesh Preprocessing");
    mesh_group->add_option("--max-point-distance-mm", config.max_point_distance_mm, 
                          "Maximum distance between sampled points in mm (default: 0.2)");
    mesh_group->add_option("--min-samples", config.min_samples_per_triangle,
                          "Min samples per triangle (default: 1)");
    mesh_group->add_option("--seed", config.random_seed,
                          "Random seed for mesh sampling (-1 = non-deterministic, default: -1)");
    
    // RANSAC segmentation  
    auto* ransac_group = app.add_option_group("RANSAC Segmentation");
    ransac_group->add_option("--ransac-distance", config.ransac_distance_threshold,
                            "Distance threshold (default: 0.01)");
    ransac_group->add_option("--ransac-iterations", config.ransac_max_iterations,
                            "Max iterations (default: 1000)");
    ransac_group->add_option("--normal-weight", config.normal_distance_weight,
                            "Normal weight 0-1 (default: 0.1)");
    ransac_group->add_option("--min-inliers", config.min_inliers,
                            "Min points per surface (default: 100)");
    ransac_group->add_option("--min-inlier-ratio", config.min_inlier_ratio,
                            "Min ratio of cloud (default: 0.01)");
    
    // Surface types
    auto* surface_group = app.add_option_group("Surface Types");
    surface_group->add_flag("--no-planes,!--fit-planes", config.fit_planes,
                           "Disable plane fitting");
    surface_group->add_flag("--no-cylinders,!--fit-cylinders", config.fit_cylinders,
                           "Disable cylinder fitting");
    surface_group->add_flag("--no-spheres,!--fit-spheres", config.fit_spheres,
                           "Disable sphere fitting");
    surface_group->add_flag("--no-cones,!--fit-cones", config.fit_cones,
                           "Disable cone fitting");
    surface_group->add_flag("--fit-tori", config.fit_tori,
                           "Enable torus fitting");
    
    // Clustering
    auto* cluster_group = app.add_option_group("Clustering");
    cluster_group->add_option("--plane-merge-angle", config.plane_merge_angle_degrees,
                             "Plane merge angle threshold (deg, default: 5.0)");
    cluster_group->add_option("--plane-merge-dist", config.plane_merge_distance,
                             "Plane merge distance (default: 0.01)");
    cluster_group->add_option("--cluster-tolerance", config.cluster_tolerance,
                             "Euclidean cluster tolerance (default: 0.02)");
    cluster_group->add_option("--min-cluster-size", config.min_cluster_size,
                             "Minimum cluster size (default: 50)");
    
    // NURBS fitting
    auto* nurbs_group = app.add_option_group("NURBS Fitting");
    nurbs_group->add_option("--nurbs-degree", config.nurbs_degree,
                           "B-spline degree (default: 3)");
    nurbs_group->add_option("--nurbs-refinement", config.nurbs_control_points,
                           "Control point density (default: 10)");
    nurbs_group->add_option("--nurbs-tolerance", config.nurbs_fitting_tolerance,
                           "Fitting tolerance (default: 0.001)");
    
    // Triangle assignment
    auto* assign_group = app.add_option_group("Triangle Assignment");
    assign_group->add_option("--assign-distance", config.assignment_distance_threshold,
                            "Max assignment distance (default: 0.02)");
    assign_group->add_option("--assign-angle", config.assignment_angle_threshold_degrees,
                            "Max normal deviation (deg, default: 15.0)");
    
    // Curve fitting  
    auto* curve_group = app.add_option_group("Curve Fitting");
    curve_group->add_option("--curve-tolerance", config.curve_fitting_tolerance,
                           "Curve fitting tolerance (default: 0.001)");
    curve_group->add_flag("--prefer-analytic,!--prefer-splines", config.prefer_analytic_curves,
                         "Prefer analytic curves (default: true)");
    
    // B-Rep construction
    auto* brep_group = app.add_option_group("B-Rep Construction");
    brep_group->add_option("--sewing-tolerance", config.sewing_tolerance,
                          "Face sewing tolerance (default: 0.001)");
    brep_group->add_option("--healing-tolerance", config.healing_tolerance,
                          "Shape healing tolerance (default: 0.001)");
    brep_group->add_option("--step-schema", config.step_schema,
                          "STEP schema: AP203|AP214|AP242 (default: AP214)")
        ->check(CLI::IsMember({"AP203", "AP214", "AP242"}));
    
    // Debug output
    auto* debug_group = app.add_option_group("Debug Output");
    debug_group->add_option("--save-point-cloud", config.save_point_cloud,
                           "Save sampled point cloud (PCD/PLY)");
    debug_group->add_option("--save-segmentation", config.save_segmentation,
                           "Save segmented mesh (PLY with colors)");
    debug_group->add_option("--save-boundaries", config.save_boundaries,
                           "Save boundary curves (PLY)");
    
    // Parse command line
    CLI11_PARSE(app, argc, argv);
    
    // Convert units string to enum
    if (units_str == "mm") config.stl_units = Units::Millimeters;
    else if (units_str == "cm") config.stl_units = Units::Centimeters;
    else if (units_str == "m") config.stl_units = Units::Meters;
    else if (units_str == "in") config.stl_units = Units::Inches;
    
    // Convert stage number to enum
    config.stop_after_stage = static_cast<PipelineStage>(stage_num);
    
    // Validate inputs
    if (!validate_files(config)) {
        return 1;
    }
    
    // Run the pipeline
    try {
        BrepperPipeline pipeline(config);
        if (!pipeline.process()) {
            LOG_ERROR("Pipeline failed");
            return 1;
        }
    } catch (const std::exception& e) {
        LOG_ERROR("Unexpected error: ", e.what());
        return 1;
    }
    
    return 0;
}