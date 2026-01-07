// Comparison table generator - run manually to see current state
// Build: included in test executable
// Run: ./build/tests/brepper_tests "[comparison_table]" -s

#include <catch2/catch_test_macros.hpp>
#include "comparison/step_comparator.hpp"
#include "brepper.hpp"
#include "common/config.hpp"
#include "test_config.hpp"
#include <filesystem>
#include <fmt/core.h>
#include <fmt/format.h>
#include <vector>

using namespace brepper;

struct ModelInfo {
    std::string filename;
    std::string display_name;
};

static std::string run_pipeline(const std::string& input_stl, 
                                const std::string& output_step) {
    Config config;
    config.verbose = false;
    config.quiet = true;
    config.input_file = input_stl;
    config.output_file = output_step;
    config.max_point_distance_mm = 0.5;
    config.min_inliers = 100;
    config.random_seed = 42;  // Use deterministic seed for reproducibility
    config.print_brep_diagnostics = false; // Hide B-Rep diagnostics in table output
    BrepperPipeline pipeline(config);
    if (pipeline.process()) {
        return output_step;
    }
    return "";
}

TEST_CASE("Generate comparison table for all models", "[comparison_table][manual]") {
    std::vector<ModelInfo> models = {
        {"cylinder_10x30_medium", "Cylinder"},
        {"sphere_25_fine", "Sphere"},
        {"cone_15x20_medium", "Cone"},
        {"stepped_block_coarse", "Stepped Block"},
        {"l_bracket_simple_medium", "L Bracket"},
        {"plate_with_hole_100x50_coarse", "Plate+Hole"},
        {"chamfered_cube_10_c1_medium", "Chamf. Cube"},
        {"rounded_cube_10_r2_fine", "Round. Cube"},
        {"pipe_elbow_10_fine", "Pipe Elbow"},
        {"dome_hemisphere_20_fine", "Hemisphere"},
    };
    
    std::string base = std::string(TEST_DATA_DIR) + "/onshape/";
    std::string temp_dir = std::filesystem::temp_directory_path().string();
    
    STEPComparator comparator;
    
    fmt::print("\n");
    fmt::print("## STL to STEP Reconstruction Quality Comparison\n\n");
    const std::string row_fmt = "| {0:<13} | {1:>10} | {2:>10} | {3:>8} | {4:>11} | {5:>11} | {6:>8} | {7:>8} | {8:>8} |\n";
    fmt::print(row_fmt, "Model", "Ref Vol", "Gen Vol", "Vol Δ%", "Ref Area", "Gen Area", "Area Δ%", "Faces", "Solids");
    // Print a blank separator row instead of dashes, as fmt::print does not support per-field fill chars
    fmt::print(row_fmt, "", "", "", "", "", "", "", "", "");
    
    for (const auto& model : models) {
        std::string stl_path = base + model.filename + ".stl";
        std::string ref_step = base + model.filename + ".step";
        std::string gen_step = temp_dir + "/brepper_cmp_" + model.filename + ".step";
        
        // Generate STEP
        std::string result = run_pipeline(stl_path, gen_step);
        
        if (result.empty() || !std::filesystem::exists(gen_step)) {
            fmt::print(row_fmt, model.display_name, "ERROR", "-", "-", "-", "-", "-", "-", "-");
            continue;
        }
        
        // Compare
        auto cmp = comparator.compare_files(ref_step, gen_step);
        
        // Format output
        std::string ref_vol = fmt::format("{:.1f}", cmp.ref_volume);
        std::string gen_vol, vol_delta;
        if (cmp.gen_volume > 0) {
            gen_vol = fmt::format("{:.1f}", cmp.gen_volume);
            double vol_err = 100.0 * (cmp.gen_volume - cmp.ref_volume) / cmp.ref_volume;
            vol_delta = fmt::format("{:.1f}", vol_err);
        } else {
            gen_vol = "N/A";
            vol_delta = "N/A";
        }
        std::string ref_area = fmt::format("{:.1f}", cmp.ref_surface_area);
        std::string gen_area, area_delta;
        if (cmp.gen_surface_area > 0) {
            gen_area = fmt::format("{:.1f}", cmp.gen_surface_area);
            double area_err = 100.0 * (cmp.gen_surface_area - cmp.ref_surface_area) / cmp.ref_surface_area;
            area_delta = fmt::format("{:.1f}", area_err);
        } else {
            gen_area = "N/A";
            area_delta = "N/A";
        }
        std::string faces = fmt::format("{}/{}", cmp.ref_faces, cmp.gen_faces);
        std::string solids = fmt::format("{}/{}", cmp.ref_solids, cmp.gen_solids);
        fmt::print(row_fmt, model.display_name, ref_vol, gen_vol, vol_delta, ref_area, gen_area, area_delta, faces, solids);
    }
    
    std::cout << "\nNotes:\n";
    std::cout << "- Vol Δ% and Area Δ% show (generated - reference) / reference * 100\n";
    std::cout << "- Faces shows ref/gen count\n";
    std::cout << "- Solids shows ref/gen count (1/1 means valid closed solid)\n";
    std::cout << "- N/A means no valid solid was generated (volume = 0)\n";
    std::cout << std::endl;
    
    // Always pass - this is just for reporting
    CHECK(true);
}
