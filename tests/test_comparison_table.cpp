// Comparison table generator with strict assertions
// Build: included in test executable
// Run: ./build/tests/brepper_tests "[comparison_table]" -s

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
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
    int expected_faces;
    int expected_solids;
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

TEST_CASE("Generate comparison table for all models", "[comparison_table]") {
    // Models with their expected face and solid counts from reference
    std::vector<ModelInfo> models = {
        {"cylinder_10x30_medium", "Cylinder", 3, 1},
        {"sphere_25_fine", "Sphere", 1, 1},
        {"cone_15x20_medium", "Cone", 2, 1},
        {"stepped_block_coarse", "Stepped Block", 16, 1},
        {"l_bracket_simple_medium", "L Bracket", 10, 1},
        {"plate_with_hole_100x50_coarse", "Plate+Hole", 7, 1},
        {"chamfered_cube_10_c1_medium", "Chamf. Cube", 26, 1},
        {"rounded_cube_10_r2_fine", "Round. Cube", 26, 1},
        {"pipe_elbow_10_fine", "Pipe Elbow", 5, 1},
        {"dome_hemisphere_20_fine", "Hemisphere", 2, 1},
    };
    
    std::string base = std::string(TEST_DATA_DIR) + "/onshape/";
    std::string temp_dir = std::filesystem::temp_directory_path().string();
    
    STEPComparator comparator;
    
    // Tolerance for geometric comparisons: 0.1%
    constexpr double TOLERANCE = 0.003;  // 0.3%
    
    fmt::print("\n");
    fmt::print("## STL to STEP Reconstruction Quality Comparison\n\n");
    const std::string row_fmt = "| {0:<13} | {1:>10} | {2:>10} | {3:>8} | {4:>11} | {5:>11} | {6:>8} | {7:>8} | {8:>8} |\n";
    fmt::print(row_fmt, "Model", "Ref Vol", "Gen Vol", "Vol Δ%", "Ref Area", "Gen Area", "Area Δ%", "Faces", "Solids");
    // Print a blank separator row instead of dashes, as fmt::print does not support per-field fill chars
    fmt::print(row_fmt, "", "", "", "", "", "", "", "", "");
    
    int total_failures = 0;
    
    for (const auto& model : models) {
        std::string stl_path = base + model.filename + ".stl";
        std::string ref_step = base + model.filename + ".step";
        std::string gen_step = temp_dir + "/brepper_cmp_" + model.filename + ".step";
        
        // Generate STEP
        std::string result = run_pipeline(stl_path, gen_step);
        
        if (result.empty() || !std::filesystem::exists(gen_step)) {
            fmt::print(row_fmt, model.display_name, "ERROR", "-", "-", "-", "-", "-", "-", "-");
            FAIL_CHECK("Pipeline failed for " << model.display_name);
            total_failures++;
            continue;
        }
        
        // Compare
        auto cmp = comparator.compare_files(ref_step, gen_step);
        
        // Format output
        std::string ref_vol = fmt::format("{:.1f}", cmp.ref_volume);
        std::string gen_vol, vol_delta;
        double vol_err_pct = 0.0;
        if (cmp.gen_volume > 0) {
            gen_vol = fmt::format("{:.1f}", cmp.gen_volume);
            vol_err_pct = 100.0 * (cmp.gen_volume - cmp.ref_volume) / cmp.ref_volume;
            vol_delta = fmt::format("{:.1f}", vol_err_pct);
        } else {
            gen_vol = "N/A";
            vol_delta = "N/A";
        }
        std::string ref_area = fmt::format("{:.1f}", cmp.ref_surface_area);
        std::string gen_area, area_delta;
        double area_err_pct = 0.0;
        if (cmp.gen_surface_area > 0) {
            gen_area = fmt::format("{:.1f}", cmp.gen_surface_area);
            area_err_pct = 100.0 * (cmp.gen_surface_area - cmp.ref_surface_area) / cmp.ref_surface_area;
            area_delta = fmt::format("{:.1f}", area_err_pct);
        } else {
            gen_area = "N/A";
            area_delta = "N/A";
        }
        std::string faces = fmt::format("{}/{}", cmp.ref_faces, cmp.gen_faces);
        std::string solids = fmt::format("{}/{}", cmp.ref_solids, cmp.gen_solids);
        fmt::print(row_fmt, model.display_name, ref_vol, gen_vol, vol_delta, ref_area, gen_area, area_delta, faces, solids);
        
        // Strict assertions for each model
        bool model_passed = true;
        
        // Must produce a valid solid
        if (cmp.gen_solids != model.expected_solids) {
            FAIL_CHECK(model.display_name << ": expected " << model.expected_solids 
                      << " solid(s), got " << cmp.gen_solids);
            model_passed = false;
        }
        
        // Face count must match reference
        if (cmp.gen_faces != model.expected_faces) {
            FAIL_CHECK(model.display_name << ": expected " << model.expected_faces 
                      << " faces, got " << cmp.gen_faces);
            model_passed = false;
        }
        
        // Volume must be within tolerance (skip if no valid solid)
        if (cmp.gen_volume > 0 && cmp.ref_volume > 0) {
            if (std::abs(vol_err_pct) > TOLERANCE * 100) {
                FAIL_CHECK(model.display_name << ": volume error " << vol_err_pct 
                          << "% exceeds tolerance of " << (TOLERANCE * 100) << "%");
                model_passed = false;
            }
        } else {
            FAIL_CHECK(model.display_name << ": no valid volume computed");
            model_passed = false;
        }
        
        // Surface area must be within tolerance
        if (cmp.gen_surface_area > 0 && cmp.ref_surface_area > 0) {
            if (std::abs(area_err_pct) > TOLERANCE * 100) {
                FAIL_CHECK(model.display_name << ": area error " << area_err_pct 
                          << "% exceeds tolerance of " << (TOLERANCE * 100) << "%");
                model_passed = false;
            }
        } else {
            FAIL_CHECK(model.display_name << ": no valid surface area computed");
            model_passed = false;
        }
        
        if (!model_passed) {
            total_failures++;
        }
    }
    
    std::cout << "\nNotes:\n";
    std::cout << "- Vol Δ% and Area Δ% show (generated - reference) / reference * 100\n";
    std::cout << "- Faces shows ref/gen count\n";
    std::cout << "- Solids shows ref/gen count (1/1 means valid closed solid)\n";
    std::cout << "- Tolerance: " << (TOLERANCE * 100) << "% for volume and area\n";
    std::cout << std::endl;
    
    // Require all models to pass
    REQUIRE(total_failures == 0);
}
