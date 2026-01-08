#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "comparison/step_comparator.hpp"
#include "brepper.hpp"
#include "common/config.hpp"
#include "test_config.hpp"
#include <filesystem>
#include <string>
#include <vector>
#include <fmt/core.h>
#include <fmt/format.h>
#include <unistd.h>

using namespace brepper;

// ============================================================================
// Model Configuration
// ============================================================================

enum class ModelSource {
    Onshape,
    CodeCAD
};

struct ModelInfo {
    std::string id;             // Filename without extension
    std::string display_name;   // Human readable name
    ModelSource source;         // Where the file comes from
    // Expected topology (optional, -1 if unknown)
    int expected_faces = -1;
    int expected_solids = -1;
};

// List of all models to test
static const std::vector<ModelInfo> ALL_MODELS = {
    // Onshape models
    {"cylinder_10x30_medium", "Cylinder", ModelSource::Onshape, 3, 1},
    {"sphere_25_fine", "Sphere", ModelSource::Onshape, 1, 1},
    {"cone_15x20_medium", "Cone", ModelSource::Onshape, 2, 1},
    {"stepped_block_coarse", "Stepped Block", ModelSource::Onshape, 16, 1},
    {"l_bracket_simple_medium", "L Bracket", ModelSource::Onshape, 10, 1},
    {"plate_with_hole_100x50_coarse", "Plate+Hole", ModelSource::Onshape, 7, 1},
    {"chamfered_cube_10_c1_medium", "Chamf. Cube", ModelSource::Onshape, 26, 1},
    {"rounded_cube_10_r2_fine", "Round. Cube", ModelSource::Onshape, 26, 1},
    {"pipe_elbow_10_fine", "Pipe Elbow", ModelSource::Onshape, 5, 1},
    {"dome_hemisphere_20_fine", "Hemisphere", ModelSource::Onshape, 2, 1},
    
    // CodeCAD models
    {"cube", "CCAD Cube", ModelSource::CodeCAD, 6, 1}
};

// ============================================================================
// Helpers
// ============================================================================

static std::pair<std::string, std::string> get_model_paths(const ModelInfo& model) {
    std::string base_dir;
    if (model.source == ModelSource::Onshape) {
        base_dir = std::string(TEST_DATA_DIR) + "/onshape/";
    } else {
        base_dir = std::string(TEST_DATA_DIR) + "/ccad/generated/";
    }
    return {base_dir + model.id + ".stl", base_dir + model.id + ".step"};
}

static std::string run_pipeline(const std::string& input_stl, 
                                const std::string& output_step,
                                bool quiet = true) {
    Config config;
    config.verbose = !quiet;
    config.quiet = quiet;
    config.input_file = input_stl;
    config.output_file = output_step;
    // Default params suitable for typical test models (mm scale)
    config.max_point_distance_mm = 0.5;
    config.min_inliers = 100;
    config.random_seed = 42; 
    config.print_brep_diagnostics = false;

    BrepperPipeline pipeline(config);
    if (pipeline.process()) {
        return output_step;
    }
    return "";
}

// ============================================================================
// Unit Tests (STEPComparator specific)
// ============================================================================

TEST_CASE("STEPComparator functionality", "[step_comparison][unit]") {
    STEPComparator comparator;
    
    SECTION("Read valid STEP file") {
        std::string step_file = std::string(TEST_DATA_DIR) + "/onshape/cylinder_10x30_medium.step";
        auto shape = comparator.read_step(step_file);
        REQUIRE(shape.has_value());
    }
    
    SECTION("Topology counting") {
        // Use sphere as simple case
        std::string step_file = std::string(TEST_DATA_DIR) + "/onshape/sphere_25_fine.step";
        auto shape = comparator.read_step(step_file);
        REQUIRE(shape.has_value());
        
        int v, e, f, sh, so;
        STEPComparator::count_topology(*shape, v, e, f, sh, so);
        CHECK(f == 1);
        CHECK(so == 1);
    }
}

// ============================================================================
// Individual Model Integration Tests
// ============================================================================

static void test_single_model(const ModelInfo& model) {
    auto [input_stl, ref_step] = get_model_paths(model);
    
    // Skip if files don't exist (e.g. CCAD generation failed)
    if (!std::filesystem::exists(input_stl) || !std::filesystem::exists(ref_step)) {
        // If CCAD, maybe we should warn? But for now, skip or fail.
        // Failing is better to alert the user.
        FAIL("Missing test files for " + model.id + ": " + input_stl);
    }
    
    std::string temp_step = std::filesystem::temp_directory_path().string() 
                          + "/brepper_test_" + std::to_string(getpid()) + "_" + model.id + ".step";
    
    // Run reconstruction
    std::string result_path = run_pipeline(input_stl, temp_step);
    REQUIRE(!result_path.empty());
    REQUIRE(std::filesystem::exists(temp_step));
    
    // Compare
    STEPComparator comparator;
    // 5% tolerance for initial pass
    comparator.set_tolerance(0.05);
    auto result = comparator.compare_files(ref_step, temp_step); 
    
    INFO("Model: " << model.display_name);
    INFO("Volume Error: " << result.volume_error_percent << "%");
    INFO("Area Error: " << result.area_error_percent << "%");
    
    CHECK(result.volume_error_percent < 5.0);
    // Area error might be higher due to triangulation differences, but keep check loose
    CHECK(result.area_error_percent < 10.0);
    
    if (std::filesystem::exists(temp_step)) {
        std::filesystem::remove(temp_step);
    }
}

// Define individual test cases for parallel execution availability
// We macro these to avoid boilerplate code
#define TEST_MODEL(index) \
    TEST_CASE("Reconstruct_" + ALL_MODELS[index].id, "[step_comparison][integration]") { \
        test_single_model(ALL_MODELS[index]); \
    }

// Manually unrolled for now as macros can't iterate. 
// If list grows large, we might move to a generator or just loop in one test case.
// For parallel execution via CTest, separate TEST_CASEs are required.
TEST_MODEL(0)   // cylinder
TEST_MODEL(1)   // sphere
TEST_MODEL(2)   // cone
TEST_MODEL(3)   // stepped_block
TEST_MODEL(4)   // l_bracket
TEST_MODEL(5)   // plate_with_hole
TEST_MODEL(6)   // chamfered_cube
TEST_MODEL(7)   // rounded_cube
TEST_MODEL(8)   // pipe_elbow
TEST_MODEL(9)   // hemisphere
TEST_MODEL(10)  // ccad cube

// ============================================================================
// Comparison Table Generator
// ============================================================================

TEST_CASE("Generate comparison table for all models", "[comparison_table]") {
    std::string temp_dir = std::filesystem::temp_directory_path().string();
    STEPComparator comparator;
    
    fmt::print("\n");
    fmt::print("## reconstruction Quality Comparison Table\n\n");
    const std::string row_fmt = "| {0:<15} | {1:>10} | {2:>10} | {3:>8} | {4:>11} | {5:>11} | {6:>8} | {7:>8} | {8:>8} |\n";
    fmt::print(row_fmt, "Model", "Ref Vol", "Gen Vol", "Vol Δ%", "Ref Area", "Gen Area", "Area Δ%", "Faces", "Solids");
    fmt::print(row_fmt, "-","-","-","-","-","-","-","-","-"); // Header separator line (markdown)

    for (const auto& model : ALL_MODELS) {
        auto [input_stl, ref_step] = get_model_paths(model);
        std::string gen_step = temp_dir + "/brepper_cmp_" + model.id + ".step";

        if (!std::filesystem::exists(input_stl)) {
            fmt::print(row_fmt, model.display_name, "MISSING", "-", "-", "-", "-", "-", "-", "-");
            continue;
        }

        std::string result = run_pipeline(input_stl, gen_step, true);
        if (result.empty()) {
             fmt::print(row_fmt, model.display_name, "FAILED", "-", "-", "-", "-", "-", "-", "-");
             continue;
        }

        auto cmp = comparator.compare_files(ref_step, gen_step);
        
        // Format fields
        std::string ref_vol = fmt::format("{:.1f}", cmp.ref_volume);
        std::string gen_vol = fmt::format("{:.1f}", cmp.gen_volume);
        
        double vol_err_pct = (cmp.ref_volume > 1e-6) ? 
            100.0 * (cmp.gen_volume - cmp.ref_volume) / cmp.ref_volume : 0.0;
            
        std::string vol_delta = fmt::format("{:+.1f}", vol_err_pct);
        
        std::string ref_area = fmt::format("{:.1f}", cmp.ref_surface_area);
        std::string gen_area = fmt::format("{:.1f}", cmp.gen_surface_area);

        double area_err_pct = (cmp.ref_surface_area > 1e-6) ?
            100.0 * (cmp.gen_surface_area - cmp.ref_surface_area) / cmp.ref_surface_area : 0.0;
            
        std::string area_delta = fmt::format("{:+.1f}", area_err_pct);

        fmt::print(row_fmt, 
            model.display_name, 
            ref_vol, gen_vol, vol_delta,
            ref_area, gen_area, area_delta,
            "?", "?" // TODO: Count faces if needed
        );
        
        // Cleanup
        std::filesystem::remove(gen_step);
    }
    fmt::print("\n");
}
