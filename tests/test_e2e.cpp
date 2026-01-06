#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "brepper.hpp"
#include "common/config.hpp"
#include "test_config.hpp"
#include <filesystem>
#include <cstdlib>

using namespace brepper;

// ============================================================================
// End-to-End Pipeline Tests
// ============================================================================

TEST_CASE("Full pipeline runs on simple primitives", "[e2e][pipeline]") {
    Config config;
    config.verbose = false;
    config.quiet = true;
    
    // Use temp directory for output
    std::string temp_dir = std::filesystem::temp_directory_path().string();
    
    SECTION("Manual cube") {
        config.input_file = std::string(TEST_DATA_DIR) + "/manual/cube.stl";
        config.output_file = temp_dir + "/test_cube_output.step";
        
        BrepperPipeline pipeline(config);
        
        // Pipeline should run without crashing (stages 2-6 are stubs)
        bool result = pipeline.process();
        
        // Currently returns true even though later stages are not implemented
        CHECK(result == true);
    }
    
    SECTION("Onshape cylinder") {
        config.input_file = std::string(TEST_DATA_DIR) + "/onshape/cylinder_10x30_medium.stl";
        config.output_file = temp_dir + "/test_cylinder_output.step";
        
        BrepperPipeline pipeline(config);
        bool result = pipeline.process();
        CHECK(result == true);
    }
    
    SECTION("Onshape sphere") {
        config.input_file = std::string(TEST_DATA_DIR) + "/onshape/sphere_25_fine.stl";
        config.output_file = temp_dir + "/test_sphere_output.step";
        
        BrepperPipeline pipeline(config);
        bool result = pipeline.process();
        CHECK(result == true);
    }
}

TEST_CASE("Pipeline handles various mesh sizes", "[e2e][pipeline][performance]") {
    Config config;
    config.verbose = false;
    config.quiet = true;
    config.max_point_distance_mm = 0.5;  // Reasonable density for testing
    
    std::string temp_dir = std::filesystem::temp_directory_path().string();
    
    SECTION("Coarse mesh (stepped_block)") {
        config.input_file = std::string(TEST_DATA_DIR) + "/onshape/stepped_block_coarse.stl";
        config.output_file = temp_dir + "/test_stepped_output.step";
        
        BrepperPipeline pipeline(config);
        CHECK(pipeline.process() == true);
    }
    
    SECTION("Fine mesh (rounded_cube with fillets)") {
        config.input_file = std::string(TEST_DATA_DIR) + "/onshape/rounded_cube_10_r2_fine.stl";
        config.output_file = temp_dir + "/test_rounded_cube_output.step";
        config.max_point_distance_mm = 1.0;  // Coarser sampling for large mesh
        
        BrepperPipeline pipeline(config);
        CHECK(pipeline.process() == true);
    }
}

TEST_CASE("Pipeline with different unit settings", "[e2e][pipeline][units]") {
    Config config;
    config.verbose = false;
    config.quiet = true;
    
    std::string temp_dir = std::filesystem::temp_directory_path().string();
    config.input_file = std::string(TEST_DATA_DIR) + "/onshape/cylinder_10x30_medium.stl";
    config.output_file = temp_dir + "/test_units_output.step";
    
    SECTION("Default units (mm)") {
        config.stl_units = Units::Millimeters;
        config.max_point_distance_mm = 1.0;
        BrepperPipeline pipeline(config);
        CHECK(pipeline.process() == true);
    }
    
    SECTION("Inches") {
        config.stl_units = Units::Inches;
        // Inches -> mm scales by 25.4, so use proportionally larger spacing
        config.max_point_distance_mm = 25.0;
        BrepperPipeline pipeline(config);
        CHECK(pipeline.process() == true);
    }
    
    SECTION("Meters") {
        config.stl_units = Units::Meters;
        // Meters -> mm scales by 1000, so use much larger spacing
        config.max_point_distance_mm = 1000.0;
        BrepperPipeline pipeline(config);
        CHECK(pipeline.process() == true);
    }
}

TEST_CASE("Pipeline with debug output options", "[e2e][pipeline][debug]") {
    Config config;
    config.verbose = false;
    config.quiet = true;
    config.debug = true;
    
    std::string temp_dir = std::filesystem::temp_directory_path().string();
    config.input_file = std::string(TEST_DATA_DIR) + "/onshape/cone_15x20_medium.stl";
    config.output_file = temp_dir + "/test_debug_output.step";
    
    // Enable debug outputs (if implemented)
    config.save_point_cloud = temp_dir + "/test_pointcloud.pcd";
    
    BrepperPipeline pipeline(config);
    CHECK(pipeline.process() == true);
}

// ============================================================================
// CLI End-to-End Tests (run actual executable)
// ============================================================================

TEST_CASE("CLI executable runs successfully", "[e2e][cli]") {
    std::string brepper_exe = std::string(TEST_DATA_DIR) + "/../build/brepper";
    std::string temp_dir = std::filesystem::temp_directory_path().string();
    
    // Check if executable exists
    if (!std::filesystem::exists(brepper_exe)) {
        SKIP("brepper executable not found at: " << brepper_exe);
    }
    
    SECTION("Help flag") {
        std::string cmd = brepper_exe + " --help > /dev/null 2>&1";
        int result = std::system(cmd.c_str());
        CHECK(result == 0);
    }
    
    SECTION("Process simple cube") {
        std::string input = std::string(TEST_DATA_DIR) + "/manual/cube.stl";
        std::string output = temp_dir + "/cli_test_cube.step";
        std::string cmd = brepper_exe + " \"" + input + "\" -o \"" + output + "\" -q 2>&1";
        
        int result = std::system(cmd.c_str());
        CHECK(result == 0);
    }
    
    SECTION("Process with verbose flag") {
        std::string input = std::string(TEST_DATA_DIR) + "/onshape/cylinder_10x30_medium.stl";
        std::string output = temp_dir + "/cli_test_cylinder.step";
        std::string cmd = brepper_exe + " \"" + input + "\" -o \"" + output + "\" -v > /dev/null 2>&1";
        
        int result = std::system(cmd.c_str());
        CHECK(result == 0);
    }
    
    SECTION("Process with custom sampling") {
        std::string input = std::string(TEST_DATA_DIR) + "/onshape/sphere_25_fine.stl";
        std::string output = temp_dir + "/cli_test_sphere.step";
        std::string cmd = brepper_exe + " \"" + input + "\" -o \"" + output + "\" "
                         "--max-point-distance-mm 2.0 -q 2>&1";
        
        int result = std::system(cmd.c_str());
        CHECK(result == 0);
    }
    
    SECTION("Process with units flag") {
        std::string input = std::string(TEST_DATA_DIR) + "/onshape/cone_15x20_medium.stl";
        std::string output = temp_dir + "/cli_test_cone.step";
        std::string cmd = brepper_exe + " \"" + input + "\" -o \"" + output + "\" "
                         "--units mm -q 2>&1";
        
        int result = std::system(cmd.c_str());
        CHECK(result == 0);
    }
}

// ============================================================================
// All Models Pipeline Smoke Test
// ============================================================================

TEST_CASE("Pipeline processes all Onshape models", "[e2e][pipeline][smoke]") {
    Config config;
    config.verbose = false;
    config.quiet = true;
    config.max_point_distance_mm = 2.0;  // Coarse for speed
    
    std::string temp_dir = std::filesystem::temp_directory_path().string();
    
    std::vector<std::string> models = {
        "chamfered_cube_10_c1_medium.stl",
        "cone_15x20_medium.stl",
        "cylinder_10x30_medium.stl",
        "dome_hemisphere_20_fine.stl",
        "l_bracket_simple_medium.stl",
        "pipe_elbow_10_fine.stl",
        "plate_with_hole_100x50_coarse.stl",
        "rounded_cube_10_r2_fine.stl",
        "sphere_25_fine.stl",
        "stepped_block_coarse.stl"
    };
    
    for (const auto& model : models) {
        SECTION(model) {
            config.input_file = std::string(TEST_DATA_DIR) + "/onshape/" + model;
            config.output_file = temp_dir + "/test_" + model + ".step";
            
            INFO("Processing: " << model);
            REQUIRE(std::filesystem::exists(config.input_file));
            
            BrepperPipeline pipeline(config);
            CHECK(pipeline.process() == true);
        }
    }
}
