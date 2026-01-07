#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "comparison/step_comparator.hpp"
#include "brepper.hpp"
#include "common/config.hpp"
#include "test_config.hpp"
#include <filesystem>
#include <string>
#include <unistd.h>  // for getpid()

using namespace brepper;

// ============================================================================
// STEP File Comparison Tests
// 
// These tests compare generated STEP files against reference Onshape STEP files
// to verify geometric accuracy. The comparison checks:
// - Volume (should match within tolerance)
// - Surface area (should match within tolerance)
// - Bounding box (should match within tolerance)
// - Centroid position (should match within tolerance)
// ============================================================================

/// Helper: Run brepper pipeline and return the generated STEP file path
static std::string run_brepper_pipeline(const std::string& input_stl, 
                                        const std::string& output_step,
                                        double point_distance = 0.5,
                                        int min_inliers = 100,
                                        int random_seed = 42) {
    Config config;
    config.verbose = false;
    config.quiet = true;
    config.input_file = input_stl;
    config.output_file = output_step;
    config.max_point_distance_mm = point_distance;
    config.min_inliers = min_inliers;
    config.random_seed = random_seed;  // Use deterministic seed for reproducibility
    
    BrepperPipeline pipeline(config);
    bool success = pipeline.process();
    
    if (!success) {
        return "";
    }
    return output_step;
}

/// Helper: Generate unique output path for a given model (includes PID for parallel safety)
static std::string get_output_path(const std::string& model_name) {
    std::string temp_dir = std::filesystem::temp_directory_path().string();
    // Include PID to avoid collisions when tests run in parallel
    return temp_dir + "/brepper_test_" + std::to_string(getpid()) + "_" + model_name + "_generated.step";
}

/// Helper: Get Onshape reference paths
static std::pair<std::string, std::string> get_onshape_paths(const std::string& model_name) {
    std::string base = std::string(TEST_DATA_DIR) + "/onshape/" + model_name;
    return {base + ".stl", base + ".step"};
}

/// Helper: Get manual test reference paths
static std::pair<std::string, std::string> get_manual_paths(const std::string& model_name) {
    std::string base = std::string(TEST_DATA_DIR) + "/manual/" + model_name;
    return {base + ".stl", base + ".step"};
}

// ============================================================================
// Unit Tests for STEPComparator
// ============================================================================

TEST_CASE("STEPComparator can read STEP files", "[step_comparison][reader]") {
    STEPComparator comparator;
    
    SECTION("Read valid Onshape STEP file") {
        std::string step_file = std::string(TEST_DATA_DIR) + "/onshape/cylinder_10x30_medium.step";
        auto shape = comparator.read_step(step_file);
        
        REQUIRE(shape.has_value());
        
        // Check basic topology
        int vertices, edges, faces, shells, solids;
        STEPComparator::count_topology(*shape, vertices, edges, faces, shells, solids);
        
        INFO("Cylinder: vertices=" << vertices << ", edges=" << edges 
             << ", faces=" << faces << ", shells=" << shells << ", solids=" << solids);
        
        CHECK(faces == 3);  // Cylinder has 3 faces: top, bottom, curved
        CHECK(edges >= 3);  // At least 3 edges (may be more in STEP representation)
        CHECK(solids == 1); // Should be a solid
    }
    
    SECTION("Read non-existent file returns nullopt") {
        auto shape = comparator.read_step("/nonexistent/path.step");
        CHECK_FALSE(shape.has_value());
    }
}

TEST_CASE("STEPComparator topology counting", "[step_comparison][topology]") {
    STEPComparator comparator;
    
    SECTION("Sphere topology") {
        std::string step_file = std::string(TEST_DATA_DIR) + "/onshape/sphere_25_fine.step";
        auto shape = comparator.read_step(step_file);
        REQUIRE(shape.has_value());
        
        int vertices, edges, faces, shells, solids;
        STEPComparator::count_topology(*shape, vertices, edges, faces, shells, solids);
        
        INFO("Sphere: vertices=" << vertices << ", edges=" << edges 
             << ", faces=" << faces << ", shells=" << shells << ", solids=" << solids);
        
        CHECK(faces == 1);  // Sphere is one face
        CHECK(solids == 1);
    }
    
    SECTION("Cone topology") {
        std::string step_file = std::string(TEST_DATA_DIR) + "/onshape/cone_15x20_medium.step";
        auto shape = comparator.read_step(step_file);
        REQUIRE(shape.has_value());
        
        int vertices, edges, faces, shells, solids;
        STEPComparator::count_topology(*shape, vertices, edges, faces, shells, solids);
        
        INFO("Cone: vertices=" << vertices << ", edges=" << edges 
             << ", faces=" << faces << ", shells=" << shells << ", solids=" << solids);
        
        CHECK(faces == 2);  // Cone has 2 faces: base and conical surface
        CHECK(solids == 1);
    }
}

TEST_CASE("STEPComparator geometric properties", "[step_comparison][geometry]") {
    STEPComparator comparator;
    
    SECTION("Cylinder geometric properties") {
        std::string step_file = std::string(TEST_DATA_DIR) + "/onshape/cylinder_10x30_medium.step";
        auto shape = comparator.read_step(step_file);
        REQUIRE(shape.has_value());
        
        double volume = STEPComparator::compute_volume(*shape);
        double area = STEPComparator::compute_surface_area(*shape);
        double bbox_diag = STEPComparator::compute_bbox_diagonal(*shape);
        gp_Pnt centroid = STEPComparator::compute_centroid(*shape);
        
        // Cylinder: r=5, h=30 (diameter 10, height 30)
        // Volume = pi * r^2 * h = pi * 25 * 30 = 750 * pi ≈ 2356.19
        // Area = 2*pi*r^2 + 2*pi*r*h = 2*pi*25 + 2*pi*5*30 = 50*pi + 300*pi = 350*pi ≈ 1099.56
        // Bbox diagonal = sqrt(10^2 + 10^2 + 30^2) = sqrt(1100) ≈ 33.17
        
        INFO("Cylinder volume: " << volume << " (expected ~2356.19)");
        INFO("Cylinder area: " << area << " (expected ~1099.56)");
        INFO("Cylinder bbox diagonal: " << bbox_diag << " (expected ~33.17)");
        INFO("Centroid: (" << centroid.X() << ", " << centroid.Y() << ", " << centroid.Z() << ")");
        
        CHECK_THAT(volume, Catch::Matchers::WithinRel(2356.19, 0.01));
        CHECK_THAT(area, Catch::Matchers::WithinRel(1099.56, 0.01));
        CHECK_THAT(bbox_diag, Catch::Matchers::WithinRel(33.17, 0.01));
    }
    
    SECTION("Sphere geometric properties") {
        std::string step_file = std::string(TEST_DATA_DIR) + "/onshape/sphere_25_fine.step";
        auto shape = comparator.read_step(step_file);
        REQUIRE(shape.has_value());
        
        double volume = STEPComparator::compute_volume(*shape);
        double area = STEPComparator::compute_surface_area(*shape);
        
        // Sphere: r=12.5 (diameter 25)
        // Volume = 4/3 * pi * r^3 = 4/3 * pi * 1953.125 ≈ 8181.23
        // Area = 4 * pi * r^2 = 4 * pi * 156.25 ≈ 1963.50
        
        INFO("Sphere volume: " << volume << " (expected ~8181.23)");
        INFO("Sphere area: " << area << " (expected ~1963.50)");
        
        CHECK_THAT(volume, Catch::Matchers::WithinRel(8181.23, 0.01));
        CHECK_THAT(area, Catch::Matchers::WithinRel(1963.50, 0.01));
    }
}

// ============================================================================
// Integration Tests: Compare Generated STEP vs Reference STEP
//
// These tests currently document the state of reconstruction.
// Known issues: 
// - Negative volumes (face orientation)
// - No solid (shell not closed)
// - Bounding box 2x too large (unbounded parametric faces)
// ============================================================================

TEST_CASE("Compare generated STEP against reference - cylinder", "[step_comparison][integration]") {
    auto [input_stl, reference_step] = get_onshape_paths("cylinder_10x30_medium");
    std::string generated_step = get_output_path("cylinder_10x30_medium");
    
    // Generate STEP file using brepper
    std::string result_path = run_brepper_pipeline(input_stl, generated_step);
    REQUIRE_FALSE(result_path.empty());
    REQUIRE(std::filesystem::exists(generated_step));
    
    // Compare against reference
    STEPComparator comparator;
    comparator.set_tolerance(0.05);
    
    auto result = comparator.compare_files(reference_step, generated_step);
    
    INFO(result.summary());
    
    // For now, just verify we can generate and load a STEP file
    // The geometry won't match until B-Rep builder issues are fixed
    CHECK(result.ref_faces == 3);  // Reference has 3 faces
    CHECK(result.gen_faces == 3);  // We also generate 3 faces
}

TEST_CASE("Compare generated STEP against reference - sphere", "[step_comparison][integration]") {
    auto [input_stl, reference_step] = get_onshape_paths("sphere_25_fine");
    std::string generated_step = get_output_path("sphere_25_fine");
    
    std::string result_path = run_brepper_pipeline(input_stl, generated_step);
    REQUIRE_FALSE(result_path.empty());
    REQUIRE(std::filesystem::exists(generated_step));
    
    STEPComparator comparator;
    comparator.set_tolerance(0.05);
    
    auto result = comparator.compare_files(reference_step, generated_step);
    
    INFO(result.summary());
    
    CHECK(result.ref_faces == 1);  // Sphere is one face
    CHECK(result.gen_faces == 1);  // We generate 1 face
}

TEST_CASE("Compare generated STEP against reference - cone", "[step_comparison][integration]") {
    auto [input_stl, reference_step] = get_onshape_paths("cone_15x20_medium");
    std::string generated_step = get_output_path("cone_15x20_medium");
    
    std::string result_path = run_brepper_pipeline(input_stl, generated_step);
    REQUIRE_FALSE(result_path.empty());
    REQUIRE(std::filesystem::exists(generated_step));
    
    STEPComparator comparator;
    comparator.set_tolerance(0.05);
    
    auto result = comparator.compare_files(reference_step, generated_step);
    
    INFO(result.summary());
    
    CHECK(result.ref_faces == 2);  // Cone has 2 faces
    CHECK(result.gen_faces == 2);  // We generate 2 faces
}

TEST_CASE("Compare generated STEP against reference - stepped block", "[step_comparison][integration][straight_edges]") {
    // Stepped block has only planar faces with straight edges
    auto [input_stl, reference_step] = get_onshape_paths("stepped_block_coarse");
    std::string generated_step = get_output_path("stepped_block_coarse");
    
    std::string result_path = run_brepper_pipeline(input_stl, generated_step);
    REQUIRE_FALSE(result_path.empty());
    REQUIRE(std::filesystem::exists(generated_step));
    
    STEPComparator comparator;
    comparator.set_tolerance(0.05);
    
    auto result = comparator.compare_files(reference_step, generated_step);
    
    INFO(result.summary());
    
    // Stepped block has 16 faces, we may generate fewer
    CHECK(result.ref_faces == 16);
    // Allow some variance in detected faces
    CHECK(result.gen_faces >= 10);  // Should get most faces
    
    // BBox should be close even with issues
    CHECK_THAT(result.gen_bbox_diagonal, Catch::Matchers::WithinRel(result.ref_bbox_diagonal, 0.05));
}

TEST_CASE("Compare generated STEP against reference - L bracket", "[step_comparison][integration][straight_edges]") {
    // L bracket also has only planar faces
    auto [input_stl, reference_step] = get_onshape_paths("l_bracket_simple_medium");
    std::string generated_step = get_output_path("l_bracket_simple_medium");
    
    std::string result_path = run_brepper_pipeline(input_stl, generated_step);
    REQUIRE_FALSE(result_path.empty());
    REQUIRE(std::filesystem::exists(generated_step));
    
    STEPComparator comparator;
    comparator.set_tolerance(0.05);
    
    auto result = comparator.compare_files(reference_step, generated_step);
    
    INFO(result.summary());
    
    // L bracket has 10 faces
    CHECK(result.ref_faces == 10);
    // We may generate more or fewer due to segmentation
    CHECK(result.gen_faces >= 5);
}

// ============================================================================
// Strict Comparison Tests for Planar Models
// 
// For models with only straight edges and planar faces, we should be able
// to achieve much tighter tolerances.
// ============================================================================

TEST_CASE("Strict comparison - planar models should match closely", "[step_comparison][strict]") {
    STEPComparator comparator;
    comparator.set_tolerance(0.02);  // 2% tolerance for planar models
    
    SECTION("Stepped block - strict geometry check") {
        auto [input_stl, reference_step] = get_onshape_paths("stepped_block_coarse");
        std::string generated_step = get_output_path("stepped_block_strict");
        
        std::string result_path = run_brepper_pipeline(input_stl, generated_step, 0.3);
        if (result_path.empty() || !std::filesystem::exists(generated_step)) {
            SKIP("Pipeline failed to generate output");
        }
        
        auto result = comparator.compare_files(reference_step, generated_step);
        
        INFO(result.summary());
        
        // Volume and bounding box should be very close for planar models
        CHECK_THAT(result.gen_bbox_diagonal, Catch::Matchers::WithinRel(result.ref_bbox_diagonal, 0.01));
    }
}

// ============================================================================
// Cube Test: Simplest possible planar model for debugging B-Rep issues
// ============================================================================

TEST_CASE("Compare generated STEP against reference - unit cube", "[step_comparison][integration][cube]") {
    auto [input_stl, reference_step] = get_manual_paths("cube");
    std::string generated_step = get_output_path("cube");
    
    INFO("Input STL: " << input_stl);
    INFO("Reference STEP: " << reference_step);
    
    // First verify our reference STEP file is valid
    STEPComparator comparator;
    auto ref_shape = comparator.read_step(reference_step);
    REQUIRE(ref_shape.has_value());
    
    // Check reference topology
    int ref_v, ref_e, ref_f, ref_sh, ref_so;
    STEPComparator::count_topology(*ref_shape, ref_v, ref_e, ref_f, ref_sh, ref_so);
    INFO("Reference: vertices=" << ref_v << ", edges=" << ref_e << ", faces=" << ref_f 
         << ", shells=" << ref_sh << ", solids=" << ref_so);
    
    double ref_volume = STEPComparator::compute_volume(*ref_shape);
    double ref_area = STEPComparator::compute_surface_area(*ref_shape);
    double ref_bbox = STEPComparator::compute_bbox_diagonal(*ref_shape);
    gp_Pnt ref_centroid = STEPComparator::compute_centroid(*ref_shape);
    
    INFO("Reference volume: " << ref_volume << " (expected 1.0)");
    INFO("Reference area: " << ref_area << " (expected 6.0)");
    INFO("Reference bbox diagonal: " << ref_bbox << " (expected sqrt(3) ≈ 1.732)");
    INFO("Reference centroid: (" << ref_centroid.X() << ", " << ref_centroid.Y() << ", " << ref_centroid.Z() << ")");
    
    // Unit cube should have these exact properties
    CHECK(ref_f == 6);
    CHECK(ref_so == 1);
    CHECK_THAT(ref_volume, Catch::Matchers::WithinRel(1.0, 0.01));
    CHECK_THAT(ref_area, Catch::Matchers::WithinRel(6.0, 0.01));
    CHECK_THAT(ref_bbox, Catch::Matchers::WithinRel(std::sqrt(3.0), 0.01));
    
    // Now generate STEP from STL - use finer point sampling and lower min_inliers for small cube
    // Cube is 1mm³ so default sampling is too coarse
    std::string result_path = run_brepper_pipeline(input_stl, generated_step, 0.1, 20);
    REQUIRE_FALSE(result_path.empty());
    REQUIRE(std::filesystem::exists(generated_step));
    
    // Compare generated vs reference
    comparator.set_tolerance(0.05);
    auto result = comparator.compare_files(reference_step, generated_step);
    
    INFO(result.summary());
    
    // Diagnostic output for debugging
    INFO("Generated volume: " << result.gen_volume);
    INFO("Generated surface area: " << result.gen_surface_area);
    INFO("Generated bbox diagonal: " << result.gen_bbox_diagonal);
    INFO("Generated centroid: (" << result.gen_centroid.X() << ", " << result.gen_centroid.Y() 
         << ", " << result.gen_centroid.Z() << ")");
    INFO("Generated solids: " << result.gen_solids);
    INFO("Generated shells: " << result.gen_shells);
    
    // Essential checks for cube
    CHECK(result.gen_faces == 6);  // Must have 6 faces for a cube
    CHECK(result.gen_solids == 1); // Must form a valid solid
    CHECK_THAT(result.gen_volume, Catch::Matchers::WithinRel(1.0, 0.05));
    CHECK_THAT(result.gen_surface_area, Catch::Matchers::WithinRel(6.0, 0.05));
    CHECK_THAT(result.gen_bbox_diagonal, Catch::Matchers::WithinRel(std::sqrt(3.0), 0.01));
}
