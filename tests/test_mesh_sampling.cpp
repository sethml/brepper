#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "mesh/stl_reader.hpp"
#include "mesh/mesh_sampling.hpp"
#include "common/config.hpp"
#include "test_config.hpp"
#include <cmath>

using namespace brepper;
using Catch::Matchers::WithinRel;

// Helper to load mesh for sampling tests
pcl::PolygonMesh load_cube() {
    Config config;
    STLReader reader(config);
    pcl::PolygonMesh mesh;
    std::string test_file = std::string(TEST_DATA_DIR) + "/manual/cube.stl";
    reader.load(test_file, mesh);
    return mesh;
}

TEST_CASE("MeshSampler generates points from cube", "[mesh_sampling]") {
    Config config;
    config.max_point_distance_mm = 0.5;  // Sample at 0.5mm spacing
    config.min_samples_per_triangle = 1;
    
    pcl::PolygonMesh mesh = load_cube();
    MeshSampler sampler(config);
    PointCloudNormalPtr cloud;
    
    REQUIRE(sampler.sample(mesh, cloud));
    
    SECTION("Points are generated") {
        CHECK(cloud->size() > 0);
        CHECK(sampler.getNumSampledPoints() > 0);
    }
    
    SECTION("All points have valid normals") {
        for (const auto& pt : *cloud) {
            float normal_len = std::sqrt(
                pt.normal_x * pt.normal_x +
                pt.normal_y * pt.normal_y +
                pt.normal_z * pt.normal_z
            );
            CHECK_THAT(static_cast<double>(normal_len), WithinRel(1.0, 0.01));
        }
    }
    
    SECTION("Normals are consistent for each face") {
        // For a cube, each face should have all normals pointing the same direction
        // STL winding order determines if they point inward or outward
        // We just verify they are unit vectors and axis-aligned for a cube
        
        for (const auto& pt : *cloud) {
            Eigen::Vector3f normal(pt.normal_x, pt.normal_y, pt.normal_z);
            float norm = normal.norm();
            
            // Normal should be unit length
            CHECK_THAT(static_cast<double>(norm), WithinRel(1.0, 0.01));
            
            // For axis-aligned cube, normal should be mostly along one axis
            float max_component = std::max({std::abs(normal.x()), 
                                            std::abs(normal.y()), 
                                            std::abs(normal.z())});
            CHECK_THAT(static_cast<double>(max_component), WithinRel(1.0, 0.01));
        }
    }
}

TEST_CASE("MeshSampler respects max_point_distance", "[mesh_sampling]") {
    pcl::PolygonMesh mesh = load_cube();
    
    SECTION("Smaller max_point_distance produces more samples") {
        Config config_coarse;
        config_coarse.max_point_distance_mm = 1.0;  // Coarse: 1mm spacing
        config_coarse.min_samples_per_triangle = 1;
        
        Config config_fine;
        config_fine.max_point_distance_mm = 0.25;  // Fine: 0.25mm spacing
        config_fine.min_samples_per_triangle = 1;
        
        MeshSampler sampler_coarse(config_coarse);
        MeshSampler sampler_fine(config_fine);
        
        PointCloudNormalPtr cloud_coarse, cloud_fine;
        REQUIRE(sampler_coarse.sample(mesh, cloud_coarse));
        REQUIRE(sampler_fine.sample(mesh, cloud_fine));
        
        // Fine sampling should produce significantly more points
        CHECK(cloud_fine->size() > cloud_coarse->size() * 2);
    }
    
    SECTION("Very fine sampling on large face") {
        // Unit cube face has area 1.0
        // With max_point_distance_mm = 0.1, we expect roughly (1/0.1)^2 = 100 samples per face
        // Plus vertices and centroid = 4 per triangle, 2 triangles per face = 8 base
        // 6 faces * ~100 = ~600+ samples total
        Config config;
        config.max_point_distance_mm = 0.1;
        config.min_samples_per_triangle = 1;
        
        MeshSampler sampler(config);
        PointCloudNormalPtr cloud;
        REQUIRE(sampler.sample(mesh, cloud));
        
        // Should have many points
        CHECK(cloud->size() > 100);
    }
}

TEST_CASE("MeshSampler handles min_samples_per_triangle", "[mesh_sampling]") {
    pcl::PolygonMesh mesh = load_cube();
    
    Config config;
    config.max_point_distance_mm = 10.0;  // Very coarse, relies on minimum
    config.min_samples_per_triangle = 10;
    
    MeshSampler sampler(config);
    PointCloudNormalPtr cloud;
    REQUIRE(sampler.sample(mesh, cloud));
    
    // 12 triangles * at least min_samples_per_triangle = 120+ points
    // Actually we add centroid + 3 vertices + (min_samples - 4) additional
    // So 12 * 10 = 120 minimum
    CHECK(cloud->size() >= 12 * config.min_samples_per_triangle);
}

TEST_CASE("MeshSampler computes reasonable triangle area", "[mesh_sampling]") {
    pcl::PolygonMesh mesh = load_cube();
    
    Config config;
    config.max_point_distance_mm = 1.0;
    
    MeshSampler sampler(config);
    PointCloudNormalPtr cloud;
    REQUIRE(sampler.sample(mesh, cloud));
    
    // Unit cube has 6 faces, each face is 2 triangles
    // Each face has area 1.0, each triangle has area 0.5
    // Average triangle area should be 0.5
    double avg_area = sampler.getAverageTriangleArea();
    CHECK_THAT(avg_area, WithinRel(0.5, 0.01));
}

TEST_CASE("MeshSampler points stay within mesh bounds", "[mesh_sampling]") {
    pcl::PolygonMesh mesh = load_cube();
    
    Config config;
    config.max_point_distance_mm = 0.1;  // Dense sampling
    config.min_samples_per_triangle = 5;
    
    MeshSampler sampler(config);
    PointCloudNormalPtr cloud;
    REQUIRE(sampler.sample(mesh, cloud));
    
    // All points should be within the unit cube bounds (with small epsilon)
    const float eps = 0.001f;
    for (const auto& pt : *cloud) {
        CHECK(pt.x >= -eps);
        CHECK(pt.x <= 1.0f + eps);
        CHECK(pt.y >= -eps);
        CHECK(pt.y <= 1.0f + eps);
        CHECK(pt.z >= -eps);
        CHECK(pt.z <= 1.0f + eps);
    }
}
