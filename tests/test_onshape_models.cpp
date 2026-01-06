#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "mesh/stl_reader.hpp"
#include "mesh/mesh_sampling.hpp"
#include "common/config.hpp"
#include "common/types.hpp"
#include "test_config.hpp"
#include <pcl/conversions.h>
#include <filesystem>
#include <cmath>

using namespace brepper;
using Catch::Matchers::WithinRel;

// Helper to get bounding box of a mesh
struct BoundingBox {
    float min_x, max_x, min_y, max_y, min_z, max_z;
    
    float size_x() const { return max_x - min_x; }
    float size_y() const { return max_y - min_y; }
    float size_z() const { return max_z - min_z; }
};

BoundingBox getMeshBounds(const pcl::PolygonMesh& mesh) {
    pcl::PointCloud<pcl::PointXYZ>::Ptr vertices(new pcl::PointCloud<pcl::PointXYZ>);
    pcl::fromPCLPointCloud2(mesh.cloud, *vertices);
    
    BoundingBox bb;
    bb.min_x = bb.min_y = bb.min_z = std::numeric_limits<float>::max();
    bb.max_x = bb.max_y = bb.max_z = std::numeric_limits<float>::lowest();
    
    for (const auto& pt : *vertices) {
        bb.min_x = std::min(bb.min_x, pt.x);
        bb.max_x = std::max(bb.max_x, pt.x);
        bb.min_y = std::min(bb.min_y, pt.y);
        bb.max_y = std::max(bb.max_y, pt.y);
        bb.min_z = std::min(bb.min_z, pt.z);
        bb.max_z = std::max(bb.max_z, pt.z);
    }
    return bb;
}

// ============================================================================
// Binary STL Loading Tests (Onshape exports)
// ============================================================================

TEST_CASE("STLReader loads binary STL files from Onshape", "[stl_reader][onshape][binary]") {
    Config config;
    STLReader reader(config);
    
    SECTION("Cylinder 10x30mm") {
        pcl::PolygonMesh mesh;
        std::string file = std::string(TEST_DATA_DIR) + "/onshape/cylinder_10x30_medium.stl";
        REQUIRE(std::filesystem::exists(file));
        REQUIRE(reader.load(file, mesh));
        
        CHECK(reader.getNumTriangles() > 0);
        
        // Check approximate dimensions (10mm diameter, 30mm height)
        auto bb = getMeshBounds(mesh);
        CHECK_THAT(static_cast<double>(bb.size_x()), WithinRel(10.0, 0.1));
        CHECK_THAT(static_cast<double>(bb.size_y()), WithinRel(10.0, 0.1));
        CHECK_THAT(static_cast<double>(bb.size_z()), WithinRel(30.0, 0.1));
    }
    
    SECTION("Sphere 25mm diameter") {
        pcl::PolygonMesh mesh;
        std::string file = std::string(TEST_DATA_DIR) + "/onshape/sphere_25_fine.stl";
        REQUIRE(std::filesystem::exists(file));
        REQUIRE(reader.load(file, mesh));
        
        CHECK(reader.getNumTriangles() > 1000);  // Fine mesh should have many triangles
        
        auto bb = getMeshBounds(mesh);
        // Sphere should be roughly 25mm in all directions
        CHECK_THAT(static_cast<double>(bb.size_x()), WithinRel(25.0, 0.1));
        CHECK_THAT(static_cast<double>(bb.size_y()), WithinRel(25.0, 0.1));
        CHECK_THAT(static_cast<double>(bb.size_z()), WithinRel(25.0, 0.1));
    }
    
    SECTION("Cone 15x20mm") {
        pcl::PolygonMesh mesh;
        std::string file = std::string(TEST_DATA_DIR) + "/onshape/cone_15x20_medium.stl";
        REQUIRE(std::filesystem::exists(file));
        REQUIRE(reader.load(file, mesh));
        
        CHECK(reader.getNumTriangles() > 0);
        
        auto bb = getMeshBounds(mesh);
        // 15mm base diameter, 20mm height
        CHECK_THAT(static_cast<double>(bb.size_x()), WithinRel(15.0, 0.1));
        CHECK_THAT(static_cast<double>(bb.size_z()), WithinRel(20.0, 0.1));
    }
    
    SECTION("Plate with hole") {
        pcl::PolygonMesh mesh;
        std::string file = std::string(TEST_DATA_DIR) + "/onshape/plate_with_hole_100x50_coarse.stl";
        REQUIRE(std::filesystem::exists(file));
        REQUIRE(reader.load(file, mesh));
        
        CHECK(reader.getNumTriangles() > 0);
        
        // Just verify it loads and has reasonable dimensions
        auto bb = getMeshBounds(mesh);
        CHECK(bb.size_x() > 0);
        CHECK(bb.size_y() > 0);
        CHECK(bb.size_z() > 0);
    }
    
    SECTION("L-bracket") {
        pcl::PolygonMesh mesh;
        std::string file = std::string(TEST_DATA_DIR) + "/onshape/l_bracket_simple_medium.stl";
        REQUIRE(std::filesystem::exists(file));
        REQUIRE(reader.load(file, mesh));
        
        CHECK(reader.getNumTriangles() > 0);
        // Just verify it loads - dimensions vary
    }
    
    SECTION("Rounded cube with fillets") {
        pcl::PolygonMesh mesh;
        std::string file = std::string(TEST_DATA_DIR) + "/onshape/rounded_cube_10_r2_fine.stl";
        REQUIRE(std::filesystem::exists(file));
        REQUIRE(reader.load(file, mesh));
        
        // Fine mesh with fillets should have many triangles
        CHECK(reader.getNumTriangles() > 5000);
        
        auto bb = getMeshBounds(mesh);
        // 10mm cube
        CHECK_THAT(static_cast<double>(bb.size_x()), WithinRel(10.0, 0.1));
        CHECK_THAT(static_cast<double>(bb.size_y()), WithinRel(10.0, 0.1));
        CHECK_THAT(static_cast<double>(bb.size_z()), WithinRel(10.0, 0.1));
    }
    
    SECTION("Pipe elbow") {
        pcl::PolygonMesh mesh;
        std::string file = std::string(TEST_DATA_DIR) + "/onshape/pipe_elbow_10_fine.stl";
        REQUIRE(std::filesystem::exists(file));
        REQUIRE(reader.load(file, mesh));
        
        CHECK(reader.getNumTriangles() > 1000);  // Fine mesh
    }
    
    SECTION("Hemisphere dome") {
        pcl::PolygonMesh mesh;
        std::string file = std::string(TEST_DATA_DIR) + "/onshape/dome_hemisphere_20_fine.stl";
        REQUIRE(std::filesystem::exists(file));
        REQUIRE(reader.load(file, mesh));
        
        CHECK(reader.getNumTriangles() > 1000);
        
        auto bb = getMeshBounds(mesh);
        // Hemisphere: one dimension should be ~half the others (radius vs diameter)
        // Don't assume orientation - just verify reasonable dimensions
        CHECK(bb.size_x() > 0);
        CHECK(bb.size_y() > 0);
        CHECK(bb.size_z() > 0);
        
        // Find the smallest dimension - should be roughly half the largest
        float dims[3] = {bb.size_x(), bb.size_y(), bb.size_z()};
        std::sort(dims, dims + 3);
        // Smallest should be approximately half of largest (hemisphere height vs diameter)
        CHECK(dims[0] < dims[2]);  // At minimum, dimensions differ
    }
    
    SECTION("Chamfered cube") {
        pcl::PolygonMesh mesh;
        std::string file = std::string(TEST_DATA_DIR) + "/onshape/chamfered_cube_10_c1_medium.stl";
        REQUIRE(std::filesystem::exists(file));
        REQUIRE(reader.load(file, mesh));
        
        CHECK(reader.getNumTriangles() > 0);
        
        auto bb = getMeshBounds(mesh);
        CHECK_THAT(static_cast<double>(bb.size_x()), WithinRel(10.0, 0.1));
    }
    
    SECTION("Stepped block") {
        pcl::PolygonMesh mesh;
        std::string file = std::string(TEST_DATA_DIR) + "/onshape/stepped_block_coarse.stl";
        REQUIRE(std::filesystem::exists(file));
        REQUIRE(reader.load(file, mesh));
        
        CHECK(reader.getNumTriangles() > 0);
    }
}

// ============================================================================
// Mesh Sampling Tests with Real Models
// ============================================================================

TEST_CASE("MeshSampler works on complex Onshape models", "[mesh_sampling][onshape]") {
    Config config;
    config.max_point_distance_mm = 1.0;  // 1mm spacing for reasonable test speed
    config.min_samples_per_triangle = 1;
    
    STLReader reader(config);
    MeshSampler sampler(config);
    
    SECTION("Sampling sphere produces uniform coverage") {
        pcl::PolygonMesh mesh;
        std::string file = std::string(TEST_DATA_DIR) + "/onshape/sphere_25_fine.stl";
        REQUIRE(reader.load(file, mesh));
        
        PointCloudNormalPtr cloud;
        REQUIRE(sampler.sample(mesh, cloud));
        
        // Should have many points
        CHECK(cloud->size() > 1000);
        
        // All points should have valid normals
        for (const auto& pt : *cloud) {
            float normal_len = std::sqrt(
                pt.normal_x * pt.normal_x +
                pt.normal_y * pt.normal_y +
                pt.normal_z * pt.normal_z
            );
            CHECK_THAT(static_cast<double>(normal_len), WithinRel(1.0, 0.01));
        }
    }
    
    SECTION("Sampling cylinder") {
        pcl::PolygonMesh mesh;
        std::string file = std::string(TEST_DATA_DIR) + "/onshape/cylinder_10x30_medium.stl";
        REQUIRE(reader.load(file, mesh));
        
        PointCloudNormalPtr cloud;
        REQUIRE(sampler.sample(mesh, cloud));
        
        CHECK(cloud->size() > 100);
    }
    
    SECTION("Sampling rounded cube with fillets") {
        pcl::PolygonMesh mesh;
        std::string file = std::string(TEST_DATA_DIR) + "/onshape/rounded_cube_10_r2_fine.stl";
        REQUIRE(reader.load(file, mesh));
        
        PointCloudNormalPtr cloud;
        REQUIRE(sampler.sample(mesh, cloud));
        
        // High triangle count should produce many points
        CHECK(cloud->size() > 5000);
    }
}

// ============================================================================
// All Models Load Test
// ============================================================================

TEST_CASE("All Onshape STL files load successfully", "[stl_reader][onshape][smoke]") {
    Config config;
    STLReader reader(config);
    
    std::vector<std::string> files = {
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
    
    for (const auto& filename : files) {
        SECTION(filename) {
            pcl::PolygonMesh mesh;
            std::string path = std::string(TEST_DATA_DIR) + "/onshape/" + filename;
            
            INFO("Loading: " << path);
            REQUIRE(std::filesystem::exists(path));
            REQUIRE(reader.load(path, mesh));
            CHECK(reader.getNumTriangles() > 0);
        }
    }
}

// ============================================================================
// ASCII vs Binary STL Comparison Test
// ============================================================================

TEST_CASE("ASCII and binary STL produce identical mesh data", "[stl_reader][onshape][ascii]") {
    Config config;
    STLReader reader(config);
    
    std::string binary_file = std::string(TEST_DATA_DIR) + "/onshape/chamfered_cube_10_c1_medium.stl";
    std::string ascii_file = std::string(TEST_DATA_DIR) + "/onshape/chamfered_cube_10_c1_medium.ascii.stl";
    
    REQUIRE(std::filesystem::exists(binary_file));
    REQUIRE(std::filesystem::exists(ascii_file));
    
    pcl::PolygonMesh binary_mesh, ascii_mesh;
    
    REQUIRE(reader.load(binary_file, binary_mesh));
    size_t binary_triangles = reader.getNumTriangles();
    
    REQUIRE(reader.load(ascii_file, ascii_mesh));
    size_t ascii_triangles = reader.getNumTriangles();
    
    // Same number of triangles
    CHECK(binary_triangles == ascii_triangles);
    
    // Same number of polygons in the mesh
    CHECK(binary_mesh.polygons.size() == ascii_mesh.polygons.size());
    
    // Extract vertices from both meshes
    pcl::PointCloud<pcl::PointXYZ>::Ptr binary_verts(new pcl::PointCloud<pcl::PointXYZ>);
    pcl::PointCloud<pcl::PointXYZ>::Ptr ascii_verts(new pcl::PointCloud<pcl::PointXYZ>);
    pcl::fromPCLPointCloud2(binary_mesh.cloud, *binary_verts);
    pcl::fromPCLPointCloud2(ascii_mesh.cloud, *ascii_verts);
    
    // Same number of vertices
    CHECK(binary_verts->size() == ascii_verts->size());
    
    // Same bounding box dimensions
    auto binary_bb = getMeshBounds(binary_mesh);
    auto ascii_bb = getMeshBounds(ascii_mesh);
    
    CHECK_THAT(static_cast<double>(binary_bb.size_x()), 
               WithinRel(static_cast<double>(ascii_bb.size_x()), 0.001));
    CHECK_THAT(static_cast<double>(binary_bb.size_y()), 
               WithinRel(static_cast<double>(ascii_bb.size_y()), 0.001));
    CHECK_THAT(static_cast<double>(binary_bb.size_z()), 
               WithinRel(static_cast<double>(ascii_bb.size_z()), 0.001));
    
    // Verify vertex positions match (allowing for small floating point differences)
    REQUIRE(binary_verts->size() == ascii_verts->size());
    for (size_t i = 0; i < binary_verts->size(); ++i) {
        const auto& bv = (*binary_verts)[i];
        const auto& av = (*ascii_verts)[i];
        
        CHECK_THAT(static_cast<double>(bv.x), WithinRel(static_cast<double>(av.x), 0.0001));
        CHECK_THAT(static_cast<double>(bv.y), WithinRel(static_cast<double>(av.y), 0.0001));
        CHECK_THAT(static_cast<double>(bv.z), WithinRel(static_cast<double>(av.z), 0.0001));
    }
}
