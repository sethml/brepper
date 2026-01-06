#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "mesh/stl_reader.hpp"
#include "common/config.hpp"
#include "test_config.hpp"
#include <pcl/conversions.h>
#include <filesystem>

using namespace brepper;
using Catch::Matchers::WithinRel;

TEST_CASE("STLReader loads valid ASCII STL file", "[stl_reader]") {
    Config config;
    STLReader reader(config);
    pcl::PolygonMesh mesh;
    
    std::string test_file = std::string(TEST_DATA_DIR) + "/manual/cube.stl";
    
    REQUIRE(std::filesystem::exists(test_file));
    REQUIRE(reader.load(test_file, mesh));
    
    SECTION("Mesh statistics are correct") {
        CHECK(reader.getNumVertices() == 8);
        CHECK(reader.getNumTriangles() == 12);
    }
    
    SECTION("All faces are triangles") {
        for (const auto& polygon : mesh.polygons) {
            CHECK(polygon.vertices.size() == 3);
        }
    }
    
    SECTION("Vertices are in expected range (unit cube)") {
        pcl::PointCloud<pcl::PointXYZ>::Ptr vertices(new pcl::PointCloud<pcl::PointXYZ>);
        pcl::fromPCLPointCloud2(mesh.cloud, *vertices);
        
        for (const auto& pt : *vertices) {
            CHECK(pt.x >= 0.0f);
            CHECK(pt.x <= 1.0f);
            CHECK(pt.y >= 0.0f);
            CHECK(pt.y <= 1.0f);
            CHECK(pt.z >= 0.0f);
            CHECK(pt.z <= 1.0f);
        }
    }
}

TEST_CASE("STLReader handles missing file", "[stl_reader]") {
    Config config;
    STLReader reader(config);
    pcl::PolygonMesh mesh;
    
    REQUIRE_FALSE(reader.load("/nonexistent/path/file.stl", mesh));
}

TEST_CASE("STLReader converts coordinates to mm based on units setting", "[stl_reader][units]") {
    // OpenCASCADE works internally in mm, so coordinates are converted on load
    
    Config config_mm;
    config_mm.stl_units = Units::Millimeters;
    
    Config config_m;
    config_m.stl_units = Units::Meters;
    
    Config config_in;
    config_in.stl_units = Units::Inches;
    
    std::string test_file = std::string(TEST_DATA_DIR) + "/manual/cube.stl";
    
    pcl::PolygonMesh mesh_mm, mesh_m, mesh_in;
    
    STLReader reader_mm(config_mm);
    STLReader reader_m(config_m);
    STLReader reader_in(config_in);
    
    REQUIRE(reader_mm.load(test_file, mesh_mm));
    REQUIRE(reader_m.load(test_file, mesh_m));
    REQUIRE(reader_in.load(test_file, mesh_in));
    
    pcl::PointCloud<pcl::PointXYZ>::Ptr vertices_mm(new pcl::PointCloud<pcl::PointXYZ>);
    pcl::PointCloud<pcl::PointXYZ>::Ptr vertices_m(new pcl::PointCloud<pcl::PointXYZ>);
    pcl::PointCloud<pcl::PointXYZ>::Ptr vertices_in(new pcl::PointCloud<pcl::PointXYZ>);
    
    pcl::fromPCLPointCloud2(mesh_mm.cloud, *vertices_mm);
    pcl::fromPCLPointCloud2(mesh_m.cloud, *vertices_m);
    pcl::fromPCLPointCloud2(mesh_in.cloud, *vertices_in);
    
    auto get_max_coord = [](const pcl::PointCloud<pcl::PointXYZ>::Ptr& verts) {
        float max_coord = 0.0f;
        for (const auto& pt : *verts) {
            max_coord = std::max({max_coord, pt.x, pt.y, pt.z});
        }
        return max_coord;
    };
    
    // Unit cube (0-1): mm stays 1.0, meters scaled to 1000.0, inches to 25.4
    CHECK_THAT(static_cast<double>(get_max_coord(vertices_mm)), WithinRel(1.0, 0.01));
    CHECK_THAT(static_cast<double>(get_max_coord(vertices_m)), WithinRel(1000.0, 0.01));
    CHECK_THAT(static_cast<double>(get_max_coord(vertices_in)), WithinRel(25.4, 0.01));
}
