#pragma once

#include "common/types.hpp"
#include "common/config.hpp"
#include <map>
#include <set>

namespace brepper {

// Represents a mesh edge as an ordered pair of vertex indices
struct MeshEdge {
    int v0, v1;  // Vertex indices (v0 < v1 for canonical form)
    
    MeshEdge(int a, int b) : v0(std::min(a, b)), v1(std::max(a, b)) {}
    
    bool operator<(const MeshEdge& other) const {
        if (v0 != other.v0) return v0 < other.v0;
        return v1 < other.v1;
    }
    
    bool operator==(const MeshEdge& other) const {
        return v0 == other.v0 && v1 == other.v1;
    }
};

// Information about a boundary edge
struct BoundaryEdge {
    MeshEdge edge;
    int surface_id_1;  // Surface on one side
    int surface_id_2;  // Surface on other side (-1 if mesh boundary)
    Eigen::Vector3f midpoint;
    
    BoundaryEdge(const MeshEdge& e, int s1, int s2, const Eigen::Vector3f& mid)
        : edge(e), surface_id_1(s1), surface_id_2(s2), midpoint(mid) {}
};

class BoundaryDetector {
public:
    explicit BoundaryDetector(const Config& config);
    
    // Detect boundary edges between different surfaces
    bool detect(
        const pcl::PolygonMesh& mesh,
        const std::vector<TriangleAssignment>& assignments,
        std::vector<BoundaryEdge>& boundary_edges
    );
    
    // Group boundary edges into connected curves
    bool extract_curves(
        const pcl::PolygonMesh& mesh,
        const std::vector<BoundaryEdge>& boundary_edges,
        std::vector<BoundaryCurve>& curves
    );
    
private:
    const Config& config_;
    
    // Build map from edge to adjacent triangles
    void build_edge_triangle_map(
        const pcl::PolygonMesh& mesh,
        std::map<MeshEdge, std::vector<int>>& edge_to_triangles
    );
    
    // Extract vertex positions from mesh
    void extract_vertices(
        const pcl::PolygonMesh& mesh,
        std::vector<Eigen::Vector3f>& vertices
    );
    
    // Find connected components of edges sharing surface pair
    void find_edge_chains(
        const std::vector<BoundaryEdge>& edges,
        int surface_id_1,
        int surface_id_2,
        std::vector<std::vector<int>>& chains
    );
    
    // Order vertices along an edge chain
    std::vector<int> order_chain_vertices(
        const std::vector<BoundaryEdge>& edges,
        const std::vector<int>& chain_edge_indices
    );
};

} // namespace brepper
