#include "boundary_detector.hpp"
#include "common/logging.hpp"
#include <pcl/conversions.h>
#include <queue>
#include <algorithm>

namespace brepper {

BoundaryDetector::BoundaryDetector(const Config& config) : config_(config) {}

void BoundaryDetector::extract_vertices(
    const pcl::PolygonMesh& mesh,
    std::vector<Eigen::Vector3f>& vertices
) {
    pcl::PointCloud<pcl::PointXYZ> cloud;
    pcl::fromPCLPointCloud2(mesh.cloud, cloud);
    
    vertices.clear();
    vertices.reserve(cloud.size());
    for (const auto& pt : cloud) {
        vertices.emplace_back(pt.x, pt.y, pt.z);
    }
}

void BoundaryDetector::build_edge_triangle_map(
    const pcl::PolygonMesh& mesh,
    std::map<MeshEdge, std::vector<int>>& edge_to_triangles
) {
    edge_to_triangles.clear();
    
    for (size_t tri_idx = 0; tri_idx < mesh.polygons.size(); ++tri_idx) {
        const auto& polygon = mesh.polygons[tri_idx];
        if (polygon.vertices.size() < 3) continue;
        
        // Add all three edges of the triangle
        for (size_t i = 0; i < polygon.vertices.size(); ++i) {
            int v0 = polygon.vertices[i];
            int v1 = polygon.vertices[(i + 1) % polygon.vertices.size()];
            MeshEdge edge(v0, v1);
            edge_to_triangles[edge].push_back(static_cast<int>(tri_idx));
        }
    }
}

bool BoundaryDetector::detect(
    const pcl::PolygonMesh& mesh,
    const std::vector<TriangleAssignment>& assignments,
    std::vector<BoundaryEdge>& boundary_edges
) {
    LOG_DEBUG("Detecting boundary edges");
    
    if (assignments.empty()) {
        LOG_WARN("No triangle assignments provided");
        return true;
    }
    
    // Build surface ID lookup from assignments
    std::vector<int> triangle_surface(mesh.polygons.size(), -1);
    for (const auto& assignment : assignments) {
        if (assignment.triangle_id >= 0 && 
            static_cast<size_t>(assignment.triangle_id) < triangle_surface.size()) {
            triangle_surface[assignment.triangle_id] = assignment.surface_id;
        }
    }
    
    // Extract vertices
    std::vector<Eigen::Vector3f> vertices;
    extract_vertices(mesh, vertices);
    
    // Build edge-to-triangle map
    std::map<MeshEdge, std::vector<int>> edge_to_triangles;
    build_edge_triangle_map(mesh, edge_to_triangles);
    
    boundary_edges.clear();
    int interior_count = 0;
    int mesh_boundary_count = 0;
    int surface_boundary_count = 0;
    
    for (const auto& [edge, triangles] : edge_to_triangles) {
        if (triangles.size() == 1) {
            // Mesh boundary edge (only one adjacent triangle)
            int surface_id = triangle_surface[triangles[0]];
            Eigen::Vector3f midpoint = (vertices[edge.v0] + vertices[edge.v1]) / 2.0f;
            boundary_edges.emplace_back(edge, surface_id, -1, midpoint);
            ++mesh_boundary_count;
        } else if (triangles.size() == 2) {
            int surface1 = triangle_surface[triangles[0]];
            int surface2 = triangle_surface[triangles[1]];
            
            if (surface1 != surface2) {
                // Boundary between two different surfaces
                Eigen::Vector3f midpoint = (vertices[edge.v0] + vertices[edge.v1]) / 2.0f;
                boundary_edges.emplace_back(edge, surface1, surface2, midpoint);
                ++surface_boundary_count;
            } else {
                ++interior_count;
            }
        }
        // Edges with >2 triangles are non-manifold; ignore for now
    }
    
    LOG_INFO("Boundary detection: ", surface_boundary_count, " surface boundaries, ",
             mesh_boundary_count, " mesh boundaries, ", interior_count, " interior edges");
    
    return true;
}

void BoundaryDetector::find_edge_chains(
    const std::vector<BoundaryEdge>& edges,
    int surface_id_1,
    int surface_id_2,
    std::vector<std::vector<int>>& chains
) {
    chains.clear();
    
    // Filter edges for this surface pair
    std::vector<int> pair_edge_indices;
    for (size_t i = 0; i < edges.size(); ++i) {
        const auto& e = edges[i];
        if ((e.surface_id_1 == surface_id_1 && e.surface_id_2 == surface_id_2) ||
            (e.surface_id_1 == surface_id_2 && e.surface_id_2 == surface_id_1)) {
            pair_edge_indices.push_back(static_cast<int>(i));
        }
    }
    
    if (pair_edge_indices.empty()) return;
    
    // Build vertex adjacency for these edges
    std::map<int, std::vector<int>> vertex_to_edges;  // vertex -> edge indices in pair_edge_indices
    for (size_t i = 0; i < pair_edge_indices.size(); ++i) {
        const auto& edge = edges[pair_edge_indices[i]].edge;
        vertex_to_edges[edge.v0].push_back(static_cast<int>(i));
        vertex_to_edges[edge.v1].push_back(static_cast<int>(i));
    }
    
    // Find connected components using BFS
    std::vector<bool> visited(pair_edge_indices.size(), false);
    
    for (size_t start = 0; start < pair_edge_indices.size(); ++start) {
        if (visited[start]) continue;
        
        std::vector<int> chain;
        std::queue<int> queue;
        queue.push(static_cast<int>(start));
        visited[start] = true;
        
        while (!queue.empty()) {
            int edge_idx = queue.front();
            queue.pop();
            chain.push_back(pair_edge_indices[edge_idx]);
            
            const auto& edge = edges[pair_edge_indices[edge_idx]].edge;
            
            // Check neighbors through both vertices
            for (int v : {edge.v0, edge.v1}) {
                for (int neighbor_idx : vertex_to_edges[v]) {
                    if (!visited[neighbor_idx]) {
                        visited[neighbor_idx] = true;
                        queue.push(neighbor_idx);
                    }
                }
            }
        }
        
        chains.push_back(std::move(chain));
    }
}

std::vector<int> BoundaryDetector::order_chain_vertices(
    const std::vector<BoundaryEdge>& edges,
    const std::vector<int>& chain_edge_indices
) {
    if (chain_edge_indices.empty()) return {};
    if (chain_edge_indices.size() == 1) {
        const auto& e = edges[chain_edge_indices[0]].edge;
        return {e.v0, e.v1};
    }
    
    // Build vertex degree map for this chain
    std::map<int, std::vector<int>> vertex_to_edges;
    for (int edge_idx : chain_edge_indices) {
        const auto& edge = edges[edge_idx].edge;
        vertex_to_edges[edge.v0].push_back(edge_idx);
        vertex_to_edges[edge.v1].push_back(edge_idx);
    }
    
    // Find an endpoint (vertex with degree 1) or any vertex for closed loops
    int start_vertex = -1;
    for (const auto& [v, adj_edges] : vertex_to_edges) {
        if (adj_edges.size() == 1) {
            start_vertex = v;
            break;
        }
    }
    
    // If no endpoint found, it's a closed loop - start anywhere
    if (start_vertex < 0) {
        start_vertex = edges[chain_edge_indices[0]].edge.v0;
    }
    
    // Walk the chain
    std::vector<int> ordered_vertices;
    std::set<int> used_edges;
    int current_vertex = start_vertex;
    
    ordered_vertices.push_back(current_vertex);
    
    while (used_edges.size() < chain_edge_indices.size()) {
        bool found_next = false;
        
        for (int edge_idx : vertex_to_edges[current_vertex]) {
            if (used_edges.count(edge_idx)) continue;
            
            const auto& edge = edges[edge_idx].edge;
            int next_vertex = (edge.v0 == current_vertex) ? edge.v1 : edge.v0;
            
            used_edges.insert(edge_idx);
            ordered_vertices.push_back(next_vertex);
            current_vertex = next_vertex;
            found_next = true;
            break;
        }
        
        if (!found_next) break;  // Disconnected or stuck
    }
    
    return ordered_vertices;
}

bool BoundaryDetector::extract_curves(
    const pcl::PolygonMesh& mesh,
    const std::vector<BoundaryEdge>& boundary_edges,
    std::vector<BoundaryCurve>& curves
) {
    LOG_DEBUG("Extracting boundary curves");
    
    curves.clear();
    
    if (boundary_edges.empty()) {
        LOG_DEBUG("No boundary edges to process");
        return true;
    }
    
    // Extract vertices
    std::vector<Eigen::Vector3f> vertices;
    extract_vertices(mesh, vertices);
    
    // Find all unique surface pairs
    std::set<std::pair<int, int>> surface_pairs;
    for (const auto& edge : boundary_edges) {
        int s1 = std::min(edge.surface_id_1, edge.surface_id_2);
        int s2 = std::max(edge.surface_id_1, edge.surface_id_2);
        surface_pairs.insert({s1, s2});
    }
    
    // Process each surface pair
    for (const auto& [s1, s2] : surface_pairs) {
        std::vector<std::vector<int>> chains;
        find_edge_chains(boundary_edges, s1, s2, chains);
        
        for (const auto& chain_edges : chains) {
            std::vector<int> ordered_verts = order_chain_vertices(boundary_edges, chain_edges);
            
            if (ordered_verts.size() < 2) continue;
            
            BoundaryCurve curve;
            curve.surface_id_left = s1;
            curve.surface_id_right = s2;
            
            for (int v : ordered_verts) {
                if (static_cast<size_t>(v) < vertices.size()) {
                    const auto& pt = vertices[v];
                    curve.points.emplace_back(pt.x(), pt.y(), pt.z());
                }
            }
            
            // Store edge IDs
            for (int edge_idx : chain_edges) {
                curve.edge_ids.push_back(edge_idx);
            }
            
            if (!curve.points.empty()) {
                curves.push_back(std::move(curve));
            }
        }
    }
    
    LOG_INFO("Extracted ", curves.size(), " boundary curves from ", surface_pairs.size(), " surface pairs");
    
    return true;
}

} // namespace brepper
