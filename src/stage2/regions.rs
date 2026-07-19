//! Stage 2.2: region extraction — connected components of the face graph
//! over INTERIOR edges (union-find). Regions are maximal, mesh-connected by
//! construction, and every face belongs to exactly one region (a single
//! facet is a valid region — the terminal per-facet-plane fallback of the
//! classification ladder depends on that). Region numbering is
//! deterministic: regions are ordered by their smallest member face index.

use super::edges::EdgeClassification;
use crate::stage1::ConnectedMesh;

/// One region of the partition: a maximal set of faces connected across
/// INTERIOR edges.
#[derive(Debug, Default)]
pub struct Region {
    /// Member mesh face indices, ascending.
    pub faces: Vec<usize>,
    /// Member mesh vertex indices (union of face vertices), sorted, deduped.
    pub vertices: Vec<usize>,
}

fn find(parent: &mut [usize], mut x: usize) -> usize {
    while parent[x] != x {
        parent[x] = parent[parent[x]]; // path halving
        x = parent[x];
    }
    x
}

/// Extract regions from the edge classification. Returns
/// `(face_region, regions)` where `face_region[fi]` is the region id of
/// face `fi`.
pub fn extract_regions(
    mesh: &ConnectedMesh,
    edge_classes: &EdgeClassification,
) -> (Vec<usize>, Vec<Region>) {
    let n = mesh.faces.len();
    let mut parent: Vec<usize> = (0..n).collect();

    for fi in 0..n {
        let vc = mesh.faces[fi].vertex_count as usize;
        for e in 0..vc {
            let ni = mesh.faces[fi].neighbors[e];
            if ni >= 0 && !edge_classes.is_boundary(fi, e) {
                let ra = find(&mut parent, fi);
                let rb = find(&mut parent, ni as usize);
                if ra != rb {
                    // Union by smaller root index keeps numbering stable.
                    parent[ra.max(rb)] = ra.min(rb);
                }
            }
        }
    }

    let mut face_region = vec![usize::MAX; n];
    let mut region_of_root = vec![usize::MAX; n];
    let mut regions: Vec<Region> = Vec::new();
    for fi in 0..n {
        let root = find(&mut parent, fi);
        let id = if region_of_root[root] == usize::MAX {
            regions.push(Region::default());
            region_of_root[root] = regions.len() - 1;
            regions.len() - 1
        } else {
            region_of_root[root]
        };
        face_region[fi] = id;
        regions[id].faces.push(fi);
    }

    for region in &mut regions {
        let mut vertices: Vec<usize> = region
            .faces
            .iter()
            .flat_map(|&fi| {
                let vc = mesh.faces[fi].vertex_count as usize;
                mesh.faces[fi].vertex_indices[..vc].iter().copied()
            })
            .collect();
        vertices.sort_unstable();
        vertices.dedup();
        region.vertices = vertices;
    }

    (face_region, regions)
}
