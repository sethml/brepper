/// Compute maximum distance between STL mesh vertices/centroids and STEP model surfaces,
/// and maximum angular deviation between adjacent STL faces on the same STEP surface.
///
/// Reads an STL file to extract vertex positions and triangle centroids, reads a
/// STEP file to extract surfaces, and computes the minimum distance from each
/// point to any STEP surface. Also builds triangle adjacency and maps each
/// triangle to its nearest STEP surface, then measures the dihedral angle
/// between adjacent triangles sharing the same STEP surface.
///
/// Units: all distances are in mm. OCCT's STEPControl_Reader converts STEP
/// file units (typically meters for Onshape exports) to mm internally.
/// STL files have no unit metadata; coordinates are assumed to be in mm.
use brepper::stage1::{self, VertexWeldOptions};
use opencascade_sys::{
    b_rep, b_rep_builder_api, b_rep_extrema, b_rep_g_prop, extrema, g_prop, gp, message, rw_stl,
    step_control, top_abs, top_exp, topo_ds,
};
use std::collections::HashMap;
use std::env;
use std::process;

fn main() {
    let args: Vec<String> = env::args().collect();

    // Parse --quad flag
    let mut positional = Vec::new();
    let mut use_quad = false;
    for arg in &args[1..] {
        if arg == "--quad" {
            use_quad = true;
        } else {
            positional.push(arg.as_str());
        }
    }
    if positional.len() != 2 {
        eprintln!("Usage: {} [--quad] <input.stl> <input.step>", args[0]);
        process::exit(1);
    }
    let stl_path = positional[0];
    let step_path = positional[1];

    // Read STL → Poly_Triangulation
    if !std::path::Path::new(stl_path).exists() {
        eprintln!("Error: STL file not found: {stl_path}");
        process::exit(1);
    }
    let progress = message::ProgressRange::new();
    let tri_handle = rw_stl::read_file_charptr_progressrange_2(stl_path, &progress);
    let tri = tri_handle.get();
    let num_nodes = tri.nb_nodes();
    let num_triangles = tri.nb_triangles();
    eprintln!("STL: {} nodes, {} triangles", num_nodes, num_triangles);

    if num_nodes == 0 {
        eprintln!("Error: STL file has no vertices");
        process::exit(1);
    }

    // Read STEP
    if !std::path::Path::new(step_path).exists() {
        eprintln!("Error: STEP file not found: {step_path}");
        process::exit(1);
    }
    let mut reader = step_control::StepReader::new();
    let read_status = reader.read_file_charptr(step_path);
    if read_status != opencascade_sys::if_select::ReturnStatus::Retdone {
        eprintln!("Error: Failed to read STEP file {step_path}: {read_status:?}");
        process::exit(1);
    }
    let progress2 = message::ProgressRange::new();
    reader.transfer_roots(&progress2);
    let step_shape = reader.one_shape();

    // Collect STEP faces and their surfaces
    let mut step_surfaces: Vec<opencascade_sys::OwnedPtr<opencascade_sys::shape_extend::HandleGeomSurface>> = Vec::new();
    {
        let mut explorer = top_exp::Explorer::new_shape_shapeenum2(
            &step_shape,
            top_abs::ShapeEnum::Face,
            top_abs::ShapeEnum::Shape,
        );
        while explorer.more() {
            let face = topo_ds::face_shape(explorer.value());
            let surface = b_rep::Tool::surface_face(face);
            step_surfaces.push(surface);
            explorer.next();
        }
    }
    let face_count = step_surfaces.len();
    eprintln!("STEP: {} faces", face_count);

    if face_count == 0 {
        eprintln!("Error: STEP file has no faces");
        process::exit(1);
    }

    // Compute STEP surface area and volume
    let mut step_area_props = g_prop::GProps::new();
    b_rep_g_prop::surface_properties_shape_gprops_bool2(
        &step_shape,
        &mut step_area_props,
        false,
        false,
    );
    let step_area = step_area_props.mass();

    let mut step_vol_props = g_prop::GProps::new();
    b_rep_g_prop::volume_properties_shape_gprops_bool3(
        &step_shape,
        &mut step_vol_props,
        true,
        false,
        false,
    );
    let step_volume = step_vol_props.mass();

    eprintln!("STEP area: {:.6} mm²  volume: {:.6} mm³", step_area, step_volume);

    // Compute STL bounding box and max dimension
    let mut bb_min = [f64::MAX; 3];
    let mut bb_max = [f64::MIN; 3];
    for i in 1..=num_nodes {
        let pt = tri.node(i);
        let coords = [pt.x(), pt.y(), pt.z()];
        for d in 0..3 {
            bb_min[d] = bb_min[d].min(coords[d]);
            bb_max[d] = bb_max[d].max(coords[d]);
        }
    }
    let extents = [
        bb_max[0] - bb_min[0],
        bb_max[1] - bb_min[1],
        bb_max[2] - bb_min[2],
    ];
    let max_dimension = extents[0].max(extents[1]).max(extents[2]);
    eprintln!("Bounding box: {:.6} x {:.6} x {:.6}", extents[0], extents[1], extents[2]);
    eprintln!("Max dimension: {:.6}", max_dimension);

    // Compute STL surface area and volume
    // Volume via signed tetrahedra method (each triangle forms a tet with the origin)
    let mut stl_area = 0.0_f64;
    let mut stl_volume = 0.0_f64;

    for i in 1..=num_triangles {
        let triangle = tri.triangle(i);
        let mut n1 = 0_i32;
        let mut n2 = 0_i32;
        let mut n3 = 0_i32;
        triangle.get(&mut n1, &mut n2, &mut n3);
        let p1 = tri.node(n1);
        let p2 = tri.node(n2);
        let p3 = tri.node(n3);

        // Triangle area = 0.5 * |cross(p2-p1, p3-p1)|
        let e1 = [p2.x() - p1.x(), p2.y() - p1.y(), p2.z() - p1.z()];
        let e2 = [p3.x() - p1.x(), p3.y() - p1.y(), p3.z() - p1.z()];
        let cx = e1[1] * e2[2] - e1[2] * e2[1];
        let cy = e1[2] * e2[0] - e1[0] * e2[2];
        let cz = e1[0] * e2[1] - e1[1] * e2[0];
        stl_area += (cx * cx + cy * cy + cz * cz).sqrt() * 0.5;

        // Signed volume contribution = (1/6) * p1 . (p2 x p3)
        let cross_x = p2.y() * p3.z() - p2.z() * p3.y();
        let cross_y = p2.z() * p3.x() - p2.x() * p3.z();
        let cross_z = p2.x() * p3.y() - p2.y() * p3.x();
        stl_volume += (p1.x() * cross_x + p1.y() * cross_y + p1.z() * cross_z) / 6.0;
    }
    stl_volume = stl_volume.abs();

    eprintln!("STL  area: {:.6} mm\u{00b2}  volume: {:.6} mm\u{00b3}", stl_area, stl_volume);

    let min_distance_to_step = |pt: &gp::Pnt| -> f64 {
        let mut vertex = b_rep_builder_api::MakeVertex::new_pnt(pt);
        let vtx_shape = vertex.vertex();
        let progress = message::ProgressRange::new();
        let dist_calc = b_rep_extrema::DistShapeShape::new_shape2_extflag_extalgo_progressrange(
            vtx_shape.as_shape(),
            &step_shape,
            0,
            extrema::ExtAlgo::Grad,
            &progress,
        );
        if dist_calc.is_done() && dist_calc.nb_solution() > 0 {
            dist_calc.value()
        } else {
            f64::MAX
        }
    };

    // Vertex distances
    let mut vtx_max_dist = 0.0_f64;
    let mut vtx_max_idx = 1_i32;
    let mut vtx_dist_sum = 0.0_f64;

    for i in 1..=num_nodes {
        let pt = tri.node(i);
        let min_dist = min_distance_to_step(&pt);
        if min_dist < f64::MAX {
            vtx_dist_sum += min_dist;
            if min_dist > vtx_max_dist {
                vtx_max_dist = min_dist;
                vtx_max_idx = i;
            }
        }
    }

    // Centroid distances — either from raw triangles, or from quad-fused faces
    let mut ctr_max_dist = 0.0_f64;
    let mut ctr_dist_sum = 0.0_f64;
    let mut ctr_count = 0_usize;

    if use_quad {
        // Use brepper's mesh reader + quad fusing for centroid computation
        let mut mesh = stage1::read_connected_mesh_from_stl(
            stl_path,
            VertexWeldOptions { tolerance: 1.0e-9 },
        )
        .unwrap_or_else(|e| {
            eprintln!("Error reading STL for quad fusing: {e:?}");
            process::exit(1);
        });
        mesh.validate_and_populate_topology().unwrap_or_else(|e| {
            eprintln!("Error validating mesh for quad fusing: {e:?}");
            process::exit(1);
        });
        stage1::fuse_coplanar_triangles(&mut mesh, 1.0e-5);
        let num_quads = mesh.faces.iter().filter(|f| f.vertex_count == 4).count();
        let num_active = mesh.faces.iter().filter(|f| f.vertex_count > 0).count();
        eprintln!(
            "Quad fusing: {} faces ({} quads, {} tris)",
            num_active,
            num_quads,
            num_active - num_quads
        );
        for face in &mesh.faces {
            if face.vertex_count == 0 {
                continue;
            }
            let n = face.vertex_count as usize;
            let mut cx = 0.0_f64;
            let mut cy = 0.0_f64;
            let mut cz = 0.0_f64;
            for vi_idx in 0..n {
                let v = &mesh.vertices[face.vertex_indices[vi_idx]];
                cx += v.x;
                cy += v.y;
                cz += v.z;
            }
            let inv_n = 1.0 / n as f64;
            let centroid = gp::Pnt::new_real3(cx * inv_n, cy * inv_n, cz * inv_n);
            let min_dist = min_distance_to_step(&centroid);
            if min_dist < f64::MAX {
                ctr_dist_sum += min_dist;
                ctr_count += 1;
                if min_dist > ctr_max_dist {
                    ctr_max_dist = min_dist;
                }
            }
        }
    } else {
        ctr_count = num_triangles as usize;
        for i in 1..=num_triangles {
            let triangle = tri.triangle(i);
            let mut n1 = 0_i32;
            let mut n2 = 0_i32;
            let mut n3 = 0_i32;
            triangle.get(&mut n1, &mut n2, &mut n3);
            let p1 = tri.node(n1);
            let p2 = tri.node(n2);
            let p3 = tri.node(n3);
            let centroid = gp::Pnt::new_real3(
                (p1.x() + p2.x() + p3.x()) / 3.0,
                (p1.y() + p2.y() + p3.y()) / 3.0,
                (p1.z() + p2.z() + p3.z()) / 3.0,
            );
            let min_dist = min_distance_to_step(&centroid);
            if min_dist < f64::MAX {
                ctr_dist_sum += min_dist;
                if min_dist > ctr_max_dist {
                    ctr_max_dist = min_dist;
                }
            }
        }
    }

    // --- Angular tolerance measurement ---
    // Build triangle adjacency: for each edge (min_node, max_node), collect triangle indices
    let mut edge_to_triangles: HashMap<(i32, i32), Vec<i32>> = HashMap::new();
    // Compute face normals and centroids for each triangle
    let mut tri_normals: Vec<[f64; 3]> = Vec::with_capacity(num_triangles as usize);
    let mut tri_centroids: Vec<[f64; 3]> = Vec::with_capacity(num_triangles as usize);

    for i in 1..=num_triangles {
        let triangle = tri.triangle(i);
        let mut n1 = 0_i32;
        let mut n2 = 0_i32;
        let mut n3 = 0_i32;
        triangle.get(&mut n1, &mut n2, &mut n3);

        let p1 = tri.node(n1);
        let p2 = tri.node(n2);
        let p3 = tri.node(n3);

        // Centroid
        let cx = (p1.x() + p2.x() + p3.x()) / 3.0;
        let cy = (p1.y() + p2.y() + p3.y()) / 3.0;
        let cz = (p1.z() + p2.z() + p3.z()) / 3.0;
        tri_centroids.push([cx, cy, cz]);

        // Face normal via cross product
        let e1 = [p2.x() - p1.x(), p2.y() - p1.y(), p2.z() - p1.z()];
        let e2 = [p3.x() - p1.x(), p3.y() - p1.y(), p3.z() - p1.z()];
        let nx = e1[1] * e2[2] - e1[2] * e2[1];
        let ny = e1[2] * e2[0] - e1[0] * e2[2];
        let nz = e1[0] * e2[1] - e1[1] * e2[0];
        let len = (nx * nx + ny * ny + nz * nz).sqrt();
        if len > 1e-30 {
            tri_normals.push([nx / len, ny / len, nz / len]);
        } else {
            tri_normals.push([0.0, 0.0, 0.0]); // degenerate
        }

        // Register edges (sorted node indices)
        let edges = [(n1, n2), (n2, n3), (n1, n3)];
        for (a, b) in edges {
            let key = (a.min(b), a.max(b));
            edge_to_triangles.entry(key).or_default().push(i);
        }
    }

    // Map each triangle to its nearest STEP surface using GeomAPI_ProjectPointOnSurf
    let mut tri_surface_idx: Vec<usize> = Vec::with_capacity(num_triangles as usize);
    use opencascade_sys::geom_api;
    for tc in &tri_centroids {
        let pt = gp::Pnt::new_real3(tc[0], tc[1], tc[2]);
        let mut best_face = 0_usize;
        let mut best_dist = f64::MAX;
        for (fi, surface) in step_surfaces.iter().enumerate() {
            let projector = geom_api::ProjectPointOnSurf::new_pnt_handlegeomsurface_extalgo(
                &pt, surface, extrema::ExtAlgo::Grad,
            );
            if projector.is_done() && projector.nb_points() > 0 {
                let d = projector.lower_distance();
                if d < best_dist {
                    best_dist = d;
                    best_face = fi;
                }
            }
        }
        tri_surface_idx.push(best_face);
    }

    // For adjacent triangle pairs on the same STEP surface, compute dihedral angle
    let mut same_surface_angles: Vec<f64> = Vec::new();
    for tris in edge_to_triangles.values() {
        if tris.len() == 2 {
            let t1 = (tris[0] - 1) as usize; // convert 1-based to 0-based
            let t2 = (tris[1] - 1) as usize;
            if tri_surface_idx[t1] == tri_surface_idx[t2] {
                let n1 = &tri_normals[t1];
                let n2 = &tri_normals[t2];
                // Skip degenerate triangles
                let len1 = n1[0] * n1[0] + n1[1] * n1[1] + n1[2] * n1[2];
                let len2 = n2[0] * n2[0] + n2[1] * n2[1] + n2[2] * n2[2];
                if len1 < 0.5 || len2 < 0.5 {
                    continue;
                }
                let dot = n1[0] * n2[0] + n1[1] * n2[1] + n1[2] * n2[2];
                let angle_rad = dot.clamp(-1.0, 1.0).acos();
                let angle_deg = angle_rad.to_degrees();
                same_surface_angles.push(angle_deg);
            }
        }
    }

    let (ang_max, ang_avg, ang_count) = if same_surface_angles.is_empty() {
        (0.0_f64, 0.0_f64, 0_usize)
    } else {
        let max = same_surface_angles.iter().cloned().fold(0.0_f64, f64::max);
        let avg = same_surface_angles.iter().sum::<f64>() / same_surface_angles.len() as f64;
        (max, avg, same_surface_angles.len())
    };

    let vtx_avg_dist = vtx_dist_sum / num_nodes as f64;
    let ctr_avg_dist = if ctr_count > 0 {
        ctr_dist_sum / ctr_count as f64
    } else {
        0.0
    };

    let worst_pt = tri.node(vtx_max_idx);
    eprintln!(
        "Worst vertex #{}: ({:.6}, {:.6}, {:.6})",
        vtx_max_idx,
        worst_pt.x(),
        worst_pt.y(),
        worst_pt.z()
    );
    eprintln!("Vertex   avg: {:.10} mm  max: {:.10} mm", vtx_avg_dist, vtx_max_dist);
    eprintln!("Centroid avg: {:.10} mm  max: {:.10} mm", ctr_avg_dist, ctr_max_dist);
    eprintln!(
        "Angular  avg: {:.4} deg  max: {:.4} deg  ({} same-surface edges)",
        ang_avg, ang_max, ang_count
    );

    // Print tab-separated values to stdout for scripting:
    // vtx_max  vtx_avg  ctr_max  ctr_avg  max_dimension  ang_max  ang_avg  step_area  step_vol  stl_area  stl_vol
    println!(
        "{:.6e}\t{:.6e}\t{:.6e}\t{:.6e}\t{:.6}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.4}\t{:.4}",
        vtx_max_dist, vtx_avg_dist, ctr_max_dist, ctr_avg_dist, max_dimension, ang_max, ang_avg,
        step_area, step_volume, stl_area, stl_volume
    );
}
