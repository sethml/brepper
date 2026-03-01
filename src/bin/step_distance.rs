/// Compute maximum distance between STL mesh vertices/centroids and STEP model surfaces.
///
/// Reads an STL file to extract vertex positions and triangle centroids, reads a
/// STEP file to extract surfaces, and computes the minimum distance from each
/// point to any STEP surface. Reports max distance for both vertices and centroids.
///
/// Units: all distances are in mm. OCCT's STEPControl_Reader converts STEP
/// file units (typically meters for Onshape exports) to mm internally.
/// STL files have no unit metadata; coordinates are assumed to be in mm.
use opencascade_sys::{
    b_rep, extrema, geom_api, gp, message, rw_stl, step_control, top_abs, top_exp, topo_ds,
};
use std::env;
use std::process;

fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() != 3 {
        eprintln!("Usage: {} <input.stl> <input.step>", args[0]);
        process::exit(1);
    }
    let stl_path = &args[1];
    let step_path = &args[2];

    // Read STL → Poly_Triangulation
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
    let mut reader = step_control::Reader::new();
    reader.read_file_charptr(step_path);
    let progress2 = message::ProgressRange::new();
    reader.transfer_roots(&progress2);
    let step_shape = reader.one_shape();

    // Collect surfaces from all STEP faces
    let mut surfaces = Vec::new();
    let mut explorer = top_exp::Explorer::new_shape_shapeenum2(
        &step_shape,
        top_abs::ShapeEnum::Face,
        top_abs::ShapeEnum::Shape,
    );
    while explorer.more() {
        let shape_ref = explorer.value();
        let face = topo_ds::face_shape(shape_ref);
        surfaces.push(b_rep::Tool::surface_face(face));
        explorer.next();
    }
    eprintln!("STEP: {} faces", surfaces.len());

    if surfaces.is_empty() {
        eprintln!("Error: STEP file has no faces");
        process::exit(1);
    }

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

    let min_distance_to_surfaces = |pt: &gp::Pnt| -> f64 {
        let mut min_dist = f64::MAX;
        for surface in &surfaces {
            let projector = geom_api::ProjectPointOnSurf::new_pnt_handlegeomsurface_extalgo(
                pt,
                surface,
                extrema::ExtAlgo::Grad,
            );
            if projector.is_done() && projector.nb_points() > 0 {
                let d = projector.lower_distance();
                if d < min_dist {
                    min_dist = d;
                }
            }
        }
        min_dist
    };

    // Vertex distances
    let mut vtx_max_dist = 0.0_f64;
    let mut vtx_max_idx = 1_i32;
    let mut vtx_dist_sum = 0.0_f64;

    for i in 1..=num_nodes {
        let pt = tri.node(i);
        let min_dist = min_distance_to_surfaces(&pt);
        if min_dist < f64::MAX {
            vtx_dist_sum += min_dist;
            if min_dist > vtx_max_dist {
                vtx_max_dist = min_dist;
                vtx_max_idx = i;
            }
        }
    }

    // Centroid distances
    let mut ctr_max_dist = 0.0_f64;
    let mut _ctr_max_idx = 1_i32;
    let mut ctr_dist_sum = 0.0_f64;

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
        let min_dist = min_distance_to_surfaces(&centroid);
        if min_dist < f64::MAX {
            ctr_dist_sum += min_dist;
            if min_dist > ctr_max_dist {
                ctr_max_dist = min_dist;
                _ctr_max_idx = i;
            }
        }
    }

    let vtx_avg_dist = vtx_dist_sum / num_nodes as f64;
    let ctr_avg_dist = ctr_dist_sum / num_triangles as f64;

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

    // Print tab-separated values to stdout for scripting:
    // vtx_max  vtx_avg  ctr_max  ctr_avg  max_dimension
    println!(
        "{:.6e}\t{:.6e}\t{:.6e}\t{:.6e}\t{:.6}",
        vtx_max_dist, vtx_avg_dist, ctr_max_dist, ctr_avg_dist, max_dimension
    );
}
