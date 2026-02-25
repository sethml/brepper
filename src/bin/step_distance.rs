/// Compute maximum distance between STL mesh vertices and STEP model surfaces.
///
/// Reads an STL file to extract vertex positions, reads a STEP file to extract
/// surfaces, and computes the minimum distance from each STL vertex to any STEP
/// surface. Reports the maximum such distance (i.e., the worst-case vertex error).
use opencascade_sys::{
    b_rep, extrema, geom_api, message, rw_stl, step_control, top_abs, top_exp, topo_ds,
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

    // For each STL vertex, find minimum distance to any STEP surface
    let mut max_dist = 0.0_f64;
    let mut max_vertex_idx = 1_i32;
    let mut dist_sum = 0.0_f64;

    for i in 1..=num_nodes {
        let pt = tri.node(i);
        let mut min_dist = f64::MAX;

        for surface in &surfaces {
            let projector = geom_api::ProjectPointOnSurf::new_pnt_handlegeomsurface_extalgo(
                &pt,
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

        if min_dist < f64::MAX {
            dist_sum += min_dist;
            if min_dist > max_dist {
                max_dist = min_dist;
                max_vertex_idx = i;
            }
        }
    }

    let avg_dist = dist_sum / num_nodes as f64;
    let worst_pt = tri.node(max_vertex_idx);
    eprintln!(
        "Worst vertex #{}: ({:.6}, {:.6}, {:.6})",
        max_vertex_idx,
        worst_pt.x(),
        worst_pt.y(),
        worst_pt.z()
    );
    eprintln!("Average distance: {:.10}", avg_dist);
    eprintln!("Maximum distance: {:.10}", max_dist);

    // Print max distance to stdout for scripting
    println!("{:.10}", max_dist);
}
