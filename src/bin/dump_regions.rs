//! Debug dump for stage 2.1/2.2 segmentation (scratch tool, not shipped).

use brepper::config;
use brepper::stage1::ConnectedMesh;
use brepper::stage2::edges::{classify_edges, compute_jet_field, EdgeClass};
use brepper::stage2::regions::extract_regions;
use brepper::stage1;

fn centroid(mesh: &ConnectedMesh, fi: usize) -> [f64; 3] {
    let f = &mesh.faces[fi];
    let vc = f.vertex_count as usize;
    let mut c = [0.0; 3];
    for &vi in &f.vertex_indices[..vc] {
        let v = &mesh.vertices[vi];
        c[0] += v.x;
        c[1] += v.y;
        c[2] += v.z;
    }
    let n = vc as f64;
    [c[0] / n, c[1] / n, c[2] / n]
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let path = &args[1];
    let mode = args.get(2).map(|s| s.as_str()).unwrap_or("summary");
    let config = config::parse_config_from(["brepper", path, "-q"]).unwrap();
    let mesh = stage1::stage1(&config).unwrap();
    let ec = classify_edges(&mesh, config.angular_tolerance_rad);
    let jet = compute_jet_field(&mesh, &ec);
    let (face_region, regions) = extract_regions(&mesh, &ec);
    eprintln!(
        "{} faces, {} regions: sharp={} interior={}",
        mesh.faces.len(),
        regions.len(),
        ec.sharp_count,
        ec.interior_count,
    );

    match mode {
        // Regions sorted by size, with bounding boxes of face centroids.
        "summary" => {
            let mut order: Vec<usize> = (0..regions.len()).collect();
            order.sort_by_key(|&r| std::cmp::Reverse(regions[r].faces.len()));
            for &r in order.iter().take(80) {
                let reg = &regions[r];
                let mut min = [f64::INFINITY; 3];
                let mut max = [f64::NEG_INFINITY; 3];
                for &fi in &reg.faces {
                    let c = centroid(&mesh, fi);
                    for k in 0..3 {
                        min[k] = min[k].min(c[k]);
                        max[k] = max[k].max(c[k]);
                    }
                }
                println!(
                    "region {r:4} faces={:5} ctr=[{:7.3},{:7.3},{:7.3}] ext=[{:6.3},{:6.3},{:6.3}]",
                    reg.faces.len(),
                    (min[0] + max[0]) / 2.0,
                    (min[1] + max[1]) / 2.0,
                    (min[2] + max[2]) / 2.0,
                    max[0] - min[0],
                    max[1] - min[1],
                    max[2] - min[2]
                );
            }
        }
        // Centroids of single-face regions (locate fragmented features).
        "singletons" => {
            for (r, reg) in regions.iter().enumerate() {
                if reg.faces.len() == 1 {
                    let c = centroid(&mesh, reg.faces[0]);
                    println!(
                        "region {r:4} face={:5} ctr=[{:7.3},{:7.3},{:7.3}]",
                        reg.faces[0], c[0], c[1], c[2]
                    );
                }
            }
        }
        // Region adjacency: edge counts by class between region pairs.
        "adj" => {
            let mut pairs: std::collections::BTreeMap<(usize, usize), [usize; 2]> =
                std::collections::BTreeMap::new();
            for fi in 0..mesh.faces.len() {
                let f = &mesh.faces[fi];
                let vc = f.vertex_count as usize;
                for e in 0..vc {
                    let ni = f.neighbors[e];
                    if ni < 0 {
                        continue;
                    }
                    let (ra, rb) = (face_region[fi], face_region[ni as usize]);
                    if ra == rb {
                        continue;
                    }
                    let key = (ra.min(rb), ra.max(rb));
                    let entry = pairs.entry(key).or_default();
                    match ec.class[fi][e] {
                        EdgeClass::Interior => entry[0] += 1,
                        EdgeClass::Sharp => entry[1] += 1,
                    }
                }
            }
            for ((ra, rb), [interior, sharp]) in pairs {
                println!("regions {ra:4} <-> {rb:4}: interior={interior:4} sharp={sharp:4}");
            }
        }
        // probe <x> <y> <z> [radius]: dump faces near a point with their
        // descriptors and per-edge kappa/class.
        "probe" => {
            let p: [f64; 3] = [
                args[3].parse().unwrap(),
                args[4].parse().unwrap(),
                args[5].parse().unwrap(),
            ];
            let radius: f64 = args.get(6).map(|s| s.parse().unwrap()).unwrap_or(0.6);
            let mut todo: Vec<usize> = (0..mesh.faces.len())
                .filter(|&fi| {
                    let c = centroid(&mesh, fi);
                    (c[0] - p[0]).powi(2) + (c[1] - p[1]).powi(2) + (c[2] - p[2]).powi(2)
                        < radius * radius
                })
                .collect();
            todo.sort_unstable();
            for fi in todo {
                let c = centroid(&mesh, fi);
                let f = &mesh.faces[fi];
                println!(
                    "face {fi} region={:3} vc={} ctr=[{:.3},{:.3},{:.3}] sigma=[{:+.4},{:+.4}]",
                    face_region[fi],
                    f.vertex_count,
                    c[0],
                    c[1],
                    c[2],
                    jet[fi].k1,
                    jet[fi].k2
                );
                let vc = f.vertex_count as usize;
                for e in 0..vc {
                    let ni = f.neighbors[e];
                    if ni < 0 {
                        continue;
                    }
                    let ni = ni as usize;
                    println!(
                        "   e{e} -> f{ni:5} region={:3} kappa={:+.4} class={:?}",
                        face_region[ni],
                        ec.kappa[fi][e],
                        ec.class[fi][e]
                    );
                }
            }
        }
        // leaks: Interior edges whose face descriptors differ a lot — true
        // junctions that the classifier missed. Prints midpoint + values.
        "leaks" => {
            let mut count = 0;
            for fi in 0..mesh.faces.len() {
                let f = &mesh.faces[fi];
                let vc = f.vertex_count as usize;
                for e in 0..vc {
                    let ni = f.neighbors[e];
                    if ni < 0 || (ni as usize) < fi {
                        continue;
                    }
                    let ni = ni as usize;
                    if ec.class[fi][e] != EdgeClass::Interior {
                        continue;
                    }
                    let jmax = (jet[fi].k1 - jet[ni].k1).abs();
                    let jmin = (jet[fi].k2 - jet[ni].k2).abs();
                    let jump = jmax.max(jmin);
                    if jump > 0.05 {
                        count += 1;
                        let ca = centroid(&mesh, fi);
                        let cb = centroid(&mesh, ni);
                        println!(
                            "LEAK f{fi}-f{ni} mid=[{:7.3},{:7.3},{:7.3}] kappa={:+.4} sigmaA=[{:+.4},{:+.4}] sigmaB=[{:+.4},{:+.4}]",
                            (ca[0] + cb[0]) / 2.0,
                            (ca[1] + cb[1]) / 2.0,
                            (ca[2] + cb[2]) / 2.0,
                            ec.kappa[fi][e],
                            jet[fi].k1,
                            jet[fi].k2,
                            jet[ni].k1,
                            jet[ni].k2
                        );
                    }
                }
            }
            eprintln!("{count} leak edges (Interior with descriptor jump > 0.05)");
        }
        // geom <fi>...: dump vertex positions + normals for given faces.
        "geom" => {
            for arg in &args[3..] {
                let fi: usize = arg.parse().unwrap();
                let f = &mesh.faces[fi];
                let vc = f.vertex_count as usize;
                println!(
                    "face {fi} vc={vc} region={} normal={:?}",
                    face_region[fi],
                    f.normal.map(|n| [n[0], n[1], n[2]])
                );
                for (i, &vi) in f.vertex_indices[..vc].iter().enumerate() {
                    let v = &mesh.vertices[vi];
                    println!("   v{i} (idx {vi:5}) = [{:9.5}, {:9.5}, {:9.5}]", v.x, v.y, v.z);
                }
            }
        }
        // rdesc <region_id>: descriptor histogram of a region's faces.
        "rdesc" => {
            let rid: usize = args[3].parse().unwrap();
            let mut hist: std::collections::BTreeMap<(i64, i64), usize> =
                std::collections::BTreeMap::new();
            for &fi in &regions[rid].faces {
                let kmax = (jet[fi].k1 * 50.0).round() as i64;
                let kmin = (jet[fi].k2 * 50.0).round() as i64;
                *hist.entry((kmax, kmin)).or_default() += 1;
            }
            for ((kmax, kmin), count) in hist {
                println!(
                    "  sigma=[{:+.3},{:+.3}]: {} faces",
                    kmax as f64 / 50.0,
                    kmin as f64 / 50.0,
                    count
                );
            }
        }
        // jet <fi>: replicate the jet fit for one face with diagnostics.
        "jet" => {
            let fi: usize = args[3].parse().unwrap();
            let f = &mesh.faces[fi];
            let vc = f.vertex_count as usize;
            let mut pts: Vec<([f64; 3], f64)> = Vec::new();
            let mut seen = std::collections::HashSet::new();
            for &vi in &f.vertex_indices[..vc] {
                if seen.insert(vi) {
                    let v = &mesh.vertices[vi];
                    pts.push(([v.x, v.y, v.z], 50.0));
                }
            }
            let mut ring1 = Vec::new();
            for e in 0..vc {
                let ni = f.neighbors[e];
                if ni >= 0 {
                    ring1.push(ni as usize);
                }
            }
            let push = |idx: usize, pts: &mut Vec<([f64; 3], f64)>, seen: &mut std::collections::HashSet<usize>| {
                let ff = &mesh.faces[idx];
                let fvc = ff.vertex_count as usize;
                for &vi in &ff.vertex_indices[..fvc] {
                    if seen.insert(vi) {
                        let v = &mesh.vertices[vi];
                        pts.push(([v.x, v.y, v.z], 1.0));
                    }
                }
            };
            for &ni in &ring1 {
                push(ni, &mut pts, &mut seen);
            }
            if pts.len() < 12 {
                for &ni in &ring1 {
                    let ff = &mesh.faces[ni];
                    let fvc = ff.vertex_count as usize;
                    for e2 in 0..fvc {
                        let n2 = ff.neighbors[e2];
                        if n2 >= 0 {
                            push(n2 as usize, &mut pts, &mut seen);
                        }
                    }
                }
            }
            println!("face {fi}: {} jet points (ring1 = {:?})", pts.len(), ring1);
            for (p, w) in &pts {
                println!("   [{:9.4},{:9.4},{:9.4}] w={}", p[0], p[1], p[2], w);
            }
            let c = centroid(&mesh, fi);
            let n = f.normal.unwrap();
            println!(
                "descriptor: sigma=[{:+.4},{:+.4}] centroid=[{:.3},{:.3},{:.3}] normal=[{:.4},{:.4},{:.4}]",
                jet[fi].k1, jet[fi].k2, c[0], c[1], c[2], n[0], n[1], n[2]
            );
        }
        _ => eprintln!("unknown mode {mode}"),
    }
}
