//! Debug dump for stage 2 (scratch tool, not shipped).
//!
//! `dump_regions <model.stl> <mode> [args]`
//!
//! Stage 2.1/2.2 segmentation: `summary`, `singletons`, `adj`, `leaks`,
//! `rdesc <region>`, `geom <face>...`, `probe <x> <y> <z> [radius]`, `jet
//! <face>`.
//!
//! Stage 2.3 classification (runs the whole of stage 2):
//! - `surfaces` — every committed surface, largest area first, with its
//!   parameters, to line up against a reference STEP's surface list.
//! - `at <x> <y> <z> [radius]` — which surfaces own the faces near a point;
//!   for tracing a stage-3 compare failure back to stage 2.
//! - `owner <face>...` — the surface owning a face and how well its
//!   neighbours' surfaces explain it (the stage-2.4 relabel criterion).
//! - `locus [a b]` — per adjacent surface pair, what the stage-2.4 junction
//!   locus concludes (tangent/transversal, locus built or not, faces on the
//!   wrong side of it); problems first. With two surface indices, the
//!   per-face detail for that one pair.

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

/// Area of a mesh face (triangle, or a quad as two triangles).
fn tri_area(mesh: &ConnectedMesh, fi: usize) -> f64 {
    let f = &mesh.faces[fi];
    let p = |k: usize| {
        let v = &mesh.vertices[f.vertex_indices[k]];
        [v.x, v.y, v.z]
    };
    let tri = |a: [f64; 3], b: [f64; 3], c: [f64; 3]| {
        let u = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
        let w = [c[0] - a[0], c[1] - a[1], c[2] - a[2]];
        let x = [u[1] * w[2] - u[2] * w[1], u[2] * w[0] - u[0] * w[2], u[0] * w[1] - u[1] * w[0]];
        0.5 * (x[0] * x[0] + x[1] * x[1] + x[2] * x[2]).sqrt()
    };
    let mut a = tri(p(0), p(1), p(2));
    if f.vertex_count == 4 {
        a += tri(p(0), p(2), p(3));
    }
    a
}

/// Distance from a point to a selected surface (analytic types only).
fn dist_to(
    out: &brepper::stage2::Stage2Output,
    sel: &brepper::stage2::SelectedSurface,
    q: [f64; 3],
) -> f64 {
    use brepper::stage2::SelectedSurface as S;
    match sel {
        S::Planar(i) => {
            let h = &out.planar_hypotheses[*i];
            (h.normal[0] * q[0] + h.normal[1] * q[1] + h.normal[2] * q[2] - h.distance).abs()
        }
        S::Conical(i) => {
            let h = &out.conical_hypotheses[*i];
            let d = [q[0] - h.apex[0], q[1] - h.apex[1], q[2] - h.apex[2]];
            let a = h.axis_direction;
            let t = d[0] * a[0] + d[1] * a[1] + d[2] * a[2];
            let perp = [d[0] - t * a[0], d[1] - t * a[1], d[2] - t * a[2]];
            let r = (perp[0].powi(2) + perp[1].powi(2) + perp[2].powi(2)).sqrt();
            if h.half_angle == 0.0 {
                (r - h.radius_at_apex).abs()
            } else {
                (r * h.half_angle.cos() - t * h.half_angle.sin()).abs()
            }
        }
        S::Toroidal(i) => {
            let h = &out.toroidal_hypotheses[*i];
            let d = [q[0] - h.center[0], q[1] - h.center[1], q[2] - h.center[2]];
            let a = h.axis_direction;
            let t = d[0] * a[0] + d[1] * a[1] + d[2] * a[2];
            let perp = [d[0] - t * a[0], d[1] - t * a[1], d[2] - t * a[2]];
            let r = (perp[0].powi(2) + perp[1].powi(2) + perp[2].powi(2)).sqrt();
            (((r - h.major_radius).powi(2) + t * t).sqrt() - h.minor_radius).abs()
        }
        S::Spherical(i) => {
            let h = &out.spherical_hypotheses[*i];
            let d = [q[0] - h.center[0], q[1] - h.center[1], q[2] - h.center[2]];
            ((d[0].powi(2) + d[1].powi(2) + d[2].powi(2)).sqrt() - h.radius).abs()
        }
        S::BSpline(_) => f64::NAN,
    }
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
        // What the stage-2.4 boundary locus concludes about each adjacent
        // pair of committed surfaces: tangent or transversal, whether a
        // junction locus exists, and how many faces sit on the wrong side of
        // it. Problems first. Read-only — it does not run snap.
        "locus" => {
            let out = brepper::stage2::stage2(&config, mesh, None).unwrap();
            let lines = match (args.get(3), args.get(4)) {
                (Some(a), Some(b)) => brepper::stage2::locus_pair_detail(
                    &out,
                    &config,
                    a.parse().unwrap(),
                    b.parse().unwrap(),
                ),
                _ => brepper::stage2::locus_report(&out, &config),
            };
            for line in lines {
                println!("{line}");
            }
        }
        // Every surface stage 2 committed, largest first: type, extent, and
        // the parameters to line up against a reference STEP's surface list.
        "surfaces" => {
            let out = brepper::stage2::stage2(&config, mesh, None).unwrap();
            let mut rows: Vec<(usize, f64, String)> = Vec::new();
            for s in &out.selected_surfaces {
                use brepper::stage2::SelectedSurface as S;
                let (faces, desc) = match s {
                    S::Planar(i) => {
                        let h = &out.planar_hypotheses[*i];
                        let n = h.normal;
                        (&h.faces, format!("plane   n=[{:+.3},{:+.3},{:+.3}] d={:+.4}", n[0], n[1], n[2], h.distance))
                    }
                    S::Conical(i) => {
                        let h = &out.conical_hypotheses[*i];
                        let a = h.axis_direction;
                        if h.is_cylinder() {
                            (&h.faces, format!("cyl     r={:.4} axis=[{:+.3},{:+.3},{:+.3}]", h.radius_at_apex, a[0], a[1], a[2]))
                        } else {
                            (&h.faces, format!("cone    ha={:.4} axis=[{:+.3},{:+.3},{:+.3}]", h.half_angle, a[0], a[1], a[2]))
                        }
                    }
                    S::Spherical(i) => {
                        let h = &out.spherical_hypotheses[*i];
                        (&h.faces, format!("sphere  r={:.4}", h.radius))
                    }
                    S::Toroidal(i) => {
                        let h = &out.toroidal_hypotheses[*i];
                        let a = h.axis_direction;
                        (&h.faces, format!("torus   R={:.4} r={:.4} axis=[{:+.3},{:+.3},{:+.3}]", h.major_radius, h.minor_radius, a[0], a[1], a[2]))
                    }
                    S::BSpline(i) => {
                        let h = &out.bspline_hypotheses[*i];
                        (&h.faces, "bspline".to_string())
                    }
                };
                let area: f64 = faces.iter().map(|&fi| tri_area(&out.mesh, fi)).sum();
                let ids = if faces.len() <= 40 { format!("  {faces:?}") } else { String::new() };
                let mut vset: Vec<usize> = faces.iter().flat_map(|&fi| {
                    let f = &out.mesh.faces[fi];
                    f.vertex_indices[..f.vertex_count as usize].to_vec()
                }).collect();
                vset.sort_unstable();
                vset.dedup();
                let mut resid: Vec<f64> = vset.iter().map(|&vi| {
                    let v = &out.mesh.vertices[vi];
                    dist_to(&out, s, [v.x, v.y, v.z])
                }).collect();
                // A B-spline has no closed-form distance here; report nothing.
                resid.retain(|d| d.is_finite());
                if resid.is_empty() {
                    resid.push(f64::NAN);
                }
                resid.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
                let q = |f: f64| resid[((resid.len() - 1) as f64 * f) as usize];
                let own = match s {
                    S::Planar(i) => out.planar_hypotheses[*i].error_max.abs(),
                    S::Conical(i) => out.conical_hypotheses[*i].error_max,
                    S::Spherical(i) => out.spherical_hypotheses[*i].error_max,
                    S::Toroidal(i) => out.toroidal_hypotheses[*i].error_max,
                    S::BSpline(_) => f64::NAN,
                };
                let stats = format!(
                    " vresid p10 {:.1e} med {:.1e} max {:.1e} | hyp.error_max {own:.1e}",
                    q(0.1), q(0.5), resid[resid.len() - 1]
                );
                rows.push((faces.len(), area, format!("{desc}{stats}{ids}")));
            }
            rows.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap());
            for (nf, area, desc) in &rows {
                println!("{nf:6} faces  area {area:10.3}  {desc}");
            }
            println!("{} surfaces", rows.len());
        }
        // at <x> <y> <z> [radius]: which committed surfaces own the faces near
        // a point — for tracing a stage-3 compare failure back to stage 2.
        "at" => {
            let p: [f64; 3] = [
                args[3].parse().unwrap(),
                args[4].parse().unwrap(),
                args[5].parse().unwrap(),
            ];
            let radius: f64 = args.get(6).map(|s| s.parse().unwrap()).unwrap_or(0.3);
            let out = brepper::stage2::stage2(&config, mesh, None).unwrap();
            use brepper::stage2::SelectedSurface as S;
            for (k, s) in out.selected_surfaces.iter().enumerate() {
                let (faces, desc) = match s {
                    S::Planar(i) => {
                        let h = &out.planar_hypotheses[*i];
                        (&h.faces, format!("plane n=[{:+.3},{:+.3},{:+.3}] d={:+.4}", h.normal[0], h.normal[1], h.normal[2], h.distance))
                    }
                    S::Conical(i) => {
                        let h = &out.conical_hypotheses[*i];
                        (&h.faces, format!("cyl/cone r={:.4} ha={:.4}", h.radius_at_apex, h.half_angle))
                    }
                    S::Spherical(i) => (&out.spherical_hypotheses[*i].faces, "sphere".to_string()),
                    S::Toroidal(i) => {
                        let h = &out.toroidal_hypotheses[*i];
                        (&h.faces, format!("torus R={:.4} r={:.4}", h.major_radius, h.minor_radius))
                    }
                    S::BSpline(i) => (&out.bspline_hypotheses[*i].faces, "bspline".to_string()),
                };
                let near: Vec<usize> = faces
                    .iter()
                    .copied()
                    .filter(|&fi| {
                        let c = centroid(&out.mesh, fi);
                        (c[0] - p[0]).powi(2) + (c[1] - p[1]).powi(2) + (c[2] - p[2]).powi(2)
                            < radius * radius
                    })
                    .collect();
                if !near.is_empty() {
                    println!(
                        "surface {k:3} ({:5} faces) {desc}  — {} nearby: {:?}",
                        faces.len(),
                        near.len(),
                        &near[..near.len().min(8)]
                    );
                }
            }
        }
        // owner <face>...: which surface owns each face, and how well its
        // neighbours' surfaces explain it (the stage-2.4 relabel criterion).
        "owner" => {
            let out = brepper::stage2::stage2(&config, mesh, None).unwrap();
            let mut f2s = vec![usize::MAX; out.mesh.faces.len()];
            use brepper::stage2::SelectedSurface as S;
            let faces_of = |k: usize| -> &Vec<usize> {
                match &out.selected_surfaces[k] {
                    S::Planar(i) => &out.planar_hypotheses[*i].faces,
                    S::Conical(i) => &out.conical_hypotheses[*i].faces,
                    S::Spherical(i) => &out.spherical_hypotheses[*i].faces,
                    S::Toroidal(i) => &out.toroidal_hypotheses[*i].faces,
                    S::BSpline(i) => &out.bspline_hypotheses[*i].faces,
                }
            };
            for k in 0..out.selected_surfaces.len() {
                for &fi in faces_of(k) {
                    f2s[fi] = k;
                }
            }
            let kind = |k: usize| match &out.selected_surfaces[k] {
                S::Planar(i) => {
                    let h = &out.planar_hypotheses[*i];
                    format!("plane d={:+.5} n=[{:+.3},{:+.3},{:+.3}]", h.distance, h.normal[0], h.normal[1], h.normal[2])
                }
                S::Conical(i) => format!("cyl r={:.4}", out.conical_hypotheses[*i].radius_at_apex),
                S::Spherical(_) => "sphere".to_string(),
                S::Toroidal(i) => {
                    let h = &out.toroidal_hypotheses[*i];
                    format!("torus R={:.4} r={:.4}", h.major_radius, h.minor_radius)
                }
                S::BSpline(_) => "bspline".to_string(),
            };
            for arg in &args[3..] {
                let fi: usize = arg.parse().unwrap();
                let vdist = |k: usize| -> f64 {
                    let f = &out.mesh.faces[fi];
                    f.vertex_indices[..f.vertex_count as usize]
                        .iter()
                        .map(|&vi| {
                            let v = &out.mesh.vertices[vi];
                            dist_to(&out, &out.selected_surfaces[k], [v.x, v.y, v.z])
                        })
                        .fold(0.0_f64, f64::max)
                };
                println!("face {fi}: owner surface {} = {} (max vertex dist {:.3e})",
                    f2s[fi], kind(f2s[fi]), vdist(f2s[fi]));
                let f = &out.mesh.faces[fi];
                for &ni in &f.neighbors[..f.vertex_count as usize] {
                    if ni >= 0 {
                        let sj = f2s[ni as usize];
                        println!("   neighbour face {ni} -> surface {sj} = {} (this face's max vertex dist to it {:.3e})",
                            kind(sj), vdist(sj));
                    }
                }
            }
        }
        _ => eprintln!("unknown mode {mode}"),
    }
}
