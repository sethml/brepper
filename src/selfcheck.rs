//! Self-check: is the solid we are about to hand back the mesh we were given?
//!
//! Runs at the end of every conversion, with no reference involved. It is what
//! turns "exited 0" into a statement about the output. Before it existed a run
//! with no `--compare` could write a shell with no solid, a face of negative
//! area, or a solid 2 mm off its own input and report success — measured on the
//! wild corpus, 4 of 32 verified parts did.
//!
//! Three questions, each answered against the INPUT MESH:
//!
//! 1. **Is it a valid closed solid?** At least one solid, and `BRepCheck`
//!    passes.
//! 2. **Is every mesh vertex on the output?** Each vertex is measured against
//!    the trimmed faces of the solid near it.
//! 3. **Does every mesh facet have output under it?** This is the direction a
//!    vertex test cannot see: a face trimmed too small, or missing, can still
//!    leave every vertex on some neighbouring face. Each facet's centroid is
//!    measured against the trimmed output and against the untrimmed surface
//!    the facet was assigned to. Over a correctly trimmed face the two agree
//!    (both are the facet's chord sag); where the output is farther than the
//!    surface, there is no output under that facet.
//!
//! The faces measured are those of the solid being written, not the ones
//! stage 3.4 built: stage 3.6 re-trims some, and a check against the earlier
//! faces fails solids that are correct (see `self_check`).
//!
//! The only allowance is the B-rep's own: OCCT stores a tolerance on every
//! face, edge and vertex saying how far the geometry there may be trusted, and
//! a point is allowed `--vertex-tolerance` plus the largest such tolerance on
//! the face it is measured against. Nothing here is tuned.
//!
//! That allowance cannot be open-ended, or a solid could pass any test by
//! declaring itself loose (measured: one output carried a 4.7e5 mm tolerance).
//! So the largest declared tolerance (`tol_max`) is itself held to
//! `--surface-tolerance`, the bound the user already gives for how far the
//! output may sit from the mesh's facets, and it is reported, because how
//! loose a solid declares itself is a quality number in its own right.

use crate::config::Config;
use crate::stage1::ConnectedMesh;
use crate::stage3::Stage3Output;
use opencascade_sys::{
    b_rep, b_rep_bnd_lib, b_rep_check, b_rep_g_prop, bnd, extrema, g_prop, geom, geom_api, gp,
    top_abs, top_exp, topo_ds, OwnedPtr,
};

/// What the self-check found.
#[derive(Debug, Clone)]
pub struct SelfCheck {
    pub solids: usize,
    pub faces: usize,
    /// At least one solid, `BRepCheck` valid, no face lost in assembly.
    pub valid: bool,
    /// Largest distance from a mesh vertex to a trimmed face it belongs to.
    pub vertex_max: f64,
    /// Mesh vertices farther than the allowance from one of their faces.
    pub vertices_off: usize,
    /// Mesh facets whose own output face does not reach them.
    pub facets_uncovered: usize,
    /// Largest amount by which a facet centroid's distance to its trimmed face
    /// exceeds its distance to the untrimmed surface.
    pub cover_worst: f64,
    /// `(solid - mesh) / mesh`. Reported, not gated: a mesh legitimately
    /// encloses a little less than the solid it was tessellated from.
    pub volume_rel: f64,
    /// Largest OCCT tolerance on any face, edge or vertex of the output.
    pub tolerance_max: f64,
    /// `tolerance_max` exceeds `--surface-tolerance`.
    pub too_loose: bool,
    /// The vertex and facet measurements were made. They are skipped on an
    /// invalid solid: its status is already decided, and a broken face can make
    /// every query against it slow (measured: 32 s on one such output).
    pub measured: bool,
    /// Why `valid` is false, when it is.
    pub problems: Vec<String>,
    /// The worst few offending vertices and facets, described for a person.
    pub examples: Vec<String>,
}

impl SelfCheck {
    /// `good`, `unfaithful` (valid solid, but not the mesh) or `invalid`.
    pub fn status(&self) -> &'static str {
        if !self.valid {
            "invalid"
        } else if self.vertices_off > 0 || self.facets_uncovered > 0 || self.too_loose {
            "unfaithful"
        } else {
            "good"
        }
    }

    /// One `key=value` line: read by people and by `tools/scoreboard.py`.
    pub fn line(&self) -> String {
        let fidelity = if self.measured {
            format!(
                "vtx_max={:.3e} vtx_off={} uncovered={} cover_worst={:.3e}",
                self.vertex_max, self.vertices_off, self.facets_uncovered, self.cover_worst
            )
        } else {
            "vtx_max=- vtx_off=- uncovered=- cover_worst=-".to_string()
        };
        format!(
            "Self-check: status={} solids={} faces={} valid={} {fidelity} vol_rel={:+.3e} \
             tol_max={:.3e}",
            self.status(),
            self.solids,
            self.faces,
            self.valid,
            self.volume_rel,
            self.tolerance_max,
        )
    }
}

/// Signed volume enclosed by the mesh (signed tetrahedra against the origin).
fn mesh_volume(mesh: &ConnectedMesh) -> f64 {
    let mut vol = 0.0;
    for f in &mesh.faces {
        let p = |k: usize| {
            let v = &mesh.vertices[f.vertex_indices[k]];
            [v.x, v.y, v.z]
        };
        for k in 1..(f.vertex_count as usize).saturating_sub(1) {
            let (a, b, c) = (p(0), p(k), p(k + 1));
            vol += (a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0])
                + a[2] * (b[0] * c[1] - b[1] * c[0]))
                / 6.0;
        }
    }
    vol
}

/// The largest OCCT tolerance on a face or on any edge or vertex bounding it.
fn face_tolerance(face: &topo_ds::Face) -> f64 {
    let mut tol = b_rep::Tool::tolerance_face(face);
    for kind in [top_abs::ShapeEnum::Edge, top_abs::ShapeEnum::Vertex] {
        let mut ex = top_exp::Explorer::new_shape_shapeenum2(
            face.as_shape(),
            kind,
            top_abs::ShapeEnum::Shape,
        );
        while ex.more() {
            let t = if kind == top_abs::ShapeEnum::Edge {
                b_rep::Tool::tolerance_edge(topo_ds::edge(ex.value()))
            } else {
                b_rep::Tool::tolerance_vertex(topo_ds::vertex(ex.value()))
            };
            tol = tol.max(t);
            ex.next();
        }
    }
    tol
}

/// One face of the solid about to be written, with what a point query needs.
struct OutFace {
    face: OwnedPtr<topo_ds::Face>,
    /// Largest OCCT tolerance on the face or its boundary, capped at
    /// `--surface-tolerance`: a face declaring more than that already fails
    /// the check as too loose, and must not also buy itself a wider allowance
    /// (or, through its grown bounding box, a query against every point).
    tol: f64,
    lo: [f64; 3],
    hi: [f64; 3],
}

/// Every face of `compound`, with its tolerance and bounding box, and the
/// largest tolerance any of them declares (before the cap).
fn out_faces(compound: &topo_ds::Shape, tol_cap: f64) -> (Vec<OutFace>, f64) {
    let mut declared_max = 0.0_f64;
    let mut out = Vec::new();
    let mut ex = top_exp::Explorer::new_shape_shapeenum2(
        compound,
        top_abs::ShapeEnum::Face,
        top_abs::ShapeEnum::Shape,
    );
    while ex.more() {
        let face = topo_ds::face_shape(ex.value()).to_owned();
        let mut bbox = bnd::Box::new();
        // No triangulation exists, so this bounds the geometry (for a B-spline,
        // its control hull). Too large is the safe direction: the box only
        // decides which faces a point is measured against.
        b_rep_bnd_lib::add_shape_box_bool(face.as_shape(), &mut bbox, false);
        let (lo, hi) = if bbox.is_void() {
            ([f64::MIN; 3], [f64::MAX; 3])
        } else {
            (
                [bbox.get_x_min(), bbox.get_y_min(), bbox.get_z_min()],
                [bbox.get_x_max(), bbox.get_y_max(), bbox.get_z_max()],
            )
        };
        let declared = face_tolerance(&face);
        declared_max = declared_max.max(declared);
        out.push(OutFace { tol: declared.min(tol_cap), face, lo, hi });
        ex.next();
    }
    (out, declared_max)
}

/// A point's closest approach to the output.
#[derive(Clone, Copy)]
struct Nearest {
    /// Distance to the nearest trimmed face among those whose box, grown by
    /// the query radius, contains the point. `f64::MAX` if there is none.
    distance: f64,
    /// That face's tolerance (the solid's largest, if there is none).
    tol: f64,
}

/// For each `(point, radius)`, the nearest output face within `radius`.
///
/// `radius` says how far from the point its face could legitimately be: a face
/// is only measured if its bounding box, grown by `radius` plus the face's own
/// tolerance, contains the point. That is what makes this cheap — a point is
/// measured against the one to three faces near it, not against the solid.
///
/// Safety: the same argument as [`crate::stage1::distances_to_shape`] —
/// `BRepExtrema_DistShapeShape` does not mutate its target, so read-only faces
/// may be shared across threads.
fn nearest_faces(tasks: &[([f64; 3], f64)], faces: &[OutFace], tol_max: f64) -> Vec<Nearest> {
    struct Shared<'a>(&'a [OutFace]);
    unsafe impl Send for Shared<'_> {}
    unsafe impl Sync for Shared<'_> {}
    let shared = Shared(faces);
    let one = |(p, radius): &([f64; 3], f64), faces: &[OutFace]| -> Nearest {
        let pt = gp::Pnt::new_real3(p[0], p[1], p[2]);
        let mut best = Nearest { distance: f64::MAX, tol: tol_max };
        for f in faces {
            let grow = radius + f.tol;
            if (0..3).any(|k| p[k] < f.lo[k] - grow || p[k] > f.hi[k] + grow) {
                continue;
            }
            let d = crate::stage1::min_distance_to_shape(&pt, f.face.as_shape());
            // Nearest relative to what each face allows, so a loose face next
            // to a tight one does not hide a point the tight one holds.
            if d - f.tol < best.distance - best.tol {
                best = Nearest { distance: d, tol: f.tol };
            }
        }
        best
    };

    let n = tasks.len();
    let threads = std::thread::available_parallelism()
        .map(|p| p.get())
        .unwrap_or(1)
        .min(8)
        .min(n.div_ceil(64).max(1));
    if threads <= 1 {
        return tasks.iter().map(|t| one(t, shared.0)).collect();
    }
    let chunk = n.div_ceil(threads);
    let mut out = vec![Nearest { distance: f64::MAX, tol: tol_max }; n];
    std::thread::scope(|s| {
        for (task_chunk, out_chunk) in tasks.chunks(chunk).zip(out.chunks_mut(chunk)) {
            let shared = &shared;
            s.spawn(move || {
                for (t, o) in task_chunk.iter().zip(out_chunk.iter_mut()) {
                    *o = one(t, shared.0);
                }
            });
        }
    });
    out
}

/// Check `compound` (the solids about to be written) against the mesh in
/// `output`.
pub fn self_check(
    config: &Config,
    output: &Stage3Output,
    compound: &topo_ds::Shape,
) -> SelfCheck {
    let stage2 = &output.stage2;
    let mesh = &stage2.mesh;
    let mut problems = Vec::new();

    // --- 1. A valid closed solid. ---
    if output.solids.is_empty() {
        problems.push("no solid was built".to_string());
    } else if !b_rep_check::Analyzer::new_shape_bool(compound, true).is_valid() {
        problems.push("BRepCheck reports the solid invalid".to_string());
    }

    let solid_volume = {
        let mut props = g_prop::GProps::new();
        b_rep_g_prop::volume_properties_shape_gprops_bool3(compound, &mut props, true, false, false);
        props.mass()
    };
    let mesh_vol = mesh_volume(mesh);
    let volume_rel = if mesh_vol.abs() > 0.0 {
        (solid_volume.abs() - mesh_vol.abs()) / mesh_vol.abs()
    } else {
        0.0
    };

    // The faces of the solid itself, not `output.faces`. Stage 3.6 does not
    // merely collect the faces stage 3.4 built: `ShapeFix_Solid` re-trims some
    // (measured on remesh/cone_tangent_planes_offset: a cone face built over
    // the complementary sector, 14 mm from its own facets, comes out of 3.6
    // correct). What is written is what must be checked.
    let (faces, declared_max) = out_faces(compound, config.surface_tolerance_mm);
    let tol_max = faces.iter().map(|f| f.tol).fold(0.0, f64::max);
    let vtol = config.vertex_tolerance_mm;
    let mut result = SelfCheck {
        solids: output.solids.len(),
        faces: faces.len(),
        valid: problems.is_empty(),
        vertex_max: 0.0,
        vertices_off: 0,
        facets_uncovered: 0,
        cover_worst: 0.0,
        volume_rel,
        tolerance_max: declared_max,
        too_loose: declared_max > config.surface_tolerance_mm,
        measured: false,
        problems,
        examples: Vec::new(),
    };
    if !result.valid || faces.is_empty() {
        return result;
    }
    result.measured = true;

    // --- 2. Every mesh vertex on the output. ---
    let tasks: Vec<([f64; 3], f64)> =
        mesh.vertices.iter().map(|v| ([v.x, v.y, v.z], vtol)).collect();
    let near = nearest_faces(&tasks, &faces, tol_max);
    let mut offenders: Vec<(f64, usize)> = Vec::new();
    for (vi, n) in near.iter().enumerate() {
        if n.distance > vtol + n.tol {
            offenders.push((n.distance, vi));
        } else {
            result.vertex_max = result.vertex_max.max(n.distance);
        }
    }
    result.vertices_off = offenders.len();
    if !offenders.is_empty() {
        // A vertex with no face inside its query radius has no distance yet;
        // measure the offenders against the whole solid so the report is a
        // number. (Few points, so the slow query is affordable.)
        for (d, vi) in offenders.iter_mut() {
            let v = &mesh.vertices[*vi];
            *d = crate::stage1::min_distance_to_shape(&gp::Pnt::new_real3(v.x, v.y, v.z), compound);
        }
        offenders.sort_by(|a, b| b.0.total_cmp(&a.0));
        result.vertex_max = result.vertex_max.max(offenders[0].0);
        for &(d, vi) in offenders.iter().take(3) {
            let v = &mesh.vertices[vi];
            result.examples.push(format!(
                "vertex {vi} ({:.4}, {:.4}, {:.4}) is {d:.3e} mm from the output",
                v.x, v.y, v.z
            ));
        }
    }

    // --- 3. Every mesh facet has output under it. ---
    // A facet's centroid is within the facet's own radius of its vertices,
    // which are on its face, so that radius is how far its face can be.
    let face_surface = crate::stage3::build_face_to_surface_map(stage2);
    let mut centroid_tasks: Vec<([f64; 3], f64)> = Vec::with_capacity(mesh.faces.len());
    let mut centroid_facet: Vec<usize> = Vec::with_capacity(mesh.faces.len());
    for (fi, f) in mesh.faces.iter().enumerate() {
        let n = f.vertex_count as usize;
        if n == 0 {
            continue;
        }
        let p = |k: usize| {
            let v = &mesh.vertices[f.vertex_indices[k]];
            [v.x, v.y, v.z]
        };
        let mut c = [0.0_f64; 3];
        for k in 0..n {
            let v = p(k);
            for i in 0..3 {
                c[i] += v[i] / n as f64;
            }
        }
        let radius = (0..n)
            .map(|k| {
                let v = p(k);
                ((v[0] - c[0]).powi(2) + (v[1] - c[1]).powi(2) + (v[2] - c[2]).powi(2)).sqrt()
            })
            .fold(0.0_f64, f64::max);
        centroid_tasks.push((c, radius));
        centroid_facet.push(fi);
    }
    let near = nearest_faces(&centroid_tasks, &faces, tol_max);
    // The untrimmed surface each facet was assigned to. Stage 3.4's faces are
    // good for this even where their trim is not: the surface is the fitted
    // one either way. Fetched once per surface, and only where needed — a
    // centroid already within the allowance of the output (every facet of a
    // planar face) is covered without asking.
    let have_surfaces = output.faces.len() == stage2.selected_surfaces.len();
    let mut surfaces: Vec<Option<OwnedPtr<geom::Handle_Geom_Surface>>> =
        (0..output.faces.len()).map(|_| None).collect();
    let mut uncovered: Vec<(f64, usize)> = Vec::new();
    for ((&(p, _), n), &fi) in centroid_tasks.iter().zip(&near).zip(&centroid_facet) {
        let allow = vtol + n.tol;
        if n.distance <= allow {
            continue;
        }
        let si = face_surface[fi];
        // How far the centroid is from its surface: the facet's chord sag,
        // which is not the output's doing.
        let sag = if have_surfaces && si != usize::MAX {
            let surface =
                surfaces[si].get_or_insert_with(|| b_rep::Tool::surface_face(&output.faces[si]));
            let proj = geom_api::ProjectPointOnSurf::new_pnt_handlegeomsurface_extalgo(
                &gp::Pnt::new_real3(p[0], p[1], p[2]),
                surface,
                extrema::ExtAlgo::Grad,
            );
            if proj.is_done() && proj.nb_points() > 0 {
                proj.lower_distance()
            } else {
                0.0
            }
        } else {
            0.0
        };
        if n.distance - sag > allow {
            uncovered.push((n.distance - sag, fi));
        }
    }
    result.facets_uncovered = uncovered.len();
    if !uncovered.is_empty() {
        uncovered.sort_by(|a, b| b.0.total_cmp(&a.0));
        // Same as for vertices: put a real distance on the worst few.
        for (excess, fi) in uncovered.iter_mut().take(3) {
            let k = centroid_facet.iter().position(|f| f == fi).unwrap();
            let c = centroid_tasks[k].0;
            if *excess > 1.0e300 {
                *excess = crate::stage1::min_distance_to_shape(
                    &gp::Pnt::new_real3(c[0], c[1], c[2]),
                    compound,
                );
            }
            result.examples.push(format!(
                "facet {fi} (centroid {:.4}, {:.4}, {:.4}) has no output within {:.3e} mm \
                 beyond its own sag",
                c[0], c[1], c[2], excess
            ));
        }
        result.cover_worst = uncovered
            .iter()
            .map(|u| u.0)
            .filter(|e| *e < 1.0e300)
            .fold(0.0, f64::max);
    }
    result
}
