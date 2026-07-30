//! Geometric comparison against a reference STEP — and the trap that made three
//! of these gates read zero for every input, right and wrong alike.
//!
//! # Never measure a distance against a solid
//!
//! `BRepExtrema_DistShapeShape` computes the **minimum** distance between two
//! shapes, and it reports **zero whenever one lies inside the other**. Point
//! against a solid therefore answers "is it inside?", not "how far is it from
//! the surface?", and shape against shape answers "do they overlap?", not "do
//! they agree?". Both are useless as accuracy gates, and neither announces
//! itself: a wrong answer and a perfect one both print `0.000000e0`.
//!
//! Measured 2026-07-29, all three found at once:
//!
//! * stage 1's vertex check passed a mesh written 1000× too small — a 25 mm
//!   sphere spanning 0.025 mm sat wholly inside its own reference and reported
//!   a distance of exactly 0. Its centroid check was worse than useless: an
//!   inscribed tessellation puts every facet centroid *inside* the solid, which
//!   is precisely the sag the check exists to measure, so it read 0 on convex
//!   geometry everywhere. One of the two recorded sag truth values was 4.3×
//!   understated as a result.
//! * stage 3.4 sampled a point per face against the reference solid, same trap.
//! * stages 3.6 and 4.1 compared our solid to the reference solid and gated on
//!   the result. Two solids that overlap at all are at distance zero, so those
//!   gates never fired. Measured on a decomposition with 9 surfaces against a
//!   reference 26 and *no planes at all*: both printed `0.000000e0`, exactly as
//!   for a perfect model.
//!
//! # What to use instead
//!
//! * A point against a boundary: [`face_compound`], then the point-distance
//!   helpers in [`crate::stage1`]. A compound of faces has no inside.
//! * Two shapes against each other: [`boundary_deviation`]. Note that swapping
//!   in face compounds is **not** sufficient here — the minimum distance
//!   between two boundaries that touch anywhere is still zero. It samples the
//!   reference's surface and measures each point to ours, which is also the one
//!   direction nothing else checks.
//!
//! If you add a comparison, ask what your metric returns for a *grossly* wrong
//! input. If that is indistinguishable from a right one, the metric is wrong.

use opencascade_sys::{b_rep, b_rep_mesh, top_abs, top_exp, top_loc, topo_ds, OwnedPtr};

/// A compound of every face of `shape` — the correct target for a
/// point-to-boundary distance.
///
/// A solid is not: see the module docs. Stage 3 builds edge and vertex
/// compounds for its own compares for the same reason.
pub fn face_compound(shape: &topo_ds::Shape) -> OwnedPtr<topo_ds::Compound> {
    let builder = topo_ds::Builder::new();
    let mut compound = topo_ds::Compound::new();
    builder.make_compound(&mut compound);
    let mut explorer = top_exp::Explorer::new_shape_shapeenum2(
        shape,
        top_abs::ShapeEnum::Face,
        top_abs::ShapeEnum::Shape,
    );
    while explorer.more() {
        builder.add(compound.as_shape_mut(), explorer.value());
        explorer.next();
    }
    compound
}

/// Triangulation deflection used by [`boundary_deviation`], as a fraction of
/// the tolerance being tested against.
///
/// The deviation is sampled at triangulation nodes, so its resolution is the
/// deflection it meshes at. A tenth of the tolerance keeps the discretization
/// an order below the quantity being gated. This is a numerical-accuracy choice
/// for a diagnostic, not a tuned threshold: it trades runtime for resolution and
/// changes no decision except by making the measurement sharper. It also tracks
/// `--surface-tolerance`, so tightening that in Phase 3 sharpens this for free.
const DEVIATION_MESH_FRACTION: f64 = 0.1;

/// Sample points: every triangulation node of every face of `shape`.
///
/// The shape must already be meshed. Node coordinates are in the face's local
/// frame, so the face location has to be applied — see `viz.rs`, which reads
/// triangulations the same way.
fn triangulation_nodes(shape: &topo_ds::Shape) -> Vec<[f64; 3]> {
    let mut points = Vec::new();
    let mut explorer = top_exp::Explorer::new_shape_shapeenum2(
        shape,
        top_abs::ShapeEnum::Face,
        top_abs::ShapeEnum::Shape,
    );
    while explorer.more() {
        let face = topo_ds::face_shape(explorer.value());
        let mut location = top_loc::Location::new();
        let tri_handle = b_rep::Tool::triangulation(face, &mut location, 0);
        if tri_handle.is_null() {
            explorer.next();
            continue;
        }
        let tri = tri_handle.get();
        let has_location = !location.is_identity();
        for i in 1..=tri.nb_nodes() {
            let node = tri.node(i);
            let (mut x, mut y, mut z) = (node.x(), node.y(), node.z());
            if has_location {
                location
                    .transformation()
                    .transforms_real3(&mut x, &mut y, &mut z);
            }
            points.push([x, y, z]);
        }
        explorer.next();
    }
    points
}

/// How far the reference's boundary strays from ours: the maximum distance from
/// any point on the reference surface to our nearest surface.
///
/// This is the direction nothing else in the pipeline measures. Stage 1 and
/// stage 3.4 both go ours-to-reference, which cannot notice a reference feature
/// we failed to reproduce — our points can all sit perfectly on the reference
/// while whole regions of the reference have nothing of ours near them. Running
/// it the other way catches exactly that.
///
/// Built from sampled points against a face compound rather than from
/// `DistShapeShape` on the two shapes, for the reason in the module docs: a
/// minimum distance between overlapping shapes is zero, and between two
/// boundaries that touch anywhere it is also zero. `BRepExtrema_ShapeProximity`
/// was tried first and rejected — its "minimal thickness of a layer containing
/// both shapes" reported 13.2 mm for a hemisphere whose volume matched its
/// reference to 1.2e-7, so whatever it measures, it is not this.
///
/// Returns `None` if the reference cannot be sampled, which callers must treat
/// as "not measured" rather than as agreement.
pub fn boundary_deviation(
    ours: &topo_ds::Shape,
    reference: &topo_ds::Shape,
    tolerance: f64,
) -> Option<f64> {
    let deflection = (tolerance * DEVIATION_MESH_FRACTION).max(1.0e-6);
    let _mesh = b_rep_mesh::IncrementalMesh::new_shape_real_bool_real_bool(
        reference,
        deflection,
        false,
        0.5, // ~28.6 degrees, matching the corpus generator's tessellation angle
        false,
    );
    let samples = triangulation_nodes(reference);
    if samples.is_empty() {
        return None;
    }
    let ours_faces = face_compound(ours);
    let (max_dist, _) = crate::stage1::max_distance_to_shape(&samples, ours_faces.as_shape());
    Some(max_dist)
}
