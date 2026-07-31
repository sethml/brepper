//! Reference labeling: ground truth per mesh face, read back from the compare
//! call the pipeline already makes.
//!
//! **Diagnostics only.** Nothing downstream reads any of this; it is printed
//! under `--reference-report` (which requires `--compare`) and changes no
//! decision. Like [`crate::estimate`], its job is to make a claim measurable
//! before anything is built on it.
//!
//! # What it is for
//!
//! Phase 1 replaces the blind seed patch with a patch grown under validation,
//! initialized in closed form from facet normals — "a planar facet spanning two
//! rulings of a cylinder has its normal exactly perpendicular to the axis, so
//! two facets give the axis exactly". That is exact on strip/quad tessellation
//! and *not* exact on a general Delaunay triangle, and the plan says to measure
//! the error per provenance before relying on it. This module is that
//! measurement, plus the two others that fall out of the same labeling for
//! nearly nothing.
//!
//! # Why it lives here and not in a probe
//!
//! Stage 1's `--compare` already runs `BRepExtrema_DistShapeShape` from every
//! mesh face centroid to a compound of the reference's faces. That call already
//! determines **which reference face won**, and
//! [`crate::stage1::min_distance_to_shape`] throws it away and returns only the
//! distance. `SupportOnShape2` hands it back. So labeling every mesh face with
//! its ground-truth surface — type and exact parameters — costs no additional
//! geometric work, and it runs on the pipeline's own mesh, after repair and
//! quad fusion, which is the mesh Phase 1's earn loop will see. A standalone
//! probe reading raw STL triangles measures a different object:
//! `step_measure_tolerance` needed a `--quad` flag bolted on for exactly that
//! reason.
//!
//! # The three things it reports
//!
//! 1. **Facet-normal exactness**, per reference surface type: the angle between
//!    the mesh facet normal and the exact analytic surface normal. This is the
//!    quantity every closed-form initializer rests on.
//! 2. **Initializer error against ground truth** — cylinder axis from
//!    `n₁ × n₂`, sphere centre from circumcenter normal lines, and the torus
//!    meridian-coplanarity residual — scored per provenance.
//! 3. **The label-ambiguity rate**, which is the untested half of the seam
//!    assumption. Phase 0.5b measured that mesh *vertices* lie on the reference
//!    surfaces; nothing has ever measured that a *facet* does not straddle a
//!    reference face boundary, which is what Design Principle #2 actually
//!    claims and what the earned patch depends on.
//!
//! # Two honest limits, stated because they bound what the numbers mean
//!
//! **The ambiguity rate errs in BOTH directions, so it is an estimate and not
//! a bound.** It is tempting to call it an upper bound on straddling; that is
//! wrong, and the two error directions are different:
//!
//! * *Over-reports.* A facet lying wholly inside reference face A can still be
//!   nearer to face B at some sample point, where B curves back toward it —
//!   near a concave junction especially. Nearest-face labeling is weakest
//!   exactly at junctions, which is exactly where this is measured.
//! * *Under-reports.* The samples cover the middle of a facet, not its rim (see
//!   [`sample_points`] for why they must). A straddle confined to a narrow
//!   sliver along one edge is invisible to them.
//!
//! What makes the number useful despite that is its size: measured at 0.03% of
//! facets over the corpus, both error terms are bounded by something tiny, so
//! the *conclusion* — that facets essentially never cross reference face
//! boundaries — is robust even though the count is not exact. A large rate
//! would be a prompt to look, not a result.
//!
//! **Normals are compared unoriented.** An OCCT face carries an orientation
//! flag that flips its natural normal, and a reference read from STEP may
//! disagree with the mesh's own outward convention face by face. Every angle
//! here is therefore the acute one, `min(θ, 180° − θ)`. That is exactly right
//! for what is being asked — perpendicularity to an axis, concurrence of normal
//! lines, coplanarity of a meridian — none of which has a sign.
//!
//! # What the first corpus run measured (2026-07-30, 79 meshes)
//!
//! **The closed-form axis error obeys `axis error ≈ perpendicularity error /
//! |n₁ × n₂|`.** Verified to 0.1% on `pipe_elbow_10` Parasolid medium, where a
//! perpendicularity max of 0.0280° over a conditioning p50 of 0.0271 predicts
//! 1.033° against a measured 1.0342°, constant across all 348 pairs. That one
//! relation explains every cylinder outlier in the sweep, and it separates the
//! two causes that a raw error number confounds:
//!
//! * *Amplification.* Adjacent facets on a smooth surface differ in normal by
//!   the exporter's angular step, so `|n₁ × n₂| ≈ sin(step)` — measured p50
//!   0.0432 on `stepped_block_toroidal_fillets` fine, which is sin(2.48°) and
//!   recovers Onshape's 2.5° setting. The amplification is therefore ~23× on a
//!   fine mesh and **gets worse the finer the mesh is**. Where the normals are
//!   exact this is the entire error, and a 90° reading — the max on three rows
//!   — is the fully degenerate case of two facets on one ruling strip, which
//!   carries no axis information at all rather than a wrong answer.
//! * *Genuinely non-perpendicular normals*, which is the provenance question
//!   Phase 1.1 asked. Measured only once in the sweep:
//!   `l_bracket_simple` Fusion medium, perpendicularity p50 0.2740°, on a small
//!   hole carrying 96 facets over 75 mm².
//!
//! The practical consequence for Phase 1.2 is that the initializer should pair
//! the seed with the angularly **most distant** admitted facet rather than a
//! neighbour: at 90° separation the amplification factor is 1, and the axis
//! error falls to the perpendicularity error itself — under 0.032° on every
//! mesh in the corpus.
//!
//! **The Delaunay caveat did not reproduce.** The plan's stated worry was that
//! "a general Delaunay triangle's normal is not exactly perpendicular to the
//! axis". On the OCCT re-meshes it reads perpendicularity p50 = p95 = 0.0000°
//! and axis error 0.0000° at every density on five of six models — because
//! OCCT's Delaunay is *constrained* by the pre-discretized face boundaries, so
//! a cylinder's triangles still span two rulings. Fusion's ShapeManager agrees.
//!
//! # Two ways to get the true normal, and why both are here
//!
//! For all five analytic types the surface normal at a sample point equals the
//! normal at that point's footpoint, because projecting onto the surface moves
//! along the normal and leaves the meridian plane (or the plane's constant
//! normal) unchanged. So those need no parameters and no evaluation: the normal
//! comes in closed form from the extracted parameters, exactly.
//!
//! That identity does **not** hold for a freeform surface, so B-spline, Bezier,
//! revolution, extrusion and offset faces are evaluated through OCCT at the
//! footpoint's (u,v), which `ParOnFaceS2` reports from the same distance query
//! that produced the label. `GeomLProp_SLProps` is used rather than a raw
//! `D1U × D1V` because it answers `IsNormalDefined`: a B-spline face has poles
//! and seams where the normal genuinely does not exist, and those facets must be
//! dropped from the statistics rather than counted as an error of some size.
//! The two paths are cross-validated against each other in the tests, on a
//! cylinder where both apply — they agree to under 1e-14.
//!
//! **What freeform coverage showed (2026-07-30).** Over the 8-model NURBS
//! corpus every freeform facet evaluated and **none** landed on a pole or seam.
//! Facet normals there are 2–3 orders of magnitude less exact than on analytic
//! faces — `bspline_terrain` p50 0.9396°, `hill_with_hole` 0.9014° against
//! 0.0000° for planes — but the freeform-vs-analytic comparison is confounded by
//! density, since the freeform models are ~100× coarser per facet. Ruledness
//! separates them better: `bspline_wave`'s *extrusion* faces read p50 0.0507° at
//! 16.7 mm² per facet against `bspline_terrain`'s 0.9396° at 0.53 mm², so the
//! 30× coarser surface is 19× more exact.
//!
//! # Ruledness is necessary but NOT sufficient — alignment is the real condition
//!
//! An earlier version of this note said "ruledness governs facet-normal
//! exactness". That is too strong, and [`min_ruling_gap_rad`] carries the exact
//! statement: a facet's normal is perpendicular to the axis **iff two of its
//! vertices share a ruling**. A tessellator that lays triangles across the
//! ruling family produces skewed normals on a geometrically perfect cylinder.
//!
//! Measured over the corpus's 12 cylinder-bearing models, 12 275 cylinder
//! facets:
//!
//! | group | count | worst perpendicularity |
//! |---|---|---|
//! | edge on a ruling | 9 127 | **0.0006°** |
//! | straddling 3+ rulings | 3 148 (**25.6%**) | **1.0528°** |
//!
//! So alignment is not a rare nicety — a quarter of cylinder facets lack it —
//! and it is the *whole* story: every non-zero perpendicularity in the corpus
//! comes from a straddling facet. Straddling is necessary for error but not
//! sufficient for *large* error, because the tilt grows with a facet's angular
//! spread and shrinks with its axial one; straddling facets run from ~0° to
//! 1.05°. That also explains the extrusion reading above: 0.0507° rather than
//! 0.0000° is itself evidence those facets are not ruling-aligned, since
//! alignment would force exactness.
//!
//! Where alignment comes from, measured: **quad fusion is a 23× discriminator
//! but not a guarantee** — 2.7% of fused quads straddle against 62.3% of
//! triangles. Simple geometry gets it for free (a plain cylinder under OCCT is
//! 120 quads / 0 triangles, 100% aligned, which is the real reason the OCCT
//! re-meshes read 0.0000°); the straddling concentrates where a tessellator
//! cannot lay a clean strip — `stepped_block_toroidal`, whose narrow flat
//! regions between fillets are its known pathology, straddles 49–53% at every
//! OCCT density.
//!
//! **For Phase 1.2 this is better news than the aggregate suggested**, because
//! the test needs no ground truth: given a candidate axis, the minimum angular
//! gap between a facet's vertices says whether that facet's normal can be
//! trusted, so the earn loop can *select* the facets it initializes from rather
//! than hoping. Combined with the amplification law, the rule is: pair two
//! **ruling-aligned** facets that are **angularly distant**.

use crate::stage1::ConnectedMesh;
use opencascade_sys::{
    b_rep, b_rep_adaptor, b_rep_builder_api, b_rep_extrema, extrema, geom, geom_abs, geom_l_prop,
    gp, message, n_collection, top_abs, top_exp, topo_ds, OwnedPtr,
};

// ---------------------------------------------------------------------------
// Small vector helpers (kept local; this module is diagnostics-only)
// ---------------------------------------------------------------------------

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

fn scale(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn norm(a: [f64; 3]) -> f64 {
    dot(a, a).sqrt()
}

/// Unit vector, or `None` if the input is too short to normalize meaningfully.
///
/// The floor is a numeric-conditioning guard, not a tolerance: below it the
/// direction carries no information at all because it is pure rounding noise.
fn unit(a: [f64; 3]) -> Option<[f64; 3]> {
    let n = norm(a);
    if n < 1.0e-30 {
        None
    } else {
        Some(scale(a, 1.0 / n))
    }
}

/// The acute angle between two directions, in degrees. Sign-free by
/// construction — see the module docs on unoriented normals.
fn acute_angle_deg(a: [f64; 3], b: [f64; 3]) -> f64 {
    match (unit(a), unit(b)) {
        (Some(u), Some(v)) => dot(u, v).abs().clamp(0.0, 1.0).acos().to_degrees(),
        _ => f64::NAN,
    }
}

// ---------------------------------------------------------------------------
// Ground truth: the reference's surfaces, as plain data
// ---------------------------------------------------------------------------

/// One reference face's surface, with exact parameters read from the STEP.
///
/// Parameters are copied out as plain `f64` rather than held as OCCT handles so
/// the whole ground truth is `Send + Sync` and needs no lifetime tied to the
/// reader.
#[derive(Debug, Clone, PartialEq)]
pub enum RefSurface {
    Plane {
        point: [f64; 3],
        normal: [f64; 3],
    },
    Cylinder {
        axis_pos: [f64; 3],
        axis_dir: [f64; 3],
        radius: f64,
    },
    Cone {
        apex: [f64; 3],
        axis_dir: [f64; 3],
        half_angle: f64,
    },
    Sphere {
        centre: [f64; 3],
        radius: f64,
    },
    Torus {
        centre: [f64; 3],
        axis_dir: [f64; 3],
        major: f64,
        minor: f64,
    },
    /// A B-spline, Bezier, revolution, extrusion or offset surface — anything
    /// with no closed-form normal. The kind is carried so the report can name
    /// it, and so coverage can be reported honestly: a corpus where half the
    /// area is freeform should not look like one where the *initializers* were
    /// measured on everything. Normals on these are evaluated through OCCT
    /// rather than in closed form — see [`ReferenceModel::normal_at_uv`].
    Other { kind: &'static str },
}

impl RefSurface {
    /// A short name for grouping in the report.
    pub fn type_name(&self) -> &'static str {
        match self {
            RefSurface::Plane { .. } => "plane",
            RefSurface::Cylinder { .. } => "cylinder",
            RefSurface::Cone { .. } => "cone",
            RefSurface::Sphere { .. } => "sphere",
            RefSurface::Torus { .. } => "torus",
            RefSurface::Other { kind } => kind,
        }
    }

    /// The exact unit surface normal at `p`, unoriented.
    ///
    /// Evaluated at `p` itself rather than at its footpoint, which is
    /// equivalent for every analytic type here — see the module docs.
    ///
    /// Returns `None` on the surface's own degeneracies (a point on a
    /// cylinder's axis, a sphere's centre, a cone's apex), where the normal is
    /// genuinely undefined, and for [`RefSurface::Other`], which has no closed
    /// form — use [`ReferenceModel::normal_at_uv`] for those.
    pub fn normal_at(&self, p: [f64; 3]) -> Option<[f64; 3]> {
        match *self {
            RefSurface::Plane { normal, .. } => unit(normal),
            RefSurface::Cylinder {
                axis_pos, axis_dir, ..
            } => {
                let d = unit(axis_dir)?;
                let v = sub(p, axis_pos);
                unit(sub(v, scale(d, dot(v, d))))
            }
            RefSurface::Cone {
                apex,
                axis_dir,
                half_angle,
            } => {
                let d = unit(axis_dir)?;
                let v = sub(p, apex);
                let radial = unit(sub(v, scale(d, dot(v, d))))?;
                // Outward normal of a cone opening along +d tilts back toward
                // -d by the half-angle; at half_angle = 0 this is the cylinder
                // case and reduces to the radial direction.
                unit(add(
                    scale(radial, half_angle.cos()),
                    scale(d, -half_angle.sin()),
                ))
            }
            RefSurface::Sphere { centre, .. } => unit(sub(p, centre)),
            RefSurface::Torus {
                centre,
                axis_dir,
                major,
                ..
            } => {
                let d = unit(axis_dir)?;
                let v = sub(p, centre);
                let meridian = unit(sub(v, scale(d, dot(v, d))))?;
                // The centre of the tube circle through p's meridian.
                let tube_centre = add(centre, scale(meridian, major));
                unit(sub(p, tube_centre))
            }
            RefSurface::Other { .. } => None,
        }
    }

    /// The surface's axis, for the types that have one.
    pub fn axis(&self) -> Option<([f64; 3], [f64; 3])> {
        match *self {
            RefSurface::Cylinder {
                axis_pos, axis_dir, ..
            } => Some((axis_pos, axis_dir)),
            RefSurface::Cone {
                apex, axis_dir, ..
            } => Some((apex, axis_dir)),
            RefSurface::Torus {
                centre, axis_dir, ..
            } => Some((centre, axis_dir)),
            _ => None,
        }
    }
}

/// Read one reference face's surface parameters.
fn read_surface(face: &topo_ds::Face) -> RefSurface {
    let adaptor = b_rep_adaptor::Surface::new_face(face);
    match adaptor.get_type() {
        geom_abs::SurfaceType::Plane => {
            let pln = adaptor.plane();
            let ax = pln.axis();
            let loc = ax.location();
            let dir = ax.direction();
            RefSurface::Plane {
                point: [loc.x(), loc.y(), loc.z()],
                normal: [dir.x(), dir.y(), dir.z()],
            }
        }
        geom_abs::SurfaceType::Cylinder => {
            let cyl = adaptor.cylinder();
            let ax = cyl.axis();
            let loc = ax.location();
            let dir = ax.direction();
            RefSurface::Cylinder {
                axis_pos: [loc.x(), loc.y(), loc.z()],
                axis_dir: [dir.x(), dir.y(), dir.z()],
                radius: cyl.radius(),
            }
        }
        geom_abs::SurfaceType::Cone => {
            let cone = adaptor.cone();
            let ax = cone.axis();
            let dir = ax.direction();
            let apex = cone.apex();
            RefSurface::Cone {
                apex: [apex.x(), apex.y(), apex.z()],
                axis_dir: [dir.x(), dir.y(), dir.z()],
                half_angle: cone.semi_angle(),
            }
        }
        geom_abs::SurfaceType::Sphere => {
            let sph = adaptor.sphere();
            let c = sph.location();
            RefSurface::Sphere {
                centre: [c.x(), c.y(), c.z()],
                radius: sph.radius(),
            }
        }
        geom_abs::SurfaceType::Torus => {
            let tor = adaptor.torus();
            let ax = tor.axis();
            let loc = ax.location();
            let dir = ax.direction();
            RefSurface::Torus {
                centre: [loc.x(), loc.y(), loc.z()],
                axis_dir: [dir.x(), dir.y(), dir.z()],
                major: tor.major_radius(),
                minor: tor.minor_radius(),
            }
        }
        other => RefSurface::Other {
            kind: match other {
                geom_abs::SurfaceType::Bsplinesurface => "bspline",
                geom_abs::SurfaceType::Beziersurface => "bezier",
                geom_abs::SurfaceType::Surfaceofrevolution => "revolution",
                geom_abs::SurfaceType::Surfaceofextrusion => "extrusion",
                geom_abs::SurfaceType::Offsetsurface => "offset",
                _ => "other",
            },
        },
    }
}

/// The reference model: a face compound to measure against, plus one
/// [`RefSurface`] per face **in compound-explorer order**.
///
/// The order matters: [`ReferenceModel::support_index`] recovers a face index
/// by walking the same compound with the same explorer, so the two orders are
/// the same by construction rather than by assumption.
pub struct ReferenceModel {
    compound: OwnedPtr<topo_ds::Compound>,
    /// Face → index, for turning a `SupportOnShape2` result back into a
    /// [`Self::surfaces`] subscript in O(1).
    ///
    /// This is OCCT's `TopTools_IndexedMapOfShape`. That name is a *deprecated
    /// typedef* in OCCT 8.0 and does not appear in the bindings; the canonical
    /// type does, spelled out as its template instantiation. `ShapeMapHasher`
    /// is the variant whose equality is `TopoDS_Shape::IsSame`, which is the
    /// comparison this needs.
    face_index: OwnedPtr<n_collection::IndexedMap_TopoDS_Shape_TopTools_ShapeMapHasher>,
    /// The underlying geometry of each face, for the freeform types that have
    /// no closed-form normal. `BRep_Tool::Surface(F)` applies the face's
    /// location, so these are in world coordinates and share the (u,v)
    /// parameterization that `ParOnFaceS2` reports.
    ///
    /// Only read by [`Self::normal_at_uv`], which runs single-threaded in
    /// [`report`]; labeling never touches these, which is what keeps the
    /// `unsafe impl Sync` on the sharing wrapper in [`label_faces`] honest.
    geometry: Vec<OwnedPtr<geom::Handle_Geom_Surface>>,
    pub surfaces: Vec<RefSurface>,
}

impl ReferenceModel {
    /// Build the ground truth from a `--compare` shape.
    pub fn load(shape: &topo_ds::Shape) -> Self {
        let compound = crate::compare::face_compound(shape);
        let mut surfaces = Vec::new();
        let mut face_index =
            n_collection::IndexedMap_TopoDS_Shape_TopTools_ShapeMapHasher::new();
        let mut explorer = top_exp::Explorer::new_shape_shapeenum2(
            compound.as_shape(),
            top_abs::ShapeEnum::Face,
            top_abs::ShapeEnum::Shape,
        );
        let mut geometry = Vec::new();
        while explorer.more() {
            let face = topo_ds::face_shape(explorer.value());
            surfaces.push(read_surface(face));
            geometry.push(b_rep::Tool::surface_face(face));
            face_index.add(explorer.value());
            explorer.next();
        }
        Self {
            compound,
            face_index,
            geometry,
            surfaces,
        }
    }

    /// The exact unit surface normal of reference face `i` at parameters
    /// `(u, v)`, unoriented — the freeform counterpart of
    /// [`RefSurface::normal_at`].
    ///
    /// Uses `GeomLProp_SLProps`, which reports `IsNormalDefined` rather than
    /// silently returning a zero-length `D1U × D1V`. That matters here: a
    /// B-spline reference face has poles and seams where the normal genuinely
    /// does not exist, and those must be dropped from the statistics rather
    /// than counted as an error of some size.
    ///
    /// `RESOLUTION` is the linear tolerance below which SLProps calls a
    /// derivative vector null. It is a numeric-conditioning floor, not a
    /// threshold: it decides only whether a quantity exists, never whether a
    /// fit is good.
    pub fn normal_at_uv(&self, i: usize, u: f64, v: f64) -> Option<[f64; 3]> {
        const RESOLUTION: f64 = 1.0e-12;
        let surface = self.geometry.get(i)?;
        let mut props =
            geom_l_prop::SLProps::new_handlegeomsurface_real2_int_real(surface, u, v, 1, RESOLUTION);
        if !props.is_normal_defined() {
            return None;
        }
        let n = props.normal();
        unit([n.x(), n.y(), n.z()])
    }

    /// The compound of reference faces — the correct distance target. Never
    /// measure against the solid: see [`crate::compare`].
    pub fn shape(&self) -> &topo_ds::Shape {
        self.compound.as_shape()
    }

    pub fn len(&self) -> usize {
        self.surfaces.len()
    }

    pub fn is_empty(&self) -> bool {
        self.surfaces.is_empty()
    }

    /// Map a `SupportOnShape2` result back to an index into [`Self::surfaces`].
    ///
    /// `FindIndex` is a const hash lookup returning OCCT's 1-based index, or 0
    /// when the shape is absent. The map is built in the same explorer order as
    /// `surfaces`, so the two agree by construction rather than by assumption.
    fn support_index(&self, support: &topo_ds::Shape) -> Option<usize> {
        match self.face_index.find_index(support) {
            0 => None,
            i => Some((i - 1) as usize),
        }
    }
}

// ---------------------------------------------------------------------------
// Labeling
// ---------------------------------------------------------------------------

/// Where a sample point's closest approach to the reference landed.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Support {
    /// Inside a reference face: an unambiguous label.
    Face(usize),
    /// On a reference edge or vertex. The sample is at a seam, which for a
    /// facet interior point means the facet reaches past a face boundary.
    Seam,
    /// The distance query returned nothing.
    None,
}

/// One mesh face's ground truth.
#[derive(Debug, Clone, Copy)]
pub struct FaceLabel {
    /// The reference face nearest the centroid, if the centroid landed in one.
    pub support: Support,
    /// Distance from the centroid to the reference surface, mm.
    pub distance: f64,
    /// The centroid's footpoint parameters on that reference face, when it
    /// landed in one. Needed only for freeform normals; the analytic types are
    /// evaluated from the sample point directly.
    pub uv: Option<(f64, f64)>,
    /// True when the facet's interior sample points do not all agree with the
    /// centroid's label. An estimate of straddling a seam, not a bound: it both
    /// over- and under-reports, for different reasons — see the module docs.
    pub ambiguous: bool,
}

impl FaceLabel {
    /// The reference face index, for faces that got a clean one.
    pub fn face(&self) -> Option<usize> {
        match self.support {
            Support::Face(i) if !self.ambiguous => Some(i),
            _ => None,
        }
    }
}

/// One distance query, returning the winning support as well as the distance.
///
/// This is [`crate::stage1::min_distance_to_shape`] with the support read back
/// instead of discarded — the whole reason this module is cheap.
fn nearest(pt: [f64; 3], model: &ReferenceModel) -> (f64, Support, Option<(f64, f64)>) {
    let p = gp::Pnt::new_real3(pt[0], pt[1], pt[2]);
    let mut vertex = b_rep_builder_api::MakeVertex::new_pnt(&p);
    let vtx_shape = vertex.vertex();
    let progress = message::ProgressRange::new();
    let dist = b_rep_extrema::DistShapeShape::new_shape2_extflag_extalgo_progressrange(
        vtx_shape.as_shape(),
        model.shape(),
        0, // F (unused/obsolete)
        extrema::ExtAlgo::Grad,
        &progress,
    );
    if !dist.is_done() || dist.nb_solution() < 1 {
        return (f64::MAX, Support::None, None);
    }
    let support = match dist.support_type_shape2(1) {
        b_rep_extrema::SupportType::Isinface => {
            let shape = dist.support_on_shape2(1);
            match model.support_index(&shape) {
                Some(i) => Support::Face(i),
                None => Support::None,
            }
        }
        // The closest point on the reference lies on an edge or a vertex of it,
        // so the sample is off the end of any single face's trimmed domain.
        _ => Support::Seam,
    };
    // The footpoint's parameters on the winning face, for the freeform normals
    // that cannot be evaluated from the sample point. Only meaningful when the
    // support is a face.
    let uv = match support {
        Support::Face(_) => {
            let (mut u, mut v) = (0.0_f64, 0.0_f64);
            dist.par_on_face_s2(1, &mut u, &mut v);
            Some((u, v))
        }
        _ => None,
    };
    (dist.value(), support, uv)
}

/// The interior sample points used to test a facet's label for agreement.
///
/// The centroid plus one point halfway from the centroid to each vertex. All
/// are strictly interior, so a facet that lies wholly within one reference face
/// gives the same label at every one of them; a facet reaching across a seam
/// does not. Vertices themselves are deliberately not used: a mesh vertex *on*
/// a seam is expected and says nothing.
fn sample_points(mesh: &ConnectedMesh, fi: usize) -> Vec<[f64; 3]> {
    let face = &mesh.faces[fi];
    let n = face.vertex_count as usize;
    let mut centroid = [0.0_f64; 3];
    for k in 0..n {
        let v = &mesh.vertices[face.vertex_indices[k]];
        centroid[0] += v.x;
        centroid[1] += v.y;
        centroid[2] += v.z;
    }
    let inv = 1.0 / n as f64;
    let centroid = scale(centroid, inv);

    let mut pts = Vec::with_capacity(n + 1);
    pts.push(centroid);
    for k in 0..n {
        let v = &mesh.vertices[face.vertex_indices[k]];
        pts.push([
            0.5 * (centroid[0] + v.x),
            0.5 * (centroid[1] + v.y),
            0.5 * (centroid[2] + v.z),
        ]);
    }
    pts
}

/// Label every mesh face with its ground-truth reference face.
///
/// Parallel across faces using the same pattern and the same safety argument as
/// [`crate::stage1::distances_to_shape`]: `BRepExtrema_DistShapeShape`,
/// `BRepBuilderAPI_MakeVertex`, `TopExp_Explorer` and the `gp` value types have
/// no global mutable state and do not mutate the target shape, so sharing a
/// read-only `&Shape` across threads is sound (LEARNINGS.md, "OCCT Thread
/// Safety").
pub fn label_faces(mesh: &ConnectedMesh, model: &ReferenceModel) -> Vec<FaceLabel> {
    let n = mesh.faces.len();
    let unlabeled = FaceLabel {
        support: Support::None,
        distance: f64::MAX,
        uv: None,
        ambiguous: false,
    };
    if n == 0 || model.is_empty() {
        return vec![unlabeled; n];
    }

    let label_one = |fi: usize| -> FaceLabel {
        if mesh.faces[fi].vertex_count == 0 {
            return unlabeled;
        }
        let pts = sample_points(mesh, fi);
        let (distance, support, uv) = nearest(pts[0], model);
        let mut ambiguous = false;
        for &p in &pts[1..] {
            let (_, s, _) = nearest(p, model);
            if s != support {
                ambiguous = true;
                break;
            }
        }
        FaceLabel {
            support,
            distance,
            uv,
            ambiguous,
        }
    };

    let threads = std::thread::available_parallelism()
        .map(|p| p.get())
        .unwrap_or(1)
        .min(8)
        .min(n.div_ceil(16).max(1));
    if threads <= 1 {
        return (0..n).map(label_one).collect();
    }

    struct SharedModel(*const ReferenceModel);
    unsafe impl Send for SharedModel {}
    unsafe impl Sync for SharedModel {}
    let shared = SharedModel(model as *const _);

    let chunk = n.div_ceil(threads);
    let mut out = vec![unlabeled; n];
    std::thread::scope(|s| {
        for (ci, out_chunk) in out.chunks_mut(chunk).enumerate() {
            let shared = &shared;
            s.spawn(move || {
                let model = unsafe { &*shared.0 };
                let base = ci * chunk;
                for (k, slot) in out_chunk.iter_mut().enumerate() {
                    let fi = base + k;
                    if mesh.faces[fi].vertex_count == 0 {
                        continue;
                    }
                    let pts = sample_points(mesh, fi);
                    let (distance, support, uv) = nearest(pts[0], model);
                    let mut ambiguous = false;
                    for &p in &pts[1..] {
                        let (_, s, _) = nearest(p, model);
                        if s != support {
                            ambiguous = true;
                            break;
                        }
                    }
                    *slot = FaceLabel {
                        support,
                        distance,
                        uv,
                        ambiguous,
                    };
                }
            });
        }
    });
    out
}

/// Below this angular gap (radians) two of a facet's vertices are taken to lie
/// on the same ruling.
///
/// Not a tuned threshold: the corpus is bimodal with five orders of magnitude
/// of empty space between the two modes — vertices on one ruling differ only by
/// float32 mesh noise (~1e-7 rad), while distinct rulings are separated by the
/// exporter's angular step (≥1e-2 rad). Any cut in that gap gives the same
/// split. The gap is reported beside the split so a mesh that lands inside it
/// announces itself.
const RULING_SAME_RAD: f64 = 1.0e-6;

// ---------------------------------------------------------------------------
// Distributions
// ---------------------------------------------------------------------------

/// p50 / p95 / max of a sample, with its size. The shape of a distribution
/// names the cause on sight where a max alone cannot — the lesson from the
/// wild-pair corpus in Phase 0.5b, applied here too.
#[derive(Debug, Clone, Copy, Default)]
pub struct Dist {
    pub p50: f64,
    pub p95: f64,
    pub max: f64,
    pub n: usize,
}

impl Dist {
    fn of(mut v: Vec<f64>) -> Self {
        v.retain(|x| x.is_finite());
        if v.is_empty() {
            return Self {
                p50: f64::NAN,
                p95: f64::NAN,
                max: f64::NAN,
                n: 0,
            };
        }
        v.sort_by(f64::total_cmp);
        let at = |q: f64| v[(((v.len() - 1) as f64) * q).round() as usize];
        Self {
            p50: at(0.50),
            p95: at(0.95),
            max: at(1.0),
            n: v.len(),
        }
    }

    fn is_empty(&self) -> bool {
        self.n == 0
    }
}

// ---------------------------------------------------------------------------
// Geometry the initializers are built from
// ---------------------------------------------------------------------------

/// Area of a mesh face (a quad as two triangles).
fn face_area(mesh: &ConnectedMesh, fi: usize) -> f64 {
    let f = &mesh.faces[fi];
    let p = |k: usize| {
        let v = &mesh.vertices[f.vertex_indices[k]];
        [v.x, v.y, v.z]
    };
    let tri = |a: [f64; 3], b: [f64; 3], c: [f64; 3]| 0.5 * norm(cross(sub(b, a), sub(c, a)));
    let mut area = tri(p(0), p(1), p(2));
    if f.vertex_count == 4 {
        area += tri(p(0), p(2), p(3));
    }
    area
}

fn face_centroid(mesh: &ConnectedMesh, fi: usize) -> [f64; 3] {
    let f = &mesh.faces[fi];
    let n = f.vertex_count as usize;
    let mut c = [0.0_f64; 3];
    for k in 0..n {
        let v = &mesh.vertices[f.vertex_indices[k]];
        c = add(c, [v.x, v.y, v.z]);
    }
    scale(c, 1.0 / n as f64)
}

/// The circumcenter of a facet's first three vertices.
///
/// For a facet inscribed in a sphere this is the point whose normal line passes
/// **exactly** through the centre, which is what makes the sphere initializer
/// closed-form. A quad here is planar by construction (stage 1 fuses only
/// coplanar triangles), so its four vertices are concyclic and any three give
/// the same circumcenter.
fn circumcenter(mesh: &ConnectedMesh, fi: usize) -> Option<[f64; 3]> {
    let f = &mesh.faces[fi];
    if (f.vertex_count as usize) < 3 {
        return None;
    }
    let p = |k: usize| {
        let v = &mesh.vertices[f.vertex_indices[k]];
        [v.x, v.y, v.z]
    };
    let (a, b, c) = (p(0), p(1), p(2));
    let ab = sub(b, a);
    let ac = sub(c, a);
    let n = cross(ab, ac);
    let denom = 2.0 * dot(n, n);
    if denom < 1.0e-30 {
        return None;
    }
    // Standard barycentric form: a + (|ac|² (ab×ac)×ab + |ab|² ac×(ab×ac)) / 2|ab×ac|²
    let term = add(
        scale(cross(n, ab), dot(ac, ac)),
        scale(cross(ac, n), dot(ab, ab)),
    );
    Some(add(a, scale(term, 1.0 / denom)))
}

/// Distance between two lines given as (point, direction).
///
/// Returns `None` when the directions are too nearly parallel for the common
/// perpendicular to be meaningful — a numeric-conditioning guard, not a
/// tolerance.
fn line_line_distance(p1: [f64; 3], d1: [f64; 3], p2: [f64; 3], d2: [f64; 3]) -> Option<f64> {
    let n = cross(d1, d2);
    let len = norm(n);
    if len < 1.0e-12 {
        return None;
    }
    Some((dot(sub(p2, p1), n) / len).abs())
}

/// Midpoint of the common perpendicular between two lines — the standard
/// near-intersection of two normal lines, and the sphere initializer's centre.
fn common_perpendicular_midpoint(
    p1: [f64; 3],
    d1: [f64; 3],
    p2: [f64; 3],
    d2: [f64; 3],
) -> Option<[f64; 3]> {
    let (u, v) = (unit(d1)?, unit(d2)?);
    let w = sub(p1, p2);
    let b = dot(u, v);
    let denom = 1.0 - b * b;
    if denom < 1.0e-12 {
        return None;
    }
    let d = dot(u, w);
    let e = dot(v, w);
    let s = (b * e - d) / denom;
    let t = (e - b * d) / denom;
    let q1 = add(p1, scale(u, s));
    let q2 = add(p2, scale(v, t));
    Some(scale(add(q1, q2), 0.5))
}

// ---------------------------------------------------------------------------
// The report
// ---------------------------------------------------------------------------

/// How exactly facet normals reproduce the true surface normal, by type.
#[derive(Debug, Clone, Default)]
pub struct NormalScore {
    /// (type name, angle error in degrees, area on that type)
    pub by_type: Vec<(&'static str, Dist, f64)>,
}

/// The cylinder initializer, scored against ground truth.
#[derive(Debug, Clone, Default)]
pub struct CylinderScore {
    /// `|90° − ∠(facet normal, true axis)|`, degrees. The property the
    /// initializer rests on: zero means the normal is exactly perpendicular to
    /// the axis, which is what makes `n₁ × n₂` the axis exactly.
    pub perpendicularity_deg: Dist,
    /// Angle between `normalize(n₁ × n₂)` over an adjacent same-surface pair
    /// and the true axis, degrees.
    pub pair_axis_deg: Dist,
    /// `|n₁ × n₂|` for those pairs — the conditioning of the cross product.
    pub conditioning: Dist,
    /// Smallest angular gap between any two of a facet's vertices around the
    /// axis. Zero means the facet has an edge along a ruling.
    pub ruling_gap_rad: Dist,
    /// Perpendicularity for facets **with** an edge along a ruling, where the
    /// normal is perpendicular to the axis exactly.
    pub perp_ruling_aligned: Dist,
    /// Perpendicularity for facets straddling three or more rulings, where it
    /// is not — see [`min_ruling_gap_rad`].
    pub perp_ruling_skew: Dist,
    /// Quad/triangle split of the two groups above. A stage-1 fused quad is a
    /// coplanar triangle pair, and on a cylinder that all but forces a
    /// two-ruling span; whether it forces it is what these count.
    pub aligned_quads: usize,
    pub aligned_tris: usize,
    pub skew_quads: usize,
    pub skew_tris: usize,
    /// Median pair-axis error binned by conditioning, worst-conditioned first:
    /// (bin label, count, median error in degrees).
    pub by_conditioning: Vec<(&'static str, usize, f64)>,
}

/// The sphere initializer, scored against ground truth.
#[derive(Debug, Clone, Default)]
pub struct SphereScore {
    /// Distance from the true centre to the facet's circumcenter normal line,
    /// mm. Exactly zero for a facet inscribed in the sphere.
    pub normal_line_miss_mm: Dist,
    /// Error of the two-facet centre estimate, mm.
    pub pair_centre_mm: Dist,
    /// Error of the radius implied by that centre, mm.
    pub pair_radius_mm: Dist,
    /// `|n₁ × n₂|` for those pairs. Two normal lines that are nearly parallel
    /// meet nowhere useful, so this separates "the normals are wrong" from
    /// "the near-intersection was ill-conditioned" — the same split the
    /// cylinder score makes.
    pub conditioning: Dist,
}

/// The torus initializer, scored against ground truth.
#[derive(Debug, Clone, Default)]
pub struct TorusScore {
    /// Distance between the facet's normal line and the true axis, mm — the
    /// residual of the meridian-coplanarity constraint
    /// `(pᵢ − a)·(nᵢ × d) = 0` that the axis solve is built on.
    pub meridian_miss_mm: Dist,
    /// `|nᵢ × d|`, the conditioning of that constraint. Near zero where the
    /// normal is nearly parallel to the axis and the row carries no
    /// information.
    pub conditioning: Dist,
}

/// Our stage-2 partition against the reference partition.
#[derive(Debug, Clone, Default)]
pub struct SegmentationScore {
    pub ref_faces: usize,
    pub our_surfaces: usize,
    /// Distinct (our surface, reference face) pairs. Equals both counts above
    /// when the partitions agree exactly.
    pub pairs: usize,
    /// Area-weighted: `Σ_R max_S area(R,S) / Σ area`. Below 1 means a reference
    /// face was split across several of our surfaces — **over-segmentation**,
    /// which no existing gate can see.
    pub ref_coverage: f64,
    /// Area-weighted: `Σ_S max_R area(S,R) / Σ area`. Below 1 means one of our
    /// surfaces spans several reference faces — **under-segmentation**.
    pub surf_purity: f64,
}

/// Everything measured against the reference for one model.
#[derive(Debug, Clone, Default)]
pub struct Report {
    pub ref_faces: usize,
    pub mesh_faces: usize,
    pub labeled: usize,
    pub ambiguous: usize,
    pub seam: usize,
    pub unlabeled: usize,
    /// Fraction of mesh area whose reference face is an analytic type. The
    /// **initializer** scores below say nothing about the rest — but the normal
    /// scores now cover freeform faces too.
    pub analytic_area_fraction: f64,
    /// Facets on a freeform reference face whose normal OCCT could evaluate.
    pub freeform_evaluated: usize,
    /// Facets on a freeform reference face at a pole or seam, where the surface
    /// normal does not exist and the facet is excluded rather than scored.
    pub freeform_undefined: usize,
    pub normals: NormalScore,
    pub cylinder: Option<CylinderScore>,
    pub sphere: Option<SphereScore>,
    pub torus: Option<TorusScore>,
    pub segmentation: Option<SegmentationScore>,
}

/// Score a labeled mesh against its ground truth.
pub fn report(mesh: &ConnectedMesh, model: &ReferenceModel, labels: &[FaceLabel]) -> Report {
    let mut rep = Report {
        ref_faces: model.len(),
        mesh_faces: mesh.faces.iter().filter(|f| f.vertex_count > 0).count(),
        ..Default::default()
    };

    // --- label accounting -------------------------------------------------
    let mut area_total = 0.0_f64;
    let mut area_analytic = 0.0_f64;
    let mut area_by_type: Vec<(&'static str, f64)> = Vec::new();
    let mut normals_by_type: Vec<(&'static str, Vec<f64>)> = Vec::new();
    let mut freeform_evaluated = 0_usize;
    let mut freeform_undefined = 0_usize;

    for (fi, label) in labels.iter().enumerate() {
        if mesh.faces[fi].vertex_count == 0 {
            continue;
        }
        match label.support {
            Support::Seam => rep.seam += 1,
            Support::None => rep.unlabeled += 1,
            Support::Face(_) if label.ambiguous => rep.ambiguous += 1,
            Support::Face(_) => rep.labeled += 1,
        }

        let Some(ri) = label.face() else { continue };
        let surf = &model.surfaces[ri];
        let area = face_area(mesh, fi);
        area_total += area;
        let name = surf.type_name();
        match area_by_type.iter_mut().find(|(n, _)| *n == name) {
            Some(slot) => slot.1 += area,
            None => area_by_type.push((name, area)),
        }
        if !matches!(surf, RefSurface::Other { .. }) {
            area_analytic += area;
        }

        // The true surface normal under this facet. Analytic types are
        // evaluated in closed form at the facet's own centroid, which is exact
        // and needs no parameters (see the module docs); freeform types have no
        // closed form, so they are evaluated through OCCT at the centroid's
        // footpoint, where `None` means the normal is genuinely undefined —
        // a pole or a seam — and the facet is dropped rather than scored.
        let true_normal = match surf {
            RefSurface::Other { .. } => label
                .uv
                .and_then(|(u, v)| model.normal_at_uv(ri, u, v))
                .inspect(|_| freeform_evaluated += 1)
                .or_else(|| {
                    freeform_undefined += 1;
                    None
                }),
            _ => surf.normal_at(face_centroid(mesh, fi)),
        };

        // Facet normal against the exact surface normal.
        let (Some(fn_), Some(sn)) = (mesh.faces[fi].normal, true_normal) else {
            continue;
        };
        let err = acute_angle_deg(fn_, sn);
        match normals_by_type.iter_mut().find(|(n, _)| *n == name) {
            Some(slot) => slot.1.push(err),
            None => normals_by_type.push((name, vec![err])),
        }
    }
    rep.analytic_area_fraction = if area_total > 0.0 {
        area_analytic / area_total
    } else {
        0.0
    };
    rep.freeform_evaluated = freeform_evaluated;
    rep.freeform_undefined = freeform_undefined;
    rep.normals.by_type = normals_by_type
        .into_iter()
        .map(|(name, errs)| {
            let area = area_by_type
                .iter()
                .find(|(n, _)| *n == name)
                .map(|(_, a)| *a)
                .unwrap_or(0.0);
            (name, Dist::of(errs), area)
        })
        .collect();
    rep.normals
        .by_type
        .sort_by(|a, b| b.2.total_cmp(&a.2));

    // --- initializer scores ----------------------------------------------
    rep.cylinder = score_cylinders(mesh, model, labels);
    rep.sphere = score_spheres(mesh, model, labels);
    rep.torus = score_tori(mesh, model, labels);
    rep
}

/// Adjacent mesh-face pairs that share a clean label, each reported once.
fn same_label_pairs(mesh: &ConnectedMesh, labels: &[FaceLabel]) -> Vec<(usize, usize, usize)> {
    let mut pairs = Vec::new();
    for (fi, face) in mesh.faces.iter().enumerate() {
        let Some(ri) = labels[fi].face() else { continue };
        for k in 0..face.vertex_count as usize {
            let nb = face.neighbors[k];
            if nb < 0 {
                continue;
            }
            let fj = nb as usize;
            if fj <= fi {
                continue; // report each pair once
            }
            if labels[fj].face() == Some(ri) {
                pairs.push((fi, fj, ri));
            }
        }
    }
    pairs
}

/// The smallest angular gap, in radians, between any two of a facet's vertices
/// measured around `axis` — zero exactly when two vertices share a ruling.
///
/// This is the quantity that decides whether a facet normal on a ruled surface
/// is perpendicular to the axis at all. Parametrize the cylinder as
/// `p(θ, z) = a + r(cos θ · u + sin θ · v) + z · d`. The facet's plane contains
/// the direction `d` — equivalently `n ⟂ d` — iff the vertices, projected along
/// `d` onto the plane perpendicular to it, are collinear. They land on a circle
/// of radius `r`, and three points on a circle are collinear only when two of
/// them coincide. **So `n ⟂ d` exactly when two vertices share a ruling**, and a
/// triangle straddling three distinct rulings tilts its normal off
/// perpendicular however perfect the underlying cylinder is.
fn min_ruling_gap_rad(mesh: &ConnectedMesh, fi: usize, axis_pos: [f64; 3], axis: [f64; 3]) -> Option<f64> {
    let d = unit(axis)?;
    let seed = if d[0].abs() < 0.9 { [1.0, 0.0, 0.0] } else { [0.0, 1.0, 0.0] };
    let u = unit(cross(d, seed))?;
    let v = cross(d, u);

    let face = &mesh.faces[fi];
    let n = face.vertex_count as usize;
    let mut thetas = Vec::with_capacity(n);
    for k in 0..n {
        let vtx = &mesh.vertices[face.vertex_indices[k]];
        let w = sub([vtx.x, vtx.y, vtx.z], axis_pos);
        let w = sub(w, scale(d, dot(w, d)));
        thetas.push(dot(w, v).atan2(dot(w, u)));
    }

    let mut min_gap = f64::MAX;
    for i in 0..n {
        for j in (i + 1)..n {
            let mut delta = (thetas[i] - thetas[j]).abs();
            if delta > std::f64::consts::PI {
                delta = std::f64::consts::TAU - delta;
            }
            min_gap = min_gap.min(delta);
        }
    }
    Some(min_gap)
}

fn score_cylinders(
    mesh: &ConnectedMesh,
    model: &ReferenceModel,
    labels: &[FaceLabel],
) -> Option<CylinderScore> {
    let mut perp = Vec::new();
    let mut gaps = Vec::new();
    let (mut perp_aligned, mut perp_skew) = (Vec::new(), Vec::new());
    let (mut aligned_quads, mut aligned_tris) = (0_usize, 0_usize);
    let (mut skew_quads, mut skew_tris) = (0_usize, 0_usize);
    for (fi, label) in labels.iter().enumerate() {
        let Some(ri) = label.face() else { continue };
        let RefSurface::Cylinder {
            axis_pos, axis_dir, ..
        } = model.surfaces[ri]
        else {
            continue;
        };
        let Some(n) = mesh.faces[fi].normal else {
            continue;
        };
        let err = (90.0 - acute_angle_deg(n, axis_dir)).abs();
        perp.push(err);
        // Split on whether the facet has an edge along a ruling. The cut sits
        // in a measured gap, not on a tuned value: two vertices on one ruling
        // differ in θ only by float32 mesh noise, ~1e-7 rad at these radii,
        // while two rulings are separated by the exporter's angular step,
        // ≥ 1e-2 rad. Any cut between those gives the same split, and the gap
        // itself is reported beside it so a corpus that lands inside it says so.
        if let Some(gap) = min_ruling_gap_rad(mesh, fi, axis_pos, axis_dir) {
            gaps.push(gap);
            let quad = mesh.faces[fi].vertex_count == 4;
            if gap < RULING_SAME_RAD {
                perp_aligned.push(err);
                if quad {
                    aligned_quads += 1;
                } else {
                    aligned_tris += 1;
                }
            } else {
                perp_skew.push(err);
                if quad {
                    skew_quads += 1;
                } else {
                    skew_tris += 1;
                }
            }
        }
    }
    if perp.is_empty() {
        return None;
    }

    let mut pair_err = Vec::new();
    let mut cond = Vec::new();
    let mut binned: Vec<(&'static str, Vec<f64>)> = vec![
        ("|n1xn2| < 0.01", Vec::new()),
        ("0.01 - 0.1", Vec::new()),
        ("0.1 - 0.5", Vec::new()),
        ("> 0.5", Vec::new()),
    ];
    for (fi, fj, ri) in same_label_pairs(mesh, labels) {
        let RefSurface::Cylinder { axis_dir, .. } = model.surfaces[ri] else {
            continue;
        };
        let (Some(n1), Some(n2)) = (mesh.faces[fi].normal, mesh.faces[fj].normal) else {
            continue;
        };
        let c = cross(n1, n2);
        let k = norm(c);
        let Some(est) = unit(c) else { continue };
        let err = acute_angle_deg(est, axis_dir);
        pair_err.push(err);
        cond.push(k);
        // Binning is a presentation choice for a diagnostic — it separates
        // "the normals are wrong" from "the cross product was ill-conditioned",
        // and no pipeline decision reads it.
        let bin = if k < 0.01 {
            0
        } else if k < 0.1 {
            1
        } else if k < 0.5 {
            2
        } else {
            3
        };
        binned[bin].1.push(err);
    }

    Some(CylinderScore {
        perpendicularity_deg: Dist::of(perp),
        ruling_gap_rad: Dist::of(gaps),
        perp_ruling_aligned: Dist::of(perp_aligned),
        perp_ruling_skew: Dist::of(perp_skew),
        aligned_quads,
        aligned_tris,
        skew_quads,
        skew_tris,
        pair_axis_deg: Dist::of(pair_err),
        conditioning: Dist::of(cond),
        by_conditioning: binned
            .into_iter()
            .map(|(name, v)| {
                let n = v.len();
                (name, n, Dist::of(v).p50)
            })
            .filter(|(_, n, _)| *n > 0)
            .collect(),
    })
}

fn score_spheres(
    mesh: &ConnectedMesh,
    model: &ReferenceModel,
    labels: &[FaceLabel],
) -> Option<SphereScore> {
    let mut miss = Vec::new();
    for (fi, label) in labels.iter().enumerate() {
        let Some(ri) = label.face() else { continue };
        let RefSurface::Sphere { centre, .. } = model.surfaces[ri] else {
            continue;
        };
        let (Some(cc), Some(n)) = (circumcenter(mesh, fi), mesh.faces[fi].normal) else {
            continue;
        };
        // Distance from the true centre to the circumcenter normal line.
        let v = sub(centre, cc);
        let Some(u) = unit(n) else { continue };
        miss.push(norm(sub(v, scale(u, dot(v, u)))));
    }
    if miss.is_empty() {
        return None;
    }

    let mut centre_err = Vec::new();
    let mut radius_err = Vec::new();
    let mut cond = Vec::new();
    for (fi, fj, ri) in same_label_pairs(mesh, labels) {
        let RefSurface::Sphere { centre, radius } = model.surfaces[ri] else {
            continue;
        };
        let (Some(c1), Some(c2)) = (circumcenter(mesh, fi), circumcenter(mesh, fj)) else {
            continue;
        };
        let (Some(n1), Some(n2)) = (mesh.faces[fi].normal, mesh.faces[fj].normal) else {
            continue;
        };
        let Some(est) = common_perpendicular_midpoint(c1, n1, c2, n2) else {
            continue;
        };
        centre_err.push(norm(sub(est, centre)));
        let r_est = 0.5 * (norm(sub(c1, est)) + norm(sub(c2, est)));
        radius_err.push((r_est - radius).abs());
        cond.push(norm(cross(n1, n2)));
    }

    Some(SphereScore {
        normal_line_miss_mm: Dist::of(miss),
        pair_centre_mm: Dist::of(centre_err),
        pair_radius_mm: Dist::of(radius_err),
        conditioning: Dist::of(cond),
    })
}

fn score_tori(
    mesh: &ConnectedMesh,
    model: &ReferenceModel,
    labels: &[FaceLabel],
) -> Option<TorusScore> {
    let mut miss = Vec::new();
    let mut cond = Vec::new();
    for (fi, label) in labels.iter().enumerate() {
        let Some(ri) = label.face() else { continue };
        let RefSurface::Torus {
            centre, axis_dir, ..
        } = model.surfaces[ri]
        else {
            continue;
        };
        let Some(n) = mesh.faces[fi].normal else {
            continue;
        };
        let p = face_centroid(mesh, fi);
        cond.push(norm(cross(n, axis_dir)));
        if let Some(d) = line_line_distance(p, n, centre, axis_dir) {
            miss.push(d);
        }
    }
    if miss.is_empty() {
        return None;
    }
    Some(TorusScore {
        meridian_miss_mm: Dist::of(miss),
        conditioning: Dist::of(cond),
    })
}

/// Score our stage-2 partition against the reference partition.
///
/// `our_partition` is one face-index list per selected surface.
pub fn score_segmentation(
    mesh: &ConnectedMesh,
    model: &ReferenceModel,
    labels: &[FaceLabel],
    our_partition: &[Vec<usize>],
) -> SegmentationScore {
    // Area of every (our surface, reference face) overlap.
    let mut overlap: std::collections::HashMap<(usize, usize), f64> =
        std::collections::HashMap::new();
    let mut total = 0.0_f64;
    for (si, faces) in our_partition.iter().enumerate() {
        for &fi in faces {
            let Some(ri) = labels[fi].face() else { continue };
            let a = face_area(mesh, fi);
            *overlap.entry((si, ri)).or_insert(0.0) += a;
            total += a;
        }
    }
    if total <= 0.0 {
        return SegmentationScore {
            ref_faces: model.len(),
            our_surfaces: our_partition.len(),
            ..Default::default()
        };
    }

    let mut best_per_ref: std::collections::HashMap<usize, f64> = std::collections::HashMap::new();
    let mut best_per_surf: std::collections::HashMap<usize, f64> = std::collections::HashMap::new();
    for (&(si, ri), &a) in &overlap {
        let e = best_per_ref.entry(ri).or_insert(0.0);
        *e = e.max(a);
        let e = best_per_surf.entry(si).or_insert(0.0);
        *e = e.max(a);
    }

    SegmentationScore {
        ref_faces: model.len(),
        our_surfaces: our_partition.len(),
        pairs: overlap.len(),
        ref_coverage: best_per_ref.values().sum::<f64>() / total,
        surf_purity: best_per_surf.values().sum::<f64>() / total,
    }
}

// ---------------------------------------------------------------------------
// Printing
// ---------------------------------------------------------------------------

fn fmt_dist_deg(d: &Dist) -> String {
    if d.is_empty() {
        return "—".to_string();
    }
    format!(
        "p50 {:.4}° p95 {:.4}° max {:.4}° (n={})",
        d.p50, d.p95, d.max, d.n
    )
}

fn fmt_dist_mm(d: &Dist) -> String {
    if d.is_empty() {
        return "—".to_string();
    }
    format!(
        "p50 {:.3e} p95 {:.3e} max {:.3e} mm (n={})",
        d.p50, d.p95, d.max, d.n
    )
}

impl Report {
    /// Print the mesh-side half: labeling, seam assumption, normals,
    /// initializers.
    pub fn print(&self) {
        eprintln!(
            "  Reference: {} faces; labeled {}/{} mesh faces, {} ambiguous, {} at a seam, {} unlabeled",
            self.ref_faces, self.labeled, self.mesh_faces, self.ambiguous, self.seam, self.unlabeled
        );
        let straddle = self.ambiguous + self.seam;
        if self.mesh_faces > 0 {
            eprintln!(
                "    seam assumption: {:.4}% of facets disagree across their own interior \
                 (estimate, not a bound; initializer scores cover the {:.1}% of \
                 labeled area that is analytic)",
                100.0 * straddle as f64 / self.mesh_faces as f64,
                100.0 * self.analytic_area_fraction,
            );
        }
        if self.freeform_evaluated + self.freeform_undefined > 0 {
            eprintln!(
                "    freeform normals: {} evaluated, {} at a pole or seam where the normal \
                 does not exist (excluded)",
                self.freeform_evaluated, self.freeform_undefined
            );
        }
        for (name, dist, area) in &self.normals.by_type {
            eprintln!(
                "    normal error on {name:<9} {}  [{:.1} mm² of mesh]",
                fmt_dist_deg(dist),
                area
            );
        }
        if let Some(c) = &self.cylinder {
            eprintln!(
                "    cylinder  perpendicularity {}",
                fmt_dist_deg(&c.perpendicularity_deg)
            );
            eprintln!(
                "              n1xn2 axis error {}",
                fmt_dist_deg(&c.pair_axis_deg)
            );
            eprintln!(
                "              ruling-aligned   {} facets with an edge on a ruling {}, \
                 {} straddling 3+ rulings {}",
                c.perp_ruling_aligned.n,
                fmt_dist_deg(&c.perp_ruling_aligned),
                c.perp_ruling_skew.n,
                fmt_dist_deg(&c.perp_ruling_skew),
            );
            eprintln!(
                "              aligned {} quad / {} tri, straddling {} quad / {} tri",
                c.aligned_quads, c.aligned_tris, c.skew_quads, c.skew_tris
            );
            eprintln!(
                "              min ruling gap   p50 {:.3e} p95 {:.3e} rad",
                c.ruling_gap_rad.p50, c.ruling_gap_rad.p95
            );
            eprintln!(
                "              |n1xn2|          p50 {:.4} p95 {:.4} max {:.4}  \
                 (axis error ≈ perpendicularity / |n1xn2|)",
                c.conditioning.p50, c.conditioning.p95, c.conditioning.max
            );
            let bins = c
                .by_conditioning
                .iter()
                .map(|(name, n, p50)| format!("{name}: n={n} p50 {p50:.4}°"))
                .collect::<Vec<_>>()
                .join("; ");
            eprintln!("              by conditioning  {bins}");
        }
        if let Some(s) = &self.sphere {
            eprintln!(
                "    sphere    normal-line miss {}",
                fmt_dist_mm(&s.normal_line_miss_mm)
            );
            eprintln!(
                "              pair centre err  {}",
                fmt_dist_mm(&s.pair_centre_mm)
            );
            eprintln!(
                "              pair radius err  {}, |n1xn2| p50 {:.4} p95 {:.4}",
                fmt_dist_mm(&s.pair_radius_mm),
                s.conditioning.p50,
                s.conditioning.p95
            );
        }
        if let Some(t) = &self.torus {
            eprintln!(
                "    torus     meridian residual {}",
                fmt_dist_mm(&t.meridian_miss_mm)
            );
            eprintln!(
                "              |n x d|          {}",
                fmt_dist_mm(&t.conditioning)
            );
        }
    }
}

impl SegmentationScore {
    pub fn print(&self) {
        eprintln!(
            "  Segmentation vs reference: {} reference faces, {} of ours, {} overlapping pairs",
            self.ref_faces, self.our_surfaces, self.pairs
        );
        eprintln!(
            "    ref coverage {:.4} (1 = no over-segmentation), surface purity {:.4} (1 = no merging)",
            self.ref_coverage, self.surf_purity
        );
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::stage1::{MeshFace, MeshVertex};

    /// An arbitrary unit direction, deliberately not axis-aligned: an axis
    /// along +Z would let a sign or an index slip through every assertion here.
    fn oblique_axis() -> [f64; 3] {
        unit([0.3, -0.5, 0.81]).unwrap()
    }

    /// `|â × b̂|` = sin of the angle between two directions — how far from
    /// parallel they are, and the *only* well-conditioned way to assert
    /// near-parallelism.
    ///
    /// [`acute_angle_deg`] cannot be used for this: it goes through `acos`,
    /// which near an argument of 1 amplifies a 1e-16 rounding error into an
    /// angle of about `sqrt(2 · 1e-16)` ≈ 1.5e-8 rad ≈ 1e-6 degrees. So no
    /// angle assertion near 0° can be tightened past ~1e-6°, however exact the
    /// underlying geometry is — measured here, where an exact strip cylinder
    /// failed a 1e-9° assertion at one sample out of sixteen. The cross product
    /// is linear in the error and resolves to machine epsilon.
    ///
    /// The same floor applies to the *reported* angles in this module. They are
    /// printed to four decimal places of a degree, far above it, so nothing in
    /// the report is affected — but do not read exactness below ~1e-6° out of
    /// a number that came from `acos`.
    fn parallel_miss(a: [f64; 3], b: [f64; 3]) -> f64 {
        match (unit(a), unit(b)) {
            (Some(u), Some(v)) => norm(cross(u, v)),
            _ => f64::NAN,
        }
    }

    /// Two unit vectors spanning the plane perpendicular to `d`.
    fn frame(d: [f64; 3]) -> ([f64; 3], [f64; 3]) {
        let seed = if d[0].abs() < 0.9 {
            [1.0, 0.0, 0.0]
        } else {
            [0.0, 1.0, 0.0]
        };
        let u = unit(cross(d, seed)).unwrap();
        let v = cross(d, u);
        (u, v)
    }

    #[test]
    fn plane_normal_is_the_axis() {
        let s = RefSurface::Plane {
            point: [1.0, 2.0, 3.0],
            normal: [0.0, 0.0, 2.0],
        };
        let n = s.normal_at([7.0, -4.0, 3.0]).unwrap();
        assert!(parallel_miss(n, [0.0, 0.0, 1.0]) < 1.0e-15);
    }

    #[test]
    fn cylinder_normal_is_radial_and_perpendicular_to_the_axis() {
        let d = oblique_axis();
        let (u, v) = frame(d);
        let s = RefSurface::Cylinder {
            axis_pos: [4.0, -1.0, 2.0],
            axis_dir: d,
            radius: 3.0,
        };
        for k in 0..16 {
            let t = k as f64 * 0.4;
            let a = k as f64 * 0.37;
            // A point on the surface: axis point + axial offset + 3 * radial.
            let radial = add(scale(u, a.cos()), scale(v, a.sin()));
            let p = add(add([4.0, -1.0, 2.0], scale(d, t)), scale(radial, 3.0));
            let n = s.normal_at(p).unwrap();
            assert!(parallel_miss(n, radial) < 1.0e-15, "not radial at k={k}");
            assert!(
                (dot(n, d)).abs() < 1.0e-12,
                "normal not perpendicular to the axis at k={k}"
            );
        }
    }

    #[test]
    fn cone_normal_tilts_by_the_half_angle_and_reduces_to_the_cylinder() {
        let d = oblique_axis();
        let (u, v) = frame(d);
        let apex = [0.5, 0.5, -1.0];
        let half = 0.4_f64;
        let s = RefSurface::Cone {
            apex,
            axis_dir: d,
            half_angle: half,
        };
        for k in 1..12 {
            let t = k as f64 * 0.9;
            let a = k as f64 * 0.53;
            let radial = add(scale(u, a.cos()), scale(v, a.sin()));
            let p = add(add(apex, scale(d, t)), scale(radial, t * half.tan()));
            let n = s.normal_at(p).unwrap();
            // The normal makes exactly (90° − half-angle) with the axis.
            let to_axis = acute_angle_deg(n, d);
            assert!(
                (to_axis - (90.0 - half.to_degrees())).abs() < 1.0e-9,
                "cone normal tilt wrong at k={k}: {to_axis}"
            );
            // And it is coplanar with the axis and the radial direction.
            assert!(dot(n, cross(d, radial)).abs() < 1.0e-12);
        }
        // half_angle = 0 is the cylinder case.
        let flat = RefSurface::Cone {
            apex,
            axis_dir: d,
            half_angle: 0.0,
        };
        let p = add(add(apex, scale(d, 2.0)), scale(u, 1.0));
        assert!(parallel_miss(flat.normal_at(p).unwrap(), u) < 1.0e-15);
    }

    #[test]
    fn sphere_normal_is_radial() {
        let s = RefSurface::Sphere {
            centre: [1.0, -2.0, 0.5],
            radius: 4.0,
        };
        let dir = unit([1.0, 2.0, -3.0]).unwrap();
        let p = add([1.0, -2.0, 0.5], scale(dir, 4.0));
        assert!(parallel_miss(s.normal_at(p).unwrap(), dir) < 1.0e-15);
    }

    #[test]
    fn torus_normal_line_meets_the_axis() {
        let d = oblique_axis();
        let (u, v) = frame(d);
        let centre = [2.0, 3.0, -1.0];
        let (major, minor) = (10.0, 2.5);
        let s = RefSurface::Torus {
            centre,
            axis_dir: d,
            major,
            minor,
        };
        for k in 0..20 {
            let a = k as f64 * 0.31; // around the axis
            let b = k as f64 * 0.77; // around the tube
            let meridian = add(scale(u, a.cos()), scale(v, a.sin()));
            let tube_centre = add(centre, scale(meridian, major));
            let p = add(
                tube_centre,
                scale(add(scale(meridian, b.cos()), scale(d, b.sin())), minor),
            );
            let n = s.normal_at(p).unwrap();
            // This is the meridian-coplanarity constraint Phase 1.1 solves for
            // the axis: the normal line and the axis are coplanar, so their
            // line-line distance is zero.
            let miss = line_line_distance(p, n, centre, d).unwrap_or(0.0);
            assert!(miss < 1.0e-9, "meridian residual {miss} at k={k}");
        }
    }

    /// Cross-validate the two normal paths against each other: the closed form
    /// used for analytic reference faces, and the `GeomLProp_SLProps`
    /// evaluation used for freeform ones. They must agree on a surface where
    /// both apply, or the freeform half is measuring something else.
    #[test]
    fn slprops_normal_agrees_with_the_closed_form() {
        use opencascade_sys::geom;
        let d = oblique_axis();
        let origin = [4.0, -1.0, 2.0];
        let radius = 3.0;

        let dir = gp::Dir::new_real3(d[0], d[1], d[2]);
        let loc = gp::Pnt::new_real3(origin[0], origin[1], origin[2]);
        let ax3 = gp::Ax3::new_pnt_dir(&loc, &dir);
        let cyl = geom::CylindricalSurface::new_ax3_real(&ax3, radius);
        let handle = geom::CylindricalSurface::to_handle(cyl).to_handle_surface();

        let closed_form = RefSurface::Cylinder {
            axis_pos: origin,
            axis_dir: d,
            radius,
        };

        for k in 0..12 {
            let (u, v) = (k as f64 * 0.5, k as f64 * 0.9 - 2.0);
            // The 3D point at those parameters, from OCCT itself.
            let mut props = geom_l_prop::SLProps::new_handlegeomsurface_real2_int_real(
                &handle, u, v, 1, 1.0e-12,
            );
            let p = props.value();
            let p = [p.x(), p.y(), p.z()];
            assert!(props.is_normal_defined(), "normal undefined at k={k}");
            let n = props.normal();
            let occt = [n.x(), n.y(), n.z()];

            let ours = closed_form.normal_at(p).unwrap();
            assert!(
                parallel_miss(occt, ours) < 1.0e-14,
                "closed form and SLProps disagree at k={k}: {:.3e}",
                parallel_miss(occt, ours)
            );
        }
    }

    /// The ruling theorem, in both directions: a facet on a cylinder has a
    /// normal exactly perpendicular to the axis **iff** two of its vertices
    /// share a ruling. Ruledness alone is not enough — a triangle straddling
    /// three rulings tilts its normal off perpendicular on a perfect cylinder.
    #[test]
    fn perpendicularity_needs_an_edge_on_a_ruling_not_just_a_ruled_surface() {
        let d = oblique_axis();
        let (u, v) = frame(d);
        let origin = [1.0, -2.0, 3.0];
        let r = 7.0;
        let on = |theta: f64, z: f64| {
            let radial = add(scale(u, theta.cos()), scale(v, theta.sin()));
            add(add(origin, scale(d, z)), scale(radial, r))
        };
        let normal_of = |a: [f64; 3], b: [f64; 3], c: [f64; 3]| {
            unit(cross(sub(b, a), sub(c, a))).unwrap()
        };

        // Two vertices on one ruling (same θ, different z): exactly perpendicular.
        for k in 0..8 {
            let t0 = k as f64 * 0.4;
            let n = normal_of(on(t0, 0.0), on(t0, 5.0), on(t0 + 0.3, 2.0));
            assert!(
                dot(n, d).abs() < 1.0e-15,
                "aligned facet {k} is not perpendicular: {:.3e}",
                dot(n, d).abs()
            );
        }

        // Three distinct rulings: NOT perpendicular, and measurably so. The
        // tilt grows with the angular spread and shrinks with the axial one,
        // which is why a straddling facet's error ranges from ~0 upward rather
        // than being large by default.
        let n = normal_of(on(0.0, 0.0), on(0.3, 5.0), on(0.6, 2.0));
        let tilt = (90.0 - acute_angle_deg(n, d)).abs();
        assert!(tilt > 0.5, "straddling facet should tilt measurably, got {tilt}°");

        // And the detector agrees with the geometry on both.
        let mesh = ConnectedMesh {
            vertices: vec![
                MeshVertex::from_xyz(on(0.0, 0.0)[0], on(0.0, 0.0)[1], on(0.0, 0.0)[2]),
                MeshVertex::from_xyz(on(0.0, 5.0)[0], on(0.0, 5.0)[1], on(0.0, 5.0)[2]),
                MeshVertex::from_xyz(on(0.3, 2.0)[0], on(0.3, 2.0)[1], on(0.3, 2.0)[2]),
                MeshVertex::from_xyz(on(0.6, 2.0)[0], on(0.6, 2.0)[1], on(0.6, 2.0)[2]),
            ],
            faces: vec![
                MeshFace {
                    vertex_count: 3,
                    vertex_indices: [0, 1, 2, 0],
                    neighbors: [-1; 4],
                    normal: None,
                },
                MeshFace {
                    vertex_count: 3,
                    vertex_indices: [0, 2, 3, 0],
                    neighbors: [-1; 4],
                    normal: None,
                },
            ],
            stats: Default::default(),
        };
        assert!(min_ruling_gap_rad(&mesh, 0, origin, d).unwrap() < RULING_SAME_RAD);
        assert!(min_ruling_gap_rad(&mesh, 1, origin, d).unwrap() > RULING_SAME_RAD);
    }

    #[test]
    fn circumcenter_of_a_right_triangle_is_the_hypotenuse_midpoint() {
        let mesh = ConnectedMesh {
            vertices: vec![
                MeshVertex::from_xyz(0.0, 0.0, 0.0),
                MeshVertex::from_xyz(4.0, 0.0, 0.0),
                MeshVertex::from_xyz(0.0, 6.0, 0.0),
            ],
            faces: vec![MeshFace {
                vertex_count: 3,
                vertex_indices: [0, 1, 2, 0],
                neighbors: [-1; 4],
                normal: Some([0.0, 0.0, 1.0]),
            }],
            stats: Default::default(),
        };
        let c = circumcenter(&mesh, 0).unwrap();
        assert!(norm(sub(c, [2.0, 3.0, 0.0])) < 1.0e-12, "got {c:?}");
    }

    #[test]
    fn common_perpendicular_midpoint_of_crossing_lines() {
        // Two lines that cross exactly at the origin.
        let m = common_perpendicular_midpoint(
            [5.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, -7.0, 0.0],
            [0.0, 1.0, 0.0],
        )
        .unwrap();
        assert!(norm(m) < 1.0e-12, "got {m:?}");
        // Two skew lines one unit apart: the midpoint sits half way.
        let m = common_perpendicular_midpoint(
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0],
            [0.0, 1.0, 0.0],
        )
        .unwrap();
        assert!(norm(sub(m, [0.0, 0.0, 0.5])) < 1.0e-12, "got {m:?}");
        // Parallel lines have no usable common perpendicular.
        assert!(common_perpendicular_midpoint(
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [1.0, 0.0, 0.0]
        )
        .is_none());
    }

    /// Build a strip-tessellated cylinder: one ring of quads spanning two
    /// rulings each, exactly as a B-rep faceter emits them.
    fn cylinder_strip(d: [f64; 3], origin: [f64; 3], radius: f64, steps: usize) -> ConnectedMesh {
        let (u, v) = frame(d);
        let mut vertices = Vec::new();
        for k in 0..steps {
            let a = 2.0 * std::f64::consts::PI * k as f64 / steps as f64;
            let radial = add(scale(u, a.cos()), scale(v, a.sin()));
            let base = add(origin, scale(radial, radius));
            vertices.push(MeshVertex::from_xyz(base[0], base[1], base[2]));
            let top = add(base, scale(d, 5.0));
            vertices.push(MeshVertex::from_xyz(top[0], top[1], top[2]));
        }
        let mut faces = Vec::new();
        for k in 0..steps {
            let k2 = (k + 1) % steps;
            let idx = [2 * k, 2 * k + 1, 2 * k2 + 1, 2 * k2];
            let p = |i: usize| {
                let vtx = &vertices[idx[i]];
                [vtx.x, vtx.y, vtx.z]
            };
            // Facet normal exactly as stage 1 computes it, from the geometry.
            let n = unit(cross(sub(p(1), p(0)), sub(p(2), p(0)))).unwrap();
            faces.push(MeshFace {
                vertex_count: 4,
                vertex_indices: idx,
                neighbors: [-1; 4],
                normal: Some(n),
            });
        }
        ConnectedMesh {
            vertices,
            faces,
            stats: Default::default(),
        }
    }

    #[test]
    fn strip_facet_normals_are_exactly_perpendicular_to_the_cylinder_axis() {
        let d = oblique_axis();
        let mesh = cylinder_strip(d, [1.0, -2.0, 3.0], 7.0, 24);
        for (fi, face) in mesh.faces.iter().enumerate() {
            let n = face.normal.unwrap();
            assert!(
                dot(n, d).abs() < 1.0e-14,
                "facet {fi} normal is not perpendicular to the axis: {:.3e}",
                dot(n, d).abs()
            );
        }
    }

    #[test]
    fn n1_cross_n2_recovers_the_cylinder_axis_exactly() {
        // The Phase 1.1 claim, on the tessellation it is claimed for.
        let d = oblique_axis();
        let mesh = cylinder_strip(d, [1.0, -2.0, 3.0], 7.0, 24);
        let mut worst = 0.0_f64;
        for i in 0..mesh.faces.len() {
            let j = (i + 1) % mesh.faces.len();
            let n1 = mesh.faces[i].normal.unwrap();
            let n2 = mesh.faces[j].normal.unwrap();
            let est = unit(cross(n1, n2)).unwrap();
            worst = worst.max(parallel_miss(est, d));
        }
        // Exact, not approximate: this is the property Phase 1.1 rests on. The
        // bar is a few ulps rather than one, because the estimate goes through
        // two cross products and two normalizations; measured worst is
        // 1.05e-15, i.e. ~5 ulps, which is 6e-14 degrees.
        assert!(worst < 1.0e-14, "worst axis miss {worst} (sin of the angle)");
    }

    #[test]
    fn sphere_circumcenter_normal_lines_are_concurrent_at_the_centre() {
        // The sphere initializer's claim: for a facet inscribed in a sphere the
        // line through its circumcenter along its normal passes exactly through
        // the centre.
        let centre = [1.0, -2.0, 0.5];
        let radius = 6.0;
        let d = oblique_axis();
        let (u, v) = frame(d);
        let mut vertices = Vec::new();
        let mut faces = Vec::new();
        for k in 0..12 {
            let a = k as f64 * 0.5;
            let b = 0.3 + k as f64 * 0.11;
            let dir = |a: f64, b: f64| {
                add(
                    scale(d, b.cos()),
                    scale(add(scale(u, a.cos()), scale(v, a.sin())), b.sin()),
                )
            };
            let base = vertices.len();
            for (da, db) in [(0.0, 0.0), (0.12, 0.0), (0.0, 0.1)] {
                let p = add(centre, scale(dir(a + da, b + db), radius));
                vertices.push(MeshVertex::from_xyz(p[0], p[1], p[2]));
            }
            let p = |i: usize| {
                let vtx = &vertices[base + i];
                [vtx.x, vtx.y, vtx.z]
            };
            let n = unit(cross(sub(p(1), p(0)), sub(p(2), p(0)))).unwrap();
            faces.push(MeshFace {
                vertex_count: 3,
                vertex_indices: [base, base + 1, base + 2, 0],
                neighbors: [-1; 4],
                normal: Some(n),
            });
        }
        let mesh = ConnectedMesh {
            vertices,
            faces,
            stats: Default::default(),
        };
        for fi in 0..mesh.faces.len() {
            let cc = circumcenter(&mesh, fi).unwrap();
            let n = mesh.faces[fi].normal.unwrap();
            let w = sub(centre, cc);
            let miss = norm(sub(w, scale(n, dot(w, n))));
            assert!(miss < 1.0e-9, "facet {fi} normal line misses by {miss}");
        }
    }
}
