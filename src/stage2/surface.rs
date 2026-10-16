//! The `SurfaceType` trait: the shared abstraction over the analytic surface
//! fitters (planar, conical, spherical, toroidal).
//!
//! Each fitter implements `fit_surface` + the validation/convexity/hypothesis
//! methods; the stage-2.3 whole-region ladder (`classify.rs`) and scoped
//! competition (`compete.rs`) drive them generically via this trait. The
//! former BFS hypothesis-generation engine that also used this trait was
//! removed with the stage-2 rewrite (see DEVELOPMENT_PLAN Rewrite Roadmap
//! step 7).

use crate::stage1::ConnectedMesh;

pub(super) trait SurfaceType: 'static {
    /// Whether [`cached_fit`] should memoize this type's fits. Worth it only
    /// where a fit is expensive enough to dwarf a hash lookup on its inputs —
    /// the cone and torus, whose fits are Levenberg–Marquardt runs.
    const CACHE_FITS: bool = false;

    /// The fitted surface parameters (e.g., axis, radius, center).
    type Fit: Clone + 'static;
    /// The output hypothesis struct.
    type Hypothesis;

    /// Fit a surface to the given faces and vertices.
    /// Returns None if fitting fails.
    fn fit_surface(
        mesh: &ConnectedMesh,
        face_list: &[usize],
        vertex_list: &[usize],
        verbosity: u8,
    ) -> Option<Self::Fit>;

    /// Compute the absolute distance from a vertex position to the fitted surface.
    fn vertex_distance(fit: &Self::Fit, pos: &[f64; 3]) -> f64;

    /// Compute the absolute distance from a centroid to the fitted surface.
    /// Default: same as vertex_distance.
    fn centroid_distance(fit: &Self::Fit, centroid: &[f64; 3]) -> f64 {
        Self::vertex_distance(fit, centroid)
    }

    /// Check convexity consistency: does `fi` match the convexity of the fit?
    /// Return true to accept, false to reject.
    /// Default implementation always accepts (e.g., for planar surfaces).
    fn convexity_check(_fit: &Self::Fit, _mesh: &ConnectedMesh, _fi: usize) -> bool {
        true
    }

    /// Determine convexity from a set of faces and update the fit.
    /// Default: no-op (for surface types without convexity).
    fn determine_convexity(_mesh: &ConnectedMesh, _face_list: &[usize], _fit: &mut Self::Fit) {}

    /// Construct the final hypothesis from the fitted surface and its faces.
    fn make_hypothesis(
        mesh: &ConnectedMesh,
        faces: Vec<usize>,
        vertices: Vec<usize>,
        fit: &Self::Fit,
    ) -> Self::Hypothesis;
}

// ---------------------------------------------------------------------------
// Run-scoped fit memoization
// ---------------------------------------------------------------------------

/// Cache key: the surface type and its exact fit inputs.
type FitKey = (std::any::TypeId, Vec<usize>, Vec<usize>);

thread_local! {
    /// `Some` while [`with_fit_cache`] is active on this thread.
    static FIT_CACHE: std::cell::RefCell<Option<std::collections::HashMap<FitKey, Box<dyn std::any::Any>>>> =
        const { std::cell::RefCell::new(None) };
}

/// Run `f` with fit memoization enabled for this thread, and drop the cache
/// afterwards (also on panic).
///
/// **Why it is safe:** `fit_surface` is a deterministic function of the mesh,
/// the face list and the vertex list (verbosity only gates printing), and the
/// mesh is borrowed immutably for the whole of stage 2.3–2.4. So a stored fit
/// for an identical `(type, faces, vertices)` key is exactly what recomputing
/// would return — every result is bit-identical. The cache is scoped to one
/// stage-2 run so face indices from different meshes can never collide, and it
/// is thread-local because concurrent runs (parallel tests) each own a mesh.
///
/// **Why it pays:** competition reaches the same patch from different seeds and
/// flood starts. Measured with a probe, repeated torus fits were 33 s of
/// `rounded_cube_10_r2_fine`'s 68 s of torus fitting and 30 s of
/// `pipe_elbow_10_fine`'s 124 s.
pub(super) fn with_fit_cache<R>(enabled: bool, f: impl FnOnce() -> R) -> R {
    struct Reset;
    impl Drop for Reset {
        fn drop(&mut self) {
            FIT_CACHE.with(|c| *c.borrow_mut() = None);
        }
    }
    if !enabled {
        return f();
    }
    FIT_CACHE.with(|c| *c.borrow_mut() = Some(std::collections::HashMap::new()));
    let _reset = Reset;
    f()
}

/// `S::fit_surface`, memoized while [`with_fit_cache`] is active and
/// `S::CACHE_FITS` is set; otherwise a plain call.
pub(super) fn cached_fit<S>(
    mesh: &ConnectedMesh,
    face_list: &[usize],
    vertex_list: &[usize],
    verbosity: u8,
) -> Option<S::Fit>
where
    S: SurfaceType,
{
    if !S::CACHE_FITS {
        return S::fit_surface(mesh, face_list, vertex_list, verbosity);
    }
    let key: FitKey = (std::any::TypeId::of::<S>(), face_list.to_vec(), vertex_list.to_vec());
    let hit = FIT_CACHE.with(|c| {
        c.borrow().as_ref().map(|m| {
            m.get(&key)
                .map(|v| v.downcast_ref::<Option<S::Fit>>().expect("keyed by type").clone())
        })
    });
    match hit {
        Some(Some(found)) => found,
        // Cache active, miss: compute and store (failures too — a failed fit
        // is as deterministic as a successful one).
        Some(None) => {
            let fit = S::fit_surface(mesh, face_list, vertex_list, verbosity);
            FIT_CACHE.with(|c| {
                if let Some(m) = c.borrow_mut().as_mut() {
                    m.insert(key, Box::new(fit.clone()));
                }
            });
            fit
        }
        // Cache inactive.
        None => S::fit_surface(mesh, face_list, vertex_list, verbosity),
    }
}
