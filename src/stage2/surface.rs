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
    type Fit: Clone + Send + 'static;
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
type FitMap = std::collections::HashMap<FitKey, Box<dyn std::any::Any + Send>>;

/// One stage-2 run's fit memo, shared by every thread working on that run.
#[derive(Clone)]
pub(super) struct FitCache(std::sync::Arc<std::sync::Mutex<FitMap>>);

thread_local! {
    /// The memo this thread fills and reads, while one is installed.
    static FIT_CACHE: std::cell::RefCell<Option<FitCache>> = const { std::cell::RefCell::new(None) };
}

/// Install `cache` on this thread for the duration of `f` (also on panic).
fn with_installed<R>(cache: Option<FitCache>, f: impl FnOnce() -> R) -> R {
    struct Restore(Option<FitCache>);
    impl Drop for Restore {
        fn drop(&mut self) {
            let prev = self.0.take();
            FIT_CACHE.with(|c| *c.borrow_mut() = prev);
        }
    }
    let prev = FIT_CACHE.with(|c| c.replace(cache));
    let _restore = Restore(prev);
    f()
}

/// Run `f` with a fresh fit memo installed on this thread (if `enabled`), and
/// drop it afterwards.
///
/// **Why it is safe:** `fit_surface` is a deterministic function of the mesh,
/// the face list and the vertex list (verbosity only gates printing), and the
/// mesh is borrowed immutably for the whole of stage 2.3–2.4. So a stored fit
/// for an identical `(type, faces, vertices)` key is exactly what recomputing
/// would return — every result is bit-identical, whichever thread computed it.
/// The memo is scoped to one stage-2 run so face indices from different meshes
/// can never collide.
///
/// **Why it pays:** competition reaches the same patch from different seeds and
/// flood starts. Measured with a probe, repeated torus fits were 33 s of
/// `rounded_cube_10_r2_fine`'s 68 s of torus fitting and 30 s of
/// `pipe_elbow_10_fine`'s 124 s.
pub(super) fn with_fit_cache<R>(enabled: bool, f: impl FnOnce() -> R) -> R {
    let cache = enabled.then(|| FitCache(Default::default()));
    with_installed(cache, f)
}

/// The memo installed on this thread, to hand to worker threads.
pub(super) fn current_fit_cache() -> Option<FitCache> {
    FIT_CACHE.with(|c| c.borrow().clone())
}

/// Run `f` on a worker thread with the run's memo (from [`current_fit_cache`]).
pub(super) fn with_shared_fit_cache<R>(cache: Option<FitCache>, f: impl FnOnce() -> R) -> R {
    with_installed(cache, f)
}

/// `S::fit_surface`, memoized while a fit memo is installed and
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
    let Some(cache) = current_fit_cache() else {
        return S::fit_surface(mesh, face_list, vertex_list, verbosity);
    };
    let key: FitKey = (std::any::TypeId::of::<S>(), face_list.to_vec(), vertex_list.to_vec());
    let hit = cache.0.lock().expect("fit memo").get(&key).map(|v| {
        v.downcast_ref::<Option<S::Fit>>().expect("keyed by type").clone()
    });
    if let Some(found) = hit {
        return found;
    }
    // Computed outside the lock, so threads never wait on each other's fits.
    // Two threads may compute the same key at once; they produce the same value.
    let fit = S::fit_surface(mesh, face_list, vertex_list, verbosity);
    cache.0.lock().expect("fit memo").insert(key, Box::new(fit.clone()));
    fit
}
