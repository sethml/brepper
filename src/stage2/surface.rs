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

pub(super) trait SurfaceType {
    /// The fitted surface parameters (e.g., axis, radius, center).
    type Fit: Clone;
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
