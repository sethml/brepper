//! Stage 2: Surface Fitting
//!
//! 2.1: Deduce planar hypotheses
//! 2.2: Deduce cylindrical hypotheses
//! 2.3: Deduce spherical hypotheses
//! 2.4: Deduce ruled surface hypotheses
//! 2.5: Deduce NURBS hypotheses
//! 2.6: Select surfaces for reconstruction

use crate::config::Config;
use crate::stage1::ConnectedMesh;
use std::error::Error;
use std::fmt::{Display, Formatter};

// ---------------------------------------------------------------------------
// Hypothesis data structures
// ---------------------------------------------------------------------------

/// A planar surface hypothesis fitted to a set of mesh faces.
#[derive(Debug, Clone)]
pub struct PlanarHypothesis {
    /// Unit normal vector pointing outward from the shell/solid.
    pub normal: [f64; 3],
    /// Signed distance from origin to plane along the normal.
    pub distance: f64,
    /// Set of mesh face indices that fit this hypothesis.
    pub faces: Vec<usize>,
    /// Set of mesh vertex indices on this plane.
    pub vertices: Vec<usize>,
    /// Maximum (positive) distance from any vertex to the plane.
    pub error_max: f64,
    /// Minimum (most negative) distance from any vertex to the plane.
    pub error_min: f64,
    /// Sum of absolute vertex-to-plane distances.
    pub error_abs_sum: f64,
}

/// A cylindrical surface hypothesis fitted to a set of mesh faces.
#[derive(Debug, Clone)]
pub struct CylindricalHypothesis {
    /// A point on the cylinder axis.
    pub axis_origin: [f64; 3],
    /// Unit direction vector along the cylinder axis.
    pub axis_direction: [f64; 3],
    /// Radius of the cylinder (always positive).
    pub radius: f64,
    /// Whether the surface normal points away from the axis (convex=true) or toward it (concave=false).
    pub convex: bool,
    /// Set of mesh face indices that fit this hypothesis.
    pub faces: Vec<usize>,
    /// Set of mesh vertex indices on this cylinder.
    pub vertices: Vec<usize>,
    /// Maximum absolute distance from any vertex to the cylinder surface.
    pub error_max: f64,
    /// Sum of absolute vertex-to-surface distances.
    pub error_abs_sum: f64,
}

/// A spherical surface hypothesis fitted to a set of mesh faces.
#[derive(Debug, Clone)]
pub struct SphericalHypothesis {
    /// Center of the sphere.
    pub center: [f64; 3],
    /// Radius of the sphere (always positive).
    pub radius: f64,
    /// Whether the surface normal points away from center (convex) or toward it (concave).
    pub convex: bool,
    /// Set of mesh face indices that fit this hypothesis.
    pub faces: Vec<usize>,
    /// Set of mesh vertex indices on this sphere.
    pub vertices: Vec<usize>,
    /// Maximum absolute distance from any vertex to the sphere surface.
    pub error_max: f64,
    /// Sum of absolute vertex-to-surface distances.
    pub error_abs_sum: f64,
}

// TODO: Stage 2.4 - Ruled surface hypothesis
// TODO: Stage 2.5 - NURBS hypothesis

/// Identifies a selected surface hypothesis by type and index.
#[derive(Debug, Clone, Copy)]
pub enum SelectedSurface {
    Planar(usize),
    Cylindrical(usize),
    Spherical(usize),
    // TODO: RuledSurface(usize),
    // TODO: Nurbs(usize),
}

// ---------------------------------------------------------------------------
// Stage 2 output
// ---------------------------------------------------------------------------

/// The output of Stage 2: the mesh with all hypotheses populated and surfaces selected.
#[derive(Debug)]
pub struct Stage2Output {
    /// The mesh from stage 1 with per-face hypothesis indices populated.
    pub mesh: ConnectedMesh,
    /// All planar hypotheses deduced in stage 2.1.
    pub planar_hypotheses: Vec<PlanarHypothesis>,
    /// All cylindrical hypotheses deduced in stage 2.2.
    pub cylindrical_hypotheses: Vec<CylindricalHypothesis>,
    /// All spherical hypotheses deduced in stage 2.3.
    pub spherical_hypotheses: Vec<SphericalHypothesis>,
    /// Surfaces selected in stage 2.6 for reconstruction. Each face should be
    /// covered by exactly one selected surface.
    pub selected_surfaces: Vec<SelectedSurface>,
}

// ---------------------------------------------------------------------------
// Errors
// ---------------------------------------------------------------------------

#[derive(Debug)]
pub enum Stage2Error {
    // TODO: specific error variants as stages are implemented
    NotImplemented(String),
}

impl Display for Stage2Error {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        match self {
            Stage2Error::NotImplemented(stage) => write!(f, "stage {stage} not yet implemented"),
        }
    }
}

impl Error for Stage2Error {}

// ---------------------------------------------------------------------------
// Stage 2 entry point
// ---------------------------------------------------------------------------

/// Run stage 2: fit surface hypotheses to mesh faces and select surfaces.
pub fn stage2(config: &Config, mesh: ConnectedMesh) -> Result<Stage2Output, Stage2Error> {
    let output = Stage2Output {
        mesh,
        planar_hypotheses: Vec::new(),
        cylindrical_hypotheses: Vec::new(),
        spherical_hypotheses: Vec::new(),
        selected_surfaces: Vec::new(),
    };

    // Stage 2.1: Deduce planar hypotheses
    // TODO: implement plane fitting via region growing from seed faces
    if !config.quiet {
        eprintln!("Stage 2.1: Deduce planar hypotheses (not yet implemented)");
    }

    if !config.stage.at_least(2, 2) {
        return Ok(output);
    }

    // Stage 2.2: Deduce cylindrical hypotheses
    // TODO: implement cylinder fitting
    if !config.quiet {
        eprintln!("Stage 2.2: Deduce cylindrical hypotheses (not yet implemented)");
    }

    if !config.stage.at_least(2, 3) {
        return Ok(output);
    }

    // Stage 2.3: Deduce spherical hypotheses
    // TODO: implement sphere fitting
    if !config.quiet {
        eprintln!("Stage 2.3: Deduce spherical hypotheses (not yet implemented)");
    }

    if !config.stage.at_least(2, 4) {
        return Ok(output);
    }

    // Stage 2.4: Deduce ruled surface hypotheses
    // TODO: optional - detect extruded curve surfaces
    if !config.quiet {
        eprintln!("Stage 2.4: Deduce ruled surface hypotheses (not yet implemented)");
    }

    if !config.stage.at_least(2, 5) {
        return Ok(output);
    }

    // Stage 2.5: Deduce NURBS hypotheses
    // TODO: fit NURBS to remaining ungrouped faces
    if !config.quiet {
        eprintln!("Stage 2.5: Deduce NURBS hypotheses (not yet implemented)");
    }

    if !config.stage.at_least(2, 6) {
        return Ok(output);
    }

    // Stage 2.6: Select surfaces for reconstruction
    // TODO: greedy selection of best-fitting hypotheses covering all faces
    if !config.quiet {
        eprintln!("Stage 2.6: Select surfaces (not yet implemented)");
    }

    Ok(output)
}
