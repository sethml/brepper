//! Stage 3: Surface Reconstruction
//!
//! 3.1: Create OCCT surface objects from hypotheses
//! 3.2: Detect and create tangency relationships
//! 3.3: Create OCCT edge wires (surface-surface intersections)
//! 3.4: Create OCCT faces (surface bounded by wires)
//! 3.5: Construct shells from connected faces
//! 3.6: Construct solids from shells

use crate::config::Config;
use crate::stage2::Stage2Output;
use std::error::Error;
use std::fmt::{Display, Formatter};

// ---------------------------------------------------------------------------
// Data structures — Stage 3 output
// ---------------------------------------------------------------------------

/// Descriptor for a reconstructed B-Rep face, corresponding to one selected surface hypothesis.
#[derive(Debug)]
pub struct FaceDescriptor {
    /// Index of the selected surface (into Stage2Output.selected_surfaces).
    pub selected_surface_idx: usize,
    // TODO: OCCT surface object (Geom_Plane, Geom_CylindricalSurface, etc.)
    // TODO: OCCT face object (TopoDS_Face)
    /// Indices of adjacent FaceDescriptors, ordered topologically (consecutive faces
    /// are adjacent to each other). A face connecting to this one twice will appear
    /// more than once.
    pub adjacent_faces: Vec<usize>,
    /// Edge wire indices, one per adjacent face. Edge wire i connects this face
    /// to adjacent face i.
    pub edge_wires: Vec<usize>,
    /// Vertex point indices, one per pair of adjacent faces. Vertex i is the
    /// intersection of this face, adjacent face i, and adjacent face (i+1)%N.
    /// Empty if there is only one adjacent face.
    pub vertex_points: Vec<usize>,
}

/// Descriptor for an edge wire between two adjacent faces.
#[derive(Debug)]
pub struct EdgeWire {
    /// Indices of the two adjacent FaceDescriptors.
    pub face_indices: [usize; 2],
    /// Indices of the two adjacent BRepVertices at the endpoints.
    pub vertex_indices: [usize; 2],
    // TODO: Vec<Geom_Curve> for 3D intersection curves (likely just one per edge)
    // TODO: For each adjacent face, Vec<Geom2d_Curve> for pcurves in UV-space
    /// Whether there is a tangency relationship between the two faces at this edge.
    pub tangent: bool,
}

/// Vertex where three or more faces meet.
#[derive(Debug)]
pub struct BRepVertex {
    /// Indices of adjacent FaceDescriptors, in topological order.
    /// May contain duplicates if a face connects to this vertex multiple ways.
    pub adjacent_faces: Vec<usize>,
    /// Indices of adjacent EdgeWires, in topological order.
    pub adjacent_wires: Vec<usize>,
    /// 3D point location of this vertex.
    pub point: [f64; 3],
    // TODO: Vec of 2D points in UV-space for each corresponding adjacent face
}

/// The output of Stage 3: the fully reconstructed B-Rep topology.
#[derive(Debug)]
pub struct Stage3Output {
    /// Reconstructed face descriptors, one per selected surface.
    pub face_descriptors: Vec<FaceDescriptor>,
    /// Edge wires connecting adjacent faces.
    pub edge_wires: Vec<EdgeWire>,
    /// Vertices where faces and edges meet.
    pub vertices: Vec<BRepVertex>,
    // TODO: OCCT shell and solid objects
}

// ---------------------------------------------------------------------------
// Errors
// ---------------------------------------------------------------------------

#[derive(Debug)]
pub enum Stage3Error {
    // TODO: specific error variants as stages are implemented
    NotImplemented(String),
}

impl Display for Stage3Error {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        match self {
            Stage3Error::NotImplemented(stage) => write!(f, "stage {stage} not yet implemented"),
        }
    }
}

impl Error for Stage3Error {}

// ---------------------------------------------------------------------------
// Stage 3 entry point
// ---------------------------------------------------------------------------

/// Run stage 3: reconstruct B-Rep surfaces, edges, and topology from fitted surfaces.
pub fn stage3(config: &Config, _input: Stage2Output) -> Result<Stage3Output, Stage3Error> {
    // Stage 3.1: Create OCCT surface objects
    // TODO: Convert each selected hypothesis into a Geom_* surface object
    if !config.quiet {
        eprintln!("Stage 3.1: Create OCCT surface objects (not yet implemented)");
    }

    if !config.stage.at_least(3, 2) {
        return Err(Stage3Error::NotImplemented("3.1".into()));
    }

    // Stage 3.2: Detect and create tangency relationships
    // TODO: Check surface normals along shared boundaries
    if !config.quiet {
        eprintln!("Stage 3.2: Detect tangency relationships (not yet implemented)");
    }

    if !config.stage.at_least(3, 3) {
        return Err(Stage3Error::NotImplemented("3.2".into()));
    }

    // Stage 3.3: Create OCCT edge wires
    // TODO: Compute surface-surface intersections, trim to vertices
    if !config.quiet {
        eprintln!("Stage 3.3: Create edge wires (not yet implemented)");
    }

    if !config.stage.at_least(3, 4) {
        return Err(Stage3Error::NotImplemented("3.3".into()));
    }

    // Stage 3.4: Create OCCT faces
    // TODO: Build BRepBuilderAPI_MakeFace from surface + bounding wires
    if !config.quiet {
        eprintln!("Stage 3.4: Create OCCT faces (not yet implemented)");
    }

    if !config.stage.at_least(3, 5) {
        return Err(Stage3Error::NotImplemented("3.4".into()));
    }

    // Stage 3.5: Construct shells
    // TODO: BRepBuilderAPI_Sewing to stitch faces
    if !config.quiet {
        eprintln!("Stage 3.5: Construct shells (not yet implemented)");
    }

    if !config.stage.at_least(3, 6) {
        return Err(Stage3Error::NotImplemented("3.5".into()));
    }

    // Stage 3.6: Construct solids
    // TODO: BRepBuilderAPI_MakeSolid from shells, classify voids
    if !config.quiet {
        eprintln!("Stage 3.6: Construct solids (not yet implemented)");
    }

    Err(Stage3Error::NotImplemented("3".into()))
}
