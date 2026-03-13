//! Interactive 3D visualization for debugging brepper stages.
//!
//! Shared rendering infrastructure used by both `viewer` and `brepper --viz-stages`.
//!
//! # Architecture
//!
//! The brepper pipeline runs in a background thread. The main thread owns the
//! window and drives the render loop. Communication uses channels:
//!
//! - Pipeline thread sends `VizOverlay` via `VizSender::show_and_wait()`
//! - Render loop displays the overlay and captures user input
//! - User presses Space (next step) or Shift+Space (next seed) or Q (quit)
//! - Render loop sends `VizAction` back to the pipeline thread
//!
//! When `--viz-stages` is not set, `VizSender` is `None` and stages skip all
//! visualization calls with zero overhead.

use std::collections::HashMap;
use std::sync::mpsc;
use three_d::*;

// ---------------------------------------------------------------------------
// MeshData: triangulated mesh ready for GPU rendering
// ---------------------------------------------------------------------------

/// Triangulated mesh data ready for GPU upload.
pub struct MeshData {
    pub positions: Vec<Vector3<f32>>,
    pub normals: Vec<Vector3<f32>>,
    pub indices: Vec<u32>,
    /// Wireframe edge positions as pairs [p0, p1, p0, p1, ...].
    pub edge_positions: Vec<Vector3<f32>>,
    /// Near-coplanar edges (drawn at reduced alpha).
    pub soft_edge_positions: Vec<Vector3<f32>>,
}

pub const PALETTE: &[(f32, f32, f32)] = &[
    (0.55, 0.65, 0.80),
    (0.80, 0.55, 0.55),
    (0.55, 0.78, 0.55),
    (0.78, 0.70, 0.50),
    (0.70, 0.55, 0.78),
    (0.55, 0.75, 0.75),
];

// ---------------------------------------------------------------------------
// Feature edge extraction
// ---------------------------------------------------------------------------

/// Extract boundary/sharp edges and soft (near-coplanar) edges.
/// Returns (sharp_edges, soft_edges) as position pairs.
pub fn extract_feature_edges(
    positions: &[Vector3<f32>],
    indices: &[u32],
) -> (Vec<Vector3<f32>>, Vec<Vector3<f32>>) {
    let mut edge_faces: HashMap<(u32, u32), Vec<Vector3<f32>>> = HashMap::new();

    for tri in indices.chunks(3) {
        if tri.len() < 3 {
            continue;
        }
        let (i0, i1, i2) = (tri[0], tri[1], tri[2]);
        let v0 = positions[i0 as usize];
        let v1 = positions[i1 as usize];
        let v2 = positions[i2 as usize];
        let normal = (v1 - v0).cross(v2 - v0);
        let len = normal.magnitude();
        let normal = if len > 1e-10 { normal / len } else { normal };

        for &(a, b) in &[(i0, i1), (i1, i2), (i2, i0)] {
            let key = if a < b { (a, b) } else { (b, a) };
            edge_faces.entry(key).or_default().push(normal);
        }
    }

    let mut sharp = Vec::new();
    let mut soft = Vec::new();
    let cos_threshold = 0.7_f32; // ~45 deg dihedral

    for ((a, b), normals) in &edge_faces {
        let p0 = positions[*a as usize];
        let p1 = positions[*b as usize];
        if normals.len() == 1 {
            sharp.push(p0);
            sharp.push(p1);
        } else {
            let cos_angle = normals[0].dot(normals[1]);
            if cos_angle < cos_threshold {
                sharp.push(p0);
                sharp.push(p1);
            } else {
                soft.push(p0);
                soft.push(p1);
            }
        }
    }
    (sharp, soft)
}

/// Extract feature edges using face definitions (quad-aware).
/// Quad internal diagonals are excluded, only face boundary edges are considered.
pub fn extract_feature_edges_from_faces(
    positions: &[Vector3<f32>],
    face_indices: &[[usize; 4]],
    face_vertex_counts: &[u8],
) -> (Vec<Vector3<f32>>, Vec<Vector3<f32>>) {
    // Build edge → face normals map using actual face edges (not triangulated)
    let mut edge_faces: HashMap<(usize, usize), Vec<Vector3<f32>>> = HashMap::new();

    for (fi, vis) in face_indices.iter().enumerate() {
        let vc = face_vertex_counts[fi] as usize;
        if vc < 3 {
            continue;
        }
        // Compute face normal via Newell's method
        let v0 = positions[vis[0]];
        let v1 = positions[vis[1]];
        let v2 = positions[vis[2]];
        let normal = (v1 - v0).cross(v2 - v0);
        let len = normal.magnitude();
        let normal = if len > 1e-10 { normal / len } else { normal };

        for j in 0..vc {
            let a = vis[j];
            let b = vis[(j + 1) % vc];
            let key = if a < b { (a, b) } else { (b, a) };
            edge_faces.entry(key).or_default().push(normal);
        }
    }

    let mut sharp = Vec::new();
    let mut soft = Vec::new();
    let cos_threshold = 0.7_f32; // ~45 deg dihedral

    for ((a, b), normals) in &edge_faces {
        let p0 = positions[*a];
        let p1 = positions[*b];
        if normals.len() == 1 {
            sharp.push(p0);
            sharp.push(p1);
        } else {
            let cos_angle = normals[0].dot(normals[1]);
            if cos_angle < cos_threshold {
                sharp.push(p0);
                sharp.push(p1);
            } else {
                soft.push(p0);
                soft.push(p1);
            }
        }
    }
    (sharp, soft)
}

/// Smooth step interpolation (ease in/out).
fn smooth_step(t: f32) -> f32 {
    t * t * (3.0 - 2.0 * t)
}

// ---------------------------------------------------------------------------
// Bounding box
// ---------------------------------------------------------------------------

pub fn compute_bounding_box(meshes: &[&MeshData]) -> (Vector3<f32>, Vector3<f32>) {
    let mut min = vec3(f32::MAX, f32::MAX, f32::MAX);
    let mut max = vec3(f32::MIN, f32::MIN, f32::MIN);
    for mesh in meshes {
        for p in &mesh.positions {
            min.x = min.x.min(p.x);
            min.y = min.y.min(p.y);
            min.z = min.z.min(p.z);
            max.x = max.x.max(p.x);
            max.y = max.y.max(p.y);
            max.z = max.z.max(p.z);
        }
    }
    (min, max)
}

// ---------------------------------------------------------------------------
// Raw GL wireframe renderer
// ---------------------------------------------------------------------------

const LINE_VERT: &str = r#"
#version 330 core
layout(location = 0) in vec3 position;
uniform mat4 mvp;
void main() {
    gl_Position = mvp * vec4(position, 1.0);
}
"#;

const LINE_FRAG: &str = r#"
#version 330 core
uniform vec4 lineColor;
out vec4 fragColor;
void main() {
    fragColor = lineColor;
}
"#;

pub struct WireframeRenderer {
    program: context::Program,
    vbo: context::Buffer,
    vertex_count: i32,
}

impl WireframeRenderer {
    pub fn new(ctx: &Context, positions: &[Vector3<f32>]) -> Self {
        use context::HasContext;
        unsafe {
            let program = ctx.create_program().unwrap();

            let vs = ctx.create_shader(context::VERTEX_SHADER).unwrap();
            ctx.shader_source(vs, LINE_VERT);
            ctx.compile_shader(vs);
            assert!(
                ctx.get_shader_compile_status(vs),
                "line VS compile failed: {}",
                ctx.get_shader_info_log(vs)
            );

            let fs = ctx.create_shader(context::FRAGMENT_SHADER).unwrap();
            ctx.shader_source(fs, LINE_FRAG);
            ctx.compile_shader(fs);
            assert!(
                ctx.get_shader_compile_status(fs),
                "line FS compile failed: {}",
                ctx.get_shader_info_log(fs)
            );

            ctx.attach_shader(program, vs);
            ctx.attach_shader(program, fs);
            ctx.link_program(program);
            assert!(
                ctx.get_program_link_status(program),
                "line program link failed: {}",
                ctx.get_program_info_log(program)
            );
            ctx.delete_shader(vs);
            ctx.delete_shader(fs);

            let vbo = ctx.create_buffer().unwrap();
            ctx.bind_buffer(context::ARRAY_BUFFER, Some(vbo));
            let float_data: Vec<f32> =
                positions.iter().flat_map(|v| [v.x, v.y, v.z]).collect();
            let byte_data: &[u8] = std::slice::from_raw_parts(
                float_data.as_ptr() as *const u8,
                float_data.len() * std::mem::size_of::<f32>(),
            );
            ctx.buffer_data_u8_slice(context::ARRAY_BUFFER, byte_data, context::STATIC_DRAW);
            ctx.bind_buffer(context::ARRAY_BUFFER, None);

            WireframeRenderer {
                program,
                vbo,
                vertex_count: positions.len() as i32,
            }
        }
    }

    pub fn render(&self, ctx: &Context, camera: &Camera, color: [f32; 4], depth_test: bool) {
        if self.vertex_count == 0 {
            return;
        }
        use context::HasContext;
        unsafe {
            let mvp = camera.projection() * camera.view();
            let mvp_ref: &[f32; 16] = mvp.as_ref();

            ctx.use_program(Some(self.program));

            let mvp_loc = ctx.get_uniform_location(self.program, "mvp");
            ctx.uniform_matrix_4_f32_slice(mvp_loc.as_ref(), false, mvp_ref);

            let color_loc = ctx.get_uniform_location(self.program, "lineColor");
            ctx.uniform_4_f32(
                color_loc.as_ref(),
                color[0],
                color[1],
                color[2],
                color[3],
            );

            let vao = ctx.create_vertex_array().unwrap();
            ctx.bind_vertex_array(Some(vao));

            ctx.bind_buffer(context::ARRAY_BUFFER, Some(self.vbo));
            ctx.enable_vertex_attrib_array(0);
            ctx.vertex_attrib_pointer_f32(0, 3, context::FLOAT, false, 12, 0);

            ctx.line_width(1.5);
            if depth_test {
                ctx.enable(context::DEPTH_TEST);
                ctx.depth_func(context::LEQUAL);
            } else {
                ctx.disable(context::DEPTH_TEST);
            }
            ctx.draw_arrays(context::LINES, 0, self.vertex_count);

            ctx.enable(context::DEPTH_TEST);
            ctx.depth_func(context::LESS);

            ctx.bind_vertex_array(None);
            ctx.delete_vertex_array(vao);
            ctx.use_program(None);

            ctx.bind_buffer(context::ARRAY_BUFFER, None);
        }
    }
}

// ---------------------------------------------------------------------------
// Cylinder / sphere mesh generation (for translucent overlays)
// ---------------------------------------------------------------------------

/// Generate a triangulated cylinder mesh for visualization.
fn generate_cylinder_mesh(
    origin: [f64; 3],
    direction: [f64; 3],
    radius: f64,
    half_length: f64,
    segments: usize,
) -> (Vec<Vector3<f32>>, Vec<Vector3<f32>>, Vec<u32>) {
    let d = direction;
    let arbitrary = if d[0].abs() < 0.9 {
        [1.0, 0.0, 0.0]
    } else {
        [0.0, 1.0, 0.0]
    };
    let u = [
        d[1] * arbitrary[2] - d[2] * arbitrary[1],
        d[2] * arbitrary[0] - d[0] * arbitrary[2],
        d[0] * arbitrary[1] - d[1] * arbitrary[0],
    ];
    let u_len = (u[0] * u[0] + u[1] * u[1] + u[2] * u[2]).sqrt();
    let u = [u[0] / u_len, u[1] / u_len, u[2] / u_len];
    let w = [
        d[1] * u[2] - d[2] * u[1],
        d[2] * u[0] - d[0] * u[2],
        d[0] * u[1] - d[1] * u[0],
    ];

    let mut positions = Vec::new();
    let mut normals = Vec::new();
    let mut indices = Vec::new();

    for ring in 0..2usize {
        let t = if ring == 0 { -half_length } else { half_length };
        for seg in 0..segments {
            let angle = 2.0 * std::f64::consts::PI * (seg as f64) / (segments as f64);
            let cos_a = angle.cos();
            let sin_a = angle.sin();

            let px = origin[0] + d[0] * t + radius * (u[0] * cos_a + w[0] * sin_a);
            let py = origin[1] + d[1] * t + radius * (u[1] * cos_a + w[1] * sin_a);
            let pz = origin[2] + d[2] * t + radius * (u[2] * cos_a + w[2] * sin_a);
            positions.push(vec3(px as f32, py as f32, pz as f32));

            let nx = u[0] * cos_a + w[0] * sin_a;
            let ny = u[1] * cos_a + w[1] * sin_a;
            let nz = u[2] * cos_a + w[2] * sin_a;
            normals.push(vec3(nx as f32, ny as f32, nz as f32));
        }
    }

    for seg in 0..segments {
        let next = (seg + 1) % segments;
        let i0 = seg as u32;
        let i1 = next as u32;
        let i2 = (segments + seg) as u32;
        let i3 = (segments + next) as u32;
        indices.push(i0);
        indices.push(i1);
        indices.push(i2);
        indices.push(i1);
        indices.push(i3);
        indices.push(i2);
    }

    (positions, normals, indices)
}

fn generate_sphere_mesh(
    center: [f64; 3],
    radius: f64,
    lat_segments: usize,
    lon_segments: usize,
) -> (Vec<Vector3<f32>>, Vec<Vector3<f32>>, Vec<u32>) {
    let mut positions = Vec::new();
    let mut normals = Vec::new();
    let mut indices = Vec::new();

    for lat in 0..=lat_segments {
        let theta = std::f64::consts::PI * (lat as f64) / (lat_segments as f64);
        let sin_t = theta.sin();
        let cos_t = theta.cos();
        for lon in 0..=lon_segments {
            let phi = 2.0 * std::f64::consts::PI * (lon as f64) / (lon_segments as f64);
            let nx = sin_t * phi.cos();
            let ny = sin_t * phi.sin();
            let nz = cos_t;
            positions.push(vec3(
                (center[0] + radius * nx) as f32,
                (center[1] + radius * ny) as f32,
                (center[2] + radius * nz) as f32,
            ));
            normals.push(vec3(nx as f32, ny as f32, nz as f32));
        }
    }

    for lat in 0..lat_segments {
        for lon in 0..lon_segments {
            let i0 = (lat * (lon_segments + 1) + lon) as u32;
            let i1 = i0 + 1;
            let i2 = i0 + (lon_segments + 1) as u32;
            let i3 = i2 + 1;
            indices.push(i0);
            indices.push(i2);
            indices.push(i1);
            indices.push(i1);
            indices.push(i2);
            indices.push(i3);
        }
    }

    (positions, normals, indices)
}

// ---------------------------------------------------------------------------
// VizAction / VizOverlay: communication types
// ---------------------------------------------------------------------------

/// Actions from user input during interactive visualization.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum VizAction {
    /// Advance to the next step (space bar).
    NextStep,
    /// Skip to the next seed/hypothesis (shift+space).
    NextSeed,
    /// User wants to quit.
    Quit,
}

/// A face highlight overlay: a set of mesh face indices drawn in a specific color.
pub struct FaceHighlight {
    pub face_indices: Vec<usize>,
    pub color: [f32; 4],
}

/// Edge outline overlay: edges of specified faces drawn as thick colored lines.
pub struct EdgeHighlight {
    pub face_indices: Vec<usize>,
    pub color: [f32; 4],
}

/// A translucent cylinder overlay for visualization.
pub struct CylinderOverlay {
    pub origin: [f64; 3],
    pub direction: [f64; 3],
    pub radius: f64,
    pub half_length: f64,
    pub color: [f32; 4],
}

/// A translucent sphere overlay for visualization.
pub struct SphereOverlay {
    pub center: [f64; 3],
    pub radius: f64,
    pub color: [f32; 4],
}

/// Arbitrary 3D line segments drawn in a specific color.
pub struct LineOverlay {
    /// Pairs of points: [p0, p1, p0, p1, ...]
    pub positions: Vec<[f32; 3]>,
    pub color: [f32; 4],
    /// If true, lines are always visible (no depth test).
    pub no_depth_test: bool,
}

/// A pre-tessellated 3D mesh overlay (for rendering OCCT shapes).
pub struct ShapeMeshOverlay {
    pub positions: Vec<[f32; 3]>,
    pub normals: Vec<[f32; 3]>,
    pub indices: Vec<u32>,
    pub color: [f32; 4],
    /// Edge wireframe positions (pairs).
    pub edge_positions: Vec<[f32; 3]>,
    pub edge_color: [f32; 4],
}

/// A visualization step: overlay to display while waiting for user input.
pub struct VizOverlay {
    pub face_highlights: Vec<FaceHighlight>,
    pub edge_highlights: Vec<EdgeHighlight>,
    pub cylinders: Vec<CylinderOverlay>,
    pub spheres: Vec<SphereOverlay>,
    pub lines: Vec<LineOverlay>,
    pub shape_meshes: Vec<ShapeMeshOverlay>,
    pub status_text: String,
    /// If set, the camera will smoothly move to center on this point.
    pub focus_point: Option<[f32; 3]>,
    /// If set (alongside focus_point), the camera will rotate to look along this face normal.
    pub focus_normal: Option<[f32; 3]>,
}

impl Default for VizOverlay {
    fn default() -> Self {
        Self::new()
    }
}

impl VizOverlay {
    pub fn new() -> Self {
        Self {
            face_highlights: Vec::new(),
            edge_highlights: Vec::new(),
            cylinders: Vec::new(),
            spheres: Vec::new(),
            lines: Vec::new(),
            shape_meshes: Vec::new(),
            status_text: String::new(),
            focus_point: None,
            focus_normal: None,
        }
    }
}

// ---------------------------------------------------------------------------
// VizSender: pipeline-side handle for sending overlays
// ---------------------------------------------------------------------------

/// Pipeline-side handle for sending visualization overlays and waiting for
/// the user's response (space / shift-space / quit).
pub struct VizSender {
    overlay_tx: mpsc::SyncSender<VizOverlay>,
    action_rx: mpsc::Receiver<VizAction>,
}

impl VizSender {
    /// Send an overlay to the render thread and block until the user responds.
    /// Returns the user's action (NextStep, NextSeed, or Quit).
    pub fn show_and_wait(&self, overlay: VizOverlay) -> VizAction {
        if !overlay.status_text.is_empty() {
            eprintln!("[viz] {}", overlay.status_text);
        }
        if self.overlay_tx.send(overlay).is_err() {
            return VizAction::Quit; // render thread gone
        }
        self.action_rx.recv().unwrap_or(VizAction::Quit)
    }
}

// ---------------------------------------------------------------------------
// VizReceiver: render-side handle
// ---------------------------------------------------------------------------

/// Render-side handle for receiving overlays from the pipeline thread.
struct VizReceiver {
    overlay_rx: mpsc::Receiver<VizOverlay>,
    action_tx: mpsc::SyncSender<VizAction>,
}

// ---------------------------------------------------------------------------
// Channel constructor
// ---------------------------------------------------------------------------

/// Create a viz channel pair: (sender for pipeline thread, receiver for render thread).
fn viz_channel() -> (VizSender, VizReceiver) {
    let (overlay_tx, overlay_rx) = mpsc::sync_channel(0);
    let (action_tx, action_rx) = mpsc::sync_channel(0);
    (
        VizSender {
            overlay_tx,
            action_rx,
        },
        VizReceiver {
            overlay_rx,
            action_tx,
        },
    )
}

// ---------------------------------------------------------------------------
// run_viz_window: main-thread render loop for interactive debugging
// ---------------------------------------------------------------------------

/// Data needed to set up the visualization window. Built on the main thread
/// before spawning the pipeline thread.
pub struct VizSetup {
    pub base_mesh: MeshData,
    pub compare_mesh: Option<MeshData>,
    pub face_indices: Vec<[usize; 4]>,
    pub face_vertex_counts: Vec<u8>,
}

// ---------------------------------------------------------------------------
// Free functions for building overlay GPU objects
// ---------------------------------------------------------------------------

fn build_face_edge_wireframe(
    ctx: &Context,
    face_idxs: &[usize],
    stl_positions: &[Vector3<f32>],
    face_indices_arr: &[[usize; 4]],
    face_vc_arr: &[u8],
) -> WireframeRenderer {
    let mut positions = Vec::new();
    for &fi in face_idxs {
        let vc = face_vc_arr[fi] as usize;
        let vis = &face_indices_arr[fi];
        for j in 0..vc {
            let j_next = (j + 1) % vc;
            positions.push(stl_positions[vis[j]]);
            positions.push(stl_positions[vis[j_next]]);
        }
    }
    WireframeRenderer::new(ctx, &positions)
}

fn build_face_highlight_mesh(
    ctx: &Context,
    face_idxs: &[usize],
    color: [f32; 4],
    stl_positions: &[Vector3<f32>],
    face_indices_arr: &[[usize; 4]],
    face_vc_arr: &[u8],
) -> Gm<Mesh, PhysicalMaterial> {
    let mut positions = Vec::new();
    let mut normals = Vec::new();
    let mut indices = Vec::new();

    for &fi in face_idxs {
        let vc = face_vc_arr[fi] as usize;
        let vis = &face_indices_arr[fi];
        let base = positions.len() as u32;
        for j in 0..vc {
            positions.push(stl_positions[vis[j]]);
        }
        if vc >= 3 {
            let v0 = stl_positions[vis[0]];
            let v1 = stl_positions[vis[1]];
            let v2 = stl_positions[vis[2]];
            let n = (v1 - v0).cross(v2 - v0);
            let len = n.magnitude();
            let n = if len > 1e-10 { n / len } else { vec3(0.0, 0.0, 1.0) };
            for _ in 0..vc {
                normals.push(n);
            }
            indices.push(base);
            indices.push(base + 1);
            indices.push(base + 2);
        }
        if vc >= 4 {
            indices.push(base);
            indices.push(base + 2);
            indices.push(base + 3);
        }
    }

    let cpu_mesh = CpuMesh {
        positions: Positions::F32(positions),
        indices: Indices::U32(indices),
        normals: Some(normals),
        ..Default::default()
    };
    let albedo = Srgba::new(
        (color[0] * 255.0) as u8,
        (color[1] * 255.0) as u8,
        (color[2] * 255.0) as u8,
        (color[3] * 255.0) as u8,
    );
    let is_transparent = color[3] < 0.99;
    let material = if is_transparent {
        PhysicalMaterial::new_transparent(
            ctx,
            &CpuMaterial {
                albedo,
                roughness: 0.4,
                metallic: 0.0,
                ..Default::default()
            },
        )
    } else {
        let mut m = PhysicalMaterial::new_opaque(
            ctx,
            &CpuMaterial {
                albedo,
                roughness: 0.4,
                metallic: 0.0,
                ..Default::default()
            },
        );
        m.render_states.cull = Cull::Back;
        m
    };
    Gm::new(Mesh::new(ctx, &cpu_mesh), material)
}

fn build_cylinder_overlay_mesh(
    ctx: &Context,
    ov: &CylinderOverlay,
) -> Gm<Mesh, PhysicalMaterial> {
    let (pos, nor, idx) =
        generate_cylinder_mesh(ov.origin, ov.direction, ov.radius, ov.half_length, 48);
    let cpu_mesh = CpuMesh {
        positions: Positions::F32(pos),
        indices: Indices::U32(idx),
        normals: Some(nor),
        ..Default::default()
    };
    let albedo = Srgba::new(
        (ov.color[0] * 255.0) as u8,
        (ov.color[1] * 255.0) as u8,
        (ov.color[2] * 255.0) as u8,
        (ov.color[3] * 255.0) as u8,
    );
    Gm::new(
        Mesh::new(ctx, &cpu_mesh),
        PhysicalMaterial::new_transparent(
            ctx,
            &CpuMaterial {
                albedo,
                roughness: 0.3,
                metallic: 0.0,
                ..Default::default()
            },
        ),
    )
}

fn build_sphere_overlay_mesh(
    ctx: &Context,
    ov: &SphereOverlay,
) -> Gm<Mesh, PhysicalMaterial> {
    let (pos, nor, idx) = generate_sphere_mesh(ov.center, ov.radius, 24, 48);
    let cpu_mesh = CpuMesh {
        positions: Positions::F32(pos),
        indices: Indices::U32(idx),
        normals: Some(nor),
        ..Default::default()
    };
    let albedo = Srgba::new(
        (ov.color[0] * 255.0) as u8,
        (ov.color[1] * 255.0) as u8,
        (ov.color[2] * 255.0) as u8,
        (ov.color[3] * 255.0) as u8,
    );
    Gm::new(
        Mesh::new(ctx, &cpu_mesh),
        PhysicalMaterial::new_transparent(
            ctx,
            &CpuMaterial {
                albedo,
                roughness: 0.3,
                metallic: 0.0,
                ..Default::default()
            },
        ),
    )
}

fn build_shape_mesh_overlay(
    ctx: &Context,
    sm: &ShapeMeshOverlay,
) -> Gm<Mesh, PhysicalMaterial> {
    let positions: Vec<Vector3<f32>> = sm.positions.iter().map(|p| vec3(p[0], p[1], p[2])).collect();
    let normals: Vec<Vector3<f32>> = sm.normals.iter().map(|n| vec3(n[0], n[1], n[2])).collect();
    let cpu_mesh = CpuMesh {
        positions: Positions::F32(positions),
        indices: Indices::U32(sm.indices.clone()),
        normals: Some(normals),
        ..Default::default()
    };
    let albedo = Srgba::new(
        (sm.color[0] * 255.0) as u8,
        (sm.color[1] * 255.0) as u8,
        (sm.color[2] * 255.0) as u8,
        (sm.color[3] * 255.0) as u8,
    );
    let is_transparent = sm.color[3] < 0.99;
    let material = if is_transparent {
        PhysicalMaterial::new_transparent(
            ctx,
            &CpuMaterial {
                albedo,
                roughness: 0.4,
                metallic: 0.0,
                ..Default::default()
            },
        )
    } else {
        let mut m = PhysicalMaterial::new_opaque(
            ctx,
            &CpuMaterial {
                albedo,
                roughness: 0.4,
                metallic: 0.0,
                ..Default::default()
            },
        );
        m.render_states.cull = Cull::Back;
        m
    };
    Gm::new(Mesh::new(ctx, &cpu_mesh), material)
}

/// Start the visualization system. Returns a `VizSender` that should be passed
/// to the pipeline thread. This function blocks on the main thread, running
/// the window event loop.
///
/// `pipeline_fn` is called in a background thread with the `VizSender`. When
/// the pipeline finishes, the window closes.
pub fn run_viz_window<F>(setup: VizSetup, pipeline_fn: F)
where
    F: FnOnce(VizSender) + Send + 'static,
{
    let (sender, receiver) = viz_channel();

    // Spawn pipeline in background thread
    let pipeline_handle = std::thread::spawn(move || {
        pipeline_fn(sender);
    });

    // Build window and GL objects on main thread
    let refs: Vec<&MeshData> = std::iter::once(&setup.base_mesh)
        .chain(setup.compare_mesh.as_ref())
        .collect();
    let (bb_min, bb_max) = compute_bounding_box(&refs);
    let center = (bb_min + bb_max) * 0.5;
    let extent = (bb_max - bb_min).magnitude();

    let window = Window::new(WindowSettings {
        title: "brepper viz".to_string(),
        max_size: Some((1920, 1200)),
        ..Default::default()
    })
    .unwrap();
    let context = window.gl();

    let mut camera = Camera::new_perspective(
        window.viewport(),
        center + vec3(extent * 0.5, extent * 0.3, extent * 0.5),
        center,
        vec3(0.0, 0.0, 1.0),
        degrees(45.0),
        extent * 0.001,
        extent * 100.0,
    );
    let min_distance = extent * 0.01;
    let max_distance = extent * 10.0;
    let z_near = extent * 0.001;
    let z_far = extent * 100.0;
    let mut ortho_height = 2.0 * extent * (std::f32::consts::FRAC_PI_4 / 2.0).tan();

    // Build base mesh
    let cpu_mesh = CpuMesh {
        positions: Positions::F32(setup.base_mesh.positions.clone()),
        indices: Indices::U32(setup.base_mesh.indices.clone()),
        normals: Some(setup.base_mesh.normals.clone()),
        ..Default::default()
    };
    let (cr, cg, cb) = PALETTE[0];
    let mut base_material = PhysicalMaterial::new_opaque(
        &context,
        &CpuMaterial {
            albedo: Srgba::new(
                (cr * 255.0) as u8,
                (cg * 255.0) as u8,
                (cb * 255.0) as u8,
                255,
            ),
            roughness: 0.6,
            metallic: 0.1,
            ..Default::default()
        },
    );
    base_material.render_states.cull = Cull::Back;
    let base_mesh_object = Gm::new(Mesh::new(&context, &cpu_mesh), base_material);

    let base_wireframe = if !setup.base_mesh.edge_positions.is_empty() {
        Some(WireframeRenderer::new(
            &context,
            &setup.base_mesh.edge_positions,
        ))
    } else {
        None
    };
    let base_soft_wireframe = if !setup.base_mesh.soft_edge_positions.is_empty() {
        Some(WireframeRenderer::new(
            &context,
            &setup.base_mesh.soft_edge_positions,
        ))
    } else {
        None
    };

    let stl_positions = setup.base_mesh.positions;
    let face_indices_arr = setup.face_indices;
    let face_vc_arr = setup.face_vertex_counts;

    // Build compare mesh (translucent)
    let compare_mesh_object = setup.compare_mesh.as_ref().map(|cmp| {
        let cpu_mesh = CpuMesh {
            positions: Positions::F32(cmp.positions.clone()),
            indices: Indices::U32(cmp.indices.clone()),
            normals: Some(cmp.normals.clone()),
            ..Default::default()
        };
        let (cr, cg, cb) = PALETTE[1];
        let mut material = PhysicalMaterial::new_transparent(
            &context,
            &CpuMaterial {
                albedo: Srgba::new(
                    (cr * 255.0) as u8,
                    (cg * 255.0) as u8,
                    (cb * 255.0) as u8,
                    (0.35 * 255.0) as u8,
                ),
                roughness: 0.6,
                metallic: 0.1,
                ..Default::default()
            },
        );
        material.render_states.cull = Cull::Back;
        Gm::new(Mesh::new(&context, &cpu_mesh), material)
    });
    let compare_wireframe = setup.compare_mesh.as_ref().and_then(|cmp| {
        if !cmp.edge_positions.is_empty() {
            Some(WireframeRenderer::new(&context, &cmp.edge_positions))
        } else {
            None
        }
    });

    let ambient = AmbientLight::new(&context, 0.3, Srgba::WHITE);
    let dir1 = DirectionalLight::new(&context, 2.0, Srgba::WHITE, &vec3(-1.0, -0.5, -1.0));
    let dir2 = DirectionalLight::new(&context, 1.0, Srgba::WHITE, &vec3(1.0, 0.3, 0.5));

    let mut show_edges = true;
    let mut show_solid = true;
    let mut hide_edges = true;
    let mut show_soft_edges = true;
    let mut perspective = true;
    let mut auto_center = true;
    let up = vec3(0.0f32, 0.0, 1.0);
    // Maximum |Z-component| of the camera-from-target unit vector to keep
    // the camera at least 15° from the +Z and -Z poles.  cos(15°) ≈ 0.9659.
    const COS_15_DEG: f32 = 0.9659258_f32;

    // Current overlay GPU objects
    let mut highlight_meshes: Vec<Gm<Mesh, PhysicalMaterial>> = Vec::new();
    let mut edge_highlight_renderers: Vec<(WireframeRenderer, [f32; 4])> = Vec::new();
    let mut cylinder_meshes: Vec<Gm<Mesh, PhysicalMaterial>> = Vec::new();
    let mut sphere_meshes: Vec<Gm<Mesh, PhysicalMaterial>> = Vec::new();
    let mut line_renderers: Vec<(WireframeRenderer, [f32; 4], bool)> = Vec::new(); // (renderer, color, no_depth_test)
    let mut shape_mesh_objects: Vec<Gm<Mesh, PhysicalMaterial>> = Vec::new();
    let mut shape_mesh_wireframes: Vec<(WireframeRenderer, [f32; 4])> = Vec::new();
    let mut waiting_for_input = false;

    // Camera animation state (rotation-only slew: target and distance are fixed)
    let mut anim_start_dir: Option<Vector3<f32>> = None; // unit vec from target to camera at anim start
    let mut anim_end_dir: Option<Vector3<f32>> = None;   // unit vec from target to camera at anim end
    let mut anim_target: Option<Vector3<f32>> = None;    // fixed target during slew
    let mut anim_dist: f32 = 0.0;                        // fixed distance during slew
    let mut anim_t: f32 = 1.0; // 1.0 = no animation

    window.render_loop(move |frame_input| {
        // Clear stale GL errors
        unsafe {
            use context::HasContext;
            while context.get_error() != context::NO_ERROR {}
        }

        camera.set_viewport(frame_input.viewport);

        // Camera animation: smoothly rotate to face normal direction (300 ms, rotation-only)
        if anim_t < 1.0 {
            anim_t = (anim_t + frame_input.elapsed_time as f32 / 300.0).min(1.0);
            let t = smooth_step(anim_t);
            if let (Some(start_dir), Some(end_dir), Some(tgt)) =
                (anim_start_dir, anim_end_dir, anim_target)
            {
                // Slerp between start and end view directions (camera-from-target)
                let cos_angle = start_dir.dot(end_dir).clamp(-1.0, 1.0);
                let new_dir = if cos_angle > 0.9999 {
                    end_dir
                } else {
                    let angle = cos_angle.acos();
                    (start_dir * ((1.0 - t) * angle).sin() + end_dir * (t * angle).sin())
                        / angle.sin()
                };
                // Clamp to 15° from poles: |dot(new_dir, Z)| <= cos(15°)
                let z_component = new_dir.z.clamp(-COS_15_DEG, COS_15_DEG);
                let xy_len = (new_dir.x * new_dir.x + new_dir.y * new_dir.y).sqrt();
                let new_dir = if xy_len > 1e-6 {
                    let scale = (1.0 - z_component * z_component).sqrt() / xy_len;
                    vec3(new_dir.x * scale, new_dir.y * scale, z_component).normalize()
                } else {
                    vec3(1.0, 0.0, 0.0) // degenerate: look from +X
                };
                camera.set_view(tgt + new_dir * anim_dist, tgt, up);
            }
        }
        // Check for new overlay from pipeline thread (non-blocking)
        if !waiting_for_input {
            match receiver.overlay_rx.try_recv() {
                Ok(overlay) => {
                    // Build GPU objects for new overlay
                    highlight_meshes = overlay
                        .face_highlights
                        .iter()
                        .map(|h| {
                            build_face_highlight_mesh(
                                &context,
                                &h.face_indices,
                                h.color,
                                &stl_positions,
                                &face_indices_arr,
                                &face_vc_arr,
                            )
                        })
                        .collect();
                    edge_highlight_renderers = overlay
                        .edge_highlights
                        .iter()
                        .map(|e| {
                            let wf = build_face_edge_wireframe(
                                &context,
                                &e.face_indices,
                                &stl_positions,
                                &face_indices_arr,
                                &face_vc_arr,
                            );
                            (wf, e.color)
                        })
                        .collect();
                    cylinder_meshes = overlay
                        .cylinders
                        .iter()
                        .map(|c| build_cylinder_overlay_mesh(&context, c))
                        .collect();
                    sphere_meshes = overlay
                        .spheres
                        .iter()
                        .map(|s| build_sphere_overlay_mesh(&context, s))
                        .collect();
                    line_renderers = overlay
                        .lines
                        .iter()
                        .map(|l| {
                            let positions: Vec<Vector3<f32>> = l
                                .positions
                                .iter()
                                .map(|p| vec3(p[0], p[1], p[2]))
                                .collect();
                            let wf = WireframeRenderer::new(&context, &positions);
                            (wf, l.color, l.no_depth_test)
                        })
                        .collect();
                    shape_mesh_objects = overlay
                        .shape_meshes
                        .iter()
                        .map(|sm| build_shape_mesh_overlay(&context, sm))
                        .collect();
                    shape_mesh_wireframes = overlay
                        .shape_meshes
                        .iter()
                        .map(|sm| {
                            let positions: Vec<Vector3<f32>> = sm
                                .edge_positions
                                .iter()
                                .map(|p| vec3(p[0], p[1], p[2]))
                                .collect();
                            let wf = WireframeRenderer::new(&context, &positions);
                            (wf, sm.edge_color)
                        })
                        .collect();
                    waiting_for_input = true;

                    // Start camera animation if focus_normal is set and auto-center is on.
                    // Slew rotates only: target and distance stay fixed.
                    if auto_center {
                        if let Some(fn_) = overlay.focus_normal {
                            let n = vec3(fn_[0], fn_[1], fn_[2]);
                            let n_len = n.magnitude();
                            if n_len > 1e-6 {
                                let tgt = *camera.target();
                                let dist = camera.position().distance(tgt);
                                let cur_dir = (*camera.position() - tgt).normalize();
                                // End direction: face normal (look from above the face)
                                let mut end = (n / n_len).normalize();
                                // Clamp end direction to 15° from poles
                                let z_component = end.z.clamp(-COS_15_DEG, COS_15_DEG);
                                let xy_len = (end.x * end.x + end.y * end.y).sqrt();
                                end = if xy_len > 1e-6 {
                                    let scale = (1.0 - z_component * z_component).sqrt() / xy_len;
                                    vec3(end.x * scale, end.y * scale, z_component).normalize()
                                } else {
                                    cur_dir // can't slew to a pole, keep current
                                };
                                // Only slew if the angle exceeds 45°
                                // (cos(45°) ≈ 0.7071; dot > threshold means angle < 45°)
                                if cur_dir.dot(end) <= std::f32::consts::FRAC_1_SQRT_2 {
                                    anim_start_dir = Some(cur_dir);
                                    anim_end_dir = Some(end);
                                    anim_target = Some(tgt);
                                    anim_dist = dist;
                                    anim_t = 0.0;
                                }
                            }
                        }
                    }
                }
                Err(mpsc::TryRecvError::Empty) => {}
                Err(mpsc::TryRecvError::Disconnected) => {
                    // Pipeline thread finished
                    return FrameOutput {
                        exit: true,
                        ..Default::default()
                    };
                }
            }
        }

        // Handle events
        let mut should_exit = false;
        for event in frame_input.events.iter() {
            match event {
                Event::KeyPress {
                    kind, modifiers, ..
                } => match kind {
                    Key::Q | Key::Escape => should_exit = true,
                    Key::Space => {
                        if waiting_for_input {
                            let action = if modifiers.shift {
                                VizAction::NextSeed
                            } else {
                                VizAction::NextStep
                            };
                            let _ = receiver.action_tx.send(action);
                            waiting_for_input = false;
                        }
                    }
                    Key::P => {
                        perspective = !perspective;
                        if perspective {
                            camera.set_perspective_projection(degrees(45.0), z_near, z_far);
                        } else {
                            let dist = camera.position().distance(*camera.target());
                            ortho_height =
                                2.0 * dist * (std::f32::consts::FRAC_PI_4 / 2.0).tan();
                            camera.set_orthographic_projection(ortho_height, z_near, z_far);
                        }
                    }
                    Key::E => {
                        if modifiers.shift {
                            show_soft_edges = !show_soft_edges;
                        } else {
                            show_edges = !show_edges;
                        }
                    }
                    Key::S => show_solid = !show_solid,
                    Key::H => hide_edges = !hide_edges,
                    Key::C => auto_center = !auto_center,
                    _ => {}
                },
                Event::MouseMotion {
                    button, delta, ..
                } => match button {
                    Some(MouseButton::Left) => {
                        let angle_h = -delta.0 * 0.005;
                        let angle_v = -delta.1 * 0.005;
                        let rotate =
                            |v: Vector3<f32>, axis: Vector3<f32>, angle: f32| -> Vector3<f32> {
                                let c = angle.cos();
                                let s = angle.sin();
                                v * c + axis.cross(v) * s + axis * axis.dot(v) * (1.0 - c)
                            };
                        // Orbit: rotate camera position around target, keeping Z up.
                        let tgt = *camera.target();
                        let rel_pos = *camera.position() - tgt;
                        // Horizontal rotation around Z axis
                        let rel_pos = rotate(rel_pos, up, angle_h);
                        // Vertical rotation around local right axis
                        let view_dir = (-rel_pos).normalize();
                        let right = view_dir.cross(up);
                        let right_len = right.magnitude();
                        let new_rel_pos = if right_len > 1e-6 {
                            let right = right / right_len;
                            let new_rel_pos = rotate(rel_pos, right, angle_v);
                            // Clamp to 15° from poles (|Z component| <= cos(15°))
                            if new_rel_pos.normalize().z.abs() <= COS_15_DEG {
                                new_rel_pos
                            } else {
                                rel_pos
                            }
                        } else {
                            rel_pos
                        };
                        camera.set_view(tgt + new_rel_pos, tgt, up);
                    }
                    Some(MouseButton::Right) => {
                        let pos = *camera.position();
                        let tgt = *camera.target();
                        let dist = pos.distance(tgt);
                        let view_dir = camera.view_direction();
                        let right = view_dir.cross(up);
                        let right_len = right.magnitude();
                        if right_len > 1e-6 {
                            let right = right / right_len;
                            let cam_up = right.cross(view_dir).normalize();
                            let speed = dist * 0.002;
                            let shift = -right * delta.0 * speed + cam_up * delta.1 * speed;
                            camera.set_view(pos + shift, tgt + shift, up);
                        }
                    }
                    _ => {}
                },
                Event::MouseWheel {
                    delta, position, ..
                } => {
                    let ray_dir = camera.view_direction_at_pixel(*position);
                    let view_dir = camera.view_direction();
                    let ray_origin = camera.position_at_pixel(*position);
                    let denom = ray_dir.dot(view_dir);
                    if denom.abs() > 1e-6 {
                        let t = (center - ray_origin).dot(view_dir) / denom;
                        let cursor_world = ray_origin + ray_dir * t;
                        let factor = (-delta.1 * 0.002).exp();
                        let new_pos =
                            cursor_world + (*camera.position() - cursor_world) * factor;
                        let new_target =
                            cursor_world + (*camera.target() - cursor_world) * factor;
                        let new_dist = new_pos.distance(center);
                        if new_dist >= min_distance && new_dist <= max_distance {
                            camera.set_view(new_pos, new_target, up);
                            if !perspective {
                                ortho_height *= factor;
                                camera.set_orthographic_projection(ortho_height, z_near, z_far);
                            }
                        }
                    }
                }
                _ => {}
            }
        }

        if should_exit {
            let _ = receiver.action_tx.send(VizAction::Quit);
            return FrameOutput {
                exit: true,
                ..Default::default()
            };
        }

        // Render
        let screen = frame_input.screen();
        if show_solid {
            screen.clear(ClearState::color_and_depth(0.15, 0.15, 0.18, 1.0, 1.0));
        } else {
            screen.clear(ClearState::color_and_depth(1.0, 1.0, 1.0, 1.0, 1.0));
        }

        let lights: [&dyn Light; 3] = [&ambient, &dir1, &dir2];

        if show_solid {
            unsafe {
                use context::HasContext;
                context.enable(context::POLYGON_OFFSET_FILL);
                context.polygon_offset(1.0, 1.0);
            }

            // Base mesh
            screen.render(
                &camera,
                std::iter::once(&base_mesh_object as &dyn Object),
                &lights,
            );

            // Highlight overlays (closer to camera in depth)
            unsafe {
                use context::HasContext;
                context.polygon_offset(0.5, 0.5);
            }
            for hm in &highlight_meshes {
                screen.render(&camera, std::iter::once(hm as &dyn Object), &lights);
            }

            unsafe {
                use context::HasContext;
                context.disable(context::POLYGON_OFFSET_FILL);
            }

            // Compare mesh (translucent)
            if let Some(cmp) = &compare_mesh_object {
                unsafe {
                    use context::HasContext;
                    context.enable(context::BLEND);
                    context.blend_func(context::SRC_ALPHA, context::ONE_MINUS_SRC_ALPHA);
                }
                screen.render(&camera, std::iter::once(cmp as &dyn Object), &lights);
                unsafe {
                    use context::HasContext;
                    context.disable(context::BLEND);
                }
            }

            // Translucent cylinder/sphere/shape mesh overlays
            if !cylinder_meshes.is_empty() || !sphere_meshes.is_empty() || !shape_mesh_objects.is_empty() {
                unsafe {
                    use context::HasContext;
                    context.enable(context::BLEND);
                    context.blend_func(context::SRC_ALPHA, context::ONE_MINUS_SRC_ALPHA);
                }
                for cm in &cylinder_meshes {
                    screen.render(&camera, std::iter::once(cm as &dyn Object), &lights);
                }
                for sm in &sphere_meshes {
                    screen.render(&camera, std::iter::once(sm as &dyn Object), &lights);
                }
                for smo in &shape_mesh_objects {
                    screen.render(&camera, std::iter::once(smo as &dyn Object), &lights);
                }
                unsafe {
                    use context::HasContext;
                    context.disable(context::BLEND);
                }
            }
        }

        // Wireframe edges
        if show_edges {
            if let Some(wf) = &base_wireframe {
                wf.render(&context, &camera, [0.05, 0.05, 0.05, 1.0], hide_edges);
            }
            if let Some(wf) = &compare_wireframe {
                wf.render(&context, &camera, [0.3, 0.1, 0.1, 0.5], hide_edges);
            }
        }
        if show_edges && show_soft_edges {
            unsafe {
                use context::HasContext;
                context.enable(context::BLEND);
                context.blend_func(context::SRC_ALPHA, context::ONE_MINUS_SRC_ALPHA);
            }
            if let Some(wf) = &base_soft_wireframe {
                wf.render(&context, &camera, [0.0, 0.0, 0.0, 0.5], hide_edges);
            }
            unsafe {
                use context::HasContext;
                context.disable(context::BLEND);
            }
        }

        // Edge highlight overlays (thick colored outlines on specific faces)
        for (wf, color) in &edge_highlight_renderers {
            wf.render(&context, &camera, *color, true);
        }

        // Shape mesh edge wireframes
        for (wf, color) in &shape_mesh_wireframes {
            wf.render(&context, &camera, *color, true);
        }

        // Line overlays (3D curves, always visible if no_depth_test)
        for (wf, color, no_depth) in &line_renderers {
            wf.render(&context, &camera, *color, !no_depth);
        }

        FrameOutput::default()
    });

    // Wait for pipeline thread to finish
    let _ = pipeline_handle.join();
}

// ---------------------------------------------------------------------------
// MeshData builders for STL and STEP files
// ---------------------------------------------------------------------------

/// Load an STL file into MeshData for visualization.
pub fn load_stl_meshdata(path: &std::path::Path) -> MeshData {
    use opencascade_sys::{message, rw_stl};

    let path_str = path.to_str().expect("invalid path");
    if !path.exists() {
        panic!("STL file not found: {}", path.display());
    }
    let progress = message::ProgressRange::new();
    let tri_handle = rw_stl::read_file_charptr_progressrange_2(path_str, &progress);
    let tri = tri_handle.get();

    let nb_nodes = tri.nb_nodes();
    let nb_tris = tri.nb_triangles();

    let mut positions = Vec::with_capacity(nb_nodes as usize);
    for i in 1..=nb_nodes {
        let pt = tri.node(i);
        positions.push(vec3(pt.x() as f32, pt.y() as f32, pt.z() as f32));
    }

    let mut indices = Vec::with_capacity(nb_tris as usize * 3);
    let mut normals = vec![vec3(0.0f32, 0.0, 0.0); nb_nodes as usize];

    for i in 1..=nb_tris {
        let triangle = tri.triangle(i);
        let mut n1 = 0i32;
        let mut n2 = 0i32;
        let mut n3 = 0i32;
        triangle.get(&mut n1, &mut n2, &mut n3);
        let i0 = (n1 - 1) as u32;
        let i1 = (n2 - 1) as u32;
        let i2 = (n3 - 1) as u32;
        indices.push(i0);
        indices.push(i1);
        indices.push(i2);
        let v0 = positions[i0 as usize];
        let v1 = positions[i1 as usize];
        let v2 = positions[i2 as usize];
        let face_normal = (v1 - v0).cross(v2 - v0);
        normals[i0 as usize] += face_normal;
        normals[i1 as usize] += face_normal;
        normals[i2 as usize] += face_normal;
    }

    for n in &mut normals {
        let len = n.magnitude();
        if len > 1e-10 {
            *n /= len;
        }
    }

    let (edge_positions, soft_edge_positions) = extract_feature_edges(&positions, &indices);
    MeshData {
        positions,
        normals,
        indices,
        edge_positions,
        soft_edge_positions,
    }
}

/// Load an STL file into MeshData with quad fusing.
///
/// Uses brepper's stage1 pipeline to weld vertices, build topology, and fuse
/// coplanar triangle pairs into quads. Edge extraction uses quad-aware logic
/// so that internal quad diagonals are not drawn as edges.
pub fn load_stl_meshdata_quad(path: &std::path::Path) -> MeshData {
    use crate::stage1::{self, VertexWeldOptions};

    let path_str = path.to_str().expect("invalid path");
    let mut mesh =
        stage1::read_connected_mesh_from_stl(path_str, VertexWeldOptions { tolerance: 1.0e-9 })
            .unwrap_or_else(|e| {
                panic!("Failed to read STL {}: {e:?}", path.display());
            });
    mesh.validate_and_populate_topology().unwrap_or_else(|e| {
        panic!("Failed to validate mesh {}: {e:?}", path.display());
    });
    stage1::fuse_coplanar_triangles(&mut mesh, 1.0e-5);

    // Build positions
    let positions: Vec<Vector3<f32>> = mesh
        .vertices
        .iter()
        .map(|v| vec3(v.x as f32, v.y as f32, v.z as f32))
        .collect();

    // Build triangle indices (quads split into two triangles for rendering)
    // and per-vertex normals
    let mut indices: Vec<u32> = Vec::new();
    let mut normals = vec![vec3(0.0f32, 0.0, 0.0); positions.len()];

    // Also collect face data for quad-aware edge extraction
    let mut face_indices_arr: Vec<[usize; 4]> = Vec::new();
    let mut face_vertex_counts: Vec<u8> = Vec::new();

    for face in &mesh.faces {
        if face.vertex_count == 0 {
            continue;
        }
        let vis = &face.vertex_indices;
        let vc = face.vertex_count as usize;

        // Compute face normal
        let v0 = positions[vis[0]];
        let v1 = positions[vis[1]];
        let v2 = positions[vis[2]];
        let face_normal = (v1 - v0).cross(v2 - v0);
        let len = face_normal.magnitude();
        let face_normal = if len > 1e-10 { face_normal / len } else { face_normal };

        // Emit triangles
        indices.push(vis[0] as u32);
        indices.push(vis[1] as u32);
        indices.push(vis[2] as u32);
        if vc == 4 {
            indices.push(vis[0] as u32);
            indices.push(vis[2] as u32);
            indices.push(vis[3] as u32);
        }

        // Accumulate normals
        for vi_idx in 0..vc {
            normals[vis[vi_idx]] += face_normal;
        }

        face_indices_arr.push(face.vertex_indices);
        face_vertex_counts.push(face.vertex_count);
    }

    for n in &mut normals {
        let len = n.magnitude();
        if len > 1e-10 {
            *n /= len;
        }
    }

    let (edge_positions, soft_edge_positions) =
        extract_feature_edges_from_faces(&positions, &face_indices_arr, &face_vertex_counts);
    MeshData {
        positions,
        normals,
        indices,
        edge_positions,
        soft_edge_positions,
    }
}

/// Load a STEP file into MeshData for visualization.
pub fn load_step_meshdata(
    path: &std::path::Path,
    lin_deflection: f64,
    ang_deflection: f64,
) -> MeshData {
    use opencascade_sys::{
        b_rep, b_rep_mesh, message, poly, step_control, top_abs, top_exp, top_loc, topo_ds,
    };

    let path_str = path.to_str().expect("invalid path");
    if !path.exists() {
        panic!("STEP file not found: {}", path.display());
    }
    let mut reader = step_control::StepReader::new();
    let read_status = reader.read_file_charptr(path_str);
    if read_status != opencascade_sys::if_select::ReturnStatus::Retdone {
        panic!("Failed to read STEP file {}: {read_status:?}", path.display());
    }
    let progress = message::ProgressRange::new();
    reader.transfer_roots(&progress);
    let shape = reader.one_shape();
    drop(reader); // release STEP mutex before meshing

    let _mesh = b_rep_mesh::IncrementalMesh::new_shape_real_bool_real_bool(
        &shape,
        lin_deflection,
        false,
        ang_deflection,
        false,
    );

    let mut all_positions: Vec<Vector3<f32>> = Vec::new();
    let mut all_normals: Vec<Vector3<f32>> = Vec::new();
    let mut all_indices: Vec<u32> = Vec::new();
    let mut all_edge_positions: Vec<Vector3<f32>> = Vec::new();
    let mut all_soft_edge_positions: Vec<Vector3<f32>> = Vec::new();

    let mut explorer = top_exp::Explorer::new_shape_shapeenum2(
        &shape,
        top_abs::ShapeEnum::Face,
        top_abs::ShapeEnum::Shape,
    );

    while explorer.more() {
        let face = topo_ds::face_shape(explorer.value());
        let reversed = face.orientation() == top_abs::Orientation::Reversed;
        let mut location = top_loc::Location::new();
        let tri_handle = b_rep::Tool::triangulation(face, &mut location, 0);
        if tri_handle.is_null() {
            explorer.next();
            continue;
        }
        let tri = tri_handle.get();
        if tri.nb_nodes() == 0 {
            explorer.next();
            continue;
        }

        let has_location = !location.is_identity();
        poly::compute_normals(tri_handle);

        let base_idx = all_positions.len() as u32;
        let nb_nodes = tri.nb_nodes();
        let nb_tris = tri.nb_triangles();

        for i in 1..=nb_nodes {
            let pt = tri.node(i);
            let (mut x, mut y, mut z) = (pt.x(), pt.y(), pt.z());
            if has_location {
                location
                    .transformation()
                    .transforms_real3(&mut x, &mut y, &mut z);
            }
            all_positions.push(vec3(x as f32, y as f32, z as f32));
        }

        let normal_sign = if reversed { -1.0f32 } else { 1.0f32 };
        for i in 1..=nb_nodes {
            let normal = tri.normal_int(i);
            all_normals.push(vec3(
                normal.x() as f32 * normal_sign,
                normal.y() as f32 * normal_sign,
                normal.z() as f32 * normal_sign,
            ));
        }

        let mut face_positions = Vec::with_capacity(nb_nodes as usize);
        let mut face_indices = Vec::with_capacity(nb_tris as usize * 3);

        for i in 1..=nb_tris {
            let triangle = tri.triangle(i);
            let mut n1 = 0i32;
            let mut n2 = 0i32;
            let mut n3 = 0i32;
            triangle.get(&mut n1, &mut n2, &mut n3);
            if reversed {
                all_indices.push(base_idx + (n1 - 1) as u32);
                all_indices.push(base_idx + (n3 - 1) as u32);
                all_indices.push(base_idx + (n2 - 1) as u32);
            } else {
                all_indices.push(base_idx + (n1 - 1) as u32);
                all_indices.push(base_idx + (n2 - 1) as u32);
                all_indices.push(base_idx + (n3 - 1) as u32);
            }
            face_indices.push((n1 - 1) as u32);
            face_indices.push((n2 - 1) as u32);
            face_indices.push((n3 - 1) as u32);
        }

        for i in 1..=nb_nodes {
            let pt = tri.node(i);
            let (mut x, mut y, mut z) = (pt.x(), pt.y(), pt.z());
            if has_location {
                location
                    .transformation()
                    .transforms_real3(&mut x, &mut y, &mut z);
            }
            face_positions.push(vec3(x as f32, y as f32, z as f32));
        }

        let (edges, soft_edges) = extract_feature_edges(&face_positions, &face_indices);
        all_edge_positions.extend(edges);
        all_soft_edge_positions.extend(soft_edges);
        explorer.next();
    }

    MeshData {
        positions: all_positions,
        normals: all_normals,
        indices: all_indices,
        edge_positions: all_edge_positions,
        soft_edge_positions: all_soft_edge_positions,
    }
}

// ---------------------------------------------------------------------------
// Build MeshData from ConnectedMesh
// ---------------------------------------------------------------------------

/// Build a MeshData from a ConnectedMesh (after stage 1) for visualization.
/// Also returns face index/vertex-count arrays for building overlays.
pub fn meshdata_from_connected_mesh(
    mesh: &crate::stage1::ConnectedMesh,
) -> (MeshData, Vec<[usize; 4]>, Vec<u8>) {
    let mut positions: Vec<Vector3<f32>> = Vec::with_capacity(mesh.vertices.len());
    for v in &mesh.vertices {
        positions.push(vec3(v.x as f32, v.y as f32, v.z as f32));
    }

    let mut indices = Vec::new();
    let mut normals_accum = vec![vec3(0.0f32, 0.0, 0.0); mesh.vertices.len()];
    let mut face_indices_arr = Vec::with_capacity(mesh.faces.len());
    let mut face_vc_arr = Vec::with_capacity(mesh.faces.len());

    for face in &mesh.faces {
        let vc = face.vertex_count as usize;
        face_indices_arr.push(face.vertex_indices);
        face_vc_arr.push(face.vertex_count);

        if vc >= 3 {
            let i0 = face.vertex_indices[0] as u32;
            let i1 = face.vertex_indices[1] as u32;
            let i2 = face.vertex_indices[2] as u32;
            indices.push(i0);
            indices.push(i1);
            indices.push(i2);
            let v0 = positions[i0 as usize];
            let v1 = positions[i1 as usize];
            let v2 = positions[i2 as usize];
            let n = (v1 - v0).cross(v2 - v0);
            normals_accum[i0 as usize] += n;
            normals_accum[i1 as usize] += n;
            normals_accum[i2 as usize] += n;
        }
        if vc >= 4 {
            let i0 = face.vertex_indices[0] as u32;
            let i2 = face.vertex_indices[2] as u32;
            let i3 = face.vertex_indices[3] as u32;
            indices.push(i0);
            indices.push(i2);
            indices.push(i3);
            let v0 = positions[i0 as usize];
            let v2 = positions[i2 as usize];
            let v3 = positions[i3 as usize];
            let n = (v2 - v0).cross(v3 - v0);
            normals_accum[i0 as usize] += n;
            normals_accum[i2 as usize] += n;
            normals_accum[i3 as usize] += n;
        }
    }

    for n in &mut normals_accum {
        let len = n.magnitude();
        if len > 1e-10 {
            *n /= len;
        }
    }

    let (edge_positions, soft_edge_positions) =
        extract_feature_edges_from_faces(&positions, &face_indices_arr, &face_vc_arr);

    let meshdata = MeshData {
        positions,
        normals: normals_accum,
        indices,
        edge_positions,
        soft_edge_positions,
    };

    (meshdata, face_indices_arr, face_vc_arr)
}

// ---------------------------------------------------------------------------
// Tessellate OCCT shapes for visualization
// ---------------------------------------------------------------------------

/// Tessellate an OCCT TopoDS_Shape into a ShapeMeshOverlay ready for viz.
/// Uses BRepMesh_IncrementalMesh for tessellation.
pub fn tessellate_shape(
    shape: &opencascade_sys::topo_ds::Shape,
    color: [f32; 4],
    edge_color: [f32; 4],
    lin_deflection: f64,
    ang_deflection: f64,
) -> ShapeMeshOverlay {
    use opencascade_sys::{
        b_rep, b_rep_mesh, poly, top_abs, top_exp, top_loc, topo_ds,
    };

    let _mesh = b_rep_mesh::IncrementalMesh::new_shape_real_bool_real_bool(
        shape,
        lin_deflection,
        false,
        ang_deflection,
        false,
    );

    let mut all_positions: Vec<[f32; 3]> = Vec::new();
    let mut all_normals: Vec<[f32; 3]> = Vec::new();
    let mut all_indices: Vec<u32> = Vec::new();
    let mut all_edge_positions: Vec<[f32; 3]> = Vec::new();

    let mut explorer = top_exp::Explorer::new_shape_shapeenum2(
        shape,
        top_abs::ShapeEnum::Face,
        top_abs::ShapeEnum::Shape,
    );

    while explorer.more() {
        let face = topo_ds::face_shape(explorer.value());
        let reversed = face.orientation() == top_abs::Orientation::Reversed;
        let mut location = top_loc::Location::new();
        let tri_handle = b_rep::Tool::triangulation(face, &mut location, 0);
        if tri_handle.is_null() {
            explorer.next();
            continue;
        }
        let tri = tri_handle.get();
        if tri.nb_nodes() == 0 {
            explorer.next();
            continue;
        }

        let has_location = !location.is_identity();
        poly::compute_normals(tri_handle);

        let base_idx = all_positions.len() as u32;
        let nb_nodes = tri.nb_nodes();
        let nb_tris = tri.nb_triangles();

        for i in 1..=nb_nodes {
            let pt = tri.node(i);
            let (mut x, mut y, mut z) = (pt.x(), pt.y(), pt.z());
            if has_location {
                location
                    .transformation()
                    .transforms_real3(&mut x, &mut y, &mut z);
            }
            all_positions.push([x as f32, y as f32, z as f32]);
        }

        let normal_sign = if reversed { -1.0f32 } else { 1.0f32 };
        for i in 1..=nb_nodes {
            let normal = tri.normal_int(i);
            all_normals.push([
                normal.x() as f32 * normal_sign,
                normal.y() as f32 * normal_sign,
                normal.z() as f32 * normal_sign,
            ]);
        }

        for i in 1..=nb_tris {
            let triangle = tri.triangle(i);
            let mut n1 = 0i32;
            let mut n2 = 0i32;
            let mut n3 = 0i32;
            triangle.get(&mut n1, &mut n2, &mut n3);
            if reversed {
                all_indices.push(base_idx + (n1 - 1) as u32);
                all_indices.push(base_idx + (n3 - 1) as u32);
                all_indices.push(base_idx + (n2 - 1) as u32);
            } else {
                all_indices.push(base_idx + (n1 - 1) as u32);
                all_indices.push(base_idx + (n2 - 1) as u32);
                all_indices.push(base_idx + (n3 - 1) as u32);
            }
        }

        explorer.next();
    }

    // Extract feature edges for the tessellation
    let positions_v3: Vec<Vector3<f32>> = all_positions.iter().map(|p| vec3(p[0], p[1], p[2])).collect();
    let (sharp_edges, _soft_edges) = extract_feature_edges(&positions_v3, &all_indices);
    for p in &sharp_edges {
        all_edge_positions.push([p.x, p.y, p.z]);
    }

    ShapeMeshOverlay {
        positions: all_positions,
        normals: all_normals,
        indices: all_indices,
        color,
        edge_positions: all_edge_positions,
        edge_color,
    }
}

/// Sample points along an OCCT curve for line rendering.
/// Returns pairs of points [p0, p1, p1, p2, ...] forming a polyline.
pub fn sample_curve_for_viz(
    curve: &opencascade_sys::geom::HandleGeomCurve,
    num_segments: usize,
) -> Vec<[f32; 3]> {
    let c = curve.get();
    let fp = c.first_parameter();
    let lp = c.last_parameter();
    let mut positions = Vec::with_capacity(num_segments * 2);
    for i in 0..num_segments {
        let t0 = fp + (lp - fp) * (i as f64) / (num_segments as f64);
        let t1 = fp + (lp - fp) * ((i + 1) as f64) / (num_segments as f64);
        let p0 = c.value(t0);
        let p1 = c.value(t1);
        positions.push([p0.x() as f32, p0.y() as f32, p0.z() as f32]);
        positions.push([p1.x() as f32, p1.y() as f32, p1.z() as f32]);
    }
    positions
}
