// Interactive 3D viewer for STL and STEP files.
// Uses opencascade-sys for geometry loading and three-d for rendering.
// Wireframe edges rendered via raw glow GL_LINES calls.

use clap::Parser;
use opencascade_sys::{
    b_rep, b_rep_mesh, message, poly, rw_stl, step_control, top_abs, top_exp, top_loc, topo_ds,
};
use std::collections::HashMap;
use std::path::PathBuf;
use three_d::*;

#[derive(Parser)]
#[command(
    name = "viewer",
    about = "Interactive 3D viewer for STL and STEP files"
)]
struct Cli {
    /// Files to view (STL or STEP, detected by extension)
    #[arg(required = true)]
    files: Vec<PathBuf>,

    /// Linear deflection for STEP tessellation (smaller = finer mesh)
    #[arg(long, default_value = "0.1")]
    deflection: f64,

    /// Angular deflection for STEP tessellation (radians)
    #[arg(long, default_value = "0.5")]
    angular_deflection: f64,
}

struct MeshData {
    positions: Vec<Vector3<f32>>,
    normals: Vec<Vector3<f32>>,
    indices: Vec<u32>,
    /// Wireframe edge positions as pairs [p0, p1, p0, p1, ...].
    edge_positions: Vec<Vector3<f32>>,
    /// Near-coplanar edges (drawn at reduced alpha).
    soft_edge_positions: Vec<Vector3<f32>>,
}

const PALETTE: &[(f32, f32, f32)] = &[
    (0.55, 0.65, 0.80),
    (0.80, 0.55, 0.55),
    (0.55, 0.78, 0.55),
    (0.78, 0.70, 0.50),
    (0.70, 0.55, 0.78),
    (0.55, 0.75, 0.75),
];

fn load_stl(path: &std::path::Path) -> MeshData {
    let path_str = path.to_str().expect("invalid path");
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

fn load_step(path: &std::path::Path, lin_deflection: f64, ang_deflection: f64) -> MeshData {
    let path_str = path.to_str().expect("invalid path");

    let mut reader = step_control::Reader::new();
    reader.read_file_charptr(path_str);
    let progress = message::ProgressRange::new();
    reader.transfer_roots(&progress);
    let shape = reader.one_shape();

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

        // Flip normals for reversed faces so lighting is correct from outside
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
            // Flip winding order for reversed faces so back-face culling works correctly
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

/// Extract boundary/sharp edges and soft (near-coplanar) edges.
/// Returns (sharp_edges, soft_edges) as position pairs.
fn extract_feature_edges(positions: &[Vector3<f32>], indices: &[u32]) -> (Vec<Vector3<f32>>, Vec<Vector3<f32>>) {
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
            // Boundary edge - always sharp
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

fn compute_bounding_box(meshes: &[&MeshData]) -> (Vector3<f32>, Vector3<f32>) {
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

// --- Raw GL wireframe renderer ---

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

struct WireframeRenderer {
    program: context::Program,
    vbo: context::Buffer,
    vertex_count: i32,
}

impl WireframeRenderer {
    fn new(ctx: &Context, positions: &[Vector3<f32>]) -> Self {
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

    fn render(&self, ctx: &Context, camera: &Camera, color: [f32; 4], depth_test: bool) {
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

            // Create a temporary VAO for line drawing
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

            // Restore depth test state
            ctx.enable(context::DEPTH_TEST);
            ctx.depth_func(context::LESS);

            ctx.bind_vertex_array(None);
            ctx.delete_vertex_array(vao);
            ctx.use_program(None);

            // Unbind our state; three-d will rebind its own VAO on next render
            ctx.bind_buffer(context::ARRAY_BUFFER, None);
        }
    }
}

// --- Main ---

fn main() {
    let cli = Cli::parse();

    let mut meshes: Vec<(MeshData, usize)> = Vec::new();
    for (i, file) in cli.files.iter().enumerate() {
        let file = std::fs::canonicalize(file).unwrap_or_else(|e| {
            eprintln!("Cannot resolve path '{}': {}", file.display(), e);
            std::process::exit(1);
        });
        let ext = file
            .extension()
            .and_then(|e| e.to_str())
            .unwrap_or("")
            .to_lowercase();
        let color_idx = i % PALETTE.len();
        match ext.as_str() {
            "stl" => {
                eprintln!("Loading STL: {}", file.display());
                let data = load_stl(&file);
                eprintln!(
                    "  {} verts, {} tris, {} edge segs, {} soft edge segs",
                    data.positions.len(),
                    data.indices.len() / 3,
                    data.edge_positions.len() / 2,
                    data.soft_edge_positions.len() / 2
                );
                meshes.push((data, color_idx));
            }
            "step" | "stp" => {
                eprintln!("Loading STEP: {}", file.display());
                let data = load_step(&file, cli.deflection, cli.angular_deflection);
                eprintln!(
                    "  {} verts, {} tris, {} edge segs, {} soft edge segs",
                    data.positions.len(),
                    data.indices.len() / 3,
                    data.edge_positions.len() / 2,
                    data.soft_edge_positions.len() / 2
                );
                meshes.push((data, color_idx));
            }
            _ => {
                eprintln!("Unknown extension '{}': {}", ext, file.display());
                std::process::exit(1);
            }
        }
    }

    if meshes.is_empty() {
        eprintln!("No files loaded.");
        std::process::exit(1);
    }

    let refs: Vec<&MeshData> = meshes.iter().map(|(d, _)| d).collect();
    let (bb_min, bb_max) = compute_bounding_box(&refs);
    let center = (bb_min + bb_max) * 0.5;
    let extent = (bb_max - bb_min).magnitude();

    let window = Window::new(WindowSettings {
        title: "brepper viewer".to_string(),
        max_size: Some((1920, 1200)),
        ..Default::default()
    })
    .unwrap();
    let context = window.gl();

    let mut camera = Camera::new_perspective(
        window.viewport(),
        center + vec3(extent * 0.5, extent * 0.3, extent * 0.5),
        center,
        vec3(0.0, 0.0, 1.0), // Z-up
        degrees(45.0),
        extent * 0.001,
        extent * 100.0,
    );
    let min_distance = extent * 0.01;
    let max_distance = extent * 10.0;

    let mut shaded_objects: Vec<Gm<Mesh, PhysicalMaterial>> = Vec::new();
    let mut wireframes: Vec<(WireframeRenderer, [f32; 4])> = Vec::new();
    let mut soft_wireframes: Vec<WireframeRenderer> = Vec::new();

    for (data, color_idx) in &meshes {
        let (cr, cg, cb) = PALETTE[*color_idx];

        let cpu_mesh = CpuMesh {
            positions: Positions::F32(data.positions.clone()),
            indices: Indices::U32(data.indices.clone()),
            normals: Some(data.normals.clone()),
            ..Default::default()
        };

        let mut material = PhysicalMaterial::new_opaque(
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
        material.render_states.cull = Cull::Back;
        shaded_objects.push(Gm::new(Mesh::new(&context, &cpu_mesh), material));

        if !data.edge_positions.is_empty() {
            wireframes.push((
                WireframeRenderer::new(&context, &data.edge_positions),
                [0.05, 0.05, 0.05, 1.0],
            ));
        }
        if !data.soft_edge_positions.is_empty() {
            soft_wireframes.push(WireframeRenderer::new(&context, &data.soft_edge_positions));
        }
    }

    let ambient = AmbientLight::new(&context, 0.3, Srgba::WHITE);
    let dir1 = DirectionalLight::new(&context, 2.0, Srgba::WHITE, &vec3(-1.0, -0.5, -1.0));
    let dir2 = DirectionalLight::new(&context, 1.0, Srgba::WHITE, &vec3(1.0, 0.3, 0.5));

    let mut show_edges = true;
    let mut show_solid = true;
    let mut hide_edges = true;
    let mut show_soft_edges = true;
    let mut perspective = true;
    let mut ortho_height = 2.0 * extent * (std::f32::consts::FRAC_PI_4 / 2.0).tan();
    let z_near = extent * 0.001;
    let z_far = extent * 100.0;
    let up = vec3(0.0f32, 0.0, 1.0);

    window.render_loop(move |frame_input| {
        // Clear any stale GL errors (e.g. from wireframe renderer on previous frame)
        unsafe {
            use context::HasContext;
            while context.get_error() != context::NO_ERROR {}
        }

        camera.set_viewport(frame_input.viewport);

        // Handle events: orbit, pan, zoom, keyboard
        let mut should_exit = false;
        for event in frame_input.events.iter() {
            match event {
                Event::KeyPress { kind, modifiers, .. } => match kind {
                    Key::Q | Key::Escape => should_exit = true,
                    Key::P => {
                        perspective = !perspective;
                        if perspective {
                            camera.set_perspective_projection(degrees(45.0), z_near, z_far);
                        } else {
                            let dist = camera.position().distance(*camera.target());
                            ortho_height = 2.0 * dist * (std::f32::consts::FRAC_PI_4 / 2.0).tan();
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
                    _ => {}
                },
                Event::MouseMotion { button, delta, .. } => match button {
                    Some(MouseButton::Left) => {
                        // Orbit around model center
                        let angle_h = -delta.0 * 0.005;
                        let angle_v = -delta.1 * 0.005;
                        let rotate = |v: Vector3<f32>, axis: Vector3<f32>, angle: f32| -> Vector3<f32> {
                            let c = angle.cos();
                            let s = angle.sin();
                            v * c + axis.cross(v) * s + axis * axis.dot(v) * (1.0 - c)
                        };
                        let rel_pos = *camera.position() - center;
                        let rel_target = *camera.target() - center;
                        // Horizontal rotation around Z
                        let rel_pos = rotate(rel_pos, up, angle_h);
                        let rel_target = rotate(rel_target, up, angle_h);
                        // Recompute right vector for vertical rotation
                        let view_dir = ((center + rel_target) - (center + rel_pos)).normalize();
                        let right = view_dir.cross(up);
                        let right_len = right.magnitude();
                        if right_len > 1e-6 {
                            let right = right / right_len;
                            let new_rel_pos = rotate(rel_pos, right, angle_v);
                            // Prevent flipping past poles
                            if new_rel_pos.normalize().dot(up).abs() < 0.98 {
                                let new_rel_target = rotate(rel_target, right, angle_v);
                                camera.set_view(center + new_rel_pos, center + new_rel_target, up);
                            }
                        }
                    }
                    Some(MouseButton::Right) => {
                        // Pan
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
                Event::MouseWheel { delta, position, .. } => {
                    // Zoom at cursor point
                    let ray_dir = camera.view_direction_at_pixel(*position);
                    let view_dir = camera.view_direction();
                    let ray_origin = camera.position_at_pixel(*position);
                    let denom = ray_dir.dot(view_dir);
                    if denom.abs() > 1e-6 {
                        let t = (center - ray_origin).dot(view_dir) / denom;
                        let cursor_world = ray_origin + ray_dir * t;
                        let factor = (-delta.1 * 0.002).exp();
                        let new_pos = cursor_world + (*camera.position() - cursor_world) * factor;
                        let new_target = cursor_world + (*camera.target() - cursor_world) * factor;
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
            return FrameOutput { exit: true, ..Default::default() };
        }

        let screen = frame_input.screen();
        if show_solid {
            screen.clear(ClearState::color_and_depth(0.15, 0.15, 0.18, 1.0, 1.0));
        } else {
            screen.clear(ClearState::color_and_depth(1.0, 1.0, 1.0, 1.0, 1.0));
        }

        if show_solid {
            // Push solid geometry slightly away from camera so wireframe edges
            // on coplanar faces win the depth test
            unsafe {
                use context::HasContext;
                context.enable(context::POLYGON_OFFSET_FILL);
                context.polygon_offset(1.0, 1.0);
            }
            screen.render(
                &camera,
                shaded_objects.iter().map(|o| o as &dyn Object),
                &[
                    &ambient as &dyn Light,
                    &dir1 as &dyn Light,
                    &dir2 as &dyn Light,
                ],
            );
            unsafe {
                use context::HasContext;
                context.disable(context::POLYGON_OFFSET_FILL);
            }
        }

        if show_edges {
            for (wf, color) in &wireframes {
                wf.render(&context, &camera, *color, hide_edges);
            }
        }

        if show_edges && show_soft_edges {
            unsafe {
                use context::HasContext;
                context.enable(context::BLEND);
                context.blend_func(context::SRC_ALPHA, context::ONE_MINUS_SRC_ALPHA);
            }
            for wf in &soft_wireframes {
                wf.render(&context, &camera, [0.0, 0.0, 0.0, 0.5], hide_edges);
            }
            unsafe {
                use context::HasContext;
                context.disable(context::BLEND);
            }
        }

        FrameOutput::default()
    });
}
