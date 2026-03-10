// Interactive 3D viewer for STL and STEP files.
// Uses opencascade-sys for geometry loading and three-d for rendering.
// Shared rendering infrastructure lives in brepper::viz.

use clap::Parser;
use brepper::viz::{self, MeshData, WireframeRenderer, PALETTE};
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

    /// Angular deflection for STEP tessellation in degrees (default: 28.6 = 0.5 rad)
    #[arg(long, default_value = "28.6")]
    angular_deflection: f64,
}

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
                let data = viz::load_stl_meshdata(&file);
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
                let data = viz::load_step_meshdata(&file, cli.deflection, cli.angular_deflection.to_radians());
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
    let (bb_min, bb_max) = viz::compute_bounding_box(&refs);
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
