use clap::{ArgAction, Parser, ValueEnum};
use opencascade_sys::{top_abs, top_exp};
use std::fmt::{Display, Formatter};

// ---------------------------------------------------------------------------
// Stage
// ---------------------------------------------------------------------------

/// Stage identifier as (major, minor). Minor=0 means run all sub-stages of that major stage.
/// Examples: (1, 2) = stage 1.2, (2, 0) = all of stage 2, (4, 1) = full pipeline.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Stage(pub u8, pub u8);

impl Stage {
    /// The final sub-stage for each major stage.
    fn last_minor(major: u8) -> u8 {
        match major {
            1 => 3,
            2 => 9,
            3 => 6,
            4 => 1,
            _ => 0,
        }
    }

    /// Normalize: if minor is 0, expand to last sub-stage of that major.
    pub fn normalized(self) -> (u8, u8) {
        if self.1 == 0 {
            (self.0, Self::last_minor(self.0))
        } else {
            (self.0, self.1)
        }
    }

    /// Returns true if the pipeline should run at least through the given stage.
    pub fn at_least(&self, major: u8, minor: u8) -> bool {
        let (sm, ss) = self.normalized();
        (sm, ss) >= (major, minor)
    }
}

impl Default for Stage {
    fn default() -> Self {
        Stage(4, 1)
    }
}

fn parse_stage(s: &str) -> Result<Stage, String> {
    if let Some((major_s, minor_s)) = s.split_once('.') {
        let major: u8 = major_s.parse().map_err(|_| format!("invalid stage: {s}"))?;
        let minor: u8 = minor_s.parse().map_err(|_| format!("invalid stage: {s}"))?;
        if !(1..=4).contains(&major) || !(1..=Stage::last_minor(major)).contains(&minor) {
            return Err(format!("invalid stage: {s} (valid: 1.1–4.1)"));
        }
        Ok(Stage(major, minor))
    } else {
        let major: u8 = s.parse().map_err(|_| format!("invalid stage: {s}"))?;
        if !(1..=4).contains(&major) {
            return Err(format!("invalid stage: {s} (valid: 1–4)"));
        }
        Ok(Stage(major, 0))
    }
}

// ---------------------------------------------------------------------------
// Units
// ---------------------------------------------------------------------------

/// Length units for STL/STEP files.
#[derive(Debug, Default, Clone, Copy, PartialEq, ValueEnum)]
pub enum Units {
    /// Millimeters
    #[default]
    Mm,
    /// Centimeters
    Cm,
    /// Meters
    M,
    /// Inches
    In,
    /// Feet
    Ft,
    /// Micrometers
    Um,
}

impl Units {
    /// Conversion factor: 1 of this unit = N millimeters.
    pub fn to_mm(self) -> f64 {
        match self {
            Units::Mm => 1.0,
            Units::Cm => 10.0,
            Units::M => 1000.0,
            Units::In => 25.4,
            Units::Ft => 304.8,
            Units::Um => 0.001,
        }
    }
}


impl Display for Units {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        match self {
            Units::Mm => write!(f, "mm"),
            Units::Cm => write!(f, "cm"),
            Units::M => write!(f, "m"),
            Units::In => write!(f, "in"),
            Units::Ft => write!(f, "ft"),
            Units::Um => write!(f, "μm"),
        }
    }
}

// ---------------------------------------------------------------------------
// CLI args (parsed by clap)
// ---------------------------------------------------------------------------

/// brepper — Convert STL mesh to STEP with fitted analytic and freeform surfaces.
#[derive(Parser, Debug)]
#[command(name = "brepper", version, about)]
struct CliArgs {
    /// Input STL file (binary or ASCII).
    pub input_stl: String,

    /// Output STEP file.
    #[arg(short = 'o', long = "output")]
    pub output_step: Option<String>,

    /// Units used by the STL file.
    #[arg(long, default_value = "mm", value_enum)]
    pub stl_units: Units,

    /// Units to use in exported STEP file.
    #[arg(long, default_value = "mm", value_enum)]
    pub step_units: Units,

    /// STEP file to compare against at each stage.
    #[arg(long = "compare")]
    pub compare_step: Option<String>,

    /// Score the mesh and our segmentation against the --compare reference:
    /// facet-normal exactness, closed-form initializer error, the seam
    /// assumption, and over/under-segmentation. Diagnostics only — see
    /// [`crate::reference`]. Requires --compare.
    #[arg(long = "reference-report")]
    pub reference_report: bool,

    /// Measure the DEVELOPMENT_PLAN "Invariants" against the --compare
    /// reference's ground truth: whether the stated fitting requirement is
    /// satisfiable, whether an angular tolerance can separate adjacent
    /// surfaces, whether a wrong surface type also explains a true support, and
    /// whether growth from different seeds and traversal orders converges on
    /// the same surface. Diagnostics only — see [`crate::stage2::invariants`].
    /// Requires --compare; runs after stage 2.2.
    #[arg(long = "invariant-report")]
    pub invariant_report: bool,

    /// Fitting tolerance for vertex-to-surface distance, in STL units.
    #[arg(long, default_value = "5e-5")]
    pub vertex_tolerance: f64,

    /// Tolerance for surface-to-triangle-face offset, in STL units.
    #[arg(long, default_value = "0.4")]
    pub surface_tolerance: f64,

    /// Max dihedral angle between adjacent triangles on same surface, in
    /// degrees. Medium-resolution CAD exports step small-radius fillets up
    /// to ~24° between rows (Onshape "medium" enforces chordal deviation, so
    /// small features get few segments regardless of the angular preset);
    /// the default must admit those.
    #[arg(long, default_value = "25")]
    pub angular_tolerance: f64,

    /// Increment verbosity by 1 per flag (-v = 1, -vv = 2, -vvv = 3).
    #[arg(short = 'v', action = ArgAction::Count)]
    pub verbose_short: u8,

    /// Set verbosity to an explicit level:
    /// --verbose=1 (stage summaries),
    /// --verbose=2 (per-hypothesis details),
    /// --verbose=3 (per-committed-hypothesis seed/count info, stats, profiles),
    /// --verbose=4 (per-seed-attempt details),
    /// --verbose=5 (per-BFS-face exploration details).
    /// Bare --verbose sets level 1.
    #[arg(long = "verbose", value_name = "LEVEL", default_missing_value = "1", num_args = 0..=1)]
    pub verbose_level: Option<u8>,

    /// Suppress non-error output.
    #[arg(short = 'q', long)]
    pub quiet: bool,

    /// Fail on degenerate mesh faces instead of repairing them.
    #[arg(long)]
    pub strict_mesh: bool,

    /// Worker threads for stage 2's fit competition (default: all cores;
    /// 1 = serial). The result is identical for every value.
    #[arg(long, value_parser = clap::value_parser!(u32).range(1..))]
    pub threads: Option<u32>,

    /// Enable debug output and intermediate files.
    #[arg(long)]
    pub debug: bool,

    /// Stop after this stage (e.g. "2.2" or "2"). Default: run all stages.
    #[arg(long, value_parser = parse_stage)]
    pub stage: Option<Stage>,

    /// Comma-separated list of stages for interactive visualization (e.g. "2.2,3.1").
    #[arg(long = "viz-stages", value_parser = parse_stage, value_delimiter = ',')]
    pub viz_stages: Option<Vec<Stage>>,

    /// Detail level for viz-stages visualization pauses:
    /// 1 = pause at committed hypotheses and stage 2.6 selections (default),
    /// 2 = also pause at each accepted seed,
    /// 3 = also pause at each BFS face acceptance.
    #[arg(long = "viz-stages-level", default_value = "1")]
    pub viz_stages_level: u8,

    /// Interactively visualize how SEED SETS grow under the clean
    /// grow-with-refit operator, against `--compare` ground truth. Each step
    /// shows one seed set and what it grew to: seed green, correctly covered
    /// blue, overshoot red, missed yellow, and the facets that blocked further
    /// growth orange. Requires --compare. See [`crate::stage2::invariants`].
    #[arg(long = "viz-seeds")]
    pub viz_seeds: bool,

    /// Which seed sets to stop on: `failing` (default — anything that does not
    /// cover its surface), `all`, or `covering`.
    #[arg(long = "viz-seed-show", default_value = "failing")]
    pub viz_seed_show: String,

    /// Facet counts the seed-set generators may use. Only the maximum is read,
    /// as the bound on chain length, so "2,3,4" means "chains of up to 4
    /// facets" — the sizes the DEVELOPMENT_PLAN Invariants section specifies.
    #[arg(long = "viz-seed-sizes", default_value = "2,3,4", value_delimiter = ',')]
    pub viz_seed_sizes: Vec<usize>,

    /// How many anchor facets to sample per reference face.
    #[arg(long = "viz-seed-samples", default_value = "12")]
    pub viz_seed_samples: usize,

    /// Max conforming seed sets to try per anchor facet before calling it a
    /// failure.
    #[arg(long = "viz-seed-max-attempts", default_value = "12")]
    pub viz_seed_max_attempts: usize,

    /// Enforce the section's "no more than three coplanar" clause on
    /// cylinder/cone seeds. Off by default: measured, it rejects every seed
    /// containing a fused quad, which is most of them.
    #[arg(long = "viz-seed-coplanar-clause")]
    pub viz_seed_coplanar_clause: bool,

    /// Confine seed growth to the stage-2.1 smooth region instead of letting
    /// validation be the only limit. Off by default, because the sharp-edge
    /// partition demonstrably cuts true surfaces in two (measured: 350 edges
    /// inside 68 reference faces), and a seed on the far side of such a cut can
    /// never cover its surface however good it is — which is a partitioning
    /// failure, not a seeding one. Turn it on to see that effect.
    #[arg(long = "viz-seed-confine-region")]
    pub viz_seed_confine_region: bool,

    /// Try only the seed patterns the DEVELOPMENT_PLAN "Surface Seed
    /// Requirements" name, without the grown-ball fallback. Use this to see
    /// where those rules alone are insufficient.
    #[arg(long = "viz-seed-spec-only")]
    pub viz_seed_spec_only: bool,
}

// ---------------------------------------------------------------------------
// Config (the resolved, normalized configuration used by the pipeline)
// ---------------------------------------------------------------------------

/// Resolved configuration used by the pipeline. Tolerances are stored in mm.
pub struct Config {
    /// Input STL file path.
    pub input_stl: String,
    /// Output STEP file path.
    pub output_step: Option<String>,
    /// Units of the input STL file.
    pub stl_units: Units,
    /// Units for the output STEP file.
    pub step_units: Units,
    /// STEP file path to compare against at each stage.
    pub compare_step: Option<String>,
    /// Shape loaded from the --compare STEP file, for bounded distance checks.
    pub compare_shape: Option<opencascade_sys::OwnedPtr<opencascade_sys::topo_ds::Shape>>,
    /// Number of faces in the compare shape (for status output).
    pub compare_face_count: usize,
    /// Print the reference-scoring diagnostics (see [`crate::reference`]).
    pub reference_report: bool,
    /// Print the invariant measurements (see [`crate::stage2::invariants`]).
    pub invariant_report: bool,
    /// Interactively visualize seed-set growth (see [`crate::stage2::invariants`]).
    pub viz_seeds: bool,
    /// Which seed sets to stop on in `--viz-seeds`.
    pub viz_seed_show: String,
    /// Bound on seed-set facet count in `--viz-seeds` (max is what is read).
    pub viz_seed_sizes: Vec<usize>,
    /// Anchor facets sampled per reference face in `--viz-seeds`.
    pub viz_seed_samples: usize,
    /// Max conforming seed sets tried per anchor facet.
    pub viz_seed_max_attempts: usize,
    /// Enforce the section's coplanarity clause on developable seeds.
    pub viz_seed_coplanar_clause: bool,
    /// Confine seed growth to the stage-2.1 smooth region.
    pub viz_seed_confine_region: bool,
    /// Try only the section's own seed patterns, no grown-ball fallback.
    pub viz_seed_spec_only: bool,
    /// Fitting tolerance for vertex-to-surface distance, in mm.
    pub vertex_tolerance_mm: f64,
    /// Tolerance for surface-to-triangle-face offset, in mm.
    pub surface_tolerance_mm: f64,
    /// Max dihedral angle (radians) between adjacent triangles on the same surface.
    pub angular_tolerance_rad: f64,
    /// Verbosity level (0 = quiet, 1 = summaries, 2 = hypotheses, 3 = commit details, 4 = seed details, 5 = BFS trace).
    pub verbosity: u8,
    /// Convenience: verbosity >= 1.
    pub verbose: bool,
    /// Suppress non-error output.
    pub quiet: bool,
    /// Fail on degenerate mesh faces instead of repairing them.
    pub strict_mesh: bool,
    /// Worker threads for stage 2's fit competition; `None` = all cores.
    pub threads: Option<usize>,
    /// Enable debug output and intermediate files.
    pub debug: bool,
    /// Stop after this stage.
    pub stage: Stage,
    /// Stages to visualize interactively.
    pub viz_stages: Vec<Stage>,
    /// Detail level for viz-stages visualization (1=commits, 2=seeds, 3=BFS faces).
    pub viz_stages_level: u8,
}

impl std::fmt::Debug for Config {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Config")
            .field("input_stl", &self.input_stl)
            .field("output_step", &self.output_step)
            .field("stl_units", &self.stl_units)
            .field("step_units", &self.step_units)
            .field("compare_step", &self.compare_step)
            .field("compare_shape", &if self.compare_shape.is_some() { format!("[{} faces]", self.compare_face_count) } else { "None".to_string() })
            .field("vertex_tolerance_mm", &self.vertex_tolerance_mm)
            .field("surface_tolerance_mm", &self.surface_tolerance_mm)
            .field("angular_tolerance_rad", &self.angular_tolerance_rad)
            .field("verbosity", &self.verbosity)
            .field("verbose", &self.verbose)
            .field("quiet", &self.quiet)
            .field("strict_mesh", &self.strict_mesh)
            .field("threads", &self.threads)
            .field("debug", &self.debug)
            .field("stage", &self.stage)
            .field("viz_stages", &self.viz_stages)
            .field("viz_stages_level", &self.viz_stages_level)
            .finish()
    }
}

impl Config {
    fn from_cli(args: CliArgs) -> Self {
        let scale = args.stl_units.to_mm();
        Self {
            input_stl: args.input_stl,
            output_step: args.output_step,
            stl_units: args.stl_units,
            step_units: args.step_units,
            compare_step: args.compare_step,
            compare_shape: None,
            compare_face_count: 0,
            reference_report: args.reference_report,
            invariant_report: args.invariant_report || args.viz_seeds,
            viz_seeds: args.viz_seeds,
            viz_seed_show: args.viz_seed_show,
            viz_seed_sizes: args.viz_seed_sizes,
            viz_seed_samples: args.viz_seed_samples,
            viz_seed_max_attempts: args.viz_seed_max_attempts,
            viz_seed_coplanar_clause: args.viz_seed_coplanar_clause,
            viz_seed_confine_region: args.viz_seed_confine_region,
            viz_seed_spec_only: args.viz_seed_spec_only,
            vertex_tolerance_mm: args.vertex_tolerance * scale,
            surface_tolerance_mm: args.surface_tolerance * scale,
            angular_tolerance_rad: args.angular_tolerance.to_radians(),
            verbosity: args.verbose_level.unwrap_or(0).max(args.verbose_short),
            verbose: args.verbose_level.unwrap_or(0).max(args.verbose_short) >= 1,
            quiet: args.quiet,
            strict_mesh: args.strict_mesh,
            threads: args.threads.map(|n| n as usize),
            debug: args.debug,
            stage: args.stage.unwrap_or_default(),
            viz_stages: args.viz_stages.unwrap_or_default(),
            viz_stages_level: args.viz_stages_level,
        }
    }

    /// Check if a specific stage should be visualized.
    pub fn viz_active(&self, major: u8, minor: u8) -> bool {
        self.viz_stages.iter().any(|s| {
            let (sm, ss) = s.normalized();
            sm == major && ss == minor
        })
    }

    /// If --compare was specified, load the STEP file shape for bounded distance comparison.
    pub fn load_compare_step(&mut self) -> Result<(), Box<dyn std::error::Error>> {
        let step_path = match &self.compare_step {
            Some(p) => p.clone(),
            None => return Ok(()),
        };

        let step_shape = crate::read_step_file(&step_path)
            .map_err(|e| -> Box<dyn std::error::Error> { e.into() })?;

        // Count faces for validation and status output
        let mut face_count = 0_usize;
        let mut explorer = top_exp::Explorer::new_shape_shapeenum2(
            &step_shape,
            top_abs::ShapeEnum::Face,
            top_abs::ShapeEnum::Shape,
        );
        while explorer.more() {
            face_count += 1;
            explorer.next();
        }

        if face_count == 0 {
            return Err(format!("STEP file {step_path} contains no faces").into());
        }

        if !self.quiet {
            eprintln!("Loaded {} faces from comparison STEP file", face_count);
        }

        self.compare_face_count = face_count;
        self.compare_shape = Some(step_shape);

        Ok(())
    }
}
/// Parse command-line arguments and return a resolved Config.
pub fn parse_config() -> Config {
    Config::from_cli(CliArgs::parse())
}

/// Parse from an explicit argument iterator (for testing or embedding).
pub fn parse_config_from<I, T>(args: I) -> Result<Config, clap::Error>
where
    I: IntoIterator<Item = T>,
    T: Into<std::ffi::OsString> + Clone,
{
    CliArgs::try_parse_from(args).map(Config::from_cli)
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    fn parse(args: &[&str]) -> Config {
        let mut full = vec!["brepper"];
        full.extend_from_slice(args);
        parse_config_from(full).expect("parse should succeed")
    }

    #[test]
    fn minimal_args() {
        let config = parse(&["input.stl"]);
        assert_eq!(config.input_stl, "input.stl");
        assert_eq!(config.stage, Stage(4, 1));
        assert!(!config.verbose);
    }

    #[test]
    fn full_args() {
        let config = parse(&[
            "input.stl", "-o", "output.step", "--compare", "ref.step",
            "--stage=2.3", "-v", "--vertex-tolerance=1e-6",
        ]);
        assert_eq!(config.input_stl, "input.stl");
        assert_eq!(config.output_step.as_deref(), Some("output.step"));
        assert_eq!(config.compare_step.as_deref(), Some("ref.step"));
        assert_eq!(config.stage, Stage(2, 3));
        assert!(config.verbose);
        // vertex_tolerance is in STL units (default mm), so 1e-6 mm
        assert!((config.vertex_tolerance_mm - 1e-6).abs() < 1e-15);
    }

    #[test]
    fn tolerance_converted_from_stl_units() {
        // --stl-units=in means 1e-5 inches = 1e-5 * 25.4 mm
        let config = parse(&["input.stl", "--stl-units=in", "--vertex-tolerance=1e-5"]);
        let expected = 1e-5 * 25.4;
        assert!((config.vertex_tolerance_mm - expected).abs() < 1e-15);
    }

    #[test]
    fn stage_at_least() {
        let s = Stage(2, 3);
        assert!(s.at_least(1, 1));
        assert!(s.at_least(2, 1));
        assert!(s.at_least(2, 3));
        assert!(!s.at_least(2, 4));
        assert!(!s.at_least(3, 1));

        // Stage(2, 0) means all of stage 2 = (2, 9)
        let s2 = Stage(2, 0);
        assert!(s2.at_least(2, 6));
        assert!(s2.at_least(2, 9));
        assert!(!s2.at_least(3, 1));

        // --stage=2.7 stops before the absorb (2.8) and boundary-refinement
        // (2.9) passes.
        let s27 = Stage(2, 7);
        assert!(!s27.at_least(2, 8));
    }

    #[test]
    fn strict_mesh_flag() {
        assert!(!parse(&["input.stl"]).strict_mesh);
        assert!(parse(&["input.stl", "--strict-mesh"]).strict_mesh);
    }

    #[test]
    fn all_units_parse() {
        for unit in ["mm", "cm", "m", "in", "ft", "um"] {
            let config = parse(&["input.stl", &format!("--stl-units={unit}")]);
            // Just verify it parses without error
            assert!(config.vertex_tolerance_mm > 0.0);
        }
    }

    #[test]
    fn missing_input() {
        let args = vec!["brepper", "-v"];
        let result = parse_config_from(args);
        assert!(result.is_err());
    }
}
