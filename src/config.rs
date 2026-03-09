use clap::{ArgAction, Parser, ValueEnum};
use opencascade_sys::{message, step_control, top_abs, top_exp};
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
            2 => 6,
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

    /// Fitting tolerance for vertex-to-surface distance, in STL units.
    #[arg(long, default_value = "1e-5")]
    pub vertex_tolerance: f64,

    /// Tolerance for surface-to-triangle-face offset, in STL units.
    #[arg(long, default_value = "0.4")]
    pub surface_tolerance: f64,

    /// Max dihedral angle between adjacent triangles on same surface, in degrees.
    #[arg(long, default_value = "17.5")]
    pub angular_tolerance: f64,

    /// Increment verbosity by 1 per flag (-v = 1, -vv = 2, -vvv = 3).
    #[arg(short = 'v', action = ArgAction::Count)]
    pub verbose_short: u8,

    /// Set verbosity to an explicit level: --verbose=1 (summaries),
    /// --verbose=2 (per-face details), --verbose=3 (full BFS trace).
    /// Bare --verbose sets level 1.
    #[arg(long = "verbose", value_name = "LEVEL", default_missing_value = "1", num_args = 0..=1)]
    pub verbose_level: Option<u8>,

    /// Suppress non-error output.
    #[arg(short = 'q', long)]
    pub quiet: bool,

    /// Enable debug output and intermediate files.
    #[arg(long)]
    pub debug: bool,

    /// Stop after this stage (e.g. "2.2" or "2"). Default: run all stages.
    #[arg(long, value_parser = parse_stage)]
    pub stage: Option<Stage>,

    /// Comma-separated list of stages for interactive visualization (e.g. "2.2,3.1").
    #[arg(long = "viz-stages", value_parser = parse_stage, value_delimiter = ',')]
    pub viz_stages: Option<Vec<Stage>>,
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
    /// Fitting tolerance for vertex-to-surface distance, in mm.
    pub vertex_tolerance_mm: f64,
    /// Tolerance for surface-to-triangle-face offset, in mm.
    pub surface_tolerance_mm: f64,
    /// Max dihedral angle (radians) between adjacent triangles on the same surface.
    pub angular_tolerance_rad: f64,
    /// Verbosity level (0 = quiet-ish, 1 = verbose, 2 = very verbose, 3 = BFS trace).
    pub verbosity: u8,
    /// Convenience: verbosity >= 1.
    pub verbose: bool,
    /// Suppress non-error output.
    pub quiet: bool,
    /// Enable debug output and intermediate files.
    pub debug: bool,
    /// Stop after this stage.
    pub stage: Stage,
    /// Stages to visualize interactively.
    pub viz_stages: Vec<Stage>,
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
            .field("debug", &self.debug)
            .field("stage", &self.stage)
            .field("viz_stages", &self.viz_stages)
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
            vertex_tolerance_mm: args.vertex_tolerance * scale,
            surface_tolerance_mm: args.surface_tolerance * scale,
            angular_tolerance_rad: args.angular_tolerance.to_radians(),
            verbosity: args.verbose_level.unwrap_or(0).max(args.verbose_short),
            verbose: args.verbose_level.unwrap_or(0).max(args.verbose_short) >= 1,
            quiet: args.quiet,
            debug: args.debug,
            stage: args.stage.unwrap_or_default(),
            viz_stages: args.viz_stages.unwrap_or_default(),
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

        let mut reader = step_control::Reader::new();
        reader.read_file_charptr(&step_path);
        let progress = message::ProgressRange::new();
        reader.transfer_roots(&progress);
        let step_shape = reader.one_shape();

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

        // Stage(2, 0) means all of stage 2 = (2, 6)
        let s2 = Stage(2, 0);
        assert!(s2.at_least(2, 6));
        assert!(!s2.at_least(3, 1));
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
