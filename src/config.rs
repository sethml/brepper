use clap::{Parser, ValueEnum};
use opencascade_sys::{b_rep, message, step_control, top_abs, top_exp, topo_ds};
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
            1 => 2,
            2 => 6,
            3 => 6,
            4 => 1,
            _ => 0,
        }
    }

    /// Normalize: if minor is 0, expand to last sub-stage of that major.
    fn normalized(self) -> (u8, u8) {
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
        if major < 1 || major > 4 || minor < 1 || minor > Stage::last_minor(major) {
            return Err(format!("invalid stage: {s} (valid: 1.1–4.1)"));
        }
        Ok(Stage(major, minor))
    } else {
        let major: u8 = s.parse().map_err(|_| format!("invalid stage: {s}"))?;
        if major < 1 || major > 4 {
            return Err(format!("invalid stage: {s} (valid: 1–4)"));
        }
        Ok(Stage(major, 0))
    }
}

// ---------------------------------------------------------------------------
// Units
// ---------------------------------------------------------------------------

/// Length units for STL/STEP files.
#[derive(Debug, Clone, Copy, PartialEq, ValueEnum)]
pub enum Units {
    /// Millimeters
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

impl Default for Units {
    fn default() -> Self {
        Units::Mm
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

    /// Enable verbose output.
    #[arg(short = 'v', long)]
    pub verbose: bool,

    /// Suppress non-error output.
    #[arg(short = 'q', long)]
    pub quiet: bool,

    /// Enable debug output and intermediate files.
    #[arg(long)]
    pub debug: bool,

    /// Stop after this stage (e.g. "2.2" or "2"). Default: run all stages.
    #[arg(long, value_parser = parse_stage)]
    pub stage: Option<Stage>,
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
    /// Surfaces extracted from the --compare STEP file, for comparison checks.
    pub compare_surfaces: Vec<opencascade_sys::OwnedPtr<opencascade_sys::shape_extend::HandleGeomSurface>>,
    /// Fitting tolerance for vertex-to-surface distance, in mm.
    pub vertex_tolerance_mm: f64,
    /// Tolerance for surface-to-triangle-face offset, in mm.
    pub surface_tolerance_mm: f64,
    /// Enable verbose output.
    pub verbose: bool,
    /// Suppress non-error output.
    pub quiet: bool,
    /// Enable debug output and intermediate files.
    pub debug: bool,
    /// Stop after this stage.
    pub stage: Stage,
}

impl std::fmt::Debug for Config {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Config")
            .field("input_stl", &self.input_stl)
            .field("output_step", &self.output_step)
            .field("stl_units", &self.stl_units)
            .field("step_units", &self.step_units)
            .field("compare_step", &self.compare_step)
            .field("compare_surfaces", &format!("[{} surfaces]", self.compare_surfaces.len()))
            .field("vertex_tolerance_mm", &self.vertex_tolerance_mm)
            .field("surface_tolerance_mm", &self.surface_tolerance_mm)
            .field("verbose", &self.verbose)
            .field("quiet", &self.quiet)
            .field("debug", &self.debug)
            .field("stage", &self.stage)
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
            compare_surfaces: Vec::new(),
            vertex_tolerance_mm: args.vertex_tolerance * scale,
            surface_tolerance_mm: args.surface_tolerance * scale,
            verbose: args.verbose,
            quiet: args.quiet,
            debug: args.debug,
            stage: args.stage.unwrap_or_default(),
        }
    }

    /// If --compare was specified, load the STEP file and extract surfaces for comparison.
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

        let mut explorer = top_exp::Explorer::new_shape_shapeenum2(
            &step_shape,
            top_abs::ShapeEnum::Face,
            top_abs::ShapeEnum::Shape,
        );
        while explorer.more() {
            let shape_ref = explorer.value();
            let face = topo_ds::face_shape(shape_ref);
            self.compare_surfaces.push(b_rep::Tool::surface_face(face));
            explorer.next();
        }

        if self.compare_surfaces.is_empty() {
            return Err(format!("STEP file {step_path} contains no faces").into());
        }

        if !self.quiet {
            eprintln!("Loaded {} surfaces from comparison STEP file", self.compare_surfaces.len());
        }

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
