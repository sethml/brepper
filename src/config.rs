use std::error::Error;
use std::fmt::{Display, Formatter};

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

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Units {
    Mm,
    // TODO: Cm, M, Inch, etc.
}

impl Default for Units {
    fn default() -> Self {
        Units::Mm
    }
}

/// Configuration parsed from command-line flags.
#[derive(Debug, Clone)]
pub struct Config {
    /// Input STL file path (required).
    pub input_stl: String,
    /// Output STEP file path.
    pub output_step: Option<String>,
    /// Units of the input STL file.
    pub stl_units: Units,
    /// Units for the output STEP file.
    pub step_units: Units,
    /// STEP file to compare against at each stage (--compare).
    pub compare_step: Option<String>,
    /// Fitting tolerance for vertex-to-surface distance (default 1e-5).
    pub vertex_tolerance: f64,
    /// Tolerance for surface-to-triangle-face offset (default 0.4).
    pub surface_tolerance: f64,
    /// Enable verbose output.
    pub verbose: bool,
    /// Suppress non-error output.
    pub quiet: bool,
    /// Enable debug output and intermediate files.
    pub debug: bool,
    /// Stop after this stage.
    pub stage: Stage,
}

impl Default for Config {
    fn default() -> Self {
        Self {
            input_stl: String::new(),
            stl_units: Units::default(),
            step_units: Units::default(),
            output_step: None,
            compare_step: None,
            vertex_tolerance: 1e-5,
            surface_tolerance: 0.4,
            verbose: false,
            quiet: false,
            debug: false,
            stage: Stage::default(),
        }
    }
}

#[derive(Debug)]
pub enum ConfigError {
    MissingInputStl,
    UnknownFlag(String),
    MissingValue(String),
    InvalidStage(String),
    InvalidUnits(String),
    InvalidFloat { flag: String, value: String },
}

impl Display for ConfigError {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        match self {
            ConfigError::MissingInputStl => write!(f, "missing required argument: <input.stl>"),
            ConfigError::UnknownFlag(flag) => write!(f, "unknown flag: {flag}"),
            ConfigError::MissingValue(flag) => write!(f, "flag {flag} requires a value"),
            ConfigError::InvalidStage(s) => write!(f, "invalid --stage value: {s} (expected e.g. '2' or '2.3')"),
            ConfigError::InvalidUnits(s) => write!(f, "invalid units: {s} (expected: mm)"),
            ConfigError::InvalidFloat { flag, value } => write!(f, "invalid number for {flag}: {value}"),
        }
    }
}

impl Error for ConfigError {}

fn parse_units(s: &str) -> Result<Units, ConfigError> {
    match s.to_lowercase().as_str() {
        "mm" => Ok(Units::Mm),
        _ => Err(ConfigError::InvalidUnits(s.to_string())),
    }
}

fn parse_stage(s: &str) -> Result<Stage, ConfigError> {
    if let Some((major_s, minor_s)) = s.split_once('.') {
        let major: u8 = major_s.parse().map_err(|_| ConfigError::InvalidStage(s.to_string()))?;
        let minor: u8 = minor_s.parse().map_err(|_| ConfigError::InvalidStage(s.to_string()))?;
        if major < 1 || major > 4 || minor < 1 || minor > Stage::last_minor(major) {
            return Err(ConfigError::InvalidStage(s.to_string()));
        }
        Ok(Stage(major, minor))
    } else {
        let major: u8 = s.parse().map_err(|_| ConfigError::InvalidStage(s.to_string()))?;
        if major < 1 || major > 4 {
            return Err(ConfigError::InvalidStage(s.to_string()));
        }
        Ok(Stage(major, 0))
    }
}

fn parse_float(flag: &str, value: &str) -> Result<f64, ConfigError> {
    value.parse::<f64>().map_err(|_| ConfigError::InvalidFloat {
        flag: flag.to_string(),
        value: value.to_string(),
    })
}

/// Take the next value from args iterator, or split on '=' from the flag itself.
fn take_value<'a>(
    flag: &str,
    eq_value: Option<&'a str>,
    args: &mut impl Iterator<Item = String>,
) -> Result<String, ConfigError> {
    if let Some(v) = eq_value {
        Ok(v.to_string())
    } else {
        args.next().ok_or_else(|| ConfigError::MissingValue(flag.to_string()))
    }
}

pub fn parse_args(args: impl Iterator<Item = String>) -> Result<Config, ConfigError> {
    let mut config = Config::default();
    let mut args = args.peekable();
    let mut positional_count = 0;

    while let Some(arg) = args.next() {
        // Split --flag=value
        let (flag, eq_value) = if let Some(idx) = arg.find('=') {
            (&arg[..idx], Some(&arg[idx + 1..]))
        } else {
            (arg.as_str(), None)
        };

        match flag {
            "-o" | "--output" => {
                config.output_step = Some(take_value(flag, eq_value, &mut args)?);
            }
            "--stl-units" => {
                let v = take_value(flag, eq_value, &mut args)?;
                config.stl_units = parse_units(&v)?;
            }
            "--step-units" => {
                let v = take_value(flag, eq_value, &mut args)?;
                config.step_units = parse_units(&v)?;
            }
            "--compare" => {
                config.compare_step = Some(take_value(flag, eq_value, &mut args)?);
            }
            "--vertex-tolerance" => {
                let v = take_value(flag, eq_value, &mut args)?;
                config.vertex_tolerance = parse_float(flag, &v)?;
            }
            "--surface-tolerance" => {
                let v = take_value(flag, eq_value, &mut args)?;
                config.surface_tolerance = parse_float(flag, &v)?;
            }
            "-v" | "--verbose" => config.verbose = true,
            "-q" | "--quiet" => config.quiet = true,
            "--debug" => config.debug = true,
            "--stage" => {
                let v = take_value(flag, eq_value, &mut args)?;
                config.stage = parse_stage(&v)?;
            }
            _ if flag.starts_with('-') => {
                return Err(ConfigError::UnknownFlag(flag.to_string()));
            }
            _ => {
                // Positional argument
                if positional_count == 0 {
                    config.input_stl = arg;
                } else {
                    return Err(ConfigError::UnknownFlag(arg));
                }
                positional_count += 1;
            }
        }
    }

    if config.input_stl.is_empty() {
        return Err(ConfigError::MissingInputStl);
    }

    Ok(config)
}

pub fn print_usage(program: &str) {
    eprintln!("brepper - Convert STL mesh to STEP with fitted surfaces");
    eprintln!();
    eprintln!("USAGE:");
    eprintln!("    {program} [OPTIONS] <input.stl> [-o <output.step>]");
    eprintln!();
    eprintln!("REQUIRED:");
    eprintln!("    <input.stl>              Input STL file (binary or ASCII)");
    eprintln!();
    eprintln!("GENERAL OPTIONS:");
    eprintln!("    -o <step>, --output=<step>  Output STEP file");
    eprintln!("    --stl-units=<units>      Units used by the STL file (default: mm)");
    eprintln!("    --step-units=<units>     Units to use in exported STEP file (default: mm)");
    eprintln!("    --compare=<step>         STEP file to compare to at each step");
    eprintln!("    --vertex-tolerance=<val> Fitting tolerance in STL units (default: 1e-5)");
    eprintln!("    --surface-tolerance=<val> Surface-to-face offset tolerance (default: 0.4 mm)");
    eprintln!("    -v, --verbose            Enable verbose output");
    eprintln!("    -q, --quiet              Suppress non-error output");
    eprintln!("    --debug                  Enable debug output and intermediate files");
    eprintln!("    --stage=<stage>          Stop after stage, e.g. 2.2 (default: 4.1)");
}

#[cfg(test)]
mod tests {
    use super::*;

    fn parse(args: &[&str]) -> Result<Config, ConfigError> {
        parse_args(args.iter().map(|s| s.to_string()))
    }

    #[test]
    fn minimal_args() {
        let config = parse(&["input.stl"]).unwrap();
        assert_eq!(config.input_stl, "input.stl");
        assert_eq!(config.stage, Stage(4, 1));
        assert!(!config.verbose);
    }

    #[test]
    fn full_args() {
        let config = parse(&[
            "input.stl", "-o", "output.step", "--compare", "ref.step",
            "--stage=2.3", "-v", "--vertex-tolerance=1e-6",
        ]).unwrap();
        assert_eq!(config.input_stl, "input.stl");
        assert_eq!(config.output_step.as_deref(), Some("output.step"));
        assert_eq!(config.compare_step.as_deref(), Some("ref.step"));
        assert_eq!(config.stage, Stage(2, 3));
        assert!(config.verbose);
        assert!((config.vertex_tolerance - 1e-6).abs() < 1e-15);
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
    fn missing_input() {
        let err = parse(&["-v"]).unwrap_err();
        assert!(matches!(err, ConfigError::MissingInputStl));
    }

    #[test]
    fn unknown_flag() {
        let err = parse(&["input.stl", "--bogus"]).unwrap_err();
        assert!(matches!(err, ConfigError::UnknownFlag(_)));
    }
}
