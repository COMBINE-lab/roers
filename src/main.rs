use peak_alloc::PeakAlloc;
use tracing::debug;
use tracing_subscriber::{filter::LevelFilter, fmt, prelude::*, EnvFilter};

#[global_allocator]
static PEAK_ALLOC: PeakAlloc = PeakAlloc;

use clap::{command, Parser, Subcommand};

#[derive(Debug, Subcommand)]
pub enum Commands {
    /// build the (expanded) reference index
    #[command(arg_required_else_help = true)]
    MakeRef(roers::AugRefOpts),
}

#[derive(Debug, Parser)]
#[command(author, version, about)]
#[command(propagate_version = true)]
pub struct Cli {
    #[command(subcommand)]
    command: Commands,
}

fn main() -> anyhow::Result<()> {
    // Check the `RUST_LOG` variable for the logger level and
    // respect the value found there. If this environment
    // variable is not set then set the logging level to
    // INFO.
    tracing_subscriber::registry()
        .with(fmt::layer())
        .with(
            EnvFilter::builder()
                .with_default_directive(LevelFilter::INFO.into())
                .from_env_lossy(),
        )
        .init();

    let cli_args = Cli::parse();

    match cli_args.command {
        Commands::MakeRef(aug_ref_opts) => { roers::make_ref(aug_ref_opts) }?,
    }

    // get stats
    let peak_mem = PEAK_ALLOC.peak_usage_as_gb();
    debug!("Peak Memory usage was {} GB", peak_mem);
    Ok(())
}
