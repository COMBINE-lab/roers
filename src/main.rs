use anyhow::Context;
use clap::builder::{PossibleValuesParser, TypedValueParser};
use grangers::grangers::{options, Grangers};
use noodles;
use peak_alloc::PeakAlloc;
use polars::lazy::dsl::concat_str;
use polars::prelude::*;
use serde_json::json;
use std::collections::HashSet;
use std::path::PathBuf;
use std::time::{Duration, Instant};
use tracing::{debug, info, warn};
use tracing_subscriber::{filter::LevelFilter, fmt, prelude::*, EnvFilter};
use std::ops::Add;
#[global_allocator]
static PEAK_ALLOC: PeakAlloc = PeakAlloc;

use clap::{command, Parser, Subcommand};

/// The type of sequences we might include in the output reference FASTA file
/// to map against for quantification with
/// alevin-fry.
#[derive(Clone, Debug)]
pub enum AugType {
    /// The sequence of Spliced transcripts
    Transcript,
    /// The sequence of introns of transcripts, merged by gene
    Intronic,
    /// The sequence of gene bodies (from the first to the last base of each gene)
    GeneBody,
    /// The sequence of transcript bodies (from the first to the last base of each transcript)
    TranscriptBody,
}

impl From<&str> for AugType {
    fn from(s: &str) -> Self {
        match s {
            // "transcript" | "t" => Self::Transcript,
            "intronic" | "i" => Self::Intronic,
            "gene-body" | "g" => Self::GeneBody,
            "transcript-body" | "t" => Self::TranscriptBody,
            _ => panic!("Invalid sequence type"),
        }
    }
}

impl AsRef<str> for AugType {
    fn as_ref(&self) -> &str {
        match self {
            Self::Transcript => "transcript",
            Self::Intronic => "intronic",
            Self::GeneBody => "gene-body",
            Self::TranscriptBody => "transcript-body",
        }
    }
}

#[derive(Debug, Subcommand)]
pub enum Commands {
    /// build the (expanded) reference index
    #[command(arg_required_else_help = true)]
    MakeRef {
        /// The path to a genome fasta file.
        genome: PathBuf,
        /// The path to a gene annotation gtf/gff3 file.
        genes: PathBuf,
        /// The path to the output directory (will be created if it doesn't exist).
        out_dir: PathBuf,

        /// Comma separated types of augmented sequences to include in the output FASTA file on top of spliced transcripts.
        /// Available options are `intronic` (or `i` for short), `gene-body` (or `g`), and `transcript-body` (or `t`).
        #[arg(
            short,
            long,
            display_order = 1,
            value_delimiter = ',',
            requires = "genome",
            value_parser = PossibleValuesParser::new(&["i", "g", "t", "intronic", "gene-body", "transcript-body"]).map(|s| AugType::from(s.as_str())),
            hide_possible_values = true,
        )]
        aug_type: Option<Vec<AugType>>,

        /// A flag of not including spliced transcripts in the output FASTA file. (usually there should be a good reason to do so)
        #[arg(long, display_order = 3)]
        no_transcript: bool,

        /// The read length of the single-cell experiment being processed (determines flank size).
        #[arg(
            short,
            long,
            help_heading = "Intronic Sequence Options",
            display_order = 1,
            default_value_t = 91,
            requires_if("intronic", "aug_type")
        )]
        read_length: i64,

        /// Determines the length of sequence subtracted from the read length to obtain the flank length.
        #[arg(
            long,
            help_heading = "Intronic Sequence Options",
            display_order = 2,
            default_value_t = 5,
            requires_if("intronic", "aug_type")
        )]
        flank_trim_length: i64,

        /// Indicates whether flank lengths will be considered when merging introns.
        #[arg(
            long,
            help_heading = "Intronic Sequence Options",
            display_order = 3,
            requires_if("intronic", "aug_type")
        )]
        no_flanking_merge: bool,

        /// The file name prefix of the generated output files.
        #[arg(short = 'p', long, default_value = "roers_ref", display_order = 2)]
        filename_prefix: String,

        /// Indicates whether identical sequences will be deduplicated.
        #[arg(long = "dedup", display_order = 1)]
        dedup_seqs: bool,

        /// The path to an extra spliced sequence fasta file.
        #[arg(long, help_heading = "Extra Spliced Sequence File", display_order = 3)]
        extra_spliced: Option<PathBuf>,

        /// The path to an extra unspliced sequence fasta file.
        #[arg(
            // short,
            long,
            help_heading = "Extra Unspliced Sequence File",
            display_order = 3,
        )]
        extra_unspliced: Option<PathBuf>,

        /// Denotes that the input annotation is a GFF3 (instead of GTF) file
        #[arg(long = "gff3", display_order = 4)]
        gff3: bool,
    },
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
        Commands::MakeRef {
            genome,
            genes,
            out_dir,
            aug_type,
            no_transcript,
            read_length,
            flank_trim_length,
            no_flanking_merge,
            filename_prefix,
            dedup_seqs,
            extra_spliced,
            extra_unspliced,
            gff3,
        } => {
            make_ref(
                genome,
                genes,
                out_dir,
                aug_type,
                no_transcript,
                read_length,
                flank_trim_length,
                no_flanking_merge,
                filename_prefix,
                dedup_seqs,
                extra_spliced,
                extra_unspliced,
                gff3,
            )
        }?,
    }

    // get stats
    let peak_mem = PEAK_ALLOC.peak_usage_as_gb();
    debug!("Peak Memory usage was {} GB", peak_mem);
    Ok(())
}

fn make_ref(
    genome_path: PathBuf,
    gtf_path: PathBuf,
    out_dir: PathBuf,
    aug_type: Option<Vec<AugType>>,
    no_transcript: bool,
    read_length: i64,
    flank_trim_length: i64,
    no_flanking_merge: bool,
    filename_prefix: String,
    dedup_seqs: bool,
    extra_spliced: Option<PathBuf>,
    extra_unspliced: Option<PathBuf>,
    gff3: bool,
) -> anyhow::Result<()> {
    // if nothing to build, then exit
    if no_transcript & aug_type.is_none() {
        anyhow::bail!("Nothing to build: --no-transcript is set and --aug-type is not provided. Cannot proceed");
    }

    // check if extra_spliced and extra_unspliced is valid
    if let Some(extra_spliced) = &extra_spliced {
        if !extra_spliced.exists() {
            anyhow::bail!("The extra spliced sequence file does not exist. Cannot proceed");
        }
    }

    if let Some(extra_unspliced) = &extra_unspliced {
        if !extra_unspliced.exists() {
            anyhow::bail!("The extra spliced sequence file does not exist. Cannot proceed");
        }
    }

    // create the folder if it doesn't exist
    std::fs::create_dir_all(&out_dir)?;
    // specify output file names
    let out_gid2name = out_dir.join("gene_id_to_name.tsv");
    let out_t2g_name = out_dir.join(if aug_type.is_some() {
        "t2g_3col.tsv"
    } else {
        "t2g.tsv"
    });

    if flank_trim_length > read_length {
        anyhow::bail!(
            "The read length: {} must be >= the flank trim length: {}",
            read_length,
            flank_trim_length
        );
    }
    let flank_length = (read_length - flank_trim_length) as i64;
    // let filename_prefix = format!("{}_fl{}.fa", filename_prefix, flank_length);
    let filename_prefix = format!("{}.fa", filename_prefix);
    let out_fa = out_dir.join(&filename_prefix);
    std::fs::File::create(&out_fa)?;

    // 1. we read the gtf/gff3 file as grangers. This will make sure that the eight fields are there.
    let start = Instant::now();
    let (gr, file_type) = if gff3 {
        (Grangers::from_gff(gtf_path.as_path(), true)?, "GFF3")
    } else {
        (Grangers::from_gtf(gtf_path.as_path(), true)?, "GTF")
    };

    let duration: Duration = start.elapsed();
    debug!("Built the Grangers object in {:?}", duration);
    info!(
        "Built the Grangers object for {:?} records",
        gr.df().height()
    );

    // we get the exon df and validate it
    // this will make sure that each exon has a valid transcript ID, and the exon numbers are valid
    let mut exon_gr = gr.exons(None, true)?;

    let mut fc = exon_gr.field_columns().clone();
    let df = exon_gr.df_mut();

    // we then make sure that the gene_id and gene_name fields are not both missing
    if fc.gene_id().is_none() && fc.gene_name().is_none() {
        anyhow::bail!(
            "The input {} file does not have a valid gene_id or gene_name field. Cannot proceed",
            file_type
        );
    } else if fc.gene_id().is_none() {
        warn!(
            "The input {} file does not have a valid gene_id field. Roers will use gene_name as gene_id",
            file_type
        );
        // we get gene name and rename it to gene_id
        let mut gene_id = df.column(fc.gene_name().unwrap())?.clone();
        gene_id.rename("gene_id");
        fc.update("gene_id", "gene_id")?;
        // push to the df·
        df.with_column(gene_id)?;
    } else if fc.gene_name().is_none() {
        warn!(
            "The input {} file does not have a valid gene_name field. Roers will use gene_id as gene_name.",
            file_type
        );
        // we get gene id and rename it to gene_name
        let mut gene_name = df.column(fc.gene_id().unwrap())?.clone();
        gene_name.rename("gene_name");
        fc.update("gene_name", "gene_name")?;
        // push to the df
        df.with_column(gene_name)?;
    }

    exon_gr.field_columns = fc;
    let gene_id_s = exon_gr.get_column_name("gene_id", false)?.to_string();
    let gene_id = gene_id_s.as_str();
    let gene_name_s = exon_gr.get_column_name("gene_name", false)?.to_string();
    let gene_name = gene_name_s.as_str();
    let transcript_id_s = exon_gr.get_column_name("transcript_id", false)?.to_string();
    let transcript_id = transcript_id_s.as_str();

    // Next, we fill the missing gene_id and gene_name fields
    if exon_gr.any_nulls(&[gene_id, gene_name], false, false)? {
        warn!("Found missing gene_id and/or gene_name; Imputing. If both missing, will impute using transcript_id; Otherwise, will impute using the existing one.");
        exon_gr.df = exon_gr
            .df
            .lazy()
            .with_columns([
                when(col(gene_id).is_null())
                    .then(
                        when(col(gene_name).is_null())
                            .then(col(transcript_id))
                            .otherwise(col(gene_name)),
                    )
                    .otherwise(col(gene_id))
                    .alias(gene_id),
                when(col(gene_name).is_null())
                    .then(
                        when(col(gene_id).is_null())
                            .then(col(transcript_id))
                            .otherwise(col(gene_id)),
                    )
                    .otherwise(col(gene_name))
                    .alias(gene_name),
            ])
            .collect()?;
    }

    // to this point, we have a valid exon df to work with.
    info!(
        "Proceed {} exon records from {} transcripts",
        exon_gr.df().height(),
        exon_gr.df().column("transcript_id")?.n_unique()?
    );

    // Next, we get the gene id to name mapping
    let mut gene_id_to_name =
        exon_gr
            .df()
            .select([gene_id, gene_name])?
            .unique(None, UniqueKeepStrategy::Any, None)?;

    // also, the t2g mapping for spliced transcripts
    let mut t2g_map = exon_gr.df().select([transcript_id, gene_id])?.unique(
        None,
        UniqueKeepStrategy::Any,
        None,
    )?;
    t2g_map.rename(transcript_id, "t2g_tx_id")?;

    // if we have augmented sequences, we need three columns
    if aug_type.is_some() {
        t2g_map.with_column(Series::new("splice_status", vec!["S"; t2g_map.height()]))?;
    }

    // Next, we write the transcript sequences to file if asked, otherwise create a file
    // we need to know if we want to deduplicate sequences
    let mut spliced_recs: Vec<noodles::fasta::Record> = Vec::new();
    if !no_transcript {
        if dedup_seqs {
            spliced_recs.extend(
                exon_gr
                    .get_transcript_sequences(&genome_path, None, true)?
                    .into_iter(),
            );
        } else {
            // Next, we write the transcript seuqences
            exon_gr.write_transcript_sequences(&genome_path, &out_fa, None, true, true)?;
        }
    }

    // Next, we write the augmented sequences
    let mut unspliced_recs: Vec<Option<noodles::fasta::Record>> = Vec::new();
    if let Some(aug_type) = &aug_type {
        for seq_typ in aug_type {
            match seq_typ {
                AugType::Intronic => {
                    // Then, we get the introns
                    let mut intron_gr = exon_gr.introns(None, None, None, true)?;

                    info!(
                        "Processing {} intronic records",
                        intron_gr.df().height(),
                    );

                    if !no_flanking_merge {
                        intron_gr.extend(
                            flank_length as i64,
                            &options::ExtendOption::Both,
                            false,
                        )?;
                    }

                    // Then, we merge the overlapping introns
                    intron_gr = intron_gr.merge(
                        &[intron_gr.get_column_name(gene_id, false)?],
                        false,
                        None,
                        None,
                    )?;

                    intron_gr.add_order(Some(&[gene_id]), "intron_number", Some(1), true)?;
                    intron_gr.df = intron_gr
                        .df
                        .lazy()
                        .with_column(
                            concat_str([col(gene_id), col("intron_number")], "-I")
                                .alias("t2g_tx_id"),
                        )
                        .collect()?;

                    if dedup_seqs {
                        unspliced_recs.extend(
                            intron_gr
                                .get_sequences(
                                    &genome_path,
                                    false,
                                    Some("t2g_tx_id"),
                                    options::OOBOption::Truncate,
                                )?
                                .into_iter(),
                        );
                    } else {
                        intron_gr.write_sequences(
                            &genome_path,
                            &out_fa,
                            false,
                            Some("t2g_tx_id"),
                            options::OOBOption::Truncate,
                            true,
                        )?;
                    }

                    // we need to update the t2g mapping for introns
                    let mut intron_t2g = intron_gr.df().select(["t2g_tx_id", gene_id])?.unique(
                        None,
                        UniqueKeepStrategy::Any,
                        None,
                    )?;
                    // intron_t2g.rename("intron_id", transcript_id)?;
                    intron_t2g.with_column(Series::new(
                        "splice_status",
                        vec!["U"; intron_t2g.height()],
                    ))?;

                    t2g_map.extend(&intron_t2g)?;
                }
                AugType::GeneBody => {
                    // Then, we get the introns
                    let mut gene_gr = exon_gr.genes(None, true)?;

                    gene_gr.df = gene_gr
                    .df
                    .lazy()
                    .with_column(
                        col(gene_id).add(lit("-G")).alias("t2g_tx_id"),
                    )
                    .collect()?;

                    // to this point, we have a valid exon df to work with.
                    info!(
                        "Proceed {} gene body records",
                        gene_gr.df().height(),
                    );

                    if dedup_seqs {
                        unspliced_recs.extend(
                            gene_gr
                                .get_sequences(
                                    &genome_path,
                                    false,
                                    Some("t2g_tx_id"),
                                    options::OOBOption::Truncate,
                                )?
                                .into_iter(),
                        );
                    } else {
                        gene_gr.write_sequences(
                            &genome_path,
                            &out_fa,
                            false,
                            Some("t2g_tx_id"),
                            options::OOBOption::Truncate,
                            true,
                        )?;
                    }

                    // we need to update the t2g mapping for genes
                    let mut gene_t2g = gene_gr.df().select(["t2g_tx_id", gene_id])?.unique(
                        None,
                        UniqueKeepStrategy::Any,
                        None,
                    )?;
                    
                    gene_t2g.with_column(Series::new("splice_status", vec!["U"; gene_t2g.height()]))?;

                    t2g_map.extend(&gene_t2g)?;

                }
                AugType::TranscriptBody => {
                    // Then, we get the introns
                    let mut tx_gr = exon_gr.transcripts(None, true)?;

                    tx_gr.df = tx_gr
                    .df
                    .lazy()
                    .with_column(
                        col(gene_id).add(lit("-T")).alias("t2g_tx_id"),
                    )
                    .collect()?;

                    // to this point, we have a valid exon df to work with.
                    info!(
                        "Proceed {} transcript body records",
                        tx_gr.df().height(),
                    );

                    if dedup_seqs {
                        unspliced_recs.extend(
                            tx_gr
                                .get_sequences(
                                    &genome_path,
                                    false,
                                    Some(gene_id),
                                    options::OOBOption::Truncate,
                                )?
                                .into_iter(),
                        );
                    } else {
                        tx_gr.write_sequences(
                            &genome_path,
                            &out_fa,
                            false,
                            Some(gene_id),
                            options::OOBOption::Truncate,
                            true,
                        )?;
                    }

                    // we need to update the t2g mapping for genes
                    let mut tx_t2g = tx_gr.df().select(["t2g_tx_id", gene_id])?.unique(
                        None,
                        UniqueKeepStrategy::Any,
                        None,
                    )?;
                    tx_t2g.with_column(Series::new("splice_status", vec!["U"; tx_t2g.height()]))?;
                    t2g_map.extend(&tx_t2g)?;
                }
                _ => {
                    // todo, this is to avoid initalization error
                    // we probably want to set this as an empty grangers
                    // object here.
                    anyhow::bail!("Invalid sequence type");
                }
            }
        }
    }

    // if there are extra spliced sequences, we include them in
    if let Some(path) = &extra_spliced {
        // create extra file reader
        let mut reader = std::fs::File::open(path)
            .map(std::io::BufReader::new)
            .map(noodles::fasta::Reader::new)?;

        // we crate the writer, and write if not dedup
        let mut writer = noodles::fasta::Writer::new(
            std::fs::OpenOptions::new()
                .append(true)
                .open(&out_fa)
                .with_context(|| {
                    format!("Could not open the output file {:?}", out_fa.as_os_str())
                })?,
        );

        // if dedup, we push the records into the seq vector
        // otherwise, we write the records to the output file
        let mut names = Vec::new();
        for result in reader.records() {
            let record = result?;
            names.push(record.name().to_owned().clone());

            if dedup_seqs {
                spliced_recs.push(record);
            } else {
                writer.write_record(&record).with_context(|| {
                    format!(
                        "Could not write the sequence of extra spliced sequence {} to the output file",
                        record.name()
                    )
                })?;
            }
        }

        // extend t2g_map for the custom spliced targets
        t2g_map.extend(&df!(
            "transcript_id" => &names,
            "gene_id" => &names,
            "splice_status" => &vec!["S"; names.len()]
        )?)?;

        // extend gene_id_to_name for the custom spliced targets
        gene_id_to_name.extend(&df!("gene_id" => &names, "gene_name" => &names)?)?;
    };

    if let Some(path) = &extra_unspliced {
        // create extra file reader
        let mut reader = std::fs::File::open(path)
            .map(std::io::BufReader::new)
            .map(noodles::fasta::Reader::new)?;

        let mut writer = noodles::fasta::Writer::new(
            std::fs::OpenOptions::new()
                .append(true)
                .open(&out_fa)
                .with_context(|| {
                    format!("Could not open the output file {:?}", out_fa.as_os_str())
                })?,
        );

        let mut names = Vec::new();
        for result in reader.records() {
            let record = result?;
            names.push(record.name().to_owned().clone());

            if dedup_seqs {
                unspliced_recs.push(Some(record));
            } else {
                writer.write_record(&record).with_context(|| {
                    format!(
                        "Could not write the sequence of extra spliced sequence {} to the output file",
                        record.name()
                    )
                })?;
            }
        }

        // extend a dataframe for the custom unspliced targets
        t2g_map.extend(&df!(
            "transcript_id" => &names,
            "gene_id" => &names,
            "splice_status" => &vec!["U"; names.len()]
        )?)?;

        // extend gene_id_to_name for the custom unspliced targets
        gene_id_to_name.extend(&df!("gene_id" => &names, "gene_name" => &names)?)?;
    };

    // at this point, if we don't do deduplication, the fasta file should be ready
    // if we do, we need to check if each seq is the unique one before write to the file
    if dedup_seqs {
        let mut writer = noodles::fasta::Writer::new(
            std::fs::OpenOptions::new()
                .append(true)
                .open(&out_fa)
                .with_context(|| {
                    format!("Could not open the output file {:?}", out_fa.as_os_str())
                })?,
        );

        // we need a hashset to record if we have seen a sequence
        let mut seq_hs = HashSet::new();
        // we process spliced records first
        if !spliced_recs.is_empty() {
            // for each record, if the sequence is not in the hashset, we write it to the file
            for r in spliced_recs.iter() {
                if seq_hs.insert(
                    r.sequence()
                        .get(..)
                        .with_context(|| {
                            format!("Failed getting sequence for record {}", r.name())
                        })?
                        .to_owned(),
                ) {
                    writer.write_record(r).with_context(|| {
                        format!(
                            "Could not write the sequence of spliced transcript {} to the output file",
                            r.name()
                        )
                    })?;
                }
            }
        }

        // we process unspliced records first
        if !unspliced_recs.is_empty() {
            // for each record, if the sequence is not in the hashset, we write it to the file
            for r in unspliced_recs.iter() {
                // as we might get empty sequence (oob)
                // we need to check if the record is empty
                if let Some(r) = r {
                    if seq_hs.insert(
                        r.sequence()
                            .get(..)
                            .with_context(|| {
                                format!("Failed getting sequence for record {}", r.name())
                            })?
                            .to_owned(),
                    ) {
                        writer.write_record(r).with_context(|| {
                            format!(
                                "Could not write the sequence of spliced transcript {} to the output file",
                                r.name()
                            )
                        })?;
                    }
                }
            }
        }
    } // if dedup_seq

    // Till this point, we are done with the fasta file
    // we need to write the t2g and gene_id_to_name files

    let mut file = std::fs::File::create(&out_t2g_name)?;
    CsvWriter::new(&mut file)
        .has_header(false)
        .with_delimiter(b'\t')
        .finish(&mut t2g_map)?;

    let mut file = std::fs::File::create(&out_gid2name)?;
    CsvWriter::new(&mut file)
        .has_header(false)
        .with_delimiter(b'\t')
        .finish(&mut gene_id_to_name)?;

    let info_file = out_dir.join("roers_makeref.json");

    let v = clap::crate_version!();

    let index_info = json!({
        "command" : "roers makeref",
        "roers_version" : v,
        "output_fasta": out_fa,
        "args" : {
            "genome" : genome_path,
            "genes" : gtf_path,
            "out_dir" : out_dir,
            "aug_type" : if let Some(at) = &aug_type {
                at.iter().map(|x| x.as_ref()).collect::<Vec<&str>>()
            } else {
                vec![]
            },
            "no_transcript" : no_transcript,
            "read_length" : read_length,
            "flank_trim_length" : flank_trim_length,
            "no_flanking_merge" : no_flanking_merge,
            "filename_prefix" : filename_prefix,
            "dedup_seqs" : dedup_seqs,
            "extra_spliced" : extra_spliced,
            "extra_unspliced" : extra_unspliced,
            "gff3" : gff3,
        }
    });

    std::fs::write(
        &info_file,
        serde_json::to_string_pretty(&index_info).unwrap(),
    )
    .with_context(|| format!("could not write {}", info_file.display()))?;

    Ok(())
}
