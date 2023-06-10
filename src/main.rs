use clap::Parser;
use grangers::grangers::{options, Grangers};
use peak_alloc::PeakAlloc;
use polars::lazy::dsl::concat_str;
use polars::prelude::*;
use std::env;
use std::path::PathBuf;
use std::time::{Duration, Instant};
use tracing_subscriber::{filter::LevelFilter, fmt, prelude::*, EnvFilter};
use grangers::grangers::{options, Grangers};
use peak_alloc::PeakAlloc;
use polars::prelude::*;
use tracing::{warn, debug, info};
use noodles;
use anyhow::Context;

#[global_allocator]
static PEAK_ALLOC: PeakAlloc = PeakAlloc;

use clap::{builder::ArgPredicate, ArgGroup, Args, Subcommand, command, Parser};

/// The type of sequences we might include in the output reference FASTA file
/// to map against for quantification with
/// alevin-fry.
#[derive(Clone, Debug)]
pub enum SequenceType {
    /// The sequence of Spliced transcripts
    Transcript,
    /// The sequence of introns of transcripts, merged by gene
    Intronic,
    /// The sequence of gene bodies (from the first to the last base of each gene)
    GeneBody,
    /// The sequence of transcript bodies (from the first to the last base of each transcript)
    TranscriptBody,
}

impl From<&str> for SequenceType {
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

#[derive(Debug, Subcommand)]
pub enum Commands {
    /// build the (expanded) reference index
    #[command(arg_required_else_help = true)]
    MakeRef {
        /// The path to a genome fasta file.
        genome: PathBuf,
        /// The path to a gene annotation gtf file.
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
            value_parser = clap::builder::PossibleValuesParser::new(["i", "g", "t", "intronic", "gene-body", "transcript-body"]),
            hide_possible_values = true,
        )]
        augmented_sequences: Option<Vec<SequenceType>>,

        /// A flag of not including spliced transcripts in the output FASTA file. (usually there should be a good reason to do so)
        #[arg(
            long,
            display_order = 3,
        )]
        no_transcript: bool,

        /// The read length of the single-cell experiment being processed (determines flank size).
        #[arg(
            short,
            long,
            help_heading = "Intronic Sequence Options",
            display_order = 1,
            default_value_t = 91,
            requires_if("intronic", "augmented_sequences"),
        )]
        read_length: i64,

        /// Determines the amount subtracted from the read length to get the flank length.
        #[arg(
            long,
            help_heading = "Intronic Sequence Options",
            display_order = 2,
            default_value_t = 5,
            requires_if("intronic", "augmented_sequences"),
        )]
        flank_trim_length: i64,

        /// A flag indicates whether flank lengths will be considered when merging introns.
        #[arg(
            long,
            help_heading = "Intronic Sequence Options",
            display_order = 3,
            requires_if("intronic", "augmented_sequences"),
        )]
        no_flanking_merge: bool,

        /// The file name prefix of the generated output files.
        #[arg(
            long,
            default_value = "augmented_ref",
            display_order = 2,
        )]
        filename_prefix: Option<String>,

        /// A flag indicates whether identical sequences will be deduplicated.
        #[arg(
            long = "dedup",
            display_order = 1,
        )]
        dedup_seqs: bool,

        /// The path to an extra spliced sequence fasta file.
        #[arg(
            long,
            help_heading = "Extra Sequence Files",
            display_order = 3,
        )]
        extra_spliced: Option<PathBuf>,

        /// The path to an extra unspliced sequence fasta file.
        #[arg(
            // short,
            long,
            help_heading = "Extra Sequence Files",
            display_order = 3,
        )]
        extra_unspliced: Option<PathBuf>,

    },
}

/// simplifying alevin-fry workflows
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
        Commands::MakeRef { genome, genes, out_dir, augmented_sequences, no_transcript, read_length, flank_trim_length, no_flanking_merge, filename_prefix, dedup_seqs, extra_spliced, extra_unspliced } => {
            make_ref(genome, genes, out_dir, augmented_sequences, no_transcript, read_length, flank_trim_length, no_flanking_merge, filename_prefix, dedup_seqs, extra_spliced, extra_unspliced)
        }?
    }



    // let ref_typ = ReferenceType::SplicedUnspliced;

    // let args: Vec<String> = env::args().collect();
    // let gtf_file = PathBuf::from(args.get(1).unwrap());
    // let fasta_file = PathBuf::from(args.get(2).unwrap());
    // let out_dir = PathBuf::from(args.get(3).unwrap());

    // // create the folder if it doesn't exist
    // std::fs::create_dir_all(&out_dir)?;
    // let out_fa = out_dir.join("splici_fl86.fa");

    // // 1. we read the gtf file as grangers. This will make sure that the eight fields are there.
    // let start = Instant::now();
    // let gr = Grangers::from_gtf(gtf_file.as_path(), true)?;
    // let duration: Duration = start.elapsed();
    // info!("Built Grangers in {:?}", duration);
    // info!("Grangers shape {:?}", gr.df().shape());

    // // we get the exon df and validate it
    // // this will make sure that each exon has a valid transcript ID, and the exon numbers are valid
    // let mut exon_gr = gr.exons(None, true)?;

    // let mut fc = exon_gr.field_columns().clone();
    // let df = exon_gr.df_mut();

    // // we then make sure that the gene_id and gene_name fields are not both missing
    // if fc.gene_id().is_none() && fc.gene_name().is_none() {
    //     anyhow::bail!(
    //         "The input GTF file must have either gene_id or gene_name field. Cannot proceed"
    //     );
    // } else if fc.gene_id().is_none() {
    //     warn!("The input GTF file do not have a gene_id field. We will use gene_name as gene_id");
    //     // we get gene name and rename it to gene_id
    //     let mut gene_id = df.column(fc.gene_name().unwrap())?.clone();
    //     gene_id.rename("gene_id");
    //     fc.update("gene_id", "gene_id")?;
    //     // push to the df
    //     df.with_column(gene_id)?;
    // } else if fc.gene_name().is_none() {
    //     warn!(
    //         "The input GTF file do not have a gene_name field. We will use gene_id as gene_name."
    //     );
    //     // we get gene id and rename it to gene_name
    //     let mut gene_name = df.column(fc.gene_id().unwrap())?.clone();
    //     gene_name.rename("gene_name");
    //     fc.update("gene_name", "gene_name")?;
    //     // push to the df
    //     df.with_column(gene_name)?;
    // }

    // exon_gr.field_columns = fc;
    // let gene_id_s = exon_gr.get_column_name("gene_id", false)?.to_string();
    // let gene_id = gene_id_s.as_str();
    // let gene_name_s = exon_gr.get_column_name("gene_name", false)?.to_string();
    // let gene_name = gene_name_s.as_str();
    // let transcript_id_s = exon_gr.get_column_name("transcript_id", false)?.to_string();
    // let transcript_id = transcript_id_s.as_str();

    // // Next, we fill the missing gene_id and gene_name fields
    // if exon_gr.any_nulls(&[gene_id, gene_name], false, false)? {
    //     warn!("Found missing gene_id and/or gene_name; Imputing. If both missing, will impute using transcript_id; Otherwise, will impute using the existing one.");
    //     exon_gr.df = exon_gr
    //         .df
    //         .lazy()
    //         .with_columns([
    //             when(col(gene_id).is_null())
    //                 .then(
    //                     when(col(gene_name).is_null())
    //                         .then(col(transcript_id))
    //                         .otherwise(col(gene_name)),
    //                 )
    //                 .otherwise(col(gene_id))
    //                 .alias(gene_id),
    //             when(col(gene_name).is_null())
    //                 .then(
    //                     when(col(gene_id).is_null())
    //                         .then(col(transcript_id))
    //                         .otherwise(col(gene_id)),
    //                 )
    //                 .otherwise(col(gene_name))
    //                 .alias(gene_name),
    //         ])
    //         .collect()?;
    // }

    // // Next, we get the gene_name to id mapping
    // let mut gene_id_to_name =
    //     exon_gr
    //         .df()
    //         .select([gene_id, gene_name])?
    //         .unique(None, UniqueKeepStrategy::Any, None)?;

    // info!(
    //     "Extracting the sequence of {} transcripts",
    //     exon_gr.df().column("transcript_id")?.n_unique()?
    // );

    // // Next, we write the transcript seuqences
    // exon_gr.write_transcript_sequences(&fasta_file, &out_fa, None, true, false)?;

    // match ref_typ {
    //     ReferenceType::Spliced => {
    //         // do nothing
    //     },
    //     ReferenceType::SplicedIntronic => {
    //         // Then, we get the introns
    //         let mut intron_gr = exon_gr.introns(None, None, None, true)?;
        
    //         intron_gr.extend(86, &options::ExtendOption::Both, false)?;
        
    //         // Then, we merge the overlapping introns
    //         intron_gr = intron_gr.merge(
    //             &[intron_gr.get_column_name("gene_id", false)?],
    //             false,
    //             None,
    //             None,
    //         )?;
        
    //         intron_gr.add_order(Some(&["gene_id"]), "intron_number", Some(1), true)?;
    //         intron_gr.df = intron_gr
    //             .df
    //             .lazy()
    //             .with_column(concat_str([col("gene_id"), col("intron_number")], "-I").alias("intron_id"))
    //             .collect()?;
        
    //         intron_gr.write_sequences(
    //             &fasta_file,
    //             &out_fa,
    //             false,
    //             Some("intron_id"),
    //             options::OOBOption::Truncate,
    //             true,
    //         )?;
    //     },
    //     ReferenceType::SplicedUnspliced => {
    //         // Then, we get the introns
    //         let mut intron_gr = exon_gr.genes(None, true)?;
    //         intron_gr.write_sequences(
    //             &fasta_file,
    //             &out_fa,
    //             false,
    //             Some("gene_id"),
    //             options::OOBOption::Truncate,
    //             true,
    //         )?;
    //     },
    // };

    // let mut file = std::fs::File::create(out_dir.join("gene_id_to_name.tsv"))?;
    // CsvWriter::new(&mut file)
    //     .has_header(false)
    //     .with_delimiter(b'\t')
    //     .finish(&mut gene_id_to_name)?;



    // next, we write transcripts and unspliced/itrons
    // 2. we quit if the required attributes are not valid:
    // - if transcript_id field dosn't exist
    // - if transcript_id field exists but is null for some exon features
    // - if gene_name and gene_id both are missing,
    // 3. we warn if one of gene_name and gene_id fields is missing, and copy the existing one as the missing one.
    // 4. If gene_name and/or gene_id fields contain null values, we imputet them:
    //      - If gene_name is missing, we impute it with gene_id
    //      - If gene_id is missing, we impute it with gene_name

    // let start = Instant::now();
    // gr.get_transcript_sequences(&fasta_file, None, true)?;
    // let duration: Duration = start.elapsed();
    // info!("extract transcript sequences in {:?}", duration);

    // gr.df = gr.df.head(Some(100000));
    // let start = Instant::now();
    // gr.get_sequences(&fasta_file, false, None, options::OOBOption::Skip)?;
    // let duration: Duration = start.elapsed();
    // info!("extract first 100,000 sequences in {:?}", duration);

    // let mo = options::MergeOptions::new(&["seqname", "gene_id", "transcript_id"], false, 1)?;
    // let start = Instant::now();
    // gr.merge(&["seqname", "gene_id", "transcript_id"], false, None)?;
    // let duration: Duration = start.elapsed();
    // info!("Merged overlapping ranges in {:?}", duration);

    // let start = Instant::now();
    // gr.flank(10, options::FlankOptions::default())?;
    // let duration: Duration = start.elapsed();
    // info!("Flanked ranges in {:?}", duration);

    // let start = Instant::now();
    // gr.introns(options::IntronsBy::Gene, None, true)?;
    // let duration: Duration = start.elapsed();
    // info!("Built intron's Grangers in {:?}", duration);

    // let start = Instant::now();
    // gr.extend(10, &options::ExtendOption::Both, false)?;
    // let duration: Duration = start.elapsed();
    // info!("Extended ranges in {:?}", duration);

    // let start = Instant::now();
    // gr.build_lapper(&mo.by)?;
    // let duration: Duration = start.elapsed();
    // info!("Built interval tree in {:?}", duration);

    // 2. we quit if the required attributes are not valid:
    // - if transcript_id field dosn't exist
    // - if transcript_id field exists but is null for some exon features
    // - if gene_name and gene_id both are missing,
    // 3. we warn if one of gene_name and gene_id fields is missing, and copy the existing one as the missing one.
    // 4. If gene_name and/or gene_id fields contain null values, we imputet them:
    //      - If gene_name is missing, we impute it with gene_id
    //     - If gene_id is missing, we impute it with gene_name
    //     - If both are missing, we impute them with transcript_id

    // let df = df!(
    //     "seqname" => ["chr1", "chr1", "chr1", "chr1", "chr1", "chr2", "chr2", "chr3"],
    //     "feature_type" => ["exon", "exon", "exon", "exon", "exon", "exon", "exon", "exon"],
    //     "start" => [1i64, 12, 1, 5, 22, 1, 5, 111],
    //     "end" => [10i64, 20, 10, 20, 30, 10, 30, 1111],
    //     "strand"=> ["+", "+", "+", "+", "+", "+", "-", "."],
    //     "gene_id" => [Some("g1"), Some("g1"), Some("g2"), Some("g2"), Some("g2"), Some("g3"), Some("g4"), None],
    // )?;

    // let mut gr = Grangers::new(df, None, None, None, IntervalType::Inclusive(1), FieldColumns::default(), false).unwrap();

    // // df.column("name")?.utf8()?.set(&df.column("name")?.is_null(), Some("."))?;
    // // println!("df: {:?}", df);

    // info!("gr: {:?}", gr.df());
    // let mo = options::MergeOptions {
    //     by: vec![
    //         "seqname".to_string(),
    //         "gene_id".to_string(),
    //         "strand".to_string(),
    //     ],
    //     ignore_strand: false,
    //     slack: 1,
    // };

    // gr.drop_nulls(None)?;
    // info!("drop_nulls' gr: \n{:?}", gr.df());

    // let merged_gr = gr.merge(&["seqname", "gene_id", "strand"], false, None)?;

    // info!("merged gr: \n{:?}", merged_gr.df());

    // let mut flanked_gr = gr.flank(10, options::FlankOptions::default())?;
    // info!("flanked gr: \n{:?}", flanked_gr.df());

    // flanked_gr.extend(10, &options::ExtendOption::Both, false)?;
    // info!("extended and flanked gr: \n{:?}", flanked_gr.df());

    // library(GenomicRanges)
    // gr <- GRanges(
    //     seqname = Rle(rep("chr1", 9)),
    //     ranges = IRanges(c(101, 101, 101, 121,141, 201, 201, 201, 221), end = c(150, 150, 110, 130, 150, 250, 250, 210,250), names = head(letters, 9)),
    //     strand = Rle(strand(c(rep("+",5), rep("-", 4)))),
    //     score = 1:9,
    //     GC = seq(1, 0, length=9))
    // gr

    // pl.when().then().otherwise()
    // use schema to define data type
    // let mut schema = Schema::new();
    // Filling missing data

    // interval operation
    // merger overlapping intervals
    // find introns
    //

    // get stats
    let peak_mem = PEAK_ALLOC.peak_usage_as_gb();
    info!("Peak Memory usage was {} GB", peak_mem);
    Ok(())
}

fn make_ref(
    genome_path: PathBuf,
    gtf_path: PathBuf, 
    out_dir: PathBuf, 
    augmented_sequences: Option<Vec<SequenceType>>, 
    no_transcript: bool, 
    read_length: i64, 
    flank_trim_length: i64, 
    no_flanking_merge: bool, 
    filename_prefix: Option<String>, 
    dedup_seqs: bool, 
    extra_spliced: Option<PathBuf>, 
    extra_unspliced: Option<PathBuf>) -> anyhow::Result<()> 
{
    // if nothing to build, then exit
    if no_transcript & augmented_sequences.is_none() {
        anyhow::bail!("Nothing to build: --no-transcript is set and no augmented_sequences is provided. Cannot proceed");
    }

    // check if extra_spliced and unspliced valid
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
    let out_fa = out_dir.join("roers_ref.fa");
    let out_gid2name = out_dir.join("gene_id_to_name.tsv");
    let flank_length = read_length - flank_trim_length;

    // 1. we read the gtf file as grangers. This will make sure that the eight fields are there.
    let start = Instant::now();
    let gr = Grangers::from_gtf(gtf_path.as_path(), true)?;
    let duration: Duration = start.elapsed();
    debug!("Built Grangers in {:?}", duration);
    info!("Built Grangers object for {:?} records", gr.df().height());

    // we get the exon df and validate it
    // this will make sure that each exon has a valid transcript ID, and the exon numbers are valid
    let mut exon_gr = gr.exons(None, true)?;

    let mut fc = exon_gr.field_columns().clone();
    let df = exon_gr.df_mut();

    // we then make sure that the gene_id and gene_name fields are not both missing
    if fc.gene_id().is_none() && fc.gene_name().is_none() {
        anyhow::bail!(
            "The input GTF file must have either gene_id or gene_name field. Cannot proceed"
        );
    } else if fc.gene_id().is_none() {
        warn!("The input GTF file do not have a gene_id field. We will use gene_name as gene_id");
        // we get gene name and rename it to gene_id
        let mut gene_id = df.column(fc.gene_name().unwrap())?.clone();
        gene_id.rename("gene_id");
        fc.update("gene_id", "gene_id")?;
        // push to the df
        df.with_column(gene_id)?;
    } else if fc.gene_name().is_none() {
        warn!(
            "The input GTF file do not have a gene_name field. We will use gene_id as gene_name."
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

    // Next, we get the gene_name to id mapping
    let mut gene_id_to_name =
        exon_gr
            .df()
            .select([gene_id, gene_name])?
            .unique(None, UniqueKeepStrategy::Any, None)?;

    info!(
        "Proceed {} exon records from {} transcripts",
        exon_gr.df().height(),
        exon_gr.df().column("transcript_id")?.n_unique()?
    );

    if no_transcript {
        std::fs::File::create(&out_fa)?;

    } else {
        // Next, we write the transcript seuqences
        exon_gr.write_transcript_sequences(&genome_path, &out_fa, None, true, false)?;
    }

    // Next, we write the augmented sequences
    if let Some(augmented_sequences) = augmented_sequences {
        for seq_typ in augmented_sequences {
            match seq_typ {
                SequenceType::Intronic => {
                    // Then, we get the introns
                    let mut intron_gr = exon_gr.introns(None, None, None, true)?;
                
                    if !no_flanking_merge {
                        intron_gr.extend(flank_length, &options::ExtendOption::Both, false)?;
                    }
                
                    // Then, we merge the overlapping introns
                    intron_gr = intron_gr.merge(
                        &[intron_gr.get_column_name("gene_id", false)?],
                        false,
                        None,
                        None,
                    )?;
                
                    intron_gr.add_order(Some(&["gene_id"]), "intron_number", Some(1), true)?;
                    intron_gr.df = intron_gr
                        .df
                        .lazy()
                        .with_column(concat_str([col("gene_id"), col("intron_number")], "-I").alias("intron_id"))
                        .collect()?;
                
                    intron_gr.write_sequences(
                        &genome_path,
                        &out_fa,
                        false,
                        Some("intron_id"),
                        options::OOBOption::Truncate,
                        true,
                    )?;
                },
                SequenceType::GeneBody => {
                    // Then, we get the introns
                    let mut intron_gr = exon_gr.genes(None, true)?;
                    intron_gr.write_sequences(
                        &genome_path,
                        &out_fa,
                        false,
                        Some("gene_id"),
                        options::OOBOption::Truncate,
                        true,
                    )?;
                },
                SequenceType::TranscriptBody => {
                    // Then, we get the introns
                    let mut intron_gr = exon_gr.transcripts(None, true)?;
                    intron_gr.write_sequences(
                        &genome_path,
                        &out_fa,
                        false,
                        Some("transcript_id"),
                        options::OOBOption::Truncate,
                        true,
                    )?;
                },
                _ => {
                    anyhow::bail!("Invalid sequence type");
                }
            }
        }
    }

    let mut file = std::fs::File::create(out_gid2name)?;
    CsvWriter::new(&mut file)
        .has_header(false)
        .with_delimiter(b'\t')
        .finish(&mut gene_id_to_name)?;

    if let Some(path) = extra_spliced {
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
            })?
        );


        for result in reader.records() {
            let record = result?;
            writer
                .write_record(&record)
                .with_context(|| {
                    format!(
                        "Could not write the sequence of extra spliced sequence {} to the output file",
                        record.name()
                    )
                })?;
        }
    };

    if let Some(path) = extra_unspliced {
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
            })?
        );


        for result in reader.records() {
            let record = result?;
            writer
                .write_record(&record)
                .with_context(|| {
                    format!(
                        "Could not write the sequence of extra spliced sequence {} to the output file",
                        record.name()
                    )
                })?;
        }
    };

    Ok(())
}