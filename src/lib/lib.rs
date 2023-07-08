use anyhow::Context;
use clap::builder::{PossibleValuesParser, TypedValueParser};
use grangers::{options, Grangers};
use polars::lazy::dsl::concat_str;
use polars::prelude::*;
use itertools::Itertools;
use serde::Serialize;
use serde_json::json;
use std::collections::{HashMap, hash_map::Entry};
use std::io::{BufWriter, Write};
use std::ops::Add;
use std::path::{Path, PathBuf};
use std::time::{Duration, Instant};
use tracing::{debug, info, warn};
use oomfi::*;

//#[global_allocator]
//static PEAK_ALLOC: PeakAlloc = PeakAlloc;

use clap::Args;

/// The type of sequences we might include in the output reference FASTA file
/// to map against for quantification with
/// alevin-fry.
#[derive(Clone, Debug, Serialize)]
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

#[derive(Args, Clone, Debug, Serialize)]
pub struct AugRefOpts {
    /// The path to a genome fasta file.
    pub genome: PathBuf,
    /// The path to a gene annotation gtf/gff3 file.
    pub genes: PathBuf,
    /// The path to the output directory (will be created if it doesn't exist).
    pub out_dir: PathBuf,

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
    pub aug_type: Option<Vec<AugType>>,

    /// A flag of not including spliced transcripts in the output FASTA file. (usually there should be a good reason to do so)
    #[arg(long, display_order = 3)]
    pub no_transcript: bool,

    /// The read length of the single-cell experiment being processed (determines flank size).
    #[arg(
        short,
        long,
        help_heading = "Intronic Sequence Options",
        display_order = 1,
        default_value_t = 91,
        requires_if("intronic", "aug_type")
    )]
    pub read_length: i64,

    /// Determines the length of sequence subtracted from the read length to obtain the flank length.
    #[arg(
        long,
        help_heading = "Intronic Sequence Options",
        display_order = 2,
        default_value_t = 5,
        requires_if("intronic", "aug_type")
    )]
    pub flank_trim_length: i64,

    /// Indicates whether flank lengths will be considered when merging introns.
    #[arg(
        long,
        help_heading = "Intronic Sequence Options",
        display_order = 3,
        requires_if("intronic", "aug_type")
    )]
    pub no_flanking_merge: bool,

    /// The file name prefix of the generated output files.
    #[arg(short = 'p', long, default_value = "roers_ref", display_order = 2)]
    pub filename_prefix: String,

    /// Indicates whether identical sequences will be deduplicated.
    #[arg(long = "dedup", display_order = 1)]
    pub dedup_seqs: bool,

    /// The path to an extra spliced sequence fasta file.
    #[arg(long, help_heading = "Extra Spliced Sequence File", display_order = 3)]
    pub extra_spliced: Option<PathBuf>,

    /// The path to an extra unspliced sequence fasta file.
    #[arg(
            // short,
            long,
            help_heading = "Extra Unspliced Sequence File",
            display_order = 3,
        )]
    pub extra_unspliced: Option<PathBuf>,

    /// Denotes that the input annotation is a GFF3 (instead of GTF) file
    #[arg(long = "gff3", display_order = 4)]
    pub gff3: bool,
}

struct SeqDedup {
    seq_bloom: Bloom,
    seq_hs: HashMap::<Vec<u8>, String>,
    collisions: Vec<(String, String)>, 
}

impl SeqDedup {
    
    fn new() -> Self {
        Self {
            seq_bloom: Bloom::new(5_000_000, 0.01),
            seq_hs: HashMap::<Vec<u8>, String>::new(),
            collisions: vec![]
        }
    }

    fn callback(&mut self, rec: &noodles::fasta::Record) -> bool {
        let prob_dup = self.seq_bloom.query(rec.sequence().as_ref());
        if prob_dup {
            match self.seq_hs.entry(
                rec.sequence()
                    .get(..)
                    .expect(
                        &format!("Failed getting sequence for record {}", rec.name())
                    )
                    .to_owned()) {
                    // if we have already seen this key then add this to the list of collisions
                    Entry::Occupied(e) => {
                        self.collisions.push((e.get().to_owned(), rec.name().to_owned()));
                        return false;
                    },
                    // otherwise, associate this sequence with the given name, and write the 
                    // sequence to file
                    Entry::Vacant(ve) => {
                        ve.insert(rec.name().to_owned());
                        return true;
                    }
                }

        } else {
            self.seq_bloom.insert(rec.sequence().as_ref());
            true 
        }
    }

    fn write_duplicate_info<P: AsRef<Path>>(&mut self, out_dir: P) -> anyhow::Result<()> {

        let dupfile = std::fs::File::create(out_dir.as_ref().join("duplicate_entries.tsv"))?;
        let mut dup_writer = BufWriter::new(dupfile);

        // sort so all of the values with the same
        // retained key are adjacent
        self.collisions.sort();

        writeln!(dup_writer, "RetainedRef\tDuplicateRef\n")?; 

        for (key, group) in &self.collisions.iter().group_by(|&x| &x.0) {
            for d in group.into_iter() {
                writeln!(dup_writer, "{}\t{}", key, d.1)?;
            }
        }
        Ok(())
    }
}


pub fn make_ref(aug_ref_opts: AugRefOpts) -> anyhow::Result<()> {
    // clean this up
    let genome_path: PathBuf = aug_ref_opts.genome;
    let gtf_path: PathBuf = aug_ref_opts.genes;
    let out_dir: PathBuf = aug_ref_opts.out_dir;
    let aug_type: Option<Vec<AugType>> = aug_ref_opts.aug_type;
    let no_transcript: bool = aug_ref_opts.no_transcript;
    let read_length: i64 = aug_ref_opts.read_length;
    let flank_trim_length: i64 = aug_ref_opts.flank_trim_length;
    let no_flanking_merge: bool = aug_ref_opts.no_flanking_merge;
    let filename_prefix: String = aug_ref_opts.filename_prefix;
    let dedup_seqs: bool = aug_ref_opts.dedup_seqs;
    let extra_spliced: Option<PathBuf> = aug_ref_opts.extra_spliced;
    let extra_unspliced: Option<PathBuf> = aug_ref_opts.extra_unspliced;
    let gff3: bool = aug_ref_opts.gff3;

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
            anyhow::bail!("The extra unspliced sequence file does not exist. Cannot proceed");
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
    let flank_length = read_length - flank_trim_length;
    // let filename_prefix = format!("{}_fl{}.fa", filename_prefix, flank_length);
    let filename_prefix = format!("{}.fa", filename_prefix);
    let out_fa_path = out_dir.join(&filename_prefix);

    // we create the writer
    // we will not append because this should be the start of the file
    let fa_out_file = 
        std::fs::OpenOptions::new()
            .write(true)
            .truncate(true)
            .create(true)
            .open(&out_fa_path)
            .with_context(|| {
                format!("Could not open the output file {:?}", out_fa_path.as_os_str())
            })?;

    let mut sd = SeqDedup::new();
    let mut sd_callback = if dedup_seqs {
        Some(
            |r: &noodles::fasta::Record| {
                sd.callback(r)
            }
        )
    } else {
        None
    };

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
        // we update the field_columns
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

    // TODO: Getting the string and then making the &str are annoying, maybe we should have a better way to do this
    exon_gr.field_columns = fc;
    let gene_id_s = exon_gr.get_column_name("gene_id", false)?;
    let gene_id = gene_id_s.as_str();
    let gene_name_s = exon_gr.get_column_name("gene_name", false)?;
    let gene_name = gene_name_s.as_str();
    let transcript_id_s = exon_gr.get_column_name("transcript_id", false)?;
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

    // Next, we write the transcript sequences to file if asked
    if !no_transcript {
        // Next, we write the transcript seuqences
        exon_gr.write_transcript_sequences_with_filter(&genome_path, &fa_out_file, None, true, &mut sd_callback)?;
    }

    // Next, we write the augmented sequences
    if let Some(aug_type) = &aug_type {
        for seq_typ in aug_type {
            match seq_typ {
                AugType::Intronic => {
                    // Then, we get the introns
                    let mut intron_gr = exon_gr.introns(None, None, None, true)?;

                    info!("Processing {} intronic records", intron_gr.df().height(),);

                    // no_flanking_merge decides when the order of merge and extend
                    // if no_flanking_merge is true, we merge first, then extend
                    if no_flanking_merge {
                        // Then, we merge the overlapping introns
                        intron_gr = intron_gr.merge(
                            &[intron_gr.get_column_name(gene_id, false)?],
                            false,
                            None,
                            None,
                        )?;
                    }

                    // add flanking end to both side of each intron
                    intron_gr.extend(flank_length, &options::ExtendOption::Both, false)?;

                    // if no_flanking_merge is false, we merge after extend
                    if !no_flanking_merge {
                        // Then, we merge the overlapping introns
                        intron_gr = intron_gr.merge(
                            &[intron_gr.get_column_name(gene_id, false)?],
                            false,
                            None,
                            None,
                        )?;
                    }

                    // Then, we give a unique id to each intron
                    intron_gr.add_order(Some(&[gene_id]), "intron_number", Some(1), true)?;
                    intron_gr.df = intron_gr
                        .df
                        .lazy()
                        .with_column(
                            concat_str([col(gene_id), col("intron_number")], "-I")
                                .alias("t2g_tx_id"),
                        )
                        .collect()?;


                    intron_gr.write_sequences_with_filter(
                        &genome_path,
                        &fa_out_file,
                        false,
                        Some("t2g_tx_id"),
                        options::OOBOption::Truncate,
                        &mut sd_callback,
                    )?;

                    // we need to update the t2g mapping for introns
                    let mut intron_t2g = intron_gr.df().select(["t2g_tx_id", gene_id])?.unique(
                        None,
                        UniqueKeepStrategy::Any,
                        None,
                    )?;

                    // if we are here, we need to add a column for splice_status cuz it is an augmented reference
                    intron_t2g.with_column(Series::new(
                        "splice_status",
                        vec!["U"; intron_t2g.height()],
                    ))?;

                    // we extend the t2g mapping for introns
                    t2g_map.extend(&intron_t2g)?;
                }
                AugType::GeneBody => {
                    // first we get the range (gene body) of each gene
                    let mut gene_gr = exon_gr.genes(None, true)?;

                    // we append a -G to each sequence here to mark that they are gene-body seuqences
                    gene_gr.df = gene_gr
                        .df
                        .lazy()
                        .with_column(col(gene_id).add(lit("-G")).alias("t2g_tx_id"))
                        .collect()?;

                    // say something
                    info!("Proceed {} gene body records", gene_gr.df().height(),);

                    gene_gr.write_sequences_with_filter(
                        &genome_path,
                        &fa_out_file,
                        false,
                        Some("t2g_tx_id"),
                        options::OOBOption::Truncate,
                        &mut sd_callback,
                    )?;

                    // we need to update the t2g mapping for genes
                    let mut gene_t2g = gene_gr.df().select(["t2g_tx_id", gene_id])?.unique(
                        None,
                        UniqueKeepStrategy::Any,
                        None,
                    )?;

                    // add the splice status column
                    gene_t2g
                        .with_column(Series::new("splice_status", vec!["U"; gene_t2g.height()]))?;

                    // extend t2g_map for genes
                    t2g_map.extend(&gene_t2g)?;
                }
                AugType::TranscriptBody => {
                    // we get the range (transcript body) of each transcript
                    let mut tx_gr = exon_gr.transcripts(None, true)?;

                    // Then we append a -T to mark the sequence as transcript body
                    tx_gr.df = tx_gr
                        .df
                        .lazy()
                        .with_column(col(gene_id).add(lit("-T")).alias("t2g_tx_id"))
                        .collect()?;

                    // say something
                    info!("Proceed {} transcript body records", tx_gr.df().height(),);

                    tx_gr.write_sequences_with_filter(
                        &genome_path,
                        &fa_out_file,
                        false,
                        Some(gene_id),
                        options::OOBOption::Truncate,
                        &mut sd_callback,
                    )?;

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

    // if there are extra spliced sequences, we include them in if dedup is set, otherwise we write them to the file
    if let Some(path) = &extra_spliced {
        // create extra file reader
        let mut reader = std::fs::File::open(path)
            .map(std::io::BufReader::new)
            .map(noodles::fasta::Reader::new)?;

        // we crate the writer, and write if not dedup
        let mut writer = noodles::fasta::Writer::new(&fa_out_file);

        // if dedup, we push the records into the seq vector
        // otherwise, we write the records to the output file
        let mut names = Vec::new();
        for result in reader.records() {
            let record = result?;

            let write_record = if let Some(ref mut dup_filt) = sd_callback {
                dup_filt(&record)
            } else { 
                true
            };

            // TODO : @DongzeHe - do we want to add the name here
            // regardless of whether or not the sequence is a duplicate
            // or do we only want to push the name if the sequence is 
            // written to the file.
            names.push(record.name().to_owned().clone());

            if write_record {
                writer.write_record(&record).with_context(|| {
                    format!(
                        "Could not write the sequence of extra spliced sequence {} to the output file",
                        record.name()
                    )
                })?;
            }
        }

        // extend t2g_map for the custom spliced targets
        let mut extra_t2g = df!(
            "t2g_tx_id" => &names,
            "gene_id" => &names,
        )?;

        if aug_type.is_some() {
            // if we are here, we need to add a column for splice_status cuz it is an augmented reference
            extra_t2g.with_column(Series::new("splice_status", vec!["S"; extra_t2g.height()]))?;
        }

        t2g_map.extend(&extra_t2g)?;

        // extend gene_id_to_name for the custom spliced targets
        gene_id_to_name.extend(&df!(gene_id => &names, gene_name => &names)?)?;
    };

    if let Some(path) = &extra_unspliced {
        // create extra file reader
        let mut reader = std::fs::File::open(path)
            .map(std::io::BufReader::new)
            .map(noodles::fasta::Reader::new)?;

        let mut writer = noodles::fasta::Writer::new(&fa_out_file);

        let mut names = Vec::new();
        for result in reader.records() {
            let record = result?;

            let write_record = if let Some(ref mut dup_filt) = sd_callback {
                dup_filt(&record)
            } else { 
                true
            };

            // TODO : @DongzeHe - do we want to add the name here
            // regardless of whether or not the sequence is a duplicate
            // or do we only want to push the name if the sequence is 
            // written to the file.
            names.push(record.name().to_owned().clone());

            if write_record {
                writer.write_record(&record).with_context(|| {
                    format!(
                        "Could not write the sequence of extra spliced sequence {} to the output file",
                        record.name()
                    )
                })?;
            }
        }

        // extend t2g_map for the custom spliced targets
        let mut extra_t2g = df!(
            "t2g_tx_id" => &names,
            "gene_id" => &names,
        )?;

        if aug_type.is_some() {
            // if we are here, we need to add a column for splice_status cuz it is an augmented reference
            extra_t2g.with_column(Series::new("splice_status", vec!["U"; extra_t2g.height()]))?;
        }

        t2g_map.extend(&extra_t2g)?;

        // extend gene_id_to_name for the custom spliced targets
        gene_id_to_name.extend(&df!(gene_id => &names, gene_name => &names)?)?;
    };

    // at this point, the fasta file should be ready
    // if we are doing deduplication, we need to write out 
    // the duplicate entry information
    if dedup_seqs { sd.write_duplicate_info(&out_dir)?; } 

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

    let info_file = out_dir.join("roers_make-ref.json");

    let v = clap::crate_version!();

    let index_info = json!({
        "command" : "roers makeref",
        "roers_version" : v,
        "output_fasta": out_fa_path,
        "output_t2g": out_t2g_name,
        "output_gid2name": out_gid2name,
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

    info!("Done!");

    Ok(())
}
