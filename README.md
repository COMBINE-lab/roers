# roers

A rust library for preparing expanded transcriptome references for quantification with [`alevin-fry`](https://alevin-fry.readthedocs.io/en/latest/).

To use outside of simpleaf, follow the following help message:

```bash
build the (expanded) reference index

Usage: roers make-ref [OPTIONS] <GENOME> <GENES> <OUT_DIR>

Arguments:
  <GENOME>   The path to a genome fasta file
  <GENES>    The path to a gene annotation gtf/gff3 file
  <OUT_DIR>  The path to the output directory (will be created if it doesn't exist)

Options:
  -a, --aug-type <AUG_TYPE>
          Comma separated types of augmented sequences to include in the output FASTA file on
          top of spliced transcripts. Available options are `intronic` (or `i` for short),
          `gene-body` (or `g`), and `transcript-body` (or `t`)
      --dedup
          Indicates whether identical sequences will be deduplicated
  -p, --filename-prefix <FILENAME_PREFIX>
          The file name prefix of the generated output files [default: roers_ref]
      --no-transcript
          A flag of not including spliced transcripts in the output FASTA file. (usually there
          should be a good reason to do so)
      --gff3
          Denotes that the input annotation is a GFF3 (instead of GTF) file
  -h, --help
          Print help
  -V, --version
          Print version

Intronic Sequence Options:
  -r, --read-length <READ_LENGTH>
          The read length of the single-cell experiment being processed (determines flank size)
          [default: 91]
      --flank-trim-length <FLANK_TRIM_LENGTH>
          Determines the length of sequence subtracted from the read length to obtain the flank
          length [default: 5]
      --no-flanking-merge
          Indicates whether flank lengths will be considered when merging introns

Extra Spliced Sequence File:
      --extra-spliced <EXTRA_SPLICED>  The path to an extra spliced sequence fasta file

Extra Unspliced Sequence File:
      --extra-unspliced <EXTRA_UNSPLICED>  The path to an extra unspliced sequence fasta file

```