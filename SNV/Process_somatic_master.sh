#!/usr/bin/bash

#Shell script for processing somatic variants

# Source data is a tab-delimited file converted from a .vcf output of GATK variant calling.
# Required columns:
# CHROM, POS, REF, ALT, ID, FILTER   - Obligatory fields of a vcf-file (https://samtools.github.io/hts-specs/VCFv4.1.pdf)
# patient                            - Patient ID
# samples                            - Samples IDs (separated by ";")
# readCounts                         - Sample-level readCounts (sample-specific values separated by ";", REF and ALT values separated by ",")
# Func.MANE                          - Gene-based annotation by a MANE transcript: Position in relation to protein-coding sequence (e.g. "exonic", "intronic", "UTR3")
# Gene.MANE                          - Gene-based annotation by a MANE transcript: Gene symbol
# GeneDetail.MANE                    - Gene-based annotation by a MANE transcript: Transcript, cDNA change, and protein change (separated with ":")
# ExonicFunc.MANE                    - Gene-based annotation by a MANE transcript: Effect on protein-coding sequence (e.g. "frameshift_deletion", "nonsynonymous_SNV", "stopgain")
# AAChange.MANE                      - Gene-based annotation by a MANE transcript: Gene symbol, transcript, cDNA change, and protein change (separated with ":")
# Func.refGene                       - Gene-based annotation by a refGene transcript: Position in relation to protein-coding sequence (e.g. "exonic", "intronic", "UTR3")
# ExonicFunc.refGene                 - Gene-based annotation by a refGene transcript: Effect on protein-coding sequence (e.g. "frameshift_deletion", "nonsynonymous_SNV", "stopgain")
# AAChange.refGene                   - Gene-based annotation by a MANE transcript: Gene symbol, transcript, cDNA change, and protein change (separated with ":")
# CLNSIG                             - ClinVar annotation: Clinical significance (e.g. "Pathogenic", "Uncertain_significance", "Benign")
# CLNSIGCONF                         - ClinVar annotation: Certainty of clinical significance (numerical summary)
# CLNALLELEID                        - ClinVar annotation: Variant ID
# CADD_phred                         - Normalized CADD prediction value
# dbscSNV_ADA_SCORE                  - dbscSNV ADA splice prediction score
# dbscSNV_RF_SCORE                   - dbscSNV RF splice prediction score
# PolyPhenCat                        - PolyPhen2 prediction for missenses: category (e.g. "probably_damaging", "benign")
# PolyPhenVal                        - PolyPhen2 prediction for missenses: numeric score
# SIFTcat                            - SIFT prediction for missenses: category (e.g. "deleterious", "tolerated")
# SIFTval                            - SIFT prediction for missenses: numeric score
# AM_class                           - AlphaMissense prediction for missenses: category (e.g. "likely_pathogenic", "likely_benign", "ambiguous")
# AM_score                           - AlphaMissense prediction for missenses: numeric score
SOMATIC_VARIANTS="Source_variant_table_from_vcf.csv"

# Paths to sample-specific output files from GATK ASEReadCounter for calling sample-specific, somatic, protein-coding missenses in RNA-seq data.
# The path file needs to have the following columns:
# sample_RNA                         - RNA sample ID
# file                               - Path to sample-specific output file from GATK ASEReadCounter
# The AseReadCounter output file has to have the following columns
# contig                             - Chromosome
# position                           - Position in chromosome
# refCount                           - Number of reads from the reference allele
# altCount                           - Number of reads from the alternative allele
GEX_PATH_FILE="File_with_paths_to_ASEReadCounter_output.csv"

# Unique list of symbols of the target genes
# tab-delimited file with header row and the first column containing gene symbols.
GENE_LIST="Actionable/DRUGS/Actionable_gene_list.csv"
# Gene list converted into a regular expression
GENE_LIST_REGEX="(^|;)($(tail -n +2 $GENE_LIST | cut -f1 | perl -pe 's/\n/|/g' | perl -pe 's/\|$//'))(;|$)"

# Database with gene- and sample-level copy number calls from ASCAT copy-number analysis.
# Required columns:
# sample                             - Sample ID
# Gene                               - Gene symbol
# ID                                 - Ensembl gene ID
# CNstatus                           - Gene copy number status ("AMP", "DEL", "Normal")
# nMajor                             - Major allele copy number for the gene
# nMinor                             - Minor allele copy number for the gene
CN_DB="cn_db.parquet"

# Metadata file, including e.g. sample tumor fraction, i.e. purity, from ASCAT copy-number analysis.
# Required columns:
# sample                             - Sample ID
# aberrant                           - Whether sample contains tumor cells (TRUE/FALSE)
# purity                             - Tumor cell fraction in the sample
METADATA="metadata.csv"

# Perl script for filtering
FILTER_SCRIPT="Actionable/SNV/filter.somatic.pl"

# R script for processing
PROCESSING_SCRIPT="Actionable/SNV/process_somatic_variant_table.R"

# Output directory: needs an existing path
OUT_DIR="Output_directory_to_be_defined"


cd $OUT_DIR

#Get candidate gene variants
$FILTER_SCRIPT $SOMATIC_VARIANTS |
    awk -v gene_regex="$GENE_LIST_REGEX" '$10 == "Gene.MANE" || $10 ~ gene_regex' \
        > somatic.snv.coding.csv

#Process the table with purity and gene copy number status, outputting a row for each sample-gene pair
Rscript $PROCESSING_SCRIPT $GENE_LIST $GEX_PATH_FILE $CN_DB $METADATA $OUT_DIR

# CGI predictions need to be retrieved from https://www.cancergenomeinterpreter.org/home, using dataToCGI.csv as an input file.
# The CGI output needs to have the following columns:
# CHROMOSOME
# POSITION
# REF
# ALT
# CGI.Protein.Change
# CGI.Oncogenic.Summary
# CGI.Consequence
# Merge CGI predictions with "Rscript summarize_predictions.R somatic.expressed.csv <CGI_output_file.csv>"