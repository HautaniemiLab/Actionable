#!/bin/bash

#Shell script for processing somatic variants

# Source data is a .csv file converted from a .vcf output of GATK variant calling.
SOMATIC_VARIANTS="Source_variant_table_from_vcf.csv"

# Paths to sample-specific output files from GATK ASEReadCounter for calling sample-specific, somatic, protein-coding missenses in RNA-seq data.
GEX_PATH_FILE="File_including_paths_to_ASEReadCounter_output.csv"

# Unique list of symbols of the target genes
GENE_LIST="Actionable/DRUGS/Actionable_gene_list.csv"

# Gene list converted into a regular expression
GENE_LIST_REGEX="(^|;)($(tail -n +2 $GENE_LIST | cut -f1 | perl -pe 's/\n/|/g' | perl -pe 's/\|$//'))(;|$)"

# Perl script for filtering
FILTER_SCRIPT="Actionable/SNV/filter.somatic.pl"

# R script for processing
PROCESSING_SCRIPT="Actionable/SNV/process_somatic_variant_table.R"

# Output directory: an existing path
OUT_DIR="Output_directory_to_be_defined"


cd $OUT_DIR

#Get candidate gene variants
$FILTER_SCRIPT $SOMATIC_VARIANTS |
    awk -v gene_regex="$GENE_LIST_REGEX" '$10 == "Gene.MANE" || $10 ~ gene_regex' \
        > somatic.actionable.all.csv

#Process the table with purity and gene copy number status, outputting a row for each sample-gene pair
Rscript $PROCESSING_SCRIPT $GEX_PATH_FILE $OUT_DIR

