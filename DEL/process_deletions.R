#!/usr/bin/env Rscript

library("tidyverse")
library("arrow")

setwd("<path_to_desired_working_directory>")

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Read in drug, protein target, and biomarker gene data, as created in Actionable/DRUGS/Create_network.R

network        <- read.delim("<appropriate_path>/NETWORK.csv", header=T, as.is=T, sep="\t")

lossgenes            <- sort(unique(network$gene[which(network$GoF_LoF=="LoF")]))
gaingenes            <- sort(unique(network$gene[which(network$GoF_LoF=="GoF")]))
lossgenes            <- sort(c(lossgenes, "TP53"))
gaingenes            <- sort(c(gaingenes, unique(network$gene[which(network$GoF_LoF=="GoF_LoF")])))
gaingenes            <- gaingenes[-which(gaingenes %in% lossgenes)]

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Upload copy number data from the ASCAT analysis

# Metadata file, including e.g. sample tumor fraction, i.e. purity, from ASCAT copy-number analysis.
# Required columns:
# patient                            - Patient ID
# sample                             - Sample ID
# purity                             - Tumor cell fraction in the sample
metadata_table       <- read.delim("<appropriate_path>/metadata.csv>", header=T, sep="\t", as.is=T)

# Database with gene- and sample-level copy number calls from ASCAT copy-number analysis.
# Required columns:
# sample                             - Sample ID
# Gene                               - Gene symbol
# ID                                 - Ensembl gene ID
# purifiedLogR                       - log-ratio of copy number of the gene-associated segment and reference copy-number (2), normalized for cancer cell fraction in the sample
# nMajor                             - Major allele copy number for the gene
# nMinor                             - Minor allele copy number for the gene
cngenes              <- open_dataset("<appropriate_path>/cn_db.parquet") |>
                            filter(Gene %in% lossgenes) |>
                            collect()

delgenes             <- cngenes %>%
                           merge(metadata_table %>% select(sample, purity, patient),
                                 by="sample") %>%
                           mutate(mkey = paste(Gene, patient, sep="_")) %>%
                           group_by(mkey) %>%
                           mutate(suspicious = any(purity<0.2) & all(which(nMajor==0) %in% which(purity<0.2))) %>% 
                           ungroup() %>%
                           filter(!suspicious) %>%
                           filter((is.na(nMajor) & purifiedLogR<(-2)) |
                                  nMajor==0) %>%
                           select(-patient, -purity, -mkey, -suspicious)
                           
write.table(delgenes, file="deletions.csv", sep="\t", col.names=T, row.names=F, quote=F)
