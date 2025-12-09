#!/usr/bin/env Rscript

library("tidyverse")
library("arrow")
library("vsn")

setwd("<path_to_desired_working_directory>")

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Read in drug, protein target, and biomarker gene data, as created in Actionable/DRUGS/Create_network.R

network              <- read.delim("<appropriate_path>/NETWORK.csv", header=T, as.is=T, sep="\t")

lossgenes            <- sort(unique(network$gene[which(network$GoF_LoF=="LoF")]))
gaingenes            <- sort(unique(network$gene[which(network$GoF_LoF=="GoF")]))
lossgenes            <- sort(c(lossgenes, "TP53"))
gaingenes            <- sort(c(gaingenes, unique(network$gene[which(network$GoF_LoF=="GoF_LoF")])))
gaingenes            <- gaingenes[-which(gaingenes %in% lossgenes)]

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Read in OncoKB v4.1 data with the following columns.
# Level                                       - FDA therapeutic level for the Gene-treatment combination
# Gene                                        - HUGO gene symbol
# Alterations                                 - Gene mutations, specific mutations or mutation categories
# Cancer.Types                                - Cancer histologies
# Drugs..for.therapeutic.implications.only.   - Drug combinations

OKBdata              <- read.delim("OKBdata_v4.1.csv", header=T, as.is=T, sep="\t", na.strings=c("", "NA")) 

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Read in copy-number impact data
# Supplementary table 2 from Jamalzadeh et al. BMC Cancer, 2024. https://doi.org/10.1186/s12885-024-11895-6
# Required column names on row 6
# Gene                                        - Gene symbol
# CNI                                         - Copy-number impact (numerical correlation value)
# CNA transition point                        - Copy number transition point
# Transition.bandwidth.left                   - Lower bound of CNA transition point confidence interval
# Transition.bandwidth.right                  - Upper bound of CNA transition point confidence interval

cni                  <- read.delim("<jamalzadeh_data.csv>", header=T, as.is=T, sep=",", skip=5)

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Metadata file, including e.g. sample tumor fraction, i.e. purity, from ASCAT copy-number analysis.
# Required columns:
# patient                                     - Patient ID
# sample                                      - Sample ID
# purity                                      - Tumor cell fraction in the sample

metadata_table       <- read.delim("<appropriate_path>/metadata.csv>", header=T, sep="\t", as.is=T)

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Upload gene copy number data from a database with gene- and sample-level copy number calls from ASCAT copy-number analysis.
# Required columns:
# sample                                      - Sample ID
# Gene                                        - Gene symbol
# ID                                          - Ensembl gene ID
# purifiedLogR                                - log-ratio of copy number of the gene-associated segment and reference copy-number (2), normalized for cancer cell fraction in the sample
# nMajor                                      - Major allele copy number for the gene
# nMinor                                      - Minor allele copy number for the gene

cngenes              <- open_dataset("<appropriate_path>/cn_db.parquet")  |>
                            filter(Gene %in% gaingenes) |>
                            collect()
               

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Read in normalized gene expression matrix with gene symbols as row names and sample IDs as column names

samplegex           <- read.delim("<expression_matrix.csv>", header=T, stringsAsFactors=F, sep="\t", row.names=1)

gexanno             <- data.frame(
                           sample = str_extract(colnames(samplegex), "^[:alnum:]+(?=_)"),
                           temp   = str_extract(colnames(samplegex), "(?<=_)[poir][:alnum:]+"),
                           row.names = colnames(samplegex)) %>%
                           mutate(stime = str_extract(temp, "^[poir][:digit:]?"),
                                  ssite = str_extract(temp, "(?<=[poir][:digit:]?)([:upper:][:lower:]{2}|LN)"),
                                  sside = str_extract(temp, "(?<=([:upper:][:lower:]{2}|LN))[LR]"),
                                  sord  = str_extract(temp, "[:digit:]$")) %>%
                           mutate(prim  = str_detect(stime, "p")) %>%
                           filter(!is.na(ssite) & !is.na(prim)) %>%
                           select(-temp)

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# PROCESS DATA

ampdata <- cngenes %>%
               merge(metadata_table %>% select(sample, ploidy), all.x=T, by="sample") %>%
               mutate(AMP_ploidy = (nMajor + nMinor)/ploidy > 2.5 & purifiedLogR > log2(2.5) & (nMajor + nMinor)>=8) %>%
               mutate(patient = str_extract(sample, "^[:alnum:]+"),
                      temp = str_extract(sample, "(?<=_)[poir][:alnum:]+")) %>%
               mutate(stime = str_extract(temp, "^[poir][:digit:]?"),
                      ssite = str_extract(temp, "(?<=[poir][:digit:]?)([:upper:][:lower:]{2}|LN)"),
                      sside = str_extract(temp, "(?<=([:upper:][:lower:]{2}|LN))[LR]"),
                      sord  = str_extract(temp, "[:digit:]$")) %>%
               select(-temp) %>%
               mutate(patient_gene = paste(patient, Gene, sep="_")) %>%
               group_by(patient_gene) %>%
               filter(any(AMP_ploidy)) %>%
               ungroup()
               
               
AMP_predicted_genes <- cni %>%
                           filter(Gene %in% gaingenes & (CNA.transition.point>2 | (CNA.transition.point!=0 & Transition.bandwidth.right>5))) %>%
                           filter(CNI > quantile(CNI, 0.75, na.rm=T)) %>%
                           mutate(Rationale="CNI") %>%
                           select(Gene, Rationale) %>%
                           add_row(OKBdata %>%
                                       filter(Alterations=="Amplification") %>%
                                       select(Gene) %>%
                                       distinct() %>%
                                       mutate(Rationale="OncoKB")) %>%
                           mutate(Rationale=ifelse(duplicated(Gene, fromLast=T), "CNI|OncoKB", Rationale)) %>%
                           filter(!duplicated(Gene))

ampdata             <- ampdata %>%
                           mutate(Pathogenic = Gene %in% AMP_predicted_genes$Gene)

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# EXPRESSION QUANTILES FOR AMPLIFICATIONS

quantlist           <- lapply(AMP_predicted_genes$Gene, function(mygene) ecdf(samplegex[mygene,]))
names(quantlist)    <- AMP_predicted_genes$Gene

genesamplequant     <- do.call(rbind, 
                               lapply(1:length(quantlist), function(i) {
                                   temp <- as.data.frame(matrix(quantlist[[i]](samplegex[names(quantlist)[i],]), nrow=1))
                                   row.names(temp) <- names(quantlist)[i]
                                   colnames(temp)  <- colnames(samplegex)
                                   return(temp)
                               }))


RNAtemp             <- do.call(rbind, lapply(which(ampdata$Pathogenic), function(i) {
                                   rnarefs <- which(gexanno$sample == ampdata$patient[i] & 
                                                    gexanno$stime  == ampdata$stime[i] & 
                                                    gexanno$ssite  == ampdata$ssite[i] & 
                                                   (gexanno$sside  == ampdata$sside[i] | is.na(ampdata$sside[i])) & 
                                                   (gexanno$sord   == ampdata$sord[i] | is.na(ampdata$sord[i])))
                                   if(length(rnarefs)==0) {
                                   rnarefs <- which(gexanno$sample == ampdata$patient[i] & 
                                                    gexanno$stime  == ampdata$stime[i] & 
                                                    gexanno$ssite  == ampdata$ssite[i] & 
                                                    is.na(gexanno$sside[i]) & 
                                                    is.na(gexanno$sord[i]))
                                   }                                                                
                                   return(data.frame(DNA_sample    = ampdata$sample[i],
                                                     Gene          = ampdata$Gene[i],
                                                     RNA_sample    = ifelse(length(rnarefs)==0, NA, 
                                                                     ifelse(length(rnarefs)==1, row.names(gexanno)[rnarefs], 
                                                                     names(which.max(genesamplequant[ampdata$Gene[i],row.names(gexanno)[rnarefs]])))),
                                                     RNA_quant     = ifelse(length(rnarefs)==0, NA,
                                                                     ifelse(length(rnarefs)==1, genesamplequant[ampdata$Gene[i], row.names(gexanno)[rnarefs]],
                                                                     max(genesamplequant[ampdata$Gene[i], row.names(gexanno)[rnarefs]])))))
                                                                }))

cnadata             <- ampdata %>%
                           mutate(mkey=paste(sample, Gene, sep="_")) %>%
                           merge(RNAtemp %>% mutate(mkey=paste(DNA_sample, Gene, sep="_")) %>% select(mkey, RNA_sample, Expression_quantile=RNA_quant),
                                 by="mkey", all.x=T) %>%
                           select(-mkey) %>%
                           group_by(sample) %>% 
                           mutate(RNA.sample=ifelse(is.na(RNA_sample), unique(c(na.omit(RNA_sample)))[1], RNA_sample)) %>%
                           ungroup() %>%
                           select(-RNA_sample)

                                     
write.table(cnadata, file="amplifications.csv", col.names=T, row.names=F, quote=F, sep="\t")














