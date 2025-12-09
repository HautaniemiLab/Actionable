#!/usr/bin/env Rscript

library("tidyverse")

args = commandArgs(trailingOnly=TRUE)

if (length(args)!=2) {
  stop("Provide two arguments \n Somatic short mutation file \n CGI output file \n", call.=FALSE)
} 

mutations        <- read.delim(args[1], header=T, as.is=T, sep="\t")
CGIdata          <- read.delim(args[2], header=T, as.is=T, sep="\t", na.strings=c("NA", "", " "))

# -----------------------------------------------------------------------------
# Clean and merge predictions and annotations for pathogenicity
# -----------------------------------------------------------------------------

mutations_annotated <- mutations %>%
                                mutate(PolyPhenCat=str_replace(str_replace(PolyPhenCat, ",NA", ""), "NA,", ""),
                                       PolyPhenVal=str_replace(str_replace(PolyPhenVal, ",NA", ""), "NA,", ""),
                                       SIFTcat=str_replace(str_replace(SIFTcat, ",NA", ""), "NA,", ""),
                                       SIFTval=str_replace(str_replace(SIFTval, ",NA", ""), "NA,", ""),
                                       AM_class=ifelse(str_detect(AM_class, "ambiguous"), "ambiguous", AM_class)) %>%
                                mutate(mergekey=paste(CHROM, POS, REF, ALT, sep="_")) %>%
                                merge(CGIdata %>% mutate(mergekey=paste(paste("chr", CHROMOSOME, sep=""), POSITION, REF, ALT, sep="_")) %>% 
                                                  select(CGI.Protein.Change, CGI.Oncogenic.Summary, CGI.Consequence, mergekey),
                                      by="mergekey", all.x=T) %>%
                               filter(!((Type=="intronic" & (is.na(CGI.Consequence) | CGI.Consequence=="intron_variant")) | 
                                        (Type %in% c("synonymous SNV", "upstream"))))

write.table(mutations_annotated, file="somatic.expressed.annotated.csv", col.names=T, row.names=F, quote=F, sep="\t")

