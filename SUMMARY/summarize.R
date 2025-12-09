#!/usr/bin/env Rscript

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Summarize data from directories DRUGS, SNV, SV, DEL, and AMP 
library("tidyverse")

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
# Read in matched DNA-RNA sample file with columns. The file includes only valid samples from eligible patients.
# patient            - patient ID
# DNA.sample         - DNA sample ID
# RNA.sample         - RNA sample ID
# stime              - sample time (e.g. 'p' for PDS, 'i' for IDS, 'r' for relapse, "p2" and "o" for operations between PDS and IDS)
# ssite              - tumor site (e.g. "Ova" for ovaries, "Ome" for omentum, "Per" for peritoneum)
# sside              - Right (R) or Left (L) if applicable
# sord               - Order of the sample, if multiple biopsies were taken from the same anatomical site at the same operation
cancersamples        <- read.delim("<appropriate_path/cancersamples.csv", header=T, as.is=T, sep="\t")

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Read in short somatic mutations
snvdata              <- read.delim("<appropriate_path>/somatic.expressed.annotated.csv", header=T, as.is=T, sep="\t", na.strings=c("NA", "", ".")) %>%
                            filter(sample %in% cancersamples$DNA.sample)

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Read in structural mutations
brkdata              <- read.delim("<appropriate_path>/breakevents_expressed.csv", header=T, as.is=T, sep="\t") %>%
                            filter(sample %in% cancersamples$DNA.sample)

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Read in gene-level deletions
deldata              <- read.delim("<appropriate_path>/deletions.csv", header=T, as.is=T, sep="\t") %>%
                            filter(sample %in% cancersamples$DNA.sample)

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Read in gene amplifications
ampdata              <- read.delim("<appropriate_path>/amplifications.csv", header=T, as.is=T, sep="\t") %>%
                            filter(sample %in% cancersamples$DNA.sample)

# --------------------------------------------------------------------------------------------------------------------------------------------------
# PROCESS DATA
# Modify column names, so that the output from different analysis pipelines provide similar information
# If allele-specific expression was not detected, but an informative sample was available, RNA-reads is imputed 0

deltemp      <- deldata %>%
                    mutate(patient = str_extract(sample, "^[:alnum:]+"),
                           temp    = str_extract(sample, "(?<=_)[poir][:alnum:]+")) %>%
                    mutate(stime   = str_extract(temp, "^[poir][:digit:]?"),
                           ssite   = str_extract(temp, "(?<=[poir][:digit:]?)([:upper:][:lower:]{2}|LN)"),
                           sside   = str_extract(temp, "(?<=([:upper:][:lower:]{2}|LN))[LR]"),
                           sord    = str_extract(temp, "[:digit:]$")) %>%
                    select(-temp) %>%
                    mutate(Gene.type = ifelse(Gene %in% lossgenes, "loss", ifelse(Gene %in% gaingenes, "gain", NA)))

snvtemp      <- snvdata %>%
                    mutate(mkey        = paste(patient, stime, ssite, sside, sord, sep="_")) %>%
                    merge(cancersamples %>% mutate(mkey= paste(patient, stime, ssite, sside, sord, sep="_")) %>% select(mkey, RNA.sample),
                          by="mkey", all.x=TRUE) %>%
                    mutate(RNAname     = ifelse(is.na(RNAname), RNA.sample, RNAname)) %>%
                    select(-mkey, -RNA.sample) %>%
                    mutate(refCount    = ifelse(is.na(refCount) & !is.na(RNAname), 0, refCount),
                           altCount    = ifelse(is.na(altCount) & !is.na(RNAname), 0, altCount)) %>%
                    mutate(Gene.type   = ifelse(Gene %in% gaingenes, "gain", ifelse(Gene %in% lossgenes, "loss", NA))) %>%
                    mutate(ExonicFunc  = ifelse((!is.na(dbscSNV_ADA_SCORE) & dbscSNV_ADA_SCORE>0.95 & !is.na(dbscSNV_RF_SCORE) & dbscSNV_RF_SCORE>0.95) | str_detect(Func.MANE, "splic"), "splicing", ExonicFunc.MANE),
                           SIFTnum     = ifelse(!is.na(SIFTcat) & str_detect(SIFTcat, "deleterious"), 1, 0),
                           PP2num      = ifelse(!is.na(PolyPhenCat) & str_detect(PolyPhenCat, "probably"), 1, 0),
                           AMnum       = ifelse(!is.na(AM_class) & str_detect(AM_class, "pathogenic"), 1, 0)) %>%
                    mutate(Deleterious = (SIFTnum + PP2num + AMnum) > 1) %>%
                    filter(!is.na(ExonicFunc) & ExonicFunc!="synonymous_SNV" & ExonicFunc!="startloss") %>%
                    mutate(Pathogenic  = ifelse(!is.na(ClinVar) & ClinVar %in% c("Benign", "Benign/Likely_benign", "Likely_benign"), FALSE,
                                         ifelse(!is.na(ClinVar) & ClinVar %in% c("Pathogenic", "Pathogenic/Likely_pathogenic", "Likely_pathogenic"), TRUE,
                                         ifelse((ExonicFunc=="nonsynonymous_SNV" | ExonicFunc=="nonframeshift_substitution") & is.na(CGI.Oncogenic.Summary) & is.na(Deleterious), FALSE,
                                         ifelse((ExonicFunc=="nonsynonymous_SNV" | ExonicFunc=="nonframeshift_substitution") & !is.na(CGI.Oncogenic.Summary) & str_detect(CGI.Oncogenic.Summary, "^onco"), TRUE,
                                         ifelse((ExonicFunc=="nonsynonymous_SNV" | ExonicFunc=="nonframeshift_substitution") & !is.na(Deleterious) & Deleterious, TRUE,
                                         ifelse(str_detect(ExonicFunc, "nonframeshift"), TRUE,
                                         ifelse(str_detect(ExonicFunc, "frameshift|stop|splic") & Gene %in% c("KMT2A", "EZH2"), TRUE,
                                         ifelse(Gene.type=="gain", FALSE,
                                         ifelse(Gene.type=="loss", str_detect(ExonicFunc, "^frameshift") | ExonicFunc %in% c("splicing", "stopgain"), NA)))))))))) %>%
                    mutate(Expressed   = altCount>=5,
                           Homogeneous = expHom.pbinom.lower>0.01 & AD.1>0)

brktemp      <- brkdata %>%
                    mutate(Gene.type   = ifelse(primary_gene %in% gaingenes, "gain", ifelse(primary_gene %in% lossgenes, "loss", NA)),
                           Homogeneous = undisruptedCopyNumber<0.5)


# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Summation to patient level, separately for oncogenes (Gene.type=="gain") and tumor-suppressor genes (Gene.type=="loss")
# The following data objects were used for producing the manuscript Figure 2.

snv_gain_per_patient    <- snvtemp %>%
                               filter(Gene.type=="gain") %>%
                               mutate(gkey=paste(patient, Mutation, sep="_")) %>%
                               group_by(gkey) %>%
                               summarize(patient=unique(patient),
                                         Mutation=unique(Mutation),
                                         Gene=unique(Gene),
                                         HGVS.change=unique(HGVS.change),
                                         ExonicFunc=unique(ExonicFunc),
                                         CredOme=any(AD.1>0 & stime %in% c("p", "p2", "o", "i") & ssite=="Ome" & Pathogenic & !is.na(Expressed) & Expressed, na.rm=TRUE),
                                         CredPer=any(AD.1>0 & stime %in% c("p", "p2", "o", "i") & (ssite %in% c("Per", "Mes", "Bow", "Oth", "Ute", "Vag")) & Pathogenic & !is.na(Expressed) & Expressed, na.rm=TRUE),
                                         CredAdn=any(AD.1>0 & stime %in% c("p", "p2", "o", "i") & (ssite %in% c("Ova", "Tub", "Adn")) & Pathogenic & !is.na(Expressed) & Expressed, na.rm=TRUE),
                                         Relapse=ifelse(!any(str_detect(stime, "r")), NA, any(AD.1[which(str_detect(stime, "r"))]>0)),
                                         Allsamples=length(which(AD.1>0))==n(),
                                         Pathogenic=any(Pathogenic),
                                         Expressed=ifelse(all(is.na(Expressed)), NA, any(Expressed, na.rm=TRUE)))
snv_loss_per_patient    <- snvtemp %>%
                               filter(Gene.type=="loss") %>%
                               mutate(gkey=paste(patient, Mutation, sep="_")) %>%
                               group_by(gkey) %>%
                               summarize(patient=unique(patient),
                                         Mutation=unique(Mutation),
                                         Gene=unique(Gene),
                                         HGVS.change=unique(HGVS.change),
                                         ExonicFunc=unique(ExonicFunc),
                                         CredOme=any(AD.1>0 & stime %in% c("p", "p2", "o", "i") & ssite=="Ome" & Pathogenic & !is.na(Homogeneous) & Homogeneous, na.rm=TRUE),
                                         CredPer=any(AD.1>0 & stime %in% c("p", "p2", "o", "i") & (ssite %in% c("Per", "Mes", "Bow", "Oth", "Ute", "Vag")) & Pathogenic & !is.na(Homogeneous) & Homogeneous, na.rm=TRUE),
                                         CredAdn=any(AD.1>0 & stime %in% c("p", "p2", "o", "i") & (ssite %in% c("Ova", "Tub", "Adn")) & Pathogenic & !is.na(Homogeneous) & Homogeneous, na.rm=TRUE),
                                         Relapse=ifelse(!any(str_detect(stime, "r")), NA, any(AD.1[which(str_detect(stime, "r"))]>0)),
                                         Allsamples=length(which(AD.1>0))==n(),
                                         Pathogenic=any(Pathogenic),
                                         Homogeneous=any(Homogeneous, na.rm=TRUE))

brk_gain_per_patient    <- brktemp %>%
                               filter(Gene.type=="gain") %>%
                               mutate(gkey=paste(patient, primary_gene, sep="_")) %>%
                               mutate(undisruptedCopyNumber=ifelse(is.na(undisruptedCopyNumber), NA, ifelse(undisruptedCopyNumber<0, 0, undisruptedCopyNumber)),
                                      affectedCopyNumber=ifelse(is.na(affectedCopyNumber), NA, ifelse(affectedCopyNumber<0, 0, affectedCopyNumber))) %>%
                               group_by(gkey) %>%
                               summarize(patient=unique(patient),
                                         Gene=unique(primary_gene),
                                         DNA.consequence=ifelse(length(unique(effect))>1, "multiple", unique(effect)),
                                         Mutation=ifelse(any(expressed, na.rm=T), paste(unique(mutation[which(expressed)], "|")), paste(unique(mutation))),
                                         DNA.change=ifelse(any(expressed, na.rm=T), paste(unique(DNA_change[which(expressed)], "|")), paste(unique(DNA_change))),
                                         RNA.change=ifelse(any(expressed, na.rm=T), paste(unique(RNA.change[which(expressed)], "|")), paste(unique(RNA.change))),
                                         RNA.consequence=ifelse(any(expressed, na.rm=T), paste(unique(RNA.consequence[which(expressed)], "|")), paste(unique(RNA.consequence))),
                                         cn_normal_range=ifelse(all(is.na(undisruptedCopyNumber)), NA, paste(unique(round(range(undisruptedCopyNumber, na.rm=T),0)), collapse="-")),
                                         cn_affected_range=ifelse(all(is.na(affectedCopyNumber)), NA, paste(unique(round(range(affectedCopyNumber, na.rm=T),0)), collapse="-")),
                                         CredOme=any(stime %in% c("p", "p2", "o", "i") & ssite=="Ome" & pathogenic & !is.na(expressed) & expressed, na.rm=TRUE),
                                         CredPer=any(stime %in% c("p", "p2", "o", "i") & (ssite %in% c("Per", "Mes", "Bow", "Oth", "Ute", "Vag")) & pathogenic & !is.na(expressed) & expressed, na.rm=TRUE),
                                         CredAdn=any(stime %in% c("p", "p2", "o", "i") & (ssite %in% c("Ova", "Tub", "Adn")) & pathogenic & !is.na(expressed) & expressed, na.rm=TRUE),
                                         Relapse=ifelse(!any(str_detect(cancersamples$stime[which(cancersamples$patient==unique(patient))], "r")), NA, any(str_detect(stime, "r"))),
                                         Allsamples=length(which(cancersamples$patient==unique(patient)))==n(), 
                                         Pathogenic=any(pathogenic),
                                         Expressed=any(expressed, na.rm=TRUE))

brk_loss_per_patient    <- brktemp %>%
                               filter(Gene.type=="loss") %>%
                               mutate(gkey=paste(patient, primary_gene, sep="_")) %>%
                               mutate(undisruptedCopyNumber=ifelse(is.na(undisruptedCopyNumber), NA, ifelse(undisruptedCopyNumber<0, 0, undisruptedCopyNumber)),
                                      affectedCopyNumber=ifelse(is.na(affectedCopyNumber), NA, ifelse(affectedCopyNumber<0, 0, affectedCopyNumber))) %>%
                               group_by(gkey) %>%
                               summarize(patient=unique(patient),
                                         Gene=unique(primary_gene),
                                         DNA.consequence=ifelse(length(unique(effect))>1, "multiple", unique(effect)),
                                         Mutation=ifelse(any(Homogeneous, na.rm=T), paste(unique(mutation[which(Homogeneous)], "|")), paste(unique(mutation))),
                                         DNA.change=ifelse(any(Homogeneous, na.rm=T), paste(unique(DNA_change[which(Homogeneous)], "|")), paste(unique(DNA_change))),
                                         RNA.change=ifelse(any(Homogeneous, na.rm=T), paste(unique(RNA.change[which(Homogeneous)], "|")), paste(unique(RNA.change))),
                                         RNA.consequence=ifelse(any(Homogeneous, na.rm=T), paste(unique(RNA.consequence[which(Homogeneous)], "|")), paste(unique(RNA.consequence))),
                                         cn_normal_range=ifelse(all(is.na(undisruptedCopyNumber)), NA, paste(unique(round(range(undisruptedCopyNumber, na.rm=T),0)), collapse="-")),
                                         cn_affected_range=ifelse(all(is.na(affectedCopyNumber)), NA, paste(unique(round(range(affectedCopyNumber, na.rm=T),0)), collapse="-")),
                                         CredOme=any(stime %in% c("p", "p2", "o", "i") & ssite=="Ome" & pathogenic & !is.na(Homogeneous) & Homogeneous, na.rm=TRUE),
                                         CredPer=any(stime %in% c("p", "p2", "o", "i") & (ssite %in% c("Per", "Mes", "Bow", "Oth", "Ute", "Vag")) & pathogenic & !is.na(Homogeneous) & Homogeneous, na.rm=TRUE),
                                         CredAdn=any(stime %in% c("p", "p2", "o", "i") & (ssite %in% c("Ova", "Tub", "Adn")) & pathogenic & !is.na(Homogeneous) & Homogeneous, na.rm=TRUE),
                                         Relapse=ifelse(!any(str_detect(cancersamples$stime[which(cancersamples$patient==unique(patient))], "r")), NA, any(str_detect(stime, "r"))),
                                         Allsamples=length(which(cancersamples$patient==unique(patient)))==n(), 
                                         Pathogenic=any(pathogenic),
                                         Homogeneous=any(Homogeneous, na.rm=T))

brk_loss_per_patient    <- brk_loss_per_patient %>%
                               add_row(deltemp  %>%
                                   filter(Gene.type=="loss") %>%
                                   mutate(gkey=paste(patient, Gene, sep="_")) %>%
                                   merge(brk_loss_per_patient %>% select(gkey, cn_normal_range), by="gkey", all.x=T) %>%
                                   filter(is.na(cn_normal_range)) %>%
                                   group_by(gkey) %>%
                                   summarize(patient=unique(patient),
                                             Gene=unique(Gene),
                                             DNA.consequence="deletion",
                                             Mutation=NA,
                                             DNA.change=paste(unique(Gene), "del", sep=":"),
                                             RNA.change=NA,
                                             RNA.consequence=NA,
                                             cn_normal_range="0",
                                             cn_affected_range=NA,
                                             CredOme=any(stime %in% c("p", "p2", "o", "i") & ssite=="Ome", na.rm=TRUE),
                                             CredPer=any(stime %in% c("p", "p2", "o", "i") & (ssite %in% c("Per", "Mes", "Bow", "Oth", "Ute", "Vag")), na.rm=TRUE),
                                             CredAdn=any(stime %in% c("p", "p2", "o", "i") & (ssite %in% c("Ova", "Tub", "Adn")), na.rm=TRUE),
                                             Relapse=ifelse(!any(str_detect(cancersamples$stime[which(cancersamples$patient==unique(patient))], "r")), NA, any(str_detect(stime, "r"))),
                                             Allsamples=length(which(cancersamples$patient==unique(patient)))==n(), 
                                             Pathogenic=TRUE,
                                             Homogeneous=TRUE))

amp_gain_per_patient     <- ampdata %>%
                                mutate(gkey=paste(patient, Gene, sep="_")) %>%
                                group_by(gkey) %>%
                                summarize(patient=unique(patient),
                                             Gene=unique(Gene),
                                             consequence="amplification",
                                             region=NA,
                                             cn_normal_range=paste(unique(round(range(nMinor, na.rm=T),0)), collapse="-"),
                                             cn_affected_range=paste(unique(round(range(nMajor, na.rm=T),0)), collapse="-"),
                                             CredOme=any(AMP_ploidy & stime %in% c("p", "p2", "o", "i") & ssite=="Ome" & Expression_quantile>=0.8, na.rm=TRUE),
                                             CredPer=any(AMP_ploidy & stime %in% c("p", "p2", "o", "i") & (ssite %in% c("Per", "Mes", "Bow", "Oth", "Ute", "Vag")) & Expression_quantile>=0.8, na.rm=TRUE),
                                             CredAdn=any(AMP_ploidy & stime %in% c("p", "p2", "o", "i") & (ssite %in% c("Ova", "Tub", "Adn")) & Expression_quantile>=0.8, na.rm=TRUE),
                                             Relapse = ifelse(!any(str_detect(stime, "r")), NA, any(AMP_ploidy & str_detect(stime, "r"))),
                                             Allsamples=length(which(AMP_ploidy))==n(),
                                             Pathogenic=any(Pathogenic),
                                             Expressed=any(Expression_quantile>=0.8))

# --------------------------------------------------------------------------------------------------------------- #
# Manually curated somatic mutation data for genes TP53, BRCA1, BRCA2, RAD51C, RAD51D, CDK12, and NF1,            #
# as well as germline risk variant data was added to their respective aberration data-frames at this point        #
# Overlaps between SNV and SV data were identified and resolved, so that every aberration was counted only once.  #                                                  
# --------------------------------------------------------------------------------------------------------------- #

# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Combining data layers, assigning ESCAT Tiers

aberrations_per_patient      <- amp_gain_per_patient %>%
                                    mutate(Patient=patient,
                                           Gene.type="gain",
                                           Aberration.name=paste(Gene, "amp", sep=":"),
                                           Genomic.position=NA,
                                           Consequence=consequence,
                                           CN.normal.range=cn_normal_range,
                                           CN.affected.range=cn_affected_range,
                                           Homogeneous=NA) %>%
                                    select(gkey, Patient, Gene, Gene.type, Aberration.name, Genomic.position, Consequence, CN.normal.range, CN.affected.range,
                                           CredOme, CredPer, CredAdn, Relapse, Allsamples, 
                                           Pathogenic, Expressed, Homogeneous) %>%
                        add_row(brk_gain_per_patient %>%
                                    mutate(Patient=patient,
                                           Gene.type="gain",
                                           Aberration.name=DNA.change,
                                           Genomic.position=Mutation,
                                           Consequence=DNA.consequence,
                                           CN.normal.range=cn_normal_range,
                                           CN.affected.range=cn_affected_range,
                                           Homogeneous=NA)  %>%
                                    select(gkey, Patient, Gene, Gene.type, Aberration.name, Genomic.position, Consequence, CN.normal.range, CN.affected.range,
                                           CredOme, CredPer, CredAdn, Relapse, Allsamples, 
                                           Pathogenic, Expressed, Homogeneous)) %>%
                        add_row(snv_gain_per_patient %>%
                                    mutate(gkey=paste(patient, Gene, sep="_"),
                                           Patient=patient,
                                           Gene.type="gain",
                                           Aberration.name=HGVS.change,
                                           Genomic.position=Mutation,
                                           Consequence=str_replace(ExonicFunc, "nonframeshift", "inframe"),
                                           CN.normal.range=NA,
                                           CN.affected.range=NA,
                                           Homogeneous=NA)  %>%
                                    select(gkey, Patient, Gene, Gene.type, Aberration.name, Genomic.position, Consequence, CN.normal.range, CN.affected.range,
                                           CredOme, CredPer, CredAdn, Relapse, Allsamples, 
                                           Pathogenic, Expressed, Homogeneous)) %>%
                        add_row(brk_loss_per_patient %>%
                                    mutate(Patient=patient,
                                           Gene.type="loss",
                                           Aberration.name=DNA.change,
                                           Genomic.position=Mutation,
                                           Consequence=DNA.consequence,
                                           CN.normal.range=cn_normal_range,
                                           CN.affected.range=cn_affected_range,
                                           Expressed=NA)  %>%
                                    select(gkey, Patient, Gene, Gene.type, Aberration.name, Genomic.position, Consequence, CN.normal.range, CN.affected.range,
                                           CredOme, CredPer, CredAdn, Relapse, Allsamples, 
                                           Pathogenic, Expressed, Homogeneous)) %>%
                        add_row(snv_loss_per_patient %>%
                                    mutate(Patient=patient,
                                           Gene.type="loss",
                                           Aberration.name=HGVS.change,
                                           Genomic.position=Mutation,
                                           Consequence=str_replace(ExonicFunc, "nonframeshift", "inframe"),
                                           CN.normal.range=NA,
                                           CN.affected.range=NA,
                                           Expressed=NA)  %>%
                                    select(gkey, Patient, Gene, Gene.type, Aberration.name, Genomic.position, Consequence, CN.normal.range, CN.affected.range,
                                           CredOme, CredPer, CredAdn, Relapse, Allsamples, 
                                           Pathogenic, Expressed, Homogeneous))



# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Read in OncoKB v4.1 data with the following columns.
# Level                                       - FDA therapeutic level for the Gene-treatment combination
# Gene                                        - HUGO gene symbol
# Alterations                                 - Gene mutations, specific mutations or mutation categories
# Cancer.Types                                - Cancer histologies
# Drugs..for.therapeutic.implications.only.   - Drug combinations

OKBdata              <- read.delim("OKBdata_v4.1.csv", header=T, as.is=T, sep="\t", na.strings=c("", "NA")) 

# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ESCAT Tiers were derived from OncoKB FDA therapeutic levels
# Tier IA - FDA 1 or 2 in ovarian cancer
# Tier IC - FDA 1 or 2 in solid tumors

mytierI          <-  Onco_drugmut %>% filter(str_detect(Cancer.Types, "[Oo]va")   & Level %in% c("1", "2")) %>% 
                                      select(-Cancer.Types, -Level) %>% 
                                      distinct() %>%
                                      mutate(Tier="IA",
                                             type=ifelse(Alterations=="Amplification", "amp", 
                                                      ifelse(Alterations=="Fusions", "fus",
                                                          ifelse(str_detect(Alterations, "^[:alpha:]*\\s?Mutations$"), "mut",
                                                              ifelse(str_detect(Alterations, "^[:upper:][:digit:]+[:alpha:]+$"), "spec",
                                                                  ifelse(str_detect(Alterations, "^[:upper:][:digit:]+$"), "pos", "man"))))),
                                             reserved=paste(Gene, Alterations)) %>%
             add_row(Onco_drugmut %>% filter(str_detect(Cancer.Types, "[Ss]olid") & Level %in% c("1", "2") & !str_detect(Gene, "\\s")) %>% 
                                      select(-Cancer.Types, -Level) %>% 
                                      distinct() %>%
                                      mutate(Tier="IC",
                                             type=ifelse(Alterations=="Amplification", "amp", 
                                                      ifelse(Alterations=="Fusions", "fus",
                                                          ifelse(str_detect(Alterations, "^[:alpha:]*\\s?Mutations$"), "mut",
                                                              ifelse(str_detect(Alterations, "^[:upper:][:digit:]+[:alpha:]+$"), "spec",
                                                                  ifelse(str_detect(Alterations, "^[:upper:][:digit:]+$"), "pos", "man"))))),
                                             reserved=paste(Gene, Alterations))) %>%
                                      group_by(reserved) %>%
                                      summarize(Gene=unique(Gene),
                                                Alterations=unique(Alterations),
                                                Drugs=paste(Drugs..for.therapeutic.implications.only., collapse=", "),
                                                Tier=Tier[1],
                                                type=unique(type))


# Tier II - FDA 3 in ovca or solid tumors - some clinical evidence in the same cancer type & must be an approved drug

mytierII         <-  Onco_drugmut %>% 
                                      group_by(Drugs..for.therapeutic.implications.only.) %>%
                                      mutate(approved = Drugs..for.therapeutic.implications.only. %in% c("Zenocutuzumab", "Selitrectinib") | any(Level %in% c("1", "2"))) %>%
                                      ungroup() %>%
                                      filter(str_detect(Cancer.Types, "([Oo]va)|([Ss]olid)")   & Level=="3" & approved) %>% 
                                      select(-Cancer.Types, -Level) %>% 
                                      distinct() %>% 
                                      mutate(Tier="II",
                                              type=ifelse(Alterations=="Amplification", "amp", 
                                                       ifelse(Alterations=="Fusions", "fus",
                                                           ifelse(str_detect(Alterations, "^[:alpha:]*\\s?Mutations$"), "mut",
                                                               ifelse(str_detect(Alterations, "^[:upper:][:digit:]+[:alpha:]+$"), "spec",
                                                                   ifelse(str_detect(Alterations, "^[:upper:][:digit:]+$"), "pos", "man"))))),
                                              reserved=paste(Gene, Alterations)) %>%
                                      filter(!(reserved %in% mytierI$reserved)) %>%
                                      group_by(reserved) %>%
                                      summarize(Gene=unique(Gene),
                                                Alterations=unique(Alterations),
                                                Drugs=paste(Drugs..for.therapeutic.implications.only., collapse=", "),
                                                Tier=Tier[1],
                                                type=unique(type)) %>%
                                      filter(Drugs != "Trastuzumab Deruxtecan") # ADC drugs not included


# Tier IIIA - FDA 1 or 2 in any other cancer

mytierIIIA       <-  Onco_drugmut %>% mutate(Alterations= str_replace(Alterations, "\\s\\(.+\\)", "")) %>%
                                      group_by(Drugs..for.therapeutic.implications.only.) %>%
                                      mutate(approved = Drugs..for.therapeutic.implications.only. %in% c("Zenocutuzumab", "Selitrectinib") | any(Level %in% c("1", "2"))) %>%
                                      ungroup() %>%
                                      filter(Level %in% c("1", "2", "3") & !str_detect(Gene, "\\s") & approved)  %>% 
                                      select(-Cancer.Types, -Level) %>%
                                      mutate(Alterations = strsplit(Alterations, split=", ")) %>% 
                                      unnest(Alterations) %>% 
                                      distinct() %>%
                                      mutate(Tier="IIIA",
                                             type=ifelse(Alterations=="Amplification", "amp", 
                                                      ifelse(Alterations=="Fusions", "fus",
                                                          ifelse(str_detect(Alterations, "^[:alpha:]*\\s?Mutations$"), "mut",
                                                              ifelse(str_detect(Alterations, "^[:upper:][:digit:]+[:alpha:]+$"), "spec",
                                                                  ifelse(str_detect(Alterations, "^[:upper:][:digit:]+$"), "pos", "man"))))),
                                             reserved=paste(Gene, Alterations)) %>%
                                      filter(!(reserved %in% mytierI$reserved) & !(reserved %in% mytierII$reserved)) %>%
                                      group_by(reserved) %>%
                                      summarize(Gene=unique(Gene),
                                                Alterations=unique(Alterations),
                                                Drugs=paste(Drugs..for.therapeutic.implications.only., collapse=", "),
                                                Tier=Tier[1],
                                                type=unique(type))

# Tier IIIB
# Any pathogenic aberration (but not amplification) in the Tier I-IIIA genes, but without match to annotated aberrations

mytierIIIB      <- mytierI %>%
                   add_row(mytierII) %>%
                   add_row(mytierIIIA) %>%
                   group_by(Gene) %>%
                   mutate(mutstat = any(str_detect(Alterations, "^[:alpha:]*\\s?Mutations$")),
                          fusstat = any(Alterations=="Fusions")) %>%
                   select(Gene, Drugs, mutstat, fusstat) %>%
                   filter(!mutstat | !fusstat) %>%
                   filter(!(Gene %in% lossgenes)) %>%
                   distinct() %>%
                   summarize(Gene=unique(Gene),
                             Drugs=paste(Drugs, collapse=", "),
                             mutstat=any(mutstat),
                             fusstat=any(fusstat)) %>%
                   mutate(Drugs=unlist(lapply(Drugs, function(dr) paste(unique(unlist(strsplit(dr, split=", "))), collapse=", ")))) %>%
                   mutate(Alterations=ifelse(!mutstat & !fusstat, "Oncogenic Mutations, Fusions",
                                          ifelse(!mutstat, "Oncogenic Mutations",
                                              ifelse(!fusstat, "Fusions", NA)))) %>%
                   mutate(Alterations=strsplit(Alterations, split=", ")) %>%
                   unnest(Alterations) %>%
                   mutate(type=ifelse(Alterations=="Fusions", "fus",
                                   ifelse(str_detect(Alterations, "^[:alpha:]*\\s?Mutations$"), "mut", NA)),
                          reserved=paste(Gene, Alterations),
                          Tier="IIIB") %>%
                   select(reserved, Gene, Alterations, Drugs, Tier, type) %>%
                   filter(!(reserved %in% mytierI$reserved) & !(reserved %in% mytierII$reserved) & !(reserved %in% mytierIIIA$reserved))

# Tier IVA - FDA 3 in ovca or solid tumors

mytierIVA1         <-  Onco_drugmut %>% filter(str_detect(Cancer.Types, "([Oo]va)|([Ss]olid)")   & Level=="3") %>% 
                                      select(-Cancer.Types, -Level) %>% 
                                      distinct() %>% 
                                      mutate(Tier="IVA",
                                              type=ifelse(Alterations=="Amplification", "amp", 
                                                       ifelse(Alterations=="Fusions", "fus",
                                                           ifelse(str_detect(Alterations, "^[:alpha:]*\\s?Mutations$"), "mut",
                                                               ifelse(str_detect(Alterations, "^[:upper:][:digit:]+[:alpha:]+$"), "spec",
                                                                   ifelse(str_detect(Alterations, "^[:upper:][:digit:]+$"), "pos", "man"))))),
                                              reserved=paste(Gene, Alterations)) %>%
                                      filter(!(reserved %in% mytierI$reserved) & !(reserved %in% mytierII$reserved) & !(reserved %in% mytierIIIA$reserved) & !(reserved %in% mytierIIIB$reserved)) %>%
                                      group_by(reserved) %>%
                                      summarize(Gene=unique(Gene),
                                                Alterations=unique(Alterations),
                                                Drugs=paste(Drugs..for.therapeutic.implications.only., collapse=", "),
                                                Tier=Tier[1],
                                                type=unique(type)) %>%
                                      filter(Drugs != "Trastuzumab Deruxtecan")


# Tier IVA
# FDA level 4 matching the aberrations

mytierIVA2        <-  Onco_drugmut %>% mutate(Alterations=str_replace(Alterations, "\\s\\(.+\\)", "")) %>%
                                      filter(Level %in% as.character(1:4) & !str_detect(Gene, "\\s")) %>% 
                                      mutate(Level=as.numeric(Level)) %>% 
                                      group_by(Gene) %>% 
                                      mutate(Level=min(Level)) %>% 
                                      filter(Level==4) %>% 
                                      select(-Cancer.Types, -Level) %>%
                                      distinct() %>%
                                      mutate(Tier="IVA",
                                             type=ifelse(Alterations=="Amplification", "amp", 
                                                      ifelse(Alterations=="Fusions", "fus",
                                                          ifelse(str_detect(Alterations, "^[:alpha:]*\\s?Mutations$"), "mut",
                                                              ifelse(str_detect(Alterations, "^[:upper:][:digit:]+[:alpha:]+$"), "spec",
                                                                  ifelse(str_detect(Alterations, "^[:upper:][:digit:]+$"), "pos", "man"))))),
                                             reserved=paste(Gene, Alterations)) %>%
                                      filter(!(reserved %in% mytierI$reserved) & !(reserved %in% mytierIIIA$reserved) & !(reserved %in% mytierIIIB$reserved)) %>%
                                      group_by(reserved) %>%
                                      summarize(Gene=unique(Gene),
                                                Alterations=unique(Alterations),
                                                Drugs=paste(Drugs..for.therapeutic.implications.only., collapse=", "),
                                                Tier=Tier[1],
                                                type=unique(type))

# Tier IVB
# Amplifications of the oncogenes go here
# Aberrations in secondary drug targets of the Tier I-IIIA drugs

mytierIVB        <- network %>%
                        filter(proteinTarget %in% gaingenes & !is.na(okb.drug)) %>%
                        select(Gene=proteinTarget, Drugs=okb.drug) %>%
                        group_by(Gene) %>%
                        summarize(Drugs=paste(unique(Drugs), collapse=", ")) %>%
                        mutate(Alterations=strsplit("Amplification, Fusions, Oncogenic Mutations", split=", ")) %>%
                        unnest(Alterations) %>%
                        mutate(reserved=paste(Gene, Alterations)) %>%
                        filter(!(reserved %in% mytierI$reserved) & !(reserved %in% mytierIIIA$reserved) & !(reserved %in% mytierIIIB$reserved) & !(reserved %in% mytierIVA1$reserved) & !(reserved %in% mytierIVA2$reserved)) %>%
                        mutate(Tier="IVB",
                               type=ifelse(Alterations=="Amplification", "amp", 
                                        ifelse(Alterations=="Fusions", "fus",
                                            ifelse(str_detect(Alterations, "^[:alpha:]*\\s?Mutations$"), "mut", NA)))) %>%
                        select(reserved, Gene, Alterations, Drugs, Tier, type)

mytiers         <- rbind(mytierI, mytierII, mytierIIIA, mytierIIIB, mytierIVA1, mytierIVA2, mytierIVB)

# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Tiers to patient-aberration data.frame

aberrations_per_patient$Tier       <- "X"
aberrations_per_patient$Drugs      <- NA
aberrations_per_patient$okb.link   <- NA

for(i in which(mytiers$type %in% c("spec", "aapos"))){
         
         matches                                   <- with(aberrations_per_patient, which(Gene==mytiers$Gene[i] & str_detect(Aberration.name, mytiers$Alterations[i])))
         aberrations_per_patient$Tier[matches]     <- mytiers$Tier[i]
         aberrations_per_patient$Drugs[matches]    <- mytiers$Drugs[i]
         aberrations_per_patient$okb.link[matches] <- mytiers$Alterations[i]
         
         }
         
for(i in which(mytiers$type=="amp")){
         
         matches                                   <- with(aberrations_per_patient, which(Gene==mytiers$Gene[i] & Consequence=="amplification" & Pathogenic & Tier=="X"))
         aberrations_per_patient$Tier[matches]     <- mytiers$Tier[i]
         aberrations_per_patient$Drugs[matches]    <- mytiers$Drugs[i]
         aberrations_per_patient$okb.link[matches] <- mytiers$Alterations[i]
         
         }

for(i in which(mytiers$type=="fus")){
         
         matches                                   <- with(aberrations_per_patient, which(Gene==mytiers$Gene[i] & str_detect(Consequence, "fus|transloca|join|multiple") & Pathogenic & Tier=="X"))
         aberrations_per_patient$Tier[matches]     <- mytiers$Tier[i]
         aberrations_per_patient$Drugs[matches]    <- mytiers$Drugs[i]
         aberrations_per_patient$okb.link[matches] <- mytiers$Alterations[i]
         
         }
         
for(i in which(mytiers$type=="mut")){
         
         matches                                   <- c(with(aberrations_per_patient, which(Gene==mytiers$Gene[i] & Gene.type=="loss" & Pathogenic & Tier=="X")),
                                                        with(aberrations_per_patient, which(Gene==mytiers$Gene[i] & Gene.type=="gain" & 
                                                                                            Consequence %in% c("inframe_deletion", "inframe_duplication", "inframe_insertion", "inframe_substitution", "nonsynonymous_SNV") &
                                                                                            Pathogenic & Tier=="X")))  
         aberrations_per_patient$Tier[matches]     <- mytiers$Tier[i]
         aberrations_per_patient$Drugs[matches]    <- mytiers$Drugs[i]
         aberrations_per_patient$okb.link[matches] <- mytiers$Alterations[i]
         
         }



# Force TP53 mutations to Tier IVA

aberrations_per_patient$Drugs[which(aberrations_per_patient$Tier=="X" & aberrations_per_patient$Gene=="TP53" & aberrations_per_patient$Pathogenic)] <- "Adavosertib"
aberrations_per_patient$Tier[which(aberrations_per_patient$Tier=="X" & aberrations_per_patient$Gene=="TP53" & aberrations_per_patient$Pathogenic)]  <- "IVA"

for(i in which(is.na(aberrations_per_patient$Drugs))) {

                   tempdrugs <- subset(network, proteinTarget==aberrations_per_patient$Gene[i] | gene==aberrations_per_patient$Gene[i], select="drug")
                   tempclen  <- unlist(apply(tempdrugs, 1, function(t){
                                                     if(str_detect(t, "[:digit:]")){
                                                         return(t)
                                                     } else {
                                                         return(str_to_title(t))
                                                     }}))

                  if(length(tempclen)>0){
                      aberrations_per_patient$Drugs[i] <- paste(unique(tempclen), collapse=", ")
                  } else {
                      next()
                  }
}  

with(subset(aberrations_per_patient, Tier=="X" & is.na(Drugs)), table(Pathogenic))
# < table of extent 0 >

aberrations_per_patient <- aberrations_per_patient %>%
                               mutate(Credible = (Gene.type=="gain" & Pathogenic & Expressed) |
                                                 (Gene.type=="loss" & Pathogenic & Homogeneous))



# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Create an annotation matrix for PROGENy and CollecTRI signature analyses
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Prioritize samples
# primary samples from metastases with RNA match - highest purity
# primary samples from adnex with RNA highest purity
# other samples - highest purity

pwsamples  <- cancersamples %>%
                  group_by(patient) %>%
                  mutate(prior1 = !(is.na(RNA.sample)) & str_detect(stime, "p") & ssite %in% c("Bow", "Mes", "Ome", "Oth", "Per", "Ute", "Vag"),
                         prior2 = !(is.na(RNA.sample)) & str_detect(stime, "p"), 
                         prior3 = !(is.na(RNA.sample))) %>%
                  filter(prior1 | prior2 | prior3) %>%
                  summarize(DNA.sample=ifelse(any(prior1), DNA.sample[which(purity==max(purity[which(prior1)]) & prior1)],
                                       ifelse(any(prior2), DNA.sample[which(purity==max(purity[which(prior2)]) & prior2)],
                                                           DNA.sample[which(purity==max(purity))])), 
                            RNA.sample=ifelse(any(prior1), RNA.sample[which(purity==max(purity[which(prior1)]) & prior1)],
                                       ifelse(any(prior2), RNA.sample[which(purity==max(purity[which(prior2)]) & prior2)],
                                                           RNA.sample[which(purity==max(purity))]))) %>%
                 merge(snvtemp %>% filter(Gene=="NF1" & Pathogenic & Homogeneous) %>% mutate(NF1homsnv=TRUE) %>% select(sample, NF1homsnv),
                       by.x="DNA.sample", by.y="sample", all.x=T) %>%
                 mutate(NF1homsnv = ifelse(is.na(NF1homsnv), FALSE, NF1homsnv)) %>%
                 merge(snvtemp %>% filter(Gene=="NF1" & Pathogenic & !Homogeneous & AD.1>0) %>% mutate(NF1hetsnv=TRUE) %>% select(sample, NF1hetsnv),
                       by.x="DNA.sample", by.y="sample", all.x=T) %>%
                 mutate(NF1hetsnv = ifelse(is.na(NF1hetsnv), FALSE, NF1hetsnv)) %>%
                 merge(brktemp %>% filter(Gene=="NF1" & Pathogenic & Homogeneous) %>% mutate(NF1hombrk=TRUE) %>% select(sample, NF1hombrk),
                       by.x="DNA.sample", by.y="sample", all.x=T) %>%
                 mutate(NF1hombrk = ifelse(is.na(NF1hombrk), FALSE, NF1hombrk)) %>%
                 merge(brktemp %>% filter(Gene=="NF1" & Pathogenic & !Homogeneous) %>% mutate(NF1hetbrk=TRUE) %>% select(sample, NF1hetbrk),
                       by.x="DNA.sample", by.y="sample", all.x=T) %>%
                 mutate(NF1hetbrk = ifelse(is.na(NF1hetbrk), FALSE, NF1hetbrk)) %>%
                 merge(deltemp %>% filter(Gene=="NF1") %>% mutate(NF1homdel=TRUE) %>% select(sample, NF1homdel),
                       by.x="DNA.sample", by.y="sample", all.x=T) %>%
                 mutate(NF1homdel = ifelse(is.na(NF1homdel), FALSE, NF1homdel)) %>%
                 merge(snvtemp %>% filter(Gene=="PTEN" & Pathogenic & Homogeneous) %>% mutate(PTENhomsnv=TRUE) %>% select(sample, PTENhomsnv),
                       by.x="DNA.sample", by.y="sample", all.x=T) %>%
                 mutate(PTENhomsnv = ifelse(is.na(PTENhomsnv), FALSE, PTENhomsnv)) %>%
                 merge(snvtemp %>% filter(Gene=="PTEN" & Pathogenic & !Homogeneous & AD.1>0) %>% mutate(PTENhetsnv=TRUE) %>% select(sample, PTENhetsnv),
                       by.x="DNA.sample", by.y="sample", all.x=T) %>%
                 mutate(PTENhetsnv = ifelse(is.na(PTENhetsnv), FALSE, PTENhetsnv)) %>%
                 merge(brktemp %>% filter(Gene=="PTEN" & Pathogenic & Homogeneous) %>% mutate(PTENhombrk=TRUE) %>% select(sample, PTENhombrk),
                       by.x="DNA.sample", by.y="sample", all.x=T) %>%
                 mutate(PTENhombrk = ifelse(is.na(PTENhombrk), FALSE, PTENhombrk)) %>%
                 merge(brktemp %>% filter(Gene=="PTEN" & Pathogenic & !Homogeneous) %>% mutate(PTENhetbrk=TRUE) %>% select(sample, PTENhetbrk),
                       by.x="DNA.sample", by.y="sample", all.x=T) %>%
                 mutate(PTENhetbrk = ifelse(is.na(PTENhetbrk), FALSE, PTENhetbrk)) %>%
                 merge(deltemp %>% filter(Gene=="PTEN") %>% mutate(PTENhomdel=TRUE) %>% select(sample, PTENhomdel),
                       by.x="DNA.sample", by.y="sample", all.x=T) %>%
                 mutate(PTENhomdel = ifelse(is.na(PTENhomdel), FALSE, PTENhomdel)) %>%
                 merge(ampdata %>% filter(Gene=="AKT2" & AMP_ploidy & Expression_quantile>=0.8) %>% mutate(AKT2exp=TRUE) %>% select(sample, AKT2exp),
                       by.x="DNA.sample", by.y="sample", all.x=T) %>%
                 mutate(AKT2exp = ifelse(is.na(AKT2exp), FALSE, AKT2exp)) %>%
                 merge(ampdata %>% filter(Gene=="AKT2" & AMP_ploidy & Expression_quantile<0.8) %>% mutate(AKT2low=TRUE) %>% select(sample, AKT2low),
                       by.x="DNA.sample", by.y="sample", all.x=T) %>%
                 mutate(AKT2low = ifelse(is.na(AKT2low), FALSE, AKT2low)) %>%
                 merge(ampdata %>% filter(Gene=="PTK2" & AMP_ploidy & Expression_quantile>=0.8) %>% mutate(PTK2exp=TRUE) %>% select(sample, PTK2exp),
                       by.x="DNA.sample", by.y="sample", all.x=T) %>%
                 mutate(PTK2exp = ifelse(is.na(PTK2exp), FALSE, PTK2exp)) %>%
                 merge(ampdata %>% filter(Gene=="PTK2" & AMP_ploidy & Expression_quantile<0.8) %>% mutate(PTK2low=TRUE) %>% select(sample, PTK2low),
                       by.x="DNA.sample", by.y="sample", all.x=T) %>%
                 mutate(PTK2low = ifelse(is.na(PTK2low), FALSE, PTK2low)) %>%
                 merge(ampdata %>% filter(Gene=="KRAS" & AMP_ploidy & Expression_quantile>=0.8) %>% mutate(KRASexp=TRUE) %>% select(sample, KRASexp),
                       by.x="DNA.sample", by.y="sample", all.x=T) %>%
                 mutate(KRASexp = ifelse(is.na(KRASexp), FALSE, KRASexp)) %>%
                 merge(ampdata %>% filter(Gene=="KRAS" & AMP_ploidy & Expression_quantile<0.8) %>% mutate(KRASlow=TRUE) %>% select(sample, KRASlow),
                       by.x="DNA.sample", by.y="sample", all.x=T) %>%
                 mutate(KRASlow = ifelse(is.na(KRASlow), FALSE, KRASlow)) %>%
                 merge(ampdata %>% filter(Gene=="CCNE1" & AMP_ploidy & Expression_quantile>=0.8) %>% mutate(CCNE1exp=TRUE) %>% select(sample, CCNE1exp),
                       by.x="DNA.sample", by.y="sample", all.x=T) %>%
                 mutate(CCNE1exp = ifelse(is.na(CCNE1exp), FALSE, CCNE1exp)) %>%
                 merge(ampdata %>% filter(Gene=="CCNE1" & AMP_ploidy & Expression_quantile<0.8) %>% mutate(CCNE1low=TRUE) %>% select(sample, CCNE1low),
                       by.x="DNA.sample", by.y="sample", all.x=T) %>%
                 mutate(CCNE1low = ifelse(is.na(CCNE1low), FALSE, CCNE1low)) %>%

                 mutate(NF1loss = ifelse(NF1homsnv | NF1hombrk | NF1homdel, "full", ifelse(NF1hetsnv | NF1hetbrk, "marginal", "normal")),
                        PTENloss = ifelse(PTENhomsnv | PTENhombrk | PTENhomdel, "full", ifelse(PTENhetsnv | PTENhetbrk, "marginal", "normal")),
                        AKT2amp = ifelse(AKT2exp, "hyperexpressed", ifelse(AKT2low, "amp_normal", "normal")),
                        PTK2amp = ifelse(PTK2exp, "hyperexpressed", ifelse(PTK2low, "amp_normal", "normal")),
                        KRASamp = ifelse(KRASexp, "hyperexpressed", ifelse(KRASlow, "amp_normal", "normal")),
                        CCNE1amp = ifelse(CCNE1exp, "hyperexpressed", ifelse(CCNE1low, "amp_normal", "normal"))) %>%
                select(patient, DNA.sample, RNA.sample, NF1loss, PTENloss, AKT2amp, PTK2amp, KRASamp, CCNE1amp)

write.table(pwsamples, file="<desired_path>/sample_summaries_driver_events.tsv", col.names=T, row.names=F, quote=F, sep="\t")

# Output Supplementary data file with patient-level mutation data

Supplementary_data2 <- aberrations_per_patient %>%
                       merge(patdata %>% select(Patient, Publication_code), by="Patient", all.x=T) %>%
                       mutate(Patient=Publication_code) %>%
                       select(Patient, Gene, Gene.type, Aberration.name, Genomic.position, Consequence, CN.normal.range, CN.affected.range, CredOme, CredPer, CredAdn, Relapse, Allsamples, Pathogenic, Expressed, Homogeneous, Tier, Drugs, Credible) %>%
                       arrange(Patient, Gene)

write.table(Supplementary_data2, file="<desired_path>/Supplementary_data2.csv", col.names=T, row.names=F, quote=T, sep=",")
                     
                       
