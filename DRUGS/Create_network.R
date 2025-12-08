#!/usr/bin/env R
library("tidyverse")

setwd("path_to_network")

# The drug-target-biomarker-mutation network is created by combining information from the Open targets platform and OncoKB database.
# This code should be run after #### for filtering the Open targets data for cancer-targeting small-molecule drugs with a protein target, excluding generic chemotherapies and anti-microbiological drugs (lists provided in the repository).
# The expected consequence types for oncogenic aberrations in a gene are provided via a separate, curated file. (provided in the repository)
# Some OncoKB drugs did not pass the filtering criteria for Open targets and these were salvaged from the data separately, and added as an external curated file. (not provided in the repository)

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Read in a comma-separated file from processing of the Open targets data with the following columns. Add the salvaged information.
# chemblId                                    - Drug identifier in the Chembl database ("CHEMBL[:digit:]+")               
# drug                                        - Drug generic name ([:upper:]+)
# drugType                                    - Drug type
# action                                      - Type of drug action
# targetENSG                                  - Ensembl gene id of the gene encoding the drug target protein
# proteinTarget                               - HUGO symbol of the gene encoding the target protein
# gene                                        - HUGO symbol of the biomarker gene linked to the drug-protein pair
# GoF_LoF                                     - Consequence of an oncogenic mutation in the biomarker gene. GoF for gain-of-function and LoF for loss-of-function.
# maximumClinicalTrialPhase                   - Maximum phase of clinical trials for the drug (1-4)
# hasBeenWithdrawn                            - Has the drug been withdrawn from clinical trials or from clinical use (TRUE/FALSE)
# yearOfFirstApproval                         - First year of approval for the drug. N.B. drugs approved before 1990 were systematically excluded with an assumption that the anti-cancer action is generic.
# isApproved                                  - Is the drug approved for clinical use (TRUE/FALSE)
# ApprovedAfter1989                           - Was the drug approved after 1989 (TRUE/FALSE)
# rightDrugType                               - Does the drug type match the selection criteria (TRUE/FALSE)
# havingProteinTarget                         - Does the drug have a physical protein target (TRUE/FALSE)
# notAntimicrobiologicalDrug                  - Is the drug not an anti-microbiological agent (TRUE/FALSE)
# noGenericChemo                              - Is the drug not a generic chemotherapeutic compound (TRUE/FALSE)
# terminated                                  - Was a phase I clinical trial for the drug terminated (TRUE/FALSE)
# potentialEffect                             - Is there a credible functional connection between the drug and the target protein (TRUE/FALSE)

OTcancer  <- read.delim("OTcancer.csv", header=T, as.is=T, sep=",") %>%
             mutate(ApprovedAfter1989=ifelse(is.na(ApprovedAfter1989), T, ApprovedAfter1989)) %>%
             add_row(read.delim("Salvaged_OT_cancerdrug_data.csv", header=T, as.is=T, sep=","))  %>%
             mutate(drugProteinRelation = action,
                    geneProteinRelation = ifelse(!hasBeenWithdrawn & ApprovedAfter1989 & rightDrugType & havingProteinTarget & 
                                                  notAntimicrobiologicalDrug & noGenericChemo & !terminated & potentialEffect, "activation", NA)) %>%
             select(-action) %>%
             select(chemblId, drug, drugType, drugProteinRelation, proteinTarget, targetENSG, geneProteinRelation, gene, GoF_LoF, maximumClinicalTrialPhase, yearOfFirstApproval,
                    isApproved, hasBeenWithdrawn, terminated, noGenericChemo, notAntimicrobiologicalDrug, ApprovedAfter1989, rightDrugType, havingProteinTarget, potentialEffect)


# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Read in biomarker gene functional classification with column names
# Gene                                        - HUGO gene symbol
# GoF_LoF                                     - Consequence of an oncogenic mutation in the biomarker gene. GoF for gain-of-function and LoF for loss-of-function.
# Is.Oncogene                                 - Annotated as an oncogene in OncoKB v4.1.
# Is.Tumor.Suppressor.Gene                    - Annotated as a tumor suppressor gene in OncoKB v4.1.

genfunc   <- read.delim("genefunctionkey.csv", header=T, as.is=T, sep="\t", na.strings=c("", "NA")))

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Read in OncoKB v4.1 data with the following columns.
# Level                                       - FDA therapeutic level for the Gene-treatment combination
# Gene                                        - HUGO gene symbol
# Alterations                                 - Gene mutations, specific mutations or mutation categories
# Cancer.Types                                - Cancer histologies
# Drugs..for.therapeutic.implications.only.   - Drug combinations

OKBdata   <- read.delim("OKBdata_v4.1.csv", header=T, as.is=T, sep="\t", na.strings=c("", "NA"))  %>%
             select(Level, Gene, Alterations, Cancer.Types, Drugs=Drugs..for.therapeutic.implications.only.) %>%
             mutate(Drug = str_split(Drugs, pattern="\\s\\+\\s")) %>%
             unnest(c(Drug)) %>%
             merge(genfunc %>% select(Gene, GoF_LoF), by="Gene", all.x=T)  %>%
             filter(Level %in% as.character(1:4)) %>%
             mutate(groupkey=paste(Drug, Gene, sep="_"),
                    Level=as.numeric(Level)) %>%
             group_by(groupkey) %>%
             summarize(Gene=unique(Gene),
                       GoF_LoF=unique(GoF_LoF),
                       Drug=unique(Drug),
                       Level=max(Level),
                       Alterations=paste(unique(Alterations[which.max(Level)]), collapse=","),
                       Cancer.Types=paste(unique(Cancer.Types[which.max(Level)]), collapse=","),
                       Drugs=paste(unique(Drugs[which.max(Level)]), collapse=","))  %>%
             mutate(OTkey=ifelse(Drug=="ABI-009", "SIROLIMUS", 
                          ifelse(Drug=="GSK2636771", "GSK-2636771",
                          ifelse(Drug=="AZD8186", "AZD-8186", str_to_upper(Drug)))))

# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Combine data

NETWORKraw      <- do.call(rbind, lapply(unique(OTcancer$drug), function(md) {

                     ottemp <- subset(OKBdata, drug==md)
                      
                      if(md %in% OKBdata$OTkey) {
                      
                          okbtemp <- subset(OKBdata, OTkey==md)
                          
                          if(any(ottemp$gene %in% okbtemp$Gene)){
                          
                              comb1 <- merge(ottemp, okbtemp, by.x="gene", by.y="Gene", all.x=T) %>%
                                          mutate(GoF_LoF=ifelse(!is.na(GoF_LoF.y), GoF_LoF.y,
                                                              ifelse(GoF_LoF.x=="gain-of-function", "GoF", 
                                                                   ifelse(GoF_LoF.x=="loss-of-function", "LoF", GoF_LoF.x))),
                                                   okb.drug=as.character(na.omit(unique(okbtemp$Drug)))) %>%
                                          select(chemblId, drug, drugType, drugProteinRelation, proteinTarget, targetENSG, geneProteinRelation, gene, GoF_LoF, maximumClinicalTrialPhase, 
                                                 yearOfFirstApproval, isApproved,hasBeenWithdrawn, terminated, GenericChemo, AntimicrobiologicalDrug, ApprovedBefore1990, rightDrugType, 
                                                 havingProteinTarget, potentialEffect, okb.drug, fda.level=Level, Alterations, Cancer.Types, Drugs)
                         
                               if(any(!(okbtemp$Gene %in% ottemp$gene))){
                               
                                   okbtemp <- subset(okbtemp, !(Gene %in% ottemp$gene))
                         
                                   comb2  <- do.call(rbind, lapply(1:nrow(okbtemp), function(i) {
                                                 return(cbind(comb1[,1:7], gene=okbtemp$Gene[i], GoF_LoF=okbtemp$GoF_LoF[i], comb1[,10:20],
                                                             okb.drug=okbtemp$Drug[i], fda.level=okbtemp$Level[i], Alterations=okbtemp$Alterations[i],
                                                             Cancer.Types=okbtemp$Cancer.Types[i], Drugs=okbtemp$Drugs[i]))}))
                         
                                   return(rbind(comb1, comb2))
 
                               } else {
                               
                                   return(comb1)
                               }
                               
                          } else {
                          
                               comb  <-  rbind(cbind(ottemp, okb.drug=NA, fda.level=NA, Alterations=NA,
                                                         Cancer.Types=NA, Drugs=NA),
                                         do.call(rbind, lapply(1:nrow(okbtemp), function(i) {
                                             return(cbind(ottemp[,1:6], geneProteinRelation=NA, gene=okbtemp$Gene[i], GoF_LoF=okbtemp$GoF_LoF[i], ottemp[,10:20],
                                                         okb.drug=okbtemp$Drug[i], fda.level=okbtemp$Level[i], Alterations=okbtemp$Alterations[i],
                                                         Cancer.Types=okbtemp$Cancer.Types[i], Drugs=okbtemp$Drugs[i]))})))
                                                         
                                return(comb)
                                                         
                          }
                          
                     } else {
                     
                         return(cbind(ottemp, okb.drug=NA, fda.level=NA, Alterations=NA, Cancer.Types=NA, Drugs=NA))

                     }
                  }
                  ))
# Some hard-coded curation for specific cases

NETWORK        <- NETWORKraw %>%
                      mutate(proteinTarget = ifelse(drug=="LUNRESERTIB", "PKMYT1", proteinTarget),
                             targetENSG = ifelse(drug=="LUNRESERTIB", "ENSG00000127564", targetENSG),
                             geneProteinRelation = ifelse(drug=="LUNRESERTIB", "dependency", geneProteinRelation),
                             havingProteinTarget = ifelse(drug=="LUNRESERTIB", TRUE, havingProteinTarget),
                             GenericChemo = ifelse(str_detect(drug, "PLATIN|RACRIL"), TRUE, GenericChemo)) %>%
                      filter(!(chemblId %in% c("CHEMBL3545096", "CHEMBL5314415", "CHEMBL5095024"))) %>%
                      mutate(drugProteinRelation=ifelse((drug=="ABIRATERONE" & proteinTarget=="AR") | (drug=="PBF-999" & proteinTarget=="ADORA2A"), "ANTAGONIST",
                                                 ifelse(drug=="ABIRATERONE" | drug=="PBF-999", "INHIBITOR", drugProteinRelation))) %>%
                      mutate(drugProteinRelation=ifelse(is.na(drugProteinRelation) & drug!="REZATAPOPT", "INHIBITOR",
                                                 ifelse(drug=="REZATAPOPT", "MODULATOR", drugProteinRelation))) %>%
                      mutate(geneProteinRelation=ifelse(is.na(geneProteinRelation) & str_detect(drug, "PARIB"), "dependency", geneProteinRelation)) %>%
                      mutate(GoF_LoF=ifelse(GoF_LoF %in% c("flagged out", "gain-of-function"), "GoF", GoF_LoF)) %>%
                      filter(GoF_LoF!="ambiguous") %>%
                      mutate(geneProteinRelation=ifelse(is.na(geneProteinRelation) & drug %in% c("ABIRATERONE", "ABIRATERONE ACETATE", "ATEZOLIZUMAB", "BEVACIZUMAB", "DOSTARLIMAB", "ENZALUTAMIDE", "IPILIMUMAB", "NIVOLUMAB", "PEMBROLIZUMAB", "RAMUCIRUMAB"), "redundant",
                                                 ifelse(is.na(geneProteinRelation) & drug %in% c("CERALASERTIB", "EMAVUSERTIB"), "dependency",
                                                 ifelse(is.na(geneProteinRelation), "activation", geneProteinRelation)))) %>%
                      filter(geneProteinRelation!="redundant")  %>%
                      mutate(yearOfFirstApproval=ifelse(drug=="ABIRATERONE", 2011,
                                                 ifelse(drug=="IXAZOMIB", 2015,
                                                 ifelse(drug=="TISLELIZUMAB", 2024, yearOfFirstApproval)))) %>%
                      filter(!hasBeenWithdrawn & !terminated & !GenericChemo & notAntimicrobiologicalDrug & ApprovedAfter1989 & rightDrugType & havingProteinTarget & potentialEffect & !is.na(proteinTarget)) %>%
                      select(-hasBeenWithdrawn, -terminated, -GenericChemo, -notAntimicrobiologicalDrug, -ApprovedAfter1989, -rightDrugType, -havingProteinTarget, -potentialEffect) %>%
                      distinct()

# Output

write.table(NETWORK, file="NETWORK.csv", col.names=T, row.names=F, quote=F, sep="\t")