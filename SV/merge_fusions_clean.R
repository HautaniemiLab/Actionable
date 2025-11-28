#!/usr/bin/env R

library("tidyverse")

# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Read in sample data file to manage data accessibility
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

cancersamples        <- read.delim("cancersamplesfile.csv", header=T, as.is=T, sep="\t")

# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Read in a file with target genes, classified according to the effect that a mutation has to have in that gene for being oncogenic, either gain or loss.
# While two genes may be affected by the same break event and any event is counted only once, the genes with fda-level indication for targeted therapy are prioritized.
# 
# The file must have the following columns
# gene         Hugo symbol
# GoF_LoF      Either GoF for gain-of-function, LoF for loss-of-function, or GoF_LoF for either.
# fda.level    Therapeutic level 1-4 as in OncoKB (https://www.oncokb.org/therapeutic-levels), or NA.
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

network              <- read.delim("networkfile.csv", header=T, as.is=T, sep="\t")

lossgenes            <- sort(unique(network$gene[which(network$GoF_LoF=="LoF")]))
gaingenes            <- sort(unique(network$gene[which(network$GoF_LoF=="GoF")]))
lossgenes            <- sort(c(lossgenes, "TP53"))
gaingenes            <- sort(c(gaingenes, unique(network$gene[which(network$GoF_LoF=="GoF_LoF")])))
gaingenes            <- gaingenes[-which(gaingenes %in% lossgenes)]
gainprim             <- with(subset(network, GoF_LoF=="GoF" & !is.na(fda.level)), unique(gene))

# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Read in a file with annotated gene breakpoints from DNA analysis
# 
# The file must have the following columns
# sample                    Sample ID
# patient                   Patient ID
# svId                      Structural variant ID from Purple/GRIDSS analysis
# resolvedType              Complex type of the structural variant
# Type                      Simple type of the structural variant
# gene1                     Hugo symbol of the gene at the earlier genomic location of the break
# transcriptId1             Transcript ID of the gene at the earlier genomic location of the break
# regionType1               DOWNSTREAM|EXONIC|INTRONIC|UPSTREAM position in relation to the gene at the earlier genomic location of the break
# strand1                   Coding strand of the gene at the earlier genomic location of the break, -1|1, at the earlier genomic location of the break
# undisruptedCopyNumber1    Estimated number of unaffected gene copies at the break, at the earlier genomic location of the break
# affectedCopyNumber1       Estimated number of affected gene copies at the break, at the earlier genomic location of the break
# exonicBasePhase1          Phase of the exon, which is assumed to take part in a potential fusion of coding exons, at the earlier genomic location of the break
# exonJoined1               Identifier of the exon, which is assumed to take part in a potential fusion of coding exons, at the earlier genomic location of the break
# chromosome1               Chromosome at the earlier genomic location of the break
# pos1                      Chromosomal position at the earlier genomic location of the break
# gene2                     Hugo symbol of the gene at the later genomic location of the break
# transcriptId2             Transcript ID of the gene at the later genomic location of the break
# regionType2               DOWNSTREAM|EXONIC|INTRONIC|UPSTREAM position in relation to the gene at the later genomic location of the break
# strand2                   Coding strand of the gene at the earlier genomic location of the break, -1|1, at the later genomic location of the break
# undisruptedCopyNumber2    Estimated number of unaffected gene copies at the break, at the later genomic location of the break
# affectedCopyNumber2       Estimated number of affected gene copies at the break, at the later genomic location of the break
# exonicBasePhase2          Phase of the exon, which is assumed to take part in a potential fusion of coding exons, at the later genomic location of the break
# exonJoined2               Identifier of the exon, which is assumed to take part in a potential fusion of coding exons, at the later genomic location of the break
# chromosome2               Chromosome at the later genomic location of the break
# pos2                      Chromosomal position at the later genomic location of the break
# affected_gene             Hugo symbol of the gene, whose coding sequence is affected by the break event. N.B. If several genes are affected, the event is multiplied and the annotation is provided for each of the genes on separate rows. These can be pooled with svId.
# effect                    The predicted effect for the affected gene. E.g. fusion, truncated_transcript, fusion_from_unknown, promoter_hijacking, simple_duplication, simple_deletion
# frame                     frameshift|in_frame
# DNA_change                Gene-based name for the structural variant
# mutation                  Coordinate-based name for the structural variant

# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

brks                 <- read.delim("breakpointfile.tsv", header=T, as.is=T, sep="\t") %>%
                            filter((gene1 %in% c(gaingenes, lossgenes) & str_detect(affected_gene, "gene1")) | (gene2 %in% c(gaingenes, lossgenes) & str_detect(affected_gene, "gene2"))) %>%
                            distinct() %>%
                            mutate(patient = str_extract(sample, "^[:alnum:]+"),
                                   temp = str_extract(sample, "(?<=_)[poir][:alnum:]+")) %>%
                            mutate(stime = str_extract(temp, "^[poir][:digit:]?"),
                                   ssite = str_extract(temp, "(?<=[poir][:digit:]?)([:upper:][:lower:]{2}|LN)"),
                                   sside = str_extract(temp, "(?<=([:upper:][:lower:]{2}|LN))[LR]"),
                                   sord  = str_extract(temp, "[:digit:]$")) %>%
                            select(-temp) %>%
                            mutate(temp5prime = str_extract(DNA_change, "^.+(?=\\:\\:)"),
                                   temp3prime = str_extract(DNA_change, "(?<=\\:\\:).+$")) %>%
                            mutate(gene5prime = str_extract(temp5prime, "^[[:alnum:]\\-\\.]+(?=\\:)"),
                                   gene3prime = str_extract(temp3prime, "^[[:alnum:]\\-\\.]+(?=\\:)")) %>%
                            select(-temp5prime, -temp3prime) %>%
                            mutate(flag = ifelse(gene1==gene2 & transcriptId1!=transcriptId2, "cross-transcript", NA)) %>%
                            filter(is.na(flag)) %>%
                            select(-flag)
              

# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Read in a file with gene fusion data from RNA analysis
# 
# The file must have the following columns
# Sample                    Sample ID
# Gene1                     5' gene of the fusion transcript
# Gene2                     3' gene of the fusion transcript
# Chromosome1               Chromosome of the 5' gene
# Chromosome2               Chromosome of the 3' gene
# Breakpoint1               Genomic coordinate of the fusion in the 5' gene
# Breakpoint2               Genomic coordinate of the fusion in the 5' gene
# FusionStrand1             Always the same as GeneStrand1, the strand of the 5' gene, -|+
# FusionStrand2             Genomic strand of the 3' sequence, not necessarily same as GeneStrand2 -|+
# GeneStrand1               Genomic strand of the 5' gene, -|+
# GeneStrand2               Genomic strand of the 3' gene, -|+
# EncompassingReads         Number of read pairs, where the read and mate are aligned to different fusion partner sequences
# JunctionReads             Number of reads crossing the fusion junction
# SupportingReads           Total number of supporting reads, sum of JunctionReads and EncompassingReads
# type                      Estimated type of the genomic event causing the fusion, deletion|deletion/read-through|duplication|inversion|translocation
# site1                     Description of the fusion site in the 5' gene: 3'UTR|5'UTR/splice-site|CDS|CDS/splice-site|exon|exon/splice-site|intergenic|intron 
# site2                     Description of the fusion site in the 3' gene: 3'UTR|5'UTR/splice-site|CDS|CDS/splice-site|exon|exon/splice-site|intergenic|intron 
# Ensembl1                  Ensembl gene id for the 5' gene
# Ensembl2                  Ensembl gene id for the 3' gene
# transcript_id1            Ensembl transcript id for the 5' gene
# transcript_id2            Ensembl transcript id for the 3' gene
# reading_frame             Preferred estimate on consequence of the fusion for reading frame, in-frame|out-of-frame|stop-codon or NA
# PROT_FUSION_TYPE          Alternative estimate on cConsequence of the fusion for reading frame, INFRAME|FRAMESHIFT|stop-codon or NA
# peptide_sequence          Predicted sequence until of the fusion protein from the start codon of the 5' gene until the first stop codon
# FUSION_TRANSL             Full-length sequence of the fusion protein, including multiple stop codons

# Create automatic interpretation of the fusions
# 1. In-frame events are interpreted fusions, no matter whether the 5' or 3' gene is the gene of interest
# 2. Frameshift events are interpreted as truncations of the 5'gene, if the 5' gene is the gene of interest
# 3. UTR events are ignored
# 4. If the 3' gene is the gene of interest and the fusion causes a frameshift, but the 5' sequence is shorter than 25 amino acids, the 3' is interpreted as a potential low-molecular-weight protein with novel ORF start codon, LMW_Cterm.
# 5. If the 3' gene is the gene of interest and the coding sequence of the 5' partner is undefined, the 3' is interpreted as a potential low-molecular-weight protein with novel ORF start codon.
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

fusions              <- read.delim("fusionfile.csv", header=T, as.is=T, na.strings=c(".", "", "NA", " "), sep=",") %>%
                            mutate(patient = str_extract(Sample, "^[:alnum:]+"),
                                   temp = str_extract(Sample, "(?<=_)[poir][:alnum:]+")) %>%
                            mutate(stime = str_extract(temp, "^[poir][:digit:]?"),
                                   ssite = str_extract(temp, "(?<=[poir][:digit:]?)([:upper:][:lower:]{2}|LN)"),
                                   sside = str_extract(temp, "(?<=([:upper:][:lower:]{2}|LN))[LR]"),
                                   sord  = str_extract(temp, "[:digit:]$")) %>%
                            select(-temp) %>%
                            mutate(gene1act = ifelse(Gene1 %in% gaingenes, "gain",
                                              ifelse(Gene1 %in% lossgenes, "loss", NA)),
                                   gene2act = ifelse(Gene2 %in% gaingenes, "gain",
                                              ifelse(Gene2 %in% lossgenes, "loss", NA))) %>%
                            filter((!is.na(gene1act) | !is.na(gene2act))) %>%
                            filter((str_detect(site1, "CDS|exon|intron") & !is.na(gene1act)) |
                                   (str_detect(site2, "CDS|exon|intron") & !is.na(gene2act))) %>%
                            mutate(group_key=paste(Chromosome1, Breakpoint1, Chromosome2, Breakpoint2, sep="_")) %>%
                            group_by(group_key) %>%
                            mutate(reading_frame=ifelse(is.na(reading_frame), paste(c(na.omit(unique(reading_frame))), collapse=";"), reading_frame),
                                   peptide_sequence=ifelse(is.na(peptide_sequence), paste(c(na.omit(unique(peptide_sequence))), collapse=";"), peptide_sequence)) %>%
                            ungroup %>%
                            mutate(reading_frame=ifelse(reading_frame=="", NA, reading_frame),
                                   peptide_sequence=ifelse(peptide_sequence=="", NA, peptide_sequence)) %>%
                            mutate(consequence = ifelse(reading_frame=="in-frame" & !is.na(reading_frame), "fusion_transcript",
                                                 ifelse(!is.na(gene1act) & site1!="5'UTR/splice-site", "truncation", 
                                                 ifelse(!is.na(gene1act), "redundant",
                                                 ifelse(!is.na(gene2act) & !is.na(reading_frame) & str_detect(FUSION_TRANSL, "^[:upper:]{1,24}\\*"), "LMW_Cterm",
                                                 ifelse(!is.na(gene2act) & !is.na(reading_frame), "redundant",
                                                 ifelse(!is.na(gene2act) & site1=="3'UTR", "redundant",
                                                 ifelse(!is.na(gene2act), "LMW_Cterm", NA)))))))) %>%
                            mutate(consequence = ifelse(!is.na(peptide_sequence) & consequence!="redundant" & str_detect(peptide_sequence, "\\|[:alpha:]{0,9}\\*"), "truncation", consequence)) %>%
                            mutate(acttarget = ifelse(!is.na(gene1act), Gene1, Gene2))
fusions              <- fusions %>%
                            add_row(fusions %>% filter(!is.na(gene1act) & !is.na(gene2act)) %>%
                                    mutate(consequence = ifelse(str_detect(site2, "CDS"), "LMW_Cterm", "redundant"),
                                           acttarget = Gene2))



# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Read in a file with intragenic fusion transcripts - expected consequence of intragenic structural events
# 
# The colnames are the same as for intergenic fusions, above
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

splices              <- read.delim("splicefile.csv", header=T, as.is=T, na.strings=c(".", "", "NA", " "), sep=",") %>%
                            mutate(patient = str_extract(Sample, "^[:alnum:]+"),
                                   temp = str_extract(Sample, "(?<=_)[poir][:alnum:]+")) %>%
                            mutate(stime = str_extract(temp, "^[poir][:digit:]?"),
                                   ssite = str_extract(temp, "(?<=[poir][:digit:]?)([:upper:][:lower:]{2}|LN)"),
                                   sside = str_extract(temp, "(?<=([:upper:][:lower:]{2}|LN))[LR]"),
                                   sord  = str_extract(temp, "[:digit:]$")) %>%
                            mutate(geneact=ifelse(Gene1 %in% gaingenes, "gain", ifelse(Gene1 %in% lossgenes, "loss", NA))) %>%
                            mutate(consequence=ifelse(!is.na(reading_frame) & reading_frame=="in-frame", "neo_isoform", 
                                               ifelse(!is.na(reading_frame), "truncation", 
                                               ifelse(is.na(GeneStrand2), "truncation", NA)))) %>%
                            select(-temp) %>%
                            distinct()

# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Match RNA-fusion data to DNA-breakpoint data based on sample site, time, and affected gene


combsamples          <- cbind(brks,
                              do.call(rbind, lapply(1:nrow(brks), function(i) {

                                  event <- brks[i,]
                                  goi   <- ifelse(event$gene1 %in% network$gene, event$gene1, event$gene2)

                                  # Ignore within-gene indels, keep within-gene inversions
                                  if(str_detect(event$DNA_change, "del|dup")){
                                            return(data.frame(RNA.sample=NA, RNA.change=NA, RNA.type=NA, RNA.consequence=NA, RNA.reads=NA, RNA.direct=NA))
                                  # Ignore if no match for the DNA sample is in fusions table
                                  } else if(is.na(event$sside) & nrow(subset(fusions, patient==event$patient & stime==event$stime & ssite==event$ssite))==0) {
                                            return(data.frame(RNA.sample=NA, RNA.change=NA, RNA.type=NA, RNA.consequence=NA, RNA.reads=NA, RNA.direct=NA))
                                  } else if(!is.na(event$sside) & nrow(subset(fusions, patient==event$patient & stime==event$stime & ssite==event$ssite & sside==event$sside))==0) {
                                            return(data.frame(RNA.sample=NA, RNA.change=NA, RNA.type=NA, RNA.consequence=NA, RNA.reads=NA, RNA.direct=NA))
                                  } else {

                                       # Match by patient, sample time and site, and the gene of interest.
                                       if(is.na(event$sside)){
                                           fusion =  subset(fusions, patient==event$patient & (acttarget==event$gene1 | acttarget==event$gene2) & stime==event$stime & ssite==event$ssite & consequence!="redundant")
                                       } else {
                                           fusion =  subset(fusions, patient==event$patient & (acttarget==event$gene1 | acttarget==event$gene2) & stime==event$stime & ssite==event$ssite & sside==event$sside & consequence!="redundant")
                                       }
                              
                                       if(nrow(fusion)==0){
                                           # Return empty row, if no matches
                                           return(data.frame(RNA.sample=NA, RNA.change=NA, RNA.type=NA, RNA.consequence=NA, RNA.reads=0, RNA.direct=NA))
                                       } else { 

                                          # Check that the orientation of the genes, which gene is the 5' and which the 3', matches in DNA and RNA                                  
                                           if(!is.na(event$gene5prime) & event$gene5prime==goi){
                                               orientationmatch  <- fusion$Gene1==event$gene5prime
                                           } else if(!is.na(event$gene3prime) & event$gene3prime==goi) {
                                               orientationmatch  <- fusion$Gene2==event$gene3prime
                                           }
                                           # Flag events with a better match
                                           dopduplicate      <- unlist(lapply(1:nrow(fusion), function(j) {
                                                                  with(brks, any(patient==fusion$patient[j] & stime==fusion$stime[j] & ssite==fusion$ssite[j] & str_detect(DNA_change, paste(fusion$Gene1[j], ".+", fusion$Gene2[j], sep="")))) &
                                                                  !str_detect(event$DNA_change, paste(fusion$Gene1[j], ".+", fusion$Gene2[j], sep=""))
                                                                  }))
                                           # Option to check that the other end of the genomic break and fusion transcript match
                                           # However, complex genomic events involving multiple locations would be filtered, and this option is not applied, below.
                                           if(any(fusion$Gene1==goi) & !is.na(event$gene1) & event$gene1==goi){
                                               positionmatch     <- paste("chr", fusion$Chromosome2, sep="")==event$chromosome2 & fusion$Breakpoint2>event$pos2-500000 & fusion$Breakpoint2<event$pos2+500000
                                           } else if(any(fusion$Gene1==goi) & !is.na(event$gene2) & event$gene2==goi){
                                               positionmatch     <- paste("chr", fusion$Chromosome2, sep="")==event$chromosome1 & fusion$Breakpoint2>event$pos1-500000 & fusion$Breakpoint2<event$pos1+500000
                                           } else if(any(fusion$Gene2==goi) & !is.na(event$gene1) & event$gene1==goi){
                                               positionmatch     <- paste("chr", fusion$Chromosome1, sep="")==event$chromosome2 & fusion$Breakpoint1>event$pos2-500000 & fusion$Breakpoint1<event$pos2+500000
                                           } else if(any(fusion$Gene2==goi) & !is.na(event$gene2) & event$gene2==goi){
                                               positionmatch     <- paste("chr", fusion$Chromosome1, sep="")==event$chromosome1 & fusion$Breakpoint1>event$pos1-500000 & fusion$Breakpoint1<event$pos1+500000
                                           } else {
                                               positionmatch==rep(FALSE, nrow(fusion))
                                           }
                                           
                                           if(any(orientationmatch & !dopduplicate & positionmatch)){       
                                               fusion            <- fusion[which(orientationmatch & !dopduplicate & positionmatch),]
                                               fusion$name       <- with(fusion, paste(Gene1, ":chr", Chromosome1, ":", Breakpoint1, "::", Gene2, ":chr", Chromosome2, ":", Breakpoint2, sep=""))
                                               return(data.frame(RNA.sample      = paste(unique(fusion$Sample),collapse=""), 
                                                                 RNA.change      = paste(unique(fusion$name),collapse="-"),
                                                                 RNA.type        = paste(unique(fusion$type),collapse="-"),
                                                                 RNA.consequence = paste(unique(fusion$consequence),collapse="-"),
                                                                 RNA.reads       = with(fusion %>% mutate(gkey = paste(Chromosome1, Chromosome2, Breakpoint1, Breakpoint2, sep="_")) %>% group_by(gkey) %>% summarize(nreads=max(SupportingReads)), sum(nreads)),
                                                                 RNA.direct      = paste(unique(positionmatch), collapse="-")))

                                           } else if(any(orientationmatch & !dopduplicate)){
                                               fusion            <- fusion[which(orientationmatch & !dopduplicate),]
                                               fusion$name       <- with(fusion, paste(Gene1, ":chr", Chromosome1, ":", Breakpoint1, "::", Gene2, ":chr", Chromosome2, ":", Breakpoint2, sep=""))
                                               return(data.frame(RNA.sample      = paste(unique(fusion$Sample),collapse=""), 
                                                                 RNA.change      = paste(unique(fusion$name),collapse="-"),
                                                                 RNA.type        = paste(unique(fusion$type),collapse="-"),
                                                                 RNA.consequence = paste(unique(fusion$consequence),collapse="-"),
                                                                 RNA.reads       = with(fusion %>% mutate(gkey = paste(Chromosome1, Chromosome2, Breakpoint1, Breakpoint2, sep="_")) %>% group_by(gkey) %>% summarize(nreads=max(SupportingReads)), sum(nreads)),
                                                                 RNA.direct      = paste(unique(positionmatch), collapse="-")))

                                            } else {
                                               return(data.frame(RNA.sample=NA, RNA.change=NA, RNA.type=NA, RNA.consequence=NA, RNA.reads=0, RNA.direct=NA))
                                            }
                                       }
                                   }
                              })))

# ---------------------------------------------------------------------------------
# Add the intragenic events

WIGsupport           <- do.call(rbind, lapply(1:nrow(brks), function(i) {
                        
                            event <- brks[i,]
                            # Match only intragenic events, skip others by returning an empty row.
                            if(is.na(event$gene1) | is.na(event$gene2) | event$gene1!=event$gene2){
                                return(data.frame(WIG.sample=NA, WIG.change=NA, WIG.type=NA, WIG.consequence=NA, WIG.reads=NA))

                            # If no matches, return an empty row.
                            } else if(is.na(event$sside) & nrow(subset(splices, patient==event$patient & stime==event$stime & ssite==event$ssite & Gene1==event$gene1))==0) {
                                return(data.frame(WIG.sample=NA, WIG.change=NA, WIG.type=NA, WIG.consequence=NA, WIG.reads=NA))
                            } else if(!is.na(event$sside) & nrow(subset(splices, patient==event$patient & stime==event$stime & ssite==event$ssite & sside==event$sside & Gene1==event$gene1))==0) {
                                return(data.frame(WIG.sample=NA, WIG.change=NA, WIG.type=NA, WIG.consequence=NA, WIG.reads=NA))
                            } else {
    
                                if(is.na(event$sside)){
                                    wig =  subset(splices, patient==event$patient & Gene1==event$gene1 & stime==event$stime & ssite==event$ssite)
                                } else {
                                    wig =  subset(splices, patient==event$patient & Gene1==event$gene1 & stime==event$stime & ssite==event$ssite & sside==event$sside)
                                }
                                wig$name  <- with(wig, paste(Gene1, ":chr", Chromosome1, ":", Breakpoint1, "_", Breakpoint2, sep=""))
                                return(data.frame(WIG.sample=paste(unique(wig$Sample), collapse="|"), 
                                                  WIG.change=paste(unique(wig$name), collapse="|"), 
                                                  WIG.type=paste(unique(wig$type), collapse="|"), 
                                                  WIG.consequence=paste(unique(wig$consequence), collapse="|"), 
                                                  WIG.reads=sum(unique(wig$SupportingReads))))
                            }
                        }))


# Merge and if no RNA support RNA.reads as 0.
combsamples          <- cbind(combsamples, WIGsupport) %>%
                            mutate(RNA.sample      = ifelse(is.na(RNA.sample), WIG.sample, RNA.sample),
                                   RNA.change      = ifelse(is.na(RNA.change), WIG.change, RNA.change),
                                   RNA.type        = ifelse(is.na(RNA.type), WIG.type, RNA.type),
                                   RNA.consequence = ifelse(is.na(RNA.consequence), WIG.consequence, RNA.consequence),
                                   RNA.reads       = ifelse(is.na(RNA.reads), WIG.reads, RNA.reads)) %>%
                            select(-WIG.sample, -WIG.change, -WIG.type, -WIG.consequence, -WIG.reads) %>%
                            mutate(RNA.reads       = ifelse(is.na(RNA.reads), 0, RNA.reads)) %>%
                            mutate(RNA.reads       = ifelse(is.na(RNA.reads) & (paste(patient, stime, ssite, sside, sep="_") %in% with(cancersamples, paste(patient, stime, ssite, sside, sep="_"))), 0, RNA.reads))

# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Group status per sample and gene
# Filter events outside coding sequence, if no fusion is detected. This is done after matching the fusions, because an upstream fusion may cause a fusion to exon2, the first splice-acceptor, skipping the coding sequence in exon1.
# If two opposite breaks are observed within a gene, these are interpreted as balanced translocation, and the number of affected and intact gene copies calculated separately.
# For simplicity, if multiple events are observed within the same gene, only the most relevants are retained and then pooled.
# Events are annotated non-pathogenic only in case of truncation after exon 2. All other events occurring on gene coding regions are considered potentially pathogenic.
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Prioritize different consequences of a single breakevent and keep only the most relevant
combsvid             <- combsamples %>%
                            filter(!(effect=="noncoding_duplication" & exonicBasePhase1<0 & exonicBasePhase2<0)) %>%          
                            filter(!(effect=="truncated_transcripts" & ((gene1==gene5prime & exonJoined1==0) | (gene2==gene5prime & exonJoined2==0)))) %>%                                                                       
                            # drop non-coding, some expressed fragments lost
                            filter(!(!is.na(gene1) & gene1 %in% c(gaingenes, lossgenes) & (is.na(gene2) | !(gene2 %in% c(gaingenes, lossgenes))) & (is.na(exonicBasePhase1) | exonicBasePhase1<0))) %>%
                            filter(!(!is.na(gene2) & gene2 %in% c(gaingenes, lossgenes) & (is.na(gene1) | !(gene1 %in% c(gaingenes, lossgenes))) & (is.na(exonicBasePhase2) | exonicBasePhase2<0))) %>%
                            mutate(sampevent=paste(sample, svId, sep="_")) %>%
                            arrange(sampevent) %>%
                            group_by(sampevent) %>%
                            mutate(nannots=n()) %>%
                            ungroup() %>%
                                             # include if only one annotation
                           mutate(priority= ifelse(nannots==1, 1,                                                                                                                                                               
                                             # prioritize expressed sense-sense fusions, gene-gene and intragenic events
                                             ifelse(!is.na(RNA.reads) & RNA.reads>0 & !is.na(gene5prime) & !is.na(gene3prime) & gene5prime==gene1 & !is.na(exonJoined1) & !is.na(exonJoined2) & exonJoined1!="antisense" & exonJoined2!="antisense", 1,        
                                             ifelse(!is.na(RNA.reads) & RNA.reads>0 & !is.na(gene5prime) & !is.na(gene3prime) & gene5prime==gene2 & !is.na(exonJoined1) & !is.na(exonJoined2) & exonJoined1!="antisense" & exonJoined2!="antisense", 1,
                                             # prioritize events causing full loss of gene copies
                                             ifelse(!is.na(gene1) & gene1 %in% lossgenes & undisruptedCopyNumber1<0.5,1,
                                             ifelse(!is.na(gene2) & gene2 %in% lossgenes & undisruptedCopyNumber2<0.5,1,
                                             # prioritize coding intragenic events
                                             ifelse(!is.na(gene1) & !is.na(gene2) & gene1==gene2 & !is.na(exonicBasePhase1) & !is.na(exonicBasePhase2) & exonicBasePhase1!=-1 & exonicBasePhase2!=-1, 2,                                                       
                                             # prioritize sense-sense fusions
                                             ifelse(effect=="fusion", 3,                                                                                                                                                                                       
                                             # prioritize simple intragenic indels
                                             ifelse(str_detect(effect, "simple"), 4,                                                                                                                                                                           
                                             # prioritize events truncating the 5' gene
                                             ifelse(str_detect(effect, "truncated"), 5, 6)))))))))) %>%                                                                                                                                                          
                           group_by(sampevent) %>%
                           filter(priority==min(priority)) %>%
                           mutate(nannots=n()) %>%
                           ungroup() %>%
                                             # Another round: include if only one annotation
                           mutate(priority = ifelse(nannots==1, 1,                                                                                                                                                                                                
                                             # prioritize events causing full loss of gene copies
                                             ifelse(!is.na(gene1) & gene1 %in% lossgenes & undisruptedCopyNumber1<0.5,1,
                                             ifelse(!is.na(gene2) & gene2 %in% lossgenes & undisruptedCopyNumber2<0.5,1,
                                             # Prioritize intragenic events, where both ends occur on the coding region
                                             ifelse(!is.na(gene1) & !is.na(gene2) & gene1==gene2 & !is.na(exonicBasePhase1) & !is.na(exonicBasePhase2) & exonicBasePhase1>=0 & exonicBasePhase2>=0, 2,
                                             # Prioritize gene-gene events, where both ends occur on a gene coding region
                                             ifelse(!is.na(exonicBasePhase1) & !is.na(exonicBasePhase2) & exonicBasePhase1>=0 & exonicBasePhase2>=0, 3,
                                             # Prioritize events, where the break occurs on coding region of the gene of interest
                                             ifelse(str_detect(affected_gene, "gene1") & !is.na(exonicBasePhase1) & exonicBasePhase1>=0 & gene1 %in% c(gaingenes, lossgenes), 4,
                                             ifelse(str_detect(affected_gene, "gene2") & !is.na(exonicBasePhase2) & exonicBasePhase2>=0 & gene2 %in% c(gaingenes, lossgenes), 4, 5)))))))) %>%
                           group_by(sampevent) %>%
                           mutate(allsame = length(unique(DNA_change))==1) %>%
                           filter(priority==min(priority)) %>%
                           mutate(nselect=length(which(priority==min(priority)))) %>%
                           filter(effect!="noncoding_duplication")

# Pool equally important consequences of a single breakevent
combbreaks           <- subset(combsvid, nselect==1) %>%
                            ungroup() %>%
                            mutate(across(setdiff(colnames(combsvid), c("RNA.reads", "undisruptedCopyNumber1", "affectedCopyNumber1", "undisruptedCopyNumber2", "affectedCopyNumber2")), ~ as.character(.x))) %>%
                            add_row(subset(combsvid, nselect %in% 2:6) %>%
                                    mutate(allnonc=all(exonJoined1=="antisense" | exonJoined1=="intergenic" | exonJoined2=="antisense" | exonJoined2=="intergenic")) %>%
                                    filter(!((!is.na(exonJoined1) & exonJoined1=="antisense") | (!is.na(exonJoined2) & exonJoined2=="antisense")) | allnonc) %>%
                                    summarise(across(setdiff(colnames(combsvid), c("RNA.reads", "undisruptedCopyNumber1", "affectedCopyNumber1", "undisruptedCopyNumber2", "affectedCopyNumber2", "sampevent")), ~ paste(unique(.x), collapse=";")),
                                              across(c("RNA.reads"), ~ ifelse(any(!is.na(.x)), max(.x, na.rm=T), NA)),
                                              across(c("undisruptedCopyNumber1", "undisruptedCopyNumber2"), ~ ifelse(any(!is.na(.x)), min(.x, na.rm=T), NA)),
                                              across(c("affectedCopyNumber1", "affectedCopyNumber2"), ~ ifelse(any(!is.na(.x)), max(.x, na.rm=T), NA)))) %>%
                            select(-nselect, -priority, -allsame, -sampevent, -nannots)


# Summarize to the level of gene-of-interest, row per sample-gene. Modify the undisruptedCopyNumber for balanced translocations.
brk_sample_gene <- combbreaks %>%
                            mutate(undisruptedCopyNumber1=ifelse(str_detect(DNA_change, "dup$"), undisruptedCopyNumber1-affectedCopyNumber1, undisruptedCopyNumber1),
                                   undisruptedCopyNumber2=ifelse(str_detect(DNA_change, "dup$"), undisruptedCopyNumber2-affectedCopyNumber2, undisruptedCopyNumber2)) %>%
                            mutate(affected_gene =        ifelse(effect %in% c("loss_transcription", "truncated_transcripts") & !is.na(gene1) & !is.na(gene2) & gene1==gene2, "gene1",
                                                          ifelse(effect %in% c("loss_transcription", "truncated_transcripts"), "gene1 gene2",
                                                          ifelse(effect=="promoter_hijacking" & !is.na(gene1) & !is.na(gene3prime) & gene1==gene3prime, "gene1",
                                                          ifelse(effect=="promoter_hijacking" & !is.na(gene2) & !is.na(gene3prime) & gene2==gene3prime, "gene2", affected_gene))))) %>%
                            mutate(primary_select =       ifelse(affected_gene=="gene1", "gene1",
                                                          ifelse(affected_gene=="gene2", "gene2",
                                                          ifelse(!is.na(gene1) & gene1 %in% lossgenes, "gene1",
                                                          ifelse(!is.na(gene2) & gene2 %in% lossgenes, "gene2",
                                                          ifelse(!is.na(gene1) & gene1=="ERBB2", "gene1",
                                                          ifelse(!is.na(gene2) & gene2=="ERBB2", "gene2",
                                                          ifelse(!is.na(gene1) & gene1 %in% gainprim, "gene1",
                                                          ifelse(!is.na(gene2) & gene2 %in% gainprim, "gene2", 
                                                          ifelse(!is.na(gene1) & gene1 %in% gaingenes, "gene1",
                                                          ifelse(!is.na(gene2) & gene2 %in% gaingenes, "gene2", 
                                                          ifelse(!is.na(undisruptedCopyNumber1) & undisruptedCopyNumber1<0.5, "gene1",
                                                          ifelse(!is.na(undisruptedCopyNumber2) & undisruptedCopyNumber2<0.5, "gene2", 
                                                          ifelse(!is.na(gene1), "gene1", "gene2")))))))))))))) %>%
                            mutate(primary_gene =         ifelse(primary_select=="gene1", gene1, gene2),
                                   primary_strand =       ifelse(primary_select=="gene1", strand1, strand2),
                                   primary_pos =          ifelse(primary_select=="gene1", pos1, pos2)) %>%
                            filter(primary_gene %in% c(gaingenes, lossgenes)) %>%
                            mutate(effect=                ifelse(effect %in% c("loss_transcription", "no_impact"), "fusion_from_unknown",
                                                          ifelse(effect=="promoter_deletion", "5prime_deletion", effect))) %>%
                            mutate(pathogenic=                   effect %in% c("fusion_from_unknown", "fusion", "5prime_deletion", "promoter_hijacking", "simple_deletion", "simple_duplication") |
                                                                (effect=="truncated_duplication" & !str_detect(DNA_change, "1_2dup")) |
                                                                (effect=="truncated_transcripts" & !str_detect(DNA_change, "e\\.2\\:\\:"))) %>%
                            mutate(groupkey =              paste(sample, primary_gene, sep="_")) %>%
                            arrange(groupkey) %>%
                            group_by(groupkey) %>%
                            mutate(nevents=                n(), 
                                   expressed=              any(!is.na(RNA.reads) & RNA.reads>0)) %>%
                            mutate(across(everything(), ~ ifelse(.x=="NA", NA, .x))) %>%
                            summarize(
                                   sample=                 unique(sample),
                                   patient=                unique(patient),
                                   stime=                  unique(stime),
                                   ssite=                  unique(ssite),
                                   sside=                  unique(sside),
                                   sord=                   unique(sord),
                                   svId=                   paste(unique(svId), collapse=";"),
                                   resolvedType=           paste(unique(resolvedType), collapse=";"),
                                   effect=                 ifelse(n()==1, unique(effect),
                                                           ifelse(n()==2 & any(str_detect(DNA_change, paste(primary_gene, ".+\\:\\:", sep=""))) & any(str_detect(DNA_change, paste("\\:\\:", primary_gene, sep=""))) & all(gene1!=gene2 | is.na(gene1) | is.na(gene2)), "balanced_translocation", "multiple")),
                                   frame=                  ifelse(any(!is.na(frame)), paste(na.omit(unique(frame)), collapse=";"), NA),
                                   undisruptedCopyNumber=  ifelse(n()==1 & all(!is.na(gene1)) & all(!is.na(gene2)) & all(gene1==gene2), min(undisruptedCopyNumber1, undisruptedCopyNumber2),
                                                           ifelse(n()==1 & all(!is.na(gene1)) & all(primary_gene==gene1), undisruptedCopyNumber1,
                                                           ifelse(n()==1, undisruptedCopyNumber2,
                                                           ifelse(n()>2, min(c(undisruptedCopyNumber1[which(gene1==primary_gene)], undisruptedCopyNumber2[which(gene2==primary_gene)])),
                                                           ifelse(n()==2 & !(any(str_detect(DNA_change, paste(primary_gene, ".+\\:\\:", sep=""))) & any(str_detect(DNA_change, paste("\\:\\:", primary_gene, sep=""))) & all(gene1!=gene2 | is.na(gene1) | is.na(gene2))), min(c(undisruptedCopyNumber1[which(gene1==primary_gene)], undisruptedCopyNumber2[which(gene2==primary_gene)])),
                                                           ifelse((all(primary_strand==1) & primary_pos[which(str_detect(DNA_change, paste(primary_gene, ".+\\:\\:", sep="")))] > primary_pos[which(str_detect(DNA_change, paste("\\:\\:", primary_gene, sep="")))]) |
                                                                  (all(primary_strand==-1)& primary_pos[which(str_detect(DNA_change, paste(primary_gene, ".+\\:\\:", sep="")))] < primary_pos[which(str_detect(DNA_change, paste("\\:\\:", primary_gene, sep="")))]),
                                                                   min(undisruptedCopyNumber1[which(gene1==primary_gene & !is.na(gene1))]-affectedCopyNumber1[which(gene1==primary_gene & !is.na(gene1))], undisruptedCopyNumber2[which(gene2==primary_gene & !is.na(gene2))]-affectedCopyNumber2[which(gene2==primary_gene & !is.na(gene2))]), 
                                                                   min(c(undisruptedCopyNumber1[which(gene1==primary_gene)], undisruptedCopyNumber2[which(gene2==primary_gene)])))))))),
                                   affectedCopyNumber=     ifelse(n()<3, max(c(affectedCopyNumber1[which(gene1==primary_gene)], affectedCopyNumber2[which(gene2==primary_gene)])), NA),
                                   primary_gene=           unique(primary_gene),
                                   DNA_change=             paste(unique(DNA_change), collapse=";"),
                                   mutation=               paste(unique(mutation), collapse=";"),
                                   pathogenic=             any(pathogenic),
                                   RNA.sample=             ifelse(any(!is.na(RNA.sample)), paste(na.omit(unique(RNA.sample)), collapse=";"), NA),
                                   RNA.change=             ifelse(any(!is.na(RNA.change)), paste(na.omit(unique(RNA.change)), collapse=";"), NA),
                                   RNA.type=               ifelse(any(!is.na(RNA.type)), paste(na.omit(unique(RNA.type)), collapse=";"), NA),
                                   RNA.consequence=        ifelse(any(!is.na(RNA.consequence)), paste(na.omit(unique(RNA.consequence)), collapse=";"), NA),
                                   RNA.reads=              ifelse(any(!is.na(RNA.reads)), sum(na.omit(RNA.reads)), NA),
                                   nevents=                n(), 
                                   expressed=              any(!is.na(RNA.reads) & RNA.reads>0))

# Write out results table
write.table(brk_sample_gene, file="breakevents_expressed.csv", col.names=T, row.names=F, quote=F, sep="\t")
                              


