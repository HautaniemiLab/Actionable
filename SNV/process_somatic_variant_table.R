#!/usr/bin/env script

#Compute VAF vs purity of samples for select mutations of specific patients

library("arrow")
library("tidyverse")

args = commandArgs(trailingOnly=TRUE)

if (length(args)!=5) {
  stop("Provide five arguments \n Gene list file \n Path to a file listing ASEReadCounter output files \n Copy number database \n Path to metadata file \n Path to output directory \n", call.=FALSE)
} 

# --------------------------------------------------
# Functions
# --------------------------------------------------

#Combine two lists into one by pasting element-wise
pasteByList = function(a_list, b_list, sep=":") {
    lapply(
        1:length(a_list),
        function(i) {
            paste(a_list[[i]], b_list[[i]], sep=sep)
        }
    )
}

#Function for expected AF of somatic mutation given purity and tumour mutation and total copy number
expectedAF = function(mut_cn, total_cn, purity) {
    (mut_cn * purity) / (total_cn * purity + 2 * (1 - purity))
}


#Strip patient code prefix and DNA# suffix
trimSampleName = function(sample, patient) {
    sub(
        paste0("^", patient, "_"), "",
        sub(
            "_DNA\\d*$", "",
            sample
        )
    )
}

# --------------------------------------------------
# Inputs
# --------------------------------------------------

#Inputs
setwd(args[5])

mutation_table = read.delim("somatic.snv.coding.csv", header=T, sep="\t", as.is=T) %>%
                     filter(FILTER=="PASS")
gene_list      = read.delim(args[1], header=T, sep="\t", as.is=T)
cnv_table      = open_dataset(args[3]) |>
                     filter(Gene %in% gene_list$gene) |>
                     collect()
metadata_table = read.delim(args[4], header=T, sep="\t", as.is=T)
rnasereads     = read.delim(args[2], header=T, as.is=T, sep="\t") %>%
                             mutate(temp = str_extract(sample_RNA, "(?<=_)[poir][:alnum:]+")) %>%
                             mutate(stime = str_extract(temp, "^[poir][:digit:]?"),
                                    ssite = str_extract(temp, "(?<=[poir][:digit:]?)([:upper:][:lower:]{2}|LN)"),
                                    sside = str_extract(temp, "(?<=([:upper:][:lower:]{2}|LN))[LR]"),
                                    sord  = str_extract(temp, "[:digit:]$")) %>%
                             select(-temp) 


# --------------------------------------------------
# Processing with metadata and copy number data
# --------------------------------------------------

#Display mutations in desired format with AD mapped to AF and samples separated to rows
mutations_with_af_purity = mutation_table %>%
    select(
        # variant-calling data
        patient, CHROM, POS, REF, ALT, ID, readCounts, samples,
        # external annotations
        # MANE and refGene annotations
        Func.MANE, Gene.MANE, GeneDetail.MANE, ExonicFunc.MANE, AAChange.MANE, AAChange.refGene,
        # ClinVar annotations
        CLNSIG, CLNSIGCONF, CLNALLELEID,
        # Prediction algorithms (AM: AlphaMissense)
        CADD_phred, dbscSNV_ADA_SCORE, dbscSNV_RF_SCORE, 
        PolyPhenCat, PolyPhenVal, SIFTcat, SIFTval, AM_class, AM_score
    ) %>%
    mutate(Mutation = paste0(CHROM, ":", POS, REF, ">", ALT)) %>%
    mutate(AD = strsplit(readCounts, ";")) %>%
    mutate(sample = strsplit(samples, ";")) %>%
    mutate(sample.AD = pasteByList(sample, AD)) %>%
    unnest(sample.AD) %>%
    select(-sample, -AD, -readCounts, -samples) %>%
    separate(sample.AD, into=c("sample", "AD"), ":") %>%
    separate(AD, into=c("AD.0", "AD.1"), ",", convert=T) %>%
    mutate(DP = as.numeric(AD.0) + as.numeric(AD.1)) %>%
    mutate(AF = round(AD.1 / (AD.0 + AD.1), 2))

#Merge purities and round numbers
mutations_with_af_purity = mutations_with_af_purity %>%
    merge(metadata_table %>% filter(aberrant) %>% select(sample, purity)) %>%
    mutate(purity = round(purity, 2)) %>%
    mutate(afPurityRatio = AF / purity) %>%
    mutate(afPurityStr = paste(sprintf("%.2f", AF), sprintf("%.2f", purity), sep="/"))

#Add copy number
mutations_with_cnv = mutations_with_af_purity %>%
    mutate(key = paste(Gene.MANE, sample, "-")) %>%
    merge(
        cnv_table %>%
            mutate(key = paste(Gene, sample, "-")) %>%
            select(-Gene, -sample, -CNstatus, -ID),
        by="key", all.x=T
    ) %>%
    select(-key) %>%
    mutate(totalCN = nMajor + nMinor) %>%
    mutate(expHomAF = expectedAF(totalCN, totalCN, purity)) %>%
    mutate(expHomCI.lo = qbinom(0.025, DP, expHomAF)) %>%
    mutate(expHomCI.hi = qbinom(0.975, DP, expHomAF)) %>%
    mutate(expHomCI.cover = expHomCI.lo <= AD.1) %>% 
    mutate(expHom.pbinom.lower = pbinom(AD.1, DP, expHomAF))

#Extract and sort protein changes based on transcript order
#Determine potentially benign conflicting pathogenicity variants from the breakdown
#Finally modify the table for GenomeSpy viewing
mutations_out = mutations_with_cnv %>%
    mutate(HGVS.change = ifelse(AAChange.MANE==".", paste(Gene.MANE, GeneDetail.MANE, sep=":"), AAChange.MANE)) %>%
    mutate(HGVS.change.alt = AAChange.refGene) %>%
    select(
        sample, patient, CHROM, POS, REF, ALT, Gene=Gene.MANE, Mutation,
        HGVS.change, HGVS.change.alt,
        ClinVar=CLNSIG, ClinVar.conf=CLNSIGCONF, ClinVar.ID=CLNALLELEID,
        Func.MANE, ExonicFunc.MANE, rsID=ID, 
        CADD_phred,dbscSNV_ADA_SCORE, dbscSNV_RF_SCORE, 
        PolyPhenCat, PolyPhenVal, SIFTcat, SIFTval, AM_class, AM_score,
        AD.0, AD.1, DP, AF, purity, afPurityRatio, afPurityStr,
        nMajor, nMinor, totalCN, LOHstatus,
        expHomAF, expHomCI.lo, expHomCI.hi, expHomCI.cover, expHom.pbinom.lower
    ) %>%
    mutate(
        conflBenign = ifelse(
            is.na(ClinVar.conf),
            F,
            sapply(ClinVar.conf, function(x) sum(as.integer(strsplit(x, ",")[[1]][4:5])) == 0)
        )
    ) %>%
    mutate(ClinVar = ifelse(ClinVar == ".", NA, ClinVar)) %>%
    mutate(rsID = ifelse(rsID == ".", "", rsID)) %>%
    mutate(across(c("ExonicFunc.MANE", "Func.MANE", "ClinVar"), ~gsub("_", " ", .))) %>%
    mutate(Type = ifelse(ExonicFunc.MANE==".", Func.MANE, ExonicFunc.MANE)) %>%
    mutate(dbscSNV_ADA_SCORE = ifelse(dbscSNV_ADA_SCORE==".", "0", dbscSNV_ADA_SCORE),
           dbscSNV_RF_SCORE = ifelse(dbscSNV_RF_SCORE==".", "0", dbscSNV_RF_SCORE)) %>%           
    mutate(dbscSNV_ADA_SCORE = as.numeric(dbscSNV_ADA_SCORE),
           dbscSNV_RF_SCORE = as.numeric(dbscSNV_RF_SCORE)) %>%
    mutate(dbscSNVmax = ifelse(dbscSNV_ADA_SCORE > dbscSNV_RF_SCORE, dbscSNV_ADA_SCORE, dbscSNV_RF_SCORE)) %>%
    mutate(Type = ifelse(dbscSNVmax>0.6, "splicing", Type)) %>%
    mutate(ClinVar.conf = gsub(",", ", ", ClinVar.conf)) %>%
    mutate(ClinVar = gsub("/", " / ", ClinVar)) %>%
    mutate(ClinVar = ifelse(grepl("^Conflicting", ClinVar), paste0(ClinVar, " (", ClinVar.conf, ")"), ClinVar)) %>%
    arrange(patient, Gene, Mutation) %>%
    mutate(index = row_number()) %>%
    select(index, everything())


# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ADD data on allele-specific expression


mutations_expresssed = do.call(rbind, lapply(1:nrow(mutations_out), function(i) {

                             temp = mutations_out[i,] %>%
                                        mutate(temp = str_extract(sample, "(?<=_)[poir][:alnum:]+")) %>%
                                        mutate(stime = str_extract(temp, "^[poir][:digit:]?"),
                                               ssite = str_extract(temp, "(?<=[poir][:digit:]?)([:upper:][:lower:]{2}|LN)"),
                                               sside = str_extract(temp, "(?<=([:upper:][:lower:]{2}|LN))[LR]"),
                                               sord  = str_extract(temp, "[:digit:]$")) %>%
                                        select(-temp) 

                             if(with(temp, nchar(REF)==1 & nchar(ALT)==1 & is.na(refCount))){
                                 if(any(with(temp, paste(patient, stime, ssite, sside, sord, sep="_")) %in% with(rnasereads, , paste(patient, stime, ssite, sside, sord, sep="_")))) {
                                     tempfile <- unlist(subset(rnasereads, paste(patient, stime, ssite, sside, sord, sep="_") == with(temp, paste(patient, stime, ssite, sside, sord, sep="_")), select="file"))[1]
                                     tempdata <- read.delim(tempfile, header=T, as.is=T, sep="\t") %>%
                                                     filter(contig==temp$CHROM & position==temp$POS) 
                                     if(nrow(tempdata==0)) {
                                        temp <- temp %>%
                                                    mutate(refCount = 0,
                                                           altCount = 0)
                                        return(temp)
                                     } else {
                                        temp <- temp %>%
                                                    mutate(refCount = tempdata$refCount[1],
                                                           altCount = tempdata$altCount[1])
                                        return(temp)        
                                     }                                   
                                 } else {
                                        return(temp)
                                 }
                             } else{
                                    return(temp)
                             }
                         }))
                             
# -------------------------------------------------------------------------------------------------
# CGI predictions were retrieved via the CGI web interface, using a variant matrix as an input file
# -------------------------------------------------------------------------------------------------

dataToCGI   <- mutations_expressed %>%
                   select(chr=CHROM, pos=POS, ref=REF, alt=ALT) %>%
                   distinct()

write.table(dataToCGI, file="dataToCGI.csv", sep="\t", col.names=T, row.names=F, quote=F)
write.table(mutations_expressed, file="somatic.expressed.csv", col.names=T, row.names=F, quote=F, sep="\t")
