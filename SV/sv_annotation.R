######################################## SV annotations ########################################
library(GenomicRanges)
library(tidyverse)

source("sv_annotation_functions.R")

# List samples and patients
# Load a dataframe with sample and patient information
samples <- read.table("purity_ploidy_estimates.tsv", sep="\t", header=T) |>
  select(sample, patient)  
# Load copy-number segmentation
segs <- read.table("segmentation_info_final.tsv", sep="\t", header=T)

# Extract annotation
lst <- apply(samples, 1, function(x) {
  process_sample(patient = x[["patient"]], sample = x[["sample"]], cnv = segs)
})
saveRDS(lst, "pre_sv_annotation.rds")

sv_annotation <- lst |>
  bind_rows() |>
  mutate(affectedCopyNumber1 = copyNumber1 - undisruptedCopyNumber1,
         affectedCopyNumber2 = copyNumber2 - undisruptedCopyNumber2) |>
  #select(sample, patient, clusterId, svId, resolvedType, Type, 4:14, 16:26)
  select(sample, patient, clusterId, svId, resolvedType, Type, 4:9, affectedCopyNumber1, 10:14, 16:21, affectedCopyNumber2, 22:26) |>
  distinct()
  
write.table(sv_annotation, "pre_sv_annotation.tsv", sep="\t", col.names = T, row.names = F)

# Add SV effect on the transcript 
cases_effect <- read.table(file.path("sv_annotation_effects.tsv"), sep="\t", header=T)

ann_vars <- sv_annotation |>
  mutate(Type_simple = case_when(
    Type %in% c("DEL", "DUP", "INV") ~ Type,
    orientation1 == 1 & orientation2 == -1 ~ "DEL",   # Transform BND cases in known cases
    orientation1 == -1 & orientation2 == 1 ~ "DUP",
    orientation1 == orientation2 ~ "INV",
    .default = Type),
    # create fake genes to be able to join the cases_effect table
    Gene1 = if_else(is.na(gene1) | is.na(strand1),
                         NA_character_, "gene1"),
    # Only label Gene2 when BOTH gene2 and strand2 are present
    Gene2 = case_when(
      is.na(gene2) | is.na(strand2) ~ NA_character_,
       gene1 == gene2                 ~ "gene1",
       TRUE                           ~ "gene2"
     )) |>
  left_join(cases_effect |> 
              filter( !Type_simple %in% c("DEL", "DUP") ), by = colnames(cases_effect)[-8], relationship = "many-to-many") |>
  left_join(cases_effect |> 
              filter( Type_simple %in% c("DEL", "DUP") ) |>
              select( -orientation1, -orientation2 ), by = colnames(cases_effect)[c(1:3, 5:6)], relationship = "many-to-many") |>
  mutate(effect = coalesce(effect.x, effect.y)) |>
  select(-effect.x, -effect.y)

# Get annotations
columns_to_swap <- c("gene", "transcriptId", "regionType",
  "strand", "orientation", "undisruptedCopyNumber", "affectedCopyNumber",
  "exonicBasePhase", "chromosome", "pos", "exonJoined")
final_ann <- get_variants(ann_vars, fusion_bases = columns_to_swap)
write.table(final_ann, "sv_annotation.tsv", sep="\t", col.names = T, row.names = F)

# Further filtering for effect preference
low <- c("loss_transcription", "truncated_transcripts", "noncoding_duplication", "promoter_deletion")
ann_filtered <- final_ann %>%
  filter(effect != "no_impact") %>%
  group_by(patient, svId) %>%
  filter(!(effect %in% low & any(!effect %in% low))) %>%
  filter(
    if (all(effect %in% low)) {
      !(replace_na(regionType1 == "UPSTREAM", FALSE) |
          replace_na(regionType2 == "UPSTREAM", FALSE))
    } else TRUE
  ) %>%
  ungroup()
write.table(ann_filtered, "customized_sv_annotation.tsv", sep="\t", col.names = T, row.names = F)

















