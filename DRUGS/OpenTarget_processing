# OpenTargets processing script.
#
#Written by Andreas Hainari.

suppressPackageStartupMessages({
  library(magrittr)
  library(dplyr)
  library(tidyr)
  library(stringr)
})

# ---- Inputs ----
antimicrobiological <- read.csv("list of antimicrobiological drugs", stringsAsFactors = FALSE)
Phase_1_followup   <- read.csv("Follow up data for phase 1 medicine", stringsAsFactors = FALSE)
generic_chemotherapy_genes_2_0 <- read.csv("List of generic chemo genes", stringsAsFactors = FALSE)
gene_name_and_id <- read.csv("List of ENSG codes for all genes", stringsAsFactors = FALSE)

mechanism <- read.csv("Opentarget_mechanism_raw.csv", stringsAsFactors = FALSE)
drug      <- read.csv("Opentarget_drug_raw.csv", stringsAsFactors = FALSE)

# ---- Prepare ENSG -> gene mapping ----
colnames(gene_name_and_id) <- c("ENSG", "Gene")
ensg_map <- gene_name_and_id %>%
  mutate(ENSG = sub("\\..*$", "", ENSG))  # strip ENSG suffixes like ".1"

# ---- Merge mechanism and drug ----
mechanism_expanded <- mechanism %>%
  mutate(chemblIds = strsplit(as.character(chemblIds), ",")) %>%
  unnest(chemblIds) %>%
  mutate(chemblIds = trimws(chemblIds))

combined <- drug %>%
  left_join(mechanism_expanded, by = c("id" = "chemblIds"))

# ---- Build 'all_targets' and explode into rows ----
clean_and_split_targets <- function(linkedTargets, targets) {
  combined_text <- paste0(linkedTargets, ",", targets)
  combined_text <- str_replace_all(combined_text, 'c\\(|\\)|"', "")
  combined_text <- str_replace_all(combined_text, "\\s+", "")
  vals <- unlist(strsplit(combined_text, ","))
  vals <- trimws(vals)
  vals <- vals[vals != ""]
  vals
}

combined_targets <- combined %>%
  rowwise() %>%
  mutate(all_targets = list(clean_and_split_targets(linkedTargets, targets))) %>%
  ungroup() %>%
  unnest(all_targets) %>%
  mutate(all_targets = trimws(all_targets)) %>%
  filter(all_targets != "",
         !grepl("^[0-9]+$", all_targets),
         !grepl("list\\(", all_targets, ignore.case = TRUE))

# ---- Split cancer vs non-cancer ----
opentarget_CANCER_raw <- combined_targets %>%
  filter(grepl("MONDO_0004992", linkedDiseases, fixed = TRUE)) %>%
  select(id, name, drugType, actionType, all_targets,
         maximumClinicalTrialPhase, hasBeenWithdrawn,
         yearOfFirstApproval, isApproved)

# Drop "NA" target rows when other targets exist for the same chembl id
cancer_clean <- opentarget_CANCER_raw %>%
  group_by(id) %>%
  filter(!(all_targets == "NA" & any(all_targets != "NA"))) %>%
  ungroup() %>%
  distinct()

opentarget_NONCANCER_raw <- combined_targets %>%
  filter(!grepl("MONDO_0004992", linkedDiseases, fixed = TRUE)) %>%
  select(id, name, drugType, actionType, all_targets)

# ---- Process cancer dataset with clearer variable names ----
cancer_df <- cancer_clean

# Map ENSG to gene symbol and compute simple flags
mapped_df <- cancer_df %>%
  mutate(
    protein_target = ensg_map$Gene[match(all_targets, ensg_map$ENSG)],
    ApprovedAfter1989 = yearOfFirstApproval > 1989,
    noGenericChemo = !(protein_target %in% generic_chemotherapy_genes_2_0$genes),
    gene = protein_target
  )

# Termination status from Phase 1 follow-up
annotated_df <- mapped_df %>%
  mutate(
    terminated = ifelse(id %in% Phase_1_followup$ChemblId[Phase_1_followup$Terminated == "TRUE"],
                        "TRUE", "FALSE")
  )

# ---- Derive functional annotations from actionType ----
inhibitory_actions <- c("ANTAGONIST", "INHIBITOR", "NEGATIVE ALLOSTERIC MODULATOR", "CROSS-LINKING AGENT")
activating_actions  <- c("AGONIST", "PARTIAL AGONIST", "POSITIVE ALLOSTERIC MODULATOR", "OPENER", "POSITIVE MODULATOR", "ACTIVATOR")
flagged_actions     <- c("VACCINE ANTIGEN", "ANTISENSE INHIBITOR", "DISRUPTING AGENT", "EXOGENOUS PROTEIN", "BINDING AGENT", "NA", "OTHER", "STABILISER")

annotated_df <- annotated_df %>%
  mutate(
    `function` = case_when(
      actionType == "MODULATOR" ~ "modulates",
      actionType %in% inhibitory_actions ~ "whose hyperactivation is caused by",
      actionType == "BLOCKER" ~ "blocking",
      actionType %in% activating_actions ~ "stimulating",
      actionType %in% flagged_actions ~ "flagged out",
      TRUE ~ NA_character_
    ),
    potentialEffect = case_when(
      `function` %in% c("blocking", "stimulating") ~ FALSE,
      TRUE ~ TRUE
    ),
    GoF_LoF = case_when(
      actionType == "MODULATOR" ~ "ambiguous",
      actionType %in% c(inhibitory_actions, "BINDING AGENT", "CROSS-LINKING AGENT") ~ "gain-of-function",
      actionType == "BLOCKER" ~ "gain-of-function",
      actionType %in% activating_actions ~ "loss-of-function",
      actionType %in% flagged_actions ~ "flagged out",
      TRUE ~ NA_character_
    )
  )

# ---- Simple boolean flags and final renaming/selection ----
flagged_df <- annotated_df %>%
  mutate(
    rightDrugType = drugType %in% c("Antibody", "Small molecule"),
    havingProteinTarget = !is.na(protein_target),
    notAntimicrobiologicalDrug = !(name %in% antimicrobiological$drug)
  )

renamed_df <- flagged_df %>%
  rename(
    chemblId      = id,
    drug          = name,
    action        = actionType,
    targetENSG    = all_targets,
    proteinTarget = protein_target
  )

final_df <- renamed_df %>%
  select(
    chemblId,
    drug,
    drugType,
    action,
    targetENSG,
    proteinTarget,
    `function`,
    gene,
    GoF_LoF,
    maximumClinicalTrialPhase,
    hasBeenWithdrawn,
    yearOfFirstApproval,
    isApproved,
    ApprovedAfter1989,
    rightDrugType,
    havingProteinTarget,
    notAntimicrobiologicalDrug,
    noGenericChemo,
    terminated,
    potentialEffect
  )

# ---- Outputs ----
write.csv(opentarget_NONCANCER_raw, "your output file for non cancerous connections", row.names = FALSE)
write.csv(opentarget_CANCER_raw, "your output file for cancerous connections", row.names = FALSE)
write.csv(final_df, "your output for formatted OpenTargets table", row.names = FALSE)
