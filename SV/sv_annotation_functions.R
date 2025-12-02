#' Annotate breakends for a set of SV ids
#'
#' @param be        data.frame of LINX breakends (columns like: svId, isStart, gene, transcriptId, disruptive, ...).
#' @param sv_vis    data.frame with SV metadata (must contain SvId, PosStart, PosEnd, Type).
#' @param svIds     vector of svId values to process.
#' @param sv        data.frame of LINX svs
#' @param clusters  data.frame of LINX clusters
#' @param protein_coding_only logical; keep only biotype == "protein_coding" before pairing (default TRUE).
#' @param keep_cols_range logical; if TRUE, keeps the final select(svId, resolvedType, Type, 2:16); if FALSE, keep all columns
#' @return tibble with paired start/end rows (full cartesian when both sides exist; keeps singletons if only one side).
annotate_breakends <- function(be, sv_vis, svIds, sv, clusters,
                                   protein_coding_only = TRUE,
                                   keep_cols_range = TRUE) {
  
  process_one_sv <- function(s) {
    sv_vis_s <- sv_vis %>% filter(SvId == s)
    if (nrow(sv_vis_s) == 0) return(tibble())
    
    sv_genes <- sv |>
      filter(svId == s) |>
      select(geneStart, geneEnd) |> 
      mutate(genes = paste(geneStart, geneEnd, sep =";")) |> 
      pull(genes)
    sv_genes <- unlist(str_split(sv_genes, ";"))
    
    ex <- be %>%
      filter(svId == s) |>
      filter(gene %in% sv_genes) %>%
      { if (protein_coding_only && "biotype" %in% names(.)) filter(., biotype == "protein_coding") else . } %>%
      # keep key cols explicitly (safer than numeric indices)
      select(id, svId, isStart, gene, transcriptId,
             disruptive, regionType, strand, orientation, 
             undisruptedCopyNumber, exonUp, exonDown, exonicBasePhase, chromosome) %>%
      mutate(
        pos = ifelse(isStart == "true", sv_vis_s$PosStart, sv_vis_s$PosEnd),
        .is_start = if (is.logical(isStart)) isStart else isStart == "true"#,
        #.disr     = if (is.logical(disruptive)) disruptive else disruptive == "true"
      )
    
    # skip if no disruptive at all
    #if (!any(ex$.disr)) return(tibble())
    
    # apply your rules for keeping rows
    # TODO: modify this part. It is complicated because I removed the condition "disruptive"
    ex <- ex %>%
      mutate(
        has_true_start = any(.is_start), 
        has_true_end   = any( !.is_start)
      ) %>%
      filter(
        (has_true_start & has_true_end) |
          (has_true_start & !has_true_end) |
          (!has_true_start & has_true_end)
      ) %>%
      select(-.is_start,  -has_true_start, -has_true_end) 
    
    # starts / ends with suffixes; keep svId for each side
    start_be <- ex %>%
      filter(isStart == "true") %>%
      select(-isStart, -disruptive) %>%
      rename_with(~ paste0(.x, "1")) 
    
    end_be <- ex %>%
      filter(isStart == "false") %>%
      select(-isStart, -disruptive) %>%
      rename_with(~ paste0(.x, "2"))   
    
    # combine (cartesian if both sides; keep singletons otherwise)
    pairs <- if (nrow(start_be) == 0 && nrow(end_be) == 0) {
              tibble()
            } else {
              full_join(
                mutate(start_be, .k = 1),
                mutate(end_be,   .k = 1),
                by = ".k"
              ) %>%
                select(-.k)
            } |>
      mutate(svId = ifelse(is.na(svId1), svId2, svId1)) |>
      left_join(
        sv_vis_s[, c("SvId","Type","ChrStart","PosStart","ChrEnd","PosEnd", "OrientStart", "OrientEnd")],
        by = c("svId" = "SvId")
      ) %>%
      transmute( # Force typeof elements in the case of NAs
        svId,
        id1,
        gene1 = as.character(gene1),
        transcriptId1 = as.character(transcriptId1),
        regionType1 = as.character(regionType1),
        strand1 = as.numeric(strand1),
        orientation1 = as.numeric(ifelse(is.na(orientation1), OrientStart, orientation1)),
        undisruptedCopyNumber1 = as.numeric(undisruptedCopyNumber1),
        exonUp1 = as.numeric(exonUp1),
        exonDown1 = as.numeric(exonDown1),
        exonicBasePhase1 = as.numeric(exonicBasePhase1),
        chromosome1 = as.character(ifelse(is.na(gene1), ChrStart, chromosome1)),
        pos1        = ifelse(is.na(gene1), PosStart, pos1),
        id2,
        gene2 = as.character(gene2),
        transcriptId2 = as.character(transcriptId2),
        regionType2 = as.character(regionType2),
        strand2 = as.numeric(strand2),
        orientation2 = as.numeric(ifelse(is.na(orientation2), OrientEnd, orientation2)),
        undisruptedCopyNumber2 = as.numeric(undisruptedCopyNumber2),
        exonUp2 = as.numeric(exonUp2),
        exonDown2 = as.numeric(exonDown2),
        exonicBasePhase2 = as.numeric(exonicBasePhase2),
        chromosome2 = as.character(ifelse(is.na(gene2), ChrEnd, chromosome2)),
        pos2        = ifelse(is.na(gene2), PosEnd,   pos2),
        Type        = as.character(Type)
        ) %>%
      filter(
        Type != "SGL",
        !(replace_na(regionType1 == "UPSTREAM", FALSE) &
            replace_na(regionType2 == "UPSTREAM", FALSE))
      ) %>%
      filter(
        if (any(!is.na(gene1) & !is.na(gene2) & gene1 == gene2 & (regionType1 != "UPSTREAM" & regionType2 != "UPSTREAM"))) {
          !is.na(gene1) & !is.na(gene2) & gene1 == gene2 & (regionType1 != "UPSTREAM" & regionType2 != "UPSTREAM")
        } else TRUE
      )
    
    # Eliminate intron-intron duplication-deletions
    if (any(
      replace_na(pairs$regionType1 == "INTRONIC", FALSE) &
      replace_na(pairs$regionType2 == "INTRONIC", FALSE) &
      (
        (replace_na(pairs$exonUp1 == pairs$exonUp2, FALSE) &
         replace_na(pairs$exonDown1 == pairs$exonDown2, FALSE)) |
        (replace_na(pairs$exonUp1 == pairs$exonDown2, FALSE) &
         replace_na(pairs$exonDown1 == pairs$exonUp2, FALSE))
      ) & replace_na(pairs$orientation1 != pairs$orientation2, FALSE) & # keep intra-intronic inversions 
      replace_na(pairs$gene1 == pairs$gene2, FALSE)
    )) {
      return(pairs[0,])  # return empty tibble if condition met
    } else { 
      return(pairs) 
    }
  }
  # run over all svIds and bind rows
  ann <- svIds %>% map(process_one_sv) %>% list_rbind() %>% as_tibble()
  
  if (nrow(ann) == 0) return(ann)
  
  ann <- ann %>%
    left_join(sv[, c("vcfId", "svId", "clusterId")], by = "svId") %>%
    left_join(clusters[, c("clusterId", "resolvedType")],   by = "clusterId") %>%
    select(-svId) %>%
    dplyr::rename(svId = vcfId) %>%
    filter(!grepl("purple|unbalanced", svId))
  
  # final select: mimic your `select(svId, resolvedType, Type, 2:16)`
  if (keep_cols_range) {
    # defensively select by names if available
    base_cols <- c("svId", "clusterId", "resolvedType", "Type")
    other_cols <- setdiff(names(ann), base_cols)
    # match the numeric slice 2:16 if it exists, else keep all others
    idx <- seq_len(ncol(ann))
    take_idx <- intersect(2:17, idx)
    if (length(take_idx) > 0) {
      ann <- ann %>% select(all_of(base_cols), all_of(names(ann)[take_idx]))
    } else {
      ann <- ann %>% select(all_of(base_cols), all_of(other_cols))
    }
  }
  
  as.data.frame(ann)
}

#' Annotate structural variant breakends with nearest CNV copy number
#'
#' Given a table of structural variant breakends (`ann`) and
#' a segmentation file of copy-number segments, this function finds
#' the nearest CNV junction to each breakpoint, selects the copy number
#' from the segment *before* or *after* the junction based on the
#' breakpoint orientation, and merges the results back into the
#' original SV data in wide format.
#'
#' @param dev_sample A data frame containing structural variant breakends.
#'   Must include columns:
#'   \describe{
#'     \item{id1, id2}{Identifiers for paired breakends (used to form `id`)}
#'     \item{chromosome1, pos1, orientation1}{Chromosome, position, and orientation of first breakend}
#'     \item{chromosome2, pos2, orientation2}{Chromosome, position, and orientation of second breakend}
#'     \item{sample}{Sample identifier (used to filter segmentation data)}
#'   }
#' @param segs_path Path to the segmentation file (e.g. `"segmentation_info_final.tsv"`).
#'   The file must contain columns: `sample`, `chromosome`, `start`, `end`, `copyNumber`.
#'
#' @return A data frame identical to `ann`, with 2 extra columns: {copyNumber1}{Copy number near the first breakend}
#'   and {copyNumber2}{Copy number near the second breakend}
#'   Breakends for which no CNV junction is found are assigned `NA`.

annotate_cnv <- function(ann_df,
                         samp,
                         cnv_file) {
  
  # -------- 1) Load CNV segments for this sample --------
  segs <- cnv_file %>%
    filter(sample == samp) %>%
    select(chromosome, start, end, copyNumber)
  
  # -------- 2) SV breakends (long) --------
  sv_be <- as.data.frame(ann_df) %>%
    transmute(
      id = paste(id1, id2, sep = "_"),
      chromosome1, pos1, orientation1,
      chromosome2, pos2, orientation2
    ) %>%
    pivot_longer(
      cols = -id,
      names_to      = c(".value", "end"),
      names_pattern = "^(chromosome|pos|orientation)([12])$"
    ) %>%
    distinct() %>%
    arrange(id, end)
  
  # -------- 3) Build CNV junctions with CN before/after --------
  cnv_junc <- segs %>%
    arrange(chromosome, start, end) %>%
    group_by(chromosome) %>%
    mutate(
      junc_pos  = end,
      cn_before = copyNumber,
      cn_after  = lead(copyNumber)
    ) %>%
    filter(!is.na(cn_after)) %>%
    ungroup() %>%
    select(chromosome, junc_pos, cn_before, cn_after)
  
  if (nrow(cnv_junc) == 0L) {
    warning("No CNV junctions found for sample ", samp, ". Returning NA copy numbers.")
    return(
      ann_df %>%
        mutate(copyNumber1 = NA_real_, copyNumber2 = NA_real_)
    )
  }
  
  # -------- 4) Nearest junction per breakend (GRanges) --------
  j_gr <- GenomicRanges::GRanges(
    seqnames = cnv_junc$chromosome,
    ranges   = IRanges(start = cnv_junc$junc_pos, end = cnv_junc$junc_pos),
    cn_before = cnv_junc$cn_before,
    cn_after  = cnv_junc$cn_after
  )
  
  sv_gr <- GRanges(
    seqnames = sv_be$chromosome,
    ranges   = IRanges(start = sv_be$pos, end = sv_be$pos)
  )
  
  # Nearest junction for each SV breakpoint
  hit_idx <- nearest(sv_gr, j_gr, ignore.strand = TRUE)  # integer index into j_gr (or NA)
  
  # Pull CN based on orientation: +1 -> before, -1 -> after
  cn_before_vec <- mcols(j_gr)$cn_before
  cn_after_vec  <- mcols(j_gr)$cn_after
  
  picked_cn <- ifelse(
    sv_be$orientation > 0,          # +1
    cn_before_vec[hit_idx],
    cn_after_vec[hit_idx]           # -1
  )
  
  sv_be_with_cn <- sv_be %>%
    mutate(copyNumber = picked_cn) %>%
    mutate(end = paste0(end)) %>%   # ensure end is a character key (1 / 2)
    pivot_wider(
      id_cols = id,
      names_from = end,
      values_from = c(chromosome, pos, orientation, copyNumber),
      names_glue = "{.value}{end}"
    ) %>%
    select(id,
           chromosome1, pos1, orientation1, copyNumber1,
           chromosome2, pos2, orientation2, copyNumber2)
  
  # -------- 5) Join back to the original annotation --------
  out <- ann_df %>%
    left_join(sv_be_with_cn)
  
  return(out)
}

#' Process and annotate structural variant data for a single sample
#'
#' This function reads LINX output tables for a given patient and sample, 
#' including structural variant (`.linx.svs.tsv`), breakend (`.linx.breakend.tsv`), 
#' cluster (`.linx.clusters.tsv`), and visualization (`.linx.vis_sv_data.tsv`) files. 
#' It then applies the [`annotate_breakends_all()`] function to generate 
#' annotated structural variant pairs, filtering for protein-coding breakends 
#' and merging start–end breakpoints as defined by the LINX data model.
#'
#' @param patient Character string specifying the patient identifier (corresponding to the directory containing LINX output files).
#' @param sample Character string specifying the sample identifier (used as the file prefix for LINX output tables).
#' @param base_dir Character string specifying the base directory containing patient folders. Defaults (`"."`).
#'
#' @return A tibble containing annotated SV pairs for the sample.
process_sample <- function(patient, sample, cnv, base_dir = ".") {
  # Build paths
  path_sv_vis <- file.path(base_dir, patient, paste0(sample, ".linx.vis_sv_data.tsv"))
  path_be     <- file.path(base_dir, patient, paste0(sample, ".linx.breakend.tsv"))
  path_cl     <- file.path(base_dir, patient, paste0(sample, ".linx.clusters.tsv"))
  path_sv     <- file.path(base_dir, patient, paste0(sample, ".linx.svs.tsv"))
  path_cnv    <- file.path(base_dir, "segmentation_info_final.tsv")
  
  # Fast & robust reads (fail nicely if a file is missing)
  if (!file.exists(path_sv_vis) || !file.exists(path_be) ||
      !file.exists(path_cl)     || !file.exists(path_sv)) {
    warning(sprintf("Missing LINX files for %s / %s — skipping", patient, sample))
    return(invisible(NULL))
  }
  
  sv_vis <- read.table(path_sv_vis, sep="\t", header=T)
  be     <- read.table(path_be, sep="\t", header=T)
  cl     <- read.table(path_cl, sep="\t", header=T)
  sv     <- read.table(path_sv, sep="\t", header=T)
  
  # Skip if no svId in breakends
  if (!"svId" %in% names(be) || nrow(be) == 0) return(tibble())
  
  svIds <- unique(be$svId)
  
  ann <- annotate_breakends(
    be      = be,
    sv_vis  = sv_vis,
    svIds   = svIds,
    sv      = sv,
    cl      = cl,
    protein_coding_only = TRUE,
    keep_cols_range     = FALSE
  )
  
  if (nrow(ann) == 0) return(ann)
  
  # annotate copy number
  ann_cnv <- annotate_cnv(ann, sample, cnv)
  
  # tag with patient/sample for downstream convenience
  ann_cnv %>% mutate(sample = sample, patient = patient, .before = 1)
}

#' Swap paired columns in a data frame
#' Internal helper function that swaps any number of column pairs with suffixes
#' `1` and `2` (e.g. `gene1` / `gene2`) when a logical swap flag column is `TRUE`.
#'
#' @param .df A data frame or tibble.
#' @param bases A character vector of base column names (e.g. `"gene"`, `"chromosome"`, `"pos"`).
#'   For each base name, columns `<base>1` and `<base>2` must exist in `.df`.
#' @param swap_flag_col Name of the logical column in `.df` indicating which rows to swap.
#'
#' @return A modified version of `.df` where, for each base name, the corresponding `<base>1`
#'  and `<base>2` columns have been swapped wherever `swap_flag_col` is `TRUE`.
swap_pairs <- function(.df, bases, swap_flag_col = "swap") {
  cols1 <- paste0(bases, "1")
  cols2 <- paste0(bases, "2")
  
  # build programmatic assignments for both sides
  lhs_updates <-
    setNames(
      lapply(seq_along(cols1), function(k)
        if_else(.df[[swap_flag_col]], .df[[cols2[k]]], .df[[cols1[k]]])
      ),
      cols1
    )
  
  rhs_updates <-
    setNames(
      lapply(seq_along(cols2), function(k)
        if_else(.df[[swap_flag_col]], .df[[cols1[k]]], .df[[cols2[k]]])
      ),
      cols2
    )
  
  mutate(.df, !!!lhs_updates, !!!rhs_updates)
}

#' Get chromosome rank to compare the position in the genome between two chromosomes
chr_rank <- function(x) {
  # Strip "chr" if present
  x <- sub("^chr", "", x, ignore.case = TRUE)
  
  # Map X/Y/MT/M
  out <- suppressWarnings(as.integer(x))
  out[is.na(out) & toupper(x) == "X"]  <- 23L
  out[is.na(out) & toupper(x) %in% c("M", "MT")] <- 25L
  out[is.na(out)] <- Inf
  out
}

#' Generate fusion 
#' Filters annotated structural variant data for fusion events and  create gene annotations.
#' It swaps gene1/gene2 (and related metadata) when the transcription affects antisense genes 
#' or in specific cases of duplication and deletion, so according to strand and orientation logic.
#'
#' @param df A data frame or tibble containing annotated SV calls.
#'   Must include columns such as `effect`, `Type_simple`, `strand1`, `orientation1`,
#'   and the paired columns `<base>1` / `<base>2` for each entry in `bases`.
#' @param bases Character vector of base names for paired columns to be swapped.
#'   Default includes:
#'   `"chromosome"`, `"pos"`, `"exonJoined"`.
#'
#' @return A tibble filtered to fusion events, with gene1/gene2 (and associated fields)
#'  swapped as needed. Includes derived columns `DNA_change` and `mutation`.
get_fusions <- function(
    df,
    bases = c("chromosome", "pos", "exonJoined")) {
  
  df %>%
    dplyr::filter(effect == "fusion") %>%
    dplyr::mutate(
      exonJoined1.tmp = dplyr::case_when(
        Type_simple == "DEL" & strand1 ==  1 ~ "Up",
        Type_simple == "DEL" & strand1 == -1 ~ "Down",
        Type_simple == "DUP" & strand1 ==  1 ~ "Down",
        Type_simple == "DUP" & strand1 == -1 ~ "Up",
        Type_simple == "INV" & strand1 == -1 & orientation1 ==  1 ~ "Down",
        Type_simple == "INV" & strand1 == -1 & orientation1 == -1 ~ "Up",
        Type_simple == "INV" & strand1 ==  1 & orientation1 ==  1 ~ "Up",
        Type_simple == "INV" & strand1 ==  1 & orientation1 == -1 ~ "Down"
      ),
      # swap when "Down"
      swap = exonJoined1.tmp == "Down",
      
      # choose exons first (pre-swap)
      exonJoined1 = dplyr::if_else(exonJoined1.tmp == "Up",
                                   as.character(exonUp1),
                                   as.character(exonDown1)),
      exonJoined2 = dplyr::if_else(exonJoined1.tmp == "Up",
                                   as.character(exonDown2),
                                   as.character(exonUp2))
    ) %>%
    { swap_pairs(., bases = bases, swap_flag_col = "swap") } %>%
    dplyr::mutate(
      affected_gene = "gene1 gene2",
      frame = case_when(
        is.na(gene1) | is.na(gene2) ~ NA,
        exonicBasePhase1 < 0 ~ NA,
        exonicBasePhase2 < 0 ~ NA,
        exonicBasePhase1 == exonicBasePhase2 ~ "in_frame",
        exonicBasePhase1 != exonicBasePhase2 ~ "frameshift"
      ),
      DNA_change = paste0(gene1, ":", transcriptId1, ":e.", exonJoined1, "::",
                          gene2, ":", transcriptId2, ":e.", exonJoined2),
      chr1r = chr_rank(chromosome1),
      chr2r = chr_rank(chromosome2),
      swap = chr2r < chr1r | (chr2r == chr1r & pos2 < pos1),
      mutation = if_else(
        swap,
        paste0(chromosome2, ":", pos2, "::", chromosome1, ":", pos1),
        paste0(chromosome1, ":", pos1, "::", chromosome2, ":", pos2)
      )
    ) %>%
    dplyr::select(-exonJoined1.tmp, -swap, -chr1r, -chr2r, -exonUp1, -exonUp2, -exonDown1, -exonDown2)
}

#' Generate simple deletion and duplications annotations
#' Filters annotated structural variant data for simple deletions and duplications 
#' and  create gene annotations.
#'
#' @param df A data frame or tibble containing annotated SV calls.
#'   Must include columns such as `effect`, `Type_simple`, `strand1`, `orientation1`,
#'   and the paired columns `<base>1` / `<base>2` for each entry in `bases`.
#'
#' @return A tibble filtered to simple deletions and duplications events.
#'  Includes derived columns `DNA_change` and `mutation`.
get_simple_deldup <- function(
    df,
    bases = c("chromosome", "pos")) {
  
  df |>
    filter(effect %in% c("simple_deletion", "simple_duplication")) |>
    mutate(swap = pos2 < pos1 ) %>%
    { swap_pairs(., bases = bases[bases!= "exonJoined"], swap_flag_col = "swap") } %>%
    mutate(
      exonUp1 = case_when(
        regionType1 == "INTRONIC" & exonUp1 == 0 & exonDown1==2 ~ 1,
        TRUE ~ exonUp1
      ),
      exonDown1 = case_when(
        regionType1 == "UPSTREAM" & exonUp1 == 0 & exonDown1==2 ~ 1,
        TRUE ~ exonDown1
      ),
      exonUp2 = case_when(
        regionType2 == "INTRONIC" & exonUp2 == 0 & exonDown2 == 2 ~ 1,
        TRUE ~ exonUp2
      ),
      exonDown2 = case_when(
        regionType2 == "UPSTREAM" & exonUp2 == 0 & exonDown2 == 2 ~ 1,
        TRUE ~ exonDown2
      )
    ) %>%
    mutate(
      exonJoined1.tmp = case_when(
        strand1 == 1 ~ "Down",
        strand1 == -1 ~ "Up"
      ),
      exonJoined1 = ifelse(exonJoined1.tmp == "Up", exonUp1, exonDown1),
      exonJoined2 = ifelse(exonJoined1.tmp == "Up", exonDown2, exonUp2),
      affected_gene = "gene1",
      frame = case_when(
        is.na(gene1) | is.na(gene2) ~ NA,
        exonicBasePhase1 < 0 ~ NA,
        exonicBasePhase2 < 0 ~ NA,
        exonicBasePhase1 == exonicBasePhase2 ~ "in_frame",
        exonicBasePhase1 != exonicBasePhase2 ~ "frameshift"
      ),
      DNA_change = case_when(
        exonJoined1 <= exonJoined2 ~ paste0(gene1, ":", transcriptId1, ":e.", exonJoined1, ifelse(exonJoined1 == exonJoined2, "", paste0("_", exonJoined2)), tolower(Type)),
        exonJoined1 > exonJoined2  ~ paste0(gene1, ":", transcriptId1, ":e.", exonJoined2, ifelse(exonJoined1 == exonJoined2, "", paste0("_", exonJoined1)), tolower(Type))
      ),
      mutation = paste0(chromosome1, ":", pos1, "_", pos2, tolower(Type)),
      exonJoined1 = as.character(exonJoined1),
      exonJoined2 = as.character(exonJoined2)) |>
    select(-exonJoined1.tmp) |>
    relocate(affected_gene, frame, DNA_change, mutation, .after = last_col())
}

#' Generate annotations for truncated transcripts
#' Filters annotated structural variant data for truncations
#' and  create gene annotations.
#'
#' @param df A data frame or tibble containing annotated SV calls.
#'   Must include columns such as `effect`, `Type_simple`, `strand1`, `orientation1`,
#'   and the paired columns `<base>1` / `<base>2` for each entry in `bases`.
#'
#' @return A tibble filtered to truncations events.
#'  Includes derived columns `DNA_change` and `mutation`.
get_truncations <- function(df) { 
  
  base <- df %>%
    filter(effect == "truncated_transcripts") %>%
    mutate(
      exonJoined1 = if_else(!is.na(gene1), as.character(exonUp1), "intergenic"),
      exonJoined2 = if_else(!is.na(gene2), as.character(exonUp2), "intergenic"),
      chr1r = chr_rank(chromosome1),
      chr2r = chr_rank(chromosome2),
      swap = chr2r < chr1r | (chr2r == chr1r & pos2 < pos1),
      mutation = if_else(
        swap,
        paste0(chromosome2, ":", pos2, "::", chromosome1, ":", pos1),
        paste0(chromosome1, ":", pos1, "::", chromosome2, ":", pos2)
      )
    ) %>%
    select(-chr1r, -chr2r, -swap)
  
  # rows with both genes present → produce left/right
  both <- base %>% 
    filter(!is.na(gene1) & !is.na(gene2))
  
  left <- both %>%
    mutate(
      exonJoined2   = "antisense",
      affected_gene = "gene1"
    )
  
  right <- both %>%
    mutate(
      exonJoined1   = "antisense",
      affected_gene = "gene2"
    )
  
  # rows where one side is intergenic → keep as single row, mark affected_gene
  intergenic <- base %>%
    filter(is.na(gene1) | is.na(gene2)) %>%
    mutate(affected_gene = if_else(is.na(gene1), "gene2", "gene1"))
  
  out <- bind_rows(left, right, intergenic)
  
  # Build DNA_change once, based on exonJoined* states
  out %>%
    mutate(
      frame = NA,
      DNA_change = case_when(
        exonJoined1 == "intergenic" ~ paste0(gene2, ":", transcriptId2, ":e.", exonJoined2, "::", "intergenic"),
        exonJoined2 == "intergenic" ~ paste0(gene1, ":", transcriptId1, ":e.", exonJoined1, "::", "intergenic"),
        
        exonJoined1 == "antisense" ~ paste0(gene2, ":", transcriptId2, ":e.", exonJoined2, "::", gene1, ":", transcriptId1, ":antisense"),
        exonJoined2 == "antisense" ~ paste0(gene1, ":", transcriptId1, ":e.", exonJoined1, "::", gene2, ":", transcriptId2, ":antisense")
      )
    ) %>%
    relocate(affected_gene, frame, DNA_change, mutation, .after = last_col()) %>%
    arrange(sample, svId)
}

#' Generate annotations for loss of transcription cases
#' Filters annotated structural variant data for loss of transcription
#' and  create gene annotations.
#'
#' @param df A data frame or tibble containing annotated SV calls.
#'   Must include columns such as `effect`, `Type_simple`, `strand1`, `orientation1`,
#'   and the paired columns `<base>1` / `<base>2` for each entry in `bases`.
#'
#' @return A tibble filtered to loss of transcription events.
#'  Includes derived columns `DNA_change` and `mutation`.
get_loss_of_transcription <- function(df) {

  base <- df %>%
    filter(effect == "loss_transcription") %>%
    mutate(
      # for loss of transcription use exonDown*
      exonJoined1 = if_else(!is.na(gene1), as.character(exonDown1), "intergenic"),
      exonJoined2 = if_else(!is.na(gene2), as.character(exonDown2), "intergenic"),
      chr1r = chr_rank(chromosome1),
      chr2r = chr_rank(chromosome2),
      swap = chr2r < chr1r | (chr2r == chr1r & pos2 < pos1),
      mutation = if_else(
        swap,
        paste0(chromosome2, ":", pos2, "::", chromosome1, ":", pos1),
        paste0(chromosome1, ":", pos1, "::", chromosome2, ":", pos2)
      )
    ) %>%
    select(-chr1r, -chr2r, -swap)
  
  # rows with both genes → create left/right antisense rows
  both <- base %>% filter(!is.na(gene1) & !is.na(gene2))
  
  left <- both %>%
    mutate(
      exonJoined2   = "antisense",
      affected_gene = "gene1"
    )
  
  right <- both %>%
    mutate(
      exonJoined1   = "antisense",
      affected_gene = "gene2"
    )
  
  # rows where one side is intergenic → keep single row, mark affected_gene
  intergenic <- base %>%
    filter(is.na(gene1) | is.na(gene2)) %>%
    mutate(affected_gene = if_else(is.na(gene1), "gene2", "gene1"))
  
  out <- bind_rows(left, right, intergenic)
  
  # Build DNA_change once, using current exonJoined* (which may be "antisense"/"intergenic")
  out %>%
    mutate(
      frame = NA,
      DNA_change = case_when(
        exonJoined1 == "intergenic" ~ paste0("intergenic", "::", gene2, ":", transcriptId2, ":e.", exonJoined2),
        exonJoined2 == "intergenic" ~ paste0("intergenic", "::", gene1, ":", transcriptId1, ":e.", exonJoined1),
        exonJoined1 == "antisense" ~ paste0(gene1, ":", transcriptId1, ":antisense", "::", gene2, ":", transcriptId2, ":e.", exonJoined2),
        exonJoined2 == "antisense" ~ paste0(gene2, ":", transcriptId2, ":antisense", "::", gene1, ":", transcriptId1, ":e.", exonJoined1)
      )
    ) %>%
    relocate(affected_gene, frame, DNA_change, mutation, .after = last_col()) %>%
    arrange(sample, svId)
}

#' Generate annotations for no-impact events
#' Filters annotated structural variant data for no-impact events
#' and  create gene annotations.
#'
#' @param df A data frame or tibble containing annotated SV calls.
#'   Must include columns such as `effect`, `Type_simple`, `strand1`, `orientation1`,
#'   and the paired columns `<base>1` / `<base>2` for each entry in `bases`.
#'
#' @return A tibble filtered to no-impact events.
#'  Includes derived columns `DNA_change` and `mutation`.
get_no_impact <- function(df) {

  base <- df |>
    filter(effect == "no_impact") |>
    mutate(chr1r = chr_rank(chromosome1),
           chr2r = chr_rank(chromosome2),
           swap = chr2r < chr1r | (chr2r == chr1r & pos2 < pos1),
           mutation = if_else(
             swap,
             paste0(chromosome2, ":", pos2, "::", chromosome1, ":", pos1),
             paste0(chromosome1, ":", pos1, "::", chromosome2, ":", pos2)
           )) |>
    select(-chr1r, -chr2r, -swap)

  
  # rows with both genes → create left/right antisense rows
  both <- base %>% filter(!is.na(gene1) & !is.na(gene2))
  
  left <- both %>%
    mutate(
      exonJoined1 = as.character(exonDown1),
      exonJoined2   = "antisense",
      affected_gene = "gene1",
      DNA_change = paste0(gene2, ":", transcriptId2, ":antisense", "::", gene1, ":", transcriptId1, ":e.", exonJoined1)
    )
  
  right <- both %>%
    mutate(
      exonJoined1   = "antisense",
      exonJoined2 = as.character(exonDown2),
      affected_gene = "gene2",
      DNA_change = paste0(gene1, ":", transcriptId1, ":antisense", "::", gene2, ":", transcriptId2, ":e.", exonJoined2)
    )
  
  # rows where one side is intergenic
  intergenic <- base %>%
    filter(is.na(gene1) | is.na(gene2)) %>%
    mutate(exonJoined1 = case_when(
             is.na(gene1) ~ "intergenic",
             strand1 == 1 ~ as.character(exonDown1),
             strand1 == -1 ~ as.character(exonUp1)
           ),
           exonJoined2 = case_when(
             is.na(gene2) ~ "intergenic",
             strand2 == 1 ~ as.character(exonUp2),
             strand2 == -1 ~ as.character(exonDown2)
           ),
           affected_gene = if_else(is.na(gene1), "gene2", "gene1"),
           DNA_change = case_when(
             is.na(gene1) & strand2 == -1 ~ paste0("intergenic::", gene2, ":", transcriptId2, ":e.", exonJoined2),
             is.na(gene1) & strand2 == 1 ~ paste0(gene2, ":", transcriptId2, ":e.", exonJoined2, "::intergenic"),
             is.na(gene2) & strand1 == -1 ~ paste0(gene1, ":", transcriptId1, ":e.", exonJoined1, "::intergenic"),
             is.na(gene2) & strand1 == 1 ~ paste0("intergenic::", gene1, ":", transcriptId1, ":e.", exonJoined1)
           ))
  
  out <- bind_rows(left, right, intergenic)
  
  # Build DNA_change once, using current exonJoined* (which may be "antisense"/"intergenic")
  out %>%
    mutate(frame = NA) %>%
    relocate(affected_gene, frame, DNA_change, mutation, .after = last_col()) %>%
    arrange(sample, svId)
}

#' Refine promoters annotation
#' Cleans and refines the final annotated SV table by:
#' \itemize{
#'   \item Removing \code{truncated_transcripts} events where one exon side
#'         is labelled \code{"antisense"} and the other \code{"0"} (or vice versa).
#'   \item Reclassifying events with exon \code{"0"} into specific promoter-related
#'         categories based on the original \code{effect}:
#'         \itemize{
#'           \item \code{fusion} → \code{promoter_hijacking}
#'           \item \code{simple_deletion} → \code{promoter_deletion}
#'           \item \code{simple_duplication} → \code{truncated_duplication}
#'         }
#' }
#'
#' @param df A data frame or tibble containing the final annotated SVs.
#'   Must include the columns \code{effect}, \code{exonJoined1}, and \code{exonJoined2}.
#'
#' @return a tibble identical to the input but with filtered rows and possibly updated
#' \code{effect} values according to the described rules.
refine_promoters <- function(df) {
  df %>%
    mutate(
      .e1 = as.character(exonJoined1),
      .e2 = as.character(exonJoined2)
    ) %>%
    # drop truncated_transcripts with one "antisense" and the other "0"
    filter(!(effect == "truncated_transcripts" &
               ((.e1 == "antisense" & .e2 == "0") | 
                (.e1 == "0"         & .e2 == "antisense")))) %>%
    # recode effect when either exonJoined == "0"
    mutate(effect = case_when(
        (.e1 == "0" | .e2 == "0") & effect == "fusion"             ~ "promoter_hijacking",
        (regionType1 == "UPSTREAM" | regionType2 == "UPSTREAM") & effect == "simple_deletion"    ~ "promoter_deletion",
        (.e1 == "0" | .e2 == "0") & effect == "simple_duplication" ~ "truncated_duplication",
        (.e1 != "0" & .e2 != "0" & .e1 != .e2) & (regionType1 == "UPSTREAM" | regionType2 == "UPSTREAM") & effect == "simple_duplication" ~ "truncated_duplication",
        (.e1 != "0" & .e2 != "0" & .e1 == .e2) & (regionType1 == "UPSTREAM" | regionType2 == "UPSTREAM") & effect == "simple_duplication" ~ "noncoding_duplication",
        TRUE ~ effect)) %>%
    select(-.e1, -.e2)
}

#' Collect multiple variant event types in one call
#' Convenience wrapper around your per-event functions (`get_fusions()`,
#' `get_simple_deldup()`, `get_truncations()`, `get_loss_of_transcription()`,
#' `get_no_impact()`). It runs the selected generators on the same input
#' data frame and returns a single, schema-aligned tibble ready for downstream use.
#'
#' @param df A data frame/tibble of annotated SV calls.
#' @param include Character vector of event groups to include. Any of:
#'   `"fusion"`, `"simple_deldup"`, `"truncated_transcripts"`,
#'   `"loss_transcription"`, `"no_impact"`. Defaults to all.
#' @param fusion_bases Character vector of base column names (e.g.
#'   `"gene"`, `"transcriptId"`, `"chromosome"`, `"pos"`, `"exonJoined"`)
#'   to be swapped inside `get_fusions()`.  
#'   **Required** if `"fusion"` is in `include`.
#' @param as_list Logical; if `TRUE`, return a named list with each event table
#'   separately. If `FALSE` (default), return a single tibble via `bind_rows()`.
#'
#' @return A tibble (default) with all selected events bound together, or a list
#' of tibbles if `as_list = TRUE`. Columns are standardized to a common schema.
get_variants <- function(
    df,
    include = c("fusion", "simple_deldup", "truncated_transcripts",
                "loss_transcription", "no_impact"),
    fusion_bases = NULL,
    as_list = FALSE
) {
  include <- intersect(
    include,
    c("fusion","simple_deldup","truncated_transcripts","loss_transcription","no_impact")
  )
  
  # Check that if fusion is extracted, then the column names to swap are provided
  if ("fusion" %in% include && is.null(fusion_bases)) {
    stop("Argument 'fusion_bases' must be provided when 'fusion' is in include.")
  }
  
  pieces <- list()
  
  if ("fusion" %in% include) {
    pieces$fusion <- get_fusions(df, fusion_bases)
  }
  if ("simple_deldup" %in% include) {
    pieces$simple_deldup <- get_simple_deldup(df, fusion_bases)
  }
  if ("truncated_transcripts" %in% include) {
    pieces$truncated_transcripts <- get_truncations(df)
  }
  if ("loss_transcription" %in% include) {
    pieces$loss_transcription <- get_loss_of_transcription(df)
  }
  if ("no_impact" %in% include) {
    pieces$no_impact <- get_no_impact(df)
  }
  
  if (as_list) return(pieces)
  
  # drop NULLs (if any were skipped) and bind
  pieces <- Filter(Negate(is.null), pieces)
  if (!length(pieces)) {
    # return an empty, schema-aligned tibble if nothing selected
    return(finalize_schema(df[0, , drop = FALSE]))
  }

  combined <- dplyr::bind_rows(pieces) %>%
    dplyr::select(
      1:5, Type_simple, 7:14,
      exonJoined1, 15:24, exonJoined2, 
      25:26, 30, 34,
      affected_gene, DNA_change, mutation
    ) %>%
    dplyr::rename(Type = Type_simple)

  refined <- refine_promoters(combined)
  return(refined)
}
