#!/usr/bin/env Rscript
# ============================================================
# compare_cluster.R
#
# Compares subregion-based transmission clusters (ClusterPicker and
# HIV-TRACE) against full-length-genome clusters, within each method:
#   ClusterPicker subregion  vs  ClusterPicker full-length reference
#   HIV-TRACE subregion      vs  HIV-TRACE full-length reference
#
# Metrics reported per region x method x GD threshold x bootstrap:
#   - V_measure_refclustered (HEADLINE): V-measure scored only over
#     sequences whose full-length reference cluster has >=2 members,
#     so the large mutually-singleton background typical of HIV
#     clustering doesn't saturate the score near 1.
#   - ARI: Adjusted Rand Index (chance-corrected, not singleton-inflated).
#   - pw_precision / pw_recall / pw_f1: pairwise precision/recall/F1 on
#     co-clustered ("linked") pairs, i.e. how much of the full-length
#     linkage structure a subregion recovers (recall) and how much of
#     what it links is also linked in full-length (precision).
#
# Usage:
#   Rscript compare_cluster.R regions.txt
#   Rscript compare_cluster.R regions.txt bs90            # one bootstrap, all GDs
#   Rscript compare_cluster.R regions.txt bs90 0.03        # one bootstrap, one GD
#   Rscript compare_cluster.R regions.txt "" 0.03          # all bootstraps, one GD
#
# `regions.txt` is a plain-text file with one region name per line
# (blank lines and lines starting with # are ignored). Region names
# must match the FASTA/ClusterPicker/HIV-TRACE file naming convention
# described in README.md.
# ============================================================

suppressPackageStartupMessages({
  library(jsonlite)
  library(dplyr)
  library(stringr)
  library(tibble)
  library(clevr)
  library(readr)
})

# ─── CONFIGURATION — edit these paths before running ──────────────────

PROJECT_ROOT <- "path/to/HIVphyloSeq_V2"     # repo/analysis root
COHORT       <- "2005-2024_US_CA"            # cohort subfolder name

CP_BASE    <- file.path(PROJECT_ROOT, "Cluster", "ClusterPicker", COHORT)
TRACE_BASE <- file.path(PROJECT_ROOT, "Cluster", "HIV_trace", COHORT)
ALN_BASE   <- file.path(PROJECT_ROOT, "Alignment", COHORT)  # sequence universe (FASTA headers)
OUTPUT_DIR <- file.path(PROJECT_ROOT, "Cluster", "V_measure_Comparisons", COHORT)

GD_THRESHOLDS    <- seq(from = 0.005, to = 0.045, by = 0.005)  # genetic distance thresholds swept
BS_TAGS          <- c("bs70", "bs90", "bs99")                  # ClusterPicker bootstrap thresholds
BS_NUMBER_FOLDER <- "bs_1000"                                  # bootstrap-replicate folder name
CP_LOG_FILE      <- "input_clusterPicks_log.txt"

# Full-length reference clustering (fixed comparison point for every region)
REF_REGION <- "fulllength"
REF_GD     <- 0.030
REF_BS_TAG <- "bs90"

dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ─── HELPERS ───────────────────────────────────────────────────────────

gd_tag <- function(gd) sprintf("gd%03d", as.integer(round(gd * 1000)))

# Sequence-ID universe for a region, read from its alignment FASTA headers.
# IDs are normalised to match parse_cp()/parse_trace() (drop "Seq_" prefix).
read_fasta_ids <- function(region) {
  path <- file.path(ALN_BASE, paste0("DECIPHER_aligned_", region, "_sequences.fasta"))
  if (!file.exists(path)) {
    message("  FASTA missing (cannot build sequence universe): ", path)
    return(character(0))
  }
  ids <- sub("^>", "", grep("^>", readLines(path, warn = FALSE), value = TRUE))
  ids <- sub("\\s.*$", "", ids)
  sub("^Seq_", "", ids)
}

# Parse a ClusterPicker log into a SequenceID -> cluster table.
parse_cp <- function(log_path) {
  if (!file.exists(log_path)) {
    message("ClusterPicker log missing: ", log_path)
    return(NULL)
  }
  rows <- grep("^\\d+\\t", readLines(log_path, warn = FALSE), value = TRUE)
  if (!length(rows)) {
    message("No data rows in ClusterPicker log: ", log_path)
    return(NULL)
  }
  bind_rows(lapply(rows, function(r) {
    fields <- str_split(r, "\\t", n = 6)[[1]]
    tibble(
      SequenceID = str_remove_all(fields[4], "\\[|\\]|\\s") |>
        str_split(",") |> unlist() |> str_replace("^Seq_", ""),
      CP_cluster = as.integer(fields[1])
    )
  }))
}

# Parse an HIV-TRACE JSON into a SequenceID -> cluster table.
# Handles both the flat (compact_json) and nested ($values) node layouts.
parse_trace <- function(json_path) {
  if (!file.exists(json_path)) {
    message("HIV-TRACE JSON missing: ", json_path)
    return(NULL)
  }
  j <- tryCatch(fromJSON(json_path, simplifyVector = FALSE), error = function(e) NULL)
  nodes <- j$trace_results$Nodes
  if (is.null(nodes)) {
    message("Invalid HIV-TRACE JSON: ", json_path)
    return(NULL)
  }

  cluster_values <- nodes$cluster
  if (is.list(cluster_values) && !is.null(cluster_values$values)) cluster_values <- cluster_values$values
  cluster_ids <- nodes$id
  if (is.list(cluster_ids) && !is.null(cluster_ids$values)) cluster_ids <- cluster_ids$values

  if (is.null(cluster_values) || is.null(cluster_ids) || length(cluster_values) != length(cluster_ids)) {
    message("Could not extract matching Nodes$cluster / Nodes$id from: ", json_path)
    return(NULL)
  }
  tibble(
    SequenceID    = unlist(cluster_ids) |> as.character() |> str_replace("^Seq_", ""),
    TRACE_cluster = unlist(cluster_values) |> as.integer()
  )
}

# Expand a cluster table to the full sequence universe: unclustered IDs
# each get their own unique singleton label.
build_labels <- function(df, universe_ids, cluster_col) {
  lab <- setNames(rep(NA_integer_, length(universe_ids)), universe_ids)
  if (!is.null(df) && nrow(df) > 0) {
    keep <- df$SequenceID %in% universe_ids
    lab[df$SequenceID[keep]] <- as.integer(df[[cluster_col]][keep])
  }
  max_label <- suppressWarnings(max(lab, na.rm = TRUE))
  if (!is.finite(max_label)) max_label <- 0L
  singleton_idx <- which(is.na(lab))
  lab[singleton_idx] <- max_label + seq_along(singleton_idx)
  lab
}

# Pairwise precision/recall/F1 on co-clustered pairs, via contingency
# counts (avoids enumerating all N-choose-2 pairs).
pairwise_prf <- function(pred, ref) {
  n_pairs <- function(x) x * (x - 1) / 2
  joint <- dplyr::count(tibble(pred = pred, ref = ref), pred, ref, name = "n")
  shared_pairs <- sum(n_pairs(joint$n))
  pred_pairs   <- sum(n_pairs(as.integer(table(pred))))
  ref_pairs    <- sum(n_pairs(as.integer(table(ref))))
  list(
    pw_precision = if (pred_pairs > 0) shared_pairs / pred_pairs else NA_real_,
    pw_recall    = if (ref_pairs  > 0) shared_pairs / ref_pairs  else NA_real_,
    pw_f1        = if (pred_pairs + ref_pairs > 0) 2 * shared_pairs / (pred_pairs + ref_pairs) else NA_real_
  )
}

# Sequence IDs whose full-length reference cluster has >= 2 members
# present in this universe (used for the headline V-measure).
ref_clustered_ids <- function(ref_df, universe_ids, ref_col) {
  in_universe <- ref_df[ref_df$SequenceID %in% universe_ids, , drop = FALSE]
  if (nrow(in_universe) == 0) return(character(0))
  sizes <- table(in_universe[[ref_col]])
  multi_member <- names(sizes)[sizes >= 2]
  unique(in_universe$SequenceID[as.character(in_universe[[ref_col]]) %in% multi_member])
}

calc_metrics <- function(pred_df, ref_df, universe_ids, pred_col, ref_col) {
  if (length(universe_ids) == 0) {
    return(list(V_measure_refclustered = NA_real_, ARI = NA_real_,
                pw_precision = NA_real_, pw_recall = NA_real_, pw_f1 = NA_real_, n_ref_clustered = 0L))
  }

  pred_labels <- build_labels(pred_df, universe_ids, pred_col)
  ref_labels  <- build_labels(ref_df,  universe_ids, ref_col)

  eval_ids <- ref_clustered_ids(ref_df, universe_ids, ref_col)
  vm_headline <- NA_real_
  if (length(eval_ids) >= 2) {
    vm_headline <- clevr::v_measure(
      build_labels(pred_df, eval_ids, pred_col),
      build_labels(ref_df,  eval_ids, ref_col)
    )
  }

  c(
    list(V_measure_refclustered = vm_headline,
         ARI                    = clevr::adj_rand_index(pred_labels, ref_labels)),
    pairwise_prf(pred_labels, ref_labels),
    list(n_ref_clustered = length(eval_ids))
  )
}

# One result row for a subregion vs. its full-length reference, for a
# given clustering method (ClusterPicker or HIV-TRACE).
compare_region <- function(region, method, pred_df, ref_df, universe_ids,
                            pred_col, ref_col, gd, bs_tag = NA_character_) {
  m <- calc_metrics(pred_df, ref_df, universe_ids, pred_col, ref_col)
  tibble(
    Region = region, Method = method, GD_percent = gd * 100, Bootstrap = bs_tag,
    Reference_Method = method, Reference_GD = REF_GD * 100, Reference_Bootstrap = REF_BS_TAG,
    V_measure_refclustered = round(m$V_measure_refclustered, 6),
    ARI          = round(m$ARI, 6),
    pw_recall    = round(m$pw_recall, 6),
    pw_precision = round(m$pw_precision, 6),
    pw_f1        = round(m$pw_f1, 6),
    n_ref_clustered = m$n_ref_clustered,
    n_universe      = length(universe_ids)
  )
}

# ─── MAIN ────────────────────────────────────────────────────────────

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript compare_cluster.R regions.txt [bs_tag] [gd_fraction]")

regions <- readLines(args[1], warn = FALSE)
regions <- trimws(regions)
regions <- regions[nzchar(regions) & !startsWith(regions, "#")]

if (length(args) >= 2 && nzchar(args[2])) BS_TAGS <- args[2]
if (length(args) >= 3 && nzchar(args[3])) {
  GD_THRESHOLDS <- as.numeric(args[3])
  if (is.na(GD_THRESHOLDS)) stop("Could not parse GD '", args[3], "' as a number.")
}

message("Loading full-length reference clusters...")
ref_cp    <- parse_cp(file.path(CP_BASE, paste0(REF_REGION, "_CP"), BS_NUMBER_FOLDER,
                                 paste0(REF_BS_TAG, "_", gd_tag(REF_GD)), CP_LOG_FILE))
ref_trace <- parse_trace(file.path(TRACE_BASE, REF_REGION,
                                    paste0(REF_REGION, "_hivtrace_GD", sprintf("%03d", as.integer(round(REF_GD * 1000))), ".json")))
if (is.null(ref_cp))    stop("Full-length ClusterPicker reference not found.")
if (is.null(ref_trace)) stop("Full-length HIV-TRACE reference not found.")

ref_universe <- read_fasta_ids(REF_REGION)
results <- list()

for (region in regions) {
  message(sprintf("=== %s ===", region))
  universe_ids <- intersect(read_fasta_ids(region), ref_universe)
  if (length(universe_ids) == 0) {
    message("  No shared sequence universe; skipping region.")
    next
  }

  for (gd in GD_THRESHOLDS) {

    # ClusterPicker
    for (bs_tag in BS_TAGS) {
      cp_path <- file.path(CP_BASE, paste0(region, "_CP"), BS_NUMBER_FOLDER,
                           paste0(bs_tag, "_", gd_tag(gd)), CP_LOG_FILE)
      cp_df <- parse_cp(cp_path)
      if (!is.null(cp_df)) {
        results[[length(results) + 1]] <- compare_region(
          region, "ClusterPicker", cp_df, ref_cp, universe_ids, "CP_cluster", "CP_cluster", gd, bs_tag
        )
      }
    }

    # HIV-TRACE
    trace_path <- file.path(TRACE_BASE, region,
                            paste0(region, "_hivtrace_GD", sprintf("%03d", as.integer(round(gd * 1000))), ".json"))
    trace_df <- parse_trace(trace_path)
    if (!is.null(trace_df)) {
      results[[length(results) + 1]] <- compare_region(
        region, "HIV-TRACE", trace_df, ref_trace, universe_ids, "TRACE_cluster", "TRACE_cluster", gd
      )
    }
  }
}

if (length(results) > 0) {
  out_path <- file.path(OUTPUT_DIR, "vmeasure_results.tsv")
  bind_rows(results) |>
    arrange(Region, GD_percent, Method, Bootstrap) |>
    write_tsv(out_path)
  message("\nWrote results to ", out_path)
} else {
  message("No results generated.")
}
