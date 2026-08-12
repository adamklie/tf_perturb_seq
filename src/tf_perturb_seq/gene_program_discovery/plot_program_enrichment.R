#!/usr/bin/env Rscript
# =============================================================================
# plot_program_enrichment.R
#
# Appended pipeline - cNMF program analysis using iPSC vs hepatocyte DEGs:
#
#   1. DESeq2: iPSC vs mature hepatocyte DEGs (pseudobulk)
#   2. KS test enrichment: top 5 programs enriched for DEGs
#   3. Fisher's exact test: top TF regulators per program
#   4. UMAP (perturbation dataset): top 5 programs, annotated with top TF
#   5. UMAP (mature hepatocyte reference): inferred NNLS usages, top 5 programs
#
# Dependencies:
#   DESeq2, Seurat, MuData (via MuDataSeurat or anndata), ggplot2, dplyr,
#   tidyr, readr, forcats, patchwork, ggsci, ggrepel
# =============================================================================

.libPaths("R-bak/x86_64-pc-linux-gnu-library/4.4/")

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(forcats)
  library(patchwork)
  library(ggsci)
  library(ggrepel)
  library(Seurat)
  library(DESeq2)
  library(Matrix)
  library(zellkonverter)
})

# ---- Paths ------------------------------------------------------------------
IPSC_RDS       <- "/hpc/group/gersbachlab/seg95/scSALSA/ipscs.rds"
HEP_RDS        <- "/hpc/group/gersbachlab/seg95/helen_data/HumanLiverSeurat_QC_donorBatchC.rds"
SPECTRA_FILE   <- "/hpc/group/gersbachlab/seg95/tf_perturb_seq/datasets/Gersbach_WTC11-hepatocyte-differentiation_TF-Perturb-seq/PerturbNMF/Result/051026_gersbach_iHep_torch_minibatch/Inference/Inference.spectra.k_80.dt_2_0.consensus.txt"
#TRANS_FILE     <- "/hpc/group/gersbachlab/seg95/tf_perturb_seq/datasets/Gersbach_WTC11-hepatocyte-differentiation_TF-Perturb-seq/calibrated_outs/Gersbach_WTC11-hepatocyte-differentiation_TF-Perturb-seq_sara_synapse_syn74842722_calibrated_trans_results.tsv"
TRANS_FILE     <- "/hpc/group/gersbachlab/seg95/CRISPR_Pipeline/helen_output_poolabcd_prod_v10/pipeline_outputs/gersbach_ihep_calibrated_trans_results.tsv"
MUDATA_FILE    <- "/hpc/group/gersbachlab/seg95/tf_perturb_seq/datasets/Gersbach_WTC11-hepatocyte-differentiation_TF-Perturb-seq/post_pipeline_processing/inference_mudata_outlier_guide_filt.h5mu"
PERTURB_USAGE  <- "/hpc/group/gersbachlab/seg95/tf_perturb_seq/datasets/Gersbach_WTC11-hepatocyte-differentiation_TF-Perturb-seq/PerturbNMF/Result/051026_gersbach_iHep_torch_minibatch/Inference/Inference.usages.k_80.dt_2_0.consensus.txt"
HEP_DIR_UMAP   <- "/hpc/group/gersbachlab/seg95/tf_perturb_seq/datasets/Gersbach_WTC11-hepatocyte-differentiation_TF-Perturb-seq/calibrated_outs_manual/hepatocyte_marker_figures"

HEP_DIR        <- "/hpc/group/gersbachlab/seg95/tf_perturb_seq/datasets/Gersbach_WTC11-hepatocyte-differentiation_TF-Perturb-seq/NNLS"  
NNLS_USAGE     <- file.path(HEP_DIR, "usage_matrix_donorBatchC.csv")   # output from run_nnls.py
OUT_DIR        <- file.path(HEP_DIR, "perturbo_manual") 


# Thresholds
FDR_THRESH     <- 0.05
LFC_THRESH     <- 1.0     # |log2FC| for DESeq2
TOP_N_GENES    <- 200     # genes per program for Fisher's test
TOP_N_PROGRAMS <- 5
TOP_N_LABEL    <- 1       # TF labels per UMAP panel

# ---- Shared helpers ---------------------------------------------------------
iterm <- pal_d3("category10")(10)

# Named palette for direction categories
DIR_COLORS <- c(
  "up"   = iterm[1],
  "down" = iterm[2],
  "ns"   = "grey80"
)

theme_hep <- function(base_size = 12) {
  theme_classic(base_size = base_size) %+replace%
    theme(
      strip.background  = element_blank(),
      strip.text        = element_text(face = "bold", size = base_size),
      axis.line         = element_line(colour = "grey30"),
      axis.ticks        = element_line(colour = "grey30", linewidth = 0.5),
      axis.ticks.length = unit(0.2, "cm"),
      axis.text         = element_text(size = base_size + 1),
      axis.title        = element_text(size = base_size + 2),
      plot.title        = element_text(face = "bold", size = base_size + 1,
                                       hjust = 0),
      plot.subtitle     = element_text(size = base_size - 1, colour = "grey40",
                                       hjust = 0),
      legend.key.size   = unit(0.4, "cm"),
      panel.spacing     = unit(0.8, "lines")
    )
}

save_plot <- function(p, path, width = 8, height = 6, dpi = 200) {
  ggsave(path, p, width = width, height = height, dpi = dpi, bg = "white")
  message("  Saved: ", path)
}

read_tsv_safe <- function(path, ...) {
  if (!file.exists(path)) {
    message("  SKIP (not found): ", path)
    return(NULL)
  }
  read_tsv(path, show_col_types = FALSE, ...)
}

# Glob helper that handles spaces and special characters in filenames
list_tsv <- function(dir, prefix, suffix = ".tsv") {
  all_files <- list.files(dir, full.names = TRUE)
  all_files[startsWith(basename(all_files), prefix) &
              endsWith(basename(all_files), suffix)]
}


if (!exists("theme_hep")) {
  iterm <- pal_d3("category10")(10)
  theme_hep <- function(base_size = 12) {
    theme_classic(base_size = base_size) %+replace%
      theme(
        strip.background  = element_blank(),
        strip.text        = element_text(face = "bold", size = base_size),
        axis.line         = element_line(colour = "grey30"),
        axis.ticks        = element_line(colour = "grey30", linewidth = 0.5),
        axis.ticks.length = unit(0.2, "cm"),
        axis.text         = element_text(size = base_size + 1),
        axis.title        = element_text(size = base_size + 2),
        plot.title        = element_text(face = "bold", size = base_size + 1, hjust = 0),
        plot.subtitle     = element_text(size = base_size - 1, colour = "grey40", hjust = 0),
        legend.key.size   = unit(0.4, "cm"),
        panel.spacing     = unit(0.8, "lines")
      )
  }
  save_plot <- function(p, path, width = 8, height = 6, dpi = 200) {
    ggsave(path, p, width = width, height = height, dpi = dpi, bg = "white")
    message("  Saved: ", path)
  }
}


# =============================================================================
# 1. DESeq2: iPSC vs mature hepatocyte
# =============================================================================
message("\n--- Step 1: DESeq2 iPSC vs mature hepatocyte ---")

ipsc <- readRDS(IPSC_RDS)
hep  <- readRDS(HEP_RDS)

# -- Pseudobulk helper ---------------------------------------------------------
# Aggregates raw counts per donor/sample into a single pseudobulk sample.
# Expects a Seurat object with raw counts in the "RNA" assay.
make_pseudobulk <- function(seu, group_col, label) {
  # Use RNA counts layer; fall back to @counts slot for older Seurat versions
  counts_mat <- tryCatch(
    LayerData(seu, assay = "RNA", layer = "counts"),
    error = function(e) GetAssayData(seu, assay = "RNA", slot = "counts")
  )

  meta <- seu@meta.data

  if (!group_col %in% colnames(meta)) {
    # No donor column ? treat all cells as a single pseudobulk sample
    message("  '", group_col, "' not found in ", label,
            " metadata; treating all cells as one pseudobulk sample.")
    pb <- Matrix::rowSums(counts_mat)
    pb_mat <- matrix(pb, ncol = 1,
                     dimnames = list(names(pb), paste0(label, "_all")))
    coldata <- data.frame(
      sample    = paste0(label, "_all"),
      condition = label,
      row.names = paste0(label, "_all")
    )
  } else {
    donors  <- unique(meta[[group_col]])
    pb_list <- lapply(donors, function(d) {
      cells <- rownames(meta)[meta[[group_col]] == d]
      Matrix::rowSums(counts_mat[, cells, drop = FALSE])
    })
    pb_mat  <- do.call(cbind, pb_list)
    colnames(pb_mat) <- paste0(label, "_", donors)
    coldata <- data.frame(
      sample    = colnames(pb_mat),
      condition = label,
      row.names = colnames(pb_mat)
    )
  }
  list(counts = pb_mat, coldata = coldata)
}

# Detect donor column (common naming conventions)
detect_donor_col <- function(seu) {
  candidates <- c("donor", "donor_id", "orig.ident", "sample", "patient",
                  "individual", "batch")
  found <- intersect(candidates, colnames(seu@meta.data))
  if (length(found) > 0) found[1] else NULL
}

hep_donor  <- detect_donor_col(hep)

message("  Hepatocyte donor column: ", ifelse(is.null(hep_donor),  "none (single pseudobulk)", hep_donor))

pb_ipsc <- make_pseudobulk(ipsc, "__none__", "iPSC")
pb_hep  <- make_pseudobulk(hep,  ifelse(is.null(hep_donor),  "__none__", hep_donor),  "hepatocyte")

# Intersect genes and merge
shared_genes <- intersect(rownames(pb_ipsc$counts), rownames(pb_hep$counts))
message("  Shared genes for DESeq2: ", length(shared_genes))

counts_combined <- cbind(
  pb_ipsc$counts[shared_genes, , drop = FALSE],
  pb_hep$counts[shared_genes,  , drop = FALSE]
)
coldata_combined <- rbind(pb_ipsc$coldata, pb_hep$coldata)
coldata_combined$condition <- factor(coldata_combined$condition,
                                     levels = c("iPSC", "hepatocyte"))

# Run DESeq2
dds <- DESeqDataSetFromMatrix(
  countData = round(as.matrix(counts_combined)),
  colData   = coldata_combined,
  design    = ~ condition
)
# dplyr::filter low-count genes
keep <- rowSums(counts(dds) >= 10) >= max(1, floor(ncol(dds) * 0.25))
dds  <- dds[keep, ]
message("  Genes after low-count dplyr::filter: ", sum(keep))

dds     <- DESeq(dds, quiet = TRUE)
res_raw <- results(dds,
                   contrast        = c("condition", "hepatocyte", "iPSC"),
                   alpha           = FDR_THRESH,
                   independentdplyr::filtering = TRUE)
res_df  <- as.data.frame(res_raw) %>%
  tibble::rownames_to_column("gene") %>%
  dplyr::filter(!is.na(padj)) %>%
  arrange(padj)

# Up in hepatocyte = positive log2FC
degs_up   <- res_df %>% dplyr::filter(padj < FDR_THRESH, log2FoldChange >  LFC_THRESH)
degs_down <- res_df %>% dplyr::filter(padj < FDR_THRESH, log2FoldChange < -LFC_THRESH)

message("  DEGs up in hepatocyte: ",   nrow(degs_up))
message("  DEGs up in iPSC (down): ", nrow(degs_down))

write_tsv(res_df,    file.path(OUT_DIR, "deseq2_iPSC_vs_hepatocyte_all.tsv"))
write_tsv(degs_up,   file.path(OUT_DIR, "deseq2_iPSC_vs_hepatocyte_up_in_hep.tsv"))
write_tsv(degs_down, file.path(OUT_DIR, "deseq2_iPSC_vs_hepatocyte_up_in_iPSC.tsv"))
message("  DESeq2 results saved.")


# =============================================================================
# 2. KS enrichment: top programs for hepatocyte-upregulated DEGs
# =============================================================================
message("\n--- Step 2: KS enrichment of programs for hepatocyte DEGs ---")

spectra <- read_tsv(SPECTRA_FILE, show_col_types = FALSE)

# cNMF spectra: rows = programs, first column = program name
program_col <- colnames(spectra)[1]
program_ids <- spectra[[program_col]]
gene_cols   <- colnames(spectra)[-1]

# Tidy: long format (program x gene x weight)
spectra_long <- spectra %>%
  pivot_longer(cols = all_of(gene_cols),
               names_to  = "gene",
               values_to = "weight")

deg_set <- degs_up$gene   # hepatocyte-upregulated DEGs

# KS test per program: are DEGs ranked higher than background?
ks_results <- spectra_long %>%
  group_by(.data[[program_col]]) %>%
  group_modify(function(df, key) {
    # Rank genes by program weight (descending: rank 1 = highest weight)
    df_ranked <- df %>% arrange(desc(weight))
    ranked_genes <- df_ranked$gene
    
    # Positions (ranks) of DEGs within this program's gene ranking
    deg_positions <- which(ranked_genes %in% deg_set)
    
    if (length(deg_positions) < 5) {
      return(data.frame(ks_stat = NA_real_, p_value = NA_real_,
                        n_deg_in_program = 0L, n_program_genes = nrow(df)))
    }
    
    # One-sided KS test: are DEGs at lower rank positions (= higher weight)?
    ks <- ks.test(deg_positions,
                  seq_along(ranked_genes),
                  alternative = "less")   # "less" = DEG positions skew toward rank 1
    
    data.frame(
      ks_stat          = ks$statistic,
      p_value          = ks$p.value,
      n_deg_in_program = length(deg_positions),
      n_program_genes  = nrow(df)
    )
  }) %>%
  ungroup() %>%
  rename(program = .data[[program_col]]) %>%
  dplyr::filter(!is.na(p_value)) %>%
  mutate(padj = p.adjust(p_value, method = "BH")) %>%
  arrange(padj, desc(ks_stat))

top_programs <- ks_results %>%
  slice_head(n = TOP_N_PROGRAMS) %>%
  pull(program)

message("  Top ", TOP_N_PROGRAMS, " programs enriched for hepatocyte DEGs:")
message(paste0("    ", top_programs, collapse = "\n"))

write_tsv(ks_results, file.path(OUT_DIR, "ks_program_deg_enrichment.tsv"))


# =============================================================================
# 3. Fisher's exact test: top TF regulators per top program
# =============================================================================
message("\n--- Step 3: Fisher's exact test - TF regulators per program ---")

trans <- read_tsv(TRANS_FILE, show_col_types = FALSE)

# Universe = all genes tested in trans results
universe_genes <- unique(trans$tested_gene_symbol)
n_universe     <- length(universe_genes)
message("  Trans result universe: ", n_universe, " genes")

# Top 200 genes per program (by spectra weight)
program_top_genes <- spectra_long %>%
  group_by(.data[[program_col]]) %>%
  slice_max(order_by = weight, n = TOP_N_GENES) %>%
  ungroup() %>%
  rename(program = .data[[program_col]])

# Significant trans genes (FDR < 0.05 for that TF)
sig_trans <- trans %>%
  dplyr::filter(empirical_pval_adj < FDR_THRESH) %>%
  dplyr::select(tf = element_symbol, gene = tested_gene_symbol, log2fc)

fisher_results <- lapply(top_programs, function(prog) {
  prog_genes <- program_top_genes %>%
    dplyr::filter(program == prog) %>%
    pull(gene)
  prog_genes <- intersect(prog_genes, universe_genes)

  if (length(prog_genes) < 10) return(NULL)

  tfs <- unique(sig_trans$tf)

  tf_res <- lapply(tfs, function(tf) {
    tf_genes <- sig_trans %>%
      dplyr::filter(tf == !!tf) %>%
      pull(gene) %>%
      intersect(universe_genes)

    if (length(tf_genes) < 3) return(NULL)

    # 2x2 contingency table
    a <- length(intersect(prog_genes, tf_genes))           # in prog & TF
    b <- length(setdiff(tf_genes, prog_genes))             # in TF, not prog
    c <- length(setdiff(prog_genes, tf_genes))             # in prog, not TF
    d <- n_universe - a - b - c                            # neither

    ft <- fisher.test(matrix(c(a, b, c, d), nrow = 2),
                      alternative = "greater")

    data.frame(
      program   = prog,
      tf        = tf,
      odds_ratio = ft$estimate,
      p_value   = ft$p.value,
      n_overlap = a,
      n_prog    = length(prog_genes),
      n_tf      = length(tf_genes)
    )
  })

  return(bind_rows(tf_res))
})

fisher_df <- bind_rows(fisher_results) %>%
  group_by(program) %>%
  mutate(padj = p.adjust(p_value, method = "BH")) %>%
  ungroup() %>%
  arrange(program, padj)

# Top TF per program (for UMAP annotation)
top_tf_per_program <- fisher_df %>%
  dplyr::filter(padj < FDR_THRESH) %>%
  group_by(program) %>%
  slice_min(order_by = padj, n = TOP_N_LABEL) %>%
  ungroup()

message("  Top TF per program:")
for (prog in top_programs) {
  tf_row <- top_tf_per_program %>% dplyr::filter(program == prog) %>% slice(1)
  if (nrow(tf_row) > 0) {
    message("    ", prog, " -> ", tf_row$tf,
            "  (OR=", round(tf_row$odds_ratio, 2),
            ", FDR=", formatC(tf_row$padj, format = "e", digits = 2), ")")
  } else {
    message("    ", prog, " -> no significant TF")
  }
}

write_tsv(fisher_df,          file.path(OUT_DIR, "fisher_tf_program_regulators.tsv"))
write_tsv(top_tf_per_program, file.path(OUT_DIR, "top_tf_per_program.tsv"))


# =============================================================================
# Shared UMAP plotting helper
# =============================================================================

# gradient palette: white -> orange -> red (program usage)
usage_gradient <- function() {
  scale_colour_gradientn(
    colours = c("grey95", iterm[5], iterm[1], iterm[2]),
    name    = "Usage",
    guide   = guide_colorbar(barheight = 3, barwidth = 0.8)
  )
}

# Build one UMAP panel for a single program
make_usage_umap <- function(df,           # data.frame with UMAP1, UMAP2, usage
                            program_name,
                            top_tf = NULL,
                            subtitle = NULL) {
  # Sort so high-usage cells are on top
  df <- df %>% arrange(usage)

  p <- ggplot(df, aes(x = UMAP1, y = UMAP2, colour = usage)) +
    geom_point(size = 0.25, alpha = 0.6, stroke = 0, shape = 16) +
    usage_gradient() +
    labs(
      title    = program_name,
      subtitle = if (!is.null(top_tf)) paste0("Top regulator: ", top_tf) else subtitle
    ) +
    theme_hep() +
    theme(
      axis.text  = element_blank(),
      axis.ticks = element_blank(),
      axis.line  = element_blank(),
      plot.title    = element_text(face = "bold", size = 10, hjust = 0.5),
      plot.subtitle = element_text(size  = 8,     colour = "grey30", hjust = 0.5)
    )
  p
}

# Combine 5 panels side-by-side with a shared annotation
assemble_usage_panels <- function(panels, title, subtitle = NULL) {
  wrap_plots(panels, nrow = 1) +
    plot_annotation(
      title    = title,
      subtitle = subtitle,
      theme    = theme(
        plot.title    = element_text(face = "bold", size = 13, hjust = 0),
        plot.subtitle = element_text(size = 10, colour = "grey40", hjust = 0)
      )
    )
}


# =============================================================================
# 4. UMAP - perturbation dataset (MuData), top 5 programs
# =============================================================================
message("\n--- Step 4: UMAP - perturbation dataset ---")

tryCatch({
  #perturb_umap_df <- read_csv(file.path(dirname(MUDATA_FILE), "perturb_umap.csv"),
  #  show_col_types = FALSE) %>% rename(cell = 1)
  umap_long <- read_tsv_safe(
    file.path(HEP_DIR_UMAP, "umap_cells_marker_expression_long.tsv")
  )
  
  perturb_umap_df <- unique(umap_long[, c("cell_barcode", "UMAP1", "UMAP2")])
  colnames(perturb_umap_df)[which(colnames(perturb_umap_df) == "cell_barcode")] <- "cell"

  # Load usage matrix (cells x programs)
  perturb_usage <- read_tsv(PERTURB_USAGE, show_col_types = FALSE)
  colnames(perturb_usage)[1] <- "cell"

  # Join UMAP and usage
  perturb_df <- inner_join(perturb_umap_df, perturb_usage, by = "cell")
  message("  Perturbation cells with UMAP + usage: ", nrow(perturb_df))

  # Build panels for top 5 programs
  panels_perturb <- lapply(top_programs, function(prog) {
    top_tf <- top_tf_per_program %>%
      dplyr::filter(program == prog) %>%
      dplyr::slice(1) %>%
      pull(tf)
    top_tf <- if (length(top_tf) == 0) "n.s." else top_tf

    df_prog <- perturb_df %>%
      dplyr::select(UMAP1, UMAP2, usage = all_of(as.character(prog)))

    make_usage_umap(df_prog, paste0("Program ", prog), top_tf = top_tf)
  })

  p_perturb <- assemble_usage_panels(
    panels_perturb,
    title    = "Top 5 hepatocyte-enriched cNMF programs - perturbation dataset",
    subtitle = "Colour: program usage  |  Subtitle: top Fisher's TF regulator"
  )

  save_plot(p_perturb,
            file.path(OUT_DIR, "ggplot_umap_perturb_top5_programs.png"),
            width = 20, height = 5)

}, error = function(e) {
  message("  ERROR in Step 4 (MuData): ", conditionMessage(e))
})

# =============================================================================
# 5. UMAP - mature hepatocyte reference, NNLS inferred usages
# =============================================================================
message("\n--- Step 5: UMAP - mature hepatocyte reference (NNLS usages) ---")

if (!file.exists(NNLS_USAGE)) {
  message("  NNLS usage file not found: ", NNLS_USAGE, " - skipping Step 5.")
} else {
  nnls_usage <- read_csv(NNLS_USAGE, show_col_types = FALSE)

  # First column is cell barcode
  usage_col1 <- colnames(nnls_usage)[1]
  nnls_usage <- nnls_usage %>% rename(cell = all_of(usage_col1))

  # Extract UMAP from the hepatocyte Seurat object
  if (!"umap" %in% names(hep@reductions) &&
      !"UMAP"  %in% names(hep@reductions)) {
    message("  No UMAP found in hepatocyte Seurat object - computing one now...")
    hep <- NormalizeData(hep, verbose = FALSE)
    hep <- FindVariableFeatures(hep, verbose = FALSE)
    hep <- ScaleData(hep, verbose = FALSE)
    hep <- RunPCA(hep, verbose = FALSE)
    hep <- RunUMAP(hep, dims = 1:30, verbose = FALSE)
  }

  umap_key     <- if ("umap" %in% names(hep@reductions)) "umap" else "UMAP"
  hep_umap_mat <- Embeddings(hep, reduction = umap_key)[, 1:2]
  hep_umap_df  <- as.data.frame(hep_umap_mat) %>%
    setNames(c("UMAP1", "UMAP2")) %>%
    tibble::rownames_to_column("cell")

  # Join UMAP with NNLS usages
  hep_usage_df <- inner_join(hep_umap_df, nnls_usage, by = "cell")
  message("  Hepatocyte cells with UMAP + NNLS usage: ", nrow(hep_usage_df))

  if (nrow(hep_usage_df) == 0) {
    message("  WARNING: cell barcode join returned 0 rows - check barcode formatting.")
  } else {
    panels_hep <- lapply(top_programs, function(prog) {
      top_tf <- top_tf_per_program %>%
        dplyr::filter(program == prog) %>%
        slice(1) %>%
        pull(tf)
      top_tf <- if (length(top_tf) == 0) "n.s." else top_tf

      prog_col <- as.character(prog)
      if (!prog_col %in% colnames(hep_usage_df)) {
        message("  Program column '", prog_col, "' not in NNLS usage ? skipping.")
        return(NULL)
      }

      df_prog <- hep_usage_df %>%
        dplyr::select(UMAP1, UMAP2, usage = all_of(prog_col))

      make_usage_umap(df_prog, paste0("Program ", prog), top_tf = top_tf)
    })

    panels_hep <- dplyr::filter(Negate(is.null), panels_hep)

    if (length(panels_hep) > 0) {
      p_hep <- assemble_usage_panels(
        panels_hep,
        title    = "Top 5 hepatocyte-enriched cNMF programs - mature hepatocyte reference (NNLS)",
        subtitle = "Colour: NNLS-inferred usage  |  Subtitle: top Fisher's TF regulator"
      )
      save_plot(p_hep,
                file.path(OUT_DIR, "ggplot_umap_hep_nnls_top5_programs.png"),
                width = 20, height = 5)
    }
  }
}

message("\n=== Program enrichment analysis complete ===")
message("Output -> ", OUT_DIR)