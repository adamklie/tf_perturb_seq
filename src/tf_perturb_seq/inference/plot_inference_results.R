#!/usr/bin/env Rscript
# =============================================================================
# plot_inference_results.R
#
# Re-creates all plots from investigate_cis_inference_results.py and
# investigate_trans_inference_results.py using ggplot2 + ggsci iTerm palette.
#
# Usage (source interactively in RStudio, or run from command line):
#   Rscript plot_inference_results.R
#
# Set the two directory paths below before running. All intermediate TSV/CSV
# files produced by the Python scripts are read from those directories.
# Output PNGs are written to the same directories unless overridden.
# =============================================================================

.libPaths("R-bak/x86_64-pc-linux-gnu-library/4.4/")
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(forcats)
  library(scales)
  library(ggsci)       # iTerm palette: scale_color_d3("category10") or pal_d3()
  library(ggrepel)
  library(patchwork)
})

# ---- User-defined paths -----------------------------------------------------
#CIS_DIR   <- "/hpc/group/gersbachlab/seg95/tf_perturb_seq/datasets/Gersbach_WTC11-hepatocyte-differentiation_TF-Perturb-seq/calibrated_outs/cis_inference_visualizations"
#TRANS_DIR <- "/hpc/group/gersbachlab/seg95/tf_perturb_seq/datasets/Gersbach_WTC11-hepatocyte-differentiation_TF-Perturb-seq/calibrated_outs/trans_inference_visualizations"

CIS_DIR   <- "/hpc/group/gersbachlab/seg95/tf_perturb_seq/datasets/Gersbach_WTC11-hepatocyte-differentiation_TF-Perturb-seq/calibrated_outs_manual/cis_inference_visualizations"
TRANS_DIR <- "/hpc/group/gersbachlab/seg95/tf_perturb_seq/datasets/Gersbach_WTC11-hepatocyte-differentiation_TF-Perturb-seq/calibrated_outs_manual/trans_inference_visualizations"

# Output directories 
CIS_OUT   <- CIS_DIR
TRANS_OUT <- TRANS_DIR

# Analysis parameters (should match what was passed to the Python scripts)
P_THRESH  <- 0.05   # BH FDR threshold used in Python (for reference lines)
FC_THRESH <- 0.5    # |log2FC| threshold for volcano threshold lines
TOP_N     <- 20     # number of top hits to show in ranked bar charts
#[1] "#1F77B4FF" "#FF7F0EFF" "#2CA02CFF" "#D62728FF" "#9467BDFF" "#8C564BFF" "#E377C2FF" "#7F7F7FFF"
#[9] "#BCBD22FF" "#17BECFFF"
# ---- Helpers ----------------------------------------------------------------
iterm <- pal_d3("category10")(10)   # iTerm / D3 category10 palette from ggsci

LABEL_COLORS <- c(
  "targeting"        = "#F6C177FF",   # orange
  "positive control" = "#D62728FF",   # red
  "negative control" = "#31748FFF",   # blue
  "non-targeting"    = "grey60"
)

theme_inference <- function(base_size = 11) {
  theme_classic(base_size = base_size) %+replace%
    theme(
      strip.background = element_blank(),
      strip.text       = element_text(face = "bold", size = base_size),
      axis.line        = element_line(colour = "grey30"),
      axis.ticks       = element_line(colour = "grey30", size = 0.5),
      axis.ticks.length = unit(0.2, "cm"),
      axis.text        = element_text(size = base_size + 1),
      axis.title       = element_text(size = base_size + 2),
      plot.title       = element_text(face = "bold", size = base_size + 1,
                                      hjust = 0),
      plot.subtitle    = element_text(size = base_size - 1, colour = "grey40",
                                      hjust = 0),
      legend.key.size  = unit(0.4, "cm"),
      panel.spacing    = unit(0.8, "lines")
    )
}

save_plot <- function(p, path, width = 8, height = 6, dpi = 200) {
  ggsave(path, p, width = width, height = height, dpi = dpi,
         bg = "white")
  message("  Saved: ", path)
}

read_tsv_safe <- function(path, ...) {
  if (is.null(path) || !file.exists(path)) {
    message("  SKIP (not found): ", path)
    return(NULL)
  }
  read_tsv(path, show_col_types = FALSE, ...)
}

# list.files() with fixed=TRUE-style glob that handles spaces in filenames
list_tsv <- function(dir, prefix, suffix = ".tsv") {
  all_files <- list.files(dir, full.names = TRUE)
  all_files[startsWith(basename(all_files), prefix) &
              endsWith(basename(all_files), suffix)]
}

# =============================================================================
# CIS PLOTS
# =============================================================================

# --- 1. Guide efficiency (only present if per-guide file existed) -------------
message("\n--- CIS: Guide efficiency ---")

ge_path <- file.path(CIS_DIR, "guide_efficiency_all_methods.tsv")
ge      <- read_tsv_safe(ge_path)

if (!is.null(ge)) {
  means <- ge %>% group_by(method) %>% summarise(mean_sig = mean(sig_guides))
  
  p_ge_hist <- ge %>%
    ggplot(aes(x = sig_guides)) +
    geom_histogram(binwidth = 1, fill = iterm[1], colour = "white",
                   linewidth = 0.3) +
    geom_vline(data = means, aes(xintercept = mean_sig),
               colour = iterm[2], linetype = "dashed", linewidth = 0.8) +
    facet_wrap(~method, scales = "free_y") +
    labs(
      title    = "Guide efficiency: significant cis knockdown per target",
      subtitle = paste0("BH FDR < ", P_THRESH),
      x        = "Number of significant guides per target",
      y        = "Number of targets"
    ) +
    theme_inference()
  
  save_plot(p_ge_hist,
            file.path(CIS_OUT, "ggplot_guide_efficiency_histogram.png"),
            width = 9, height = 5)
  
  p_ge_dot <- ge %>%
    arrange(desc(sig_guides)) %>%
    slice_head(n = 40) %>%
    mutate(target_symbol = fct_reorder(target_symbol, sig_guides)) %>%
    ggplot(aes(x = sig_guides, y = target_symbol, colour = frac_sig)) +
    geom_point(size = 3) +
    scale_colour_gradient(low = "grey80", high = iterm[1],
                          name  = "Fraction\nsignificant",
                          labels = percent_format(accuracy = 1)) +
    facet_wrap(~method, scales = "free_x") +
    labs(title = "Top targets by number of significant guides",
         x = "Significant guides", y = NULL) +
    theme_inference()
  
  save_plot(p_ge_dot,
            file.path(CIS_OUT, "ggplot_guide_efficiency_dotplot.png"),
            width = 9, height = 8)
} else {
  message("  Guide efficiency skipped (no per-guide file for this dataset)")
}

# --- 2. Top downregulated TFs ------------------------------------------------
message("\n--- CIS: Top downregulated TFs ---")

top_tfs <- read_tsv_safe(file.path(CIS_DIR, "top_downregulated_tfs_all_methods.tsv"))

if (!is.null(top_tfs)) {
  # The lfc column is standardised to "log2fc" in the combined file; fall back
  # to whatever numeric column is present if renamed differently
  lfc_col <- intersect(c("log2fc", "log2_fc", "sceptre_log2_fc",
                         "perturbo_log2_fc"), names(top_tfs))[1]
  
  p_top_tfs <- top_tfs %>%
    mutate(symbol = fct_reorder(symbol, neg_log10_fdr)) %>%
    ggplot(aes(x = neg_log10_fdr, y = symbol, fill = .data[[lfc_col]])) +
    geom_col(colour = "white", linewidth = 0.3) +
    scale_fill_gradient2(
      low = "#9CCFD8FF", high = "white",
      midpoint = 0, name = "log2 FC"
    ) +
    facet_wrap(~method, scales = "free") +
    labs(
      title    = "Top downregulated TFs (cis knockdown)",
      subtitle = paste0("BH FDR < ", P_THRESH, ", ranked by FDR"),
      x        = expression(-log[10](BH~FDR)),
      y        = NULL
    ) +
    theme_inference()
  
  n_methods <- n_distinct(top_tfs$method)
  save_plot(p_top_tfs,
            file.path(CIS_OUT, "ggplot_top_downregulated_tfs.png"),
            width = 6 * n_methods, height = max(6, TOP_N * 0.35))
}

# --- 3. All significant cis hits: lollipop -----------------------------------
message("\n--- CIS: All significant hits ---")

cis_sig <- read_tsv_safe(file.path(CIS_DIR, "cis_significant_all_methods.tsv"))

if (!is.null(cis_sig)) {
  plot_data <- cis_sig %>%
    group_by(method) %>%
    slice_max(order_by = neg_log10_fdr, n = TOP_N) %>%
    ungroup() %>%
    mutate(symbol = fct_reorder(symbol, neg_log10_fdr))
  
  p_lollipop <- plot_data %>%
    ggplot(aes(x = neg_log10_fdr, y = symbol, colour = log2fc)) +
    geom_segment(aes(x = 0, xend = neg_log10_fdr, yend = symbol),
                 colour = "grey75", linewidth = 0.4) +
    geom_point(size = 2.5) +
    scale_colour_gradient2(
      low = iterm[2], mid = "grey85", high = iterm[1],
      midpoint = 0, name = "log2 FC"
    ) +
    facet_wrap(~method, scales = "free") +
    labs(
      title = paste0("Top ", TOP_N, " significant cis hits per method"),
      x     = expression(-log[10](BH~FDR)),
      y     = NULL
    ) +
    theme_inference()
  
  save_plot(p_lollipop,
            file.path(CIS_OUT, "ggplot_cis_significant_lollipop.png"),
            width = 7 * n_distinct(cis_sig$method), height = 8)
}

# --- 4. UMAP: hepatocyte markers ---------------------------------------------
message("\n--- CIS: UMAP hepatocyte markers ---")

umap_long <- read_tsv_safe(
  file.path(CIS_DIR, "umap_marker_expression_long.tsv")
)

if (!is.null(umap_long)) {
  umap_long <- umap_long %>%
    group_by(marker) %>%
    arrange(log_norm_expr) %>%
    ungroup()
  
  p_umap <- umap_long %>%
    ggplot(aes(x = UMAP1, y = UMAP2, colour = log_norm_expr)) +
    geom_point(size = 0.4, alpha = 0.6, stroke = 0, shape = 16) +
    scale_colour_gradientn(
      colours = c("grey92", iterm[5], iterm[1], iterm[2]),
      name    = "log-norm\nexpr",
      guide   = guide_colorbar(barheight = 4)
    ) +
    facet_wrap(~marker, ncol = 3) +
    labs(
      title    = "Hepatocyte marker expression",
      subtitle = paste0("n = ",
                        format(n_distinct(umap_long$cell_barcode), big.mark = ","),
                        " cells")
    ) +
    theme_inference() +
    theme(
      axis.text  = element_blank(),
      axis.ticks = element_blank(),
      axis.line  = element_blank(),
      axis.title = element_text(size = 8, colour = "grey50")
    )
  
  nrows <- ceiling(n_distinct(umap_long$marker) / 3)
  save_plot(p_umap,
            file.path(CIS_OUT, "ggplot_umap_hepatocyte_markers.png"),
            width = 12, height = 5 * nrows)
}

# --- 5 & 6. Volcano helpers (shared by cis and trans) -------------------------
plot_label_volcano <- function(df, title, p_thresh = 0.01, fc_thresh = 0.5) {
  df <- df %>%
    mutate(
      neg_log10p = pmin(-log10(pmax(p_value, 1e-350)), 350),
      direction  = case_when(
        p_value < p_thresh & log2fc >  fc_thresh ~ "up",
        p_value < p_thresh & log2fc < -fc_thresh ~ "down",
        TRUE ~ "ns"
      ),
      type = replace_na(as.character(type), "unlabelled")
    )
  
  type_counts <- df %>%
    dplyr::filter(type != "unlabelled") %>%
    count(type, direction) %>%
    pivot_wider(names_from = direction, values_from = n, values_fill = 0L) %>%
    mutate(
      up   = dplyr::coalesce(up,   0L),
      down = dplyr::coalesce(down, 0L),
      ns   = dplyr::coalesce(ns,   0L),
      label = paste0(type, "  (Up: ", up, ", Down:", down, ", NS:", ns, ")")
    )
  label_map <- setNames(type_counts$label, type_counts$type)
  
  present_labels <- label_map[intersect(names(LABEL_COLORS), names(label_map))]
  colour_map     <- setNames(LABEL_COLORS[names(present_labels)],
                             present_labels)
  
  df <- df %>%
    mutate(type_label = if_else(type == "unlabelled",
                                NA_character_, label_map[type]))
  
  # Render layers in LABEL_ORDER so positive controls are always drawn last
  # (on top of targeting and other groups).
  layer_order <- c("targeting", "negative control", "non-targeting",
                   "positive control")
  
  p <- ggplot() +
    geom_point(
      data = dplyr::filter(df, type == "unlabelled"),
      aes(x = log2fc, y = neg_log10p),
      colour = "grey88", size = 0.8, alpha = 0.3, stroke = 0, shape = 16
    )
  for (.type in layer_order) {
    .sub <- dplyr::filter(df, type == .type)
    if (nrow(.sub) == 0) next
    p <- p + geom_point(
      data = .sub,
      aes(x = log2fc, y = neg_log10p, colour = type_label),
      size = 1.8, alpha = 0.75, stroke = 0, shape = 16
    )
  }
  p +
    geom_hline(yintercept = -log10(p_thresh),
               linetype = "dashed", colour = "grey40", linewidth = 0.6) +
    geom_vline(xintercept = c(-fc_thresh, fc_thresh),
               linetype = "dashed", colour = "grey40", linewidth = 0.6) +
    geom_vline(xintercept = 0, colour = "grey30", linewidth = 0.4) +
    scale_colour_manual(values = colour_map, na.value = "grey88",
                        name = NULL, na.translate = FALSE) +
    labs(
      title    = title,
      subtitle = paste0("nominal p < ", p_thresh,
                        "  |  |log2FC| > ", fc_thresh),
      x        = expression(log[2]~"Fold Change"),
      y        = expression(-log[10](p))
    ) +
    theme_inference() +
    theme(legend.position = "right", legend.text = element_text(size = 8))
}


plot_window_volcano <- function(df, title, p_thresh = 0.01, fc_thresh = 0.5) {
  df <- df %>%
    mutate(
      neg_log10p     = pmin(-log10(pmax(p_value, 1e-350)), 350),
      direction      = case_when(
        p_value < p_thresh & log2fc >  fc_thresh ~ "up",
        p_value < p_thresh & log2fc < -fc_thresh ~ "down",
        TRUE ~ "ns"
      ),
      is_self_target = as.logical(is_self_target),
      type           = replace_na(as.character(type), "unlabelled")
    )
  
  bg   <- dplyr::filter(df, !is_self_target)
  self <- dplyr::filter(df, is_self_target)
  
  # Build legend label for background "other cis pairs" group
  bg_label <- paste0(
    "other cis pairs  (Up: ", sum(bg$direction == "up"),
    " Up: ", sum(bg$direction == "down"),
    " Down:", sum(bg$direction == "ns"), ")"
  )
  
  # Build legend labels for self-target groups
  type_counts <- self %>%
    dplyr::filter(type != "unlabelled") %>%
    count(type, direction) %>%
    pivot_wider(names_from = direction, values_from = n, values_fill = 0L) %>%
    mutate(
      up    = dplyr::coalesce(up,   0L),
      down  = dplyr::coalesce(down, 0L),
      ns    = dplyr::coalesce(ns,   0L),
      label = paste0(type, " (self)  (Up: ", up, ", Down:", down, ", NS:", ns, ")")
    )
  label_map  <- setNames(type_counts$label, type_counts$type)
  present    <- label_map[intersect(names(LABEL_COLORS), names(label_map))]
  
  # Include "other cis pairs" as a named colour entry so it appears in the legend
  colour_map <- c(
    setNames(LABEL_COLORS[names(present)], present),
    setNames("#2ca02c", bg_label)
  )
  
  # Tag background rows with the legend label
  bg   <- mutate(bg,   type_label = bg_label)
  self <- self %>%
    mutate(type_label = if_else(type == "unlabelled",
                                NA_character_, label_map[type]))
  
  # Plot background first (smaller, more transparent), self-target on top
  ggplot() +
    geom_point(data = bg,
               aes(x = log2fc, y = neg_log10p, colour = type_label),
               size = 1.4, alpha = 0.25, stroke = 0, shape = 16) +
    geom_point(data = self,
               aes(x = log2fc, y = neg_log10p, colour = type_label),
               size = 2, alpha = 0.85, stroke = 0, shape = 16) +
    geom_hline(yintercept = -log10(p_thresh),
               linetype = "dashed", colour = "grey40", linewidth = 0.6) +
    geom_vline(xintercept = c(-fc_thresh, fc_thresh),
               linetype = "dashed", colour = "grey40", linewidth = 0.6) +
    geom_vline(xintercept = 0, colour = "grey30", linewidth = 0.4) +
    scale_colour_manual(values = colour_map, na.value = "grey88",
                        name = NULL, na.translate = FALSE) +
    labs(
      title    = title,
      subtitle = paste0("nominal p < ", p_thresh,
                        "  |  |log2FC| > ", fc_thresh),
      x        = expression(log[2]~"Fold Change"),
      y        = expression(-log[10](p))
    ) +
    theme_inference() +
    theme(legend.position = "right", legend.text = element_text(size = 8))
}

# --- 5. Cis label volcanos ---------------------------------------------------
message("\n--- CIS: Label volcano plots ---")

for (f in list_tsv(CIS_DIR, "volcano_calibrated_element_")) {
  df_v <- read_tsv_safe(f)
  if (is.null(df_v)) next
  print(head(df_v))
  tag   <- sub("^volcano_calibrated_element_(.+)\\.tsv$", "\\1", basename(f))
  title <- paste0("Cis volcano  (", gsub("_", " ", tag), ")")
  p_v   <- plot_label_volcano(df_v, title)
  print(p_v)
  save_plot(p_v,
            file.path(CIS_OUT, paste0("ggplot_volcano_cis_element_", tag, ".png")),
            width = 8, height = 6)
}

# Also handle the non-"calibrated" element volcanos if present
for (f in list_tsv(CIS_DIR, "volcano_cis_element_")) {
  # Only pick up the TSV versions that weren't already covered above
  tag <- sub("^volcano_cis_element_(.+)\\.tsv$", "\\1", basename(f))
  if (file.exists(file.path(CIS_DIR,
                            paste0("ggplot_volcano_cis_element_", tag, ".png")))) next
  df_v <- read_tsv_safe(f)
  if (is.null(df_v) || !all(c("p_value", "log2fc") %in% names(df_v))) next
  title <- paste0("Cis volcano  (", gsub("_", " ", tag), ")")
  p_v   <- plot_label_volcano(df_v, title)
  print(p_v)
  save_plot(p_v,
            file.path(CIS_OUT, paste0("ggplot_volcano_cis_element_", tag, ".png")),
            width = 8, height = 6)
}

# --- 6. Cis window volcanos --------------------------------------------------
message("\n--- CIS: Window volcano plots ---")

for (f in list_tsv(CIS_DIR, "volcano_cis_window_")) {
  df_w <- read_tsv_safe(f)
  if (is.null(df_w) || !("is_self_target" %in% names(df_w))) next
  tag   <- sub("^volcano_cis_window_(.+)\\.tsv$", "\\1", basename(f))
  title <- paste0("Cis window volcano  (", gsub("_", " ", tag), ")")
  p_w   <- plot_window_volcano(df_w, title)
  print(p_w)
  save_plot(p_w,
            file.path(CIS_OUT, paste0("ggplot_volcano_cis_window_", tag, ".png")),
            width = 9, height = 6)
}


# =============================================================================
# TRANS PLOTS
# =============================================================================

# --- 7. Per-element trans distribution ---------------------------------------
message("\n--- TRANS: Per-element distribution ---")

trans_elem <- read_tsv_safe(
  file.path(TRANS_DIR, "trans_per_element_all_methods.tsv")
)

if (!is.null(trans_elem)) {
  means_trans <- trans_elem %>%
    group_by(method) %>%
    summarise(mean_sig = mean(n_trans_sig))
  
  p_te_hist <- trans_elem %>%
    ggplot(aes(x = n_trans_sig)) +
    geom_histogram(bins = 60, fill = iterm[1], colour = "white",
                   linewidth = 0.3) +
    geom_vline(data = means_trans, aes(xintercept = mean_sig),
               colour = iterm[2], linetype = "dashed", linewidth = 0.8) +
    facet_wrap(~method, scales = "free") +
    labs(
      title    = "Trans hits per targeting element",
      subtitle = paste0("BH FDR < ", P_THRESH),
      x        = "Significant trans genes",
      y        = "Number of elements"
    ) +
    theme_inference()
  
  save_plot(p_te_hist,
            file.path(TRANS_OUT, "ggplot_trans_per_element_distribution.png"),
            width = 8, height = 5)
}

# --- 8. Top trans elements ---------------------------------------------------
message("\n--- TRANS: Top elements ---")

for (f in list_tsv(TRANS_DIR, "trans_top_elements_")) {
  df_te <- read_tsv_safe(f)
  if (is.null(df_te)) next
  tag <- sub("^trans_top_elements_(.+)\\.tsv$", "\\1", basename(f))
  
  p_te <- df_te %>%
    slice_max(order_by = n_trans_sig, n = TOP_N) %>%
    mutate(target_symbol = fct_reorder(target_symbol, n_trans_sig)) %>%
    ggplot(aes(x = n_trans_sig, y = target_symbol)) +
    geom_col(fill = iterm[1], colour = "white", linewidth = 0.3) +
    labs(
      title    = paste0("Top ", TOP_N, " targeting elements by trans hits"),
      subtitle = paste0(tag, "  |  BH FDR < ", P_THRESH),
      x        = "Significant trans genes",
      y        = NULL
    ) +
    theme_inference()
  
  save_plot(p_te,
            file.path(TRANS_OUT,
                      paste0("ggplot_trans_top_elements_", tag, ".png")),
            width = 7, height = max(4, TOP_N * 0.35))
}

# --- 9. Top trans genes ------------------------------------------------------
message("\n--- TRANS: Top genes ---")

for (f in list_tsv(TRANS_DIR, "trans_top_genes_")) {
  df_tg <- read_tsv_safe(f)
  if (is.null(df_tg)) next
  tag     <- sub("^trans_top_genes_(.+)\\.tsv$", "\\1", basename(f))
  #lfc_abs <- max(abs(range(df_tg$median_lfc, na.rm = TRUE)))
  lfc_abs <- quantile(abs(df_tg$median_lfc), 0.95, na.rm = TRUE)
  
  
  p_tg <- df_tg %>%
    slice_max(order_by = n_targeting_elements, n = TOP_N) %>%
    mutate(gene_symbol = fct_reorder(gene_symbol, n_targeting_elements)) %>%
    ggplot(aes(x = n_targeting_elements, y = gene_symbol, fill = median_lfc)) +
    geom_col(colour = "white", linewidth = 0.3) +
    scale_fill_gradient2(
      low = iterm[1], mid = "grey90", high = iterm[2],
      midpoint = 0, limits = c(-lfc_abs, lfc_abs),
      oob  = squish,
      name = "Median\nlog2 FC"
    ) +
    labs(
      title    = paste0("Top ", TOP_N, " trans-regulated genes"),
      subtitle = paste0(tag, "  |  ranked by breadth across targeting elements"),
      x        = "Targeting elements for which gene is significant",
      y        = NULL
    ) +
    theme_inference()
  
  save_plot(p_tg,
            file.path(TRANS_OUT,
                      paste0("ggplot_trans_top_genes_", tag, ".png")),
            width = 8, height = max(4, TOP_N * 0.35))
}

# --- 10. Trans volcano -------------------------------------------------------
message("\n--- TRANS: Volcano plots ---")

for (f in list_tsv(TRANS_DIR, "volcano_trans_element_")) {
  df_tv <- read_tsv_safe(f)
  if (is.null(df_tv) || !all(c("p_value", "log2fc") %in% names(df_tv))) next
  tag   <- sub("^volcano_trans_element_(.+)\\.tsv$", "\\1", basename(f))
  title <- paste0("Trans volcano  (", gsub("_", " ", tag), ")")
  p_tv  <- plot_label_volcano(df_tv, title, p_thresh = 0.01,
                              fc_thresh = FC_THRESH)
  save_plot(p_tv,
            file.path(TRANS_OUT,
                      paste0("ggplot_volcano_trans_element_", tag, ".png")),
            width = 9, height = 6)
}

message("\n=== All plots complete ===")
message("CIS plots   -> ", CIS_OUT)
message("TRANS plots -> ", TRANS_OUT)