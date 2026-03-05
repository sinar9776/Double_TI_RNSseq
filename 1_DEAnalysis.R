#============================================================
# 1_DEAnalysis.R — Modular DE (gene/isoform)
# - RAW + SHRUNK tables (tT_*, up_*, down_*)
# - Adds extra contrasts:
#   (1) drug_OHT vs drug        => tT_<drug>_OHT__<drug>
#   (2) drug+OHT vs DMSO (combined) => tT_<drug>_OHT_vs_DMSO
# - Figures: PDF + TIFF only (NO png)
#============================================================

suppressPackageStartupMessages({
  library(BiocParallel)
  library(parallel)
  library(tidyverse)
  library(DESeq2)
  library(apeglm)
  library(ggplot2)
})

set.seed(1)

#========================
# Parameters
#========================
analysis_level <- "gene"  # "gene" or "isoform"
n_cores <- 1L
padj_cutoff <- 0.01
lfc_cutoff  <- 0
shrink_type <- "apeglm"

#========================
# Parallel backend (RStudio-safe)
#========================
register(SnowParam(workers = n_cores, type = "SOCK"))
options(mc.cores = n_cores)

#========================
# Helpers
#========================
strip_version <- function(x) sub("\\..*$", "", x)

save_plot <- function(filename_base, plot_obj, width = 6, height = 5) {
  ggsave(paste0(filename_base, ".pdf"), plot_obj, width = width, height = height)
  ggsave(paste0(filename_base, ".tiff"), plot_obj, width = width, height = height,
         dpi = 600, device = "tiff", compression = "lzw")
}

calc_volcano_limits <- function(files, q = 0.995) {
  if (length(files) == 0) return(list(x_lim = NULL, y_lim = NULL))
  all_lfc <- c(); all_y <- c()
  
  for (fp in files) {
    df <- tryCatch(read.delim(fp, check.names = FALSE), error = function(e) NULL)
    if (is.null(df)) next
    if (!all(c("log2FoldChange", "padj") %in% names(df))) next
    
    lfc <- df$log2FoldChange
    p   <- df$padj
    ok <- is.finite(lfc) & !is.na(lfc) & is.finite(p) & !is.na(p)
    if (!any(ok)) next
    
    all_lfc <- c(all_lfc, lfc[ok])
    all_y   <- c(all_y, -log10(pmax(p[ok], .Machine$double.xmin)))
  }
  
  if (length(all_lfc) == 0 || length(all_y) == 0) return(list(x_lim = NULL, y_lim = NULL))
  
  x_max <- as.numeric(quantile(abs(all_lfc), probs = q, na.rm = TRUE))
  if (!is.finite(x_max) || x_max <= 0) x_max <- max(abs(all_lfc), na.rm = TRUE)
  
  y_ok <- all_y[is.finite(all_y)]
  y_max <- as.numeric(quantile(y_ok, probs = q, na.rm = TRUE))
  if (!is.finite(y_max) || y_max <= 0) y_max <- max(y_ok, na.rm = TRUE)
  
  list(x_lim = c(-x_max, x_max), y_lim = c(0, y_max))
}

#========================
# Paths
#========================
setwd()
base_dir <- getwd()
data_dir <- file.path(base_dir, "data")
res_dir  <- file.path(base_dir, "results/1_DEAnalysis")

tbl_dir  <- file.path(res_dir, "tables/raw")
fig_dir  <- file.path(res_dir, "figures/raw")
shr_tbl  <- file.path(res_dir, "tables/shrunk")
shr_fig  <- file.path(res_dir, "figures/shrunk")

dir.create(tbl_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(shr_tbl, recursive = TRUE, showWarnings = FALSE)
dir.create(shr_fig, recursive = TRUE, showWarnings = FALSE)

#========================
# Load data
#========================
metadata <- read.delim(file.path(base_dir, "metadata.txt"), check.names = FALSE)
stopifnot(all(c("SampleID", "group") %in% colnames(metadata)))

annotation <- read.delim(file.path(base_dir, "annotation.txt"), check.names = FALSE)

if (analysis_level == "isoform") {
  df <- readr::read_tsv(file.path(data_dir, "isoform_counts.tsv"), show_col_types = FALSE) %>% as.data.frame()
  id_col <- "isoform_id"
} else if (analysis_level == "gene") {
  df <- readr::read_tsv(file.path(base_dir, "gene_counts.tsv"), show_col_types = FALSE) %>% as.data.frame()
  id_col <- "gene_id"
} else {
  stop("analysis_level must be 'gene' or 'isoform'")
}

rownames(df) <- df[[1]]
count_mat <- df[, intersect(colnames(df), metadata$SampleID), drop = FALSE]
metadata  <- metadata[match(colnames(count_mat), metadata$SampleID), ]
stopifnot(all(colnames(count_mat) == metadata$SampleID))

#========================
# Annotation (version-stripped IDs for safe merges)
#========================
sym_col <- if ("symbol" %in% names(annotation)) "symbol" else if ("gene_name" %in% names(annotation)) "gene_name" else NA

ann_keep_cols <- unique(c("gene_id", id_col, sym_col))
ann_keep_cols <- ann_keep_cols[ann_keep_cols %in% names(annotation)]
annotation2 <- annotation[, ann_keep_cols, drop = FALSE]
if (!is.na(sym_col)) names(annotation2)[names(annotation2) == sym_col] <- "SYMBOL"

if ("gene_id" %in% names(annotation2))    annotation2$gene_id    <- strip_version(annotation2$gene_id)
if ("isoform_id" %in% names(annotation2)) annotation2$isoform_id <- strip_version(annotation2$isoform_id)

if (analysis_level == "gene") {
  annotation2 <- annotation2 %>%
    group_by(gene_id) %>%
    summarise(SYMBOL = dplyr::first(na.omit(SYMBOL)), .groups = "drop")
}

if (!"gene_id" %in% colnames(annotation2)) {
  stop("annotation.txt must include a 'gene_id' column for downstream compatibility.")
}

#========================
# Filter low counts
#========================
keep <- rowSums(count_mat >= 10) >= min(table(metadata$group))
count_mat <- count_mat[keep, , drop = FALSE]

#========================
# DESeq objects
#========================
dds1 <- DESeqDataSetFromMatrix(round(as.matrix(count_mat)), metadata, design = ~ group)
dds1$group <- relevel(dds1$group, ref = "DMSO")

dds2 <- dds1
if ("DMSO_OHT" %in% levels(dds2$group)) dds2$group <- relevel(dds2$group, ref = "DMSO_OHT")

dds1 <- DESeq(dds1, parallel = TRUE)
dds2 <- DESeq(dds2, parallel = TRUE)

#========================
# RAW DE (default contrasts)
#========================
counts_summary_raw <- data.frame()

for (g in levels(dds1$group)) {
  
  ref <- if (grepl("_OHT$", g)) "DMSO_OHT" else "DMSO"
  if (g == "DMSO_OHT") ref <- "DMSO"
  if (g == ref || !(ref %in% levels(dds1$group))) next
  
  nm <- if (g == "DMSO_OHT") "DMSO_OHT_vs_DMSO" else g
  message(" [RAW] ", nm)
  
  res_raw <- results(dds1, contrast = c("group", g, ref))
  res_df  <- as.data.frame(res_raw)
  res_df[[id_col]] <- strip_version(rownames(res_df))
  res_df <- merge(res_df, annotation2, by = id_col, all.x = TRUE, sort = FALSE)
  
  nc1 <- counts(dds1, normalized = TRUE)
  idx_t <- which(metadata$group == g)
  idx_c <- which(metadata$group == ref)
  res_df$avg_treated <- rowMeans(nc1[, idx_t, drop = FALSE], na.rm = TRUE)
  res_df$avg_control <- rowMeans(nc1[, idx_c, drop = FALSE], na.rm = TRUE)
  
  res_df$sig_raw <- ifelse(!is.na(res_df$padj) & res_df$padj < padj_cutoff & res_df$log2FoldChange >  lfc_cutoff, "up",
                           ifelse(!is.na(res_df$padj) & res_df$padj < padj_cutoff & res_df$log2FoldChange < -lfc_cutoff, "down", "ns"))
  
  write.table(res_df, file.path(tbl_dir, paste0("tT_", nm, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(subset(res_df, sig_raw == "up"),   file.path(tbl_dir, paste0("up_", nm, ".tsv")),   sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(subset(res_df, sig_raw == "down"), file.path(tbl_dir, paste0("down_", nm, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
  
  counts_summary_raw <- rbind(counts_summary_raw, data.frame(
    contrast = nm, n_up = sum(res_df$sig_raw == "up"), n_down = sum(res_df$sig_raw == "down")
  ))
}

#========================
# RAW EXTRA 1: drug_OHT vs drug  => tT_<drug>_OHT__<drug>
#========================
all_groups <- levels(dds1$group)
drug_groups <- setdiff(all_groups[!grepl("_OHT$", all_groups)], c("DMSO", "DMSO_OHT"))

for (drug in drug_groups) {
  g_oht <- paste0(drug, "_OHT")
  if (!(g_oht %in% all_groups)) next
  
  nm <- paste0(drug, "_OHT__", drug)
  message(" [RAW] ", nm)
  
  res_raw <- results(dds1, contrast = c("group", g_oht, drug))
  res_df  <- as.data.frame(res_raw)
  res_df[[id_col]] <- strip_version(rownames(res_df))
  res_df <- merge(res_df, annotation2, by = id_col, all.x = TRUE, sort = FALSE)
  
  nc1 <- counts(dds1, normalized = TRUE)
  idx_t <- which(metadata$group == g_oht)
  idx_c <- which(metadata$group == drug)
  res_df$avg_treated <- rowMeans(nc1[, idx_t, drop = FALSE], na.rm = TRUE)
  res_df$avg_control <- rowMeans(nc1[, idx_c, drop = FALSE], na.rm = TRUE)
  
  res_df$sig_raw <- ifelse(!is.na(res_df$padj) & res_df$padj < padj_cutoff & res_df$log2FoldChange >  lfc_cutoff, "up",
                           ifelse(!is.na(res_df$padj) & res_df$padj < padj_cutoff & res_df$log2FoldChange < -lfc_cutoff, "down", "ns"))
  
  write.table(res_df, file.path(tbl_dir, paste0("tT_", nm, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(subset(res_df, sig_raw == "up"),   file.path(tbl_dir, paste0("up_", nm, ".tsv")),   sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(subset(res_df, sig_raw == "down"), file.path(tbl_dir, paste0("down_", nm, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
  
  counts_summary_raw <- rbind(counts_summary_raw, data.frame(
    contrast = nm, n_up = sum(res_df$sig_raw == "up"), n_down = sum(res_df$sig_raw == "down")
  ))
}

#========================
# RAW EXTRA 2: drug+OHT vs DMSO (combined effect) => tT_<drug>_OHT_vs_DMSO
#========================
for (drug in drug_groups) {
  g_oht <- paste0(drug, "_OHT")
  if (!(g_oht %in% all_groups)) next
  if (!("DMSO" %in% all_groups)) next
  
  nm <- paste0(drug, "_OHT_vs_DMSO")
  message(" [RAW] ", nm)
  
  res_raw <- results(dds1, contrast = c("group", g_oht, "DMSO"))
  res_df  <- as.data.frame(res_raw)
  res_df[[id_col]] <- strip_version(rownames(res_df))
  res_df <- merge(res_df, annotation2, by = id_col, all.x = TRUE, sort = FALSE)
  
  nc1 <- counts(dds1, normalized = TRUE)
  idx_t <- which(metadata$group == g_oht)
  idx_c <- which(metadata$group == "DMSO")
  res_df$avg_treated <- rowMeans(nc1[, idx_t, drop = FALSE], na.rm = TRUE)
  res_df$avg_control <- rowMeans(nc1[, idx_c, drop = FALSE], na.rm = TRUE)
  
  res_df$sig_raw <- ifelse(!is.na(res_df$padj) & res_df$padj < padj_cutoff & res_df$log2FoldChange >  lfc_cutoff, "up",
                           ifelse(!is.na(res_df$padj) & res_df$padj < padj_cutoff & res_df$log2FoldChange < -lfc_cutoff, "down", "ns"))
  
  write.table(res_df, file.path(tbl_dir, paste0("tT_", nm, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(subset(res_df, sig_raw == "up"),   file.path(tbl_dir, paste0("up_", nm, ".tsv")),   sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(subset(res_df, sig_raw == "down"), file.path(tbl_dir, paste0("down_", nm, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
  
  counts_summary_raw <- rbind(counts_summary_raw, data.frame(
    contrast = nm, n_up = sum(res_df$sig_raw == "up"), n_down = sum(res_df$sig_raw == "down")
  ))
}

write.table(counts_summary_raw, file.path(tbl_dir, "DE_counts_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

#========================
# RAW plots (Dispersion + PCA + MA/Volcano)
#========================
pdf(file.path(fig_dir, "Dispersion_dds1.pdf"), width = 7, height = 6); plotDispEsts(dds1); dev.off()
tiff(file.path(fig_dir, "Dispersion_dds1.tiff"), width = 7, height = 6, units = "in", res = 600, compression = "lzw"); plotDispEsts(dds1); dev.off()

pdf(file.path(fig_dir, "Dispersion_dds2.pdf"), width = 7, height = 6); plotDispEsts(dds2); dev.off()
tiff(file.path(fig_dir, "Dispersion_dds2.tiff"), width = 7, height = 6, units = "in", res = 600, compression = "lzw"); plotDispEsts(dds2); dev.off()

#---------- PCA per comparison (raw design: each g vs its ref) ----------
pca_colors <- c("treated" = "#D62728", "control" = "#1F77B4")

for (g in levels(dds1$group)) {
  
  ref <- if (grepl("_OHT$", g)) "DMSO_OHT" else "DMSO"
  if (g == "DMSO_OHT") ref <- "DMSO"
  if (g == ref || !(ref %in% levels(dds1$group))) next
  
  keep_samples <- metadata$group %in% c(ref, g)
  sub_metadata <- metadata[keep_samples, ]
  sub_counts   <- count_mat[, sub_metadata$SampleID, drop = FALSE]
  
  dds_sub <- DESeqDataSetFromMatrix(round(sub_counts), sub_metadata, design = ~ group)
  dds_sub$group <- droplevels(dds_sub$group)
  vsd_sub <- vst(dds_sub, blind = TRUE)
  
  pca_data <- plotPCA(vsd_sub, intgroup = "group", returnData = TRUE)
  percentVar <- round(100 * attr(pca_data, "percentVar"))
  
  pca_data$group <- factor(pca_data$group, levels = c(ref, g))
  color_map <- setNames(c(pca_colors["control"], pca_colors["treated"]), c(ref, g))
  
  p_pca <- ggplot(pca_data, aes(PC1, PC2, color = group)) +
    geom_point(size = 0.6, alpha = 0.8) +
    scale_color_manual(values = color_map, name = "Condition", labels = c(ref, g)) +
    labs(x = paste0("PC1 (", percentVar[1], "%)"),
         y = paste0("PC2 (", percentVar[2], "%)"),
         title = paste("PCA:", g, "vs", ref)) +
    theme_bw(base_size = 10) +
    theme(legend.position = "right")
  
  save_plot(file.path(fig_dir, paste0("PCA_", g, "_vs_", ref)), p_pca, width = 6, height = 5)
}

#---------- MA & Volcano (RAW) ----------
col_scale <- c(up = "#D62728", down = "#1F77B4", ns = "gray80")
files_raw <- list.files(tbl_dir, pattern = "^tT_.*\\.tsv$", full.names = TRUE)

lim_raw <- calc_volcano_limits(files_raw, q = 0.995)
x_lim_raw <- lim_raw$x_lim
y_lim_raw <- lim_raw$y_lim

for (fp in files_raw) {
  
  nm_raw <- tools::file_path_sans_ext(basename(fp))
  nm <- sub("^tT_", "", nm_raw)
  
  res_df <- tryCatch(read.delim(fp, check.names = FALSE), error = function(e) NULL)
  if (is.null(res_df) || nrow(res_df) == 0) next
  if (!all(c("baseMean", "log2FoldChange", "padj") %in% names(res_df))) next
  
  if (!"sig_raw" %in% names(res_df)) {
    res_df$sig_raw <- ifelse(!is.na(res_df$padj) & res_df$padj < padj_cutoff & res_df$log2FoldChange >  lfc_cutoff, "up",
                             ifelse(!is.na(res_df$padj) & res_df$padj < padj_cutoff & res_df$log2FoldChange < -lfc_cutoff, "down", "ns"))
  }
  
  res_df$cat <- factor(res_df$sig_raw, levels = c("ns", "down", "up"))
  plot_df <- res_df %>% arrange(cat)
  
  p_ma <- ggplot(plot_df, aes(x = log10(baseMean + 1), y = log2FoldChange, color = cat)) +
    geom_point(size = 0.6, alpha = 0.85) +
    scale_color_manual(values = col_scale,
                       breaks = c("up","down","ns"),
                       labels = c("Up","Down","NS"),
                       name = "DE category") +
    geom_hline(yintercept = c(-lfc_cutoff, lfc_cutoff), linetype = "dashed") +
    labs(title = paste("MA:", nm), x = "log10(baseMean + 1)", y = "log2(Fold Change)") +
    theme_bw(10)
  
  p_vol <- ggplot(plot_df, aes(x = log2FoldChange, y = -log10(pmax(padj, .Machine$double.xmin)), color = cat)) +
    geom_point(size = 0.6, alpha = 0.8) +
    scale_color_manual(values = col_scale,
                       breaks = c("up","down","ns"),
                       labels = c("Up","Down","NS"),
                       name = "DE category") +
    geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), linetype = "dashed") +
    geom_hline(yintercept = -log10(padj_cutoff), linetype = "dashed") +
    labs(title = paste("Volcano:", nm), x = "log2(Fold Change)", y = "-log10(padj)") +
    theme_bw(10)
  
  if (!is.null(x_lim_raw) && !is.null(y_lim_raw)) p_vol <- p_vol + coord_cartesian(xlim = x_lim_raw, ylim = y_lim_raw)
  
  save_plot(file.path(fig_dir, paste0("MA_", nm)), p_ma, width = 6, height = 5)
  save_plot(file.path(fig_dir, paste0("Volcano_", nm)), p_vol, width = 6.5, height = 5)
}

#========================
# SHRUNK DE (default contrasts)
#========================
counts_summary_shr <- data.frame()

for (g in levels(dds1$group)) {
  
  ref <- if (grepl("_OHT$", g)) "DMSO_OHT" else "DMSO"
  if (g == "DMSO_OHT") ref <- "DMSO"
  if (g == ref || !(ref %in% levels(dds1$group))) next
  
  nm <- if (g == "DMSO_OHT") "DMSO_OHT_vs_DMSO" else g
  message(" [SHRUNK] ", nm)
  
  dds_use <- if (ref == "DMSO_OHT") dds2 else dds1
  coef_name <- paste0("group_", g, "_vs_", ref)
  if (!(coef_name %in% resultsNames(dds_use))) next
  
  res_raw <- results(dds_use, contrast = c("group", g, ref))
  res_shr <- tryCatch(
    lfcShrink(dds_use, coef = coef_name, res = res_raw, type = shrink_type),
    error = function(e) lfcShrink(dds_use, coef = coef_name, res = res_raw, type = "normal")
  )
  
  res_df <- as.data.frame(res_shr)
  res_df[[id_col]] <- strip_version(rownames(res_df))
  res_df <- merge(res_df, annotation2, by = id_col, all.x = TRUE, sort = FALSE)
  
  nc <- counts(dds_use, normalized = TRUE)
  idx_t <- which(metadata$group == g)
  idx_c <- which(metadata$group == ref)
  res_df$avg_treated <- rowMeans(nc[, idx_t, drop = FALSE], na.rm = TRUE)
  res_df$avg_control <- rowMeans(nc[, idx_c, drop = FALSE], na.rm = TRUE)
  
  res_df$sig_shr <- ifelse(!is.na(res_df$padj) & res_df$padj < padj_cutoff & res_df$log2FoldChange >  lfc_cutoff, "up",
                           ifelse(!is.na(res_df$padj) & res_df$padj < padj_cutoff & res_df$log2FoldChange < -lfc_cutoff, "down", "ns"))
  
  write.table(res_df, file.path(shr_tbl, paste0("tT_", nm, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(subset(res_df, sig_shr == "up"),   file.path(shr_tbl, paste0("up_", nm, ".tsv")),   sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(subset(res_df, sig_shr == "down"), file.path(shr_tbl, paste0("down_", nm, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
  
  counts_summary_shr <- rbind(counts_summary_shr, data.frame(
    contrast = nm, n_up_shr = sum(res_df$sig_shr == "up"), n_down_shr = sum(res_df$sig_shr == "down")
  ))
}

#========================
# SHRUNK EXTRA 1: drug_OHT vs drug  => tT_<drug>_OHT__<drug>
#========================
for (drug in drug_groups) {
  
  g_oht <- paste0(drug, "_OHT")
  if (!(g_oht %in% all_groups)) next
  
  nm <- paste0(drug, "_OHT__", drug)
  message(" [SHRUNK] ", nm)
  
  dds_tmp <- dds1
  dds_tmp$group <- relevel(dds_tmp$group, ref = drug)
  dds_tmp <- DESeq(dds_tmp, parallel = TRUE)
  
  coef_name <- paste0("group_", g_oht, "_vs_", drug)
  if (!(coef_name %in% resultsNames(dds_tmp))) next
  
  res_raw <- results(dds_tmp, contrast = c("group", g_oht, drug))
  res_shr <- tryCatch(
    lfcShrink(dds_tmp, coef = coef_name, res = res_raw, type = shrink_type),
    error = function(e) lfcShrink(dds_tmp, coef = coef_name, res = res_raw, type = "normal")
  )
  
  res_df <- as.data.frame(res_shr)
  res_df[[id_col]] <- strip_version(rownames(res_df))
  res_df <- merge(res_df, annotation2, by = id_col, all.x = TRUE, sort = FALSE)
  
  nc1 <- counts(dds_tmp, normalized = TRUE)
  idx_t <- which(metadata$group == g_oht)
  idx_c <- which(metadata$group == drug)
  res_df$avg_treated <- rowMeans(nc1[, idx_t, drop = FALSE], na.rm = TRUE)
  res_df$avg_control <- rowMeans(nc1[, idx_c, drop = FALSE], na.rm = TRUE)
  
  res_df$sig_shr <- ifelse(!is.na(res_df$padj) & res_df$padj < padj_cutoff & res_df$log2FoldChange >  lfc_cutoff, "up",
                           ifelse(!is.na(res_df$padj) & res_df$padj < padj_cutoff & res_df$log2FoldChange < -lfc_cutoff, "down", "ns"))
  
  write.table(res_df, file.path(shr_tbl, paste0("tT_", nm, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(subset(res_df, sig_shr == "up"),   file.path(shr_tbl, paste0("up_", nm, ".tsv")),   sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(subset(res_df, sig_shr == "down"), file.path(shr_tbl, paste0("down_", nm, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
  
  counts_summary_shr <- rbind(counts_summary_shr, data.frame(
    contrast = nm, n_up_shr = sum(res_df$sig_shr == "up"), n_down_shr = sum(res_df$sig_shr == "down")
  ))
}

#========================
# SHRUNK EXTRA 2: drug+OHT vs DMSO (combined effect) => tT_<drug>_OHT_vs_DMSO
#========================
for (drug in drug_groups) {
  
  g_oht <- paste0(drug, "_OHT")
  if (!(g_oht %in% all_groups)) next
  if (!("DMSO" %in% all_groups)) next
  
  nm <- paste0(drug, "_OHT_vs_DMSO")
  message(" [SHRUNK] ", nm)
  
  # dds1 has ref DMSO, so coef exists directly
  coef_name <- paste0("group_", g_oht, "_vs_DMSO")
  if (!(coef_name %in% resultsNames(dds1))) next
  
  res_raw <- results(dds1, contrast = c("group", g_oht, "DMSO"))
  res_shr <- tryCatch(
    lfcShrink(dds1, coef = coef_name, res = res_raw, type = shrink_type),
    error = function(e) lfcShrink(dds1, coef = coef_name, res = res_raw, type = "normal")
  )
  
  res_df <- as.data.frame(res_shr)
  res_df[[id_col]] <- strip_version(rownames(res_df))
  res_df <- merge(res_df, annotation2, by = id_col, all.x = TRUE, sort = FALSE)
  
  nc1 <- counts(dds1, normalized = TRUE)
  idx_t <- which(metadata$group == g_oht)
  idx_c <- which(metadata$group == "DMSO")
  res_df$avg_treated <- rowMeans(nc1[, idx_t, drop = FALSE], na.rm = TRUE)
  res_df$avg_control <- rowMeans(nc1[, idx_c, drop = FALSE], na.rm = TRUE)
  
  res_df$sig_shr <- ifelse(!is.na(res_df$padj) & res_df$padj < padj_cutoff & res_df$log2FoldChange >  lfc_cutoff, "up",
                           ifelse(!is.na(res_df$padj) & res_df$padj < padj_cutoff & res_df$log2FoldChange < -lfc_cutoff, "down", "ns"))
  
  write.table(res_df, file.path(shr_tbl, paste0("tT_", nm, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(subset(res_df, sig_shr == "up"),   file.path(shr_tbl, paste0("up_", nm, ".tsv")),   sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(subset(res_df, sig_shr == "down"), file.path(shr_tbl, paste0("down_", nm, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
  
  counts_summary_shr <- rbind(counts_summary_shr, data.frame(
    contrast = nm, n_up_shr = sum(res_df$sig_shr == "up"), n_down_shr = sum(res_df$sig_shr == "down")
  ))
}

write.table(counts_summary_shr, file.path(shr_tbl, "DE_counts_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

#========================
# SHRUNK plots (MA/Volcano)
#========================
files_shr <- list.files(shr_tbl, pattern = "^tT_.*\\.tsv$", full.names = TRUE)

lim_shr <- calc_volcano_limits(files_shr, q = 0.995)
x_lim_shr <- lim_shr$x_lim
y_lim_shr <- lim_shr$y_lim

for (fp in files_shr) {
  
  nm_raw <- tools::file_path_sans_ext(basename(fp))
  nm <- sub("^tT_", "", nm_raw)
  
  res_df <- tryCatch(read.delim(fp, check.names = FALSE), error = function(e) NULL)
  if (is.null(res_df) || nrow(res_df) == 0) next
  if (!"sig_shr" %in% colnames(res_df)) next
  
  res_df$cat <- factor(res_df$sig_shr, levels = c("ns", "down", "up"))
  plot_df <- res_df %>% arrange(cat)
  
  p_ma <- ggplot(plot_df, aes(x = log10(baseMean + 1), y = log2FoldChange, color = cat)) +
    geom_point(size = 0.6, alpha = 0.85) +
    scale_color_manual(values = col_scale,
                       breaks = c("up","down","ns"),
                       labels = c("Up","Down","NS"),
                       name = "DE category") +
    geom_hline(yintercept = c(-lfc_cutoff, lfc_cutoff), linetype = "dashed") +
    labs(title = paste("MA (shrunk):", nm),
         x = "log10(baseMean + 1)", y = "log2(Fold Change)") +
    theme_bw(10)
  
  p_vol <- ggplot(plot_df, aes(x = log2FoldChange, y = -log10(pmax(padj, .Machine$double.xmin)), color = cat)) +
    geom_point(size = 0.6, alpha = 0.8) +
    scale_color_manual(values = col_scale,
                       breaks = c("up","down","ns"),
                       labels = c("Up","Down","NS"),
                       name = "DE category") +
    geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), linetype = "dashed") +
    geom_hline(yintercept = -log10(padj_cutoff), linetype = "dashed") +
    labs(title = paste("Volcano (shrunk):", nm),
         x = "log2(Fold Change)", y = "-log10(padj)") +
    theme_bw(10)
  
  if (!is.null(x_lim_shr) && !is.null(y_lim_shr)) p_vol <- p_vol + coord_cartesian(xlim = x_lim_shr, ylim = y_lim_shr)
  
  save_plot(file.path(shr_fig, paste0("MA_", nm)), p_ma, width = 6, height = 5)
  save_plot(file.path(shr_fig, paste0("Volcano_", nm)), p_vol, width = 6.5, height = 5)
}

message("Done: results/1_DEAnalysis (RAW + SHRUNK, incl. drug_OHT__drug and drug_OHT_vs_DMSO).")
