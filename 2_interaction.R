## ============================================================
## 2_interaction.R — Interaction analysis (drug-OHT)
## Reference: drug = DMSO | OHT = OFF
## Outputs:
##   tables/: tT_int_<drug>.tsv, up_int_<drug>.tsv, down_int_<drug>.tsv
##   figures/: MA_int_<drug>.png, Volcano_int_<drug>.png
## ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(DESeq2)
  library(ggplot2)
})

set.seed(1)
message(">>> 2_interaction: baseline-OHT interaction mode")

## -------------------- Directories --------------------
setwd("/hpcnfs/data/BA/glucocorticoid_keep/sinaR/m.rezaei")
base_dir    <- getwd()
data_dir    <- file.path(base_dir, "data")
results_dir <- file.path(base_dir, "results", "2_interaction")
tbl_dir     <- file.path(results_dir, "tables")
fig_dir     <- file.path(results_dir, "figures")
dir.create(tbl_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

## -------------------- Parameters --------------------
padj_cutoff <- 0.05
lfc_cutoff  <- 0

## ============================================================
## PART 1 — Interaction DESeq2 Calculation 
## ============================================================

# Choose analysis level
analysis_level <- "gene"   # or "isoform"

## ---- Load metadata ----
metadata <- read.delim(file.path(base_dir, "metadata.txt"), check.names = FALSE)
stopifnot(all(c("SampleID", "drug", "OHT") %in% colnames(metadata)))

metadata$OHT  <- factor(metadata$OHT, levels = c("OFF", "ON"))
metadata$drug <- factor(metadata$drug)
metadata$drug <- relevel(metadata$drug, ref = "DMSO")

## ---- Load counts ----
if (analysis_level == "gene") {
  count_file <- file.path(base_dir, "gene_counts.tsv")
  id_col <- "gene_id"
} else {
  count_file <- file.path(base_dir, "isoform_counts.tsv")
  id_col <- "isoform_id"
}

counts <- readr::read_tsv(count_file, show_col_types = FALSE) %>% as.data.frame()
rownames(counts) <- counts[[1]]
count_mat <- counts[, intersect(colnames(counts), metadata$SampleID), drop = FALSE]
metadata  <- metadata[match(colnames(count_mat), metadata$SampleID), ]
stopifnot(all(colnames(count_mat) == metadata$SampleID))

## ---- Load annotation ----
annotation <- read.delim(file.path(base_dir, "annotation.txt"), check.names = FALSE)
sym_col <- if ("symbol" %in% names(annotation)) "symbol" else if ("gene_name" %in% names(annotation)) "gene_name" else NA
ann_keep_cols <- unique(c("gene_id", id_col, sym_col))
ann_keep_cols <- ann_keep_cols[ann_keep_cols %in% names(annotation)]
annotation2 <- annotation[, ann_keep_cols, drop = FALSE]
if (!is.na(sym_col)) names(annotation2)[names(annotation2) == sym_col] <- "SYMBOL"
if (analysis_level == "gene") {
  annotation2 <- annotation2 %>%
    group_by(gene_id) %>%
    summarise(SYMBOL = dplyr::first(na.omit(SYMBOL)), .groups = "drop")
}
if (!"gene_id" %in% names(annotation2)) stop("Annotation must include gene_id column.")

## ---- Build DESeq2 dataset ----
dds <- DESeqDataSetFromMatrix(
  countData = round(as.matrix(count_mat)),
  colData   = metadata,
  design    = ~ drug + OHT + drug:OHT
)
dds <- DESeq(dds, parallel = TRUE)

## ---- Get interaction terms ----
res_names <- resultsNames(dds)
interaction_terms <- grep("drug.*OHTON", res_names, value = TRUE)
if (length(interaction_terms) == 0) stop("No interaction terms detected.")
message("Detected interaction terms: ", paste(interaction_terms, collapse = ", "))

## ---- Save tables ----
for (term in interaction_terms) {
  drug <- gsub("^drug", "", gsub("\\.OHTON$", "", term))
  message("\n>>> Interaction: ", drug)
  
  res <- results(dds, name = term)
  res_df <- as.data.frame(res)
  res_df[[id_col]] <- sub("\\..*$", "", rownames(res_df))
  
  # Merge annotation and ensure gene_id present
  res_df <- merge(res_df, annotation2, by = id_col, all.x = TRUE, sort = FALSE)
  if (!"gene_id" %in% names(res_df) && "gene_id" %in% names(annotation2)) {
    res_df <- merge(res_df, annotation2[, c(id_col, "gene_id")], by = id_col, all.x = TRUE)
  }
  
  res_df$sig <- ifelse(!is.na(res_df$padj) & res_df$padj < padj_cutoff & res_df$log2FoldChange > lfc_cutoff, "up",
                       ifelse(!is.na(res_df$padj) & res_df$padj < padj_cutoff & res_df$log2FoldChange < -lfc_cutoff, "down", "ns"))
  
  # Save tables
  write.table(res_df, file.path(tbl_dir, paste0("tT_int_", drug, ".tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(subset(res_df, sig == "up"), file.path(tbl_dir, paste0("up_int_", drug, ".tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(subset(res_df, sig == "down"), file.path(tbl_dir, paste0("down_int_", drug, ".tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)
}

message("\n>>> Part 1 complete — Interaction DESeq2 results saved.")
message("Tables → ", tbl_dir)

## ============================================================
## PART 2 — Plotting (MA & Volcano)
## ============================================================

col_scale <- c("up" = "#D62728", "down" = "#1F77B4", "ns" = "gray80")
files <- list.files(tbl_dir, pattern = "^tT_int_.*\\.tsv$", full.names = TRUE)

for (fp in files) {
  nm <- sub("^tT_int_", "", tools::file_path_sans_ext(basename(fp)))
  res_df <- read.delim(fp, check.names = FALSE)
  res_df$cat <- factor(res_df$sig, levels = c("up", "down", "ns"))
  
  ## --- MA Plot ---
  p_ma <- ggplot() +
    geom_point(data = subset(res_df, cat == "ns"),
               aes(x = log10(baseMean + 1), y = log2FoldChange),
               color = "gray80", size = 1.1, alpha = 0.6) +
    geom_point(data = subset(res_df, cat == "down"),
               aes(x = log10(baseMean + 1), y = log2FoldChange),
               color = col_scale["down"], size = 1.3, alpha = 0.9) +
    geom_point(data = subset(res_df, cat == "up"),
               aes(x = log10(baseMean + 1), y = log2FoldChange),
               color = col_scale["up"], size = 1.3, alpha = 0.9) +
    geom_hline(yintercept = c(-lfc_cutoff, lfc_cutoff), linetype = "dashed") +
    labs(title = paste("MA (interaction):", nm),
         x = "log10(baseMean + 1)",
         y = "log2(Fold Change)") +
    theme_bw(10) +
    theme(legend.position = "right") +
    guides(color = guide_legend(title = "DE category",
                                override.aes = list(size = 3))) +
    scale_color_manual(name = "DE category",
                       values = c("Upregulated" = "#D62728", "Downregulated" = "#1F77B4"))
  
  ## --- Volcano Plot ---
  p_vol <- ggplot() +
    geom_point(data = subset(res_df, cat == "ns"),
               aes(x = log2FoldChange, y = -log10(padj)),
               color = "gray80", size = 1.1, alpha = 0.6) +
    geom_point(data = subset(res_df, cat == "down"),
               aes(x = log2FoldChange, y = -log10(padj)),
               color = col_scale["down"], size = 1.3, alpha = 0.9) +
    geom_point(data = subset(res_df, cat == "up"),
               aes(x = log2FoldChange, y = -log10(padj)),
               color = col_scale["up"], size = 1.3, alpha = 0.9) +
    geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), linetype = "dashed") +
    geom_hline(yintercept = -log10(padj_cutoff), linetype = "dashed") +
    labs(title = paste("Volcano (interaction):", nm),
         x = "log2(Fold Change)",
         y = "-log10(adjusted p-value)") +
    theme_bw(10) +
    theme(legend.position = "right") +
    guides(color = guide_legend(title = "DE category",
                                override.aes = list(size = 3))) +
    scale_color_manual(name = "DE category",
                       values = c("Upregulated" = "#D62728", "Downregulated" = "#1F77B4"))
  
  ## --- Save plots ---
  ggsave(file.path(fig_dir, paste0("MA_", nm, ".png")), p_ma, width = 6, height = 5, dpi = 300)
  ggsave(file.path(fig_dir, paste0("Volcano_", nm, ".png")), p_vol, width = 6.5, height = 5, dpi = 300)
}

message("\n>>> Part 2 complete — MA & Volcano plots saved.")
message("Figures saved in: ", fig_dir)

