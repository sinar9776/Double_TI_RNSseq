## ============================================================
## 5_ORA_Categories.R — ORA per category
## ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(clusterProfiler)
  library(msigdbr)
  library(org.Mm.eg.db)
})

set.seed(1)
message(">>> 5_ORA_Categories running")

## -------------------- Parameters --------------------
analysis_level <- "gene"
padj_cutoff    <- 0.05
species_name   <- "Mus musculus"

## Colors (KEEP CONSISTENT with 3_enrichment)
col_up    <- "#D62728"
col_down  <- "#1F77B4"
col_other <- "#7F7F7F"

strip_version <- function(x) sub("\\..*$", "", x)

## -------------------- Directories --------------------
setwd("/hpcnfs/data/BA/glucocorticoid_keep/sinaR/m.rezaei")
project_root <- getwd()
cat_tbl <- file.path(project_root, "results/4_categories/tables")
tbl_dir <- file.path(project_root, "results/5_ORA_categories/tables")
fig_dir <- file.path(project_root, "results/5_ORA_categories/figures")
dir.create(tbl_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

## -------------------- Annotation (optional symbols) --------------------
annotation_fp <- file.path(project_root, "annotation.txt")
id2sym <- NULL
if (file.exists(annotation_fp)) {
  annotation <- read.delim(annotation_fp, check.names = FALSE)
  sym_col <- if ("symbol" %in% names(annotation)) "symbol" else if ("gene_name" %in% names(annotation)) "gene_name" else NA
  if (!is.na(sym_col) && "gene_id" %in% names(annotation)) {
    ann2 <- annotation %>%
      dplyr::transmute(gene_id = strip_version(gene_id), SYMBOL = .data[[sym_col]]) %>%
      dplyr::distinct(gene_id, .keep_all = TRUE)
    id2sym <- setNames(ann2$SYMBOL, ann2$gene_id)
  }
}

## -------------------- Universe (REQUIRED) --------------------
universe_fp <- file.path(project_root, "results", "1_DEAnalysis", "tables",
                         paste0("expressed_universe_", analysis_level, ".tsv"))
if (!file.exists(universe_fp)) stop("Missing universe: ", universe_fp)

u <- read.delim(universe_fp, check.names = FALSE)
if (ncol(u) < 1) stop("Universe file has no columns: ", universe_fp)

universe_ens <- unique(na.omit(strip_version(u[[1]])))
if (length(universe_ens) < 50) stop("Universe too small (<50): ", universe_fp)

## Universe (ENTREZ) for KEGG ORA
universe_entrez <- mapIds(
  org.Mm.eg.db,
  keys      = universe_ens,
  column    = "ENTREZID",
  keytype   = "ENSEMBL",
  multiVals = "first"
) %>% unname() %>% na.omit() %>% unique()

## -------------------- MSigDB gene sets --------------------
message("Loading MSigDB...")
msigdb <- list(
  HALLMARK = msigdbr::msigdbr(species = species_name, category = "H"),
  GOBP     = msigdbr::msigdbr(species = species_name, category = "C5", subcategory = "BP"),
  C2CP     = msigdbr::msigdbr(species = species_name, category = "C2", subcategory = "CP")
)

term2gene_list <- lapply(msigdb, function(df) {
  df %>%
    dplyr::transmute(gs_name, ensembl_gene = strip_version(ensembl_gene)) %>%
    dplyr::distinct()
})

term_sizes_u <- lapply(term2gene_list, function(t2g) {
  t2g %>%
    dplyr::filter(ensembl_gene %in% universe_ens) %>%
    dplyr::count(gs_name, name = "setSize")
})

## -------------------- Theme --------------------
theme_enrich <- function() {
  theme_bw(9) +
    theme(
      panel.grid  = element_blank(),
      plot.title  = element_text(face = "bold", hjust = 0.5),
      axis.text.y = element_text(size = 8)
    )
}

## -------------------- Plot --------------------
plot_ora_bar <- function(df_sig, out_base, fill_col) {
  if (nrow(df_sig) == 0) return(invisible(NULL))
  
  dfp <- df_sig %>%
    dplyr::arrange(p.adjust) %>%
    dplyr::slice_head(n = 20) %>%
    dplyr::mutate(Description = forcats::fct_reorder(Description, GeneRatio))
  
  p <- ggplot(dfp, aes(x = GeneRatio, y = Description)) +
    geom_col(fill = fill_col) +
    geom_text(aes(x = GeneRatio / 2, label = RatioLabel),
              size = 3.2, color = "black", fontface = "bold") +
    labs(x = "GeneRatio", y = NULL) +
    theme_enrich()
  
  ggsave(paste0(out_base, ".pdf"),  p, width = 7, height = 5)
  ggsave(paste0(out_base, ".tiff"), p, width = 7, height = 5, dpi = 300,
         device = "tiff", compression = "lzw")
}

## -------------------- Fill color from tag (MATCH 3_enrichment) --------------------
pick_fill <- function(tag) {
  if (grepl("up",   tag, ignore.case = TRUE)) return(col_up)
  if (grepl("down", tag, ignore.case = TRUE)) return(col_down)
  col_other
}

## ============================================================
## ORA (MSigDB + KEGG) for each category file
## (includes shifted_baseline_* automatically)
## ============================================================
cat_files <- list.files(cat_tbl, pattern = "\\.tsv$", full.names = TRUE)
if (length(cat_files) == 0) stop("No category files found in ", cat_tbl)

for (fp in cat_files) {
  tag <- tools::file_path_sans_ext(basename(fp))
  message("ORA: ", tag)
  
  df <- read.delim(fp, check.names = FALSE)
  
  id_col <- if ("gene_id" %in% names(df)) "gene_id" else
    if ("feature_id" %in% names(df)) "feature_id" else
      if ("isoform_id" %in% names(df)) "isoform_id" else NULL
  if (is.null(id_col)) { message("  [skip] No ID column: ", tag); next }
  
  genes_ens <- unique(na.omit(strip_version(df[[id_col]])))
  if (length(genes_ens) < 1) { message("  [skip] 0 genes: ", tag); next }
  
  fill_col <- pick_fill(tag)
  inputN   <- length(genes_ens)
  
  ## ---- MSigDB ORA ----
  for (set_name in names(term2gene_list)) {
    t2g <- term2gene_list[[set_name]]
    
    res <- tryCatch(
      clusterProfiler::enricher(
        gene          = genes_ens,
        universe      = universe_ens,
        TERM2GENE     = t2g,
        pAdjustMethod = "BH",
        qvalueCutoff  = padj_cutoff
      ),
      error = function(e) NULL
    )
    if (is.null(res)) next
    
    df_out <- as.data.frame(res)
    if (nrow(df_out) == 0) next
    
    df_out <- df_out %>%
      dplyr::left_join(term_sizes_u[[set_name]], by = c("ID" = "gs_name")) %>%
      dplyr::mutate(
        setSize    = tidyr::replace_na(setSize, 0L),
        GeneRatio  = Count / pmax(inputN, 1L),
        RatioLabel = paste0(Count, "/", setSize)
      )
    
    if (!is.null(id2sym)) {
      df_out$SYMBOLS <- vapply(
        strsplit(df_out$geneID, "/"),
        function(ids) paste(na.omit(id2sym[strip_version(ids)]), collapse = ";"),
        character(1)
      )
    }
    
    out_tsv <- file.path(tbl_dir, paste0("ORA_", set_name, "_", tag, ".tsv"))
    write.table(df_out, out_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
    
    df_sig <- df_out %>% dplyr::filter(p.adjust < padj_cutoff)
    if (nrow(df_sig) == 0) next
    
    out_base <- file.path(fig_dir, paste0("ORA_", set_name, "_", tag))
    plot_ora_bar(df_sig, out_base, fill_col = fill_col)
  }
  
  ## ---- KEGG ORA ----
  message("  KEGG ORA: ", tag)
  
  genes_entrez <- mapIds(
    org.Mm.eg.db,
    keys      = genes_ens,
    column    = "ENTREZID",
    keytype   = "ENSEMBL",
    multiVals = "first"
  ) %>% unname() %>% na.omit() %>% unique()
  
  if (length(genes_entrez) < 1 || length(universe_entrez) < 50) next
  
  kegg <- tryCatch(
    clusterProfiler::enrichKEGG(
      gene          = genes_entrez,
      universe      = universe_entrez,
      organism      = "mmu",
      pAdjustMethod = "BH",
      qvalueCutoff  = padj_cutoff
    ),
    error = function(e) NULL
  )
  if (is.null(kegg)) next
  
  df_out <- as.data.frame(kegg)
  if (nrow(df_out) == 0) next
  
  df_out <- df_out %>%
    dplyr::mutate(
      GeneRatio  = Count / pmax(length(genes_entrez), 1L),
      setSize    = as.integer(sub("/.*$", "", BgRatio)),
      RatioLabel = paste0(Count, "/", setSize)
    )
  
  out_tsv <- file.path(tbl_dir, paste0("ORA_KEGG_", tag, ".tsv"))
  write.table(df_out, out_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
  
  df_sig <- df_out %>% dplyr::filter(p.adjust < padj_cutoff)
  if (nrow(df_sig) == 0) next
  
  out_base <- file.path(fig_dir, paste0("ORA_KEGG_", tag))
  plot_ora_bar(df_sig, out_base, fill_col = fill_col)
}



## ============================================================
## HALLMARK HEATMAP — one heatmap per category across all drugs
## ============================================================
message(">>> Building Hallmark heatmaps...")

library(pheatmap)

## ---- Load Hallmark gene sets (same method as barplots) ----
msigdb_h <- msigdbr::msigdbr(species = species_name, collection = "H")

term2gene_h <- msigdb_h %>%
  dplyr::transmute(
    term = gs_name,
    gene = strip_version(ensembl_gene)
  ) %>%
  dplyr::distinct()

hallmark_terms  <- unique(term2gene_h$term)

## ---- ORA helper using clusterProfiler (IDENTICAL to barplots) ----
run_ora_hallmark <- function(genes_ens) {
  genes_ens <- unique(genes_ens)
  
  if (length(genes_ens) < 1) {
    return(tibble::tibble(term = hallmark_terms, pvalue = 1, padj = 1))
  }
  
  res <- tryCatch(
    clusterProfiler::enricher(
      gene          = genes_ens,
      universe      = universe_ens,
      TERM2GENE     = term2gene_h,
      pAdjustMethod = "BH",
      qvalueCutoff  = 1,
      pvalueCutoff  = 1
    ),
    error = function(e) NULL
  )
  
  if (is.null(res)) {
    return(tibble::tibble(term = hallmark_terms, pvalue = 1, padj = 1))
  }
  
  df_res <- as.data.frame(res) %>%
    dplyr::transmute(term = ID, pvalue = pvalue, padj = p.adjust)
  
  ## fill missing terms with p=1
  missing_terms <- setdiff(hallmark_terms, df_res$term)
  if (length(missing_terms) > 0) {
    df_res <- dplyr::bind_rows(
      df_res,
      tibble::tibble(term = missing_terms, pvalue = 1, padj = 1)
    )
  }
  df_res
}

## ---- Settings ----
cat_names_heatmap <- c(
  "enhanced_up", "enhanced_down",
  "suppressed_up", "suppressed_down",
  "independent_up", "independent_down",
  "shifted_baseline_enhanced_up", "shifted_baseline_enhanced_down",
  "shifted_baseline_suppressed_up", "shifted_baseline_suppressed_down",
  "shifted_baseline_independent_up", "shifted_baseline_independent_down",
  "additive_up", "additive_down"
)

drugs_ordered <- c(
  "Budesonide_2", "Budesonide_4", "Budesonide_8", "Budesonide_16",
  "Prednisolone_2", "Prednisolone_4", "Prednisolone_8", "Prednisolone_16"
)

cap_logp       <- 10
pval_thresh    <- 0.05
top_n_fallback <- 20
pal_heatmap    <- colorRampPalette(c("white", "red", "darkred"))(100)
heatmap_breaks <- seq(0, cap_logp, length.out = 101)

## ---- Loop over categories ----
for (cat in cat_names_heatmap) {
  message("  [Heatmap] ", cat)
  
  enrich_list <- vector("list", length(drugs_ordered))
  names(enrich_list) <- drugs_ordered
  
  for (drug in drugs_ordered) {
    fp <- file.path(cat_tbl, paste0(cat, "_", drug, ".tsv"))
    
    if (!file.exists(fp)) {
      enrich_list[[drug]] <- tibble::tibble(term = hallmark_terms, pvalue = 1, padj = 1)
      next
    }
    
    df_cat <- read.delim(fp, check.names = FALSE)
    id_col_cat <- if ("gene_id" %in% names(df_cat)) "gene_id" else
      if ("feature_id" %in% names(df_cat)) "feature_id" else NULL
    
    if (is.null(id_col_cat)) {
      enrich_list[[drug]] <- tibble::tibble(term = hallmark_terms, pvalue = 1, padj = 1)
      next
    }
    
    genes <- unique(na.omit(strip_version(df_cat[[id_col_cat]])))
    
    if (length(genes) < 1) {
      enrich_list[[drug]] <- tibble::tibble(term = hallmark_terms, pvalue = 1, padj = 1)
      next
    }
    
    enrich_list[[drug]] <- run_ora_hallmark(genes)
  }
  
  ## ---- Build matrix ----
  sig_terms <- unique(unlist(lapply(enrich_list, function(df) {
    if (is.null(df)) return(NULL)
    df$term[!is.na(df$padj) & df$padj < pval_thresh]
  })))
  
  if (length(sig_terms) == 0) {
    all_pvals <- dplyr::bind_rows(lapply(names(enrich_list), function(d) {
      enrich_list[[d]] %>% dplyr::mutate(drug = d)
    }))
    sig_terms <- all_pvals %>%
      dplyr::group_by(term) %>%
      dplyr::summarise(min_p = min(pvalue, na.rm = TRUE), .groups = "drop") %>%
      dplyr::arrange(min_p) %>%
      dplyr::slice_head(n = top_n_fallback) %>%
      dplyr::pull(term)
    
    if (length(sig_terms) == 0) {
      message("    No signal for ", cat, " -> skipping heatmap.")
      next
    }
  }
  
  mat <- matrix(
    0,
    nrow = length(sig_terms),
    ncol = length(drugs_ordered),
    dimnames = list(sig_terms, drugs_ordered)
  )
  
  star_mat <- matrix(
    "",
    nrow = length(sig_terms),
    ncol = length(drugs_ordered),
    dimnames = list(sig_terms, drugs_ordered)
  )
  
  for (drug in drugs_ordered) {
    df_j <- enrich_list[[drug]]
    if (is.null(df_j)) next
    for (tr in sig_terms) {
      row_j <- df_j[df_j$term == tr, ]
      if (nrow(row_j) == 0) next
      padj_val <- row_j$padj[1]
      if (is.na(padj_val) || padj_val <= 0) next
      mat[tr, drug]  <- min(-log10(padj_val), cap_logp)
      if (padj_val < pval_thresh) star_mat[tr, drug] <- "*"
    }
  }
  
  if (all(mat == 0)) {
    message("    All zeros for ", cat, " -> skipping heatmap.")
    next
  }
  
  ## ---- Save heatmap ----
  out_base_hm <- file.path(fig_dir, paste0("HALLMARK_heatmap_", cat))
  plot_height <- max(4, nrow(mat) * 0.17 + 2)
  
  ph <- pheatmap::pheatmap(
    mat,
    color           = pal_heatmap,
    breaks          = heatmap_breaks,
    main            = paste("Hallmark ORA:", cat),
    legend_breaks   = c(0, 2.5, 5, 7.5, 10),
    legend_labels   = c("0", "2.5", "5", "7.5", "-log10(padj)\\n10"),
    cluster_rows    = nrow(mat) >= 2,
    cluster_cols    = FALSE,
    fontsize_row    = 7,
    fontsize_col    = 9,
    na_col          = "white",
    display_numbers = star_mat,
    number_color    = "black",
    fontsize_number = 10,
    silent          = TRUE
  )
  
  ggsave(paste0(out_base_hm, ".pdf"),
         ph, width = 9, height = plot_height)
  ggsave(paste0(out_base_hm, ".tiff"),
         ph, width = 9, height = plot_height,
         dpi = 300, device = "tiff", compression = "lzw")
  
  write.csv(
    as.data.frame(mat),
    file.path(tbl_dir, paste0("HALLMARK_matrix_", cat, ".csv")),
    row.names = TRUE
  )
  
  message("    Saved: HALLMARK_heatmap_", cat)
}

message(">>> Hallmark heatmaps done.")

## -------------------- Reproducibility --------------------
writeLines(capture.output(sessionInfo()), file.path(tbl_dir, "sessionInfo.txt"))

message(">>> done")
message("Tables: ", tbl_dir)
message("Figures: ", fig_dir)
