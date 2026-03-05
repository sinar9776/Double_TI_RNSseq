## ============================================================
## 3_enrichment.R — GSEA + ORA (gene_id-based)
## Color rules + ratio labels
## ============================================================


suppressPackageStartupMessages({
  library(tidyverse)
  library(fgsea)
  library(msigdbr)
  library(clusterProfiler)
  library(enrichplot)
  library(org.Mm.eg.db)
})
BiocManager::install("clusterProfiler")
set.seed(1)
message(">>> 3_enrichment running")

## -------------------- Parameters --------------------
n_cores     <- 4
padj_cutoff <- 0.05

analysis_level <- "gene"          # used only for universe filename; enrichment uses gene_id
species_name   <- "Mus musculus"  # manual

col_up   <- "#D62728"
col_down <- "#1F77B4"

strip_version <- function(x) sub("\\..*$", "", x)

## -------------------- Parallel backend --------------------
BiocParallel::register(BiocParallel::SnowParam(workers = n_cores, type = "SOCK"))
options(mc.cores = n_cores)

## -------------------- Directories --------------------
setwd("/hpcnfs/data/BA/glucocorticoid_keep/sinaR/m.rezaei")
project_root <- getwd()
data_dir     <- file.path(project_root, "data")

de_tbl   <- file.path(project_root, "results", "1_DEAnalysis", "tables", "shrunk")
int_tbl  <- file.path(project_root, "results", "2_interaction", "tables")

res_dir  <- file.path(project_root, "results", "3_enrichment")
tbl_dir  <- file.path(res_dir, "tables")
fig_dir  <- file.path(res_dir, "figures")
dir.create(tbl_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

## -------------------- Parameters --------------------
padj_cutoff <- 0.05

## -------------------- Annotation --------------------
annotation <- read.delim(file.path(project_root, "annotation.txt"), check.names = FALSE)
annotation <- annotation %>%
  dplyr::select(gene_id, symbol) %>%
  distinct(gene_id, .keep_all = TRUE)
id2sym <- setNames(annotation$symbol, annotation$gene_id)

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

## -------------------- MSigDB (for GSEA + ORA) --------------------
message("Loading MSigDB...")
msigdb <- list(
  HALLMARK = msigdbr::msigdbr(species = "Mus musculus", category = "H"),
  GOBP     = msigdbr::msigdbr(species = "Mus musculus", category = "C5", subcategory = "BP"),
  C2CP     = msigdbr::msigdbr(species = "Mus musculus", category = "C2", subcategory = "CP")
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

gmt_list <- lapply(term2gene_list, function(t2g) split(t2g$ensembl_gene, t2g$gs_name))

## -------------------- KEGG gene sets for GSEA (ENSEMBL IDs) --------------------
build_kegg_gmt <- function(universe_ens, organism = "mmu") {
  k <- clusterProfiler::download_KEGG(species = organism)
  
  # Join pathway IDs to readable names
  id2name <- k$KEGGPATHID2NAME %>%
    tibble::as_tibble() %>%
    dplyr::transmute(from, pathway_name = to)
  
  t2g_entrez <- k$KEGGPATHID2EXTID %>%
    tibble::as_tibble() %>%
    dplyr::left_join(id2name, by = "from") %>%
    dplyr::transmute(gs_name = ifelse(!is.na(pathway_name), pathway_name, from),
                     ENTREZID = as.character(to)) %>%
    dplyr::filter(!is.na(ENTREZID))
  
  # Keep names on the vector for correct lookup
  entrez_to_ens <- mapIds(
    org.Mm.eg.db,
    keys      = unique(t2g_entrez$ENTREZID),
    keytype   = "ENTREZID",
    column    = "ENSEMBL",
    multiVals = "first"
  )
  
  t2g_ens <- t2g_entrez %>%
    dplyr::mutate(ensembl_gene = strip_version(entrez_to_ens[ENTREZID])) %>%
    dplyr::filter(!is.na(ensembl_gene), ensembl_gene %in% universe_ens) %>%
    dplyr::distinct(gs_name, ensembl_gene)
  
  split(t2g_ens$ensembl_gene, t2g_ens$gs_name)
}

message("Loading KEGG gene sets for GSEA...")
kegg_gmt <- build_kegg_gmt(universe_ens, organism = "mmu")
message("KEGG pathways loaded: ", length(kegg_gmt))

## -------------------- Theme --------------------
theme_enrich <- function() {
  theme_bw(9) +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.text.y = element_text(size = 8)
    )
}

## -------------------- Plots --------------------
plot_gsea_bar <- function(gres_full, out_base, title_text) {
  dfp <- gres_full %>%
    dplyr::filter(!is.na(padj), padj < padj_cutoff)
  if (nrow(dfp) == 0) return(invisible(NULL))
  
  dfp <- dfp %>%
    dplyr::mutate(
      leadingEdgeCount = vapply(leadingEdge, length, integer(1)),
      ratio_label = paste0(leadingEdgeCount, "/", size)
    ) %>%
    dplyr::arrange(dplyr::desc(abs(NES))) %>%
    dplyr::slice_head(n = 20) %>%
    dplyr::mutate(
      pathway = forcats::fct_reorder(pathway, NES),
      dir = ifelse(NES > 0, "Up", "Down")
    )
  
  p <- ggplot(dfp, aes(x = NES, y = pathway, fill = dir)) +
    geom_col() +
    geom_text(aes(x = NES / 2, label = ratio_label),
              color = "black", size = 3.2, fontface = "bold") +
    scale_fill_manual(values = c(Up = col_up, Down = col_down),
                      name = "Direction") +
    labs(title = title_text, x = "NES", y = NULL) +
    theme_enrich()
  
  ggsave(paste0(out_base, ".pdf"),  p, width = 7, height = 5)
  ggsave(paste0(out_base, ".tiff"), p, width = 7, height = 5, dpi = 300,
         device = "tiff", compression = "lzw")
}

plot_ora_bar <- function(df_sig, out_base, fill_col) {
  if (nrow(df_sig) == 0) return(invisible(NULL))
  
  dfp <- df_sig %>%
    dplyr::arrange(p.adjust) %>%
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

## -------------------- Helpers --------------------
write_gsea_table <- function(gres, out_tsv) {
  if (is.null(gres) || nrow(gres) == 0) return(invisible(NULL))
  
  out_tbl <- gres %>%
    dplyr::arrange(padj, dplyr::desc(NES)) %>%
    dplyr::transmute(
      pathway = .data$pathway,
      size, ES, NES, pval, padj,
      leadingEdge
    )
  
  if (!is.null(id2sym)) {
    out_tbl <- out_tbl %>%
      dplyr::mutate(
        leadingEdge_symbols = vapply(
          leadingEdge,
          function(ids) paste(na.omit(id2sym[strip_version(ids)]), collapse = ";"),
          character(1)
        )
      )
  }
  
  out_tbl <- out_tbl %>%
    dplyr::mutate(leadingEdge = vapply(leadingEdge, paste, collapse = ";", character(1)))
  
  write.table(out_tbl, out_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
  out_tbl
}

## ============================================================
## 1) GSEA (MSigDB + KEGG)
## ============================================================
tt_files <- c(
  list.files(de_tbl,  pattern = "^tT_.*\\.tsv$", full.names = TRUE),
  list.files(int_tbl, pattern = "^tT_.*\\.tsv$", full.names = TRUE)
)
if (length(tt_files) == 0) stop("No tT_*.tsv files found.")

for (fp in tt_files) {
  tag_raw <- tools::file_path_sans_ext(basename(fp))
  tag <- sub("^tT_", "", tag_raw)
  message("GSEA: ", tag)
  
  df <- read.delim(fp, check.names = FALSE)
  if (!("gene_id" %in% names(df))) stop("Missing gene_id in: ", basename(fp))
  
  rank_col <- if ("stat" %in% names(df)) "stat" else if ("log2FoldChange" %in% names(df)) "log2FoldChange" else NA
  if (is.na(rank_col)) stop("Need stat or log2FoldChange in: ", basename(fp))
  
  df2 <- df %>%
    dplyr::transmute(gene_id = strip_version(gene_id), rank = .data[[rank_col]]) %>%
    dplyr::filter(!is.na(gene_id), is.finite(rank), gene_id %in% universe_ens) %>%
    dplyr::group_by(gene_id) %>%
    dplyr::summarise(rank = rank[which.max(abs(rank))], .groups = "drop")
  
  ranks <- setNames(df2$rank, df2$gene_id)
  if (length(ranks) < 50) next
  
  ## ---- MSigDB GSEA ----
  for (set_tag in names(gmt_list)) {
    gres <- tryCatch(
      fgsea::fgseaMultilevel(
        pathways = gmt_list[[set_tag]],
        stats    = ranks,
        BPPARAM  = BiocParallel::bpparam()
      ),
      error = function(e) NULL
    )
    if (is.null(gres) || nrow(gres) == 0) next
    
    out_tsv <- file.path(tbl_dir, paste0("GSEA_", set_tag, "_", tag, ".tsv"))
    write_gsea_table(gres, out_tsv)
    
    out_base <- file.path(fig_dir, paste0("GSEA_", set_tag, "_", tag))
    plot_gsea_bar(
      gres_full  = gres %>% dplyr::transmute(pathway, NES, padj, size, leadingEdge),
      out_base   = out_base,
      title_text = paste0(set_tag, " - ", tag)
    )
  }
  
  ## ---- KEGG GSEA ----
  gres_kegg <- tryCatch(
    fgsea::fgseaMultilevel(
      pathways = kegg_gmt,
      stats    = ranks,
      BPPARAM  = BiocParallel::bpparam()
    ),
    error = function(e) NULL
  )
  
  if (!is.null(gres_kegg) && nrow(gres_kegg) > 0) {
    out_tsv <- file.path(tbl_dir, paste0("GSEA_KEGG_", tag, ".tsv"))
    write_gsea_table(gres_kegg, out_tsv)
    
    out_base <- file.path(fig_dir, paste0("GSEA_KEGG_", tag))
    plot_gsea_bar(
      gres_full  = gres_kegg %>% dplyr::transmute(pathway, NES, padj, size, leadingEdge),
      out_base   = out_base,
      title_text = paste0("KEGG - ", tag)
    )
  }
}

## ============================================================
## 2) ORA (MSigDB) — x=GeneRatio, label=Count/setSize
## ============================================================
up_files <- c(
  list.files(de_tbl,  pattern = "^up_.*\\.tsv$", full.names = TRUE),
  list.files(int_tbl, pattern = "^up_.*\\.tsv$", full.names = TRUE)
)
down_files <- c(
  list.files(de_tbl,  pattern = "^down_.*\\.tsv$", full.names = TRUE),
  list.files(int_tbl, pattern = "^down_.*\\.tsv$", full.names = TRUE)
)
ora_files <- c(up_files, down_files)
if (length(ora_files) == 0) stop("No up_/down_ tables found.")

for (fp in ora_files) {
  tag <- tools::file_path_sans_ext(basename(fp))
  message("ORA: ", tag)
  
  df <- read.delim(fp, check.names = FALSE)
  if (!("gene_id" %in% names(df))) stop("Missing gene_id in: ", basename(fp))
  
  genes_ens <- unique(na.omit(strip_version(df$gene_id)))
  if (length(genes_ens) < 10) next
  
  fill_col <- if (grepl("^up_", tag)) col_up else col_down
  inputN <- length(genes_ens)
  
  for (set_name in names(term2gene_list)) {
    t2g <- term2gene_list[[set_name]]
    
    res <- tryCatch(
      clusterProfiler::enricher(
        gene      = genes_ens,
        universe  = universe_ens,
        TERM2GENE = t2g,
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
}

## ============================================================
## 3) KEGG ORA (enrichKEGG) — x=GeneRatio, label=Count/setSize
## ============================================================
for (fp in ora_files) {
  tag <- tools::file_path_sans_ext(basename(fp))
  message("KEGG ORA: ", tag)
  
  df <- read.delim(fp, check.names = FALSE)
  if (!("gene_id" %in% names(df))) stop("Missing gene_id in: ", basename(fp))
  
  genes_ens <- unique(na.omit(strip_version(df$gene_id)))
  if (length(genes_ens) < 10) next
  
  fill_col <- if (grepl("^up_", tag)) col_up else col_down
  
  genes_entrez <- mapIds(
    org.Mm.eg.db,
    keys      = genes_ens,
    column    = "ENTREZID",
    keytype   = "ENSEMBL",
    multiVals = "first"
  ) %>% unname() %>% na.omit() %>% unique()
  
  if (length(genes_entrez) < 10 || length(universe_entrez) < 50) next
  
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
  
  inputN <- length(genes_entrez)
  
  df_out <- df_out %>%
    dplyr::mutate(
      GeneRatio  = Count / pmax(inputN, 1L),
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

## -------------------- Reproducibility --------------------
writeLines(capture.output(sessionInfo()), file.path(res_dir, "sessionInfo.txt"))

message(">>> done")
message("Tables: ", tbl_dir)
message("Figures: ", fig_dir)