BDCI
Baseline-Dependent Category Identifier

BaseDepEx is a modular, end-to-end R pipeline for identifying and characterizing genes whose transcriptional response to a drug treatment depends on the activation state of a regulatory factor (e.g., MYC via OHT-inducible system). It combines DESeq2-based differential expression, interaction modeling, functional enrichment, and mechanistic gene categorization at both gene and isoform resolution.

Table of Contents
Overview
Biological Rationale
Pipeline Architecture
Gene Category Definitions
Directory Structure
Input Requirements
Installation & Dependencies
Usage
Output Summary
Key Parameters
Citation
License
Overview
BaseDepEx is designed for experimental designs that combine a drug treatment with a secondary regulatory perturbation (e.g., an inducible transcription factor). For each drug, the pipeline asks: does this gene respond to the drug differently when the regulatory factor is active versus inactive?

The pipeline produces:

Differential expression tables (raw + LFC-shrunk) for all relevant contrasts
Interaction term statistics (drug × baseline)
GSEA and ORA results across MSigDB and KEGG gene sets
A principled, multi-criterion classification of genes into mechanistic categories
Publication-quality scatter plots and enrichment bar charts (PDF + TIFF)
Biological Rationale
The pipeline is built around a 2 × N factorial design:

Axis 1 (drug): multiple drug treatments vs. DMSO control
Axis 2 (baseline): regulatory factor OFF (OHT–) vs. ON (OHT+)
Four core comparisons are derived per drug:

Contrast	Notation	Biological question
drug / DMSO	tT_<drug>	What does the drug do at baseline?
(drug+OHT) / (DMSO+OHT)	tT_<drug>_OHT	What does the drug do under activation?
(drug+OHT) / drug	tT_<drug>_OHT__<drug>	Does expression change when OHT is added on top of the drug, independent of DMSO? Used to detect a shifted baseline.
(drug+OHT) / DMSO	tT_<drug>_OHT_vs_DMSO	What is the combined effect of drug + activation?
A separate interaction model (drug × OHT) statistically tests whether the drug effect is significantly modified by the baseline state.

Pipeline Architecture
The workflow is fully modular. Each script reads from and writes to standardized result directories, allowing individual steps to be re-run or modified independently.

data/
├── gene_counts.tsv / isoform_counts.tsv
├── metadata.txt
└── annotation.txt
        │
        ▼
┌─────────────────────────────────┐
│  1_DEAnalysis.R                 │  DESeq2 DE (all contrasts, raw + shrunk)
│  → results/1_DEAnalysis/        │  MA plots, Volcano plots, PCA, Dispersion
└─────────────┬───────────────────┘
              │
              ▼
┌─────────────────────────────────┐
│  2_interaction.R                │  DESeq2 interaction model (drug × OHT)
│  → results/2_interaction/       │  Interaction MA + Volcano plots
└─────────────┬───────────────────┘
              │
        ┌─────┴──────┐
        ▼            ▼
┌───────────────┐  ┌──────────────────────────────────┐
│ 3_enrichment  │  │  4_categories.R                  │
│ GSEA + ORA    │  │  Gene categorization + UpSet plots│
│ on DE lists   │  │  → results/4_categories/          │
└───────────────┘  └──────────────┬───────────────────┘
                                  │
                        ┌─────────┴─────────┐
                        ▼                   ▼
              ┌──────────────────┐  ┌────────────────────┐
              │ 5_ORA_Categories │  │  6_scatterplots.R  │
              │ ORA per category │  │  LFC scatter plots │
              └──────────────────┘  └────────────────────┘
Step	Script	Description
1	1_DEAnalysis.R	DESeq2 DE analysis for all contrasts; raw and apeglm-shrunk LFC tables; MA, Volcano, PCA, Dispersion plots
2	2_interaction.R	Factorial DESeq2 model (~ drug + OHT + drug:OHT); tests drug × OHT interaction per drug
3	3_enrichment.R	GSEA (fgseaMultilevel) and ORA (enricher / enrichKEGG) on DE gene lists; MSigDB Hallmark, GO BP, C2:CP, KEGG
4	4_categories.R	Integrates DE + interaction results into mechanistic gene categories; UpSet plots per drug
5	5_ORA_Categories.R	ORA on each category gene list against the expressed universe
6	6_scatterplots.R	LFC scatter plots (drug/DMSO vs drug+OHT/DMSO+OHT; OHT/DMSO vs drug+OHT/DMSO) colored by category
Gene Category Definitions
Categories are assigned per drug by integrating four binary signals per gene:

drug/DMSO — is the gene significantly DE in the drug-treated condition at baseline (OHT off)?
(drug+OHT)/(DMSO+OHT) — is it significantly DE under activation (OHT on)?
Interaction term — does the drug effect significantly differ between OHT off and OHT on (from the factorial model)?
(drug+OHT)/drug — does expression change significantly when OHT is added directly on top of the drug treatment (i.e., is the baseline itself shifted by the drug)?
The fourth comparison is used to split each primary category into a true variant (expression was stable on drug alone before OHT was added, so the OHT-driven change is interpretable) and a shifted-baseline variant (expression already changed when OHT was applied to drug-treated cells, indicating the baseline was not stable and the category assignment may be confounded by a pre-existing shift).

Category	Criteria	Interpretation
enhanced_up	Up in drug±OHT; interaction up; (drug+OHT)/drug not significant	Drug-driven upregulation is amplified by baseline activation
enhanced_down	Down in drug±OHT; interaction down; (drug+OHT)/drug not significant	Drug-driven repression is amplified by baseline activation
suppressed_up	Up in drug/DMSO, interaction down; (drug+OHT)/drug not significant	Drug-driven upregulation is blunted under activation
suppressed_down	Down in drug/DMSO, interaction up; (drug+OHT)/drug not significant	Drug-driven repression is blunted under activation
switched_positive	Down in drug/DMSO, up in (drug+OHT)/(DMSO+OHT); interaction up	Drug direction reverses when the regulatory factor is active
switched_negative	Up in drug/DMSO, down in (drug+OHT)/(DMSO+OHT); interaction down	Drug direction reverses when the regulatory factor is active
independent_up	Up in both ±OHT conditions; no interaction; (drug+OHT)/drug not significant	Drug-driven upregulation is unaffected by baseline activation state
independent_down	Down in both ±OHT conditions; no interaction; (drug+OHT)/drug not significant	Drug-driven repression is unaffected by baseline activation state
shifted_baseline_*	Any of the above + (drug+OHT)/drug significant	The drug-treated baseline is itself shifted by OHT; category is assigned but flagged separately
additive_up	Subset of independent_up; also significantly up in OHT/DMSO; no interaction	Gene is independently upregulated by both the drug and baseline activation
additive_down	Subset of independent_down; also significantly down in OHT/DMSO; no interaction	Gene is independently downregulated by both the drug and baseline activation
Directory Structure
BaseDepEx/
├── RScripts/
│   ├── 1_DEAnalysis.R
│   ├── 2_interaction.R
│   ├── 3_enrichment.R
│   ├── 4_categories.R
│   ├── 5_ORA_categories.R
│   └── 6_scatterplot.R
├── data/
│   ├── gene_counts.tsv          # Gene-level count matrix
│   ├── isoform_counts.tsv       # Isoform-level count matrix (optional)
│   ├── metadata.txt             # Sample metadata
│   └── annotation.txt           # Gene/isoform annotation
├── results/
│   ├── 1_DEAnalysis/
│   │   ├── tables/
│   │   │   ├── raw/             # Raw DESeq2 result tables
│   │   │   └── shrunk/          # apeglm-shrunk LFC tables
│   │   └── figures/
│   │       ├── raw/             # MA, Volcano, PCA, Dispersion (raw)
│   │       └── shrunk/          # MA, Volcano (shrunk)
│   ├── 2_interaction/
│   │   ├── tables/              # tT_int_*, up_int_*, down_int_*
│   │   └── figures/             # MA + Volcano per interaction term
│   ├── 3_enrichment/
│   │   ├── tables/              # GSEA + ORA TSV results
│   │   └── figures/             # Enrichment bar plots (PDF + TIFF)
│   ├── 4_categories/
│   │   ├── tables/              # Category gene lists + UpSet data
│   │   └── figures/             # UpSet plots (PDF + TIFF)
│   ├── 5_ORA_categories/
│   │   ├── tables/              # ORA results per category
│   │   └── figures/             # Bar plots per category (PDF + TIFF)
│   └── 6_scatterplots/
│       └── figures/             # Scatter plots (PDF + TIFF)
├── requirements.txt
├── LICENSE
└── README.md
Input Requirements
All input files must be placed in the data/ directory before running the pipeline.

metadata.txt
Tab-delimited. Required columns:

Column	Description
SampleID	Must match column names in count matrices
group	DESeq2 group label (e.g., DMSO, 11j, 11j_OHT, DMSO_OHT)
drug	Drug identity (e.g., DMSO, 11j, KVS); used by 2_interaction.R
OHT	Baseline state: OFF or ON; used by 2_interaction.R
The group column must follow the convention <drug> and <drug>_OHT for the pipeline to correctly pair contrasts. The control group must be labeled DMSO.

gene_counts.tsv / isoform_counts.tsv
Tab-delimited count matrices. First column is the feature ID (gene_id or isoform_id); remaining columns are samples, with names matching SampleID in metadata. Raw integer counts (not normalized, not log-transformed).

annotation.txt
Tab-delimited. Required columns:

Column	Description
gene_id	Ensembl gene ID (with or without version suffix)
symbol or gene_name	HGNC/MGI gene symbol
isoform_id	Required only for isoform-level analysis
expressed_universe_<level>.tsv
Generated automatically by 1_DEAnalysis.R as a byproduct of low-count filtering. Must be present in results/1_DEAnalysis/tables/ before running scripts 3 or 5. It contains one column of Ensembl gene IDs representing all expressed genes (the ORA/GSEA background universe).

Installation & Dependencies
BaseDepEx requires R ≥ 4.2. All packages are available through CRAN or Bioconductor.

Bioconductor packages
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c(
  "DESeq2",
  "apeglm",
  "BiocParallel",
  "fgsea",
  "clusterProfiler",
  "msigdbr",
  "org.Mm.eg.db"    # Replace with org.Hs.eg.db for human data
))
CRAN packages
install.packages(c(
  "tidyverse",
  "ggplot2",
  "patchwork",
  "ComplexUpset",
  "forcats"
))
Species note: The pipeline is configured for Mus musculus (org.Mm.eg.db, organism = "mmu"). For human data, replace with org.Hs.eg.db and organism = "hsa" in scripts 3 and 5, and update species_name to "Homo sapiens".

Usage
Scripts are run sequentially from the repository root. Set analysis_level <- "gene" or "isoform" at the top of each script as needed.

# Run from the BaseDepEx root directory
Rscript RScripts/1_DEAnalysis.R
Rscript RScripts/2_interaction.R
Rscript RScripts/3_enrichment.R
Rscript RScripts/4_categories.R
Rscript RScripts/5_ORA_categories.R
Rscript RScripts/6_scatterplot.R
Alternatively, source each script from an interactive R session with getwd() returning the repository root:

setwd("/path/to/BaseDepEx")
source("RScripts/1_DEAnalysis.R")
# ... continue sequentially
Scripts 3, 5, and 6 depend on outputs from scripts 1 and 4. Scripts 4 and 5 require the expressed universe file produced by script 1. The pipeline does not need to be re-run from scratch if only downstream parameters are changed.

Output Summary
All figures are saved in both PDF (vector, publication-ready) and TIFF (600 dpi, LZW-compressed, print-ready) formats. Interaction plots (step 2) are saved as PNG.

Step	Key outputs
1	tT_<contrast>.tsv, up_<contrast>.tsv, down_<contrast>.tsv (raw + shrunk); MA, Volcano, PCA plots
2	tT_int_<drug>.tsv, up/down_int_<drug>.tsv; interaction MA + Volcano
3	GSEA_<set>_<contrast>.tsv, ORA_<set>_<contrast>.tsv; enrichment bar plots
4	<category>_<drug>.tsv (one file per category per drug); upset_data_<drug>.tsv; UpSet figures
5	ORA_<set>_<category>_<drug>.tsv; category enrichment bar plots
6	Per-drug scatter plots + all-drug panels; CORR_OHT scatter plots; separate legend files
Key Parameters
The following parameters appear at the top of each script and can be adjusted without modifying any logic:

Parameter	Default	Description
analysis_level	"gene"	"gene" or "isoform"
padj_cutoff	0.01 (DE) / 0.05 (enrichment)	Adjusted p-value threshold for significance
lfc_cutoff	0	Minimum absolute log₂FC for DE calling
shrink_type	"apeglm"	LFC shrinkage method ("apeglm", "ashr", or "normal")
n_cores	4	Parallel workers for DESeq2 and fgsea
species_name	"Mus musculus"	Species for MSigDB gene set retrieval
drug_list	c("11j","KVS","11j_PlaB","KVS_PlaB","PlaB")	Drugs to process in scripts 4–6
