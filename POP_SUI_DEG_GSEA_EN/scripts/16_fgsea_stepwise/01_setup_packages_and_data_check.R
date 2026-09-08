# =============================================================================
# STEP 1 of 13 - working directory, package installation, data file check
# =============================================================================
# This is the fgsea-only POP vs SUI pipeline, split into one file per step so
# each step can be opened and run with RStudio's "Source" on its own, in
# order - useful if your professor wants to see the pipeline step by step
# instead of one long script.
#
# RUN ORDER: this is step 1 of 13 - the first one. Nothing needs to be
# sourced before it.
#
# ALL 13 STEPS, IN ORDER (source each with "Source", top to bottom, in the
# SAME R session - do not close RStudio between steps, or you will need to
# start again from step 1):
#   01_setup_packages_and_data_check.R   (this file)
#   02_POP_data_load.R
#   03_POP_DEG_DESeq2.R
#   04_POP_fgsea_KEGG.R
#   05_SUI_data_DEG_fgsea_KEGG.R
#   06_shared_pathways_POP_SUI.R
#   07_figures_volcano_heatmaps.R
#   08_figures_KEGG_barplots.R
#   09_GO_biological_process.R
#   10_GSVA.R
#   11_ORA_clusterProfiler.R
#   12_ECM_junction_panel.R
#   13_pathview_diagrams.R
#
# WHAT THIS STEP DOES: sets the working directory, installs/loads every R
# package the pipeline needs, checks the 4 required data files are present,
# and defines kegg_label() (turns a KEGG pathway number into a readable
# name) - a small helper every later step uses.
#
# WHY fgsea IS REQUIRED, NOT OPTIONAL: every pathway-level result in this
# entire pipeline (KEGG, GO) comes from fgsea alone - there is no other GSEA
# method anywhere in these 13 files. If fgsea cannot be installed, the
# script below stops with an install instruction rather than silently
# falling back to something else.
# =============================================================================

setwd("E:/POP+SUI ONE")  # adjust this one line if your folder is different
cat("Working directory set to:", getwd(), "\n\n")

# What each package is for:
#   - limma, edgeR: normalization and statistical models for RNA-seq/microarray
#   - DESeq2: primary method for finding differentially expressed genes (POP)
#   - org.Hs.eg.db + AnnotationDbi: gene ID dictionary (Entrez <-> symbol) and
#     the gene list for each KEGG pathway / GO Biological Process term
#   - GO.db: translates a GO ID into its readable name - used in step 9
#   - fgsea: the pathway enrichment engine used throughout this pipeline
#   - readxl: to read the Wei 2020 .xls file
#   - ggplot2, ggrepel, pheatmap: figures
#   - GSVA, clusterProfiler: two additional, independent enrichment methods
#     (steps 10-11) - optional extras, each step skips itself cleanly if its
#     package fails to install (typically a Windows/Rtools compilation issue)

cat("=== Checking required packages (only installs what's missing) ===\n\n")

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  cat("Installing BiocManager (the Bioconductor package installer)...\n")
  install.packages("BiocManager")
}

# fgsea is REQUIRED - the entire pathway analysis (steps 4-9, 12) depends on
# it, so it is installed alongside the other core packages rather than
# treated as an optional extra.
bioc_packages <- c("limma", "edgeR", "DESeq2", "org.Hs.eg.db", "AnnotationDbi", "GO.db", "fgsea")
for (pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("Installing (Bioconductor):", pkg, "... (this can take a few minutes)\n")
    BiocManager::install(pkg, update = FALSE, ask = FALSE)
  } else {
    cat("OK, already installed:", pkg, "\n")
  }
}
if (!requireNamespace("fgsea", quietly = TRUE)) {
  stop(
    "\n\nfgsea could not be installed. This pipeline requires it - there is\n",
    "no fallback method. Install manually with:\n",
    "  BiocManager::install(\"fgsea\")\n",
    "then run this script again.\n"
  )
}

cran_packages <- c("readxl", "ggplot2", "ggrepel", "pheatmap")
for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("Installing (CRAN):", pkg, "...\n")
    install.packages(pkg)
  } else {
    cat("OK, already installed:", pkg, "\n")
  }
}

# GSVA and clusterProfiler: independent extra methods (steps 10-11), not
# required for the core POP-vs-SUI KEGG/GO comparison (steps 2-9, 12). Some
# of their dependencies compile C/C++/Fortran code, which on Windows needs
# Rtools matching your R version - if that is missing, the install below
# will print warnings and not succeed for that one package. Rather than
# letting that crash a later step, steps 10 and 11 each check their own
# flag and skip themselves cleanly if their package is missing.
optional_pkgs <- c("GSVA", "clusterProfiler")
for (pkg in optional_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("Installing (Bioconductor):", pkg, "... (needs internet, this one time)\n")
    tryCatch(
      BiocManager::install(pkg, update = FALSE, ask = FALSE),
      error = function(e) cat("Install of", pkg, "raised an error - will be skipped:", conditionMessage(e), "\n")
    )
  } else {
    cat("OK, already installed:", pkg, "\n")
  }
}

cat("\n=== Loading packages into the session ===\n\n")
suppressMessages({
  library(limma)
  library(edgeR)
  library(DESeq2)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(GO.db)
  library(fgsea)
  library(ggplot2)
  library(pheatmap)
  library(readxl)
  library(ggrepel)
})
cat("Core packages loaded successfully.\n\n")

has_GSVA <- requireNamespace("GSVA", quietly = TRUE)
has_clusterProfiler <- requireNamespace("clusterProfiler", quietly = TRUE)
if (has_GSVA) suppressMessages(library(GSVA))
if (has_clusterProfiler) suppressMessages(library(clusterProfiler))
cat("Optional packages - GSVA:", has_GSVA, "| clusterProfiler:", has_clusterProfiler, "\n")
cat("(FALSE means that package's step (10 or 11) will be skipped - this does\n")
cat("not affect any other step.)\n\n")

# select() exists in more than one loaded package (AnnotationDbi and, if you
# have the tidyverse installed, dplyr too) - forcing which one to use.
select <- AnnotationDbi::select

dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)
cat("'results' and 'figures' folders ready inside", getwd(), "\n\n")


## --- data file check --------------------------------------------------------
required_files <- c(
  "GSE208261_raw_counts_GRCh38.p13_NCBI.tsv",
  "GSE208261_sample_metadata.csv",
  "43032_2020_144_MOESM2_ESM.xls",
  "kegg_pathway_names.csv"
)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop(
    "\n\nMISSING DATA FILE(S). The script stopped here on purpose, BEFORE\n",
    "trying to run any analysis, to avoid confusing errors later.\n\n",
    "Current working directory: ", getwd(), "\n",
    "File(s) not found:\n  - ", paste(missing_files, collapse = "\n  - "), "\n\n",
    "Fix: copy the missing file(s) into ", getwd(), "\n",
    "using exactly that name, then run the script again.\n"
  )
}
cat("=== All 4 data files were found. Proceeding. ===\n\n")

# kegg_label(): turns a KEGG pathway number (e.g. "04510") into a readable
# label ("KEGG 04510 - Focal adhesion"), using kegg_pathway_names.csv.
kegg_tab <- read.csv("kegg_pathway_names.csv", colClasses = c("character", "character"))
kegg_names <- setNames(kegg_tab$Name, kegg_tab$PATH5)
kegg_label <- function(id) {
  nm <- kegg_names[id]
  ifelse(is.na(nm), paste0("KEGG ", id), paste0("KEGG ", id, " - ", nm))
}

cat("=== STEP 1 DONE === Now Source 02_POP_data_load.R\n")
