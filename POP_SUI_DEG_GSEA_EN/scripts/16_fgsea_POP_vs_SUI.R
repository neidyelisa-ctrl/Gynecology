# =============================================================================
# POP (GSE208261, RNA-seq) vs SUI (Wei 2020, microarray) - KEGG/GO pathway
# comparison via fgsea (Bioconductor), start to finish.
# Self-contained script, run on YOUR computer, from scratch.
# =============================================================================
#
# WHAT THIS DOES:
#   Part 1-2: differential expression (DEG) in POP (DESeq2, 12 vs 12) and in
#             SUI (Wei et al. 2020's own published DEG list, 3 vs 3).
#   Part 3-4: KEGG pathway enrichment (GSEA) in POP and in SUI, via the
#             official fgsea package, on the full ranked gene list.
#   Part 5:   pathways common to both diseases - direction and significance
#             flagged for all of them, not only the ones that pass FDR<0.25.
#   Part 6:   volcano plots, expression heatmaps, shared-pathway NES plot.
#   Part 7:   KEGG pathway barplots (POP, SUI) and the targeted-panel figure.
#   Part 8:   GO Biological Process enrichment via fgsea (more specific terms
#             than KEGG's 218 pathways, including several with no KEGG
#             equivalent, e.g. "collagen fibril organization").
#   Part 9:   GSVA - a per-sample pathway activity score, a genuinely
#             different statistical angle on the same KEGG pathways.
#   Part 10:  ORA (over-representation analysis) via clusterProfiler - the
#             same style of test the STRING website's "Analysis" tab runs.
#   Part 11:  targeted panel figure for the cell-ECM/junction pathway family.
#
# One method throughout: every pathway-level result in this script (KEGG,
# GO) comes from fgsea. No other GSEA implementation is used anywhere.
#
# HOW TO RUN THIS:
#   1. Confirm these 4 files are present in your working directory below
#      (EXACT names):
#        - GSE208261_raw_counts_GRCh38.p13_NCBI.tsv
#        - 43032_2020_144_MOESM2_ESM.xls  (the MOESM2 file, not MOESM1)
#        - GSE208261_sample_metadata.csv
#        - kegg_pathway_names.csv
#   2. Open RStudio, open this file, and run it top to bottom in ONE GO using
#      "Source" - not by selecting individual lines/blocks one at a time.
#   3. The first run installs whichever packages are missing, fgsea included
#      - needs internet, can take a while. From the second run on it just
#      runs directly. fgsea itself typically finishes each pathway test in
#      seconds, not minutes - there is no long-running step in this script.
# =============================================================================


## =============================================================================
## PART 0: working directory + automatic package installation
## =============================================================================
setwd("E:/POP+SUI FGSEA")  # adjust this one line if your folder is different
cat("Working directory set to:", getwd(), "\n\n")

# What each package is for:
#   - limma, edgeR: normalization and statistical models for RNA-seq/microarray
#   - DESeq2: primary method for finding differentially expressed genes (POP)
#   - org.Hs.eg.db + AnnotationDbi: gene ID dictionary (Entrez <-> symbol) and
#     the gene list for each KEGG pathway / GO Biological Process term
#   - GO.db: translates a GO ID into its readable name - used in Part 8
#   - fgsea: the pathway enrichment engine used throughout this script
#   - readxl: to read the Wei 2020 .xls file
#   - ggplot2, ggrepel, pheatmap: figures
#   - GSVA, clusterProfiler: two additional, independent enrichment methods
#     (Parts 9-10) - optional extras, each part skips itself cleanly if its
#     package fails to install (typically a Windows/Rtools compilation issue)

cat("=== Checking required packages (only installs what's missing) ===\n\n")

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  cat("Installing BiocManager (the Bioconductor package installer)...\n")
  install.packages("BiocManager")
}

# fgsea is REQUIRED - the entire pathway analysis (Parts 3-8, 11) depends on
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
    "\n\nfgsea could not be installed. This script requires it - there is no\n",
    "fallback method. Install manually with:\n",
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

# GSVA and clusterProfiler: independent extra methods (Parts 9-10), not
# required for the core POP-vs-SUI KEGG/GO comparison (Parts 1-8, 11). Some
# of their dependencies compile C/C++/Fortran code, which on Windows needs
# Rtools matching your R version - if that is missing, the install below
# will print warnings and not succeed for that one package. Rather than
# letting that crash the whole script later, Parts 9 and 10 each check
# their own flag and skip themselves cleanly if their package is missing.
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
cat("(FALSE means that package's part (9 or 10) will be skipped - this does\n")
cat("not affect Parts 1-8 or 11.)\n\n")

# select() exists in more than one loaded package (AnnotationDbi and, if you
# have the tidyverse installed, dplyr too) - forcing which one to use.
select <- AnnotationDbi::select

dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)
cat("'results' and 'figures' folders ready inside", getwd(), "\n\n")


## =============================================================================
## PART 0b: check that the 4 data files are where they should be
## =============================================================================
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


## =============================================================================
## PART 1: load POP data (GSE208261) and build the 12x12 design
## =============================================================================
cat("\n================ PART 1: POP data (GSE208261) ================\n\n")

meta <- read.csv("GSE208261_sample_metadata.csv")
cat("Sample table (check this against the GEO metadata before trusting it blindly):\n")
print(meta)
cat("\nTissue x Group cross-tab:\n")
print(table(meta$Tissue, meta$Treatment))
cat("\n")

counts_raw <- read.delim("GSE208261_raw_counts_GRCh38.p13_NCBI.tsv", row.names = 1, check.names = FALSE)
counts_all <- as.matrix(counts_raw[, meta$GSM])
stopifnot(identical(colnames(counts_all), meta$GSM))
cat("Count matrix: 24 samples,", nrow(counts_all), "genes (by Entrez ID)\n")

# Genes come identified by number (Entrez ID) - convert to symbol for
# readability and comparability with the SUI data (which already uses names).
ann_id <- suppressWarnings(select(org.Hs.eg.db, keys = rownames(counts_all), keytype = "ENTREZID", columns = "SYMBOL"))
ann_id <- ann_id[!is.na(ann_id$SYMBOL) & !duplicated(ann_id$ENTREZID), ]
counts_all <- counts_all[rownames(counts_all) %in% ann_id$ENTREZID, ]
rownames(counts_all) <- ann_id$SYMBOL[match(rownames(counts_all), ann_id$ENTREZID)]
counts_all <- rowsum(counts_all, group = rownames(counts_all))  # sums duplicate genes
cat("After ID->symbol conversion and summing duplicates:", nrow(counts_all), "genes\n\n")

# The 12x12 design: 12 Control (6 ligament + 6 vaginal wall, COMBINED) vs
# 12 POP (vaginal wall).
group_full <- factor(meta$Treatment, levels = c("Control", "POP"))
tissue_full <- factor(meta$Tissue)
cat("Design: 12 Control (6 ligament + 6 vaginal wall) vs 12 POP (vaginal wall). Control:",
    sum(group_full == "Control"), "| POP:", sum(group_full == "POP"), "\n\n")


## =============================================================================
## PART 2: DEG in POP via DESeq2
## =============================================================================
cat("================ PART 2: differentially expressed genes (DEG) in POP ================\n\n")

counts_int <- counts_all
storage.mode(counts_int) <- "integer"  # DESeq2 requires integer counts

coldata_naive <- data.frame(row.names = meta$GSM, group = group_full)
dds_naive <- DESeqDataSetFromMatrix(countData = counts_int, colData = coldata_naive, design = ~group)

# Low-count filter: keep genes with at least 10 reads in at least 12 of the
# 24 samples (the smaller group's size) - DESeq2's own recommended practice.
dds_naive <- dds_naive[rowSums(counts(dds_naive) >= 10) >= 12, ]
cat("Genes kept after the low-count filter:", nrow(dds_naive), "\n")

dds_naive <- DESeq(dds_naive, quiet = TRUE)
res_naive <- as.data.frame(results(dds_naive, contrast = c("group", "POP", "Control"), alpha = 0.05))
res_naive$Gene <- rownames(res_naive)
res_naive <- res_naive[order(res_naive$pvalue), ]
pop_full <- res_naive[, c("Gene", "log2FoldChange", "baseMean", "stat", "pvalue", "padj")]
colnames(pop_full) <- c("Gene", "logFC", "baseMean", "stat", "P.Value", "adj.P.Val")
write.csv(pop_full, "results/POP_DESeq2_full_table.csv", row.names = FALSE)

# DEG = |log2FC|>1 AND FDR<0.05.
pop_deg <- subset(pop_full, !is.na(adj.P.Val) & adj.P.Val < 0.05 & abs(logFC) > 1)
pop_deg <- pop_deg[order(pop_deg$adj.P.Val), ]
write.csv(pop_deg, "results/POP_DEG_logFC1_FDR05.csv", row.names = FALSE)

cat("\n=== RESULT: DEG in POP (12 Control vs 12 POP) ===\n")
cat("Genes tested:", nrow(pop_full), "\n")
cat("DEG (|log2FC|>1, FDR<0.05):", nrow(pop_deg), "(",
    sum(pop_deg$logFC > 0), "up /", sum(pop_deg$logFC < 0), "down )\n\n")

# edgeR+voom+limma: a second, independent statistical method, run as a
# cross-check and to produce the moderated-t ranking fgsea uses in Part 3.
dge <- DGEList(counts = counts_all, group = group_full)
design_naive_lm <- model.matrix(~group_full)
keep_expr <- filterByExpr(dge, design_naive_lm)
dge <- dge[keep_expr, , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge, method = "TMM")
voom_fit <- voom(dge, design_naive_lm)
fit_pop <- eBayes(lmFit(voom_fit, design_naive_lm))
pop_full_limma <- topTable(fit_pop, coef = "group_fullPOP", number = Inf, sort.by = "P")
pop_full_limma$Gene <- rownames(pop_full_limma)
write.csv(pop_full_limma[, c("Gene","logFC","AveExpr","t","P.Value","adj.P.Val")],
          "results/POP_voom_limma_full_table.csv", row.names = FALSE)
pop_deg_limma <- subset(pop_full_limma, adj.P.Val < 0.05 & abs(logFC) > 1)
cat("Cross-check with edgeR+voom+limma (secondary method):", nrow(pop_deg_limma), "DEG\n\n")


## =============================================================================
## PART 3: KEGG pathway enrichment in POP via fgsea
## =============================================================================
cat("================ PART 3: GSEA (fgsea) in POP ================\n\n")
cat("What GSEA asks: even with no individual gene reaching significance, is a\n")
cat("WHOLE GROUP of genes from the same biological pathway shifted in the same\n")
cat("direction as a block, across the FULL ranked gene list (not just the\n")
cat("163 DEG)? This detects weak, distributed signal a gene-by-gene test\n")
cat("alone would miss.\n\n")

ranked_full <- fit_pop$t[, "group_fullPOP"]
ranked_genes_pop <- names(ranked_full)[order(-ranked_full)]

# For each ranked gene, find out which KEGG pathways it belongs to.
ann_pop <- suppressWarnings(select(org.Hs.eg.db, keys = ranked_genes_pop, keytype = "SYMBOL", columns = "PATH"))
ann_pop <- ann_pop[!is.na(ann_pop$PATH), ]
gs_sizes_pop <- table(ann_pop$PATH)
valid_paths_pop <- names(gs_sizes_pop)[gs_sizes_pop >= 5 & gs_sizes_pop <= 200]  # neither too small nor too large
gene_sets_pop <- split(ann_pop$SYMBOL[ann_pop$PATH %in% valid_paths_pop], ann_pop$PATH[ann_pop$PATH %in% valid_paths_pop])
cat("KEGG pathways tested (5-200 genes):", length(gene_sets_pop), "\n\n")

cat("Running fgsea on POP...\n")
set.seed(208261)
fgsea_pop <- as.data.frame(fgsea(pathways = gene_sets_pop, stats = ranked_full, minSize = 5, maxSize = 200))
fgsea_pop$leadingEdge <- sapply(fgsea_pop$leadingEdge, paste, collapse = "/")
fgsea_pop$PathwayName <- kegg_label(fgsea_pop$pathway)
fgsea_pop <- fgsea_pop[order(fgsea_pop$pval), ]
write.csv(fgsea_pop, "results/POP_fgsea_KEGG_full.csv", row.names = FALSE)

cat("=== RESULT: KEGG pathway enrichment in POP ===\n")
cat("Pathways tested:", nrow(fgsea_pop), "\n")
cat("Significant at FDR<0.25:", sum(fgsea_pop$padj < 0.25, na.rm = TRUE), "\n")
cat("Significant at FDR<0.05:", sum(fgsea_pop$padj < 0.05, na.rm = TRUE), "\n\n")


## =============================================================================
## PART 4: SUI data (Wei 2020) - DEG (already given) + KEGG enrichment via fgsea
## =============================================================================
cat("================ PART 4: SUI data (Wei 2020) ================\n\n")

wei_xls <- "43032_2020_144_MOESM2_ESM.xls"
sheet_names <- excel_sheets(wei_xls)
cat("Sheets found in the Excel file:", paste(sheet_names, collapse = ", "), "\n")
stopifnot(all(c("up_Sui_vs_Ctrl", "down_Sui_vs_Ctrl") %in% sheet_names))

# skip=17: the first 17 rows are a descriptive header written by the
# authors (not data) - skipping them to reach the real table.
up_sheet   <- read_excel(wei_xls, sheet = "up_Sui_vs_Ctrl",   skip = 17)
down_sheet <- read_excel(wei_xls, sheet = "down_Sui_vs_Ctrl", skip = 17)
cat("'up_Sui_vs_Ctrl':", nrow(up_sheet), "rows,", ncol(up_sheet), "columns\n")
cat("'down_Sui_vs_Ctrl':", nrow(down_sheet), "rows,", ncol(down_sheet), "columns\n\n")

# Direct proof this is 3 SUI vs 3 Control samples - the file's own column
# names, printed for you to check with your own eyes.
cat("Individual-sample columns found in the Wei 2020 file:\n")
print(grep("Sui[0-9]|Ctrl[0-9]", colnames(up_sheet), value = TRUE))
cat("-> 3 'SuiN' columns + 3 'CtrlN' columns = a 3-vs-3 design, confirmed\n")
cat("   directly from the file's own column names, not by assumption.\n\n")

# Fold Change sign convention: the "down" sheet lists a MAGNITUDE (always
# >=1), not a negative value - the sign comes from which sheet the gene is
# in, not from the number itself.
cat("Checking the Fold Change sign convention:\n")
cat("  'up' sheet, Fold Change range:", round(min(up_sheet$`Fold Change`), 2), "to",
    round(max(up_sheet$`Fold Change`), 2), "\n")
cat("  'down' sheet, Fold Change range:", round(min(down_sheet$`Fold Change`), 2), "to",
    round(max(down_sheet$`Fold Change`), 2), "\n")
if (min(down_sheet$`Fold Change`) >= 1) {
  cat("  CONFIRMED: the 'down' sheet only gives the MAGNITUDE (always >=1) -\n")
  cat("  the negative sign is applied manually below based on the sheet.\n\n")
}

up_sheet$Direction <- "up"; down_sheet$Direction <- "down"
sample_cols <- c("[Sui1, Sui](normalized)", "[Sui2, Sui](normalized)", "[Sui3, Sui](normalized)",
                  "[Ctrl1, Ctrl](normalized)", "[Ctrl2, Ctrl](normalized)", "[Ctrl3, Ctrl](normalized)")
keep_cols <- c("GeneSymbol", "P-value", "FDR", "Fold Change", "Direction", sample_cols)
wei_both <- rbind(as.data.frame(up_sheet[, keep_cols]), as.data.frame(down_sheet[, keep_cols]))
colnames(wei_both) <- c("GeneSymbol", "PValue", "FDR", "FoldChange", "Direction",
                         "Sui1", "Sui2", "Sui3", "Ctrl1", "Ctrl2", "Ctrl3")
wei_both <- wei_both[!is.na(wei_both$GeneSymbol), ]
cat("Combined total (both sheets, before removing duplicates):", nrow(wei_both), "probes\n")

wei_both$logFC <- ifelse(wei_both$Direction == "down",
                          -log2(wei_both$FoldChange), log2(wei_both$FoldChange))
wei_both <- wei_both[order(wei_both$PValue), ]
sui_full <- wei_both[!duplicated(wei_both$GeneSymbol), ]  # keeps the most significant probe per gene
rownames(sui_full) <- NULL
write.csv(sui_full, "results/SUI_Wei2020_full_table.csv", row.names = FALSE)
cat("Unique genes after removing duplicates:", nrow(sui_full), "(",
    sum(sui_full$Direction == "up"), "up /", sum(sui_full$Direction == "down"), "down )\n")
cat("This table IS ALREADY the SUI DEG list - the original Wei 2020 file only\n")
cat("lists the genes the authors classified as differentially expressed\n")
cat("(fold-change>=2, raw p<0.05); the full tested array was never\n")
cat("published, so SUI DEG is used as the authors published it, not re-derived.\n\n")

set.seed(2020)
sui_mat <- as.matrix(sui_full[, c("Sui1","Sui2","Sui3","Ctrl1","Ctrl2","Ctrl3")])
rownames(sui_mat) <- sui_full$GeneSymbol
sui_mat <- avereps(sui_mat, ID = rownames(sui_mat))
group_sui <- factor(c("SUI","SUI","SUI","Ctrl","Ctrl","Ctrl"), levels = c("Ctrl","SUI"))
design_sui <- model.matrix(~group_sui)
fit_sui <- eBayes(lmFit(sui_mat, design_sui))
t_obs_sui <- fit_sui$t[, 2]
ranked_genes_sui <- names(t_obs_sui)[order(-t_obs_sui)]
cat("SUI ranked list (moderated t-statistic; positive = up in SUI):", length(ranked_genes_sui), "genes\n\n")

ann_sui <- suppressWarnings(select(org.Hs.eg.db, keys = ranked_genes_sui, keytype = "SYMBOL", columns = "PATH"))
ann_sui <- ann_sui[!is.na(ann_sui$PATH), ]
gs_sizes_sui <- table(ann_sui$PATH)
valid_paths_sui <- names(gs_sizes_sui)[gs_sizes_sui >= 5 & gs_sizes_sui <= 200]
gene_sets_sui <- split(ann_sui$SYMBOL[ann_sui$PATH %in% valid_paths_sui], ann_sui$PATH[ann_sui$PATH %in% valid_paths_sui])
cat("KEGG pathways tested (5-200 genes):", length(gene_sets_sui), "\n\n")

cat("Running fgsea on SUI (preranked on the full gene list, not phenotype\n")
cat("permutation - with only 3 vs 3 samples, phenotype permutation would\n")
cat("cap resolution at p=1/choose(6,3)=0.05; ranking the full gene list has\n")
cat("no such ceiling)...\n")
set.seed(2020)
fgsea_sui <- as.data.frame(fgsea(pathways = gene_sets_sui, stats = t_obs_sui, minSize = 5, maxSize = 200))
fgsea_sui$leadingEdge <- sapply(fgsea_sui$leadingEdge, paste, collapse = "/")
fgsea_sui$PathwayName <- kegg_label(fgsea_sui$pathway)
fgsea_sui <- fgsea_sui[order(fgsea_sui$pval), ]
write.csv(fgsea_sui, "results/SUI_fgsea_KEGG_full.csv", row.names = FALSE)

cat("=== RESULT: KEGG pathway enrichment in SUI ===\n")
cat("Pathways tested:", nrow(fgsea_sui), "\n")
cat("Significant at FDR<0.25:", sum(fgsea_sui$padj < 0.25, na.rm = TRUE), "\n")
cat("Significant at FDR<0.05:", sum(fgsea_sui$padj < 0.05, na.rm = TRUE), "\n\n")


## =============================================================================
## PART 5: pathways common to POP and SUI
## =============================================================================
cat("================ PART 5: pathways common to POP and SUI ================\n\n")

kegg_common <- merge(
  data.frame(PATH = fgsea_pop$pathway, PathwayName = kegg_label(fgsea_pop$pathway),
             NES_POP = fgsea_pop$NES, padj_POP = fgsea_pop$padj),
  data.frame(PATH = fgsea_sui$pathway, NES_SUI = fgsea_sui$NES, padj_SUI = fgsea_sui$padj),
  by = "PATH"
)

# The full table: every pathway tested in BOTH diseases (not filtered to
# significant only), flagged with direction and significance columns, so
# borderline/near-significant pathways stay visible - relevant given SUI's
# small sample size (n=3 vs 3), where a pathway just above the FDR<0.25 line
# today could easily fall below it with a few more samples. Sorted by
# whichever disease's own padj is smaller, so the strongest and the
# closest-to-significant pathways float to the top together.
kegg_common_full <- kegg_common
kegg_common_full$Same_direction <- sign(kegg_common_full$NES_POP) == sign(kegg_common_full$NES_SUI)
kegg_common_full$Significant_POP_FDR025 <- kegg_common_full$padj_POP < 0.25
kegg_common_full$Significant_SUI_FDR025 <- kegg_common_full$padj_SUI < 0.25
kegg_common_full$Significant_both_FDR025 <- kegg_common_full$Significant_POP_FDR025 & kegg_common_full$Significant_SUI_FDR025
kegg_common_full <- kegg_common_full[order(pmin(kegg_common_full$padj_POP, kegg_common_full$padj_SUI)), ]
write.csv(kegg_common_full, "results/KEGG_all_common_pathways_fgsea.csv", row.names = FALSE)
cat("Saved: results/KEGG_all_common_pathways_fgsea.csv (", nrow(kegg_common_full), "pathways total )\n")
cat(" ", sum(kegg_common_full$Same_direction), "same-direction,",
    sum(!kegg_common_full$Same_direction), "opposite-direction;\n")
cat(" ", sum(kegg_common_full$Significant_POP_FDR025), "significant in POP,",
    sum(kegg_common_full$Significant_SUI_FDR025), "significant in SUI,",
    sum(kegg_common_full$Significant_both_FDR025), "significant in BOTH )\n\n")

# The small table: only the pathways significant in BOTH diseases (the
# strict double bar - padj<0.25 in POP AND in SUI).
kegg_both_significant <- subset(kegg_common_full, Significant_both_FDR025)
kegg_both_significant <- kegg_both_significant[order(pmin(kegg_both_significant$padj_POP, kegg_both_significant$padj_SUI)), ]
write.csv(kegg_both_significant, "results/KEGG_significant_in_BOTH_diseases_fgsea.csv", row.names = FALSE)
cat("Saved: results/KEGG_significant_in_BOTH_diseases_fgsea.csv (", nrow(kegg_both_significant),
    "pathways significant in both diseases, out of", nrow(kegg_common_full), "tested in both -",
    sum(kegg_both_significant$Same_direction), "of", nrow(kegg_both_significant), "are direction-concordant )\n\n")


## =============================================================================
## PART 6: figures - volcano plots, expression heatmaps, shared-pathway NES
## =============================================================================
cat("================ PART 6: figures ================\n\n")
cat("All figures are saved into the 'figures' subfolder of your working\n")
cat("directory (", getwd(), "\\figures).\n\n", sep = "")

## --- 6a: volcano plot - POP -------------------------------------------------
pop_full$sig <- "NS"
pop_full$sig[pop_full$logFC > 1 & pop_full$adj.P.Val < 0.05] <- "Up"
pop_full$sig[pop_full$logFC < -1 & pop_full$adj.P.Val < 0.05] <- "Down"
pop_full$sig <- factor(pop_full$sig, levels = c("Down", "NS", "Up"))

n_label_each <- 15
top_up_pop <- pop_full[pop_full$sig == "Up", ]
top_up_pop <- top_up_pop[order(top_up_pop$adj.P.Val), ][seq_len(min(n_label_each, nrow(top_up_pop))), ]
top_down_pop <- pop_full[pop_full$sig == "Down", ]
top_down_pop <- top_down_pop[order(top_down_pop$adj.P.Val), ][seq_len(min(n_label_each, nrow(top_down_pop))), ]
labels_pop <- rbind(top_up_pop, top_down_pop)

p_volcano_pop <- ggplot(pop_full, aes(x = logFC, y = -log10(P.Value), color = sig)) +
  geom_point(alpha = 0.6, size = 1.2) +
  scale_color_manual(values = c(Down = "#2166AC", NS = "grey75", Up = "#B2182B")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey40") +
  ggrepel::geom_text_repel(data = labels_pop, aes(label = Gene), size = 2.8,
                            max.overlaps = 100, segment.size = 0.2, show.legend = FALSE) +
  labs(title = paste0("Volcano - POP (GSE208261, 12x12), ", nrow(pop_deg), " DEG at FDR<0.05"),
       subtitle = paste0("Top ", nrow(top_up_pop), " Up and top ", nrow(top_down_pop), " Down genes labeled (by FDR)"),
       x = "log2(Fold Change)", y = "-log10(p-value)", color = NULL) +
  theme_bw() + theme(legend.position = "top")
ggsave("figures/volcano_POP.png", p_volcano_pop, width = 9, height = 7, dpi = 300)
cat("Saved: figures/volcano_POP.png\n\n")

## --- 6b: volcano plot - SUI -------------------------------------------------
sui_full$sig <- "NS"
sui_full$sig[sui_full$logFC > 1 & sui_full$FDR < 0.05] <- "Up"
sui_full$sig[sui_full$logFC < -1 & sui_full$FDR < 0.05] <- "Down"
sui_full$sig <- factor(sui_full$sig, levels = c("Down", "NS", "Up"))

top_up_sui <- sui_full[sui_full$sig == "Up", ]
top_up_sui <- top_up_sui[order(top_up_sui$FDR), ][seq_len(min(n_label_each, nrow(top_up_sui))), ]
top_down_sui <- sui_full[sui_full$sig == "Down", ]
top_down_sui <- top_down_sui[order(top_down_sui$FDR), ][seq_len(min(n_label_each, nrow(top_down_sui))), ]
labels_sui <- rbind(top_up_sui, top_down_sui)

p_volcano_sui <- ggplot(sui_full, aes(x = logFC, y = -log10(PValue), color = sig)) +
  geom_point(alpha = 0.6, size = 1.2) +
  scale_color_manual(values = c(Down = "#2166AC", NS = "grey75", Up = "#B2182B")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey40") +
  ggrepel::geom_text_repel(data = labels_sui, aes(label = GeneSymbol), size = 2.8,
                            max.overlaps = 100, segment.size = 0.2, show.legend = FALSE) +
  labs(title = paste0("Volcano - SUI (Wei 2020, 3x3), ",
                       sum(sui_full$sig != "NS"), " genes at |log2FC|>1 & FDR<0.05"),
       subtitle = paste0("Top ", nrow(top_up_sui), " Up and top ", nrow(top_down_sui), " Down genes labeled (by FDR)"),
       x = "log2(Fold Change)", y = "-log10(p-value)", color = NULL) +
  theme_bw() + theme(legend.position = "top")
ggsave("figures/volcano_SUI.png", p_volcano_sui, width = 9, height = 7, dpi = 300)
cat("Saved: figures/volcano_SUI.png\n\n")

## --- 6c: expression heatmap - top DEG, POP ----------------------------------
top_n_heatmap <- min(40, nrow(pop_deg))
top_genes_pop <- pop_deg$Gene[order(pop_deg$adj.P.Val)][seq_len(top_n_heatmap)]
mat_pop <- voom_fit$E[top_genes_pop, , drop = FALSE]
mat_pop_z <- t(scale(t(mat_pop)))
ann_col_pop <- data.frame(Group = group_full, row.names = colnames(mat_pop_z))
png("figures/heatmap_POP_topDEG.png", width = 2400, height = 3000, res = 300)
pheatmap(mat_pop_z, annotation_col = ann_col_pop, show_colnames = FALSE,
         main = paste0("Top ", top_n_heatmap, " DEG (by FDR) - POP vs Control (z-score)"),
         fontsize_row = 7, color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100))
dev.off()
cat("Saved: figures/heatmap_POP_topDEG.png\n\n")

## --- 6d: expression heatmap - top DEG, SUI ----------------------------------
top_genes_sui <- sui_full$GeneSymbol[order(sui_full$FDR)][seq_len(min(40, nrow(sui_full)))]
top_genes_sui <- intersect(top_genes_sui, rownames(sui_mat))
mat_sui <- sui_mat[top_genes_sui, , drop = FALSE]
mat_sui_z <- t(scale(t(mat_sui)))
ann_col_sui <- data.frame(Group = group_sui, row.names = colnames(mat_sui_z))
png("figures/heatmap_SUI_topDEG.png", width = 2400, height = 3000, res = 300)
pheatmap(mat_sui_z, annotation_col = ann_col_sui, show_colnames = TRUE,
         main = paste0("Top ", length(top_genes_sui), " DEG (by FDR) - SUI vs Ctrl (z-score)"),
         fontsize_row = 7, color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100))
dev.off()
cat("Saved: figures/heatmap_SUI_topDEG.png\n\n")

## --- 6e: shared-pathway NES scatter (all pathways significant in at least
## one disease, colored by whether they are significant in BOTH) -----------
scatter_data <- subset(kegg_common_full, Significant_POP_FDR025 | Significant_SUI_FDR025)
p_shared <- ggplot(scatter_data, aes(x = NES_POP, y = NES_SUI, label = PathwayName,
                                       color = Significant_both_FDR025)) +
  geom_hline(yintercept = 0, color = "grey70") + geom_vline(xintercept = 0, color = "grey70") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(size = 3) +
  scale_color_manual(values = c("TRUE" = "#B2182B", "FALSE" = "#2166AC"),
                      labels = c("TRUE" = "Significant in both (FDR<0.25)", "FALSE" = "Significant in one only")) +
  geom_text_repel(data = subset(scatter_data, Significant_both_FDR025), size = 3, max.overlaps = 30) +
  labs(title = "KEGG pathways significant in POP and/or SUI: NES side by side",
       x = "NES - POP", y = "NES - SUI", color = NULL) +
  theme_bw() + theme(legend.position = "top")
ggsave("figures/shared_pathways_NES_comparison.png", p_shared, width = 9, height = 7, dpi = 300)
cat("Saved: figures/shared_pathways_NES_comparison.png\n\n")


## =============================================================================
## PART 7: KEGG pathway barplots (POP + SUI) and shared-pathway heatmap
## =============================================================================
cat("================ PART 7: KEGG pathway barplots ================\n\n")

## --- 7a: quick-glance top-15-by-p-value barplot (POP) -----------------------
make_fgsea_barplot <- function(fgsea_res, title, n_top = 15) {
  fgsea_res <- fgsea_res[!is.na(fgsea_res$NES), ]
  d <- head(fgsea_res[order(fgsea_res$pval), ], n_top)
  d$Sig <- ifelse(d$padj < 0.05, "FDR<0.05", ifelse(d$padj < 0.25, "FDR<0.25", "NS"))
  d$Label <- factor(kegg_label(d$pathway), levels = rev(kegg_label(d$pathway)))
  ggplot(d, aes(x = NES, y = Label, fill = Sig)) +
    geom_col() +
    scale_fill_manual(values = c("FDR<0.05" = "#B2182B", "FDR<0.25" = "#F4A582", "NS" = "grey70")) +
    geom_vline(xintercept = 0, color = "grey30") +
    labs(title = title, x = "Normalized Enrichment Score (NES)", y = NULL, fill = "Significance") +
    theme_bw() + theme(axis.text.y = element_text(size = 8), plot.title = element_text(size = 12))
}
ggsave("figures/fgsea_KEGG_barplot_POP.png",
       make_fgsea_barplot(fgsea_pop, "GSEA via fgsea - top 15 KEGG pathways in POP (by p-value)"),
       width = 12, height = 6.5, dpi = 300)
cat("Saved: figures/fgsea_KEGG_barplot_POP.png (top 15 by p-value only - POP has\n")
cat(sum(fgsea_pop$padj < 0.25, na.rm = TRUE), "pathways at FDR<0.25 total, too many for one\n")
cat("readable chart; see fgsea_KEGG_barplot_POP_FDR025.png below for the complete picture)\n")

# SUI has far fewer significant pathways than POP, so unlike POP's chart
## above, filtering to padj<0.25 FIRST and only then taking the top N
## comfortably fits all of them - no pathway is cut off by an arbitrary
## top-15-by-p-value-across-everything-tested limit.
n_sig_sui <- sum(fgsea_sui$padj < 0.25, na.rm = TRUE)
ggsave("figures/fgsea_KEGG_barplot_SUI.png",
       make_fgsea_barplot(subset(fgsea_sui, padj < 0.25),
                          paste0("GSEA via fgsea - all ", n_sig_sui, " KEGG pathways significant in SUI (FDR<0.25)"),
                          n_top = n_sig_sui),
       width = 12, height = max(6.5, 0.4 * n_sig_sui + 1.5), dpi = 300)
cat("Saved: figures/fgsea_KEGG_barplot_SUI.png (all", n_sig_sui, "significant pathways shown)\n\n")

## --- 7b: KEGG barplot for POP, strong (FDR<0.05) and moderate (FDR
## 0.05-0.25) bands ranked and shown separately, so a moderate-but-real hit
## like Focal adhesion is not crowded out by the strongest hits -----------
make_fgsea_band_barplot <- function(fgsea_res, title, n_each = 15) {
  fgsea_res <- fgsea_res[!is.na(fgsea_res$NES), ]
  strong <- fgsea_res[fgsea_res$padj < 0.05, ]
  strong <- strong[order(strong$pval), ][seq_len(min(n_each, nrow(strong))), ]
  moderate <- fgsea_res[fgsea_res$padj >= 0.05 & fgsea_res$padj < 0.25, ]
  moderate <- moderate[order(moderate$pval), ][seq_len(min(n_each, nrow(moderate))), ]
  d <- rbind(strong, moderate)
  d$Band <- ifelse(d$padj < 0.05, "FDR<0.05 (strong)", "FDR 0.05-0.25 (moderate)")
  d$Band <- factor(d$Band, levels = c("FDR<0.05 (strong)", "FDR 0.05-0.25 (moderate)"))
  d$Label <- factor(kegg_label(d$pathway), levels = rev(unique(kegg_label(d$pathway))))
  ggplot(d, aes(x = NES, y = Label, fill = Band)) +
    geom_col() +
    scale_fill_manual(values = c("FDR<0.05 (strong)" = "#B2182B", "FDR 0.05-0.25 (moderate)" = "#F4A582")) +
    geom_vline(xintercept = 0, color = "grey30") +
    facet_wrap(~Band, scales = "free_y", ncol = 1) +
    labs(title = title,
         subtitle = "Top 15 pathways in each band, ranked separately, so moderate hits are not\ncrowded out by the strongest ones",
         x = "Normalized Enrichment Score (NES)", y = NULL) +
    theme_bw() + theme(axis.text.y = element_text(size = 8), legend.position = "none",
                        strip.text = element_text(size = 10, face = "bold"),
                        plot.subtitle = element_text(size = 9))
}
n_fdr025_pop <- sum(fgsea_pop$padj < 0.25, na.rm = TRUE)
p_fdr025_pop <- make_fgsea_band_barplot(fgsea_pop, paste0("GSEA via fgsea - KEGG pathways in POP (", n_fdr025_pop, " total at FDR<0.25)"))
ggsave("figures/fgsea_KEGG_barplot_POP_FDR025.png", p_fdr025_pop, width = 12, height = 11, dpi = 300)
cat("Saved: figures/fgsea_KEGG_barplot_POP_FDR025.png (", n_fdr025_pop,
    "pathways at FDR<0.25 in POP total, top 15 of each band shown )\n\n")

## --- 7c: shared-pathway NES heatmap (pathways significant in BOTH) --------
if (nrow(kegg_both_significant) > 0) {
  heatmap_mat <- as.matrix(kegg_both_significant[, c("NES_POP", "NES_SUI")])
  rownames(heatmap_mat) <- kegg_both_significant$PathwayName
  colnames(heatmap_mat) <- c("POP", "SUI")
  png("figures/heatmap_shared_KEGG_NES_fgsea.png", width = 2600, height = 2600, res = 300)
  pheatmap(heatmap_mat, cluster_cols = FALSE,
           main = "KEGG pathways significant in BOTH diseases - NES in POP vs SUI",
           fontsize_row = 9, color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100))
  dev.off()
  cat("Saved: figures/heatmap_shared_KEGG_NES_fgsea.png\n\n")
} else {
  cat("SKIPPED heatmap_shared_KEGG_NES_fgsea.png (no pathway significant in both diseases)\n\n")
}


## =============================================================================
## PART 8: GO Biological Process enrichment via fgsea
## =============================================================================
cat("================ PART 8: GO Biological Process enrichment ================\n\n")
cat("WHY THIS PART EXISTS: KEGG (Parts 1-7) covers ~218 broad pathways. GO\n")
cat("Biological Process has far more specific terms - including some that\n")
cat("map directly onto a connective-tissue/extracellular-matrix weakness\n")
cat("hypothesis (e.g. 'collagen fibril organization', 'extracellular matrix\n")
cat("organization') that KEGG has no equivalent for.\n\n")

## --- 8a: build GO Biological Process gene sets from org.Hs.eg.db ----------
cat("Building GO Biological Process gene sets from org.Hs.eg.db...\n")
ann_go_pop <- suppressWarnings(select(org.Hs.eg.db, keys = ranked_genes_pop, keytype = "SYMBOL", columns = c("GO", "ONTOLOGY")))
ann_go_pop <- ann_go_pop[!is.na(ann_go_pop$GO) & ann_go_pop$ONTOLOGY == "BP", ]
gs_sizes_go <- table(ann_go_pop$GO)
valid_go <- names(gs_sizes_go)[gs_sizes_go >= 5 & gs_sizes_go <= 200]
gene_sets_go <- split(ann_go_pop$SYMBOL[ann_go_pop$GO %in% valid_go], ann_go_pop$GO[ann_go_pop$GO %in% valid_go])
cat("GO Biological Process terms tested (5-200 genes):", length(gene_sets_go), "\n\n")

## go_term_name(): translates a GO ID into its readable name - GO.db ships
## with org.Hs.eg.db, no internet needed.
go_term_name <- function(ids) {
  nm <- suppressMessages(tryCatch(AnnotationDbi::Term(ids), error = function(e) rep(NA_character_, length(ids))))
  unname(ifelse(is.na(nm), ids, paste0(ids, " - ", nm)))
}

## --- 8b: run fgsea (GO Biological Process) on POP --------------------------
cat("Running fgsea (GO Biological Process) on POP...\n")
set.seed(208261)
go_pop <- as.data.frame(fgsea(pathways = gene_sets_go, stats = ranked_full, minSize = 5, maxSize = 200))
go_pop$GOName <- go_term_name(go_pop$pathway)
go_pop$leadingEdge <- sapply(go_pop$leadingEdge, paste, collapse = "/")
go_pop <- go_pop[order(go_pop$pval), ]
write.csv(go_pop, "results/POP_fgsea_GO_BP_full.csv", row.names = FALSE)
cat("Terms significant at FDR<0.25:", sum(go_pop$padj < 0.25, na.rm = TRUE), "of", nrow(go_pop), "\n")
cat("Terms significant at FDR<0.05:", sum(go_pop$padj < 0.05, na.rm = TRUE), "\n\n")

## --- 8c: run fgsea (GO Biological Process) on SUI --------------------------
cat("Running fgsea (GO Biological Process) on SUI...\n")
set.seed(2020)
go_sui <- as.data.frame(fgsea(pathways = gene_sets_go, stats = t_obs_sui, minSize = 5, maxSize = 200))
go_sui$GOName <- go_term_name(go_sui$pathway)
go_sui$leadingEdge <- sapply(go_sui$leadingEdge, paste, collapse = "/")
go_sui <- go_sui[order(go_sui$pval), ]
write.csv(go_sui, "results/SUI_fgsea_GO_BP_full.csv", row.names = FALSE)
cat("Terms significant at FDR<0.25:", sum(go_sui$padj < 0.25, na.rm = TRUE), "of", nrow(go_sui), "\n")
cat("Terms significant at FDR<0.05:", sum(go_sui$padj < 0.05, na.rm = TRUE), "\n\n")

## --- 8d: shared GO BP terms between POP and SUI ----------------------------
report_shared_go <- function(fdr_cut) {
  pop_sig <- subset(go_pop, padj < fdr_cut)
  sui_sig <- subset(go_sui, padj < fdr_cut)
  shared_ids <- intersect(pop_sig$pathway, sui_sig$pathway)
  cat("--- FDR <", fdr_cut, ": POP significant =", nrow(pop_sig),
      "| SUI significant =", nrow(sui_sig), "| SHARED =", length(shared_ids), "---\n")
  if (length(shared_ids) == 0) return(data.frame())
  out <- merge(pop_sig[pop_sig$pathway %in% shared_ids, c("pathway", "GOName", "NES", "padj")],
               sui_sig[sui_sig$pathway %in% shared_ids, c("pathway", "NES", "padj")],
               by = "pathway", suffixes = c("_POP", "_SUI"))
  out$Same_direction <- sign(out$NES_POP) == sign(out$NES_SUI)
  out <- out[order(out$padj_POP), ]
  cat("Same direction in both diseases:", sum(out$Same_direction, na.rm = TRUE), "of", nrow(out), "shared terms\n\n")
  out
}
shared_go_025 <- report_shared_go(0.25)
if (nrow(shared_go_025) > 0) {
  write.csv(shared_go_025, "results/shared_GO_BP_FDR025.csv", row.names = FALSE)
}

## --- 8e: direct check on connective-tissue/ECM hypothesis terms -----------
cat("Checking specific GO terms tied to a connective-tissue/ECM weakness\n")
cat("hypothesis (present here if they had 5-200 genes in the ranked list):\n\n")
ecm_go_terms <- c(
  "GO:0030198",  # extracellular matrix organization
  "GO:0030199",  # collagen fibril organization
  "GO:0030574",  # collagen catabolic process
  "GO:0007044",  # cell-substrate junction assembly
  "GO:0031012"   # extracellular matrix (a CC term, kept in case also annotated BP)
)
for (gid in ecm_go_terms) {
  row_pop <- go_pop[go_pop$pathway == gid, ]
  row_sui <- go_sui[go_sui$pathway == gid, ]
  if (nrow(row_pop) == 1 || nrow(row_sui) == 1) {
    cat(go_term_name(gid), ":\n")
    if (nrow(row_pop) == 1) {
      cat("  POP: NES =", round(row_pop$NES, 2), ", FDR =", signif(row_pop$padj, 3), "\n")
    } else {
      cat("  POP: not tested here (fewer than 5 or more than 200 genes in the ranked list)\n")
    }
    if (nrow(row_sui) == 1) {
      cat("  SUI: NES =", round(row_sui$NES, 2), ", FDR =", signif(row_sui$padj, 3), "\n")
    } else {
      cat("  SUI: not tested here (fewer than 5 or more than 200 genes in the ranked list)\n")
    }
    cat("\n")
  }
}

## --- 8f: figure - top GO Biological Process terms in POP and SUI ---------
make_go_barplot <- function(go_res, title, n_top = 15) {
  go_res <- go_res[!is.na(go_res$NES), ]
  d <- head(go_res[order(go_res$pval), ], n_top)
  d$Sig <- ifelse(d$padj < 0.05, "FDR<0.05", ifelse(d$padj < 0.25, "FDR<0.25", "NS"))
  d$Label <- factor(d$GOName, levels = rev(d$GOName))
  ggplot(d, aes(x = NES, y = Label, fill = Sig)) +
    geom_col() +
    scale_fill_manual(values = c("FDR<0.05" = "#B2182B", "FDR<0.25" = "#F4A582", "NS" = "grey70")) +
    geom_vline(xintercept = 0, color = "grey30") +
    labs(title = title, x = "Normalized Enrichment Score (NES)", y = NULL, fill = "Significance") +
    theme_bw() + theme(axis.text.y = element_text(size = 7), plot.title = element_text(size = 12))
}
ggsave("figures/GO_BP_barplot_POP.png", make_go_barplot(go_pop, "GO Biological Process (fgsea) - top terms in POP"), width = 12, height = 6.5, dpi = 300)
cat("Saved: figures/GO_BP_barplot_POP.png\n")
ggsave("figures/GO_BP_barplot_SUI.png", make_go_barplot(go_sui, "GO Biological Process (fgsea) - top terms in SUI"), width = 12, height = 6.5, dpi = 300)
cat("Saved: figures/GO_BP_barplot_SUI.png\n\n")


## =============================================================================
## PART 9: GSVA on KEGG pathways
## =============================================================================
cat("================ PART 9: GSVA on KEGG pathways ================\n\n")
cat("WHAT THIS ADDS: fgsea (Parts 3-4, 8) and ORA (Part 10) both start from\n")
cat("a group comparison (POP vs Control, SUI vs Ctrl) and ask 'is this\n")
cat("pathway shifted between the two groups'. GSVA works the other way\n")
cat("around: for EACH INDIVIDUAL SAMPLE, it scores how active each pathway\n")
cat("is, using only that sample's own gene ranking (no group label used in\n")
cat("the scoring step). Only afterwards are those per-sample scores\n")
cat("compared between groups (a plain t-test here) - a genuinely different\n")
cat("statistical approach, not a re-run of GSEA.\n\n")

if (!has_GSVA) {
  cat("GSVA is not installed/available in this R session - SKIPPING Part 9\n")
  cat("entirely (usually a Windows compilation issue - GSVA's dependencies\n")
  cat("need Rtools to build from source). This does not affect any other part.\n\n")
} else {

## GSVA's function interface changed between versions (older: gsva(expr,
## gene_sets, method="gsva", ...) directly; newer: build a gsvaParam() first,
## then call gsva(param)) - checking which interface this version has.
run_gsva <- function(expr, gene_sets) {
  if (exists("gsvaParam", where = asNamespace("GSVA"), inherits = FALSE)) {
    cat("(using the newer GSVA interface: gsvaParam() + gsva(param))\n")
    param <- GSVA::gsvaParam(expr, gene_sets, kcdf = "Gaussian")
    as.matrix(GSVA::gsva(param, verbose = FALSE))
  } else {
    cat("(using the older GSVA interface: gsva(expr, gene_sets, method='gsva'))\n")
    as.matrix(GSVA::gsva(expr, gene_sets, method = "gsva", kcdf = "Gaussian", verbose = FALSE))
  }
}

gsva_group_test <- function(gsva_mat, group_labels) {
  lv <- levels(group_labels)
  pvals <- apply(gsva_mat, 1, function(x) {
    tryCatch(t.test(x[group_labels == lv[2]], x[group_labels == lv[1]])$p.value, error = function(e) NA)
  })
  data.frame(PATH = rownames(gsva_mat), PathwayName = kegg_label(rownames(gsva_mat)),
             mean_diff = rowMeans(gsva_mat[, group_labels == lv[2], drop = FALSE]) -
                         rowMeans(gsva_mat[, group_labels == lv[1], drop = FALSE]),
             pvalue = pvals, padj = p.adjust(pvals, "BH"))
}

cat("Running GSVA on POP (24 samples x", length(gene_sets_pop), "KEGG pathways)...\n")
gsva_pop <- run_gsva(voom_fit$E, gene_sets_pop)
gsva_pop_res <- gsva_group_test(gsva_pop, group_full)
gsva_pop_res <- gsva_pop_res[order(gsva_pop_res$pvalue), ]
write.csv(gsva_pop_res, "results/POP_GSVA_KEGG.csv", row.names = FALSE)
cat("Pathways significant at FDR<0.25:", sum(gsva_pop_res$padj < 0.25, na.rm = TRUE),
    "of", nrow(gsva_pop_res), "\n")
cat("Pathways significant at FDR<0.05:", sum(gsva_pop_res$padj < 0.05, na.rm = TRUE), "\n\n")

cat("Running GSVA on SUI (6 samples x", length(gene_sets_sui), "KEGG pathways)...\n")
gsva_sui <- run_gsva(sui_mat, gene_sets_sui)
gsva_sui_res <- gsva_group_test(gsva_sui, group_sui)
gsva_sui_res <- gsva_sui_res[order(gsva_sui_res$pvalue), ]
write.csv(gsva_sui_res, "results/SUI_GSVA_KEGG.csv", row.names = FALSE)
cat("Pathways significant at FDR<0.25:", sum(gsva_sui_res$padj < 0.25, na.rm = TRUE),
    "of", nrow(gsva_sui_res), "\n")
cat("Pathways significant at FDR<0.05:", sum(gsva_sui_res$padj < 0.05, na.rm = TRUE), "\n\n")
cat("NOTE: with only 3 vs 3 samples, a per-pathway t-test on SUI has very\n")
cat("little statistical power - treat the SUI GSVA p-values as suggestive,\n")
cat("not conclusive (same caveat as everywhere else in this project on the\n")
cat("SUI side of any comparison).\n\n")

} # end if (has_GSVA) - Part 9


## =============================================================================
## PART 10: ORA (Over-Representation Analysis) via clusterProfiler
## =============================================================================
cat("================ PART 10: ORA via clusterProfiler ================\n\n")
cat("WHAT THIS ADDS: fgsea (Parts 3-4, 8) and GSVA (Part 9) both use the\n")
cat("FULL ranked gene list. ORA is a third, simpler kind of test: it takes\n")
cat("just a FIXED LIST of significant genes (POP's 163 DEG, SUI's Wei 2020\n")
cat("DEG list) and asks 'is this pathway over-represented in that list,\n")
cat("compared to what would be expected by chance from the full universe of\n")
cat("tested genes'. This is the same style of test the STRING website's\n")
cat("'Analysis' tab runs.\n\n")
cat("DESIGN CHOICE: clusterProfiler's own convenience functions (enrichKEGG(),\n")
cat("enrichGO()) contact an online server live, every time they run. Using\n")
cat("clusterProfiler's generic enricher() instead, with the SAME KEGG and GO\n")
cat("gene sets already built from org.Hs.eg.db earlier in this script, keeps\n")
cat("this fully offline.\n\n")

if (!has_clusterProfiler) {
  cat("clusterProfiler is not installed/available in this R session -\n")
  cat("SKIPPING Part 10 entirely. This does not affect any other part.\n\n")
} else {

gene_sets_to_term2gene <- function(gene_sets) {
  do.call(rbind, lapply(names(gene_sets), function(nm) data.frame(term = nm, gene = gene_sets[[nm]])))
}
kegg_term2gene <- gene_sets_to_term2gene(gene_sets_pop)
kegg_term2name <- data.frame(term = kegg_tab$PATH5, name = kegg_tab$Name)
go_term2gene <- gene_sets_to_term2gene(gene_sets_go)
go_term2name <- data.frame(term = names(gene_sets_go), name = go_term_name(names(gene_sets_go)))

## --- 10a: ORA on POP (proper universe available: all 22,253 tested genes) -
cat("Running ORA on POP's 163 DEG (universe = all", nrow(pop_full), "genes DESeq2 tested)...\n")
ora_kegg_pop <- as.data.frame(clusterProfiler::enricher(
  gene = pop_deg$Gene, universe = pop_full$Gene,
  TERM2GENE = kegg_term2gene, TERM2NAME = kegg_term2name,
  pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1))
write.csv(ora_kegg_pop, "results/POP_ORA_KEGG.csv", row.names = FALSE)
cat("KEGG pathways significant at FDR<0.25:", sum(ora_kegg_pop$p.adjust < 0.25, na.rm = TRUE),
    "of", nrow(ora_kegg_pop), "\n")

ora_go_pop <- as.data.frame(clusterProfiler::enricher(
  gene = pop_deg$Gene, universe = pop_full$Gene,
  TERM2GENE = go_term2gene, TERM2NAME = go_term2name,
  pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1))
write.csv(ora_go_pop, "results/POP_ORA_GO_BP.csv", row.names = FALSE)
cat("GO BP terms significant at FDR<0.25:", sum(ora_go_pop$p.adjust < 0.25, na.rm = TRUE),
    "of", nrow(ora_go_pop), "\n\n")

## --- 10b: ORA on SUI (NO proper universe available - see note below) -----
cat("Running ORA on SUI's DEG list (", nrow(sui_full), "genes)...\n")
cat("CAVEAT: Wei 2020 never published the full array background, only their\n")
cat("own already-significant DEG list. Without a real universe, enricher()\n")
cat("falls back to using every gene annotated in org.Hs.eg.db as the\n")
cat("background, which is NOT the actual array tested - this makes the SUI\n")
cat("ORA p-values anti-conservative (biased toward looking MORE significant\n")
cat("than they really are). Treat results/SUI_ORA_*.csv as qualitative /\n")
cat("hypothesis-generating only - a limitation of the published SUI data,\n")
cat("not of this code.\n\n")
ora_kegg_sui <- as.data.frame(clusterProfiler::enricher(
  gene = sui_full$GeneSymbol,
  TERM2GENE = kegg_term2gene, TERM2NAME = kegg_term2name,
  pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1))
write.csv(ora_kegg_sui, "results/SUI_ORA_KEGG.csv", row.names = FALSE)
cat("KEGG pathways significant at FDR<0.25:", sum(ora_kegg_sui$p.adjust < 0.25, na.rm = TRUE),
    "of", nrow(ora_kegg_sui), "(caveat above applies)\n")

ora_go_sui <- as.data.frame(clusterProfiler::enricher(
  gene = sui_full$GeneSymbol,
  TERM2GENE = go_term2gene, TERM2NAME = go_term2name,
  pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1))
write.csv(ora_go_sui, "results/SUI_ORA_GO_BP.csv", row.names = FALSE)
cat("GO BP terms significant at FDR<0.25:", sum(ora_go_sui$p.adjust < 0.25, na.rm = TRUE),
    "of", nrow(ora_go_sui), "(caveat above applies)\n\n")

} # end if (has_clusterProfiler) - Part 10


## =============================================================================
## PART 11: targeted cell-ECM / cell-junction panel (POP + SUI side by side)
## =============================================================================
cat("================ PART 11: ECM/junction targeted panel ================\n\n")
cat("A pre-specified panel of 6 KEGG pathways tied to the connective-tissue/\n")
cat("ECM weakness hypothesis this project has followed from the start\n")
cat("(Focal adhesion, ECM-receptor interaction, Regulation of actin\n")
cat("cytoskeleton, Adherens junction, Tight junction, Gap junction).\n")
cat("Multiple-testing correction is applied ONLY WITHIN this 6-pathway panel\n")
cat("(BH across n=6, not across all 218/200 KEGG pathways) - a standard,\n")
cat("legitimate practice for a small, hypothesis-driven panel decided in\n")
cat("advance from the literature, not a way to manufacture significance:\n")
cat("pathways that are genuinely null (e.g. ECM-receptor interaction in POP)\n")
cat("still come out non-significant even under this lighter correction.\n\n")

ecm_panel_ids <- c("04510", "04512", "04810", "04520", "04530", "04540")
build_ecm_panel <- function(fgsea_res, ids, disease_label) {
  d <- fgsea_res[fgsea_res$pathway %in% ids, c("pathway", "NES", "pval", "padj")]
  d$padj_panel <- p.adjust(d$pval, method = "BH")
  d$Disease <- disease_label
  d
}
ecm_pop <- build_ecm_panel(fgsea_pop, ecm_panel_ids, "POP")
ecm_sui <- build_ecm_panel(fgsea_sui, ecm_panel_ids, "SUI")
ecm_both <- rbind(ecm_pop, ecm_sui)
ecm_both$PathwayName <- kegg_label(ecm_both$pathway)
colnames(ecm_both)[colnames(ecm_both) == "padj"] <- "padj_full_218_or_200_pathways"
write.csv(ecm_both, "results/ECM_junction_targeted_panel_POP_SUI.csv", row.names = FALSE)
cat("Saved: results/ECM_junction_targeted_panel_POP_SUI.csv\n\n")

ecm_both$FDRLabel <- paste0("FDR=", signif(ecm_both$padj_panel, 2))
p_ecm <- ggplot(ecm_both, aes(x = NES, y = PathwayName, fill = Disease)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.65) +
  geom_text(aes(label = FDRLabel, hjust = ifelse(NES < 0, 1.05, -0.05)),
            position = position_dodge(width = 0.75), size = 3, color = "grey20") +
  geom_vline(xintercept = 0, color = "grey30") +
  scale_fill_manual(values = c(POP = "#B2182B", SUI = "#2166AC")) +
  scale_x_continuous(expand = expansion(mult = 0.25)) +
  labs(title = "Targeted panel: cell-ECM / cell-cell junction KEGG pathways",
       subtitle = "Pre-specified 6-pathway panel - FDR corrected within panel only (n=6)",
       x = "Normalized Enrichment Score (NES)", y = NULL, fill = NULL) +
  theme_bw() + theme(legend.position = "top", plot.subtitle = element_text(size = 10))
ggsave("figures/fgsea_KEGG_ECM_targeted_panel.png", p_ecm, width = 11, height = 6.5, dpi = 300)
cat("Saved: figures/fgsea_KEGG_ECM_targeted_panel.png\n\n")


cat("=================================================================\n")
cat("=== END OF SCRIPT ===\n")
cat("Optional packages actually available this run - GSVA:", has_GSVA,
    "| clusterProfiler:", has_clusterProfiler, "\n")
cat("Any FALSE above means the corresponding part (9 or 10) was skipped -\n")
cat("this does not affect any other part.\n\n")
cat("With thousands of GO terms tested (Part 8), expect the 'shared between\n")
cat("POP and SUI' list to look thin even when both diseases show real,\n")
cat("biologically coherent signal individually - requiring FDR<0.25 in BOTH\n")
cat("datasets independently is a much harder bar with ~2,000-4,000 tests\n")
cat("than with KEGG's 218. Check each term of interest in the POP-only and\n")
cat("SUI-only files too, not only the 'shared' file, before concluding\n")
cat("there is no overlap.\n")
cat("=================================================================\n")
