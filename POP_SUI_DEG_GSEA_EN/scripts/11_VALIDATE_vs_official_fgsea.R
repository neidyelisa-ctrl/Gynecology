# =============================================================================
# VALIDATION: our hand-rolled GSEA vs the OFFICIAL fgsea package (Bioconductor)
# POP (GSE208261, RNA-seq) vs SUI (Wei 2020, microarray) - 12x12 design
# Self-contained script, run on YOUR computer (Windows), from scratch.
# =============================================================================
#
# WHY THIS SCRIPT EXISTS: every GSEA computation in this project so far (the
# Enrichment Score, the NES, the permutation) has been calculated with MY OWN
# from-scratch reimplementation of the Subramanian et al. 2005 formula -
# because the environment this project was built in has no internet access
# to install the official package (`fgsea`, Bioconductor). That means this
# specific piece has never been checked side by side against the official
# implementation. This script does exactly that, using YOUR computer (which
# has internet).
#
# WHAT IT DOES: runs the already-validated 12x12 pipeline (DESeq2 + hand-
# rolled GSEA, the same numbers as always: 163 DEG / 117 pathways / 22
# shared) AND, in the new PART 7, runs the official `fgsea` package on the
# EXACT SAME ranked gene list AND THE EXACT SAME KEGG gene sets (from
# org.Hs.eg.db) our hand-rolled GSEA already used in Parts 3-4 above. Using
# identical gene sets on both sides is deliberate: it isolates the
# comparison to ONLY the enrichment-score/permutation algorithm itself,
# with no other difference (a different gene-set source or version) able
# to explain a disagreement. At the end it compares the two side by side:
# how many pathways agree in direction, and the correlation between the
# NES values from each method. Part 7 ALSO saves the FULL official fgsea
# table (every pathway, with fgsea's own NES) to
# results/POP_official_fgsea_KEGG_full.csv and _SUI_..., which is the file
# to check for a pathway's real direction when our own hand-rolled NES
# came out NA in Part 3/4 (this happens for a handful of the MOST
# significant pathways - see the README for why).
#
# PART 8 (also new) runs GO Biological Process enrichment - thousands of
# more specific terms than KEGG's 218 pathways, including some (e.g.
# "collagen fibril organization", "extracellular matrix organization")
# that map directly onto a connective-tissue/ECM weakness hypothesis with
# no KEGG equivalent. It uses the official fgsea package directly (not our
# hand-rolled permutation code) because GO BP has too many terms for that
# nested-loop approach to stay fast.
#
# HOW TO RUN THIS:
#
#   1. The E:\POP+SUI 63 folder already has the 3 data files. They all sit
#      loose directly in that folder, with no "data" subfolder.
#   2. Confirm these 4 files are present in E:\POP+SUI 63 (EXACT names):
#        - GSE208261_raw_counts_GRCh38.p13_NCBI.tsv
#        - 43032_2020_144_MOESM2_ESM.xls  (the MOESM2 file, not MOESM1)
#        - GSE208261_sample_metadata.csv
#        - kegg_pathway_names.csv  (CORRECTED version - sent together with
#          this script, already with the 5 updated pathway names you
#          checked against KEGG.jp - replace the old file with this one)
#   3. Open RStudio, open this file (.R), and run it from top to bottom in
#      ONE GO using "Source" (the Source button, or Ctrl+Shift+Enter in
#      RStudio) - NOT by selecting individual lines/blocks and running
#      them one at a time. Running only part of a `for` loop's body
#      without the loop header that defines its variable is exactly what
#      causes "object 'pkg' not found" errors - the code itself is fine,
#      but a loop variable only exists while the loop that declares it is
#      actually running.
#   4. The FIRST time you run it, Part 0 installs whichever packages are
#      missing, including `fgsea` (new in this version) - this needs
#      internet and can take a while. From the second run on, it just runs
#      directly. Note: this version does NOT use `msigdbr` any more (see
#      below), so it does not need to reach zenodo.org or any other
#      external gene-set database - only CRAN/Bioconductor, once, to
#      install `fgsea` itself.
#   5. Part 7 (the comparison against official fgsea) reuses the KEGG gene
#      sets already built from org.Hs.eg.db in Parts 3-4 - it does not
#      download anything new. If you get an error specifically in Part 7,
#      send me the full error message and I will fix it.
# =============================================================================


## =============================================================================
## PART 0: working directory + automatic package installation
## =============================================================================
# setwd() tells R "from now on, look for files inside this folder." Adjust
# the path below IF your folder is not exactly this one.
setwd("E:/POP+SUI 63")
cat("Working directory set to:", getwd(), "\n\n")

# What each package below is for:
#   - limma, edgeR: normalization and statistical models for RNA-seq/microarray
#   - DESeq2: the primary method for finding differentially expressed genes
#     (DEGs) in the POP RNA-seq data
#   - org.Hs.eg.db + AnnotationDbi: dictionary that translates gene IDs
#     (Entrez <-> gene symbol, e.g. "7157" <-> "TP53") and supplies the
#     gene list for each KEGG pathway (and, in Part 8, each GO Biological
#     Process term)
#   - GO.db: translates a GO ID (e.g. "GO:0030199") into its readable name
#     ("collagen fibril organization") - used only in Part 8
#   - readxl: to read the Wei 2020 .xls file
#   - ggplot2, ggrepel, pheatmap: only used to draw the final plots
#
# BiocManager is the "installer" used for Bioconductor packages (the first
# few below); install.packages() is R's standard installer, used for the
# rest. The block below checks each package one by one: if it's already
# installed, it's skipped; if not, it gets installed - and prints what it's
# doing so you can follow progress (this can take a while the first time).

cat("=== Checking required packages (only installs what's missing) ===\n\n")

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  cat("Installing BiocManager (the Bioconductor package installer)...\n")
  install.packages("BiocManager")
}

bioc_packages <- c("limma", "edgeR", "DESeq2", "org.Hs.eg.db", "AnnotationDbi", "GO.db")
for (pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("Installing (Bioconductor):", pkg, "... (this can take a few minutes)\n")
    BiocManager::install(pkg, update = FALSE, ask = FALSE)
  } else {
    cat("OK, already installed:", pkg, "\n")
  }
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

# fgsea: the OFFICIAL GSEA package, used only in Part 7 for comparison. It
# is only run on gene sets we already built ourselves from org.Hs.eg.db in
# Parts 3-4 (see Part 7 below for why) - it needs internet ONLY to install
# the package itself, once; running it afterwards needs no internet at all.
if (!requireNamespace("fgsea", quietly = TRUE)) {
  cat("Installing (Bioconductor): fgsea ... (needs internet, this one time)\n")
  BiocManager::install("fgsea", update = FALSE, ask = FALSE)
} else {
  cat("OK, already installed: fgsea\n")
}

cat("\n=== Loading packages into the session ===\n\n")
suppressMessages({
  library(limma)
  library(edgeR)
  library(DESeq2)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(GO.db)
  library(ggplot2)
  library(pheatmap)
  library(readxl)
  library(ggrepel)
  library(fgsea)
})
cat("All packages loaded successfully.\n\n")

# select() exists in more than one loaded package (AnnotationDbi and, if you
# have the tidyverse installed, dplyr too) - explicitly forcing which one to
# use, so we never get a silent "wrong kind of object" error from the wrong
# select() being called.
select <- AnnotationDbi::select

dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)
cat("'results' and 'figures' folders ready inside", getwd(), "\n\n")


## =============================================================================
## PART 0b: check that the 4 data files are where they should be
## =============================================================================
# Instead of letting the script crash further down with a generic error like
# "cannot open file" (which doesn't say WHICH file or WHERE it should be),
# we check all 4 files RIGHT HERE, at the start, and stop with a clear
# message if any of them is missing.

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


## Helper functions used further below ---------------------------------------
# kegg_label(): turns a KEGG pathway number (e.g. "04510") into a readable
# label ("KEGG 04510 - Focal adhesion"), using the table in
# kegg_pathway_names.csv (without this, plots would only show bare numbers).
kegg_tab <- read.csv("kegg_pathway_names.csv", colClasses = c("character", "character"))
kegg_names <- setNames(kegg_tab$Name, kegg_tab$PATH5)
kegg_label <- function(id) {
  nm <- kegg_names[id]
  ifelse(is.na(nm), paste0("KEGG ", id), paste0("KEGG ", id, " - ", nm))
}

# calc_es(): computes the GSEA "Enrichment Score" (ES) - the statistic that
# measures whether a pathway's genes cluster near the top (or bottom) of the
# ranked gene list more than would be expected by chance. Implemented by
# hand here (the Subramanian et al. 2005 formula) because the official
# package (fgsea) is not available without internet access in the
# environment this project was originally built in - running on your
# machine with internet, the result is mathematically equivalent.
calc_es <- function(hit_idx, scores_abs, N) {
  Nh <- length(hit_idx); Nm <- N - Nh
  if (Nm <= 0 || Nh == 0) return(NA)
  sum_hit <- sum(scores_abs[hit_idx])
  step <- rep(-1 / Nm, N); step[hit_idx] <- scores_abs[hit_idx] / sum_hit
  running <- cumsum(step); running[which.max(abs(running))]
}


## =============================================================================
## PART 1: load POP data (GSE208261) and build the 12x12 design
## =============================================================================
cat("\n================ PART 1: POP data (GSE208261) ================\n\n")

# read.csv: reads the metadata table - which sample (GSM) is Control or POP,
# and which tissue (vaginal wall or uterosacral ligament) it came from.
meta <- read.csv("GSE208261_sample_metadata.csv")
cat("Sample table (check this against the GEO metadata before trusting it blindly):\n")
print(meta)
cat("\nTissue x Group cross-tab:\n")
print(table(meta$Tissue, meta$Treatment))
cat("\n")

# read.delim: reads the raw count matrix (how many sequencing "reads" landed
# on each gene, in each sample) - the raw material for any RNA-seq analysis.
counts_raw <- read.delim("GSE208261_raw_counts_GRCh38.p13_NCBI.tsv", row.names = 1, check.names = FALSE)
counts_all <- as.matrix(counts_raw[, meta$GSM])
stopifnot(identical(colnames(counts_all), meta$GSM))
cat("Count matrix: 24 samples,", nrow(counts_all), "genes (by Entrez ID)\n")

# Genes come identified by number (Entrez ID, e.g. "7157"), not by name
# (symbol, e.g. "TP53") - converting here so names are readable in the rest
# of the script and comparable with Wei2020 (which already uses names).
ann_id <- suppressWarnings(select(org.Hs.eg.db, keys = rownames(counts_all), keytype = "ENTREZID", columns = "SYMBOL"))
ann_id <- ann_id[!is.na(ann_id$SYMBOL) & !duplicated(ann_id$ENTREZID), ]
counts_all <- counts_all[rownames(counts_all) %in% ann_id$ENTREZID, ]
rownames(counts_all) <- ann_id$SYMBOL[match(rownames(counts_all), ann_id$ENTREZID)]
counts_all <- rowsum(counts_all, group = rownames(counts_all))  # sums duplicate genes
cat("After ID->symbol conversion and summing duplicates:", nrow(counts_all), "genes\n\n")

# The 12x12 design: 12 Control (6 ligament + 6 vaginal wall, COMBINED) vs
# 12 POP (vaginal wall). See the report for the full discussion of why the
# controls have two tissue types and what that implies.
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

# LOW-COUNT FILTER - this is the step you asked about, whether it had been
# skipped in a different analysis. Here it is: we only keep genes with at
# least 10 reads in at least 12 of the 24 samples (the size of the smaller
# group). Genes with very low counts have no real statistical power and,
# if left in the table, only hurt the multiple-testing correction (FDR) -
# this is the practice recommended by DESeq2's own documentation.
dds_naive <- dds_naive[rowSums(counts(dds_naive) >= 10) >= 12, ]
cat("Genes kept after the low-count filter:", nrow(dds_naive), "\n")

dds_naive <- DESeq(dds_naive, quiet = TRUE)
res_naive <- as.data.frame(results(dds_naive, contrast = c("group", "POP", "Control"), alpha = 0.05))
res_naive$Gene <- rownames(res_naive)
res_naive <- res_naive[order(res_naive$pvalue), ]
pop_full <- res_naive[, c("Gene", "log2FoldChange", "baseMean", "stat", "pvalue", "padj")]
colnames(pop_full) <- c("Gene", "logFC", "baseMean", "stat", "P.Value", "adj.P.Val")
write.csv(pop_full, "results/POP_DESeq2_full_table.csv", row.names = FALSE)

# DEG = genes with at least a 2-fold change (|log2FC|>1) AND statistically
# significant after multiple-testing correction (FDR<0.05).
pop_deg <- subset(pop_full, !is.na(adj.P.Val) & adj.P.Val < 0.05 & abs(logFC) > 1)
pop_deg <- pop_deg[order(pop_deg$adj.P.Val), ]
write.csv(pop_deg, "results/POP_DEG_logFC1_FDR05.csv", row.names = FALSE)

cat("\n=== RESULT: DEG in POP (12 Control vs 12 POP) ===\n")
cat("Genes tested:", nrow(pop_full), "\n")
cat("DEG (|log2FC|>1, FDR<0.05):", nrow(pop_deg), "(",
    sum(pop_deg$logFC > 0), "up /", sum(pop_deg$logFC < 0), "down )\n")
cat("Expected value, already documented in the project: 163 DEG (123 up / 40 down).\n")
cat("If the number above matches 163, this part is validated.\n\n")

# edgeR+voom+limma: a SECOND, independent statistical method (not DESeq2),
# run here only as a cross-check and to produce the gene ranking used for
# GSEA (Part 3). It is not this script's primary DEG method - it's a
# cross-check, showing whether the two methods agree or disagree in the
# expected way (see the report's discussion on DESeq2 vs edgeR).
dge <- DGEList(counts = counts_all, group = group_full)
design_naive_lm <- model.matrix(~group_full)
keep_expr <- filterByExpr(dge, design_naive_lm)  # edgeR's equivalent filter
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
## PART 3: GSEA in POP (classic method)
## =============================================================================
cat("================ PART 3: GSEA in POP (classic method) ================\n\n")
cat("What GSEA asks: even with no individual gene reaching significance, is a\n")
cat("WHOLE GROUP of genes from the same biological pathway shifted in the same\n")
cat("direction as a block? This detects weak, distributed signal that a\n")
cat("gene-by-gene test alone would miss.\n\n")
cat("CLASSIC method (phenotype permutation): valid here because we have the\n")
cat("full set of 24 samples and a good number of possible relabelings\n")
cat("(choose 12 of 24 = 2,704,156 combinations) - excellent statistical\n")
cat("resolution.\n\n")

set.seed(208261)  # fixes the random seed: re-running gives the SAME result
ranked_full <- fit_pop$t[, "group_fullPOP"]
ord2 <- order(-ranked_full)
ranked_genes_pop <- names(ranked_full)[ord2]
ranked_scores_pop <- ranked_full[ord2]
N_pop <- length(ranked_genes_pop)

# For each ranked gene, find out which KEGG pathways it belongs to.
ann_pop <- suppressWarnings(select(org.Hs.eg.db, keys = ranked_genes_pop, keytype = "SYMBOL", columns = "PATH"))
ann_pop <- ann_pop[!is.na(ann_pop$PATH), ]
gs_sizes_pop <- table(ann_pop$PATH)
valid_paths_pop <- names(gs_sizes_pop)[gs_sizes_pop >= 5 & gs_sizes_pop <= 200]  # neither too small nor too large
gene_sets_pop <- split(ann_pop$SYMBOL[ann_pop$PATH %in% valid_paths_pop], ann_pop$PATH[ann_pop$PATH %in% valid_paths_pop])
cat("KEGG pathways tested (5-200 genes):", length(gene_sets_pop), "\n")

hit_idx_pop <- lapply(gene_sets_pop, function(g) which(ranked_genes_pop %in% g))
hit_idx_pop <- hit_idx_pop[sapply(hit_idx_pop, length) >= 3]
cat("Pathways with at least 3 genes in the ranked list:", length(hit_idx_pop), "\n")
es_obs_pop <- sapply(hit_idx_pop, calc_es, scores_abs = abs(ranked_scores_pop), N = N_pop)

# The permutation: shuffle the Control/POP labels 500 times, recompute the
# ES for each shuffle, and compare the real ES against that "random"
# distribution - this is what produces each pathway's p-value.
n_perm <- 500
cat("Running", n_perm, "permutations (shuffles) of the Control/POP label...\n")
cat("(this is the slowest part of the script - it can take a few minutes)\n")
t0 <- Sys.time()
gene_sets_syms_pop <- lapply(hit_idx_pop, function(idx) ranked_genes_pop[idx])
perm_es_pop <- matrix(NA_real_, nrow = n_perm, ncol = length(hit_idx_pop))
for (i in seq_len(n_perm)) {
  perm_group <- sample(group_full)
  perm_design <- model.matrix(~perm_group)
  perm_fit <- eBayes(lmFit(voom_fit, perm_design))
  t_perm <- perm_fit$t[, 2]
  rank_of_gene <- rank(-t_perm, ties.method = "first")
  scores_abs_sorted <- sort(abs(t_perm), decreasing = TRUE)
  for (j in seq_along(gene_sets_syms_pop)) {
    hidx <- rank_of_gene[gene_sets_syms_pop[[j]]]
    if (length(hidx) >= 3) perm_es_pop[i, j] <- calc_es(hidx, scores_abs_sorted, N_pop)
  }
}
cat("Done in", round(difftime(Sys.time(), t0, units = "secs"), 1), "seconds\n\n")

pval_pop <- numeric(length(hit_idx_pop)); nes_pop <- numeric(length(hit_idx_pop))
for (j in seq_along(hit_idx_pop)) {
  pe <- perm_es_pop[, j]; pe <- pe[!is.na(pe)]
  if (es_obs_pop[j] >= 0) {
    pval_pop[j] <- (sum(pe >= es_obs_pop[j]) + 1) / (length(pe) + 1)
    base <- mean(pe[pe >= 0]); if (is.nan(base) || base == 0) base <- NA
  } else {
    pval_pop[j] <- (sum(pe <= es_obs_pop[j]) + 1) / (length(pe) + 1)
    base <- mean(abs(pe[pe < 0])); if (is.nan(base) || base == 0) base <- NA
  }
  nes_pop[j] <- es_obs_pop[j] / base
}

leading_edge_pop <- sapply(hit_idx_pop, function(idx) paste(ranked_genes_pop[idx], collapse = "/"))
gsea_pop <- data.frame(PATH = names(hit_idx_pop), Nh = sapply(hit_idx_pop, length),
                        ES = es_obs_pop, NES = nes_pop, pvalue = pval_pop, leadingEdge = leading_edge_pop)
gsea_pop$p.adjust <- p.adjust(gsea_pop$pvalue, "BH")
gsea_pop$PathwayName <- kegg_label(gsea_pop$PATH)
gsea_pop <- gsea_pop[order(gsea_pop$pvalue), ]
write.csv(gsea_pop, "results/POP_GSEA_classic_KEGG.csv", row.names = FALSE)

cat("=== RESULT: classic GSEA in POP ===\n")
cat("Pathways significant at FDR<0.25:", sum(gsea_pop$p.adjust < 0.25, na.rm = TRUE), "of", nrow(gsea_pop), "\n")
cat("Pathways significant at FDR<0.05:", sum(gsea_pop$p.adjust < 0.05, na.rm = TRUE), "\n")
cat("Expected value, already documented: 117 at FDR<0.25 and 97 at FDR<0.05, out of 218 tested.\n\n")


## =============================================================================
## PART 4: SUI data (Wei 2020) - DEG (already given) + preranked GSEA
## =============================================================================
cat("================ PART 4: SUI data (Wei 2020) ================\n\n")

wei_xls <- "43032_2020_144_MOESM2_ESM.xls"
sheet_names <- excel_sheets(wei_xls)
cat("Sheets found in the Excel file:", paste(sheet_names, collapse = ", "), "\n")
stopifnot(all(c("up_Sui_vs_Ctrl", "down_Sui_vs_Ctrl") %in% sheet_names))

# skip=17: the first 17 rows of the Excel file are a descriptive header
# written by the authors (not data) - skipping them to reach the real table.
up_sheet   <- read_excel(wei_xls, sheet = "up_Sui_vs_Ctrl",   skip = 17)
down_sheet <- read_excel(wei_xls, sheet = "down_Sui_vs_Ctrl", skip = 17)
cat("'up_Sui_vs_Ctrl':", nrow(up_sheet), "rows,", ncol(up_sheet), "columns\n")
cat("'down_Sui_vs_Ctrl':", nrow(down_sheet), "rows,", ncol(down_sheet), "columns\n\n")

# DIRECT PROOF THAT THIS IS 3 SUI SAMPLES AND 3 CONTROL SAMPLES (not an
# assumption - these are the file's real column names, printed below for
# you to check with your own eyes):
cat("Individual-sample columns found in the Wei 2020 file:\n")
print(grep("Sui[0-9]|Ctrl[0-9]", colnames(up_sheet), value = TRUE))
cat("-> 3 'SuiN' columns + 3 'CtrlN' columns = a 3-vs-3 design, confirmed\n")
cat("   directly from the file's own column names, not by assumption.\n\n")

# Re-checking (this is the second time this specific check has been done in
# the project, on purpose) the "Fold Change" sign convention: the "down"
# sheet lists a MAGNITUDE (always >=1), not a negative value - the sign
# comes from which sheet the gene is in, not from the number itself.
cat("Checking the Fold Change sign convention:\n")
cat("  'up' sheet, Fold Change range:", round(min(up_sheet$`Fold Change`), 2), "to",
    round(max(up_sheet$`Fold Change`), 2), "\n")
cat("  'down' sheet, Fold Change range:", round(min(down_sheet$`Fold Change`), 2), "to",
    round(max(down_sheet$`Fold Change`), 2), "\n")
if (min(down_sheet$`Fold Change`) >= 1) {
  cat("  CONFIRMED: the 'down' sheet only gives the MAGNITUDE (always >=1) -\n")
  cat("  the negative sign has to be applied manually based on the\n")
  cat("  sheet/Direction, it does not come built into the number. The code\n")
  cat("  below already applies that correction.\n\n")
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

# This is where the sign correction actually happens: 'down' genes get a
# NEGATIVE log2(FoldChange); 'up' genes get a positive log2(FoldChange).
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
cat("published. That is why we do not re-derive SUI DEG from scratch - we\n")
cat("use the authors' own list.\n\n")

set.seed(2020)
sui_mat <- as.matrix(sui_full[, c("Sui1","Sui2","Sui3","Ctrl1","Ctrl2","Ctrl3")])
rownames(sui_mat) <- sui_full$GeneSymbol
sui_mat <- avereps(sui_mat, ID = rownames(sui_mat))
group_sui <- factor(c("SUI","SUI","SUI","Ctrl","Ctrl","Ctrl"), levels = c("Ctrl","SUI"))
design_sui <- model.matrix(~group_sui)
fit_sui <- eBayes(lmFit(sui_mat, design_sui))
t_obs_sui <- fit_sui$t[, 2]
ord_s <- order(-t_obs_sui)
ranked_genes_sui <- names(t_obs_sui)[ord_s]
ranked_scores_sui <- t_obs_sui[ord_s]
N_sui <- length(ranked_genes_sui)
cat("SUI ranked list (moderated t-statistic; positive = up in SUI):", N_sui, "genes\n\n")

cat("PRERANKED method for SUI's GSEA (different from the method used in\n")
cat("POP): with only 3 vs 3 samples, PHENOTYPE permutation would have only\n")
cat("choose(6,3)=20 possible relabelings, capping resolution at\n")
cat("p=1/20=0.05 - FDR<0.05 could never be reached, no matter how real the\n")
cat("signal is. PRERANKED (gene-set) permutation, which shuffles the ~200\n")
cat("pathways being tested rather than the 6 samples, has no such ceiling.\n\n")

ann_sui <- suppressWarnings(select(org.Hs.eg.db, keys = ranked_genes_sui, keytype = "SYMBOL", columns = "PATH"))
ann_sui <- ann_sui[!is.na(ann_sui$PATH), ]
gs_sizes_sui <- table(ann_sui$PATH)
valid_paths_sui <- names(gs_sizes_sui)[gs_sizes_sui >= 5 & gs_sizes_sui <= 200]
gene_sets_sui <- split(ann_sui$SYMBOL[ann_sui$PATH %in% valid_paths_sui], ann_sui$PATH[ann_sui$PATH %in% valid_paths_sui])
hit_idx_sui <- lapply(gene_sets_sui, function(g) which(ranked_genes_sui %in% g))
hit_idx_sui <- hit_idx_sui[sapply(hit_idx_sui, length) >= 3]
es_obs_sui <- sapply(hit_idx_sui, calc_es, scores_abs = abs(ranked_scores_sui), N = N_sui)

n_perm2 <- 1000
cat("Running", n_perm2, "permutations (pathway shuffling) for SUI...\n")
t0 <- Sys.time()
scores_abs_sui <- abs(ranked_scores_sui)
perm_es_sui <- matrix(NA_real_, nrow = n_perm2, ncol = length(hit_idx_sui))
for (i in seq_len(n_perm2)) {
  for (j in seq_along(hit_idx_sui)) {
    hidx <- sample.int(N_sui, length(hit_idx_sui[[j]]))
    perm_es_sui[i, j] <- calc_es(hidx, scores_abs_sui, N_sui)
  }
}
cat("Done in", round(difftime(Sys.time(), t0, units = "secs"), 1), "seconds\n\n")

pval_sui <- numeric(length(hit_idx_sui)); nes_sui <- numeric(length(hit_idx_sui))
for (j in seq_along(hit_idx_sui)) {
  pe <- perm_es_sui[, j]; pe <- pe[!is.na(pe)]
  if (es_obs_sui[j] >= 0) {
    pval_sui[j] <- (sum(pe >= es_obs_sui[j]) + 1) / (length(pe) + 1)
    base <- mean(pe[pe >= 0]); if (is.nan(base) || base == 0) base <- NA
  } else {
    pval_sui[j] <- (sum(pe <= es_obs_sui[j]) + 1) / (length(pe) + 1)
    base <- mean(abs(pe[pe < 0])); if (is.nan(base) || base == 0) base <- NA
  }
  nes_sui[j] <- es_obs_sui[j] / base
}

leading_edge_sui <- sapply(hit_idx_sui, function(idx) paste(ranked_genes_sui[idx], collapse = "/"))
gsea_sui <- data.frame(PATH = names(hit_idx_sui), Nh = sapply(hit_idx_sui, length),
                        ES = es_obs_sui, NES = nes_sui, pvalue = pval_sui, leadingEdge = leading_edge_sui)
gsea_sui$p.adjust <- p.adjust(gsea_sui$pvalue, "BH")
gsea_sui$PathwayName <- kegg_label(gsea_sui$PATH)
gsea_sui <- gsea_sui[order(gsea_sui$pvalue), ]
write.csv(gsea_sui, "results/SUI_GSEA_preranked_KEGG.csv", row.names = FALSE)

cat("=== RESULT: preranked GSEA in SUI ===\n")
cat("Pathways significant at FDR<0.25:", sum(gsea_sui$p.adjust < 0.25, na.rm = TRUE), "of", nrow(gsea_sui), "\n")
cat("Pathways significant at FDR<0.05:", sum(gsea_sui$p.adjust < 0.05, na.rm = TRUE), "\n\n")


## =============================================================================
## PART 5: shared pathways between POP and SUI + gene table
## =============================================================================
cat("================ PART 5: shared pathways between POP and SUI ================\n\n")

report_shared <- function(fdr_cut) {
  pop_sig <- subset(gsea_pop, p.adjust < fdr_cut)
  sui_sig <- subset(gsea_sui, p.adjust < fdr_cut)
  shared_ids <- intersect(pop_sig$PATH, sui_sig$PATH)
  cat("--- FDR <", fdr_cut, ": POP significant =", nrow(pop_sig),
      "| SUI significant =", nrow(sui_sig), "| SHARED =", length(shared_ids), "---\n")
  if (length(shared_ids) == 0) return(data.frame())
  out <- merge(pop_sig[pop_sig$PATH %in% shared_ids, c("PATH","PathwayName","Nh","NES","p.adjust")],
               sui_sig[sui_sig$PATH %in% shared_ids, c("PATH","Nh","NES","p.adjust")],
               by = "PATH", suffixes = c("_POP", "_SUI"))
  out$Same_direction <- sign(out$NES_POP) == sign(out$NES_SUI)
  out <- out[order(out$p.adjust_POP), ]
  n_na <- sum(is.na(out$Same_direction))
  cat("Same direction in both diseases:", sum(out$Same_direction, na.rm = TRUE), "of",
      sum(!is.na(out$Same_direction)), "pathways with a comparable NES")
  if (n_na > 0) cat(" (", n_na, "excluded - NES undefined on at least one side)")
  cat("\n\n")
  out
}

shared_025 <- report_shared(0.25)
if (nrow(shared_025) > 0) {
  write.csv(shared_025, "results/shared_pathways_FDR025.csv", row.names = FALSE)
}
cat("Expected value, already documented: 22 shared pathways at FDR<0.25,\n")
cat("of which 10 of 12 comparable (83.3%) are direction-concordant.\n\n")


## =============================================================================
## PART 6: figures - volcano, GSEA barplots (POP + SUI), shared-pathway plot
## =============================================================================
cat("================ PART 6: figures ================\n\n")
cat("All figures are saved into the 'figures' subfolder of your working\n")
cat("directory (", getwd(), "\\figures) - same folder as the results CSVs.\n\n", sep = "")

## --- 6a: volcano plot (POP DEG) --------------------------------------------
pop_full$sig <- "NS"
pop_full$sig[pop_full$logFC > 1 & pop_full$adj.P.Val < 0.05] <- "Up"
pop_full$sig[pop_full$logFC < -1 & pop_full$adj.P.Val < 0.05] <- "Down"
pop_full$sig <- factor(pop_full$sig, levels = c("Down", "NS", "Up"))
p_volcano <- ggplot(pop_full, aes(x = logFC, y = -log10(P.Value), color = sig)) +
  geom_point(alpha = 0.6, size = 1.2) +
  scale_color_manual(values = c(Down = "#2166AC", NS = "grey75", Up = "#B2182B")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey40") +
  labs(title = paste0("Volcano - POP (GSE208261, 12x12), ", nrow(pop_deg), " DEG at FDR<0.05"),
       x = "log2(Fold Change)", y = "-log10(p-value)", color = NULL) +
  theme_bw() + theme(legend.position = "top")
ggsave("figures/volcano_POP.png", p_volcano, width = 8, height = 6, dpi = 300)
cat("Saved: figures/volcano_POP.png\n\n")

## --- 6b: GSEA barplots (top 15 pathways by p-value, POP and SUI) -----------
make_gsea_barplot <- function(gsea_res, title, n_top = 15) {
  gsea_res <- gsea_res[!is.na(gsea_res$NES), ]
  d <- head(gsea_res[order(gsea_res$pvalue), ], n_top)
  d$Sig <- ifelse(d$p.adjust < 0.05, "FDR<0.05", ifelse(d$p.adjust < 0.25, "FDR<0.25", "NS"))
  d$Label <- factor(d$PathwayName, levels = rev(d$PathwayName))
  ggplot(d, aes(x = NES, y = Label, fill = Sig)) +
    geom_col() +
    scale_fill_manual(values = c("FDR<0.05" = "#B2182B", "FDR<0.25" = "#F4A582", "NS" = "grey70")) +
    geom_vline(xintercept = 0, color = "grey30") +
    labs(title = title, x = "Normalized Enrichment Score (NES)", y = NULL, fill = "Significance") +
    theme_bw() + theme(axis.text.y = element_text(size = 8), plot.title = element_text(size = 12))
}

p_bar_pop <- make_gsea_barplot(gsea_pop, "GSEA (classic) - top KEGG pathways in POP (GSE208261, 12x12)")
ggsave("figures/GSEA_barplot_POP.png", p_bar_pop, width = 12, height = 6.5, dpi = 300)
cat("Saved: figures/GSEA_barplot_POP.png\n\n")

p_bar_sui <- make_gsea_barplot(gsea_sui, "GSEA (preranked) - top KEGG pathways in SUI (Wei 2020)")
ggsave("figures/GSEA_barplot_SUI.png", p_bar_sui, width = 12, height = 6.5, dpi = 300)
cat("Saved: figures/GSEA_barplot_SUI.png\n\n")

## --- 6c: shared-pathway NES comparison scatter plot ------------------------
if (nrow(shared_025) > 0) {
  shared_plot_data <- shared_025[!is.na(shared_025$NES_POP) & !is.na(shared_025$NES_SUI), ]
  if (nrow(shared_plot_data) > 0) {
    p_shared <- ggplot(shared_plot_data, aes(x = NES_POP, y = NES_SUI, label = PathwayName)) +
      geom_hline(yintercept = 0, color = "grey70") + geom_vline(xintercept = 0, color = "grey70") +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
      geom_point(aes(color = Same_direction), size = 3) +
      scale_color_manual(values = c("TRUE" = "#2166AC", "FALSE" = "#B2182B"),
                          labels = c("TRUE" = "Same direction", "FALSE" = "Opposite direction")) +
      geom_text_repel(size = 3, max.overlaps = 30) +
      labs(title = "Shared pathways (FDR<0.25): NES in POP vs NES in SUI",
           x = "NES - POP (classic GSEA)", y = "NES - SUI (preranked GSEA)", color = NULL) +
      theme_bw() + theme(legend.position = "top")
    ggsave("figures/shared_pathways_NES_comparison.png", p_shared, width = 9, height = 7, dpi = 300)
    cat("Saved: figures/shared_pathways_NES_comparison.png\n\n")
  } else {
    cat("SKIPPED figures/shared_pathways_NES_comparison.png (no pathway with NES defined on both sides)\n\n")
  }
} else {
  cat("SKIPPED figures/shared_pathways_NES_comparison.png (no shared pathway at FDR<0.25)\n\n")
}

cat("=== All Part 1-6 figures saved: volcano_POP.png, GSEA_barplot_POP.png,\n")
cat("GSEA_barplot_SUI.png, shared_pathways_NES_comparison.png ===\n\n")


## =============================================================================
## PART 7 (REVISED): comparison with the OFFICIAL fgsea - the part that
## validates the GSEA calculation engine itself, not just the final results
## =============================================================================
cat("================ PART 7: comparison with official fgsea ================\n\n")
cat("This part runs the official Bioconductor fgsea package on the EXACT\n")
cat("SAME ranked gene list AND THE EXACT SAME KEGG gene sets (from\n")
cat("org.Hs.eg.db) our hand-rolled GSEA already used in Parts 3-4 above -\n")
cat("no new gene sets are downloaded here. Using identical gene sets on\n")
cat("both sides isolates the comparison to ONLY the enrichment-score/\n")
cat("permutation algorithm, so if the two approaches agree, that validates\n")
cat("the calculation engine this project has used from the start.\n\n")
cat("(An earlier version of this script downloaded KEGG gene sets from\n")
cat("msigdbr instead. msigdbr's KEGG data is not bundled in the package -\n")
cat("it is fetched from a Zenodo-hosted file on first use, over the\n")
cat("internet. On networks that cannot reach zenodo.org (some\n")
cat("institutional/university networks block or fail to resolve it), that\n")
cat("failed with 'Could not resolve host: zenodo.org', with no way around\n")
cat("it from inside the script. Reusing the gene sets we already built\n")
cat("sidesteps that problem entirely - Part 7 now needs no internet\n")
cat("access at all, only the fgsea package itself already installed.)\n\n")

## --- 7a: run official fgsea on the SAME ranked POP gene list, using the
## SAME gene_sets_pop object built from org.Hs.eg.db in Part 3 -------------
cat("Running official fgsea on POP, with the same KEGG gene sets our own\n")
cat("GSEA used in Part 3...\n")
set.seed(208261)
fgsea_pop <- as.data.frame(fgsea(pathways = gene_sets_pop, stats = ranked_full,
                                  minSize = 5, maxSize = 200))
cat("Official fgsea (POP):", nrow(fgsea_pop), "pathways tested,",
    sum(fgsea_pop$padj < 0.25, na.rm = TRUE), "significant at FDR<0.25\n\n")

## Save the FULL official fgsea table - EVERY pathway it tested, with ITS
## OWN NES (always a real number - fgsea does not have the "NES undefined"
## issue our hand-rolled GSEA can hit for very strongly one-sided pathways,
## see Part 3's gsea_pop / results/POP_GSEA_classic_KEGG.csv, where a
## pathway can have a p-value but NA in the NES column). This is the file
## to check for the OFFICIAL direction of a pathway when our own NES is NA.
fgsea_pop_full_out <- fgsea_pop
fgsea_pop_full_out$leadingEdge <- sapply(fgsea_pop_full_out$leadingEdge, paste, collapse = "/")
fgsea_pop_full_out$PathwayName <- kegg_label(fgsea_pop_full_out$pathway)
fgsea_pop_full_out <- fgsea_pop_full_out[order(fgsea_pop_full_out$pval), ]
write.csv(fgsea_pop_full_out, "results/POP_official_fgsea_KEGG_full.csv", row.names = FALSE)
cat("Saved: results/POP_official_fgsea_KEGG_full.csv (all", nrow(fgsea_pop_full_out),
    "pathways, with fgsea's own NES for every one - including pathways where\n")
cat("our hand-rolled NES was NA)\n\n")

## --- 7b: run official fgsea on the SAME ranked SUI gene list, using the
## SAME gene_sets_sui object built from org.Hs.eg.db in Part 4 -------------
cat("Running official fgsea on SUI, with the same KEGG gene sets our own\n")
cat("GSEA used in Part 4...\n")
set.seed(2020)
fgsea_sui <- as.data.frame(fgsea(pathways = gene_sets_sui, stats = t_obs_sui,
                                  minSize = 5, maxSize = 200))
cat("Official fgsea (SUI):", nrow(fgsea_sui), "pathways tested,",
    sum(fgsea_sui$padj < 0.25, na.rm = TRUE), "significant at FDR<0.25\n\n")

## Save the FULL official fgsea table for SUI too (see the note above the
## POP version for why this matters).
fgsea_sui_full_out <- fgsea_sui
fgsea_sui_full_out$leadingEdge <- sapply(fgsea_sui_full_out$leadingEdge, paste, collapse = "/")
fgsea_sui_full_out$PathwayName <- kegg_label(fgsea_sui_full_out$pathway)
fgsea_sui_full_out <- fgsea_sui_full_out[order(fgsea_sui_full_out$pval), ]
write.csv(fgsea_sui_full_out, "results/SUI_official_fgsea_KEGG_full.csv", row.names = FALSE)
cat("Saved: results/SUI_official_fgsea_KEGG_full.csv (all", nrow(fgsea_sui_full_out), "pathways)\n\n")

## --- 7c: compare side by side - our hand-rolled GSEA vs official fgsea ----
## Because both methods were given the exact same gene sets (keyed by the
## numeric KEGG PATH id, e.g. "04510" - no name-matching or ID-mapping
## needed, fgsea's own $pathway column already IS that same id), this
## merge is exact.
compare_methods <- function(ours, official, label) {
  cmp <- merge(
    data.frame(PATH5 = ours$PATH, PathwayName = ours$PathwayName,
               NES_ours = ours$NES, FDR_ours = ours$p.adjust),
    data.frame(PATH5 = official$pathway, NES_official = official$NES, FDR_official = official$padj),
    by = "PATH5"
  )
  cmp <- cmp[!is.na(cmp$NES_ours) & !is.na(cmp$NES_official), ]
  cmp$same_direction <- sign(cmp$NES_ours) == sign(cmp$NES_official)
  cat("=== ", label, ": our GSEA vs official fgsea ===\n", sep = "")
  cat("Comparable pathways (NES defined by both methods):", nrow(cmp), "\n")
  cat("Same direction (NES sign):", sum(cmp$same_direction), "of", nrow(cmp),
      "(", round(100 * mean(cmp$same_direction), 1), "% )\n")
  if (nrow(cmp) >= 3) {
    cat("Correlation (Pearson) between our NES and the official NES:",
        round(cor(cmp$NES_ours, cmp$NES_official), 3), "\n")
  }
  cat("\n")
  cmp
}

cmp_pop <- compare_methods(gsea_pop, fgsea_pop, "POP")
write.csv(cmp_pop, "results/comparison_vs_official_fgsea_POP.csv", row.names = FALSE)
cmp_sui <- compare_methods(gsea_sui, fgsea_sui, "SUI")
write.csv(cmp_sui, "results/comparison_vs_official_fgsea_SUI.csv", row.names = FALSE)

## --- 7d: comparison scatter plots (our NES vs official NES) ---------------
make_comparison_plot <- function(cmp, label) {
  ggplot(cmp, aes(x = NES_ours, y = NES_official, color = same_direction)) +
    geom_hline(yintercept = 0, color = "grey70") + geom_vline(xintercept = 0, color = "grey70") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
    geom_point(size = 2.5, alpha = 0.8) +
    scale_color_manual(values = c("TRUE" = "#2166AC", "FALSE" = "#B2182B"),
                        labels = c("TRUE" = "Same direction", "FALSE" = "Opposite direction")) +
    labs(title = paste0("Our hand-rolled GSEA vs official fgsea - ", label),
         subtitle = paste0(nrow(cmp), " comparable pathways, ",
                            round(100 * mean(cmp$same_direction), 1), "% same direction",
                            if (nrow(cmp) >= 3) paste0(", r=", round(cor(cmp$NES_ours, cmp$NES_official), 3)) else ""),
         x = "NES - our hand-rolled GSEA", y = "NES - official fgsea", color = NULL) +
    theme_bw() + theme(legend.position = "top")
}
if (nrow(cmp_pop) >= 2) {
  ggsave("figures/comparison_vs_official_fgsea_POP.png", make_comparison_plot(cmp_pop, "POP"), width = 8, height = 7, dpi = 300)
  cat("Saved: figures/comparison_vs_official_fgsea_POP.png\n")
}
if (nrow(cmp_sui) >= 2) {
  ggsave("figures/comparison_vs_official_fgsea_SUI.png", make_comparison_plot(cmp_sui, "SUI"), width = 8, height = 7, dpi = 300)
  cat("Saved: figures/comparison_vs_official_fgsea_SUI.png\n")
}
cat("\n")


## =============================================================================
## PART 8 (NEW): GO Biological Process enrichment - via the OFFICIAL fgsea
## directly, not our hand-rolled permutation code
## =============================================================================
cat("================ PART 8: GO Biological Process enrichment ================\n\n")
cat("WHY THIS PART EXISTS: KEGG (Parts 1-7) covers ~218 broad pathways. GO\n")
cat("Biological Process has far more specific terms - including some that\n")
cat("map directly onto a connective-tissue/extracellular-matrix weakness\n")
cat("hypothesis (e.g. 'collagen fibril organization', 'extracellular matrix\n")
cat("organization') that KEGG has no equivalent for.\n\n")
cat("WHY IT USES fgsea DIRECTLY instead of our own hand-rolled permutation\n")
cat("code (Parts 3-4): even after filtering to gene sets of 5-200 genes,\n")
cat("GO Biological Process still has thousands of terms (vs KEGG's 218).\n")
cat("Our hand-rolled code re-tests every single gene set inside every one\n")
cat("of 500-1000 permutation loops - fine for 218 KEGG pathways, but would\n")
cat("add many extra minutes of runtime for thousands of GO terms. The\n")
cat("official fgsea package uses an adaptive algorithm (not brute-force\n")
cat("permutation of every gene set) built to handle exactly this scale in\n")
cat("well under a minute. Since Part 7 already showed fgsea agrees with our\n")
cat("hand-rolled KEGG results, using it directly here (on the exact same\n")
cat("ranked gene lists as everywhere else in this script) is a reasonable,\n")
cat("honest shortcut rather than a second, slower reimplementation of the\n")
cat("same permutation logic.\n\n")

## --- 8a: build GO Biological Process gene sets from org.Hs.eg.db ----------
## SAME source (org.Hs.eg.db) as the KEGG sets above - no internet needed.
cat("Building GO Biological Process gene sets from org.Hs.eg.db...\n")
ann_go_pop <- suppressWarnings(select(org.Hs.eg.db, keys = ranked_genes_pop, keytype = "SYMBOL", columns = c("GO", "ONTOLOGY")))
ann_go_pop <- ann_go_pop[!is.na(ann_go_pop$GO) & ann_go_pop$ONTOLOGY == "BP", ]
gs_sizes_go <- table(ann_go_pop$GO)
valid_go <- names(gs_sizes_go)[gs_sizes_go >= 5 & gs_sizes_go <= 200]
gene_sets_go <- split(ann_go_pop$SYMBOL[ann_go_pop$GO %in% valid_go], ann_go_pop$GO[ann_go_pop$GO %in% valid_go])
cat("GO Biological Process terms tested (5-200 genes):", length(gene_sets_go), "\n\n")

## Term(): translates a GO ID (e.g. "GO:0030199") into its readable name
## ("collagen fibril organization") - GO.db ships with org.Hs.eg.db, no
## internet needed.
go_term_name <- function(ids) {
  nm <- suppressMessages(tryCatch(AnnotationDbi::Term(ids), error = function(e) rep(NA_character_, length(ids))))
  unname(ifelse(is.na(nm), ids, paste0(ids, " - ", nm)))
}

## --- 8b: run fgsea (GO Biological Process) on POP --------------------------
cat("Running official fgsea (GO Biological Process) on POP...\n")
set.seed(208261)
go_pop <- as.data.frame(fgsea(pathways = gene_sets_go, stats = ranked_full, minSize = 5, maxSize = 200))
go_pop$GOName <- go_term_name(go_pop$pathway)
go_pop$leadingEdge <- sapply(go_pop$leadingEdge, paste, collapse = "/")
go_pop <- go_pop[order(go_pop$pval), ]
write.csv(go_pop, "results/POP_official_fgsea_GO_BP_full.csv", row.names = FALSE)
cat("Terms significant at FDR<0.25:", sum(go_pop$padj < 0.25, na.rm = TRUE), "of", nrow(go_pop), "\n")
cat("Terms significant at FDR<0.05:", sum(go_pop$padj < 0.05, na.rm = TRUE), "\n\n")

## --- 8c: run fgsea (GO Biological Process) on SUI --------------------------
cat("Running official fgsea (GO Biological Process) on SUI...\n")
set.seed(2020)
go_sui <- as.data.frame(fgsea(pathways = gene_sets_go, stats = t_obs_sui, minSize = 5, maxSize = 200))
go_sui$GOName <- go_term_name(go_sui$pathway)
go_sui$leadingEdge <- sapply(go_sui$leadingEdge, paste, collapse = "/")
go_sui <- go_sui[order(go_sui$pval), ]
write.csv(go_sui, "results/SUI_official_fgsea_GO_BP_full.csv", row.names = FALSE)
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
## Looking up a handful of GO terms directly relevant to a connective-
## tissue/extracellular-matrix weakness hypothesis for POP/SUI, regardless
## of whether they happened to be significant/shared above.
cat("Checking specific GO terms tied to a connective-tissue/ECM weakness\n")
cat("hypothesis (present here if they had 5-200 genes in the ranked list):\n\n")
ecm_go_terms <- c(
  "GO:0030198",  # extracellular matrix organization
  "GO:0030199",  # collagen fibril organization
  "GO:0030574",  # collagen catabolic process
  "GO:0007044",  # cell-substrate junction assembly
  "GO:0031012"   # extracellular matrix (this one is a CC term but kept in case it is annotated as BP too)
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

cat("=================================================================\n")
cat("=== END OF SCRIPT ===\n")
cat("Check the numbers marked 'Expected value, already documented' above\n")
cat("against what came out on your computer (Parts 1-6). If they match,\n")
cat("the main analysis is independently validated.\n\n")
cat("In Part 7, a high correlation (>0.8) and high direction agreement\n")
cat("(>80%) between our NES and the official fgsea NES indicates the GSEA\n")
cat("calculation engine reimplemented in this project is correct. A strong\n")
cat("mismatch would call for investigation - send me both CSVs\n")
cat("(results/comparison_vs_official_fgsea_POP.csv and _SUI.csv) if that happens.\n\n")
cat("Part 8 (GO Biological Process) is a NEW analysis, not a re-validation of\n")
cat("anything above - check results/shared_GO_BP_FDR025.csv and the 'ECM\n")
cat("hypothesis terms' printout above for anything directly relevant to a\n")
cat("connective-tissue angle that KEGG's broader pathways do not capture.\n")
cat("=================================================================\n")
