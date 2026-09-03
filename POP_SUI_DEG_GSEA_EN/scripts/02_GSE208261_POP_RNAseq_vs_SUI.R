# =============================================================================
# GSE208261 (POP, RNA-seq) vs SUI (Wei 2020): DEG, GSEA, shared pathways
# =============================================================================
# Extends POP_SUI_analysis.R with a THIRD, independent POP dataset - this
# time RNA-seq instead of microarray, from a different set of patients and a
# different measurement technology entirely. This script demonstrates the
# statistically valid way to combine microarray and RNA-seq evidence: each
# dataset is analyzed ENTIRELY ON ITS OWN, with the statistical method
# appropriate to its own technology (DESeq2 for RNA-seq counts here, with
# edgeR+voom+limma run alongside for comparison and reused for GSEA ranking;
# limma directly on log-intensities for the microarray SUI data below, same
# as the main script). Raw expression values from the two technologies are
# NEVER pooled or merged into one matrix - only the DERIVED, per-dataset
# results (which genes/pathways are significant, and in which direction) are
# compared at the end. This is the same architecture already used across
# every dataset in this project (3 different microarray platforms were never
# merged either) - RNA-seq is just one more independently-analyzed source.
#
# DATA-QUALITY NOTE found while preparing this script: the sample naming in
# GSE208261 (Control_D1-6/Y1-6, POP_D1-6/Y1-6) suggests a tissue split within
# BOTH groups, but the actual GEO metadata (data/GSE208261_sample_metadata.csv,
# extracted from the submitters' own family.soft file) shows this is true
# ONLY for the controls (D = uterosacral ligaments, Y = anterior vaginal
# wall) - ALL 12 POP samples (both "_D" and "_Y") are anterior vaginal wall;
# there are no POP-side ligament samples at all, and the real meaning of the
# "_D"/"_Y" suffix within the POP arm is not stated anywhere in the
# metadata. Consequence: the only tissue-matched, valid comparison here is
# 6 controls vs 12 POP, BOTH anterior vaginal wall (the same tissue as
# GSE53868, and the closest match of any POP dataset in this project to
# Wei2020's periurethral vaginal wall). The 6 ligament-control samples have
# no POP counterpart and are excluded below.
# =============================================================================

suppressMessages({
  library(limma)
  library(edgeR)
  library(DESeq2)
  library(org.Hs.eg.db)
  library(ggplot2)
  library(pheatmap)
})
if (!requireNamespace("readxl", quietly = TRUE)) install.packages("readxl")
if (!requireNamespace("ggrepel", quietly = TRUE)) install.packages("ggrepel")
library(readxl)
library(ggrepel)

dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

# Full KEGG pathway ID -> name lookup, loaded from the shared, audited table
# data/kegg_pathway_names.csv (223 IDs covering every pathway used anywhere
# in this project) instead of a per-script partial dict - one source of
# truth so every figure across every script names the same ID the same way.
kegg_tab <- read.csv("data/kegg_pathway_names.csv", colClasses = c("character", "character"))
kegg_names <- setNames(kegg_tab$Name, kegg_tab$PATH5)
kegg_label <- function(id) {
  nm <- kegg_names[id]
  ifelse(is.na(nm), paste0("KEGG ", id), paste0("KEGG ", id, " - ", nm))
}
calc_es <- function(hit_idx, scores_abs, N) {
  Nh <- length(hit_idx); Nm <- N - Nh
  if (Nm <= 0 || Nh == 0) return(NA)
  sum_hit <- sum(scores_abs[hit_idx])
  step <- rep(-1 / Nm, N); step[hit_idx] <- scores_abs[hit_idx] / sum_hit
  running <- cumsum(step); running[which.max(abs(running))]
}


## =============================================================================
## STEP 1: DEG for POP (GSE208261, RNA-seq) - tissue-matched subset only
## =============================================================================
cat("\n================ STEP 1: DEG for POP (GSE208261, RNA-seq) ================\n\n")

meta <- read.csv("data/GSE208261_sample_metadata.csv")
meta_use <- meta[meta$Tissue == "anterior vaginal wall", ]
cat("Samples used (anterior vaginal wall only):", nrow(meta_use), "(",
    sum(meta_use$Treatment == "Control"), "Control /",
    sum(meta_use$Treatment == "POP"), "POP )\n")
cat("Excluded (uterosacral ligament, no POP counterpart in this series):",
    sum(meta$Tissue != "anterior vaginal wall"), "control samples\n\n")

counts_raw <- read.delim("data/GSE208261_raw_counts.tsv", row.names = 1, check.names = FALSE)
counts <- as.matrix(counts_raw[, meta_use$GSM])
stopifnot(identical(colnames(counts), meta_use$GSM))
cat("Raw count matrix:", nrow(counts), "Entrez genes x", ncol(counts), "samples\n\n")

# Entrez GeneID -> gene symbol, collapse any duplicate symbols by summing counts
ann_id <- suppressWarnings(select(org.Hs.eg.db, keys = rownames(counts), keytype = "ENTREZID", columns = "SYMBOL"))
ann_id <- ann_id[!is.na(ann_id$SYMBOL) & !duplicated(ann_id$ENTREZID), ]
counts <- counts[rownames(counts) %in% ann_id$ENTREZID, ]
rownames(counts) <- ann_id$SYMBOL[match(rownames(counts), ann_id$ENTREZID)]
counts <- rowsum(counts, group = rownames(counts))
cat("After Entrez->symbol mapping and collapsing duplicates:", nrow(counts), "genes\n\n")

group_pop2 <- factor(meta_use$Treatment, levels = c("Control", "POP"))
design_pop2 <- model.matrix(~group_pop2)

cat("Two RNA-seq DEG methods are run and reported side by side, on the exact\n")
cat("same 6-vs-12 tissue-matched samples: DESeq2 (primary - negative-binomial\n")
cat("GLM, its own dispersion shrinkage and independent filtering) and\n")
cat("edgeR+voom+limma (secondary - used again below to rank genes for GSEA).\n")
cat("They can legitimately disagree at the margin: DESeq2's independent\n")
cat("filtering removes genes that stand no chance of significance BEFORE\n")
cat("multiple-testing correction, shrinking the correction's denominator and\n")
cat("increasing power exactly in borderline cases like this dataset. Neither\n")
cat("is 'wrong' - reporting both, rather than only the one that finds\n")
cat("something, is the transparent way to handle a real, known point of\n")
cat("disagreement between two standard, widely used tools.\n\n")

## --- DESeq2 (primary DEG method) --------------------------------------------
counts_int <- counts; storage.mode(counts_int) <- "integer"
coldata_pop2 <- data.frame(row.names = meta_use$GSM, group = group_pop2)
dds <- DESeqDataSetFromMatrix(countData = counts_int[, meta_use$GSM], colData = coldata_pop2, design = ~group)
dds <- dds[rowSums(counts(dds) >= 10) >= 6, ]
dds <- DESeq(dds, quiet = TRUE)
res_deseq2 <- as.data.frame(results(dds, contrast = c("group", "POP", "Control"), alpha = 0.05))
res_deseq2$Gene <- rownames(res_deseq2)
res_deseq2 <- res_deseq2[order(res_deseq2$pvalue), ]

pop2_full <- res_deseq2[, c("Gene", "log2FoldChange", "baseMean", "stat", "pvalue", "padj")]
colnames(pop2_full) <- c("Gene", "logFC", "baseMean", "stat", "P.Value", "adj.P.Val")
write.csv(pop2_full, "results/06_POP_GSE208261_DESeq2_full.csv", row.names = FALSE)

pop2_deg <- subset(pop2_full, !is.na(adj.P.Val) & adj.P.Val < 0.05 & abs(logFC) > 1)
pop2_deg <- pop2_deg[order(pop2_deg$adj.P.Val), ]
write.csv(pop2_deg, "results/06_POP_GSE208261_DEG_logFC1_FDR05.csv", row.names = FALSE)

cat("=== DESeq2 (primary DEG method) ===\n")
cat("Genes kept (count>=10 in >=6 samples):", nrow(dds), "of", nrow(counts), "\n")
cat("DEG (|log2FC|>1, FDR<0.05):", nrow(pop2_deg), "(", sum(pop2_deg$logFC > 0), "up /",
    sum(pop2_deg$logFC < 0), "down )\n")
cat("Minimum FDR reached:", signif(min(pop2_full$adj.P.Val, na.rm = TRUE), 3), "\n\n")
print(head(pop2_full, 10))
cat("\n")

## --- edgeR + voom + limma (secondary - reused below for GSEA ranking) ------
dge <- DGEList(counts = counts, group = group_pop2)
keep_expr <- filterByExpr(dge, design_pop2)
dge <- dge[keep_expr, , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge, method = "TMM")
cat("Library sizes range from", round(min(dge$samples$lib.size) / 1e6, 1), "to",
    round(max(dge$samples$lib.size) / 1e6, 1), "million reads (uneven depth -",
    "the 6 'POP_D' samples run notably shallower than the rest; TMM",
    "normalization corrects for this before testing).\n\n")

voom_fit <- voom(dge, design_pop2)
fit_pop2 <- eBayes(lmFit(voom_fit, design_pop2))
pop2_full_limma <- topTable(fit_pop2, coef = "group_pop2POP", number = Inf, sort.by = "P")
pop2_full_limma$Gene <- rownames(pop2_full_limma)
pop2_full_limma <- pop2_full_limma[, c("Gene", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val")]
write.csv(pop2_full_limma, "results/06_POP_GSE208261_voom_limma_full.csv", row.names = FALSE)
pop2_deg_limma <- subset(pop2_full_limma, adj.P.Val < 0.05 & abs(logFC) > 1)

cat("=== edgeR+voom+limma (secondary - only used to rank genes for GSEA below) ===\n")
cat("Genes kept after filterByExpr:", nrow(dge), "of", nrow(counts), "\n")
cat("DEG (|log2FC|>1, FDR<0.05):", nrow(pop2_deg_limma), "\n")
cat("Minimum FDR reached:", signif(min(pop2_full_limma$adj.P.Val), 3), "\n\n")

cat("SUMMARY: DESeq2 finds", nrow(pop2_deg), "DEG(s) where voom+limma finds",
    nrow(pop2_deg_limma), "- a real, modest difference at the margin between\n")
cat("two standard tools on a dataset with weak signal, not a contradiction or\n")
cat("a bug in either. The DESeq2 result is used as the primary DEG call below\n")
cat("and throughout the rest of this script and script 03.\n\n")


## =============================================================================
## STEP 2: GSEA for POP (GSE208261) - CLASSIC (phenotype/group permutation)
## =============================================================================
cat("================ STEP 2: GSEA for POP (GSE208261) - CLASSIC method ================\n\n")
cat("Method: classic GSEA, phenotype permutation (shuffling the 18 group\n")
cat("labels, unpaired). With 6 vs 12 samples there are choose(18,6)=18,564\n")
cat("possible relabelings - ample resolution, unlike the 3-vs-3 SUI design.\n")
cat("Ranking = moderated t-statistic from Step 1's secondary edgeR+voom+limma\n")
cat("fit (not DESeq2 - permuting DESeq2's full pipeline 500 times is far too\n")
cat("slow; voom+limma's precomputed-weights refit is the standard fast\n")
cat("approach for permutation testing, and GSEA needs a sensible relative\n")
cat("ranking, not a calibrated single-gene significance call).\n")
cat("For speed, the voom precision weights are computed ONCE on the true\n")
cat("design and then reused for every permutation (only the linear model fit\n")
cat("is repeated, not the mean-variance trend) - standard practice for\n")
cat("permutation testing on voom-transformed data.\n\n")

set.seed(208261)
ranked_full <- fit_pop2$t[, "group_pop2POP"]
ord2 <- order(-ranked_full)
ranked_genes_pop2 <- names(ranked_full)[ord2]
ranked_scores_pop2 <- ranked_full[ord2]
N_pop2 <- length(ranked_genes_pop2)

ann_pop2 <- suppressWarnings(select(org.Hs.eg.db, keys = ranked_genes_pop2, keytype = "SYMBOL", columns = "PATH"))
ann_pop2 <- ann_pop2[!is.na(ann_pop2$PATH), ]
gs_sizes_pop2 <- table(ann_pop2$PATH)
valid_paths_pop2 <- names(gs_sizes_pop2)[gs_sizes_pop2 >= 5 & gs_sizes_pop2 <= 200]
gene_sets_pop2 <- split(ann_pop2$SYMBOL[ann_pop2$PATH %in% valid_paths_pop2], ann_pop2$PATH[ann_pop2$PATH %in% valid_paths_pop2])
cat("KEGG pathways tested (5-200 members):", length(gene_sets_pop2), "\n")

hit_idx_pop2 <- lapply(gene_sets_pop2, function(g) which(ranked_genes_pop2 %in% g))
hit_idx_pop2 <- hit_idx_pop2[sapply(hit_idx_pop2, length) >= 3]
cat("Pathways with >=3 genes in the ranked list:", length(hit_idx_pop2), "\n")
es_obs_pop2 <- sapply(hit_idx_pop2, calc_es, scores_abs = abs(ranked_scores_pop2), N = N_pop2)

n_perm3 <- 500
cat("Running", n_perm3, "phenotype (group-label) permutations...\n")
t0 <- Sys.time()
gene_sets_syms_pop2 <- lapply(hit_idx_pop2, function(idx) ranked_genes_pop2[idx])
perm_es_pop2 <- matrix(NA_real_, nrow = n_perm3, ncol = length(hit_idx_pop2))
for (i in seq_len(n_perm3)) {
  perm_group <- sample(group_pop2)
  perm_design <- model.matrix(~perm_group)
  perm_fit <- eBayes(lmFit(voom_fit, perm_design))
  t_perm <- perm_fit$t[, 2]
  rank_of_gene <- rank(-t_perm, ties.method = "first")
  scores_abs_sorted <- sort(abs(t_perm), decreasing = TRUE)
  for (j in seq_along(gene_sets_syms_pop2)) {
    hidx <- rank_of_gene[gene_sets_syms_pop2[[j]]]
    if (length(hidx) >= 3) perm_es_pop2[i, j] <- calc_es(hidx, scores_abs_sorted, N_pop2)
  }
}
cat("Done in", round(difftime(Sys.time(), t0, units = "secs"), 1), "seconds\n\n")

pval_pop2 <- numeric(length(hit_idx_pop2)); nes_pop2 <- numeric(length(hit_idx_pop2))
for (j in seq_along(hit_idx_pop2)) {
  pe <- perm_es_pop2[, j]; pe <- pe[!is.na(pe)]
  if (es_obs_pop2[j] >= 0) {
    pval_pop2[j] <- (sum(pe >= es_obs_pop2[j]) + 1) / (length(pe) + 1)
    base <- mean(pe[pe >= 0]); if (is.nan(base) || base == 0) base <- NA
  } else {
    pval_pop2[j] <- (sum(pe <= es_obs_pop2[j]) + 1) / (length(pe) + 1)
    base <- mean(abs(pe[pe < 0])); if (is.nan(base) || base == 0) base <- NA
  }
  nes_pop2[j] <- es_obs_pop2[j] / base
}

leading_edge_pop2 <- sapply(hit_idx_pop2, function(idx) paste(ranked_genes_pop2[idx], collapse = "/"))
gsea_pop2 <- data.frame(PATH = names(hit_idx_pop2), Nh = sapply(hit_idx_pop2, length),
                         ES = es_obs_pop2, NES = nes_pop2, pvalue = pval_pop2, leadingEdge = leading_edge_pop2)
gsea_pop2$p.adjust <- p.adjust(gsea_pop2$pvalue, "BH")
gsea_pop2$PathwayName <- kegg_label(gsea_pop2$PATH)
gsea_pop2 <- gsea_pop2[order(gsea_pop2$pvalue), ]
write.csv(gsea_pop2, "results/07_GSEA_classic_POP_GSE208261_KEGG.csv", row.names = FALSE)

cat("=== GSEA classic, POP GSE208261 (KEGG) ===\n")
cat("Pathways significant at FDR<0.25:", sum(gsea_pop2$p.adjust < 0.25, na.rm = TRUE), "of", nrow(gsea_pop2), "\n")
cat("Pathways significant at FDR<0.05:", sum(gsea_pop2$p.adjust < 0.05, na.rm = TRUE), "\n\n")
print(head(gsea_pop2[, c("PathwayName","Nh","NES","pvalue","p.adjust")], 10))
cat("\n")


## =============================================================================
## STEP 3: DEG and GSEA for SUI (Wei 2020) - same method as the main script
## =============================================================================
cat("================ STEP 3: DEG and GSEA for SUI (Wei 2020) ================\n\n")
cat("Re-derived here (not just re-loaded) so this script is fully standalone -\n")
cat("identical method to POP_SUI_analysis.R Steps 2 and 4; see that script for\n")
cat("the full methodological rationale.\n\n")

wei_xls <- "data/Wei2020_TableS2_mRNA.xls"
up_sheet   <- read_excel(wei_xls, sheet = "up_Sui_vs_Ctrl",   skip = 17)
down_sheet <- read_excel(wei_xls, sheet = "down_Sui_vs_Ctrl", skip = 17)
up_sheet$Direction <- "up"; down_sheet$Direction <- "down"
sample_cols <- c("[Sui1, Sui](normalized)", "[Sui2, Sui](normalized)", "[Sui3, Sui](normalized)",
                  "[Ctrl1, Ctrl](normalized)", "[Ctrl2, Ctrl](normalized)", "[Ctrl3, Ctrl](normalized)")
keep_cols <- c("GeneSymbol", "P-value", "FDR", "Fold Change", "Direction", sample_cols)
wei_both <- rbind(as.data.frame(up_sheet[, keep_cols]), as.data.frame(down_sheet[, keep_cols]))
colnames(wei_both) <- c("GeneSymbol", "PValue", "FDR", "FoldChange", "Direction",
                         "Sui1", "Sui2", "Sui3", "Ctrl1", "Ctrl2", "Ctrl3")
wei_both <- wei_both[!is.na(wei_both$GeneSymbol), ]
wei_both$logFC <- ifelse(wei_both$Direction == "down",
                          -log2(wei_both$FoldChange), log2(wei_both$FoldChange))
wei_both <- wei_both[order(wei_both$PValue), ]
sui_full <- wei_both[!duplicated(wei_both$GeneSymbol), ]
rownames(sui_full) <- NULL
sui_deg <- sui_full
cat("DEG SUI (authors' original criterion):", nrow(sui_deg), "(",
    sum(sui_deg$Direction == "up"), "up /", sum(sui_deg$Direction == "down"), "down )\n\n")

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

ann_sui <- suppressWarnings(select(org.Hs.eg.db, keys = ranked_genes_sui, keytype = "SYMBOL", columns = "PATH"))
ann_sui <- ann_sui[!is.na(ann_sui$PATH), ]
gs_sizes_sui <- table(ann_sui$PATH)
valid_paths_sui <- names(gs_sizes_sui)[gs_sizes_sui >= 5 & gs_sizes_sui <= 200]
gene_sets_sui <- split(ann_sui$SYMBOL[ann_sui$PATH %in% valid_paths_sui], ann_sui$PATH[ann_sui$PATH %in% valid_paths_sui])
hit_idx_sui <- lapply(gene_sets_sui, function(g) which(ranked_genes_sui %in% g))
hit_idx_sui <- hit_idx_sui[sapply(hit_idx_sui, length) >= 3]
es_obs_sui <- sapply(hit_idx_sui, calc_es, scores_abs = abs(ranked_scores_sui), N = N_sui)

n_perm4 <- 1000
cat("Running", n_perm4, "gene-set-label permutations for SUI (preranked GSEA)...\n")
t0 <- Sys.time()
scores_abs_sui <- abs(ranked_scores_sui)
perm_es_sui <- matrix(NA_real_, nrow = n_perm4, ncol = length(hit_idx_sui))
for (i in seq_len(n_perm4)) {
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

cat("=== GSEA preranked, SUI (KEGG) ===\n")
cat("Pathways significant at FDR<0.25:", sum(gsea_sui$p.adjust < 0.25, na.rm = TRUE), "of", nrow(gsea_sui), "\n")
cat("Pathways significant at FDR<0.05:", sum(gsea_sui$p.adjust < 0.05, na.rm = TRUE), "\n\n")
cat("\n")


## =============================================================================
## STEP 4: Shared pathways (POP GSE208261 x SUI) + gene-level direction table
## =============================================================================
cat("================ STEP 4: Shared pathways GSE208261(POP) x SUI ================\n\n")

report_shared2 <- function(fdr_cut) {
  pop_sig <- subset(gsea_pop2, p.adjust < fdr_cut)
  sui_sig <- subset(gsea_sui, p.adjust < fdr_cut)
  shared_ids <- intersect(pop_sig$PATH, sui_sig$PATH)
  cat("--- FDR <", fdr_cut, ": POP(GSE208261) significant =", nrow(pop_sig),
      "| SUI significant =", nrow(sui_sig), "| SHARED =", length(shared_ids), "---\n")
  if (length(shared_ids) == 0) return(data.frame())
  out <- merge(pop_sig[pop_sig$PATH %in% shared_ids, c("PATH","PathwayName","Nh","NES","p.adjust")],
               sui_sig[sui_sig$PATH %in% shared_ids, c("PATH","Nh","NES","p.adjust")],
               by = "PATH", suffixes = c("_POP", "_SUI"))
  out$Same_direction <- sign(out$NES_POP) == sign(out$NES_SUI)
  out <- out[order(out$p.adjust_POP), ]
  print(out[, c("PathwayName","Nh_POP","NES_POP","p.adjust_POP","Nh_SUI","NES_SUI","p.adjust_SUI","Same_direction")])
  n_na <- sum(is.na(out$Same_direction))
  cat("Same direction in both diseases:", sum(out$Same_direction, na.rm = TRUE), "of",
      sum(!is.na(out$Same_direction)), "pathways with a comparable NES on both sides")
  if (n_na > 0) cat(" (", n_na, "excluded - NES undefined on at least one side)")
  cat("\n\n")
  out
}

shared2_025 <- report_shared2(0.25)
shared2_005 <- report_shared2(0.05)
if (nrow(shared2_025) > 0) write.csv(shared2_025, "results/08_shared_pathways_GSE208261xSUI_FDR025.csv", row.names = FALSE)
if (nrow(shared2_005) > 0) write.csv(shared2_005, "results/08_shared_pathways_GSE208261xSUI_FDR005.csv", row.names = FALSE)

shared2_for_table <- if (nrow(shared2_005) > 0) shared2_005 else shared2_025
cat("\n--- Gene-level direction table for the shared pathways ---\n")
if (nrow(shared2_for_table) == 0) {
  cat("No shared pathway at either threshold - no gene-level table to build.\n\n")
  gene_dir_table2 <- data.frame()
} else {
  pop_lookup2 <- setNames(pop2_full$logFC, pop2_full$Gene)
  sui_lookup2 <- setNames(sui_full$logFC, sui_full$GeneSymbol)
  pop2_tested <- ranked_genes_pop2
  sui_tested2 <- ranked_genes_sui

  rows2 <- lapply(seq_len(nrow(shared2_for_table)), function(i) {
    pid <- shared2_for_table$PATH[i]
    pname <- shared2_for_table$PathwayName[i]
    full_members <- suppressWarnings(select(org.Hs.eg.db, keys = pid, keytype = "PATH", columns = "SYMBOL")$SYMBOL)
    members <- intersect(unique(full_members), union(pop2_tested, sui_tested2))
    data.frame(
      Pathway = pname, KEGG_ID = pid, Gene = members,
      logFC_POP = unname(pop_lookup2[members]), logFC_SUI = unname(sui_lookup2[members]),
      stringsAsFactors = FALSE)
  })
  gene_dir_table2 <- do.call(rbind, rows2)
  gene_dir_table2$Direction_POP <- ifelse(is.na(gene_dir_table2$logFC_POP), "not tested",
                                           ifelse(gene_dir_table2$logFC_POP > 0, "up", "down"))
  gene_dir_table2$Direction_SUI <- ifelse(is.na(gene_dir_table2$logFC_SUI), "not tested",
                                           ifelse(gene_dir_table2$logFC_SUI > 0, "up", "down"))
  gene_dir_table2$Concordant <- with(gene_dir_table2,
                                      Direction_POP != "not tested" & Direction_SUI != "not tested" &
                                        Direction_POP == Direction_SUI)
  gene_dir_table2 <- gene_dir_table2[order(gene_dir_table2$Pathway, -gene_dir_table2$Concordant), ]
  write.csv(gene_dir_table2, "results/08_shared_pathways_gene_direction_table.csv", row.names = FALSE)

  both_tested2 <- subset(gene_dir_table2, Direction_POP != "not tested" & Direction_SUI != "not tested")
  cat("Genes in shared pathways tested in both diseases:", nrow(both_tested2),
      "| concordant direction:", sum(both_tested2$Concordant),
      "(", round(100 * mean(both_tested2$Concordant), 1), "% )\n\n")
  print(head(gene_dir_table2, 20))
}
cat("\n")


## =============================================================================
## STEP 5: GRAPHICS
## =============================================================================
cat("================ STEP 5: Graphics ================\n\n")

make_volcano <- function(df, logfc_col, p_col, label_col, title, fc_cut = 1, p_cut = 0.05) {
  df <- df[!is.na(df[[logfc_col]]) & !is.na(df[[p_col]]), ]
  df$negLog10P <- -log10(df[[p_col]])
  df$sig <- "NS"
  df$sig[df[[logfc_col]] > fc_cut & df[[p_col]] < p_cut] <- "Up"
  df$sig[df[[logfc_col]] < -fc_cut & df[[p_col]] < p_cut] <- "Down"
  df$sig <- factor(df$sig, levels = c("Down", "NS", "Up"))
  top_lab <- df[order(df[[p_col]]), ][1:15, ]
  ggplot(df, aes(x = .data[[logfc_col]], y = negLog10P, color = sig)) +
    geom_point(alpha = 0.6, size = 1.3) +
    scale_color_manual(values = c(Down = "#2166AC", NS = "grey75", Up = "#B2182B")) +
    geom_vline(xintercept = c(-fc_cut, fc_cut), linetype = "dashed", color = "grey40") +
    geom_hline(yintercept = -log10(p_cut), linetype = "dashed", color = "grey40") +
    geom_text_repel(data = top_lab, aes(label = .data[[label_col]]), size = 3, color = "black",
                     max.overlaps = 20, segment.size = 0.2) +
    labs(title = title, x = "log2(Fold Change)",
         y = expression(-log[10](italic(p)~value)), color = NULL) +
    theme_bw() + theme(legend.position = "top", plot.title = element_text(face = "bold"))
}

pop2_full$Gene_label <- pop2_full$Gene
deg_note <- if (nrow(pop2_deg) == 0) "none reach FDR<0.05" else
  paste0(nrow(pop2_deg), " gene(s) reach FDR<0.05")
deg_note_full <- if (nrow(pop2_deg) == 0) "none reach FDR<0.05" else
  paste0(nrow(pop2_deg), " gene(s) reach FDR<0.05: ", paste(pop2_deg$Gene, collapse = ", "))
cat("Volcano plot DEG note:", deg_note_full, "\n\n")
v_pop2 <- make_volcano(pop2_full, "logFC", "P.Value", "Gene_label",
                        paste0("Volcano plot - POP (GSE208261, RNA-seq, DESeq2)\n",
                               "POP vs control, anterior vaginal wall\n",
                               "(points shown at raw p<0.05, |log2FC|>1; ", deg_note, ")"))
ggsave("figures/10_volcano_POP_GSE208261.png", v_pop2, width = 9.5, height = 6.5, dpi = 300)
cat("Saved: figures/10_volcano_POP_GSE208261.png\n\n")

# Heatmap of the top genes by p-value (log-CPM from the voom pipeline, for a
# reliable expression-visualization scale even though gene selection/ranking
# above comes from DESeq2).
top_pop2_genes <- head(pop2_full$Gene, 40)
logcpm2 <- voom_fit$E[rownames(voom_fit$E) %in% top_pop2_genes, , drop = FALSE]
ann_col_pop2 <- data.frame(Group = meta_use$Treatment, row.names = meta_use$GSM)
ord_pop2 <- order(meta_use$Treatment)
png("figures/11_heatmap_top_POP_GSE208261.png", width = 2600, height = 3200, res = 300)
pheatmap(logcpm2[, ord_pop2], scale = "row",
         annotation_col = ann_col_pop2[ord_pop2, , drop = FALSE],
         cluster_cols = FALSE, cluster_rows = TRUE,
         main = paste0("Top 40 genes by p-value (DESeq2) - POP (GSE208261, RNA-seq)\n(row z-score, log-CPM; ", deg_note, ")"),
         fontsize = 7, show_colnames = TRUE, fontsize_col = 6)
dev.off()
cat("Saved: figures/11_heatmap_top_POP_GSE208261.png\n\n")

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
    theme_bw() + theme(axis.text.y = element_text(size = 8))
}
p_bar_pop2 <- make_gsea_barplot(gsea_pop2, "GSEA (classic) - top KEGG pathways in POP (GSE208261, RNA-seq)")
ggsave("figures/12_GSEA_barplot_POP_GSE208261.png", p_bar_pop2, width = 10, height = 6, dpi = 300)
cat("Saved: figures/12_GSEA_barplot_POP_GSE208261.png\n\n")

compare_pool2 <- if (nrow(shared2_025) > 0) shared2_025 else data.frame()
if (nrow(compare_pool2) > 0) {
  p_shared2 <- ggplot(compare_pool2, aes(x = NES_POP, y = NES_SUI, label = PathwayName)) +
    geom_hline(yintercept = 0, color = "grey70") + geom_vline(xintercept = 0, color = "grey70") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
    geom_point(aes(color = Same_direction), size = 3) +
    scale_color_manual(values = c("TRUE" = "#2166AC", "FALSE" = "#B2182B"),
                        labels = c("TRUE" = "Same direction", "FALSE" = "Opposite direction")) +
    geom_text_repel(size = 3, max.overlaps = 30) +
    labs(title = "Shared pathways (FDR<0.25): NES in POP (GSE208261, RNA-seq) vs NES in SUI",
         x = "NES - POP GSE208261 (classic GSEA)", y = "NES - SUI (preranked GSEA)", color = NULL) +
    theme_bw() + theme(legend.position = "top")
  ggsave("figures/13_shared_pathways_NES_comparison_GSE208261.png", p_shared2, width = 9, height = 7, dpi = 300)
  cat("Saved: figures/13_shared_pathways_NES_comparison_GSE208261.png\n\n")
} else {
  cat("SKIPPED: figures/13_shared_pathways_NES_comparison_GSE208261.png (no shared pathway at FDR<0.25)\n\n")
}

if (nrow(gene_dir_table2) > 0) {
  both_tested2 <- subset(gene_dir_table2, Direction_POP != "not tested" & Direction_SUI != "not tested")
  if (nrow(both_tested2) >= 2) {
    fc_mat2 <- as.matrix(both_tested2[, c("logFC_POP", "logFC_SUI")])
    dir_mat2 <- ifelse(fc_mat2 > 0, 1, -1)
    rownames(dir_mat2) <- rownames(fc_mat2) <- make.unique(both_tested2$Gene)
    ord2b <- order(both_tested2$Pathway, -both_tested2$Concordant, -abs(rowSums(fc_mat2, na.rm = TRUE)))
    dir_mat2 <- dir_mat2[ord2b, , drop = FALSE]; fc_mat2 <- fc_mat2[ord2b, , drop = FALSE]
    dir_mat2 <- head(dir_mat2, 50); fc_mat2 <- head(fc_mat2, 50)
    ann_row2 <- data.frame(Pathway = both_tested2$Pathway[ord2b][seq_len(nrow(dir_mat2))],
                            row.names = rownames(dir_mat2))
    colnames(dir_mat2) <- colnames(fc_mat2) <- c("POP", "SUI")
    png("figures/14_gene_direction_heatmap_GSE208261xSUI.png",
        width = 2600, height = max(1800, 45 * nrow(dir_mat2)), res = 300)
    pheatmap(dir_mat2, cluster_cols = FALSE, cluster_rows = FALSE,
             annotation_row = ann_row2,
             color = c("#2166AC", "#B2182B"), breaks = c(-2, 0, 2),
             legend = FALSE,
             display_numbers = matrix(sprintf("%.2f", fc_mat2), nrow(fc_mat2)),
             number_color = "white", fontsize_number = 7,
             main = "Gene direction, shared pathways (POP GSE208261 vs SUI)\n(fill = up/down; number = log2 Fold Change)",
             fontsize = 7, fontsize_row = 6)
    dev.off()
    cat("Saved: figures/14_gene_direction_heatmap_GSE208261xSUI.png\n\n")
  } else {
    cat("SKIPPED: figures/14_gene_direction_heatmap_GSE208261xSUI.png (fewer than 2 genes tested in both)\n\n")
  }
} else {
  cat("SKIPPED: figures/14_gene_direction_heatmap_GSE208261xSUI.png (no shared-pathway gene table)\n\n")
}

cat("\n=== ALL GSE208261 FIGURES SAVED TO figures/ ===\n")
cat("10_volcano_POP_GSE208261.png\n11_heatmap_top_POP_GSE208261.png\n")
cat("12_GSEA_barplot_POP_GSE208261.png\n13_shared_pathways_NES_comparison_GSE208261.png\n")
cat("14_gene_direction_heatmap_GSE208261xSUI.png\n")
cat("\n=== DONE (script 02) ===\n")
