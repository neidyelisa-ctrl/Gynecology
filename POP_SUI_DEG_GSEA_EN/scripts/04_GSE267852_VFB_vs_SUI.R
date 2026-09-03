# =============================================================================
# GSE267852 (POP, RNA-seq, vaginal fibroblasts) vs SUI (Wei 2020):
# DEG, GSEA, shared pathways
# =============================================================================
# A FOURTH, independent POP dataset - Tchoukalova/Chen (Mayo Clinic), GSE267852,
# "Cell type specific differences in the transcriptomes of adipose derived
# stem cells and vaginal fibroblasts in patients with pelvic organ prolapse."
# RNA-seq of PRIMARY CULTURED CELLS (not tissue biopsy like GSE53868/GSE208261):
# adipose-derived stem cells (ASCs) and vaginal fibroblasts (VFBs), 6 POP vs
# 6 continent controls each, PAIRED cell type per subject. This script uses
# ONLY the 12 VFB samples (6 vs 6) - the cell type actually comparable to the
# other vaginal-tissue-derived datasets in this project; the 12 ASC samples
# (adipose-derived, a different tissue origin entirely) are not used here.
#
# SET EXPECTATIONS BEFORE RUNNING: the original authors' own abstract states
# they found "no differentially expressed genes (DEG) between POP and CTRL in
# ASCs and VFBs" using DESeq2 at FDR<0.05 - the same primary method and
# threshold used below. They only found a signal (23 up / 29 down) using a
# MUCH more lenient, uncorrected criterion (raw p<0.01, no FDR correction at
# all). If this script also finds few or no DEG at FDR<0.05, THAT MATCHES
# THE ORIGINAL PUBLISHED RESULT - it is not a pipeline failure, and relaxing
# the threshold to manufacture significance would contradict the source
# study's own stated finding, not just this project's convention.
#
# Method: DESeq2 (primary DEG call, matching the original study's own method)
# alongside edgeR+voom+limma (secondary, reused to rank genes for GSEA).
# GSEA here uses CLASSIC (phenotype/group permutation) - unlike GSE208261
# (6 vs 12) and unlike SUI (3 vs 3), this design is BALANCED 6 vs 6, giving
# choose(12,6) = 924 possible relabelings, ample resolution for a real
# phenotype-permutation test.
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

kegg_names <- c(
  "00190"="Oxidative phosphorylation","01100"="Metabolic pathways",
  "03010"="Ribosome","03013"="RNA transport","03015"="mRNA surveillance pathway",
  "03018"="RNA degradation","03020"="RNA polymerase","03022"="Basal transcription factors",
  "03030"="DNA replication","03040"="Spliceosome","03050"="Proteasome",
  "03060"="Protein export","03320"="PPAR signaling pathway",
  "04010"="MAPK signaling pathway","04012"="ErbB signaling pathway",
  "04020"="Calcium signaling pathway","04060"="Cytokine-cytokine receptor interaction",
  "04062"="Chemokine signaling pathway","04066"="HIF-1 signaling pathway",
  "04068"="FoxO signaling pathway","04110"="Cell cycle",
  "04115"="p53 signaling pathway","04120"="Ubiquitin mediated proteolysis",
  "04140"="Autophagy","04141"="Protein processing in endoplasmic reticulum",
  "04142"="Lysosome","04144"="Endocytosis","04145"="Phagosome",
  "04150"="mTOR signaling pathway","04151"="PI3K-Akt signaling pathway",
  "04210"="Apoptosis","04310"="Wnt signaling pathway",
  "04330"="Notch signaling pathway","04340"="Hedgehog signaling pathway",
  "04350"="TGF-beta signaling pathway","04370"="VEGF signaling pathway",
  "04510"="Focal adhesion","04512"="ECM-receptor interaction",
  "04520"="Adherens junction","04530"="Tight junction","04540"="Gap junction",
  "04610"="Complement and coagulation cascades","04620"="Toll-like receptor signaling pathway",
  "04621"="NOD-like receptor signaling pathway","04630"="JAK-STAT signaling pathway",
  "04660"="T cell receptor signaling pathway","04662"="B cell receptor signaling pathway",
  "04664"="Fc epsilon RI signaling pathway","04670"="Leukocyte transendothelial migration",
  "04810"="Regulation of actin cytoskeleton","04910"="Insulin signaling pathway",
  "04914"="Progesterone-mediated oocyte maturation","05010"="Alzheimer disease",
  "05012"="Parkinson disease","05014"="Amyotrophic lateral sclerosis",
  "05016"="Huntington disease","05020"="Prion diseases","05140"="Leishmaniasis",
  "05142"="Chagas disease (American trypanosomiasis)","05144"="Malaria",
  "05145"="Toxoplasmosis","05146"="Amoebiasis","05160"="Hepatitis C",
  "05161"="Hepatitis B","05164"="Influenza A","05166"="HTLV-I infection",
  "05169"="Epstein-Barr virus infection","05200"="Pathways in cancer",
  "05203"="Viral carcinogenesis","05205"="Proteoglycans in cancer",
  "05219"="Bladder cancer","05222"="Small cell lung cancer",
  "05223"="Non-small cell lung cancer")
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
## STEP 1: DEG for POP (GSE267852, VFB only) - DESeq2 primary, voom+limma secondary
## =============================================================================
cat("\n================ STEP 1: DEG for POP (GSE267852, VFB, RNA-seq) ================\n\n")

meta <- read.csv("data/GSE267852_sample_metadata.csv")
meta_use <- meta[meta$CellType == "VFB", ]
cat("Samples used (vaginal fibroblasts, VFB, only):", nrow(meta_use), "(",
    sum(meta_use$Group == "Ctrl"), "Ctrl /", sum(meta_use$Group == "POP"), "POP )\n")
cat("Excluded: 12 AMSC (adipose-derived stem cell) samples - different tissue origin.\n\n")

counts_raw <- read.delim("data/GSE267852_raw_counts.tsv", check.names = FALSE)
count_cols <- paste0(meta_use$JC, "_count")
stopifnot(all(count_cols %in% colnames(counts_raw)))
counts <- as.matrix(counts_raw[, count_cols])
rownames(counts) <- counts_raw$GeneName
colnames(counts) <- meta_use$GSM
cat("Raw count matrix:", nrow(counts), "gene rows (Ensembl-annotated, symbol names) x",
    ncol(counts), "samples\n")

# Collapse duplicate gene symbols by summing counts (5,579 duplicated names in
# the full annotation - multiple Ensembl IDs sharing one symbol).
counts <- rowsum(counts, group = rownames(counts))
cat("After collapsing duplicate gene symbols:", nrow(counts), "genes\n\n")

group_vfb <- factor(meta_use$Group, levels = c("Ctrl", "POP"))
design_vfb <- model.matrix(~group_vfb)

cat("Design: 6 controls vs 6 POP, vaginal fibroblasts, PAIRED cell type per\n")
cat("subject in the original study but UNPAIRED here (patient identity across\n")
cat("VFB and ASC samples is not given in the public metadata) - a standard\n")
cat("unpaired two-group comparison. Two RNA-seq DEG methods reported side by\n")
cat("side, as in script 02: DESeq2 (primary, matches the original study's own\n")
cat("method) and edgeR+voom+limma (secondary, reused below for GSEA ranking).\n\n")

## --- DESeq2 (primary) --------------------------------------------------------
counts_int <- counts; storage.mode(counts_int) <- "integer"
coldata_vfb <- data.frame(row.names = meta_use$GSM, group = group_vfb)
dds <- DESeqDataSetFromMatrix(countData = counts_int, colData = coldata_vfb, design = ~group)
dds <- dds[rowSums(counts(dds) >= 10) >= 6, ]
dds <- DESeq(dds, quiet = TRUE)
res_deseq2 <- as.data.frame(results(dds, contrast = c("group", "POP", "Ctrl"), alpha = 0.05))
res_deseq2$Gene <- rownames(res_deseq2)
res_deseq2 <- res_deseq2[order(res_deseq2$pvalue), ]

vfb_full <- res_deseq2[, c("Gene", "log2FoldChange", "baseMean", "stat", "pvalue", "padj")]
colnames(vfb_full) <- c("Gene", "logFC", "baseMean", "stat", "P.Value", "adj.P.Val")
write.csv(vfb_full, "results/13_POP_GSE267852_VFB_DESeq2_full.csv", row.names = FALSE)

vfb_deg <- subset(vfb_full, !is.na(adj.P.Val) & adj.P.Val < 0.05 & abs(logFC) > 1)
vfb_deg <- vfb_deg[order(vfb_deg$adj.P.Val), ]
write.csv(vfb_deg, "results/13_POP_GSE267852_VFB_DEG_logFC1_FDR05.csv", row.names = FALSE)

cat("=== DESeq2 (primary DEG method) ===\n")
cat("Genes kept (count>=10 in >=6 samples):", nrow(dds), "of", nrow(counts), "\n")
cat("DEG (|log2FC|>1, FDR<0.05):", nrow(vfb_deg), "(", sum(vfb_deg$logFC > 0), "up /",
    sum(vfb_deg$logFC < 0), "down )\n")
cat("Minimum FDR reached:", signif(min(vfb_full$adj.P.Val, na.rm = TRUE), 3), "\n")
cat("Genes at raw p<0.01 (the original study's own less-stringent criterion):",
    sum(vfb_full$P.Value < 0.01, na.rm = TRUE), "of", nrow(vfb_full), "\n\n")
print(head(vfb_full, 10))
cat("\n")

## --- edgeR + voom + limma (secondary - reused below for GSEA ranking) ------
dge <- DGEList(counts = counts, group = group_vfb)
keep_expr <- filterByExpr(dge, design_vfb)
dge <- dge[keep_expr, , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge, method = "TMM")
voom_fit <- voom(dge, design_vfb)
fit_vfb <- eBayes(lmFit(voom_fit, design_vfb))
vfb_full_limma <- topTable(fit_vfb, coef = "group_vfbPOP", number = Inf, sort.by = "P")
vfb_full_limma$Gene <- rownames(vfb_full_limma)
vfb_full_limma <- vfb_full_limma[, c("Gene", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val")]
write.csv(vfb_full_limma, "results/13_POP_GSE267852_VFB_voom_limma_full.csv", row.names = FALSE)
vfb_deg_limma <- subset(vfb_full_limma, adj.P.Val < 0.05 & abs(logFC) > 1)

cat("=== edgeR+voom+limma (secondary - only used to rank genes for GSEA below) ===\n")
cat("Genes kept after filterByExpr:", nrow(dge), "of", nrow(counts), "\n")
cat("DEG (|log2FC|>1, FDR<0.05):", nrow(vfb_deg_limma), "\n")
cat("Minimum FDR reached:", signif(min(vfb_full_limma$adj.P.Val), 3), "\n\n")

cat("SUMMARY: DESeq2 finds", nrow(vfb_deg), "DEG(s) where voom+limma finds",
    nrow(vfb_deg_limma), "\n")
if (nrow(vfb_deg) == 0) {
  cat("This matches the original study's own published result (0 DEG at\n")
  cat("FDR<0.05 in VFB, their abstract) - a genuine, expected null result at\n")
  cat("this threshold, not a pipeline problem. GSEA below tests for a\n")
  cat("coordinated, whole-transcriptome signal instead.\n")
}
cat("\n")


## =============================================================================
## STEP 2: GSEA for POP (GSE267852 VFB) - CLASSIC (phenotype permutation)
## =============================================================================
cat("================ STEP 2: GSEA for POP (GSE267852 VFB) - CLASSIC method ================\n\n")
cat("Method: classic GSEA, phenotype permutation (shuffling the 12 group\n")
cat("labels, unpaired, BALANCED 6 vs 6). choose(12,6)=924 possible\n")
cat("relabelings - better resolution than GSE208261's 6-vs-12 (18,564 total\n")
cat("but more asymmetric) and vastly better than SUI's 3-vs-3 (20 total).\n")
cat("Ranking = moderated t-statistic from Step 1's secondary edgeR+voom+limma\n")
cat("fit (not DESeq2 - see script 02 for why).\n\n")

set.seed(267852)
ranked_full <- fit_vfb$t[, "group_vfbPOP"]
ord2 <- order(-ranked_full)
ranked_genes_vfb <- names(ranked_full)[ord2]
ranked_scores_vfb <- ranked_full[ord2]
N_vfb <- length(ranked_genes_vfb)

ann_vfb <- suppressWarnings(select(org.Hs.eg.db, keys = ranked_genes_vfb, keytype = "SYMBOL", columns = "PATH"))
ann_vfb <- ann_vfb[!is.na(ann_vfb$PATH), ]
gs_sizes_vfb <- table(ann_vfb$PATH)
valid_paths_vfb <- names(gs_sizes_vfb)[gs_sizes_vfb >= 5 & gs_sizes_vfb <= 200]
gene_sets_vfb <- split(ann_vfb$SYMBOL[ann_vfb$PATH %in% valid_paths_vfb], ann_vfb$PATH[ann_vfb$PATH %in% valid_paths_vfb])
cat("KEGG pathways tested (5-200 members):", length(gene_sets_vfb), "\n")

hit_idx_vfb <- lapply(gene_sets_vfb, function(g) which(ranked_genes_vfb %in% g))
hit_idx_vfb <- hit_idx_vfb[sapply(hit_idx_vfb, length) >= 3]
cat("Pathways with >=3 genes in the ranked list:", length(hit_idx_vfb), "\n")
es_obs_vfb <- sapply(hit_idx_vfb, calc_es, scores_abs = abs(ranked_scores_vfb), N = N_vfb)

n_perm <- 500
cat("Running", n_perm, "phenotype (group-label) permutations...\n")
t0 <- Sys.time()
gene_sets_syms_vfb <- lapply(hit_idx_vfb, function(idx) ranked_genes_vfb[idx])
perm_es_vfb <- matrix(NA_real_, nrow = n_perm, ncol = length(hit_idx_vfb))
for (i in seq_len(n_perm)) {
  perm_group <- sample(group_vfb)
  perm_design <- model.matrix(~perm_group)
  perm_fit <- eBayes(lmFit(voom_fit, perm_design))
  t_perm <- perm_fit$t[, 2]
  rank_of_gene <- rank(-t_perm, ties.method = "first")
  scores_abs_sorted <- sort(abs(t_perm), decreasing = TRUE)
  for (j in seq_along(gene_sets_syms_vfb)) {
    hidx <- rank_of_gene[gene_sets_syms_vfb[[j]]]
    if (length(hidx) >= 3) perm_es_vfb[i, j] <- calc_es(hidx, scores_abs_sorted, N_vfb)
  }
}
cat("Done in", round(difftime(Sys.time(), t0, units = "secs"), 1), "seconds\n\n")

pval_vfb <- numeric(length(hit_idx_vfb)); nes_vfb <- numeric(length(hit_idx_vfb))
for (j in seq_along(hit_idx_vfb)) {
  pe <- perm_es_vfb[, j]; pe <- pe[!is.na(pe)]
  if (es_obs_vfb[j] >= 0) {
    pval_vfb[j] <- (sum(pe >= es_obs_vfb[j]) + 1) / (length(pe) + 1)
    base <- mean(pe[pe >= 0]); if (is.nan(base) || base == 0) base <- NA
  } else {
    pval_vfb[j] <- (sum(pe <= es_obs_vfb[j]) + 1) / (length(pe) + 1)
    base <- mean(abs(pe[pe < 0])); if (is.nan(base) || base == 0) base <- NA
  }
  nes_vfb[j] <- es_obs_vfb[j] / base
}

leading_edge_vfb <- sapply(hit_idx_vfb, function(idx) paste(ranked_genes_vfb[idx], collapse = "/"))
gsea_vfb <- data.frame(PATH = names(hit_idx_vfb), Nh = sapply(hit_idx_vfb, length),
                        ES = es_obs_vfb, NES = nes_vfb, pvalue = pval_vfb, leadingEdge = leading_edge_vfb)
gsea_vfb$p.adjust <- p.adjust(gsea_vfb$pvalue, "BH")
gsea_vfb$PathwayName <- kegg_label(gsea_vfb$PATH)
gsea_vfb <- gsea_vfb[order(gsea_vfb$pvalue), ]
write.csv(gsea_vfb, "results/14_GSEA_classic_POP_GSE267852VFB_KEGG.csv", row.names = FALSE)

cat("=== GSEA classic, POP GSE267852 VFB (KEGG) ===\n")
cat("Pathways significant at FDR<0.25:", sum(gsea_vfb$p.adjust < 0.25, na.rm = TRUE), "of", nrow(gsea_vfb), "\n")
cat("Pathways significant at FDR<0.05:", sum(gsea_vfb$p.adjust < 0.05, na.rm = TRUE), "\n\n")
print(head(gsea_vfb[, c("PathwayName","Nh","NES","pvalue","p.adjust")], 10))
cat("\n")


## =============================================================================
## STEP 3: DEG and GSEA for SUI (Wei 2020) - same method as script 02
## =============================================================================
cat("================ STEP 3: DEG and GSEA for SUI (Wei 2020) ================\n\n")
cat("Re-derived here (fully standalone script) - identical method to\n")
cat("POP_SUI_analysis.R and script 02; see those for full rationale.\n\n")

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

n_perm2 <- 1000
cat("Running", n_perm2, "gene-set-label permutations for SUI (preranked GSEA)...\n")
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

cat("=== GSEA preranked, SUI (KEGG) ===\n")
cat("Pathways significant at FDR<0.25:", sum(gsea_sui$p.adjust < 0.25, na.rm = TRUE), "of", nrow(gsea_sui), "\n")
cat("Pathways significant at FDR<0.05:", sum(gsea_sui$p.adjust < 0.05, na.rm = TRUE), "\n\n")


## =============================================================================
## STEP 4: Shared pathways (POP GSE267852 VFB x SUI) + gene-level direction table
## =============================================================================
cat("================ STEP 4: Shared pathways GSE267852(VFB) x SUI ================\n\n")

report_shared4 <- function(fdr_cut) {
  pop_sig <- subset(gsea_vfb, p.adjust < fdr_cut)
  sui_sig <- subset(gsea_sui, p.adjust < fdr_cut)
  shared_ids <- intersect(pop_sig$PATH, sui_sig$PATH)
  cat("--- FDR <", fdr_cut, ": POP(GSE267852 VFB) significant =", nrow(pop_sig),
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

shared4_025 <- report_shared4(0.25)
shared4_005 <- report_shared4(0.05)
if (nrow(shared4_025) > 0) write.csv(shared4_025, "results/15_shared_pathways_GSE267852VFBxSUI_FDR025.csv", row.names = FALSE)
if (nrow(shared4_005) > 0) write.csv(shared4_005, "results/15_shared_pathways_GSE267852VFBxSUI_FDR005.csv", row.names = FALSE)

shared4_for_table <- if (nrow(shared4_005) > 0) shared4_005 else shared4_025
cat("\n--- Gene-level direction table for the shared pathways ---\n")
if (nrow(shared4_for_table) == 0) {
  cat("No shared pathway at either threshold - no gene-level table to build.\n\n")
  gene_dir_table4 <- data.frame()
} else {
  pop_lookup4 <- setNames(vfb_full$logFC, vfb_full$Gene)
  sui_lookup4 <- setNames(sui_full$logFC, sui_full$GeneSymbol)
  pop4_tested <- ranked_genes_vfb
  sui_tested4 <- ranked_genes_sui

  rows4 <- lapply(seq_len(nrow(shared4_for_table)), function(i) {
    pid <- shared4_for_table$PATH[i]
    pname <- shared4_for_table$PathwayName[i]
    full_members <- suppressWarnings(select(org.Hs.eg.db, keys = pid, keytype = "PATH", columns = "SYMBOL")$SYMBOL)
    members <- intersect(unique(full_members), union(pop4_tested, sui_tested4))
    data.frame(
      Pathway = pname, KEGG_ID = pid, Gene = members,
      logFC_POP = unname(pop_lookup4[members]), logFC_SUI = unname(sui_lookup4[members]),
      stringsAsFactors = FALSE)
  })
  gene_dir_table4 <- do.call(rbind, rows4)
  gene_dir_table4$Direction_POP <- ifelse(is.na(gene_dir_table4$logFC_POP), "not tested",
                                           ifelse(gene_dir_table4$logFC_POP > 0, "up", "down"))
  gene_dir_table4$Direction_SUI <- ifelse(is.na(gene_dir_table4$logFC_SUI), "not tested",
                                           ifelse(gene_dir_table4$logFC_SUI > 0, "up", "down"))
  gene_dir_table4$Concordant <- with(gene_dir_table4,
                                      Direction_POP != "not tested" & Direction_SUI != "not tested" &
                                        Direction_POP == Direction_SUI)
  gene_dir_table4 <- gene_dir_table4[order(gene_dir_table4$Pathway, -gene_dir_table4$Concordant), ]
  write.csv(gene_dir_table4, "results/15_shared_pathways_gene_direction_table.csv", row.names = FALSE)

  both_tested4 <- subset(gene_dir_table4, Direction_POP != "not tested" & Direction_SUI != "not tested")
  cat("Genes in shared pathways tested in both diseases:", nrow(both_tested4),
      "| concordant direction:", sum(both_tested4$Concordant),
      "(", round(100 * mean(both_tested4$Concordant), 1), "% )\n\n")
  print(head(gene_dir_table4, 20))
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

vfb_full$Gene_label <- vfb_full$Gene
deg_note <- if (nrow(vfb_deg) == 0) "none reach FDR<0.05" else paste0(nrow(vfb_deg), " gene(s) reach FDR<0.05")
deg_note_full <- if (nrow(vfb_deg) == 0) "none reach FDR<0.05" else
  paste0(nrow(vfb_deg), " gene(s) reach FDR<0.05: ", paste(vfb_deg$Gene, collapse = ", "))
cat("Volcano plot DEG note:", deg_note_full, "\n\n")
v_vfb <- make_volcano(vfb_full, "logFC", "P.Value", "Gene_label",
                       paste0("Volcano plot - POP (GSE267852, vaginal fibroblasts, DESeq2)\n",
                              "POP vs control\n",
                              "(points shown at raw p<0.05, |log2FC|>1; ", deg_note, ")"))
ggsave("figures/19_volcano_POP_GSE267852VFB.png", v_vfb, width = 9.5, height = 6.5, dpi = 300)
cat("Saved: figures/19_volcano_POP_GSE267852VFB.png\n\n")

top_vfb_genes <- head(vfb_full$Gene, 40)
logcpm_vfb <- voom_fit$E[rownames(voom_fit$E) %in% top_vfb_genes, , drop = FALSE]
ann_col_vfb <- data.frame(Group = meta_use$Group, row.names = meta_use$GSM)
ord_vfb <- order(meta_use$Group)
png("figures/20_heatmap_top_POP_GSE267852VFB.png", width = 2600, height = 3200, res = 300)
pheatmap(logcpm_vfb[, ord_vfb], scale = "row",
         annotation_col = ann_col_vfb[ord_vfb, , drop = FALSE],
         cluster_cols = FALSE, cluster_rows = TRUE,
         main = paste0("Top 40 genes by p-value (DESeq2) - POP (GSE267852, VFB)\n(row z-score, log-CPM; ", deg_note, ")"),
         fontsize = 7, show_colnames = TRUE, fontsize_col = 7)
dev.off()
cat("Saved: figures/20_heatmap_top_POP_GSE267852VFB.png\n\n")

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
p_bar_vfb <- make_gsea_barplot(gsea_vfb, "GSEA (classic) - top KEGG pathways in POP (GSE267852, VFB)")
ggsave("figures/21_GSEA_barplot_POP_GSE267852VFB.png", p_bar_vfb, width = 10, height = 6, dpi = 300)
cat("Saved: figures/21_GSEA_barplot_POP_GSE267852VFB.png\n\n")

compare_pool4 <- if (nrow(shared4_025) > 0) shared4_025 else data.frame()
if (nrow(compare_pool4) > 0) {
  p_shared4 <- ggplot(compare_pool4, aes(x = NES_POP, y = NES_SUI, label = PathwayName)) +
    geom_hline(yintercept = 0, color = "grey70") + geom_vline(xintercept = 0, color = "grey70") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
    geom_point(aes(color = Same_direction), size = 3) +
    scale_color_manual(values = c("TRUE" = "#2166AC", "FALSE" = "#B2182B"),
                        labels = c("TRUE" = "Same direction", "FALSE" = "Opposite direction")) +
    geom_text_repel(size = 3, max.overlaps = 30) +
    labs(title = "Shared pathways (FDR<0.25): NES in POP (GSE267852 VFB) vs NES in SUI",
         x = "NES - POP GSE267852 VFB (classic GSEA)", y = "NES - SUI (preranked GSEA)", color = NULL) +
    theme_bw() + theme(legend.position = "top")
  ggsave("figures/22_shared_pathways_NES_comparison_GSE267852VFB.png", p_shared4, width = 9, height = 7, dpi = 300)
  cat("Saved: figures/22_shared_pathways_NES_comparison_GSE267852VFB.png\n\n")
} else {
  cat("SKIPPED: figures/22_shared_pathways_NES_comparison_GSE267852VFB.png (no shared pathway at FDR<0.25)\n\n")
}

if (nrow(gene_dir_table4) > 0) {
  both_tested4 <- subset(gene_dir_table4, Direction_POP != "not tested" & Direction_SUI != "not tested")
  if (nrow(both_tested4) >= 2) {
    fc_mat4 <- as.matrix(both_tested4[, c("logFC_POP", "logFC_SUI")])
    dir_mat4 <- ifelse(fc_mat4 > 0, 1, -1)
    rownames(dir_mat4) <- rownames(fc_mat4) <- make.unique(both_tested4$Gene)
    ord4b <- order(both_tested4$Pathway, -both_tested4$Concordant, -abs(rowSums(fc_mat4, na.rm = TRUE)))
    dir_mat4 <- dir_mat4[ord4b, , drop = FALSE]; fc_mat4 <- fc_mat4[ord4b, , drop = FALSE]
    dir_mat4 <- head(dir_mat4, 50); fc_mat4 <- head(fc_mat4, 50)
    ann_row4 <- data.frame(Pathway = both_tested4$Pathway[ord4b][seq_len(nrow(dir_mat4))],
                            row.names = rownames(dir_mat4))
    colnames(dir_mat4) <- colnames(fc_mat4) <- c("POP", "SUI")
    png("figures/23_gene_direction_heatmap_GSE267852VFBxSUI.png",
        width = 2600, height = max(1800, 45 * nrow(dir_mat4)), res = 300)
    pheatmap(dir_mat4, cluster_cols = FALSE, cluster_rows = FALSE,
             annotation_row = ann_row4,
             color = c("#2166AC", "#B2182B"), breaks = c(-2, 0, 2),
             legend = FALSE,
             display_numbers = matrix(sprintf("%.2f", fc_mat4), nrow(fc_mat4)),
             number_color = "white", fontsize_number = 7,
             main = "Gene direction, shared pathways (POP GSE267852 VFB vs SUI)\n(fill = up/down; number = log2 Fold Change)",
             fontsize = 7, fontsize_row = 6)
    dev.off()
    cat("Saved: figures/23_gene_direction_heatmap_GSE267852VFBxSUI.png\n\n")
  } else {
    cat("SKIPPED: figures/23_gene_direction_heatmap_GSE267852VFBxSUI.png (fewer than 2 genes tested in both)\n\n")
  }
} else {
  cat("SKIPPED: figures/23_gene_direction_heatmap_GSE267852VFBxSUI.png (no shared-pathway gene table)\n\n")
}

cat("\n=== ALL GSE267852 VFB FIGURES SAVED TO figures/ ===\n")
cat("19_volcano_POP_GSE267852VFB.png\n20_heatmap_top_POP_GSE267852VFB.png\n")
cat("21_GSEA_barplot_POP_GSE267852VFB.png\n22_shared_pathways_NES_comparison_GSE267852VFB.png (if any shared)\n")
cat("23_gene_direction_heatmap_GSE267852VFBxSUI.png (if any shared)\n")
cat("\n=== DONE (script 04) ===\n")