# =============================================================================
# GSEA (preranked, KEGG) on the ORIGINAL POP (GSE208261) and SUI (GSE149072,
# 36h) RNA-seq datasets - the pair that produced the "only 4 common genes"
# finding (INPP4B, ECM1, BEND3, KREMEN1) the professor's suggestion refers
# to. Run separately as requested; the point is to find pathways
# significantly enriched in BOTH, since only 4 shared genes is too few for
# ORA/GO/KEGG on gene overlap directly.
#
# WHY THIS IS THE STRONGEST GSEA IN THE WHOLE PROJECT: unlike the small
# literature panels (58-90 genes, no raw data, gene-set-label permutation
# only) or even the GSE53868 microarray (real transcriptome but the SUI
# side there was still a literature panel), BOTH datasets here are real
# RNA-seq with the FULL transcriptome tested (26,810 human genes / 12,807
# rat genes) AND real per-sample raw counts available for BOTH conditions.
# This is exactly the situation GSEA was designed for, and a genuine
# phenotype permutation null (not the weaker gene-set-label permutation)
# is used for both sides.
#
# Reuses the exact filtering/model already used for the official pipeline
# results (results/POP_GSE208261_DESeq2_completo.csv,
# results/SUI_36hr_DESeq2_completo.csv - see
# R/analysis_SUI_POP_crosscomparison.R for the original derivation):
#   POP:  rowSums(counts>=10)>=6, design ~ age_group + condition
#   SUI:  rowSums(counts>=10)>=3, design ~ group (36h Treated vs Untreated)
# Both DESeq2 "stat" columns (real Wald statistics, unaffected by
# lfcShrink) are reused directly as the GSEA ranking metric - no need to
# refit DESeq2 for the observed ranking.
#
# PERMUTATION NULL (phenotype permutation, real, per dataset): refitting
# full DESeq2 hundreds of times is too slow, so a fast, standard proxy is
# used instead - log2(size-factor-normalized counts + 1), per-gene Welch
# t-statistic, recomputed under group-label permutation. This preserves
# real per-sample structure (unlike gene-set-label permutation) at a
# fraction of the cost of refitting the negative-binomial GLM each time.
# One simplification, flagged here rather than hidden: the permutation
# null for POP does not carry the age_group covariate (a plain 2-group
# Welch test), while the observed ranking (POP's real "stat" column) DOES
# adjust for it - a real, if likely small, mismatch between observed and
# null model, noted openly.
# =============================================================================

suppressMessages(library(org.Hs.eg.db))
set.seed(2024)

calc_es <- function(hit_idx, scores_abs, N) {
  Nh <- length(hit_idx); Nm <- N - Nh
  if (Nm <= 0 || Nh == 0) return(NA)
  sum_hit <- sum(scores_abs[hit_idx])
  step <- rep(-1 / Nm, N)
  step[hit_idx] <- scores_abs[hit_idx] / sum_hit
  running <- cumsum(step)
  running[which.max(abs(running))]
}

fast_welch_t <- function(logmat, group) {
  # vectorized per-gene Welch t-statistic (group2 - group1), for the
  # permutation null only
  lv <- levels(group)
  g1 <- logmat[, group == lv[1], drop = FALSE]
  g2 <- logmat[, group == lv[2], drop = FALSE]
  n1 <- ncol(g1); n2 <- ncol(g2)
  m1 <- rowMeans(g1); m2 <- rowMeans(g2)
  v1 <- rowSums((g1 - m1)^2) / (n1 - 1)
  v2 <- rowSums((g2 - m2)^2) / (n2 - 1)
  se <- sqrt(v1 / n1 + v2 / n2)
  (m2 - m1) / se
}

run_gsea_with_permutation_matrix <- function(gene_ids, observed_stat, logmat, group, gene_sets,
                                              min_gs = 5, max_gs = 200, n_perm = 500) {
  N <- length(gene_ids)
  ord <- order(-observed_stat)
  ranked_ids <- gene_ids[ord]
  ranked_scores <- observed_stat[ord]
  scores_abs <- abs(ranked_scores)

  hit_idx_list <- lapply(gene_sets, function(g) which(ranked_ids %in% g))
  keep <- sapply(hit_idx_list, length) >= min_gs & sapply(hit_idx_list, length) <= max_gs
  hit_idx_list <- hit_idx_list[keep]
  cat("  Pathways testable (", min_gs, "-", max_gs, " members present):", length(hit_idx_list), "\n", sep = "")
  if (length(hit_idx_list) == 0) return(NULL)

  es_obs <- sapply(hit_idx_list, calc_es, scores_abs = scores_abs, N = N)

  gene_sets_ids <- lapply(hit_idx_list, function(idx) ranked_ids[idx])

  perm_mat <- matrix(NA_real_, nrow = n_perm, ncol = length(hit_idx_list))
  rownames_logmat <- rownames(logmat)
  for (i in seq_len(n_perm)) {
    perm_group <- sample(group)
    t_perm <- fast_welch_t(logmat, perm_group)
    names(t_perm) <- rownames_logmat
    t_perm <- t_perm[ranked_ids]  # align to the SAME gene order as observed ranking
    rank_of_gene <- rank(-t_perm, ties.method = "first")
    scores_abs_sorted <- sort(abs(t_perm), decreasing = TRUE)
    for (j in seq_along(gene_sets_ids)) {
      hidx <- rank_of_gene[as.character(gene_sets_ids[[j]])]
      if (length(hidx) > 0 && !any(is.na(hidx))) {
        perm_mat[i, j] <- calc_es(hidx, scores_abs_sorted, N)
      }
    }
  }

  pval <- numeric(length(hit_idx_list)); nes <- numeric(length(hit_idx_list))
  for (j in seq_along(hit_idx_list)) {
    pe <- perm_mat[, j]; pe <- pe[!is.na(pe)]
    if (length(pe) == 0) { pval[j] <- NA; nes[j] <- NA; next }
    if (es_obs[j] >= 0) {
      pval[j] <- (sum(pe >= es_obs[j]) + 1) / (length(pe) + 1)
      base <- mean(pe[pe >= 0]); if (is.nan(base) || base == 0) base <- NA
    } else {
      pval[j] <- (sum(pe <= es_obs[j]) + 1) / (length(pe) + 1)
      base <- mean(abs(pe[pe < 0])); if (is.nan(base) || base == 0) base <- NA
    }
    nes[j] <- es_obs[j] / base
  }

  leading_edge <- sapply(hit_idx_list, function(idx) paste(ranked_ids[idx], collapse = "/"))
  res <- data.frame(PATH = names(hit_idx_list), Nh = sapply(hit_idx_list, length),
                     ES = es_obs, NES = nes, pvalue = pval, leadingEdge = leading_edge)
  res$p.adjust <- p.adjust(res$pvalue, "BH")
  res[order(res$pvalue), ]
}

## KEGG gene sets by Entrez ID (offline, org.Hs.eg.db)
all_entrez <- keys(org.Hs.eg.db, keytype = "ENTREZID")
ann_path <- suppressWarnings(select(org.Hs.eg.db, keys = all_entrez, keytype = "ENTREZID", columns = "PATH"))
ann_path <- ann_path[!is.na(ann_path$PATH), ]
kegg_sets <- split(ann_path$ENTREZID, ann_path$PATH)
cat("KEGG pathway definitions loaded:", length(kegg_sets), "\n\n")

## =============================================================================
## PART 1 - POP (GSE208261)
## =============================================================================
cat("============================================================\n")
cat("PART 1: POP (GSE208261) preranked GSEA\n")
cat("============================================================\n\n")

pop_counts <- read.delim("data/GSE208261_raw_counts_GRCh38.p13_NCBI.tsv", row.names = 1, check.names = FALSE)
pop_counts <- round(as.matrix(pop_counts)); mode(pop_counts) <- "integer"

sample_labels <- c(
  GSM6339911="Control_D1", GSM6339912="Control_D2", GSM6339913="Control_D3",
  GSM6339914="Control_D4", GSM6339915="Control_D5", GSM6339916="Control_D6",
  GSM6339917="Control_Y1", GSM6339918="Control_Y2", GSM6339919="Control_Y3",
  GSM6339920="Control_Y4", GSM6339921="Control_Y5", GSM6339922="Control_Y6",
  GSM6339923="POP_D1", GSM6339924="POP_D2", GSM6339925="POP_D3",
  GSM6339926="POP_D4", GSM6339927="POP_D5", GSM6339928="POP_D6",
  GSM6339929="POP_Y1", GSM6339930="POP_Y2", GSM6339931="POP_Y3",
  GSM6339932="POP_Y4", GSM6339933="POP_Y5", GSM6339934="POP_Y6"
)
label <- sample_labels[colnames(pop_counts)]
condition <- factor(ifelse(grepl("^POP", label), "POP", "Control"), levels = c("Control", "POP"))

keep_pop <- rowSums(pop_counts >= 10) >= 6
pop_counts_sub <- pop_counts[keep_pop, ]

sf_pop <- DESeq2::estimateSizeFactorsForMatrix(pop_counts_sub)
pop_logmat <- log2(sweep(pop_counts_sub, 2, sf_pop, "/") + 1)

pop_full <- read.csv("results/POP_GSE208261_DESeq2_completo.csv")
pop_full$Human_Entrez_ID <- as.character(pop_full$Human_Entrez_ID)
pop_full <- pop_full[!is.na(pop_full$stat) & pop_full$Human_Entrez_ID %in% rownames(pop_logmat), ]
pop_full <- pop_full[!duplicated(pop_full$Human_Entrez_ID), ]
pop_logmat <- pop_logmat[pop_full$Human_Entrez_ID, ]
stopifnot(identical(rownames(pop_logmat), pop_full$Human_Entrez_ID))

cat("POP ranked gene list:", nrow(pop_full), "genes\n")
pop_gsea <- run_gsea_with_permutation_matrix(
  gene_ids = pop_full$Human_Entrez_ID, observed_stat = pop_full$stat,
  logmat = pop_logmat, group = condition, gene_sets = kegg_sets
)
write.csv(pop_gsea, "results/GSEA_POP_GSE208261_RNAseq_KEGG.csv", row.names = FALSE)
cat("POP GSEA - significant FDR<0.05:", sum(pop_gsea$p.adjust < 0.05, na.rm = TRUE),
    "| FDR<0.25:", sum(pop_gsea$p.adjust < 0.25, na.rm = TRUE), "of", nrow(pop_gsea), "\n\n")

## =============================================================================
## PART 2 - SUI 36h (GSE149072, rat -> human ortholog via babelgene)
## =============================================================================
cat("============================================================\n")
cat("PART 2: SUI 36h (GSE149072) preranked GSEA\n")
cat("============================================================\n\n")

sui_counts <- read.csv("data/GSE149072_rawCounts.csv", row.names = 1, check.names = FALSE)
sui_counts <- round(as.matrix(sui_counts)); mode(sui_counts) <- "integer"
untreated <- grep("Rat_Urethra_Untreated_36hr", colnames(sui_counts), value = TRUE)
treated   <- grep("Rat_Urethra_Treated_36hr", colnames(sui_counts), value = TRUE)
sel <- c(untreated, treated)
sui_counts_sub <- sui_counts[, sel]
group_sui <- factor(c(rep("Untreated", length(untreated)), rep("Treated", length(treated))),
                     levels = c("Untreated", "Treated"))

keep_sui <- rowSums(sui_counts_sub >= 10) >= 3
sui_counts_sub <- sui_counts_sub[keep_sui, ]

sf_sui <- DESeq2::estimateSizeFactorsForMatrix(sui_counts_sub)
sui_logmat <- log2(sweep(sui_counts_sub, 2, sf_sui, "/") + 1)

load("data/babelgene_orthologs.rda")
rat_orth_all <- unique(orthologs_df[orthologs_df$taxon_id == 10116,
                                     c("symbol", "human_symbol", "human_entrez")])
rat_n <- table(rat_orth_all$symbol); hum_n <- table(rat_orth_all$human_symbol)
rat_orth <- rat_orth_all[rat_n[rat_orth_all$symbol] == 1 & hum_n[rat_orth_all$human_symbol] == 1, ]
colnames(rat_orth) <- c("Rat_Gene_Symbol", "Human_Ortholog_Symbol", "Human_Entrez_ID")

sui_full <- read.csv("results/SUI_36hr_DESeq2_completo.csv")
sui_full <- sui_full[!is.na(sui_full$stat) & sui_full$Rat_Gene_Symbol %in% rownames(sui_logmat), ]
sui_full <- merge(sui_full, rat_orth, by = "Rat_Gene_Symbol")
sui_full$Human_Entrez_ID <- as.character(sui_full$Human_Entrez_ID)
sui_full <- sui_full[!duplicated(sui_full$Human_Entrez_ID), ]

sui_logmat_orth <- sui_logmat[sui_full$Rat_Gene_Symbol, ]
rownames(sui_logmat_orth) <- sui_full$Human_Entrez_ID  # relabel rows by human ortholog for gene-set lookup
stopifnot(identical(rownames(sui_logmat_orth), sui_full$Human_Entrez_ID))

cat("SUI 36h ranked gene list (1:1 human orthologs):", nrow(sui_full), "of",
    sum(!is.na(read.csv("results/SUI_36hr_DESeq2_completo.csv")$stat)), "rat genes tested\n")
sui_gsea <- run_gsea_with_permutation_matrix(
  gene_ids = sui_full$Human_Entrez_ID, observed_stat = sui_full$stat,
  logmat = sui_logmat_orth, group = group_sui, gene_sets = kegg_sets
)
write.csv(sui_gsea, "results/GSEA_SUI36h_GSE149072_RNAseq_KEGG.csv", row.names = FALSE)
cat("SUI 36h GSEA - significant FDR<0.05:", sum(sui_gsea$p.adjust < 0.05, na.rm = TRUE),
    "| FDR<0.25:", sum(sui_gsea$p.adjust < 0.25, na.rm = TRUE), "of", nrow(sui_gsea), "\n\n")

## =============================================================================
## PART 3 - Shared pathways
## =============================================================================
cat("============================================================\n")
cat("PART 3: Pathways significant in BOTH POP and SUI 36h\n")
cat("============================================================\n\n")

pop_sig <- subset(pop_gsea, p.adjust < 0.25)
sui_sig <- subset(sui_gsea, p.adjust < 0.25)
common <- intersect(pop_sig$PATH, sui_sig$PATH)
cat("POP significant (FDR<0.25):", nrow(pop_sig), "| SUI significant (FDR<0.25):", nrow(sui_sig),
    "| SHARED:", length(common), "\n")
if (length(common) > 0) {
  out <- merge(pop_sig[pop_sig$PATH %in% common, c("PATH","Nh","ES","NES","pvalue","p.adjust","leadingEdge")],
               sui_sig[sui_sig$PATH %in% common, c("PATH","Nh","ES","NES","pvalue","p.adjust","leadingEdge")],
               by = "PATH", suffixes = c("_POP","_SUI36h"))
  # same direction check: both up (ES>0) or both down (ES<0) in their own contrast
  out$Same_direction <- sign(out$ES_POP) == sign(out$ES_SUI36h)
  print(out[, c("PATH","ES_POP","ES_SUI36h","Same_direction","p.adjust_POP","p.adjust_SUI36h")])
  write.csv(out, "results/GSEA_shared_KEGG_POP_vs_SUI36h.csv", row.names = FALSE)
}

cat("\n=== END OF SCRIPT ===\n")
