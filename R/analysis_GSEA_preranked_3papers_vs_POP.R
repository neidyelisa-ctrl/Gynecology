# =============================================================================
# GSEA (KEGG) for POP (GSE53868 microarray) and for THREE literature SUI
# gene panels (Chen et al. 2003, Chen et al. 2006, Tong et al. 2010), then
# the pathways each panel shares with POP.
#
# Requested explicitly: run GSEA preranked on the 3 papers using their own
# p-value + fold-change tables as the ranking metric, even though each is a
# short (60-90 gene) candidate list, not a full transcriptome. This section
# documents exactly what that means methodologically before showing any
# numbers, because it changes how the results should be read:
#
#   - POP (GSE53868): the FULL ranked transcriptome is available (31,013
#     genes after removing zero-variance probes, paired limma t-statistic).
#     This is genuine, standard preranked GSEA - identical to the analysis
#     already run and delivered (results/GSEA_POP_GSE53868_KEGG.csv),
#     reused here rather than recomputed (nothing about it changed).
#
#   - The three papers: only report a pre-filtered table of already-
#     significant candidate genes (60-90 rows), not the ~7,000-33,000
#     genes originally tested on their arrays. Running "preranked GSEA" on
#     such a short list is a NON-STANDARD use of the method - the ranking
#     and running-sum statistic are computed exactly the same way, but:
#       (a) a KEGG pathway can only be "tested" if >=2 of its members
#           happen to be among these 60-90 genes, so very few pathways are
#           testable at all (most KEGG pathways will have 0 or 1 member
#           present by chance in such a short list);
#       (b) with no raw per-sample expression matrix for these papers
#           (only the published summary table), there is no way to do
#           phenotype permutation as was done for POP. The null here uses
#           gene-set-label permutation instead - the same fallback the
#           official GSEA-Preranked tool itself defaults to whenever it is
#           given a ranked list without a matching raw expression matrix,
#           but it is the more liberal, less conservative permutation type
#           (documented in the previous script/README - correlated genes
#           can inflate it). It is used here because it is the only option
#           the data allow, not because it is preferred.
#     These results should be read as an exploratory, low-powered signal,
#     not put on the same footing as the real POP GSEA above them.
#
#   Ranking score for each paper (no raw test statistic available, so the
#   standard substitute is used): score = sign(direction) * -log10(p-value)
#   (p-values of exactly 0, as printed in Chen 2006's table, are floored to
#   half the smallest nonzero p-value in that paper before taking -log10,
#   to avoid an infinite score).
# =============================================================================

suppressMessages(library(org.Hs.eg.db))
set.seed(123)

## -----------------------------------------------------------------------
## 1) POP (GSE53868): DEGs at |logFC|>1 & FDR<0.05, and the already-computed
##    real preranked GSEA (KEGG, phenotype permutation) - reused as-is.
## -----------------------------------------------------------------------
pop <- read.csv("results/GSE53868_limma_completo.csv")
pop_deg <- subset(pop, adj.P.Val < 0.05 & abs(logFC) > 1)
pop_deg <- pop_deg[order(pop_deg$adj.P.Val), ]
write.csv(pop_deg, "results/GSE53868_DEG_logFC1_FDR05.csv", row.names = FALSE)
cat("POP (GSE53868) DEGs, |logFC|>1 & FDR<0.05:", nrow(pop_deg), "of", nrow(pop), "tested\n\n")

pop_gsea <- read.csv("results/GSEA_POP_GSE53868_KEGG.csv")
pop_gsea_sig <- subset(pop_gsea, p.adjust < 0.25)
cat("POP real preranked GSEA (KEGG, phenotype permutation, reused from before):\n")
cat(" ", nrow(pop_gsea_sig), "of", nrow(pop_gsea), "KEGG pathways significant at FDR<0.25;",
    sum(pop_gsea$p.adjust < 0.05, na.rm = TRUE), "at the stricter FDR<0.05\n\n")

## -----------------------------------------------------------------------
## 2) KEGG gene sets (offline, org.Hs.eg.db) - no external universe needed
##    for preranked GSEA (unlike ORA, N is just the length of the ranked
##    list itself, so the earlier universe-size bug does not apply here).
## -----------------------------------------------------------------------
all_sym <- keys(org.Hs.eg.db, keytype = "SYMBOL")
ann_path <- suppressWarnings(select(org.Hs.eg.db, keys = all_sym, keytype = "SYMBOL", columns = "PATH"))
ann_path <- ann_path[!is.na(ann_path$PATH), ]
kegg_sets_all <- split(ann_path$SYMBOL, ann_path$PATH)
cat("KEGG pathway definitions loaded:", length(kegg_sets_all), "pathways\n\n")

## -----------------------------------------------------------------------
## 3) Small-panel preranked GSEA: weighted ES + gene-set-label permutation.
## -----------------------------------------------------------------------
calc_es <- function(hit_idx, scores_abs, N) {
  Nh <- length(hit_idx); Nm <- N - Nh
  if (Nm <= 0 || Nh == 0) return(NA)
  sum_hit <- sum(scores_abs[hit_idx])
  step <- rep(-1 / Nm, N)
  step[hit_idx] <- scores_abs[hit_idx] / sum_hit
  running <- cumsum(step)
  running[which.max(abs(running))]
}

run_preranked_small <- function(ranked_genes, ranked_scores, kegg_sets_all, min_gs = 2, n_perm = 1000) {
  N <- length(ranked_genes)
  scores_abs <- abs(ranked_scores)
  # restrict KEGG sets to those with >=min_gs members actually present in this short ranked list
  testable <- lapply(kegg_sets_all, function(g) which(ranked_genes %in% g))
  testable <- testable[sapply(testable, length) >= min_gs]
  cat("  KEGG pathways testable (>=", min_gs, "members present in this", N, "-gene list):",
      length(testable), "\n")
  if (length(testable) == 0) return(NULL)

  es_obs <- sapply(testable, calc_es, scores_abs = scores_abs, N = N)

  # gene-set-label permutation: the ranking (and its scores) stays FIXED as
  # observed; on each permutation, redraw a random set of positions of the
  # same size as the real gene set, from 1:N, to stand in for the "hits"
  perm_mat <- matrix(NA_real_, nrow = n_perm, ncol = length(testable))
  for (i in seq_len(n_perm)) {
    for (j in seq_along(testable)) {
      hidx <- sample.int(N, length(testable[[j]]))
      perm_mat[i, j] <- calc_es(hidx, scores_abs, N)
    }
  }

  pval <- numeric(length(testable)); nes <- numeric(length(testable))
  for (j in seq_along(testable)) {
    pe <- perm_mat[, j]; pe <- pe[!is.na(pe)]
    if (is.na(es_obs[j])) { pval[j] <- NA; nes[j] <- NA; next }
    if (es_obs[j] >= 0) {
      pval[j] <- (sum(pe >= es_obs[j]) + 1) / (length(pe) + 1)
      base <- mean(pe[pe >= 0]); if (is.nan(base) || base == 0) base <- NA
    } else {
      pval[j] <- (sum(pe <= es_obs[j]) + 1) / (length(pe) + 1)
      base <- mean(abs(pe[pe < 0])); if (is.nan(base) || base == 0) base <- NA
    }
    nes[j] <- es_obs[j] / base
  }

  leading_edge <- sapply(testable, function(idx) paste(ranked_genes[idx], collapse = "/"))
  res <- data.frame(PATH = names(testable), Nh = sapply(testable, length),
                     ES = es_obs, NES = nes, pvalue = pval, leadingEdge = leading_edge)
  res$p.adjust <- p.adjust(res$pvalue, "BH")
  res[order(res$pvalue), ]
}

build_score <- function(direction, pvalue) {
  p <- pvalue
  floor_p <- min(p[p > 0], na.rm = TRUE) / 2
  p[p <= 0] <- floor_p
  sign <- ifelse(direction == "up", 1, -1)
  sign * -log10(p)
}

run_paper <- function(csv_path, symbol_col, direction_col, pvalue_col, tag) {
  cat("===", tag, "===\n")
  df <- read.csv(csv_path)
  df <- df[!duplicated(df[[symbol_col]]), ]
  score <- build_score(df[[direction_col]], df[[pvalue_col]])
  ord <- order(-score)
  ranked_genes <- df[[symbol_col]][ord]
  ranked_scores <- score[ord]
  cat("  Panel size:", length(ranked_genes), "genes\n")
  res <- run_preranked_small(ranked_genes, ranked_scores, kegg_sets_all)
  if (!is.null(res)) {
    cat("  Significant at FDR<0.25:", sum(res$p.adjust < 0.25, na.rm = TRUE),
        "| FDR<0.05:", sum(res$p.adjust < 0.05, na.rm = TRUE), "\n")
    print(head(res[, c("PATH","Nh","ES","NES","pvalue","p.adjust")], 10))
  }
  cat("\n")
  res
}

chen03_res  <- run_paper("data/chen2003_90genes.csv", "GeneSymbol", "Direction", "pvalue",       "Chen 2003 (69 genes, proliferative phase)")
chen06_res  <- run_paper("data/chen2006_79genes.csv", "GeneSymbol", "Direction", "RMA_pvalue",   "Chen 2006 (60 genes, secretory phase)")
tong10_res  <- run_paper("data/tong2010_75genes.csv", "GeneSymbol", "Direction", "pvalue",        "Tong 2010 (66 genes, postmenopausal)")

write.csv(chen03_res, "results/GSEA_preranked_Chen2003_KEGG.csv", row.names = FALSE)
write.csv(chen06_res, "results/GSEA_preranked_Chen2006_KEGG.csv", row.names = FALSE)
write.csv(tong10_res, "results/GSEA_preranked_Tong2010_KEGG.csv", row.names = FALSE)

## -----------------------------------------------------------------------
## 4) Pathways each paper shares with POP (same KEGG ID significant in
##    BOTH the paper's small-panel GSEA (FDR<0.25) and POP's real GSEA
##    (FDR<0.25) - the standard GSEA screening threshold both times, for a
##    fair, consistent comparison).
## -----------------------------------------------------------------------
report_shared <- function(paper_res, tag) {
  cat("=== Shared KEGG pathways: ", tag, " vs POP ===\n", sep = "")
  if (is.null(paper_res)) { cat("  No testable pathways for this panel.\n\n"); return(data.frame()) }
  paper_sig <- subset(paper_res, p.adjust < 0.25)
  common_ids <- intersect(paper_sig$PATH, pop_gsea_sig$PATH)
  cat("  ", tag, " significant pathways (FDR<0.25): ", nrow(paper_sig),
      " | POP significant pathways (FDR<0.25): ", nrow(pop_gsea_sig),
      " | SHARED: ", length(common_ids), "\n", sep = "")
  if (length(common_ids) == 0) { cat("\n"); return(data.frame()) }
  out <- merge(paper_sig[paper_sig$PATH %in% common_ids, c("PATH","Nh","NES","pvalue","p.adjust","leadingEdge")],
               pop_gsea_sig[pop_gsea_sig$PATH %in% common_ids, c("PATH","Nh","NES","pvalue","p.adjust","leadingEdge")],
               by = "PATH", suffixes = c(paste0("_", gsub(" ", "", tag)), "_POP"))
  print(out[, 1:5])
  cat("\n")
  out
}

shared_chen03 <- report_shared(chen03_res, "Chen2003")
shared_chen06 <- report_shared(chen06_res, "Chen2006")
shared_tong10 <- report_shared(tong10_res, "Tong2010")

if (nrow(shared_chen03) > 0) write.csv(shared_chen03, "results/GSEA_shared_KEGG_Chen2003_vs_POP.csv", row.names = FALSE)
if (nrow(shared_chen06) > 0) write.csv(shared_chen06, "results/GSEA_shared_KEGG_Chen2006_vs_POP.csv", row.names = FALSE)
if (nrow(shared_tong10) > 0) write.csv(shared_tong10, "results/GSEA_shared_KEGG_Tong2010_vs_POP.csv", row.names = FALSE)

cat("=== END OF SCRIPT ===\n")
