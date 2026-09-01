# =============================================================================
# Script 13 — GSEA preranked no painel do Wei 2020 (SUI), vias KEGG. Mesmo
# método dos scripts 03/09 (Chen 2006/2003), com UMA diferença: como as
# Tabelas 5/6 do artigo só trazem fold change (sem p-valor por gene), o
# escore de ranking aqui é log2(fold change) diretamente, com o sinal da
# direção (up/down) - uma métrica de ranking padrão em GSEA preranked
# quando não há estatística de teste contínua disponível.
#
# Só 40 genes (20 up/20 down) - ainda MENOR que os painéis do Chen (59-69).
# Ver nota completa no cabeçalho do script 12 sobre por que este painel não
# é, na prática, maior que os anteriores (seria, se tivéssemos a Tabela
# Suplementar S2 completa do artigo, com os 7.102 mRNAs - não incluída no
# PDF fornecido).
# =============================================================================

suppressMessages(library(org.Hs.eg.db))
set.seed(2020)

wei <- read.csv("data/wei2020_top40_mRNA.csv")
wei <- wei[!duplicated(wei$GeneSymbol), ]
cat("Wei 2020 — genes candidatos:", nrow(wei),
    "(", sum(wei$Direction == "up"), "up /", sum(wei$Direction == "down"), "down )\n\n")

score <- ifelse(wei$Direction == "up", 1, -1) * log2(wei$FoldChange)
ord <- order(-score)
ranked_genes <- wei$GeneSymbol[ord]
ranked_scores <- score[ord]
N <- length(ranked_genes)

all_sym <- keys(org.Hs.eg.db, keytype = "SYMBOL")
ann_path <- suppressWarnings(select(org.Hs.eg.db, keys = all_sym, keytype = "SYMBOL", columns = "PATH"))
ann_path <- ann_path[!is.na(ann_path$PATH), ]
kegg_sets_all <- split(ann_path$SYMBOL, ann_path$PATH)
cat("Vias KEGG carregadas (universo completo, org.Hs.eg.db):", length(kegg_sets_all), "\n\n")

calc_es <- function(hit_idx, scores_abs, N) {
  Nh <- length(hit_idx); Nm <- N - Nh
  if (Nm <= 0 || Nh == 0) return(NA)
  sum_hit <- sum(scores_abs[hit_idx])
  step <- rep(-1 / Nm, N)
  step[hit_idx] <- scores_abs[hit_idx] / sum_hit
  running <- cumsum(step)
  running[which.max(abs(running))]
}

min_gs <- 2
n_perm <- 1000
testable <- lapply(kegg_sets_all, function(g) which(ranked_genes %in% g))
testable <- testable[sapply(testable, length) >= min_gs]
cat("Vias KEGG testáveis (>=", min_gs, "membros presentes nesta lista de", N, "genes):",
    length(testable), "\n")

if (length(testable) == 0) {
  cat("\nNenhuma via testável com este painel (nenhuma via KEGG tem >=2 destes 40 genes).\n")
  cat("=== FIM ===\n")
} else {
  scores_abs <- abs(ranked_scores)
  es_obs <- sapply(testable, calc_es, scores_abs = scores_abs, N = N)

  set.seed(2020)
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
  res <- res[order(res$pvalue), ]

  dir.create("results", showWarnings = FALSE)
  write.csv(res, "results/GSEA_preranked_Wei2020_KEGG.csv", row.names = FALSE)

  cat("\n=== GSEA preranked do Wei 2020 (SUI, KEGG) ===\n")
  cat("Vias significativas a FDR<0.25:", sum(res$p.adjust < 0.25, na.rm = TRUE), "\n")
  cat("Vias significativas a FDR<0.05:", sum(res$p.adjust < 0.05, na.rm = TRUE), "\n\n")
  print(head(res[, c("PATH","Nh","ES","NES","pvalue","p.adjust","leadingEdge")], 10))
}
