# =============================================================================
# Script 15 — GSEA preranked no painel COMPLETO do Wei 2020 (6.118 genes,
# Tabela Suplementar S2), vias KEGG. Mesmo método dos scripts 03/09 (Chen),
# mas com um painel ~100x maior — a maioria das vias KEGG agora tem >=2
# membros presentes, ao contrário dos painéis pequenos.
#
# Ainda usa permutação de RÓTULO DE GENE SET (não de fenótipo) porque o
# arquivo fornecido já vem PRÉ-FILTRADO pelo estudo original (só os 7.102
# genes que já passaram no corte do artigo, não os ~20.730 genes
# codificantes testados no array inteiro) - não temos como reconstruir a
# permutação de fenótipo sem o array completo. Ver nota no cabeçalho do
# script 14.
#
# Escore de ranking: sinal(direção) * -log10(p-valor) - mesma métrica dos
# scripts 03/09 (Chen), agora com p-valor REAL por gene (não é estimado/
# ausente como no painel top-40 do script 13).
# =============================================================================

suppressMessages(library(org.Hs.eg.db))
set.seed(2020)

wei <- read.csv("data/wei2020_mRNA_full.csv")
cat("Wei 2020 (completo) — genes candidatos:", nrow(wei),
    "(", sum(wei$Direction == "up"), "up /", sum(wei$Direction == "down"), "down )\n\n")

build_score <- function(direction, pvalue) {
  p <- pvalue
  floor_p <- min(p[p > 0], na.rm = TRUE) / 2
  p[p <= 0] <- floor_p
  sign <- ifelse(direction == "up", 1, -1)
  sign * -log10(p)
}

score <- build_score(wei$Direction, wei$PValue)
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

min_gs <- 5
n_perm <- 1000
testable <- lapply(kegg_sets_all, function(g) which(ranked_genes %in% g))
testable <- testable[sapply(testable, length) >= min_gs]
cat("Vias KEGG testáveis (>=", min_gs, "membros presentes nesta lista de", N, "genes):",
    length(testable), "\n")

scores_abs <- abs(ranked_scores)
es_obs <- sapply(testable, calc_es, scores_abs = scores_abs, N = N)

cat("Rodando", n_perm, "permutações de rótulo de gene set em", length(testable), "vias...\n")
t0 <- Sys.time()
set.seed(2020)
perm_mat <- matrix(NA_real_, nrow = n_perm, ncol = length(testable))
for (i in seq_len(n_perm)) {
  for (j in seq_along(testable)) {
    hidx <- sample.int(N, length(testable[[j]]))
    perm_mat[i, j] <- calc_es(hidx, scores_abs, N)
  }
}
cat("Concluído em", round(difftime(Sys.time(), t0, units = "secs"), 1), "segundos\n\n")

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

leading_edge <- sapply(testable, function(idx) paste(head(ranked_genes[idx], 30), collapse = "/"))
res <- data.frame(PATH = names(testable), Nh = sapply(testable, length),
                   ES = es_obs, NES = nes, pvalue = pval, leadingEdge = leading_edge)
res$p.adjust <- p.adjust(res$pvalue, "BH")
res <- res[order(res$pvalue), ]

dir.create("results", showWarnings = FALSE)
write.csv(res, "results/GSEA_preranked_Wei2020full_KEGG.csv", row.names = FALSE)

cat("=== GSEA preranked do Wei 2020 COMPLETO (SUI, KEGG) ===\n")
cat("Vias significativas a FDR<0.25:", sum(res$p.adjust < 0.25, na.rm = TRUE), "de", nrow(res), "\n")
cat("Vias significativas a FDR<0.05:", sum(res$p.adjust < 0.05, na.rm = TRUE), "\n\n")
print(head(res[, c("PATH","Nh","ES","NES","pvalue","p.adjust")], 15))
