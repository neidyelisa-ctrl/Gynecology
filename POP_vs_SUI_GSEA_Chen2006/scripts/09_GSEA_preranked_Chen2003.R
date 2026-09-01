# =============================================================================
# Script 9 — GSEA preranked no Chen et al. 2003 (SUI, fase proliferativa),
# vias KEGG. Mesmo método do script 03 (Chen 2006), aplicado ao segundo
# artigo. Ver o cabeçalho do script 03 para a explicação metodológica
# completa (painel curto, permutação de rótulo de gene set, etc.) - não
# repetida aqui.
#
# Diferença de formato: o Chen 2003 reporta só 1 p-valor por gene (não
# RMA/MAS5 separados como o Chen 2006), então o escore de ranking usa
# diretamente a coluna `pvalue` do CSV.
# =============================================================================

suppressMessages(library(org.Hs.eg.db))
set.seed(2003)

chen03 <- read.csv("data/chen2003_90genes.csv")
chen03 <- chen03[!duplicated(chen03$GeneSymbol), ]
cat("Chen 2003 — DEGs curados com símbolo HGNC atual:", nrow(chen03),
    "(", sum(chen03$Direction == "up"), "up /", sum(chen03$Direction == "down"), "down )\n\n")

build_score <- function(direction, pvalue) {
  p <- pvalue
  floor_p <- min(p[p > 0], na.rm = TRUE) / 2
  p[p <= 0] <- floor_p
  sign <- ifelse(direction == "up", 1, -1)
  sign * -log10(p)
}

score <- build_score(chen03$Direction, chen03$pvalue)
ord <- order(-score)
ranked_genes <- chen03$GeneSymbol[ord]
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

scores_abs <- abs(ranked_scores)
es_obs <- sapply(testable, calc_es, scores_abs = scores_abs, N = N)

set.seed(2003)
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
write.csv(res, "results/GSEA_preranked_Chen2003_KEGG.csv", row.names = FALSE)

cat("\n=== GSEA preranked do Chen 2003 (SUI, KEGG) ===\n")
cat("Vias significativas a FDR<0.25:", sum(res$p.adjust < 0.25, na.rm = TRUE), "\n")
cat("Vias significativas a FDR<0.05:", sum(res$p.adjust < 0.05, na.rm = TRUE), "\n\n")
print(head(res[, c("PATH","Nh","ES","NES","pvalue","p.adjust","leadingEdge")], 10))
