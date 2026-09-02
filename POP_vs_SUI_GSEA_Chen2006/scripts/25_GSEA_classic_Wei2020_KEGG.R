# =============================================================================
# Script 25 — GSEA CLÁSSICO (não preranked) do SUI (Wei 2020), com
# permutação de FENÓTIPO real — agora possível porque temos a intensidade
# por amostra (3 SUI x 3 Ctrl) extraída da Tabela S2
# (`data/wei2020_persample_normalized.csv`, script 24). Antes (scripts
# 03/09/13/15) só dava pra fazer preranked com permutação de rótulo de
# gene set, porque só tínhamos p-valor/fold-change prontos, não a matriz
# de expressão bruta por amostra.
#
# Mesmo algoritmo dos scripts 02/17 (POP): ES ponderado (Subramanian et
# al. 2005) + permutação de FENÓTIPO, ranking = t moderado (limma) do
# transcriptoma inteiro (aqui, os 6.118 genes da Tabela S2 - ainda não é
# o array inteiro do estudo original, mas já não é mais um "painel
# pequeno" de 40-90 genes escolhidos a dedo).
#
# ATENÇÃO - LIMITAÇÃO ESTATÍSTICA IMPORTANTE DO DESENHO 3 x 3 (leia antes
# de interpretar os p-valores): com só 3 SUI e 3 Ctrl, existem exatamente
# choose(6,3) = 20 formas possíveis de dividir as 6 amostras em dois
# grupos de 3 - ou seja, só 20 permutações DISTINTAS existem no total
# (a observada + 19 outras). Rodar 1000 permutações aleatórias não dá mais
# resolução real - é amostragem repetida das mesmas 20 possibilidades. Por
# isso este script faz permutação EXAUSTIVA (as 20 combinações exatas, não
# uma amostra aleatória) - é o correto estatisticamente para n tão pequeno,
# e deixa claro o teto de resolução: o p-valor mínimo possível aqui é
# 1/20 = 0,05 (ou 2/20=0,10 dependendo de como a via cai nas duas caudas) -
# ou seja, NENHUMA via pode, em princípio, ficar abaixo de p=0,05 nesta
# permutação, e o FDR (BH) sobre p-valores discretos e pouco granulares
# tende a ficar conservador. Isto não é um bug do script - é uma
# consequência inescapável de só termos 3 réplicas por grupo neste dataset.
# =============================================================================

suppressMessages({
  library(limma)
  library(org.Hs.eg.db)
})

## --- 1) Ler a matriz de intensidade por amostra e colapsar por gene -------
wei_raw <- read.csv("data/wei2020_persample_normalized.csv")
wei_mat <- as.matrix(wei_raw[, c("Sui1","Sui2","Sui3","Ctrl1","Ctrl2","Ctrl3")])
rownames(wei_mat) <- wei_raw$GeneSymbol
wei_gene <- avereps(wei_mat, ID = rownames(wei_mat))
cat("Genes únicos (colapsados por sonda):", nrow(wei_gene), "x", ncol(wei_gene), "amostras (3 SUI + 3 Ctrl)\n\n")

group <- factor(c("SUI","SUI","SUI","Ctrl","Ctrl","Ctrl"), levels = c("Ctrl","SUI"))

get_ranking <- function(expr_mat, grp) {
  design <- model.matrix(~grp)
  fit <- eBayes(lmFit(expr_mat, design))
  fit$t[, 2]  # t da condicao SUI (positivo = para cima no SUI)
}

t_obs <- get_ranking(wei_gene, group)
ord <- order(-t_obs)
ranked_genes <- names(t_obs)[ord]
ranked_scores <- t_obs[ord]
N <- length(ranked_genes)
cat("Ranked list (moderado-t, positivo = para cima no SUI):", N, "genes\n\n")

## --- 2) Vias KEGG (offline, org.Hs.eg.db) -----------------------------------
ann <- suppressWarnings(select(org.Hs.eg.db, keys = ranked_genes, keytype = "SYMBOL", columns = "PATH"))
ann <- ann[!is.na(ann$PATH), ]
gs_sizes <- table(ann$PATH)
valid_paths <- names(gs_sizes)[gs_sizes >= 5 & gs_sizes <= 200]
gene_sets <- split(ann$SYMBOL[ann$PATH %in% valid_paths], ann$PATH[ann$PATH %in% valid_paths])
cat("Vias KEGG testadas (5-200 genes):", length(gene_sets), "\n\n")

calc_es <- function(hit_idx, scores_abs, N) {
  Nh <- length(hit_idx); Nm <- N - Nh
  sum_hit <- sum(scores_abs[hit_idx])
  step <- rep(-1 / Nm, N); step[hit_idx] <- scores_abs[hit_idx] / sum_hit
  running <- cumsum(step); running[which.max(abs(running))]
}
hit_idx_list <- lapply(gene_sets, function(g) which(ranked_genes %in% g))
hit_idx_list <- hit_idx_list[sapply(hit_idx_list, length) >= 3]
cat("Vias com >=3 genes na ranked list:", length(hit_idx_list), "\n\n")
es_obs <- sapply(hit_idx_list, calc_es, scores_abs = abs(ranked_scores), N = N)

## --- 3) Permutação EXAUSTIVA de fenótipo (as 20 combinações possíveis) ----
all_perms <- combn(6, 3)  # 20 colunas, cada uma = indices que viram "SUI"
n_perm <- ncol(all_perms)
cat("Permutação EXAUSTIVA:", n_perm, "combinações possíveis (todas testadas)\n")
t0 <- Sys.time()

gene_sets_syms <- lapply(hit_idx_list, function(idx) ranked_genes[idx])
perm_es_mat <- matrix(NA_real_, nrow = n_perm, ncol = length(hit_idx_list))

for (i in seq_len(n_perm)) {
  perm_group <- factor(rep("Ctrl", 6), levels = c("Ctrl", "SUI"))
  perm_group[all_perms[, i]] <- "SUI"
  t_perm <- get_ranking(wei_gene, perm_group)
  rank_of_gene <- rank(-t_perm, ties.method = "first")
  scores_abs_sorted <- sort(abs(t_perm), decreasing = TRUE)
  for (j in seq_along(gene_sets_syms)) {
    hidx <- rank_of_gene[gene_sets_syms[[j]]]
    if (length(hidx) >= 3) perm_es_mat[i, j] <- calc_es(hidx, scores_abs_sorted, N)
  }
}
cat("Concluído em", round(difftime(Sys.time(), t0, units = "secs"), 1), "segundos\n\n")

pval <- numeric(length(hit_idx_list)); nes <- numeric(length(hit_idx_list))
for (j in seq_along(hit_idx_list)) {
  pe <- perm_es_mat[, j]; pe <- pe[!is.na(pe)]
  if (es_obs[j] >= 0) {
    pval[j] <- sum(pe >= es_obs[j]) / length(pe)  # exaustivo: sem +1 artificial, a observada esta nas 20
    base <- mean(pe[pe >= 0]); if (is.nan(base) || base == 0) base <- NA
  } else {
    pval[j] <- sum(pe <= es_obs[j]) / length(pe)
    base <- mean(abs(pe[pe < 0])); if (is.nan(base) || base == 0) base <- NA
  }
  pval[j] <- max(pval[j], 1 / n_perm)  # piso: nao pode ser 0 com permutacao exaustiva
  nes[j] <- es_obs[j] / base
}

leading_edge <- sapply(hit_idx_list, function(idx) paste(ranked_genes[idx], collapse = "/"))
gsea_res <- data.frame(PATH = names(hit_idx_list), Nh = sapply(hit_idx_list, length),
                        ES = es_obs, NES = nes, pvalue = pval, leadingEdge = leading_edge)
gsea_res$p.adjust <- p.adjust(gsea_res$pvalue, "BH")
gsea_res <- gsea_res[order(gsea_res$pvalue), ]

dir.create("results", showWarnings = FALSE)
write.csv(gsea_res, "results/GSEA_classic_Wei2020_KEGG.csv", row.names = FALSE)

cat("=== GSEA CLÁSSICO do Wei2020 (SUI, KEGG, permutação exaustiva de fenótipo) ===\n")
cat("Vias significativas FDR<0.25:", sum(gsea_res$p.adjust < 0.25, na.rm = TRUE), "de", nrow(gsea_res), "\n")
cat("Vias significativas FDR<0.05:", sum(gsea_res$p.adjust < 0.05, na.rm = TRUE), "\n")
cat("(lembrete: p-valor minimo possivel = 1/20 = 0.05 - ver nota no cabecalho)\n\n")
print(head(gsea_res[, c("PATH","Nh","ES","NES","pvalue","p.adjust")], 15))
