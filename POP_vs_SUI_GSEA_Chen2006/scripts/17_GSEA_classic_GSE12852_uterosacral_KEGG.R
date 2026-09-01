# =============================================================================
# Script 17 — GSEA "clássico" (transcriptoma inteiro, permutação de
# fenótipo) no GSE12852, ligamento UTEROSSACRAL sozinho (ver diagnóstico no
# script 16 - a análise combinada com o ligamento redondo dilui o sinal;
# aqui usamos só as 17 amostras de ligamento uterossacral, 8 POP x 9
# controle, desenho não-pareado simples).
#
# Mesmo algoritmo do script 02 (POP GSE53868): ES ponderado (Subramanian et
# al. 2005) + permutação de FENÓTIPO. Aqui a permutação é mais simples que
# no GSE53868 (que era pareado, sign-flip): como o GSE12852 uterossacral é
# um desenho caso-controle não-pareado, a permutação de fenótipo é
# reembaralhar aleatoriamente qual das 17 pacientes é "POP" (8) e qual é
# "controle" (9), recalculando o teste moderado (limma) e o ranking inteiro
# a cada permutação - preserva a correlação real entre genes, igual ao
# GSE53868.
# =============================================================================

suppressMessages({
  library(limma)
  library(org.Hs.eg.db)
})
set.seed(12852)

## --- 1) Reconstruir a matriz de expressão (uterossacral, gene-level) ------
raw_lines <- readLines("data/GSE12852_series_matrix.txt")
start_row <- grep("^!series_matrix_table_begin", raw_lines) + 1
end_row   <- grep("^!series_matrix_table_end", raw_lines) - 1
expr <- as.matrix(read.delim("data/GSE12852_series_matrix.txt", skip = start_row - 1,
                              nrows = end_row - start_row, header = TRUE,
                              row.names = 1, check.names = FALSE, quote = "\""))

sample_title_line <- raw_lines[grep("^!Sample_title", raw_lines)]
titles <- gsub('"', "", strsplit(sample_title_line, "\t")[[1]][-1])
tissue <- ifelse(grepl("^round", titles), "round", "uterosacral")
condition <- ifelse(grepl("_POP_", titles), "POP", "control")

idx <- tissue == "uterosacral"
expr_u <- expr[, idx]
cond_u <- factor(condition[idx], levels = c("control", "POP"))
cat("Ligamento uterossacral:", sum(cond_u == "POP"), "POP x", sum(cond_u == "control"), "controle\n")

expr_log2 <- log2(pmax(expr_u, 1))
annot <- read.delim("data/GPL2986_annotation.tsv", colClasses = "character")
annot <- annot[!duplicated(annot$ID), ]; rownames(annot) <- annot$ID
common <- intersect(rownames(expr_log2), annot$ID)
expr_gene <- avereps(expr_log2[common, ], ID = annot[common, "GeneSymbol"])
N_genes <- nrow(expr_gene)
cat("Genes únicos (colapsados por sonda):", N_genes, "\n\n")

paired <- FALSE
get_ranking <- function(expr_mat, cond) {
  design <- model.matrix(~cond)
  fit <- eBayes(lmFit(expr_mat, design))
  fit$t[, 2]  # t da condicao POP (positivo = para cima no POP)
}

t_obs <- get_ranking(expr_gene, cond_u)
ord <- order(-t_obs)
ranked_genes <- names(t_obs)[ord]
ranked_scores <- t_obs[ord]
N <- length(ranked_genes)
cat("Ranked list (moderado-t, positivo = para cima no POP):", N, "genes\n\n")

## --- 2) Vias KEGG (offline, org.Hs.eg.db), tamanho 5-200 --------------------
ann <- suppressWarnings(select(org.Hs.eg.db, keys = ranked_genes, keytype = "SYMBOL", columns = "PATH"))
ann <- ann[!is.na(ann$PATH), ]
gs_sizes <- table(ann$PATH)
valid_paths <- names(gs_sizes)[gs_sizes >= 5 & gs_sizes <= 200]
gene_sets <- split(ann$SYMBOL[ann$PATH %in% valid_paths], ann$PATH[ann$PATH %in% valid_paths])
cat("Vias KEGG testadas (5-200 genes):", length(gene_sets), "\n\n")

calc_es <- function(hit_idx, scores_abs, N) {
  Nh <- length(hit_idx); Nm <- N - Nh
  sum_hit <- sum(scores_abs[hit_idx])
  step <- rep(-1 / Nm, N)
  step[hit_idx] <- scores_abs[hit_idx] / sum_hit
  running <- cumsum(step)
  running[which.max(abs(running))]
}

hit_idx_list <- lapply(gene_sets, function(g) which(ranked_genes %in% g))
hit_idx_list <- hit_idx_list[sapply(hit_idx_list, length) >= 3]
cat("Vias com >=3 genes na ranked list:", length(hit_idx_list), "\n")

es_obs <- sapply(hit_idx_list, calc_es, scores_abs = abs(ranked_scores), N = N)

## --- 3) Permutação de FENÓTIPO: reembaralhar POP/controle entre as 17 -----
##        pacientes, recalcular limma + ranking inteiro a cada vez
n_perm <- 500
cat("Rodando", n_perm, "permutações de fenótipo em", length(hit_idx_list), "vias...\n")
t0 <- Sys.time()
n_pop <- sum(cond_u == "POP"); n_total <- length(cond_u)

perm_es_mat <- matrix(NA_real_, nrow = n_perm, ncol = length(hit_idx_list))
gene_sets_syms <- lapply(hit_idx_list, function(idx2) ranked_genes[idx2])

for (i in seq_len(n_perm)) {
  perm_labels <- factor(rep("control", n_total), levels = c("control", "POP"))
  perm_labels[sample.int(n_total, n_pop)] <- "POP"
  t_perm <- get_ranking(expr_gene, perm_labels)
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
    pval[j] <- (sum(pe >= es_obs[j]) + 1) / (length(pe) + 1)
    base <- mean(pe[pe >= 0]); if (is.nan(base) || base == 0) base <- NA
  } else {
    pval[j] <- (sum(pe <= es_obs[j]) + 1) / (length(pe) + 1)
    base <- mean(abs(pe[pe < 0])); if (is.nan(base) || base == 0) base <- NA
  }
  nes[j] <- es_obs[j] / base
}

leading_edge <- sapply(hit_idx_list, function(idx2) paste(ranked_genes[idx2], collapse = "/"))
gsea_res <- data.frame(PATH = names(hit_idx_list), Nh = sapply(hit_idx_list, length),
                        ES = es_obs, NES = nes, pvalue = pval, leadingEdge = leading_edge)
gsea_res$p.adjust <- p.adjust(gsea_res$pvalue, "BH")
gsea_res <- gsea_res[order(gsea_res$pvalue), ]

dir.create("results", showWarnings = FALSE)
write.csv(gsea_res, "results/GSEA_classic_GSE12852_uterosacral_KEGG.csv", row.names = FALSE)

cat("=== GSEA clássico GSE12852 (uterossacral, KEGG, permutação de fenótipo) ===\n")
cat("Vias significativas FDR<0.25:", sum(gsea_res$p.adjust < 0.25, na.rm = TRUE), "de", nrow(gsea_res), "\n")
cat("Vias significativas FDR<0.05:", sum(gsea_res$p.adjust < 0.05, na.rm = TRUE), "\n\n")
print(head(gsea_res[, c("PATH","Nh","ES","NES","pvalue","p.adjust")], 15))
