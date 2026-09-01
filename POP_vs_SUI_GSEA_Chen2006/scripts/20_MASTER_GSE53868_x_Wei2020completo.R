# =============================================================================
# Script 20 (MASTER) — GSE53868 (POP, parede vaginal) x Wei 2020 completo
# (SUI, parede vaginal periuretral) — a comparação PRIORITÁRIA, porque os
# dois tecidos têm a MESMA origem anatômica (diferente do GSE12852, que é
# ligamento uterossacral/redondo). Este script junta num arquivo só o que
# estava espalhado nos scripts 01, 02, 14 e 15 - para rodar do início ao
# fim e conferir os números principais desta comparação.
#
# CRITÉRIOS USADOS (respondendo diretamente a essa pergunta):
#   - DEG do POP: |log2FC| > 1 E FDR (Benjamini-Hochberg) < 0,05
#     (o critério pedido desde o início, usado em TODOS os DEGs deste
#     projeto - GSE53868, GSE12852)
#   - Painel do Wei 2020: os próprios critérios do artigo original
#     (fold change >= 2 E P < 0,05, SEM correção de FDR - é assim que o
#     GeneSpring GX do estudo original define a lista de 7.102 genes;
#     nós NÃO aplicamos um filtro adicional de FDR/logFC>1 em cima disso)
#   - Teste de CONCORDÂNCIA DE DIREÇÃO (binomial): usa TODOS os genes
#     testáveis nos dois lados, SEM nenhum filtro de significância ou
#     magnitude - só o SINAL (para cima/para baixo). Isso é proposital:
#     é um teste de padrão agregado, não de genes individuais - inclui
#     genes fracos de propósito, porque o objetivo é detectar um viés
#     sistemático de direção, não achar genes "campeões".
#   - GSEA (preranked e clássico): sem nenhum filtro de significância -
#     usa o ranking do transcriptoma INTEIRO (ou do painel inteiro do
#     Wei2020), como o método exige.
#
# TEMPO ESTIMADO: ~10 minutos (a maior parte é a permutação de fenótipo
# do GSEA clássico do POP, ~9 min sozinha).
# =============================================================================

suppressMessages({
  library(limma)
  library(org.Hs.eg.db)
  library(GO.db)
})

dir.create("results", showWarnings = FALSE)

## =============================================================================
## PARTE 1 — DEG do POP (GSE53868), limma pareado
## =============================================================================
cat("\n############## PARTE 1: DEG do POP (GSE53868) ##############\n\n")

raw_lines <- readLines("data/GSE53868_series_matrix.txt")
start_row <- grep("^!series_matrix_table_begin", raw_lines) + 1
end_row   <- grep("^!series_matrix_table_end", raw_lines) - 1
expr <- read.delim("data/GSE53868_series_matrix.txt", skip = start_row - 1,
                    nrows = end_row - start_row, header = TRUE,
                    row.names = 1, check.names = FALSE, quote = "\"")

sample_title_line <- raw_lines[grep("^!Sample_title", raw_lines)]
sample_titles <- gsub('"', "", strsplit(sample_title_line, "\t")[[1]][-1])
individual_line <- raw_lines[grep("^!Sample_characteristics_ch1.*individual:", raw_lines)][1]
individuals <- gsub("individual: ", "", gsub('"', "", strsplit(individual_line, "\t")[[1]][-1]))
tissue <- ifelse(grepl("\\(POP site\\)", sample_titles), "POP_site", "NonPOP_site")

coldata <- data.frame(row.names = colnames(expr),
                       tissue = factor(tissue, levels = c("NonPOP_site", "POP_site")),
                       individual = factor(individuals))

design <- model.matrix(~ individual + tissue, data = coldata)
fit <- eBayes(lmFit(as.matrix(expr), design))
pop_full <- topTable(fit, coef = "tissuePOP_site", number = Inf, sort.by = "P")
pop_full$Gene <- rownames(pop_full)
pop_full <- pop_full[, c("Gene", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val")]
write.csv(pop_full, "results/GSE53868_limma_completo.csv", row.names = FALSE)

pop_deg <- subset(pop_full, adj.P.Val < 0.05 & abs(logFC) > 1)
write.csv(pop_deg, "results/GSE53868_DEG_logFC1_FDR05.csv", row.names = FALSE)
cat("DEGs do POP (|log2FC|>1, FDR<0.05):", nrow(pop_deg), "de", nrow(pop_full), "testados\n")

## =============================================================================
## PARTE 2 — GSEA CLÁSSICO do POP (permutação de FENÓTIPO, sign-flip pareado)
## =============================================================================
cat("\n############## PARTE 2: GSEA classico do POP (KEGG) ##############\n\n")

expr_m <- as.matrix(expr)
pop_cols <- colnames(expr_m)[tissue == "POP_site"]
nonpop_cols <- colnames(expr_m)[tissue == "NonPOP_site"]
pop_ind <- individuals[tissue == "POP_site"]; nonpop_ind <- individuals[tissue == "NonPOP_site"]
pop_cols <- pop_cols[match(nonpop_ind, pop_ind)]
D <- expr_m[, pop_cols] - expr_m[, nonpop_cols]
valid_rows <- stats::complete.cases(D) & (apply(D, 1, sd) > 0)
D <- D[valid_rows, , drop = FALSE]
n_pat <- ncol(D)

paired_t <- function(D) { m <- rowMeans(D); s <- apply(D, 1, sd); n <- ncol(D); m / (s / sqrt(n)) }
t_obs <- paired_t(D)
ord <- order(-t_obs)
ranked_genes <- rownames(D)[ord]; ranked_scores <- t_obs[ord]; N <- length(ranked_genes)

ann <- suppressWarnings(select(org.Hs.eg.db, keys = ranked_genes, keytype = "SYMBOL", columns = "PATH"))
ann <- ann[!is.na(ann$PATH), ]
gs_sizes <- table(ann$PATH)
valid_paths <- names(gs_sizes)[gs_sizes >= 5 & gs_sizes <= 200]
gene_sets <- split(ann$SYMBOL[ann$PATH %in% valid_paths], ann$PATH[ann$PATH %in% valid_paths])

calc_es <- function(hit_idx, scores_abs, N) {
  Nh <- length(hit_idx); Nm <- N - Nh
  sum_hit <- sum(scores_abs[hit_idx])
  step <- rep(-1 / Nm, N); step[hit_idx] <- scores_abs[hit_idx] / sum_hit
  running <- cumsum(step); running[which.max(abs(running))]
}
hit_idx_list <- lapply(gene_sets, function(g) which(ranked_genes %in% g))
hit_idx_list <- hit_idx_list[sapply(hit_idx_list, length) >= 3]
es_obs <- sapply(hit_idx_list, calc_es, scores_abs = abs(ranked_scores), N = N)

n_perm <- 500
cat("Rodando", n_perm, "permutacoes de fenotipo (~9 min)...\n")
set.seed(42)
gene_sets_syms <- lapply(hit_idx_list, function(idx) ranked_genes[idx])
perm_es_mat <- matrix(NA_real_, nrow = n_perm, ncol = length(hit_idx_list))
for (i in seq_len(n_perm)) {
  signs <- sample(c(-1, 1), n_pat, replace = TRUE)
  D_perm <- sweep(D, 2, signs, `*`)
  t_perm <- paired_t(D_perm)
  rank_of_gene <- rank(-t_perm, ties.method = "first")
  scores_abs_sorted <- sort(abs(t_perm), decreasing = TRUE)
  for (j in seq_along(gene_sets_syms)) {
    hidx <- rank_of_gene[gene_sets_syms[[j]]]
    if (length(hidx) >= 3) perm_es_mat[i, j] <- calc_es(hidx, scores_abs_sorted, N)
  }
}

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
leading_edge <- sapply(hit_idx_list, function(idx) paste(ranked_genes[idx], collapse = "/"))
pop_gsea <- data.frame(PATH = names(hit_idx_list), Nh = sapply(hit_idx_list, length),
                        ES = es_obs, NES = nes, pvalue = pval, leadingEdge = leading_edge)
pop_gsea$p.adjust <- p.adjust(pop_gsea$pvalue, "BH")
pop_gsea <- pop_gsea[order(pop_gsea$pvalue), ]
write.csv(pop_gsea, "results/GSEA_classic_POP_KEGG.csv", row.names = FALSE)
cat("POP GSEA - sig FDR<0.05:", sum(pop_gsea$p.adjust < 0.05, na.rm = TRUE),
    "| FDR<0.25:", sum(pop_gsea$p.adjust < 0.25, na.rm = TRUE), "de", nrow(pop_gsea), "\n")

## =============================================================================
## PARTE 3 — Wei 2020 completo: preranked GSEA (KEGG)
## =============================================================================
cat("\n############## PARTE 3: GSEA preranked do Wei 2020 (KEGG) ##############\n\n")

wei <- read.csv("data/wei2020_mRNA_full.csv")
build_score <- function(direction, pvalue) {
  p <- pvalue; floor_p <- min(p[p > 0], na.rm = TRUE) / 2; p[p <= 0] <- floor_p
  sign <- ifelse(direction == "up", 1, -1); sign * -log10(p)
}
score <- build_score(wei$Direction, wei$PValue)
ordw <- order(-score)
ranked_genes_w <- wei$GeneSymbol[ordw]; ranked_scores_w <- score[ordw]; Nw <- length(ranked_genes_w)

all_sym <- keys(org.Hs.eg.db, keytype = "SYMBOL")
ann_path <- suppressWarnings(select(org.Hs.eg.db, keys = all_sym, keytype = "SYMBOL", columns = "PATH"))
ann_path <- ann_path[!is.na(ann_path$PATH), ]
kegg_sets_all <- split(ann_path$SYMBOL, ann_path$PATH)

min_gs <- 5; n_perm_w <- 1000
testable <- lapply(kegg_sets_all, function(g) which(ranked_genes_w %in% g))
testable <- testable[sapply(testable, length) >= min_gs]
scores_abs_w <- abs(ranked_scores_w)
es_obs_w <- sapply(testable, calc_es, scores_abs = scores_abs_w, N = Nw)

set.seed(2020)
perm_mat_w <- matrix(NA_real_, nrow = n_perm_w, ncol = length(testable))
for (i in seq_len(n_perm_w)) {
  for (j in seq_along(testable)) {
    hidx <- sample.int(Nw, length(testable[[j]]))
    perm_mat_w[i, j] <- calc_es(hidx, scores_abs_w, Nw)
  }
}
pval_w <- numeric(length(testable)); nes_w <- numeric(length(testable))
for (j in seq_along(testable)) {
  pe <- perm_mat_w[, j]; pe <- pe[!is.na(pe)]
  if (is.na(es_obs_w[j])) { pval_w[j] <- NA; nes_w[j] <- NA; next }
  if (es_obs_w[j] >= 0) {
    pval_w[j] <- (sum(pe >= es_obs_w[j]) + 1) / (length(pe) + 1)
    base <- mean(pe[pe >= 0]); if (is.nan(base) || base == 0) base <- NA
  } else {
    pval_w[j] <- (sum(pe <= es_obs_w[j]) + 1) / (length(pe) + 1)
    base <- mean(abs(pe[pe < 0])); if (is.nan(base) || base == 0) base <- NA
  }
  nes_w[j] <- es_obs_w[j] / base
}
leading_edge_w <- sapply(testable, function(idx) paste(head(ranked_genes_w[idx], 30), collapse = "/"))
wei_gsea <- data.frame(PATH = names(testable), Nh = sapply(testable, length),
                        ES = es_obs_w, NES = nes_w, pvalue = pval_w, leadingEdge = leading_edge_w)
wei_gsea$p.adjust <- p.adjust(wei_gsea$pvalue, "BH")
wei_gsea <- wei_gsea[order(wei_gsea$pvalue), ]
write.csv(wei_gsea, "results/GSEA_preranked_Wei2020full_KEGG.csv", row.names = FALSE)
cat("Wei2020 GSEA - sig FDR<0.05:", sum(wei_gsea$p.adjust < 0.05, na.rm = TRUE),
    "| FDR<0.25:", sum(wei_gsea$p.adjust < 0.25, na.rm = TRUE), "de", nrow(wei_gsea), "\n")

## =============================================================================
## PARTE 4 — Cruzamento gene a gene (concordância de direção, SEM filtro)
## =============================================================================
cat("\n############## PARTE 4: Concordancia de direcao (genes) ##############\n\n")

cross <- merge(wei[, c("GeneSymbol", "Direction")], pop_full, by.x = "GeneSymbol", by.y = "Gene", all.x = TRUE)
cross$Testado_no_POP <- !is.na(cross$logFC)
cross$Direcao_POP <- ifelse(is.na(cross$logFC), NA, ifelse(cross$logFC > 0, "up", "down"))
cross$Concordante <- cross$Direction == cross$Direcao_POP
cross$DEG_POP_FDR05_logFC1 <- !is.na(cross$adj.P.Val) & cross$adj.P.Val < 0.05 & abs(cross$logFC) > 1
write.csv(cross, "results/20_cruzamento_genes_GSE53868_x_Wei2020.csv", row.names = FALSE)

testaveis <- subset(cross, Testado_no_POP & !is.na(Direcao_POP))
concordantes <- subset(testaveis, Concordante)
bt <- binom.test(nrow(concordantes), nrow(testaveis), p = 0.5)
cat("Interseccao estrita (DEG POP ^ Wei2020):", sum(cross$DEG_POP_FDR05_logFC1, na.rm=TRUE), "genes\n")
cat("Concordancia (TODOS os genes testaveis, sem filtro de significancia):",
    nrow(concordantes), "de", nrow(testaveis), "=",
    round(100 * nrow(concordantes) / nrow(testaveis), 1), "% | p =", format(bt$p.value, digits = 4), "\n")

## =============================================================================
## PARTE 5 — Vias compartilhadas (GSEA x GSEA), com checagem de DIRECAO
## =============================================================================
cat("\n############## PARTE 5: Vias compartilhadas + direcao ##############\n\n")

pop_gsea$PATH <- sprintf("%05d", as.integer(pop_gsea$PATH))
wei_gsea$PATH <- sprintf("%05d", as.integer(wei_gsea$PATH))

for (fdr in c(0.25, 0.05)) {
  pop_sig <- subset(pop_gsea, p.adjust < fdr)
  wei_sig <- subset(wei_gsea, p.adjust < fdr)
  common <- intersect(pop_sig$PATH, wei_sig$PATH)
  cat("--- FDR<", fdr, ": POP sig =", nrow(pop_sig), " | Wei2020 sig =", nrow(wei_sig),
      " | COMPARTILHADAS =", length(common), " ---\n", sep = "")
  if (length(common) > 0) {
    out <- merge(pop_sig[pop_sig$PATH %in% common, c("PATH","NES","p.adjust")],
                 wei_sig[wei_sig$PATH %in% common, c("PATH","NES","p.adjust")],
                 by = "PATH", suffixes = c("_POP", "_Wei2020"))
    out$Mesma_direcao <- sign(out$NES_POP) == sign(out$NES_Wei2020)
    print(out)
    write.csv(out, sprintf("results/20_shared_pathways_FDR%03d.csv", fdr * 100), row.names = FALSE)
    cat("Concordantes:", sum(out$Mesma_direcao, na.rm = TRUE), "de", sum(!is.na(out$Mesma_direcao)), "\n")
  }
  cat("\n")
}

cat("=== FIM (script 20) ===\n")
