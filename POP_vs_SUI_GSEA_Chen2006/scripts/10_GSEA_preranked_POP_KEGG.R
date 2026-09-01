# =============================================================================
# Script 10 — GSEA PRERANKED do POP (mesmo ranking do script 02: t pareado,
# transcriptoma inteiro), mas com PERMUTAÇÃO DE RÓTULO DE GENE SET, para
# comparar diretamente com o "GSEA normal/clássico" do script 02 (que usa
# permutação de FENÓTIPO). É o mesmo par de métodos já usado no Chen
# 2003/2006 (script 03/09), aplicado agora ao POP - que tem transcriptoma
# inteiro E matriz bruta por amostra, então dá pra rodar os dois métodos de
# verdade no mesmo dataset e comparar.
#
# POR QUE ISSO IMPORTA (é exatamente o ponto que motivou a correção
# documentada no script 02): a primeira versão deste projeto usou permutação
# de rótulo de gene set no POP e achou 140 de 218 vias (64%) "significativas"
# a FDR<0.25 - um número implausivelmente alto. A causa: genes de uma mesma
# via são correlacionados na expressão real do tecido; embaralhar só os
# RÓTULOS (mantendo o ranking fixo) ignora essa correlação e subestima a
# variância real sob H0, inflando falsos positivos (alerta explícito do
# próprio Subramanian et al. 2005). A permutação de FENÓTIPO (script 02,
# sign-flip pareado) re-ranqueia o transcriptoma inteiro a cada permutação,
# preservando a correlação real entre genes - é o método correto quando se
# tem a matriz bruta por amostra (que temos aqui).
#
# Rodar os dois no MESMO dataset (POP) deixa a diferença visível de forma
# direta e quantificada, não só citada - útil para discutir na tese por que
# a permutação de rótulo de gene set (usada nos scripts 03/09 para o Chen,
# por FALTA de matriz bruta) é uma limitação conhecida, não um capricho
# metodológico.
# =============================================================================

suppressMessages(library(org.Hs.eg.db))
set.seed(10)

## --- 1) Mesmo ranking do script 02 (t pareado, transcriptoma inteiro) -----
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

expr <- as.matrix(expr)
pop_cols <- colnames(expr)[tissue == "POP_site"]
nonpop_cols <- colnames(expr)[tissue == "NonPOP_site"]
pop_ind <- individuals[tissue == "POP_site"]
nonpop_ind <- individuals[tissue == "NonPOP_site"]
pop_cols <- pop_cols[match(nonpop_ind, pop_ind)]

D <- expr[, pop_cols] - expr[, nonpop_cols]
valid_rows <- stats::complete.cases(D) & (apply(D, 1, sd) > 0)
D <- D[valid_rows, , drop = FALSE]

paired_t <- function(D) { m <- rowMeans(D); s <- apply(D, 1, sd); n <- ncol(D); m / (s / sqrt(n)) }
t_obs <- paired_t(D)
ord <- order(-t_obs)
ranked_genes <- rownames(D)[ord]
ranked_scores <- t_obs[ord]
N <- length(ranked_genes)
cat("Ranked list (idêntico ao script 02):", N, "genes\n\n")

## --- 2) Mesmas vias KEGG (5-200 genes) do script 02 ------------------------
ann <- suppressWarnings(select(org.Hs.eg.db, keys = ranked_genes, keytype = "SYMBOL", columns = "PATH"))
ann <- ann[!is.na(ann$PATH), ]
gs_sizes <- table(ann$PATH)
valid_paths <- names(gs_sizes)[gs_sizes >= 5 & gs_sizes <= 200]
gene_sets <- split(ann$SYMBOL[ann$PATH %in% valid_paths], ann$PATH[ann$PATH %in% valid_paths])

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
cat("Vias KEGG testadas (mesmas 218 do script 02):", length(hit_idx_list), "\n\n")

scores_abs <- abs(ranked_scores)
es_obs <- sapply(hit_idx_list, calc_es, scores_abs = scores_abs, N = N)

## --- 3) PERMUTAÇÃO DE RÓTULO DE GENE SET (ranking fica FIXO) ---------------
n_perm <- 1000
cat("Rodando", n_perm, "permutações de RÓTULO DE GENE SET em", length(hit_idx_list), "vias...\n")
t0 <- Sys.time()
perm_mat <- matrix(NA_real_, nrow = n_perm, ncol = length(hit_idx_list))
for (i in seq_len(n_perm)) {
  for (j in seq_along(hit_idx_list)) {
    hidx <- sample.int(N, length(hit_idx_list[[j]]))
    perm_mat[i, j] <- calc_es(hidx, scores_abs, N)
  }
}
cat("Concluído em", round(difftime(Sys.time(), t0, units = "secs"), 1), "segundos\n\n")

pval <- numeric(length(hit_idx_list)); nes <- numeric(length(hit_idx_list))
for (j in seq_along(hit_idx_list)) {
  pe <- perm_mat[, j]; pe <- pe[!is.na(pe)]
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
gsea_res <- data.frame(PATH = names(hit_idx_list), Nh = sapply(hit_idx_list, length),
                        ES = es_obs, NES = nes, pvalue = pval, leadingEdge = leading_edge)
gsea_res$p.adjust <- p.adjust(gsea_res$pvalue, "BH")
gsea_res <- gsea_res[order(gsea_res$pvalue), ]

dir.create("results", showWarnings = FALSE)
write.csv(gsea_res, "results/GSEA_preranked_POP_KEGG_genesetpermutation.csv", row.names = FALSE)

cat("=== GSEA PRERANKED do POP (permutação de rótulo de gene set) ===\n")
cat("Vias significativas a FDR<0.25:", sum(gsea_res$p.adjust < 0.25, na.rm = TRUE), "de", nrow(gsea_res), "\n")
cat("Vias significativas a FDR<0.05:", sum(gsea_res$p.adjust < 0.05, na.rm = TRUE), "\n\n")

## --- 4) Comparação direta com o GSEA "normal" (fenótipo) do script 02 -----
classic <- read.csv("results/GSEA_classic_POP_KEGG.csv")
classic$PATH <- sprintf("%05d", as.integer(classic$PATH))
gsea_res$PATH <- sprintf("%05d", as.integer(gsea_res$PATH))

comp <- merge(classic[, c("PATH","Nh","NES","pvalue","p.adjust")],
              gsea_res[, c("PATH","NES","pvalue","p.adjust")],
              by = "PATH", suffixes = c("_normal_fenotipo", "_preranked_genesetperm"))
comp$Direcao_normal <- ifelse(comp$NES_normal_fenotipo > 0, "para cima no POP", "para baixo no POP")
comp$Direcao_preranked <- ifelse(comp$NES_preranked_genesetperm > 0, "para cima no POP", "para baixo no POP")
comp$Direcao_concorda <- sign(comp$NES_normal_fenotipo) == sign(comp$NES_preranked_genesetperm)
comp$Sig_normal_FDR05 <- comp$p.adjust_normal_fenotipo < 0.05
comp$Sig_preranked_FDR05 <- comp$p.adjust_preranked_genesetperm < 0.05
comp$Sig_so_no_preranked <- comp$Sig_preranked_FDR05 & !comp$Sig_normal_FDR05
comp <- comp[order(comp$pvalue_normal_fenotipo), ]
write.csv(comp, "results/10_comparacao_metodos_POP_normal_vs_preranked.csv", row.names = FALSE)

cat("=== Comparação: GSEA normal (fenótipo) vs. GSEA preranked (rótulo de gene set), mesmo dataset POP ===\n")
cat("Vias significativas (FDR<0.05) - normal (fenótipo):", sum(comp$Sig_normal_FDR05), "\n")
cat("Vias significativas (FDR<0.05) - preranked (rótulo de gene set):", sum(comp$Sig_preranked_FDR05), "\n")
cat("Significativas nos DOIS métodos:", sum(comp$Sig_normal_FDR05 & comp$Sig_preranked_FDR05), "\n")
cat("Significativas SÓ no preranked (rótulo de gene set) -- provável falso positivo por permutação mais fraca:",
    sum(comp$Sig_so_no_preranked), "\n")
either_sig <- comp$Sig_normal_FDR05 | comp$Sig_preranked_FDR05
concord_known <- comp$Direcao_concorda[either_sig]
cat("Concordância de direção (mesmo sinal de NES) nas vias significativas em pelo menos um método:",
    sum(concord_known, na.rm = TRUE), "de", sum(!is.na(concord_known)),
    "(", sum(is.na(concord_known)), "com NES=NA em algum lado - caso de borda da normalização, ver README)\n")

cat("\n=== FIM ===\n")
