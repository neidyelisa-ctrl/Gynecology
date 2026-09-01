# =============================================================================
# Script 2/4 — GSEA "clássico" (preranked, transcriptoma inteiro) no POP
# (GSE53868), vias KEGG.
#
# "Clássico" aqui = a forma padrão de GSEA (Subramanian et al. 2005, PNAS):
# TODO o transcriptoma ranqueado por um escore contínuo (estatística-t
# pareada), sem nenhum corte de significância prévio — diferente de ORA
# (que só olha os genes já significativos). Isso é o que diferencia GSEA de
# um teste hipergeométrico simples: capta efeito coordenado e pequeno
# espalhado por muitos genes de uma via, mesmo que nenhum gene individual
# passe no limiar de significância.
#
# IMPLEMENTAÇÃO: este sandbox não tem acesso a CRAN/Bioconductor/MSigDB ao
# vivo (confirmado por teste direto: CONNECT bloqueado para cloud.r-project.org,
# bioconductor.org, zenodo.org, KEGG, Reactome, Enrichr, etc. — política de
# rede do ambiente). Os pacotes padrão (fgsea, clusterProfiler, msigdbr) não
# podem ser instalados aqui. Por isso, o algoritmo do GSEA-Preranked é
# reimplementado do zero (weighted running-sum enrichment score, ES, com
# p=1) usando as vias KEGG já anotadas no pacote offline org.Hs.eg.db
# (Bioconductor, instalado via apt — não precisa de download em tempo de
# execução). O script 05 traz a versão "oficial" (fgsea + msigdbr) para você
# rodar no seu R local, com internet normal, e obter o resultado com a base
# COMPLETA e atualizada do MSigDB (Hallmark/KEGG/Reactome/GO) em vez do KEGG
# congelado do org.Hs.eg.db.
#
# PERMUTAÇÃO: usamos permutação de FENÓTIPO (sign-flip pareado — inverter
# aleatoriamente, por paciente, qual amostra é "POP site" e qual é "non-POP
# site"), não permutação de rótulo de gene set. Isso é importante: genes de
# uma mesma via são correlacionados na expressão real do tecido, e permutar
# só os rótulos do gene set (mantendo o ranking fixo) subestima a variância
# real sob a hipótese nula, inflando falsos positivos — o próprio
# Subramanian et al. 2005 alerta para isso e recomenda permutação de
# fenótipo sempre que houver matriz de expressão bruta por amostra (que
# temos aqui, já que o GSE53868 é um desenho pareado com 12 pacientes).
# =============================================================================

suppressMessages(library(org.Hs.eg.db))
set.seed(42)

## --- 1) Reconstruir a matriz de diferença pareada (mesmos 12 pacientes) ---
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
stopifnot(setequal(pop_ind, nonpop_ind))
pop_cols <- pop_cols[match(nonpop_ind, pop_ind)]  # alinhar por paciente

D <- expr[, pop_cols] - expr[, nonpop_cols]  # genes x 12 pacientes; + = para cima no sitio do prolapso

# descarta sondas com valor faltante ou variância zero entre pacientes
# (t pareado indefinido) — mesmo universo de genes em toda permutação
valid_rows <- stats::complete.cases(D) & (apply(D, 1, sd) > 0)
cat("Descartando", sum(!valid_rows), "sondas com valor faltante ou variância zero\n")
D <- D[valid_rows, , drop = FALSE]
n_pat <- ncol(D)
cat("Matriz de diferença pareada D:", nrow(D), "genes x", n_pat, "pacientes\n")

paired_t <- function(D) {
  m <- rowMeans(D); s <- apply(D, 1, sd); n <- ncol(D)
  m / (s / sqrt(n))
}

t_obs <- paired_t(D)
ord <- order(-t_obs)
ranked_genes <- rownames(D)[ord]
ranked_scores <- t_obs[ord]
N <- length(ranked_genes)
cat("Ranked list (estatística-t pareada, positivo = para cima no sítio do prolapso):", N, "genes\n\n")

## --- 2) Vias KEGG (offline, org.Hs.eg.db), tamanho 5-200 --------------------
ann <- suppressWarnings(select(org.Hs.eg.db, keys = ranked_genes, keytype = "SYMBOL", columns = "PATH"))
ann <- ann[!is.na(ann$PATH), ]
gs_sizes <- table(ann$PATH)
valid_paths <- names(gs_sizes)[gs_sizes >= 5 & gs_sizes <= 200]
gene_sets <- split(ann$SYMBOL[ann$PATH %in% valid_paths], ann$PATH[ann$PATH %in% valid_paths])
cat("Vias KEGG testadas (5-200 genes, org.Hs.eg.db offline):", length(gene_sets), "\n\n")

## --- 3) Enrichment score (weighted, p=1), padrão Subramanian et al. 2005 --
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

## --- 4) Permutação de FENÓTIPO (sign-flip pareado), 500 permutações -------
n_perm <- 500
cat("Rodando", n_perm, "permutações de fenótipo (sign-flip) em", length(hit_idx_list), "vias...\n")
t0 <- Sys.time()

gene_sets_syms <- lapply(hit_idx_list, function(idx) ranked_genes[idx])
perm_es_mat <- matrix(NA_real_, nrow = n_perm, ncol = length(hit_idx_list))
colnames(perm_es_mat) <- names(hit_idx_list)

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

leading_edge <- sapply(hit_idx_list, function(idx) paste(ranked_genes[idx], collapse = "/"))
gsea_res <- data.frame(PATH = names(hit_idx_list), Nh = sapply(hit_idx_list, length),
                        ES = es_obs, NES = nes, pvalue = pval, leadingEdge = leading_edge)
gsea_res$p.adjust <- p.adjust(gsea_res$pvalue, "BH")
gsea_res <- gsea_res[order(gsea_res$pvalue), ]

dir.create("results", showWarnings = FALSE)
write.csv(gsea_res, "results/GSEA_classic_POP_KEGG.csv", row.names = FALSE)

cat("=== GSEA clássico do POP (KEGG, permutação de fenótipo) ===\n")
cat("Vias significativas a FDR<0.25 (limiar de rastreio padrão do GSEA):",
    sum(gsea_res$p.adjust < 0.25, na.rm = TRUE), "de", nrow(gsea_res), "\n")
cat("Vias significativas a FDR<0.05:", sum(gsea_res$p.adjust < 0.05, na.rm = TRUE), "\n\n")
cat("NOTA: com 500 permutações, muitas vias empatam no p-valor mínimo possível\n")
cat("(1/501 ~= 0.001996) — não dá para refinar esse ranking sem mais permutações.\n")
print(head(gsea_res[, c("PATH", "Nh", "ES", "NES", "pvalue", "p.adjust")], 15))
