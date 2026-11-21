# =============================================================================
# GSEA (KEGG) para a tabela de microarray SUI_vs_Ctrl fornecida pela usuaria
# (43032_2020_144_MOESM2_ESM.xls, material suplementar de um artigo de 2020,
# nao identificado por titulo/DOI aqui - so o nome do arquivo supplementar),
# e as vias compartilhadas com o POP.
#
# Mesma metodologia ja estabelecida para os 3 artigos de literatura SUI
# (Chen 2003, Chen 2006, Tong 2010 - ver R/analysis_GSEA_preranked_3papers_vs_POP.R):
#   - A tabela da usuaria so tem os genes JA FILTRADOS como significativos
#     pelos autores originais (Fold Change>=2.0, p<0.05) - nao o background
#     completo do array. GSEA preranked numa lista assim e um uso
#     nao-padrao do metodo (mesma ressalva metodologica de antes), mas aqui
#     a lista e MUITO maior (4615 up + 2487 down = 7102 genes unicos, vs
#     60-90 genes dos 3 artigos anteriores) - poder estatistico bem maior,
#     muito mais vias devem ser testaveis.
#   - Score de ranking: sign(Regulation) * -log10(p-value) - igual aos 3
#     artigos anteriores, para manter a comparacao consistente (nao log
#     fold-change, que nao e comparavel entre plataformas).
#   - Sem matriz bruta por amostra PAREADA disponivel para permutacao de
#     fenotipo -> permutacao de rotulo de gene set (mesmo fallback usado
#     antes). NOTA: a tabela tem valores normalizados por amostra (Sui1-3,
#     Ctrl1-3), mas so para os genes JA FILTRADOS como significativos - nao
#     ajuda a estender o universo testavel, entao mantive o mesmo metodo
#     dos 3 artigos anteriores por consistencia.
#   - POP: reusei o GSEA real ja calculado do GSE53868 (microarray, mesma
#     plataforma tipo microarray que esta tabela de SUI - mesmo criterio
#     de pareamento ja usado para os 3 artigos anteriores, nao o RNA-seq
#     GSE208261), sem recalcular nada.
# =============================================================================

suppressMessages({
  library(org.Hs.eg.db)
  library(readxl)
})
set.seed(123)

## -----------------------------------------------------------------------
## 1) Carregar a tabela da usuaria (2 abas: up e down), linha 18 = cabecalho
## -----------------------------------------------------------------------
f <- "data/43032_2020_144_SUI_vs_Ctrl_microarray_MOESM2.xls"
up   <- read_excel(f, sheet = "up_Sui_vs_Ctrl",   skip = 17)
down <- read_excel(f, sheet = "down_Sui_vs_Ctrl", skip = 17)
sui_tab <- rbind(up[, c("EntrezID","GeneSymbol","P-value","FDR","Fold Change","Regulation")],
                  down[, c("EntrezID","GeneSymbol","P-value","FDR","Fold Change","Regulation")])
names(sui_tab) <- c("EntrezID","GeneSymbol","pvalue","FDR","FoldChange","Regulation")
sui_tab <- sui_tab[!is.na(sui_tab$EntrezID), ]
sui_tab$EntrezID <- as.character(as.integer(sui_tab$EntrezID))
cat("Tabela SUI carregada:", nrow(sui_tab), "probes (", sum(sui_tab$Regulation=="up"), "up,",
    sum(sui_tab$Regulation=="down"), "down)\n")

# Varios probes podem bater no mesmo EntrezID - mantem so o de p-valor mais
# significativo por gene (deduplicacao padrao em analises de microarray).
sui_tab <- sui_tab[order(sui_tab$pvalue), ]
sui_tab <- sui_tab[!duplicated(sui_tab$EntrezID), ]
cat("Apos deduplicar por EntrezID (melhor probe por gene):", nrow(sui_tab), "genes unicos\n\n")

write.csv(sui_tab, "results/SUI_43032_2020_144_combined_deduplicated.csv", row.names = FALSE)

## -----------------------------------------------------------------------
## 2) POP: GSEA real ja calculado (GSE208261, RNA-seq - pedido explicito da
##    usuaria, CORRIGINDO a escolha anterior de GSE53868/microarray), reusado
##    sem recalcular.
## -----------------------------------------------------------------------
pop_gsea <- read.csv("results/GSEA_POP_GSE208261_RNAseq_KEGG.csv")
pop_gsea$PATH <- sprintf("%05d", as.integer(pop_gsea$PATH))  # normaliza (zero a esquerda)
pop_gsea_sig <- subset(pop_gsea, p.adjust < 0.25)
cat("POP (GSE208261, RNA-seq) GSEA real, reusado: ", nrow(pop_gsea_sig), "de", nrow(pop_gsea),
    "vias KEGG significativas a FDR<0.25;",
    sum(pop_gsea$p.adjust < 0.05, na.rm = TRUE), "a FDR<0.05\n\n")

## -----------------------------------------------------------------------
## 3) KEGG gene sets via ENTREZID (direto, sem conversao de simbolo)
## -----------------------------------------------------------------------
all_eg <- keys(org.Hs.eg.db, keytype = "ENTREZID")
ann_path <- suppressWarnings(select(org.Hs.eg.db, keys = all_eg, keytype = "ENTREZID", columns = "PATH"))
ann_path <- ann_path[!is.na(ann_path$PATH), ]
kegg_sets_all <- split(ann_path$ENTREZID, ann_path$PATH)
cat("KEGG pathway definitions carregadas:", length(kegg_sets_all), "vias\n\n")

## -----------------------------------------------------------------------
## 4) Preranked GSEA pequeno-painel (identico ao script anterior)
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
  testable <- lapply(kegg_sets_all, function(g) which(ranked_genes %in% g))
  testable <- testable[sapply(testable, length) >= min_gs]
  cat("  KEGG vias testaveis (>=", min_gs, "membros presentes nesta lista de", N, "genes):",
      length(testable), "\n")
  if (length(testable) == 0) return(NULL)

  es_obs <- sapply(testable, calc_es, scores_abs = scores_abs, N = N)

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

score <- build_score(sui_tab$Regulation, sui_tab$pvalue)
ord <- order(-score)
ranked_genes <- sui_tab$EntrezID[ord]
ranked_scores <- score[ord]

cat("=== SUI (43032_2020_144) preranked GSEA ===\n")
sui_res <- run_preranked_small(ranked_genes, ranked_scores, kegg_sets_all)
if (!is.null(sui_res)) {
  cat("  Significativas a FDR<0.25:", sum(sui_res$p.adjust < 0.25, na.rm = TRUE),
      "| FDR<0.05:", sum(sui_res$p.adjust < 0.05, na.rm = TRUE), "\n")
  print(head(sui_res[, c("PATH","Nh","ES","NES","pvalue","p.adjust")], 15))
}
write.csv(sui_res, "results/GSEA_preranked_SUI_43032_2020_144_KEGG.csv", row.names = FALSE)

## -----------------------------------------------------------------------
## 5) Vias compartilhadas (mesmo criterio dos 3 artigos anteriores: FDR<0.25
##    dos dois lados)
## -----------------------------------------------------------------------
cat("\n=== Vias compartilhadas: SUI (43032_2020_144) vs POP (GSE208261, RNA-seq) ===\n")
sui_sig <- subset(sui_res, p.adjust < 0.25)
common_ids <- intersect(sui_sig$PATH, pop_gsea_sig$PATH)
cat("SUI significativas (FDR<0.25):", nrow(sui_sig),
    "| POP significativas (FDR<0.25):", nrow(pop_gsea_sig),
    "| COMPARTILHADAS:", length(common_ids), "\n\n")

if (length(common_ids) > 0) {
  # Usa sinal de ES (sempre disponivel), nao NES - o GSEA_POP_GSE53868_KEGG.csv
  # ja calculado anteriormente tem NES=NA em 30/218 vias (caso degenerado em
  # que todas as permutacoes caem do mesmo lado, falha o calculo da base de
  # normalizacao) - sinal de ES e equivalente e nao tem esse problema.
  shared <- merge(sui_sig[sui_sig$PATH %in% common_ids, c("PATH","Nh","ES","NES","pvalue","p.adjust","leadingEdge")],
                  pop_gsea_sig[pop_gsea_sig$PATH %in% common_ids, c("PATH","Nh","ES","NES","pvalue","p.adjust","leadingEdge")],
                  by = "PATH", suffixes = c("_SUI","_POP"))
  shared$Same_direction <- sign(shared$ES_SUI) == sign(shared$ES_POP)
  print(shared[, c("PATH","Nh_SUI","ES_SUI","Nh_POP","ES_POP","Same_direction","p.adjust_SUI","p.adjust_POP")])
  write.csv(shared, "results/GSEA_shared_KEGG_SUI43032_2020_144_vs_POP_GSE208261.csv", row.names = FALSE)
}

cat("\n=== FIM DO SCRIPT ===\n")
