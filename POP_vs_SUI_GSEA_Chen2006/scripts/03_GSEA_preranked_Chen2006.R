# =============================================================================
# Script 3/4 — GSEA preranked nos DEGs do Chen et al. 2006 (SUI), vias KEGG.
#
# Artigo: Chen B, Wen Y, Zhang Z, Guo Y, Warrington JA, Polan ML. "Microarray
# analysis of differentially expressed genes in vaginal tissues from women
# with stress urinary incontinence compared with asymptomatic women." Hum
# Reprod. 2006;21(1):22-29 (publicado online em ago/2005 — é o mesmo PDF
# enviado pela usuária como "chen2005_DEG.pdf"). 5 pares SUI x continentes,
# parede vaginal periuretral, fase secretória do ciclo menstrual, array
# Affymetrix U133A. O artigo já aplica seus próprios critérios de
# significância (MAS 5.0 + RMA, testes paramétrico/não-paramétrico + PAM,
# correção Bonferroni/Holm/Hochberg/BH/Westfall-Young, p<0.05) e reporta a
# lista final de 79 DEGs comuns às Tabelas II (39 up) e III (40 down) — é
# essa lista final (não o transcriptoma bruto do array, que o artigo não
# disponibiliza) que serve de entrada aqui.
#
# CURADORIA (data/chen2006_79genes.csv): das 79 linhas das Tabelas II/III,
# 60 foram mapeadas com confiança para o símbolo HGNC atual (34 up / 26
# down); 19 permanecem genuinamente ambíguas com a nomenclatura de 2005
# ("hypothetical protein FLJxxxxx", "KIAAxxxx protein", "Zinc finger
# protein" genérico, sondas de controle) e foram excluídas — não foram
# inventados símbolos para elas.
#
# LIMITAÇÃO METODOLÓGICA IMPORTANTE (leia antes de interpretar os números):
# GSEA foi desenhado para ranquear um TRANSCRIPTOMA INTEIRO (milhares de
# genes), não uma lista curta já pré-filtrada como significativa pelos
# próprios autores. Rodar "GSEA preranked" nestes 60 genes é um uso
# NÃO-PADRÃO do método:
#   (a) uma via KEGG só é "testável" se >=2 dos seus membros estiverem
#       nesta lista de 60 genes — a maioria das vias KEGG terá 0 ou 1
#       membro por acaso numa lista tão curta, então poucas vias são
#       testáveis;
#   (b) sem a matriz de expressão bruta por amostra do artigo original
#       (não publicada), não é possível fazer permutação de FENÓTIPO como
#       no script 02 (POP). A permutação usada aqui é de RÓTULO DE GENE SET
#       (sortear posições aleatórias do mesmo tamanho do gene set, mantendo
#       o ranking fixo) — o mesmo fallback que a própria ferramenta oficial
#       GSEA-Preranked usa quando não há matriz bruta, mas é o tipo de
#       permutação mais liberal/menos conservador (genes correlacionados
#       numa via real podem inflar a significância — documentado no próprio
#       Subramanian et al. 2005).
# Resultado: leia os números abaixo como um sinal exploratório, de baixo
# poder estatístico — não no mesmo patamar do GSEA "clássico" do POP
# (script 02), que usa o transcriptoma inteiro com permutação de fenótipo.
#
# ESCORE DE RANKING (não há estatística-t contínua publicada, só p-valor +
# fold-change por gene): score = sinal(direção) * -log10(p-valor RMA).
# P-valores impressos como 0 na Tabela II/III são "floored" para metade do
# menor p-valor não-nulo do painel, para evitar escore infinito.
# =============================================================================

suppressMessages(library(org.Hs.eg.db))
set.seed(123)

## --- 1) Ler e preparar a lista de DEGs do Chen 2006 -------------------------
chen <- read.csv("data/chen2006_79genes.csv")
chen <- chen[!duplicated(chen$GeneSymbol), ]
cat("Chen 2006 — DEGs curados com símbolo HGNC atual:", nrow(chen),
    "(", sum(chen$Direction == "up"), "up /", sum(chen$Direction == "down"), "down )\n\n")

build_score <- function(direction, pvalue) {
  p <- pvalue
  floor_p <- min(p[p > 0], na.rm = TRUE) / 2
  p[p <= 0] <- floor_p
  sign <- ifelse(direction == "up", 1, -1)
  sign * -log10(p)
}

score <- build_score(chen$Direction, chen$RMA_pvalue)
ord <- order(-score)
ranked_genes <- chen$GeneSymbol[ord]
ranked_scores <- score[ord]
N <- length(ranked_genes)

## --- 2) Vias KEGG (offline, org.Hs.eg.db) -----------------------------------
all_sym <- keys(org.Hs.eg.db, keytype = "SYMBOL")
ann_path <- suppressWarnings(select(org.Hs.eg.db, keys = all_sym, keytype = "SYMBOL", columns = "PATH"))
ann_path <- ann_path[!is.na(ann_path$PATH), ]
kegg_sets_all <- split(ann_path$SYMBOL, ann_path$PATH)
cat("Vias KEGG carregadas (universo completo, org.Hs.eg.db):", length(kegg_sets_all), "\n\n")

## --- 3) GSEA preranked em painel pequeno: ES ponderado + permutação -------
##        de rótulo de gene set (ver nota metodológica acima)
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

set.seed(123)
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
write.csv(res, "results/GSEA_preranked_Chen2006_KEGG.csv", row.names = FALSE)

cat("\n=== GSEA preranked do Chen 2006 (SUI, 60 genes, KEGG) ===\n")
cat("Vias significativas a FDR<0.25:", sum(res$p.adjust < 0.25, na.rm = TRUE), "\n")
cat("Vias significativas a FDR<0.05:", sum(res$p.adjust < 0.05, na.rm = TRUE), "\n\n")
print(head(res[, c("PATH","Nh","ES","NES","pvalue","p.adjust","leadingEdge")], 10))
