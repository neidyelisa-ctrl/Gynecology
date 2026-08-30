# =============================================================================
# Duas fontes de SUI (Tong et al. 2010 e Chen et al. 2006), cada uma cruzada
# com o POP (GSE53868), em DOIS PASSOS por artigo, como pedido:
#   1) genes em comum / concordancia de direcao
#   2) se o passo 1 nao for muito promissor, comparacao a nivel de VIAS
#      (GO/KEGG rodado separadamente em cada lista, depois procurando termos
#      em comum entre as duas listas de vias - convergencia biologica mesmo
#      sem sobreposicao de genes individuais)
#
# FONTE 1 (NOVA): Tong, Lang & Zhu (2010) Int Urogynecol J 21:1545-1551 -
#   "Microarray analysis of differentially expressed genes in vaginal
#   tissues in postmenopausal women. The role of stress urinary
#   incontinence". 3 pares de mulheres POS-MENOPAUSA, SUI vs continentes
#   (idade/IMC/paridade pareados), array Affymetrix U133 Plus 2.0. Tabela 2:
#   75 genes (Ratio = razao SUI/controle, nao log2 apesar do nome da coluna
#   no artigo - valor>1 = up, valor<1 = down; confirmado com o texto do
#   artigo, que descreve os mesmos 75 genes como 31 up / 44 down). 66 genes
#   unicos curados em data/tong2010_75genes.csv (24 up / 42 down; diferenca
#   de contagem em relacao ao artigo por causa de linhas duplicadas na
#   Tabela 2 do artigo, ex. NFKBIZ, APOE e GAS2L1 aparecem 2x - mantive so
#   1 vez por gene).
#   IMPORTANTE: e o UNICO dos artigos de SUI usados neste projeto com
#   mulheres POS-menopausa (Chen 2003/2006 e Poelmans sao pre-menopausa ou
#   nao especificam) - relevante para interpretacao caso a concordancia
#   seja fraca.
#
# FONTE 2 (RE-TESTE): Chen et al. 2006, Hum Reprod 21:22-29 (o PDF enviado
#   agora pela usuaria como "chen2005_DEG.pdf" - CONFERIDO: e o MESMO
#   artigo publicado em 2006, fase secretora, 79 genes, ja usado para
#   construir data/chen2006_79genes.csv nesta sessao anteriormente. Nao e
#   uma fonte nova - e o mesmo Chen 2006 sendo re-analisado com o mesmo
#   framework de 2 passos usado para o Tong 2010, para comparacao lado a
#   lado). Reaproveita os resultados ja calculados (crossing
#   results/Chen2006_x_GSE53868.csv, GO/KEGG standalone
#   results/GO_BP_Chen2006_58genes.csv / KEGG_Chen2006_58genes.csv).
# =============================================================================

suppressMessages({
  library(org.Hs.eg.db)
  library(GO.db)
})

pop <- read.csv("results/GSE53868_limma_completo.csv")
universe_pop <- pop$Gene
cat("POP (GSE53868):", nrow(pop), "genes testados\n\n")

run_enrichment <- function(hit_genes, universe_genes, keytype_col, ont_filter = NULL, min_gs = 2, max_gs = 2000) {
  hit_genes <- unique(intersect(hit_genes, universe_genes)); universe_genes <- unique(universe_genes)
  if (length(hit_genes) == 0) return(NULL)
  ann <- suppressWarnings(select(org.Hs.eg.db, keys = universe_genes, keytype = "SYMBOL", columns = keytype_col))
  ann <- ann[!is.na(ann[[keytype_col]]), c("SYMBOL", keytype_col)]
  if (!is.null(ont_filter)) {
    ann2 <- suppressWarnings(select(org.Hs.eg.db, keys = universe_genes, keytype = "SYMBOL", columns = c("GO","ONTOLOGY")))
    keep <- unique(ann2$GO[ann2$ONTOLOGY == ont_filter & !is.na(ann2$GO)])
    ann <- ann[ann[[keytype_col]] %in% keep, ]
  }
  colnames(ann) <- c("SYMBOL", "TERM_ID"); ann <- unique(ann)
  gs <- table(ann$TERM_ID); valid <- names(gs)[gs >= min_gs & gs <= max_gs]; ann <- ann[ann$TERM_ID %in% valid, ]
  N <- length(universe_genes); n <- length(hit_genes)
  res <- lapply(valid, function(t) {
    g <- unique(ann$SYMBOL[ann$TERM_ID == t]); M <- length(g)
    h <- intersect(g, hit_genes); k <- length(h)
    if (k == 0) return(NULL)
    p <- phyper(k - 1, M, N - M, n, lower.tail = FALSE)
    data.frame(TERM_ID = t, Count = k, M = M, N = N, n = n, pvalue = p, geneID = paste(h, collapse = "/"))
  })
  res <- do.call(rbind, res); if (is.null(res)) return(NULL)
  res$p.adjust <- p.adjust(res$pvalue, "BH"); res[order(res$pvalue), ]
}

run_go_kegg <- function(genes, universe, tag) {
  cat("--- GO/KEGG:", tag, "(n=", length(genes), "genes) ---\n")
  go <- run_enrichment(genes, universe, "GO", ont_filter = "BP")
  if (!is.null(go)) {
    terms <- suppressMessages(select(GO.db, keys = go$TERM_ID, keytype = "GOID", columns = "TERM"))
    go <- merge(go, terms, by.x = "TERM_ID", by.y = "GOID"); go <- go[order(go$pvalue), ]
    cat("GO BP sig (FDR<0.05):", sum(go$p.adjust < 0.05, na.rm = TRUE), "de", nrow(go), "\n")
  } else cat("GO: nenhum termo\n")
  kegg <- run_enrichment(genes, universe, "PATH")
  if (!is.null(kegg)) cat("KEGG sig (FDR<0.05):", sum(kegg$p.adjust < 0.05, na.rm = TRUE), "de", nrow(kegg), "\n")
  else cat("KEGG: nenhuma via\n")
  cat("\n")
  list(go = go, kegg = kegg)
}

# compara duas tabelas de enriquecimento (cada uma ja rodada em sua propria
# lista) e retorna os termos SIGNIFICATIVOS (FDR<0.05) em comum nas DUAS
compare_pathways <- function(tabA, tabB, id_col = "TERM_ID") {
  if (is.null(tabA) || is.null(tabB)) return(NULL)
  sigA <- subset(tabA, p.adjust < 0.05)
  sigB <- subset(tabB, p.adjust < 0.05)
  common_ids <- intersect(sigA[[id_col]], sigB[[id_col]])
  if (length(common_ids) == 0) return(data.frame())
  out <- merge(sigA[sigA[[id_col]] %in% common_ids, ], sigB[sigB[[id_col]] %in% common_ids, ],
               by = id_col, suffixes = c("_A", "_B"))
  out[order(out$pvalue_A), ]
}

## =============================================================================
## PARTE 1 - Tong et al. 2010 (66 genes, pos-menopausa) x GSE53868
## =============================================================================
cat("\n============================================================\n")
cat("PARTE 1: Tong et al. 2010 (SUI, pos-menopausa) x GSE53868 (POP)\n")
cat("============================================================\n\n")

tong <- read.csv("data/tong2010_75genes.csv")
tong_genes <- unique(tong[, c("GeneSymbol", "Direction", "Ratio")])
cat("Genes curados do Tong 2010:", nrow(tong_genes),
    "(", sum(tong_genes$Direction == "up"), "up /",
    sum(tong_genes$Direction == "down"), "down )\n\n")

cross_tong <- merge(tong_genes, pop, by.x = "GeneSymbol", by.y = "Gene", all.x = TRUE)
cross_tong$Testado <- !is.na(cross_tong$logFC)
cross_tong$Sig_POP_FDR05 <- !is.na(cross_tong$adj.P.Val) & cross_tong$adj.P.Val < 0.05 & abs(cross_tong$logFC) > 0.5
cross_tong$Direcao_POP <- ifelse(is.na(cross_tong$logFC), NA, ifelse(cross_tong$logFC > 0, "up", "down"))
cross_tong$Concordante <- cross_tong$Direction == cross_tong$Direcao_POP
cross_tong <- cross_tong[order(cross_tong$adj.P.Val), ]
write.csv(cross_tong, "results/Tong2010_x_GSE53868.csv", row.names = FALSE)

cat("Genes testados no GSE53868:", sum(cross_tong$Testado), "de", nrow(cross_tong), "\n")
sig_tong <- subset(cross_tong, Sig_POP_FDR05)
cat("Genes SIGNIFICATIVOS no GSE53868 entre os candidatos do Tong 2010:", nrow(sig_tong), "\n")
if (nrow(sig_tong) > 0) print(sig_tong[, c("GeneSymbol","Direction","Ratio","logFC","adj.P.Val","Direcao_POP","Concordante")])

testaveis_tong <- subset(cross_tong, Testado & !is.na(Direcao_POP))
n_conc_tong <- sum(testaveis_tong$Concordante, na.rm = TRUE)
n_tot_tong <- nrow(testaveis_tong)
cat("\nConcordancia de direcao:", n_conc_tong, "de", n_tot_tong, "(", round(100*n_conc_tong/n_tot_tong,1), "%)\n")
bt_tong <- binom.test(n_conc_tong, n_tot_tong, p = 0.5)
cat("Teste binomial de sinal: p =", format(bt_tong$p.value, digits = 4), "\n\n")

promising_tong <- (nrow(sig_tong) >= 2) || (bt_tong$p.value < 0.05)
cat("Passo 1 (genes) foi promissor?", promising_tong, "\n\n")

concordant_tong <- unique(testaveis_tong$GeneSymbol[testaveis_tong$Concordante])
res_go_tong_cross <- NULL
if (promising_tong && length(concordant_tong) >= 3) {
  res_go_tong_cross <- run_go_kegg(concordant_tong, universe_pop, "genes concordantes Tong2010 x GSE53868")
}

# --- Passo 2: comparacao a nivel de VIAS (rodado sempre, para dar o panorama
#     completo mesmo quando o passo 1 e fraco) ---
cat("--- Passo 2: GO/KEGG do painel Tong2010 sozinho (66 genes) ---\n")
res_go_tong_standalone <- run_go_kegg(tong_genes$GeneSymbol, keys(org.Hs.eg.db, keytype = "SYMBOL"), "painel Tong2010 sozinho")

go_pop_standalone <- read.csv("results/GO_BP_GSE53868.csv")
kegg_pop_standalone <- read.csv("results/KEGG_GSE53868.csv")

cat("--- Vias em comum: Tong2010 (painel) x GSE53868 (DEGs) ---\n")
go_common_tong <- compare_pathways(res_go_tong_standalone$go, go_pop_standalone)
kegg_common_tong <- compare_pathways(res_go_tong_standalone$kegg, kegg_pop_standalone)
cat("GO BP em comum:", if (is.null(go_common_tong)) 0 else nrow(go_common_tong), "\n")
cat("KEGG em comum:", if (is.null(kegg_common_tong)) 0 else nrow(kegg_common_tong), "\n\n")
if (!is.null(go_common_tong) && nrow(go_common_tong) > 0) {
  write.csv(go_common_tong, "results/GO_BP_vias_comuns_Tong2010_GSE53868.csv", row.names = FALSE)
  print(head(go_common_tong[, c("TERM_A","Count_A","geneID_A","Count_B","geneID_B")], 10))
}
if (!is.null(kegg_common_tong) && nrow(kegg_common_tong) > 0) {
  write.csv(kegg_common_tong, "results/KEGG_vias_comuns_Tong2010_GSE53868.csv", row.names = FALSE)
  print(head(kegg_common_tong[, c("TERM_ID","Count_A","geneID_A","Count_B","geneID_B")], 10))
}

## =============================================================================
## PARTE 2 - Chen et al. 2006 (58 genes) x GSE53868 - RE-TESTE, mesmo
##   framework de 2 passos, reaproveitando o que ja foi calculado antes
## =============================================================================
cat("\n============================================================\n")
cat("PARTE 2: Chen et al. 2006 (SUI, fase secretora - re-teste) x GSE53868\n")
cat("============================================================\n\n")

cross_chen06 <- read.csv("results/Chen2006_x_GSE53868.csv")
sig_chen06 <- subset(cross_chen06, Significativo_GSE53868_FDR05)
testaveis_chen06 <- subset(cross_chen06, Testado_no_GSE53868 & !is.na(Direcao_GSE53868))
n_conc_chen06 <- sum(testaveis_chen06$Concordante, na.rm = TRUE)
n_tot_chen06 <- nrow(testaveis_chen06)
bt_chen06 <- binom.test(n_conc_chen06, n_tot_chen06, p = 0.5)

cat("Genes individualmente significativos:", nrow(sig_chen06), "de", nrow(cross_chen06), "\n")
cat("Concordancia de direcao:", n_conc_chen06, "de", n_tot_chen06,
    "(", round(100*n_conc_chen06/n_tot_chen06,1), "%), binom p =",
    format(bt_chen06$p.value, digits = 4), "\n")
promising_chen06 <- (nrow(sig_chen06) >= 2) || (bt_chen06$p.value < 0.05)
cat("Passo 1 (genes) foi promissor?", promising_chen06,
    "(ja confirmado em analise anterior desta sessao)\n\n")

# Passo 2, mesmo assim, para comparacao formal lado a lado com o Tong 2010
go_chen06_standalone <- read.csv("results/GO_BP_Chen2006_58genes.csv")
kegg_chen06_standalone <- read.csv("results/KEGG_Chen2006_58genes.csv")

cat("--- Vias em comum: Chen2006 (painel) x GSE53868 (DEGs) ---\n")
go_common_chen06 <- compare_pathways(go_chen06_standalone, go_pop_standalone)
kegg_common_chen06 <- compare_pathways(kegg_chen06_standalone, kegg_pop_standalone)
cat("GO BP em comum:", if (is.null(go_common_chen06)) 0 else nrow(go_common_chen06), "\n")
cat("KEGG em comum:", if (is.null(kegg_common_chen06)) 0 else nrow(kegg_common_chen06), "\n\n")
if (!is.null(go_common_chen06) && nrow(go_common_chen06) > 0) {
  write.csv(go_common_chen06, "results/GO_BP_vias_comuns_Chen2006_GSE53868.csv", row.names = FALSE)
  print(head(go_common_chen06[, c("TERM_A","Count_A","geneID_A","Count_B","geneID_B")], 10))
}
if (!is.null(kegg_common_chen06) && nrow(kegg_common_chen06) > 0) {
  write.csv(kegg_common_chen06, "results/KEGG_vias_comuns_Chen2006_GSE53868.csv", row.names = FALSE)
  print(head(kegg_common_chen06[, c("TERM_ID","Count_A","geneID_A","Count_B","geneID_B")], 10))
}

cat("\n=== FIM DO SCRIPT ===\n")
