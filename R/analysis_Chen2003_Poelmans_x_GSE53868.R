# =============================================================================
# Duas fontes NOVAS de genes candidatos de SUI (literatura), cruzadas
# SEPARADAMENTE (nao combinadas) contra o POP (GSE53868, ja processado com
# limma em results/GSE53868_limma_completo.csv - script original:
# R/analysis_GSE53868_limma_x_Chen2006.R).
#
# FONTE 1: Chen et al. 2003, Am J Obstet Gynecol 189:89-97 - "Menstrual
#   phase-dependent gene expression differences in periurethral vaginal
#   tissue from women with stress incontinence". Fase proliferativa, 5
#   pares SUI vs continentes, array HuGeneFL (6800 genes). Tabela II: 90
#   genes candidatos (62 up / 28 down no artigo original); 69 curados com
#   simbolo HGNC confiavel em data/chen2003_90genes.csv (43 up / 26 down -
#   21 excluidos por serem entradas ambiguas da nomenclatura de 2003, ex.
#   clones/BACs sem nome, "T-cell antigen receptor" generico, transposase
#   de camundongo que nao e gene humano real).
#   -> Dado de EXPRESSAO (tem direcao), mesmo tratamento estatistico usado
#      antes para o Chen 2006: genes ja considerados significativos pelos
#      autores originais, cruzados por (a) significancia individual no
#      GSE53868 e (b) teste binomial de concordancia de direcao no painel
#      inteiro.
#
# FONTE 2: Poelmans et al. 2023, material suplementar (Tabela S1) - 188
#   genes candidatos de SUI selecionados por GWAS (gene-wide p<1e-3 em
#   pelo menos um de 4 estudos: Penney et al. 2020, Cartwright et al.,
#   HUNT, UK Biobank); 183 extraidos com sucesso em
#   data/poelmans_2023_SUI_GWAS_188genes.csv.
#   -> Dado GENETICO/ASSOCIACAO, SEM direcao de expressao (p-valor de
#      associacao GWAS nao tem "sentido para cima/baixo" equivalente a
#      fold-change). Tratamento correto e DIFERENTE do Chen: e uma busca
#      de genes-candidato (esses genes geneticamente associados a SUI
#      tambem aparecem como DEG no tecido de POP?), nao um teste de
#      concordancia de sinal.
# =============================================================================

suppressMessages({
  library(org.Hs.eg.db)
  library(GO.db)
})

pop <- read.csv("results/GSE53868_limma_completo.csv")
pop_sig <- subset(pop, adj.P.Val < 0.05 & abs(logFC) > 0.5)
cat("POP (GSE53868): ", nrow(pop_sig), " DEGs (padj<0.05, |logFC|>0.5) de ", nrow(pop), " genes testados\n\n", sep = "")

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
    cat("GO BP: ", sum(go$p.adjust < 0.05, na.rm = TRUE), " termos sig (FDR<0.05) de ", nrow(go), " testados\n", sep = "")
    print(head(go[, c("TERM","Count","pvalue","p.adjust","geneID")], 8))
  } else cat("GO: nenhum termo\n")
  kegg <- run_enrichment(genes, universe, "PATH")
  if (!is.null(kegg)) {
    cat("\nKEGG: ", sum(kegg$p.adjust < 0.05, na.rm = TRUE), " vias sig (FDR<0.05) de ", nrow(kegg), " testadas\n", sep = "")
    print(head(kegg[, c("TERM_ID","Count","pvalue","p.adjust","geneID")], 8))
  } else cat("KEGG: nenhuma via\n")
  cat("\n")
  list(go = go, kegg = kegg)
}

universe_pop <- pop$Gene

## =============================================================================
## PARTE 1 - Chen et al. 2003 (69 genes curados) x GSE53868
## =============================================================================
cat("\n============================================================\n")
cat("PARTE 1: Chen et al. 2003 (SUI, fase proliferativa) x GSE53868 (POP)\n")
cat("============================================================\n\n")

chen03 <- read.csv("data/chen2003_90genes.csv")
chen03_genes <- unique(chen03[, c("GeneSymbol", "Direction", "AvgFoldChange")])
chen03_genes <- chen03_genes[!duplicated(chen03_genes$GeneSymbol), ]
cat("Genes curados do Chen 2003:", nrow(chen03_genes),
    "(", sum(chen03_genes$Direction == "up"), "up /",
    sum(chen03_genes$Direction == "down"), "down )\n\n")

cross03 <- merge(chen03_genes, pop, by.x = "GeneSymbol", by.y = "Gene", all.x = TRUE)
cross03$Testado_no_GSE53868 <- !is.na(cross03$logFC)
cross03$Significativo_GSE53868_FDR05 <- !is.na(cross03$adj.P.Val) & cross03$adj.P.Val < 0.05 & abs(cross03$logFC) > 0.5
cross03$Direcao_GSE53868 <- ifelse(is.na(cross03$logFC), NA, ifelse(cross03$logFC > 0, "up", "down"))
cross03$Concordante <- cross03$Direction == cross03$Direcao_GSE53868
cross03 <- cross03[order(cross03$adj.P.Val), ]
write.csv(cross03, "results/Chen2003_x_GSE53868.csv", row.names = FALSE)

cat("Genes testados no GSE53868:", sum(cross03$Testado_no_GSE53868), "de", nrow(cross03), "\n")
sig03 <- subset(cross03, Significativo_GSE53868_FDR05)
cat("Genes SIGNIFICATIVOS no GSE53868 entre os candidatos do Chen 2003:", nrow(sig03), "\n")
if (nrow(sig03) > 0) print(sig03[, c("GeneSymbol","Direction","AvgFoldChange","logFC","adj.P.Val","Direcao_GSE53868","Concordante")])

testaveis03 <- subset(cross03, Testado_no_GSE53868 & !is.na(Direcao_GSE53868))
n_conc03 <- sum(testaveis03$Concordante, na.rm = TRUE)
n_tot03 <- nrow(testaveis03)
cat("\nConcordancia de direcao (Chen 2003 SUI vs GSE53868 POP):", n_conc03, "de", n_tot03, "\n")
bt03 <- binom.test(n_conc03, n_tot03, p = 0.5)
cat("Teste binomial de sinal: p =", format(bt03$p.value, digits = 4), "\n\n")

concordant03 <- unique(testaveis03$GeneSymbol[testaveis03$Concordante])
res_go_chen03 <- run_go_kegg(concordant03, universe_pop, "genes concordantes Chen 2003 x GSE53868")
if (!is.null(res_go_chen03$go)) write.csv(res_go_chen03$go, "results/GO_BP_Chen2003_x_GSE53868_concordantes.csv", row.names = FALSE)
if (!is.null(res_go_chen03$kegg)) write.csv(res_go_chen03$kegg, "results/KEGG_Chen2003_x_GSE53868_concordantes.csv", row.names = FALSE)

## =============================================================================
## PARTE 2 - Poelmans et al. 2023 (183 genes GWAS) x GSE53868
##   SEM direcao (dado genetico/associacao, nao expressao) -> tratado como
##   busca de genes-candidato, nao teste de concordancia de sinal.
## =============================================================================
cat("\n============================================================\n")
cat("PARTE 2: Poelmans et al. 2023 (SUI, GWAS) x GSE53868 (POP)\n")
cat("============================================================\n\n")

poel <- read.csv("data/poelmans_2023_SUI_GWAS_188genes.csv")
cat("Genes candidatos GWAS (Poelmans 2023):", nrow(poel), "\n")

parse_p <- function(x) as.numeric(gsub("[*†]", "", x))
sig_flag <- function(x) grepl("[*†]", x)
poel$n_cohorts_sig <- rowSums(sapply(poel[, c("Penney_p","Cartwright_p","HUNT_p","UKBiobank_p")], sig_flag))
poel$min_p <- apply(sapply(poel[, c("Penney_p","Cartwright_p","HUNT_p","UKBiobank_p")], parse_p), 1, min, na.rm = TRUE)
poel_highconf <- subset(poel, n_cohorts_sig >= 2)
cat("Subconjunto de alta confianca (sig. em >=2 dos 4 estudos GWAS):", nrow(poel_highconf), "genes\n\n")

cross_poel <- merge(poel, pop, by.x = "Gene", by.y = "Gene", all.x = TRUE)
cross_poel$Testado_no_GSE53868 <- !is.na(cross_poel$logFC)
cross_poel$DEG_no_GSE53868_FDR05 <- !is.na(cross_poel$adj.P.Val) & cross_poel$adj.P.Val < 0.05 & abs(cross_poel$logFC) > 0.5
cross_poel <- cross_poel[order(cross_poel$adj.P.Val), ]
write.csv(cross_poel, "results/Poelmans_x_GSE53868.csv", row.names = FALSE)

cat("Genes GWAS testados no array do GSE53868:", sum(cross_poel$Testado_no_GSE53868), "de", nrow(cross_poel), "\n")
overlap_poel <- subset(cross_poel, DEG_no_GSE53868_FDR05)
cat("Genes GWAS de SUI que SAO TAMBEM DEG no tecido de POP (GSE53868, padj<0.05,|logFC|>0.5):", nrow(overlap_poel), "\n")
if (nrow(overlap_poel) > 0) {
  print(overlap_poel[, c("Gene","n_cohorts_sig","min_p","logFC","adj.P.Val")])
}
write.csv(overlap_poel, "results/Poelmans_x_GSE53868_overlap_DEG.csv", row.names = FALSE)

# alta confianca (>=2 estudos GWAS) que tambem sao DEG de POP
overlap_highconf <- subset(overlap_poel, n_cohorts_sig >= 2)
cat("\n...dos quais, com corroboracao em >=2 estudos GWAS:", nrow(overlap_highconf), "\n")
if (nrow(overlap_highconf) > 0) print(overlap_highconf[, c("Gene","n_cohorts_sig","min_p","logFC","adj.P.Val")])

# genes da Tabela S2 do Poelmans (evidencia de SUI FORA de GWAS - landscape genes +
# interatores: AAT/SERPINA1, BDNF, CDH1, CTNNB1, ESR1, ESR2, GNAI3, ITGA8, ITGB1, MMP1,
# PARP1 - lista extraida diretamente do texto do docx, nao inventada) - checamos quais
# estao presentes e o status deles no array GSE53868.
landscape_genes <- c("AAT","SERPINA1","BDNF","CDH1","CTNNB1","ESR1","ESR2","GNAI3","ITGA8","ITGB1","MMP1","PARP1")
landscape_in_array <- intersect(landscape_genes, pop$Gene)
cat("\nGenes da Tabela S2 (evidencia de SUI nao-GWAS) presentes no array GSE53868:", paste(landscape_in_array, collapse=", "), "\n")
if (length(landscape_in_array) > 0) {
  print(subset(pop, Gene %in% landscape_in_array)[, c("Gene","logFC","adj.P.Val")])
}

cat("\n")
if (nrow(overlap_poel) >= 5) {
  res_go_poel <- run_go_kegg(overlap_poel$Gene, universe_pop, "genes GWAS(Poelmans) que sao DEG no POP")
} else {
  cat("Overlap pequeno demais (n=", nrow(overlap_poel), ") para GO/KEGG com significado - ",
      "rodando em vez disso no PAINEL COMPLETO de 183 genes candidatos GWAS (caracterizacao ",
      "funcional do painel, independente de ser ou nao DEG no POP).\n", sep = "")
  res_go_poel <- run_go_kegg(poel$Gene, keys(org.Hs.eg.db, keytype = "SYMBOL"), "painel completo Poelmans (183 genes)")
}
if (!is.null(res_go_poel$go)) write.csv(res_go_poel$go, "results/GO_BP_Poelmans.csv", row.names = FALSE)
if (!is.null(res_go_poel$kegg)) write.csv(res_go_poel$kegg, "results/KEGG_Poelmans.csv", row.names = FALSE)

## =============================================================================
## PARTE 3 - Export para STRING (submissao manual pelo usuario) - PPI/hub genes
##   Live STRING API bloqueada neste ambiente (confirmado repetidamente); os
##   arquivos abaixo sao para a usuaria colar em string-db.org > Multiple
##   Proteins > Homo sapiens.
## =============================================================================
writeLines(concordant03, "results/STRING_input_Chen2003_concordantes_GSE53868.txt")
writeLines(unique(overlap_poel$Gene), "results/STRING_input_Poelmans_overlap_GSE53868.txt")
cat("\nArquivos para submissao manual no STRING (string-db.org) exportados:\n")
cat(" - results/STRING_input_Chen2003_concordantes_GSE53868.txt (", length(concordant03), " genes)\n", sep="")
cat(" - results/STRING_input_Poelmans_overlap_GSE53868.txt (", length(unique(overlap_poel$Gene)), " genes)\n", sep="")

cat("\n=== FIM DO SCRIPT ===\n")
