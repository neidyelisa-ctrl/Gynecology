# =============================================================================
# Script 6 — A abordagem ORIGINAL (antes do conselho do professor de rodar
# GSEA separado): cruzar DEG(POP) com a lista de DEG do Chen 2006 (SUI),
# olhar genes em comum / concordância de direção, e rodar GO/KEGG nesse
# conjunto cruzado.
#
# Por que isso é DIFERENTE do GSEA dos scripts 02-04: GSEA usa o
# transcriptoma INTEIRO ranqueado; este script usa só os genes que já
# passaram (ou quase) no corte de significância dos dois lados — é ORA
# (over-representation analysis / teste hipergeométrico), o método clássico
# de "genes em comum -> GO/KEGG". Mais simples, mais fácil de explicar numa
# banca, mas mais sensível ao tamanho pequeno da lista do Chen 2006 (59
# genes candidatos, não o transcriptoma inteiro do estudo original deles).
#
# Roda 3 versões, da mais para a menos rigorosa:
#   (a) interseção ESTRITA: gene tem que ser DEG (FDR<0.05, |log2FC|>1) no
#       POP E estar na lista do Chen 2006 (SUI) — critério mais rígido,
#       tende a sobrar poucos genes.
#   (b) genes CONCORDANTES: qualquer gene do Chen 2006 cuja direção
#       (para cima/para baixo) bate com a direção observada no POP,
#       independente do gene ser "significativo" no POP sozinho — critério
#       mais sensível, é o que rendeu o achado de queratinização.
#   (c) teste de sinal (binomial): a fração de genes concordantes é maior
#       que os 50% esperados por acaso?
# =============================================================================

suppressMessages({
  library(org.Hs.eg.db)
  library(GO.db)
})

## --- 1) Carregar POP (tabela completa, já calculada no script 01) ---------
pop <- read.csv("results/GSE53868_limma_completo.csv")
pop_deg <- subset(pop, adj.P.Val < 0.05 & abs(logFC) > 1)  # mesmo critério pedido (logFC>1, FDR<0.05)
cat("DEGs do POP (|log2FC|>1, FDR<0.05):", nrow(pop_deg), "de", nrow(pop), "testados\n\n")

## --- 2) Carregar Chen 2006 (SUI) -------------------------------------------
chen <- read.csv("data/chen2006_79genes.csv")
chen <- chen[!duplicated(chen$GeneSymbol), ]
cat("Chen 2006 (SUI) — genes candidatos:", nrow(chen), "\n\n")

## --- 3) Cruzamento gene a gene ---------------------------------------------
cross <- merge(chen[, c("GeneSymbol", "Direction")], pop, by.x = "GeneSymbol", by.y = "Gene", all.x = TRUE)
cross$Testado_no_POP <- !is.na(cross$logFC)
cross$Direcao_POP <- ifelse(is.na(cross$logFC), NA, ifelse(cross$logFC > 0, "up", "down"))
cross$Concordante <- cross$Direction == cross$Direcao_POP
cross$DEG_POP_FDR05_logFC1 <- !is.na(cross$adj.P.Val) & cross$adj.P.Val < 0.05 & abs(cross$logFC) > 1
cross <- cross[order(cross$adj.P.Val), ]
write.csv(cross, "results/06_cruzamento_Chen2006_x_POP_genes.csv", row.names = FALSE)

## (a) interseção estrita: DEG no POP (logFC>1, FDR<0.05) E na lista do Chen 2006
intersecao_estrita <- subset(cross, DEG_POP_FDR05_logFC1)
cat("=== (a) Interseção ESTRITA: DEG(POP, logFC>1 & FDR<0.05) ∩ Chen2006 ===\n")
cat(nrow(intersecao_estrita), "genes:", paste(intersecao_estrita$GeneSymbol, collapse = ", "), "\n\n")

## (b) genes concordantes em direção (independente de significância no POP)
testaveis <- subset(cross, Testado_no_POP & !is.na(Direcao_POP))
concordantes <- subset(testaveis, Concordante)
cat("=== (b) Genes CONCORDANTES em direção (Chen2006 x POP) ===\n")
cat(nrow(concordantes), "de", nrow(testaveis), "testáveis:",
    paste(concordantes$GeneSymbol, collapse = ", "), "\n\n")

## (c) teste de sinal
bt <- binom.test(nrow(concordantes), nrow(testaveis), p = 0.5)
cat("=== (c) Teste binomial de sinal ===\n")
cat("Concordância:", nrow(concordantes), "/", nrow(testaveis),
    "=", round(100 * nrow(concordantes) / nrow(testaveis), 1), "%\n")
cat("H0: concordância = 50% (acaso). p-valor =", format(bt$p.value, digits = 4), "\n\n")

## --- 4) GO/KEGG offline (teste hipergeométrico, universo = genes do array) -
run_enrichment <- function(hit_genes, universe_genes, keytype_col, ont_filter = NULL, min_gs = 2, max_gs = 2000) {
  hit_genes <- unique(intersect(hit_genes, universe_genes)); universe_genes <- unique(universe_genes)
  if (length(hit_genes) < 2) return(NULL)
  ann <- suppressWarnings(select(org.Hs.eg.db, keys = universe_genes, keytype = "SYMBOL", columns = keytype_col))
  ann <- ann[!is.na(ann[[keytype_col]]), c("SYMBOL", keytype_col)]
  if (!is.null(ont_filter)) {
    ann2 <- suppressWarnings(select(org.Hs.eg.db, keys = universe_genes, keytype = "SYMBOL", columns = c("GO", "ONTOLOGY")))
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

universe <- pop$Gene  # os ~31 mil genes do array GSE53868 - universo correto (nao o catalogo Entrez inteiro)

cat("=== GO/KEGG na interseção ESTRITA (", nrow(intersecao_estrita), " genes) ===\n", sep = "")
if (nrow(intersecao_estrita) >= 2) {
  go_a <- run_enrichment(intersecao_estrita$GeneSymbol, universe, "GO", ont_filter = "BP")
  if (!is.null(go_a)) {
    terms <- suppressMessages(select(GO.db, keys = go_a$TERM_ID, keytype = "GOID", columns = "TERM"))
    go_a <- merge(go_a, terms, by.x = "TERM_ID", by.y = "GOID"); go_a <- go_a[order(go_a$pvalue), ]
    write.csv(go_a, "results/06_GO_BP_interseccao_estrita.csv", row.names = FALSE)
    cat("GO BP sig (FDR<0.05):", sum(go_a$p.adjust < 0.05, na.rm = TRUE), "de", nrow(go_a), "\n")
  }
} else cat("Menos de 2 genes -- GO/KEGG não é um teste estatístico válido com N=", nrow(intersecao_estrita), "\n")
cat("\n")

cat("=== GO/KEGG nos genes CONCORDANTES (", nrow(concordantes), " genes) ===\n", sep = "")
go_b <- run_enrichment(concordantes$GeneSymbol, universe, "GO", ont_filter = "BP")
if (!is.null(go_b)) {
  terms <- suppressMessages(select(GO.db, keys = go_b$TERM_ID, keytype = "GOID", columns = "TERM"))
  go_b <- merge(go_b, terms, by.x = "TERM_ID", by.y = "GOID"); go_b <- go_b[order(go_b$pvalue), ]
  write.csv(go_b, "results/06_GO_BP_concordantes.csv", row.names = FALSE)
  cat("GO BP sig (FDR<0.05):", sum(go_b$p.adjust < 0.05, na.rm = TRUE), "de", nrow(go_b), "\n")
  print(head(go_b[, c("TERM", "Count", "pvalue", "p.adjust", "geneID")], 10))
}
kegg_b <- run_enrichment(concordantes$GeneSymbol, universe, "PATH")
if (!is.null(kegg_b)) {
  write.csv(kegg_b, "results/06_KEGG_concordantes.csv", row.names = FALSE)
  cat("\nKEGG sig (FDR<0.05):", sum(kegg_b$p.adjust < 0.05, na.rm = TRUE), "de", nrow(kegg_b), "\n")
}

cat("\n=== FIM ===\n")
