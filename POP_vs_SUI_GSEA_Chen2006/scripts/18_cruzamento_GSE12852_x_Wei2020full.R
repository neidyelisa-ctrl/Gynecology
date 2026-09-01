# =============================================================================
# Script 18 — Cruza o POP do GSE12852 (ligamento uterossacral, ver
# diagnóstico do script 16) com o painel COMPLETO do Wei 2020 (SUI, 6.118
# genes) - mesmo processo dos scripts 06/08/14, agora com um SEGUNDO
# dataset de POP independente do GSE53868.
#
# Como o GSE12852 não tem nenhum DEG individual a FDR<0.05 (dataset menor
# e mais heterogêneo - ver script 16), a interseção estrita
# (DEG-POP ∩ DEG-SUI) não é um teste informativo aqui (H0: 0 DEGs no POP,
# então a interseção com qualquer lista é sempre vazia por construção). O
# teste principal é a CONCORDÂNCIA DE DIREÇÃO usando a direção do t
# moderado de TODOS os genes testados (não exige significância individual)
# - o mesmo raciocínio de GSEA: capta um efeito coordenado mesmo sem genes
# individualmente significativos.
# =============================================================================

suppressMessages({
  library(org.Hs.eg.db)
  library(GO.db)
})

pop <- read.csv("results/GSE12852_uterosacral_limma_completo.csv")
cat("GSE12852 (uterossacral) - genes testados:", nrow(pop), "\n")
cat("DEGs (|log2FC|>1, FDR<0.05):", sum(pop$adj.P.Val < 0.05 & abs(pop$logFC) > 1), "(ver script 16 - 0, dataset pequeno)\n\n")

wei <- read.csv("data/wei2020_mRNA_full.csv")
cat("Wei 2020 (completo) - genes candidatos:", nrow(wei),
    "(", sum(wei$Direction == "up"), "up /", sum(wei$Direction == "down"), "down )\n\n")

cross <- merge(wei[, c("GeneSymbol", "Direction")], pop, by.x = "GeneSymbol", by.y = "Gene", all.x = TRUE)
cross$Testado_no_POP <- !is.na(cross$logFC)
cross$Direcao_POP <- ifelse(is.na(cross$logFC), NA, ifelse(cross$logFC > 0, "up", "down"))
cross$Concordante <- cross$Direction == cross$Direcao_POP
cross <- cross[order(cross$P.Value), ]
write.csv(cross, "results/18_cruzamento_GSE12852_x_Wei2020full_genes.csv", row.names = FALSE)

testaveis <- subset(cross, Testado_no_POP & !is.na(Direcao_POP))
concordantes <- subset(testaveis, Concordante)
cat("=== Concordância de direção (Wei2020-completo x GSE12852-uterossacral) ===\n")
cat(nrow(concordantes), "de", nrow(testaveis), "testáveis =",
    round(100 * nrow(concordantes) / nrow(testaveis), 1), "%\n")
bt <- binom.test(nrow(concordantes), nrow(testaveis), p = 0.5)
cat("H0: concordância = 50% (acaso). p-valor =", format(bt$p.value, digits = 6), "\n\n")

## Genes nominalmente significativos no POP (p bruto<0.01, ainda que não
## sobrevivam FDR - ver script 16, top hits) que também estão no Wei2020
nominal_pop <- subset(cross, Testado_no_POP & P.Value < 0.01)
cat("=== Genes com p bruto<0.01 no GSE12852 E presentes na lista do Wei2020 ===\n")
cat(nrow(nominal_pop), "genes:", paste(nominal_pop$GeneSymbol, collapse = ", "), "\n\n")

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

universe <- pop$Gene
cat("=== GO/KEGG nos", nrow(concordantes), "genes concordantes (Wei2020 x GSE12852) ===\n")
go_b <- run_enrichment(concordantes$GeneSymbol, universe, "GO", ont_filter = "BP")
if (!is.null(go_b)) {
  terms <- suppressMessages(select(GO.db, keys = go_b$TERM_ID, keytype = "GOID", columns = "TERM"))
  go_b <- merge(go_b, terms, by.x = "TERM_ID", by.y = "GOID"); go_b <- go_b[order(go_b$pvalue), ]
  write.csv(go_b, "results/18_GO_BP_concordantes.csv", row.names = FALSE)
  cat("GO BP sig (FDR<0.05):", sum(go_b$p.adjust < 0.05, na.rm = TRUE), "de", nrow(go_b), "\n")
  print(head(go_b[, c("TERM","Count","pvalue","p.adjust","geneID")], 15))
}
kegg_b <- run_enrichment(concordantes$GeneSymbol, universe, "PATH")
if (!is.null(kegg_b)) {
  write.csv(kegg_b, "results/18_KEGG_concordantes.csv", row.names = FALSE)
  cat("\nKEGG sig (FDR<0.05):", sum(kegg_b$p.adjust < 0.05, na.rm = TRUE), "de", nrow(kegg_b), "\n")
  print(head(kegg_b[, c("TERM_ID","Count","pvalue","p.adjust","geneID")], 15))
}

cat("\n=== FIM ===\n")
