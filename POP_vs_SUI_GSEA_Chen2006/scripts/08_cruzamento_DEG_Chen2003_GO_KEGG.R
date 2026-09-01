# =============================================================================
# Script 8 — MESMO PROCESSO do script 06, agora com o SEGUNDO artigo do
# Chen (Chen et al. 2003, Am J Obstet Gynecol 189:89-97): cruzar DEG(POP)
# com a lista de DEG deste artigo (SUI, fase PROLIFERATIVA do ciclo
# menstrual — diferente do Chen 2006, que é fase SECRETÓRIA), olhar genes
# em comum / concordância de direção, e rodar GO/KEGG nesse conjunto.
#
# ARTIGO: Chen B, Wen Y, Zhang Z, Wang H, Warrington JA, Polan ML.
# "Menstrual phase-dependent gene expression differences in periurethral
# vaginal tissue from women with stress incontinence." Am J Obstet Gynecol.
# 2003;189(1):89-97. 5 pares SUI x continentes (os MESMOS grupos do Chen
# 2006, mas amostrados na fase proliferativa, não secretória), array
# Affymetrix HuGeneFL (6800 genes — array mais antigo/menor que o U133A do
# Chen 2006). 90 genes candidatos no artigo (62 up / 28 down); 69 mapeados
# com confiança para símbolo HGNC atual (43 up / 26 down) em
# data/chen2003_90genes.csv (conferido linha a linha contra o PDF original
# antes de rodar este script - ex.: PPIF p=.02358 FC=-3.1, ALOX12
# p=.03810 FC=-2.7, ambos batem exatamente).
# =============================================================================

suppressMessages({
  library(org.Hs.eg.db)
  library(GO.db)
})

pop <- read.csv("results/GSE53868_limma_completo.csv")
pop_deg <- subset(pop, adj.P.Val < 0.05 & abs(logFC) > 1)
cat("DEGs do POP (|log2FC|>1, FDR<0.05):", nrow(pop_deg), "de", nrow(pop), "testados\n\n")

chen03 <- read.csv("data/chen2003_90genes.csv")
chen03 <- chen03[!duplicated(chen03$GeneSymbol), ]
cat("Chen 2003 (SUI, fase proliferativa) — genes candidatos:", nrow(chen03),
    "(", sum(chen03$Direction == "up"), "up /", sum(chen03$Direction == "down"), "down )\n\n")

cross <- merge(chen03[, c("GeneSymbol", "Direction")], pop, by.x = "GeneSymbol", by.y = "Gene", all.x = TRUE)
cross$Testado_no_POP <- !is.na(cross$logFC)
cross$Direcao_POP <- ifelse(is.na(cross$logFC), NA, ifelse(cross$logFC > 0, "up", "down"))
cross$Concordante <- cross$Direction == cross$Direcao_POP
cross$DEG_POP_FDR05_logFC1 <- !is.na(cross$adj.P.Val) & cross$adj.P.Val < 0.05 & abs(cross$logFC) > 1
cross <- cross[order(cross$adj.P.Val), ]
write.csv(cross, "results/08_cruzamento_Chen2003_x_POP_genes.csv", row.names = FALSE)

intersecao_estrita <- subset(cross, DEG_POP_FDR05_logFC1)
cat("=== (a) Interseção ESTRITA: DEG(POP, logFC>1 & FDR<0.05) ∩ Chen2003 ===\n")
cat(nrow(intersecao_estrita), "genes:", paste(intersecao_estrita$GeneSymbol, collapse = ", "), "\n\n")

testaveis <- subset(cross, Testado_no_POP & !is.na(Direcao_POP))
concordantes <- subset(testaveis, Concordante)
cat("=== (b) Genes CONCORDANTES em direção (Chen2003 x POP) ===\n")
cat(nrow(concordantes), "de", nrow(testaveis), "testáveis:",
    paste(concordantes$GeneSymbol, collapse = ", "), "\n\n")

bt <- binom.test(nrow(concordantes), nrow(testaveis), p = 0.5)
cat("=== (c) Teste binomial de sinal ===\n")
cat("Concordância:", nrow(concordantes), "/", nrow(testaveis),
    "=", round(100 * nrow(concordantes) / nrow(testaveis), 1), "%\n")
cat("H0: concordância = 50% (acaso). p-valor =", format(bt$p.value, digits = 4), "\n\n")

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

cat("=== GO/KEGG na interseção ESTRITA (", nrow(intersecao_estrita), " genes) ===\n", sep = "")
if (nrow(intersecao_estrita) >= 2) {
  go_a <- run_enrichment(intersecao_estrita$GeneSymbol, universe, "GO", ont_filter = "BP")
  if (!is.null(go_a)) {
    terms <- suppressMessages(select(GO.db, keys = go_a$TERM_ID, keytype = "GOID", columns = "TERM"))
    go_a <- merge(go_a, terms, by.x = "TERM_ID", by.y = "GOID"); go_a <- go_a[order(go_a$pvalue), ]
    write.csv(go_a, "results/08_GO_BP_interseccao_estrita.csv", row.names = FALSE)
    cat("GO BP sig (FDR<0.05):", sum(go_a$p.adjust < 0.05, na.rm = TRUE), "de", nrow(go_a), "\n")
  }
} else cat("Menos de 2 genes -- GO/KEGG não é um teste estatístico válido com N=", nrow(intersecao_estrita), "\n")
cat("\n")

cat("=== GO/KEGG nos genes CONCORDANTES (", nrow(concordantes), " genes) ===\n", sep = "")
go_b <- run_enrichment(concordantes$GeneSymbol, universe, "GO", ont_filter = "BP")
if (!is.null(go_b)) {
  terms <- suppressMessages(select(GO.db, keys = go_b$TERM_ID, keytype = "GOID", columns = "TERM"))
  go_b <- merge(go_b, terms, by.x = "TERM_ID", by.y = "GOID"); go_b <- go_b[order(go_b$pvalue), ]
  write.csv(go_b, "results/08_GO_BP_concordantes.csv", row.names = FALSE)
  cat("GO BP sig (FDR<0.05):", sum(go_b$p.adjust < 0.05, na.rm = TRUE), "de", nrow(go_b), "\n")
  print(head(go_b[, c("TERM", "Count", "pvalue", "p.adjust", "geneID")], 10))
} else cat("Nenhum termo GO testável (menos de 2 genes concordantes anotados)\n")

kegg_b <- run_enrichment(concordantes$GeneSymbol, universe, "PATH")
if (!is.null(kegg_b)) {
  write.csv(kegg_b, "results/08_KEGG_concordantes.csv", row.names = FALSE)
  cat("\nKEGG sig (FDR<0.05):", sum(kegg_b$p.adjust < 0.05, na.rm = TRUE), "de", nrow(kegg_b), "\n")
  print(head(kegg_b[, c("TERM_ID", "Count", "pvalue", "p.adjust", "geneID")], 10))
}

cat("\n=== FIM ===\n")
