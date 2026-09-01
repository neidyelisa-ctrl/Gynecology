# =============================================================================
# Script 22 — Mesmo processo do script 21, agora Chen 2003 (SUI, fase
# proliferativa) x GSE12852 (POP, ligamento uterossacral). Completa a
# matriz 2 POP x 3 SUI (ver cabeçalho do script 21).
# =============================================================================

suppressMessages({
  library(org.Hs.eg.db)
  library(GO.db)
})

pop <- read.csv("results/GSE12852_uterosacral_limma_completo.csv")
cat("GSE12852 (uterossacral) - genes testados:", nrow(pop), "\n\n")

chen03 <- read.csv("data/chen2003_90genes.csv")
chen03 <- chen03[!duplicated(chen03$GeneSymbol), ]
cat("Chen 2003 - genes candidatos:", nrow(chen03), "\n\n")

cross <- merge(chen03[, c("GeneSymbol", "Direction")], pop, by.x = "GeneSymbol", by.y = "Gene", all.x = TRUE)
cross$Testado_no_POP <- !is.na(cross$logFC)
cross$Direcao_POP <- ifelse(is.na(cross$logFC), NA, ifelse(cross$logFC > 0, "up", "down"))
cross$Concordante <- cross$Direction == cross$Direcao_POP
cross$DEG_POP_FDR05_logFC1 <- !is.na(cross$adj.P.Val) & cross$adj.P.Val < 0.05 & abs(cross$logFC) > 1
write.csv(cross, "results/22_cruzamento_Chen2003_x_GSE12852_genes.csv", row.names = FALSE)

intersecao_estrita <- subset(cross, DEG_POP_FDR05_logFC1)
cat("=== Interseção ESTRITA: DEG(GSE12852) ∩ Chen2003 ===\n")
cat(nrow(intersecao_estrita), "genes\n\n")

testaveis <- subset(cross, Testado_no_POP & !is.na(Direcao_POP))
concordantes <- subset(testaveis, Concordante)
cat("=== Concordância de direção (Chen2003 x GSE12852) ===\n")
cat(nrow(concordantes), "de", nrow(testaveis), "=",
    round(100 * nrow(concordantes) / nrow(testaveis), 1), "%\n")
bt <- binom.test(nrow(concordantes), nrow(testaveis), p = 0.5)
cat("p-valor =", format(bt$p.value, digits = 4), "\n\n")

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
cat("=== GO/KEGG nos", nrow(concordantes), "genes concordantes (Chen2003 x GSE12852) ===\n")
go_b <- run_enrichment(concordantes$GeneSymbol, universe, "GO", ont_filter = "BP")
if (!is.null(go_b)) {
  terms <- suppressMessages(select(GO.db, keys = go_b$TERM_ID, keytype = "GOID", columns = "TERM"))
  go_b <- merge(go_b, terms, by.x = "TERM_ID", by.y = "GOID"); go_b <- go_b[order(go_b$pvalue), ]
  write.csv(go_b, "results/22_GO_BP_concordantes.csv", row.names = FALSE)
  cat("GO BP sig (FDR<0.05):", sum(go_b$p.adjust < 0.05, na.rm = TRUE), "de", nrow(go_b), "\n")
  print(head(go_b[, c("TERM","Count","pvalue","p.adjust","geneID")], 10))
} else cat("Nenhum termo GO testável (poucos genes concordantes)\n")

kegg_b <- run_enrichment(concordantes$GeneSymbol, universe, "PATH")
if (!is.null(kegg_b)) {
  write.csv(kegg_b, "results/22_KEGG_concordantes.csv", row.names = FALSE)
  cat("\nKEGG sig (FDR<0.05):", sum(kegg_b$p.adjust < 0.05, na.rm = TRUE), "de", nrow(kegg_b), "\n")
  print(head(kegg_b[, c("TERM_ID","Count","pvalue","p.adjust","geneID")], 10))
}

cat("\n=== FIM ===\n")
