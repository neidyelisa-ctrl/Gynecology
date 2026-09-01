# =============================================================================
# Script 12 — MESMO PROCESSO dos scripts 06/08, agora com o TERCEIRO painel
# de SUI: Wei et al. 2020 (Reproductive Sciences 27:1490-1501), mRNA.
#
# ATENÇÃO - LIMITAÇÃO IMPORTANTE DESTE PAINEL (leia antes de interpretar):
# O artigo original identificou 7.102 mRNAs diferencialmente expressos
# (fold change >=2, P<0.05) - MUITO mais que os painéis do Chen. Mas o PDF
# fornecido só traz as Tabelas 5 e 6 do artigo (top 20 up-regulados e top 20
# down-regulados, só com nome do gene e fold change, SEM p-valor individual
# por gene) - a tabela completa de 7.102 genes está no Material Suplementar
# S2, hospedado à parte no site da revista, não incluído neste PDF. Ou seja:
# este painel (40 genes) é, na prática, MENOR que os do Chen (59-69 genes),
# não maior - o "7.102" do resumo do artigo não está disponível para nós
# ainda. Se a usuária conseguir baixar a Tabela Suplementar S2 (link
# "Electronic supplementary material" na página do artigo), o ganho real de
# poder estatístico está lá, não neste PDF.
#
# Como as Tabelas 5/6 não trazem p-valor por gene (só fold change), o
# escore de ranking do script 13 usa log2(fold change) diretamente (não
# -log10(p) como nos scripts 03/09) - é uma métrica de ranking válida e
# padrão em GSEA preranked quando só se tem fold change.
#
# Nota: o gene ZMAT4 aparece 2x na Tabela 6 original (duas sondas/
# transcritos diferentes, fold change 365,95 e 269,48) - mantida só a
# primeira ocorrência em data/wei2020_top40_mRNA.csv; o script já removeria
# a duplicata de qualquer forma (`!duplicated(GeneSymbol)`).
# =============================================================================

suppressMessages({
  library(org.Hs.eg.db)
  library(GO.db)
})

pop <- read.csv("results/GSE53868_limma_completo.csv")
pop_deg <- subset(pop, adj.P.Val < 0.05 & abs(logFC) > 1)
cat("DEGs do POP (|log2FC|>1, FDR<0.05):", nrow(pop_deg), "de", nrow(pop), "testados\n\n")

wei <- read.csv("data/wei2020_top40_mRNA.csv")
wei <- wei[!duplicated(wei$GeneSymbol), ]
cat("Wei 2020 (SUI, top mRNAs do artigo) — genes candidatos:", nrow(wei),
    "(", sum(wei$Direction == "up"), "up /", sum(wei$Direction == "down"), "down )\n\n")

cross <- merge(wei[, c("GeneSymbol", "Direction")], pop, by.x = "GeneSymbol", by.y = "Gene", all.x = TRUE)
cross$Testado_no_POP <- !is.na(cross$logFC)
cross$Direcao_POP <- ifelse(is.na(cross$logFC), NA, ifelse(cross$logFC > 0, "up", "down"))
cross$Concordante <- cross$Direction == cross$Direcao_POP
cross$DEG_POP_FDR05_logFC1 <- !is.na(cross$adj.P.Val) & cross$adj.P.Val < 0.05 & abs(cross$logFC) > 1
cross <- cross[order(cross$adj.P.Val), ]
write.csv(cross, "results/12_cruzamento_Wei2020_x_POP_genes.csv", row.names = FALSE)

intersecao_estrita <- subset(cross, DEG_POP_FDR05_logFC1)
cat("=== (a) Interseção ESTRITA: DEG(POP, logFC>1 & FDR<0.05) ∩ Wei2020 ===\n")
cat(nrow(intersecao_estrita), "genes:", paste(intersecao_estrita$GeneSymbol, collapse = ", "), "\n\n")

testaveis <- subset(cross, Testado_no_POP & !is.na(Direcao_POP))
concordantes <- subset(testaveis, Concordante)
cat("=== (b) Genes CONCORDANTES em direção (Wei2020 x POP) ===\n")
cat(nrow(concordantes), "de", nrow(testaveis), "testáveis:",
    paste(concordantes$GeneSymbol, collapse = ", "), "\n\n")

if (nrow(testaveis) > 0) {
  bt <- binom.test(nrow(concordantes), nrow(testaveis), p = 0.5)
  cat("=== (c) Teste binomial de sinal ===\n")
  cat("Concordância:", nrow(concordantes), "/", nrow(testaveis),
      "=", round(100 * nrow(concordantes) / nrow(testaveis), 1), "%\n")
  cat("H0: concordância = 50% (acaso). p-valor =", format(bt$p.value, digits = 4), "\n\n")
}

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
    write.csv(go_a, "results/12_GO_BP_interseccao_estrita.csv", row.names = FALSE)
    cat("GO BP sig (FDR<0.05):", sum(go_a$p.adjust < 0.05, na.rm = TRUE), "de", nrow(go_a), "\n")
  }
} else cat("Menos de 2 genes -- GO/KEGG não é um teste estatístico válido com N=", nrow(intersecao_estrita), "\n")
cat("\n")

cat("=== GO/KEGG nos genes CONCORDANTES (", nrow(concordantes), " genes) ===\n", sep = "")
if (nrow(concordantes) >= 2) {
  go_b <- run_enrichment(concordantes$GeneSymbol, universe, "GO", ont_filter = "BP")
  if (!is.null(go_b)) {
    terms <- suppressMessages(select(GO.db, keys = go_b$TERM_ID, keytype = "GOID", columns = "TERM"))
    go_b <- merge(go_b, terms, by.x = "TERM_ID", by.y = "GOID"); go_b <- go_b[order(go_b$pvalue), ]
    write.csv(go_b, "results/12_GO_BP_concordantes.csv", row.names = FALSE)
    cat("GO BP sig (FDR<0.05):", sum(go_b$p.adjust < 0.05, na.rm = TRUE), "de", nrow(go_b), "\n")
    print(head(go_b[, c("TERM", "Count", "pvalue", "p.adjust", "geneID")], 10))
  } else cat("Nenhum termo GO testável\n")
  kegg_b <- run_enrichment(concordantes$GeneSymbol, universe, "PATH")
  if (!is.null(kegg_b)) {
    write.csv(kegg_b, "results/12_KEGG_concordantes.csv", row.names = FALSE)
    cat("\nKEGG sig (FDR<0.05):", sum(kegg_b$p.adjust < 0.05, na.rm = TRUE), "de", nrow(kegg_b), "\n")
    print(head(kegg_b[, c("TERM_ID", "Count", "pvalue", "p.adjust", "geneID")], 10))
  }
} else cat("Menos de 2 genes concordantes -- GO/KEGG não é um teste válido\n")

cat("\n=== FIM ===\n")
