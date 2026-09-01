# =============================================================================
# Script 14 — MESMO PROCESSO dos scripts 06/08/12, agora com a Tabela
# Suplementar S2 COMPLETA do Wei et al. 2020 (não mais o top-40 do artigo
# impresso). Este painel SUBSTITUI o dos scripts 12-13 como a referência
# principal do Wei 2020 - muito maior e com p-valor real por gene.
#
# ORIGEM DO ARQUIVO: data/wei2020_TableS2_mRNA_original.xls - a Tabela
# Suplementar S2 do artigo (Electronic Supplementary Material 2, hospedada
# separadamente pela revista, obtida pela usuária diretamente na página do
# artigo). Contém a saída completa do GeneSpring GX (fold change >=2,
# P<0.05): 7.102 linhas (4.615 up / 2.487 down) - bate exatamente com o
# número reportado no resumo do artigo. Cada linha tem P-value, FDR, Fold
# Change E a intensidade normalizada por amostra individual (3 SUI x 3
# Ctrl - o subconjunto de 6 amostras usado para o array, conforme os
# Métodos do artigo: "we compared three samples... to matched tissues
# without SUI using the Arraystar microarray analysis").
#
# data/wei2020_mRNA_full.csv - extraído das duas abas (up_Sui_vs_Ctrl,
# down_Sui_vs_Ctrl), removendo linhas sem GeneSymbol e MÚLTIPLAS SONDAS por
# gene (555 genes com sonda duplicada no "up", 245 no "down") mantendo só a
# sonda de menor p-valor por gene - resultado: 6.118 genes únicos (3.991
# up / 2.127 down).
#
# LIMITAÇÃO QUE CONTINUA VALENDO: esta tabela já vem PRÉ-FILTRADA pelo
# próprio software do estudo original (só os genes que passaram fold
# change>=2 & P<0.05 aparecem aqui) - não é o array inteiro (~20.730 genes
# codificantes testados). Ainda assim, 6.118 genes é uma base MUITO maior
# que os 59-69 dos dois artigos do Chen, e permite testar quase todas as
# vias KEGG (ao contrário dos painéis pequenos, onde a maioria das vias não
# tinha nem 2 genes presentes).
# =============================================================================

suppressMessages({
  library(org.Hs.eg.db)
  library(GO.db)
})

pop <- read.csv("results/GSE53868_limma_completo.csv")
pop_deg <- subset(pop, adj.P.Val < 0.05 & abs(logFC) > 1)
cat("DEGs do POP (|log2FC|>1, FDR<0.05):", nrow(pop_deg), "de", nrow(pop), "testados\n\n")

wei <- read.csv("data/wei2020_mRNA_full.csv")
cat("Wei 2020 (Tabela S2 COMPLETA) — genes candidatos:", nrow(wei),
    "(", sum(wei$Direction == "up"), "up /", sum(wei$Direction == "down"), "down )\n\n")

cross <- merge(wei[, c("GeneSymbol", "Direction")], pop, by.x = "GeneSymbol", by.y = "Gene", all.x = TRUE)
cross$Testado_no_POP <- !is.na(cross$logFC)
cross$Direcao_POP <- ifelse(is.na(cross$logFC), NA, ifelse(cross$logFC > 0, "up", "down"))
cross$Concordante <- cross$Direction == cross$Direcao_POP
cross$DEG_POP_FDR05_logFC1 <- !is.na(cross$adj.P.Val) & cross$adj.P.Val < 0.05 & abs(cross$logFC) > 1
cross <- cross[order(cross$adj.P.Val), ]
write.csv(cross, "results/14_cruzamento_Wei2020full_x_POP_genes.csv", row.names = FALSE)

intersecao_estrita <- subset(cross, DEG_POP_FDR05_logFC1)
cat("=== (a) Interseção ESTRITA: DEG(POP, logFC>1 & FDR<0.05) ∩ Wei2020(completo) ===\n")
cat(nrow(intersecao_estrita), "genes:", paste(intersecao_estrita$GeneSymbol, collapse = ", "), "\n\n")

testaveis <- subset(cross, Testado_no_POP & !is.na(Direcao_POP))
concordantes <- subset(testaveis, Concordante)
cat("=== (b) Genes CONCORDANTES em direção (Wei2020-completo x POP) ===\n")
cat(nrow(concordantes), "de", nrow(testaveis), "testáveis\n\n")

bt <- binom.test(nrow(concordantes), nrow(testaveis), p = 0.5)
cat("=== (c) Teste binomial de sinal ===\n")
cat("Concordância:", nrow(concordantes), "/", nrow(testaveis),
    "=", round(100 * nrow(concordantes) / nrow(testaveis), 1), "%\n")
cat("H0: concordância = 50% (acaso). p-valor =", format(bt$p.value, digits = 6), "\n\n")

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
    write.csv(go_a, "results/14_GO_BP_interseccao_estrita.csv", row.names = FALSE)
    cat("GO BP sig (FDR<0.05):", sum(go_a$p.adjust < 0.05, na.rm = TRUE), "de", nrow(go_a), "\n")
    print(head(go_a[, c("TERM", "Count", "pvalue", "p.adjust", "geneID")], 10))
  }
  kegg_a <- run_enrichment(intersecao_estrita$GeneSymbol, universe, "PATH")
  if (!is.null(kegg_a)) {
    write.csv(kegg_a, "results/14_KEGG_interseccao_estrita.csv", row.names = FALSE)
    cat("\nKEGG sig (FDR<0.05):", sum(kegg_a$p.adjust < 0.05, na.rm = TRUE), "de", nrow(kegg_a), "\n")
    print(head(kegg_a[, c("TERM_ID", "Count", "pvalue", "p.adjust", "geneID")], 10))
  }
} else cat("Menos de 2 genes -- GO/KEGG não é um teste estatístico válido\n")
cat("\n")

cat("=== GO/KEGG nos genes CONCORDANTES (", nrow(concordantes), " genes) ===\n", sep = "")
go_b <- run_enrichment(concordantes$GeneSymbol, universe, "GO", ont_filter = "BP")
if (!is.null(go_b)) {
  terms <- suppressMessages(select(GO.db, keys = go_b$TERM_ID, keytype = "GOID", columns = "TERM"))
  go_b <- merge(go_b, terms, by.x = "TERM_ID", by.y = "GOID"); go_b <- go_b[order(go_b$pvalue), ]
  write.csv(go_b, "results/14_GO_BP_concordantes.csv", row.names = FALSE)
  cat("GO BP sig (FDR<0.05):", sum(go_b$p.adjust < 0.05, na.rm = TRUE), "de", nrow(go_b), "\n")
  print(head(go_b[, c("TERM", "Count", "pvalue", "p.adjust", "geneID")], 15))
}
kegg_b <- run_enrichment(concordantes$GeneSymbol, universe, "PATH")
if (!is.null(kegg_b)) {
  write.csv(kegg_b, "results/14_KEGG_concordantes.csv", row.names = FALSE)
  cat("\nKEGG sig (FDR<0.05):", sum(kegg_b$p.adjust < 0.05, na.rm = TRUE), "de", nrow(kegg_b), "\n")
  print(head(kegg_b[, c("TERM_ID", "Count", "pvalue", "p.adjust", "geneID")], 15))
}

cat("\n=== FIM ===\n")
