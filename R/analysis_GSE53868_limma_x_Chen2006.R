# =============================================================================
# GSE53868 (POP, microarray Agilent 4x44K, GPL18142) - analise DEG com limma
# (pacote correto para microarray, DIFERENTE do DESeq2 usado para RNA-seq) +
# cruzamento com Chen et al. 2006 (SUI, 58 genes curados de
# data/chen2006_79genes.csv).
#
# DESENHO DO GSE53868: 12 mulheres com POP, biopsia PAREADA por paciente -
# tecido do SITIO DO PROLAPSO (POP site) vs tecido da mesma paciente fora do
# prolapso (non-POP site, precervical). NAO e POP-vs-controle-saudavel-de-
# outra-pessoa como o GSE208261 - e uma comparacao DENTRO da mesma paciente,
# sitio afetado vs sitio nao afetado. Isso e importante para a interpretacao
# (ver README/chat) - mesmo assim, ambas as direcoes ("POP site" = analogo a
# "mais doente", "non-POP site" = analogo a "mais saudavel" NA MESMA pessoa)
# sao consistentes com "positivo = para cima no lado mais afetado/doente".
#
# Dados ja vem como matrix normalizada (log2), com SYMBOLS de gene como
# ID_REF (nao precisa de anotacao de plataforma separada) - conferido
# diretamente no arquivo antes de escrever este script.
# =============================================================================

suppressMessages({
  library(limma)
  library(org.Hs.eg.db)
  library(GO.db)
})

## -----------------------------------------------------------------------
## 1) Ler o series matrix (formato texto, GEO) - pula o cabecalho de
##    metadados (linhas comecando com "!") e le so a tabela de expressao.
## -----------------------------------------------------------------------
raw_lines <- readLines("data/GSE53868_series_matrix.txt")
start_row <- grep("^!series_matrix_table_begin", raw_lines) + 1
end_row   <- grep("^!series_matrix_table_end", raw_lines) - 1

expr <- read.delim("data/GSE53868_series_matrix.txt", skip = start_row - 1,
                    nrows = end_row - start_row, header = TRUE,
                    row.names = 1, check.names = FALSE, quote = "\"")

cat("Matriz de expressao:", nrow(expr), "genes x", ncol(expr), "amostras\n")

## -----------------------------------------------------------------------
## 2) Metadados das amostras (extraidos do cabecalho do series matrix) -
##    12 pares: tecido do sitio de POP vs tecido do sitio sem POP, mesma
##    paciente ("individual").
## -----------------------------------------------------------------------
sample_title_line <- raw_lines[grep("^!Sample_title", raw_lines)]
sample_titles <- strsplit(sample_title_line, "\t")[[1]][-1]
sample_titles <- gsub('"', "", sample_titles)

individual_line <- raw_lines[grep("^!Sample_characteristics_ch1.*individual:", raw_lines)][1]
individuals <- gsub('"', "", strsplit(individual_line, "\t")[[1]][-1])
individuals <- gsub("individual: ", "", individuals)

tissue <- ifelse(grepl("\\(POP site\\)", sample_titles), "POP_site", "NonPOP_site")
individual <- factor(individuals)
tissue <- factor(tissue, levels = c("NonPOP_site", "POP_site"))

coldata <- data.frame(row.names = colnames(expr), tissue = tissue, individual = individual)
stopifnot(all(rownames(coldata) == colnames(expr)))
print(table(coldata$tissue))

## -----------------------------------------------------------------------
## 3) limma com bloco por paciente (desenho pareado) - equivalente ao
##    paired t-test, mas moderado (shrinkage de variancia, mais robusto
##    com N pequeno). design ~ individual + tissue: positivo em "tissue" =
##    para cima no POP_site (o "lado doente"), negativo = para cima no
##    NonPOP_site (o "lado saudavel" da mesma pessoa).
## -----------------------------------------------------------------------
design <- model.matrix(~ individual + tissue, data = coldata)
fit <- lmFit(as.matrix(expr), design)
fit <- eBayes(fit)

res <- topTable(fit, coef = "tissuePOP_site", number = Inf, sort.by = "P")
res$Gene <- rownames(res)
res <- res[, c("Gene", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val")]

write.csv(res, "results/GSE53868_limma_completo.csv", row.names = FALSE)

sig <- subset(res, adj.P.Val < 0.05 & abs(logFC) > 0.5)
cat("\nGSE53868 - DEGs (padj<0.05, |logFC|>0.5):", nrow(sig), "de", nrow(res), "testados\n")
write.csv(sig, "results/GSE53868_limma_DEG_sig.csv", row.names = FALSE)

## -----------------------------------------------------------------------
## 4) Cruzamento com Chen et al. 2006 (58 genes curados, SUI vs continente)
##    MESMA CONVENCAO DE DIRECAO em ambos: positivo = "para cima no lado
##    afetado/doente" (SUI no Chen; sitio de POP no GSE53868). Nao ha
##    ambiguidade de contraste tipo Treated-vs-Untreated aqui - os dois
##    datasets sao doenca-vs-referencia-saudavel(-analoga), na MESMA
##    orientacao. Confirmado explicitamente abaixo antes de qualquer
##    conclusao de concordancia.
## -----------------------------------------------------------------------
chen <- read.csv("data/chen2006_79genes.csv")
chen_genes <- unique(chen[, c("GeneSymbol", "Direction")])
chen_genes <- chen_genes[!duplicated(chen_genes$GeneSymbol), ]

cross <- merge(chen_genes, res, by.x = "GeneSymbol", by.y = "Gene", all.x = TRUE)
cross$Testado_no_GSE53868 <- !is.na(cross$logFC)
cross$Significativo_GSE53868_FDR05 <- !is.na(cross$adj.P.Val) & cross$adj.P.Val < 0.05 & abs(cross$logFC) > 0.5
cross$Direcao_GSE53868 <- ifelse(is.na(cross$logFC), NA, ifelse(cross$logFC > 0, "up", "down"))
cross$Concordante <- cross$Direction == cross$Direcao_GSE53868
cross <- cross[order(cross$adj.P.Val), ]

write.csv(cross, "results/Chen2006_x_GSE53868.csv", row.names = FALSE)

cat("\n=== Cruzamento Chen 2006 (SUI) x GSE53868 (POP, novo dataset) ===\n")
cat("Genes testados no GSE53868:", sum(cross$Testado_no_GSE53868), "de", nrow(cross), "\n")
cat("Genes SIGNIFICATIVOS no GSE53868 entre os candidatos do Chen 2006:",
    sum(cross$Significativo_GSE53868_FDR05, na.rm = TRUE), "\n")
sig_cross <- subset(cross, Significativo_GSE53868_FDR05)
print(sig_cross[, c("GeneSymbol", "Direction", "logFC", "adj.P.Val", "Direcao_GSE53868", "Concordante")])

testaveis <- subset(cross, Testado_no_GSE53868 & !is.na(Direcao_GSE53868))
n_concord <- sum(testaveis$Concordante, na.rm = TRUE)
n_total <- nrow(testaveis)
cat("\nConcordancia de direcao (Chen 2006 SUI vs GSE53868 POP):", n_concord, "de", n_total, "\n")
bt <- binom.test(n_concord, n_total, p = 0.5)
cat("Teste binomial de sinal: p =", format(bt$p.value, digits = 4), "\n")

## -----------------------------------------------------------------------
## 5) Genes em comum ENTRE OS DOIS DATASETS DE POP (GSE208261 x GSE53868)
##    que TAMBEM estao na lista do Chen 2006 - a evidencia mais forte
##    possivel com os dados que temos: SUI (literatura) + POP em DOIS
##    datasets independentes, todos na mesma direcao.
## -----------------------------------------------------------------------
pop1 <- read.csv("results/Chen2006_79genes_x_POP_completo.csv")  # GSE208261
pop1_sig <- subset(pop1, Significativo_no_POP_FDR05)

triple <- merge(sig_cross[, c("GeneSymbol","Direction","logFC","adj.P.Val","Direcao_GSE53868")],
                 pop1_sig[, c("GeneSymbol","log2FoldChange","padj","Direcao_POP")],
                 by = "GeneSymbol")
cat("\n=== Genes significativos nos DOIS datasets de POP E na lista do Chen 2006 ===\n")
if (nrow(triple) > 0) print(triple) else cat("Nenhum (nenhuma sobreposicao tripla)\n")
write.csv(triple, "results/Chen2006_x_POP_ambos_datasets.csv", row.names = FALSE)

## -----------------------------------------------------------------------
## 6) GO e KEGG offline (mesmo metodo hipergeometrico do resto do projeto)
##    nos DEGs do GSE53868 sozinho
## -----------------------------------------------------------------------
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

universe_53868 <- rownames(expr)
hits_53868 <- sig$Gene
cat("\n=== GO/KEGG do GSE53868 (DEGs sozinhos) ===\n")
go53868 <- run_enrichment(hits_53868, universe_53868, "GO", ont_filter = "BP")
if (!is.null(go53868)) {
  terms <- suppressMessages(select(GO.db, keys = go53868$TERM_ID, keytype = "GOID", columns = "TERM"))
  go53868 <- merge(go53868, terms, by.x = "TERM_ID", by.y = "GOID"); go53868 <- go53868[order(go53868$pvalue), ]
  write.csv(go53868, "results/GO_BP_GSE53868.csv", row.names = FALSE)
  cat("GO BP sig:", sum(go53868$p.adjust < 0.05, na.rm = TRUE), "de", nrow(go53868), "\n")
  print(head(go53868[, c("TERM","Count","pvalue","p.adjust","geneID")], 10))
} else cat("GO: nenhum termo (lista de DEGs pode estar vazia ou pequena demais)\n")

kegg53868 <- run_enrichment(hits_53868, universe_53868, "PATH")
if (!is.null(kegg53868)) {
  write.csv(kegg53868, "results/KEGG_GSE53868.csv", row.names = FALSE)
  cat("\nKEGG sig:", sum(kegg53868$p.adjust < 0.05, na.rm = TRUE), "de", nrow(kegg53868), "\n")
  print(head(kegg53868[, c("TERM_ID","Count","pvalue","p.adjust","geneID")], 10))
} else cat("KEGG: nenhuma via\n")

## -----------------------------------------------------------------------
## 7) GO e KEGG nos genes CRUZADOS (Chen 2006 concordantes com GSE53868)
## -----------------------------------------------------------------------
concordant_genes <- unique(testaveis$GeneSymbol[testaveis$Concordante])
cat("\n=== GO/KEGG nos", length(concordant_genes), "genes concordantes (Chen 2006 x GSE53868) ===\n")
go_cross <- run_enrichment(concordant_genes, keys(org.Hs.eg.db, keytype = "SYMBOL"), "GO", ont_filter = "BP")
if (!is.null(go_cross)) {
  terms <- suppressMessages(select(GO.db, keys = go_cross$TERM_ID, keytype = "GOID", columns = "TERM"))
  go_cross <- merge(go_cross, terms, by.x = "TERM_ID", by.y = "GOID"); go_cross <- go_cross[order(go_cross$pvalue), ]
  write.csv(go_cross, "results/GO_BP_Chen_x_GSE53868_concordantes.csv", row.names = FALSE)
  cat("GO BP sig:", sum(go_cross$p.adjust < 0.05, na.rm = TRUE), "de", nrow(go_cross), "\n")
  print(head(go_cross[, c("TERM","Count","pvalue","p.adjust","geneID")], 10))
}
kegg_cross <- run_enrichment(concordant_genes, keys(org.Hs.eg.db, keytype = "SYMBOL"), "PATH")
if (!is.null(kegg_cross)) {
  write.csv(kegg_cross, "results/KEGG_Chen_x_GSE53868_concordantes.csv", row.names = FALSE)
  cat("\nKEGG sig:", sum(kegg_cross$p.adjust < 0.05, na.rm = TRUE), "de", nrow(kegg_cross), "\n")
  print(head(kegg_cross[, c("TERM_ID","Count","pvalue","p.adjust","geneID")], 10))
}
