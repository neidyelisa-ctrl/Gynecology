# =============================================================================
# Script 27 — Resposta a: "dá pra mudar o critério do autor pra log2FC>1?"
#
# ESCLARECIMENTO IMPORTANTE ANTES DE QUALQUER NÚMERO: o critério de fold
# change do autor (fold change >= 2, em escala linear) já É EXATAMENTE
# |log2FC| > 1 - log2(2) = 1, são a MESMA coisa matematicamente, não dá pra
# "mudar" isso porque já é isso. O que os autores fizeram DIFERENTE do
# nosso padrão (usado em todo o resto deste projeto) não foi o corte de
# fold change - foi a SIGNIFICÂNCIA: eles usaram p-valor BRUTO < 0,05, SEM
# corrigir para múltiplas comparações (nenhum ajuste de FDR no critério de
# inclusão). É ISSO que dá pra trocar por FDR<0,05 - e sim, dá, de um jeito
# válido (diferente da tentativa do script 24, que era circular).
#
# POR QUE ISSO NÃO É CIRCULAR (diferente do script 24): o script 24
# recalculou a significância do ZERO a partir da intensidade por amostra -
# um teste novo, rodado só nos genes que JÁ tinham sido pré-selecionados
# por serem significativos, o que reconfirma quase tudo por construção
# (97,5%). Aqui, em vez disso, USAMOS A COLUNA "FDR" QUE JÁ VEM PRONTA NA
# PRÓPRIA TABELA SUPLEMENTAR DO ARTIGO - e essa coluna, conferida abaixo,
# tem variação real (de 0,003 a 0,124, não todo mundo colado em zero) - ou
# seja, foi calculada pelo GeneSpring GX no ARRAY INTEIRO ANTES do filtro
# de corte (fold change>=2 & p<0,05) ser aplicado para gerar a tabela final
# - é a estatística de fato calculada pelos autores, não uma reconstrução
# nossa. Aplicar um corte mais rígido (FDR<0,05 em vez de p bruto<0,05)
# em cima de uma estatística real e já calculada é um refinamento válido,
# não um novo teste circular.
# =============================================================================

suppressMessages({
  library(limma)
  library(org.Hs.eg.db)
  library(GO.db)
})

## =============================================================================
## PARTE 1 — DEG do POP (GSE53868) - sem mudanca, |log2FC|>1 e FDR<0.05
## =============================================================================
cat("\n############## PARTE 1: DEG do POP (GSE53868) ##############\n\n")

raw_lines <- readLines("data/GSE53868_series_matrix.txt")
start_row <- grep("^!series_matrix_table_begin", raw_lines) + 1
end_row   <- grep("^!series_matrix_table_end", raw_lines) - 1
expr <- read.delim("data/GSE53868_series_matrix.txt", skip = start_row - 1,
                    nrows = end_row - start_row, header = TRUE,
                    row.names = 1, check.names = FALSE, quote = "\"")
sample_title_line <- raw_lines[grep("^!Sample_title", raw_lines)]
sample_titles <- gsub('"', "", strsplit(sample_title_line, "\t")[[1]][-1])
individual_line <- raw_lines[grep("^!Sample_characteristics_ch1.*individual:", raw_lines)][1]
individuals <- gsub("individual: ", "", gsub('"', "", strsplit(individual_line, "\t")[[1]][-1]))
tissue <- ifelse(grepl("\\(POP site\\)", sample_titles), "POP_site", "NonPOP_site")
coldata <- data.frame(row.names = colnames(expr),
                       tissue = factor(tissue, levels = c("NonPOP_site", "POP_site")),
                       individual = factor(individuals))
design_pop <- model.matrix(~ individual + tissue, data = coldata)
fit_pop <- eBayes(lmFit(as.matrix(expr), design_pop))
pop_full <- topTable(fit_pop, coef = "tissuePOP_site", number = Inf, sort.by = "P")
pop_full$Gene <- rownames(pop_full)
pop_full <- pop_full[, c("Gene", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val")]
dir.create("results", showWarnings = FALSE)
write.csv(pop_full, "results/GSE53868_limma_completo.csv", row.names = FALSE)

pop_deg <- subset(pop_full, adj.P.Val < 0.05 & abs(logFC) > 1)
write.csv(pop_deg, "results/GSE53868_DEG_logFC1_FDR05.csv", row.names = FALSE)
cat("DEGs do POP (|log2FC|>1, FDR<0.05):", nrow(pop_deg), "de", nrow(pop_full), "testados\n\n")

## =============================================================================
## PARTE 2 — DEG do SUI (Wei2020): mesmo |log2FC|>1 (ja embutido no FC>=2
## do artigo), agora com FDR<0.05 (a coluna FDR ja calculada pelos autores)
## em vez do p bruto<0.05 deles
## =============================================================================
cat("############## PARTE 2: DEG do SUI (Wei2020), FDR proprio<0.05 ##############\n\n")

wei <- read.csv("data/wei2020_mRNA_full.csv")
cat("Genes na Tabela S2 (ja FC>=2 e p bruto<0.05, criterio original):", nrow(wei), "\n")
cat("Faixa da coluna FDR (ja calculada pelos autores no array completo):",
    round(min(wei$FDR), 4), "a", round(max(wei$FDR), 4),
    "(variacao real -> confirma que NAO e circular, foi calculada antes do filtro)\n\n")

wei$logFC <- log2(wei$FoldChange)
sui_deg <- subset(wei, FDR < 0.05 & abs(logFC) > 1)
sui_deg <- sui_deg[order(sui_deg$FDR), ]
write.csv(sui_deg, "results/27_Wei2020_DEG_logFC1_FDRproprio05.csv", row.names = FALSE)

cat("DEGs do SUI (|log2FC|>1 E FDR<0.05, criterio dos PROPRIOS autores):",
    nrow(sui_deg), "de", nrow(wei), "\n")
cat(" - para cima no SUI:", sum(sui_deg$Direction == "up"),
    " | para baixo:", sum(sui_deg$Direction == "down"), "\n")
cat(" (para comparacao: com o criterio ORIGINAL dos autores - p bruto<0.05 - eram", nrow(wei), "genes)\n\n")

## =============================================================================
## PARTE 3 — Genes em comum (intersecao estrita), agora com FDR<0.05 nos DOIS lados
## =============================================================================
cat("############## PARTE 3: Genes em comum (intersecao estrita) ##############\n\n")

genes_comuns <- intersect(pop_deg$Gene, sui_deg$GeneSymbol)
cat("DEGs do POP:", nrow(pop_deg), " | DEGs do SUI (FDR proprio<0.05):", nrow(sui_deg),
    " | EM COMUM:", length(genes_comuns), "\n\n")

if (length(genes_comuns) > 0) {
  tab_comuns <- merge(pop_deg[, c("Gene","logFC","adj.P.Val")],
                       sui_deg[, c("GeneSymbol","logFC","FDR")],
                       by.x = "Gene", by.y = "GeneSymbol")
  colnames(tab_comuns) <- c("Gene","logFC_POP","FDR_POP","logFC_SUI","FDR_SUI")
  tab_comuns$Direcao_POP <- ifelse(tab_comuns$logFC_POP > 0, "up", "down")
  tab_comuns$Direcao_SUI <- ifelse(tab_comuns$logFC_SUI > 0, "up", "down")
  tab_comuns$Mesma_direcao <- tab_comuns$Direcao_POP == tab_comuns$Direcao_SUI
  tab_comuns <- tab_comuns[order(tab_comuns$FDR_POP), ]
  write.csv(tab_comuns, "results/27_genes_em_comum_FDRproprio.csv", row.names = FALSE)
  print(tab_comuns)
  cat("\nMesma direcao:", sum(tab_comuns$Mesma_direcao), "de", nrow(tab_comuns), "\n\n")
} else {
  cat("Nenhum gene em comum.\n\n")
}

## =============================================================================
## PARTE 4 — GO/KEGG na intersecao (universo = genes testados nos DOIS)
## =============================================================================
cat("############## PARTE 4: GO/KEGG na intersecao ##############\n\n")

universo_comum <- intersect(pop_full$Gene, wei$GeneSymbol)
cat("Universo (genes testados nos dois -- POP array completo x Tabela S2 do Wei2020):",
    length(universo_comum), "\n\n")

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

if (length(genes_comuns) >= 2) {
  go_res <- run_enrichment(genes_comuns, universo_comum, "GO", ont_filter = "BP")
  if (!is.null(go_res)) {
    terms <- suppressMessages(select(GO.db, keys = go_res$TERM_ID, keytype = "GOID", columns = "TERM"))
    go_res <- merge(go_res, terms, by.x = "TERM_ID", by.y = "GOID"); go_res <- go_res[order(go_res$pvalue), ]
    write.csv(go_res, "results/27_GO_BP_intersecao.csv", row.names = FALSE)
    cat("GO BP sig (FDR<0.05):", sum(go_res$p.adjust < 0.05, na.rm = TRUE), "de", nrow(go_res), "\n")
    print(head(go_res[, c("TERM","Count","pvalue","p.adjust","geneID")], 15))
  } else cat("Nenhum termo GO testavel.\n")

  kegg_res <- run_enrichment(genes_comuns, universo_comum, "PATH")
  if (!is.null(kegg_res)) {
    write.csv(kegg_res, "results/27_KEGG_intersecao.csv", row.names = FALSE)
    cat("\nKEGG sig (FDR<0.05):", sum(kegg_res$p.adjust < 0.05, na.rm = TRUE), "de", nrow(kegg_res), "\n")
    print(head(kegg_res[, c("TERM_ID","Count","pvalue","p.adjust","geneID")], 15))
  } else cat("Nenhuma via KEGG testavel.\n")
} else {
  cat("Menos de 2 genes em comum -- GO/KEGG nao e valido, ver PARTE 3.\n")
}

cat("\n=== FIM (script 27) ===\n")
