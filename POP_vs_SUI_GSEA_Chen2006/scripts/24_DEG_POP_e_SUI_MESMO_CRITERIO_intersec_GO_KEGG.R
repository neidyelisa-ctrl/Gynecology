# =============================================================================
# Script 24 — Pedido explícito: achar os DEG nos DOIS arquivos (GSE53868 e
# Wei2020 Tabela S2), do MESMO jeito e com o MESMO critério nos dois lados,
# achar os genes em COMUM entre os dois grupos de DEG (interseção estrita,
# não concordância de direção), e então rodar GO/KEGG nessa interseção.
#
# POR QUE OS SCRIPTS ANTERIORES (14, 20) NÃO ERAM ISSO: eles usavam a
# classificação de significância já pronta do PRÓPRIO artigo do Wei 2020
# (P<0,05 bruto, sem correção de FDR - o critério do GeneSpring GX do
# estudo original) do lado do SUI, contra |log2FC|>1 e FDR<0,05 (o nosso
# critério, calculado por nós) do lado do POP - dois padrões estatísticos
# DIFERENTES nos dois lados. Isso não está "errado" (é comparar cada
# estudo com o critério que ele próprio usou), mas não é o que foi pedido
# agora: um critério ÚNICO e CONSISTENTE nos dois lados.
#
# O QUE ESTE SCRIPT FAZ DE DIFERENTE: em vez de usar a coluna P-value/FDR
# já pronta da Tabela S2 do Wei 2020, ele usa a INTENSIDADE NORMALIZADA
# POR AMOSTRA (3 SUI x 3 Ctrl) que também está na Tabela S2
# (`data/wei2020_persample_normalized.csv`, extraída da planilha original)
# e roda o MESMO tipo de teste (limma, moderado-t, FDR de
# Benjamini-Hochberg) e o MESMO critério de corte (|log2FC|>1 E FDR<0,05)
# usado no GSE53868 - agora os dois lados são comparáveis de verdade.
#
# LIMITAÇÃO HONESTA QUE PERMANECE (não dá pra evitar com os dados que
# temos): a Tabela S2 do Wei 2020 já vem PRÉ-FILTRADA pelo estudo original
# (só as ~7.102 sondas que já passaram no critério deles estão na planilha
# - as outras ~13.000 sondas testadas no array original, que não foram
# significativas para eles, não estão disponíveis para nós). Isso significa
# que o FDR que calculamos aqui é sobre um universo de 7.102 testes, não
# sobre o array inteiro (~20.730 genes) - um universo menor tende a dar FDR
# mais permissivo que o do array completo. Ainda assim, é a comparação mais
# correta possível com os dados disponíveis, e MUITO mais rigorosa do que
# usar o p-valor bruto do estudo original sem nenhuma correção.
# =============================================================================

suppressMessages({
  library(limma)
  library(org.Hs.eg.db)
  library(GO.db)
})

## =============================================================================
## PARTE 1 — DEG do POP (GSE53868), limma pareado - |log2FC|>1, FDR<0.05
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
cat("DEGs do POP (|log2FC|>1, FDR<0.05):", nrow(pop_deg), "de", nrow(pop_full), "testados\n")
cat(" - para cima no POP:", sum(pop_deg$logFC > 0), " | para baixo:", sum(pop_deg$logFC < 0), "\n\n")

## =============================================================================
## PARTE 2 — DEG do SUI (Wei 2020), MESMO método (limma) e MESMO critério
## =============================================================================
cat("############## PARTE 2: DEG do SUI (Wei 2020), mesmo criterio ##############\n\n")

wei_raw <- read.csv("data/wei2020_persample_normalized.csv")
cat("Sondas na Tabela S2 (com intensidade por amostra):", nrow(wei_raw), "\n")

wei_mat <- as.matrix(wei_raw[, c("Sui1","Sui2","Sui3","Ctrl1","Ctrl2","Ctrl3")])
rownames(wei_mat) <- wei_raw$GeneSymbol

# colapsa sondas duplicadas do mesmo gene (media do log2, metodo padrao limma::avereps)
wei_gene <- avereps(wei_mat, ID = rownames(wei_mat))
cat("Genes unicos apos colapsar sondas duplicadas:", nrow(wei_gene), "\n\n")

group_sui <- factor(c("SUI","SUI","SUI","Ctrl","Ctrl","Ctrl"), levels = c("Ctrl","SUI"))
design_sui <- model.matrix(~group_sui)
fit_sui <- eBayes(lmFit(wei_gene, design_sui))
sui_full <- topTable(fit_sui, coef = "group_suiSUI", number = Inf, sort.by = "P")
sui_full$Gene <- rownames(sui_full)
sui_full <- sui_full[, c("Gene", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val")]
write.csv(sui_full, "results/24_Wei2020_limma_completo_MESMOCRITERIO.csv", row.names = FALSE)

sui_deg <- subset(sui_full, adj.P.Val < 0.05 & abs(logFC) > 1)
sui_deg <- sui_deg[order(sui_deg$adj.P.Val), ]
write.csv(sui_deg, "results/24_Wei2020_DEG_logFC1_FDR05.csv", row.names = FALSE)
cat("DEGs do SUI (|log2FC|>1, FDR<0.05, recalculado por nos com limma):",
    nrow(sui_deg), "de", nrow(sui_full), "testados\n")
cat(" - para cima no SUI:", sum(sui_deg$logFC > 0), " | para baixo:", sum(sui_deg$logFC < 0), "\n\n")

## =============================================================================
## PARTE 3 — Genes em COMUM entre os dois grupos de DEG (interseção estrita)
## =============================================================================
cat("############## PARTE 3: Genes em comum (intersecao estrita) ##############\n\n")

genes_comuns <- intersect(pop_deg$Gene, sui_deg$Gene)
cat("DEGs do POP:", nrow(pop_deg), " | DEGs do SUI:", nrow(sui_deg),
    " | EM COMUM:", length(genes_comuns), "\n\n")

if (length(genes_comuns) > 0) {
  tab_comuns <- merge(pop_deg[, c("Gene","logFC","adj.P.Val")],
                       sui_deg[, c("Gene","logFC","adj.P.Val")],
                       by = "Gene", suffixes = c("_POP", "_SUI"))
  tab_comuns$Direcao_POP <- ifelse(tab_comuns$logFC_POP > 0, "up", "down")
  tab_comuns$Direcao_SUI <- ifelse(tab_comuns$logFC_SUI > 0, "up", "down")
  tab_comuns$Mesma_direcao <- tab_comuns$Direcao_POP == tab_comuns$Direcao_SUI
  tab_comuns <- tab_comuns[order(tab_comuns$adj.P.Val_POP), ]
  write.csv(tab_comuns, "results/24_genes_em_comum_POP_x_SUI.csv", row.names = FALSE)
  cat("Genes em comum:\n")
  print(tab_comuns[, c("Gene","logFC_POP","logFC_SUI","Direcao_POP","Direcao_SUI","Mesma_direcao")])
  cat("\nMesma direção:", sum(tab_comuns$Mesma_direcao), "de", nrow(tab_comuns), "\n\n")
} else {
  cat("Nenhum gene em comum entre os dois grupos de DEG.\n\n")
}

## =============================================================================
## PARTE 4 — GO e KEGG na interseção (universo = genes testados nos DOIS)
## =============================================================================
cat("############## PARTE 4: GO/KEGG na intersecao ##############\n\n")

universo_comum <- intersect(pop_full$Gene, sui_full$Gene)
cat("Universo (genes testados nos dois arrays):", length(universo_comum), "\n\n")

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
    write.csv(go_res, "results/24_GO_BP_intersecao.csv", row.names = FALSE)
    cat("GO BP sig (FDR<0.05):", sum(go_res$p.adjust < 0.05, na.rm = TRUE), "de", nrow(go_res), "\n")
    print(head(go_res[, c("TERM","Count","pvalue","p.adjust","geneID")], 15))
  } else cat("Nenhum termo GO testavel.\n")

  kegg_res <- run_enrichment(genes_comuns, universo_comum, "PATH")
  if (!is.null(kegg_res)) {
    write.csv(kegg_res, "results/24_KEGG_intersecao.csv", row.names = FALSE)
    cat("\nKEGG sig (FDR<0.05):", sum(kegg_res$p.adjust < 0.05, na.rm = TRUE), "de", nrow(kegg_res), "\n")
    print(head(kegg_res[, c("TERM_ID","Count","pvalue","p.adjust","geneID")], 15))
  } else cat("Nenhuma via KEGG testavel.\n")
} else {
  cat("Menos de 2 genes em comum -- GO/KEGG nao e um teste estatistico valido com N=",
      length(genes_comuns), "(ver PARTE 3 acima).\n")
  cat("Alternativa se isso acontecer: ver results/24_genes_em_comum_POP_x_SUI.csv\n")
  cat("com um limiar mais frouxo (ex.: |log2FC|>0.5) so no lado do POP e/ou do SUI,\n")
  cat("ou usar o teste de CONCORDANCIA DE DIRECAO (scripts 14/20), que nao exige\n")
  cat("significancia individual nos dois lados.\n")
}

cat("\n=== FIM (script 24) ===\n")
