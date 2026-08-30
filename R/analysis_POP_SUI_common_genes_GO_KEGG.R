# =============================================================================
# PIPELINE LIMPO (do jeito que a usuaria pediu, passo a passo):
#   1. DEG de POP (GSE53868, limma, pareado por paciente)
#   2. DEG de SUI (Chen et al. 2006, genes ja significativos no artigo original)
#   3. Genes em COMUM entre os dois (interseccao simples de simbolos de gene)
#      - Criterio: |log2FC|>1 e FDR<0.05 nos dois lados
#      - Se poucos genes em comum, cai para |log2FC|>0.5 e FDR<0.05
#   4. GO e KEGG dos genes em comum
#   5. PPI e hub genes dos genes em comum
#
# Script pronto para rodar na integra no proprio R da usuaria - so precisa de
# data/GSE53868_series_matrix.txt e data/chen2006_79genes.csv (os dois ja
# estao no repositorio GitHub do projeto).
# =============================================================================

suppressMessages({
  library(limma)
  library(org.Hs.eg.db)
  library(GO.db)
  library(igraph)
})

dir.create("results", showWarnings = FALSE)

## =========================================================================
## PASSO 1: DEG de POP (GSE53868)
## =========================================================================

raw_lines <- readLines("data/GSE53868_series_matrix.txt")
start_row <- grep("^!series_matrix_table_begin", raw_lines) + 1
end_row   <- grep("^!series_matrix_table_end", raw_lines) - 1

expr <- read.delim("data/GSE53868_series_matrix.txt", skip = start_row - 1,
                    nrows = end_row - start_row, header = TRUE,
                    row.names = 1, check.names = FALSE, quote = "\"")

sample_titles <- gsub('"', "", strsplit(raw_lines[grep("^!Sample_title", raw_lines)], "\t")[[1]][-1])
individual_line <- raw_lines[grep("^!Sample_characteristics_ch1.*individual:", raw_lines)][1]
individuals <- gsub("individual: ", "", gsub('"', "", strsplit(individual_line, "\t")[[1]][-1]))

tissue <- factor(ifelse(grepl("\\(POP site\\)", sample_titles), "POP_site", "NonPOP_site"),
                  levels = c("NonPOP_site", "POP_site"))
individual <- factor(individuals)
coldata <- data.frame(row.names = colnames(expr), tissue = tissue, individual = individual)

## Desenho pareado (mesma paciente, 2 sitios) - positivo = para cima no
## sitio de POP (o lado afetado)
design <- model.matrix(~ individual + tissue, data = coldata)
fit <- eBayes(lmFit(as.matrix(expr), design))
pop_res <- topTable(fit, coef = "tissuePOP_site", number = Inf, sort.by = "P")
pop_res$Gene <- rownames(pop_res)
write.csv(pop_res, "results/GSE53868_limma_completo.csv", row.names = FALSE)

get_pop_deg <- function(lfc_cut) {
  subset(pop_res, adj.P.Val < 0.05 & abs(logFC) > lfc_cut)
}

## =========================================================================
## PASSO 2: DEG de SUI (Chen et al. 2006) - ja sao os genes significativos
## do artigo original (79 DEGs testados por t-test + correcao de multiplos
## testes, ver metodologia do artigo). Convertendo o fold-change medio
## (escala linear) para log2 para aplicar o mesmo limiar.
## =========================================================================

chen <- read.csv("data/chen2006_79genes.csv")
chen <- chen[!duplicated(chen$GeneSymbol), ]
chen$log2FC <- log2(chen$RMA_FoldChange)

get_sui_deg <- function(lfc_cut) {
  subset(chen, abs(log2FC) > lfc_cut)
}

## =========================================================================
## PASSO 3: Genes em comum - tenta |log2FC|>1 primeiro, cai para 0.5 se
## poucos genes aparecerem
## =========================================================================

try_common <- function(lfc_cut) {
  pop_deg <- get_pop_deg(lfc_cut)
  sui_deg <- get_sui_deg(lfc_cut)
  common <- intersect(pop_deg$Gene, sui_deg$GeneSymbol)
  list(pop_deg = pop_deg, sui_deg = sui_deg, common = common, lfc_cut = lfc_cut)
}

cat("=== Tentativa 1: |log2FC|>1, FDR<0.05 ===\n")
r1 <- try_common(1)
cat("POP DEGs:", nrow(r1$pop_deg), "| SUI DEGs:", nrow(r1$sui_deg),
    "| Genes em comum:", length(r1$common), "\n")

if (length(r1$common) >= 3) {
  final <- r1
  cat("-> Usando limiar |log2FC|>1 (genes em comum suficientes)\n")
} else {
  cat("Poucos genes em comum com |log2FC|>1 - caindo para |log2FC|>0.5\n\n")
  cat("=== Tentativa 2: |log2FC|>0.5, FDR<0.05 ===\n")
  r2 <- try_common(0.5)
  cat("POP DEGs:", nrow(r2$pop_deg), "| SUI DEGs:", nrow(r2$sui_deg),
      "| Genes em comum:", length(r2$common), "\n")
  final <- r2
}

cat("\n=== GENES EM COMUM (limiar final: |log2FC|>", final$lfc_cut, ") ===\n")
common_table <- merge(
  final$sui_deg[, c("GeneSymbol", "Direction", "log2FC")],
  final$pop_deg[, c("Gene", "logFC", "adj.P.Val")],
  by.x = "GeneSymbol", by.y = "Gene"
)
colnames(common_table) <- c("Gene", "SUI_Direction_Chen2006", "SUI_log2FC_Chen2006",
                             "POP_logFC_GSE53868", "POP_padj_GSE53868")
common_table$Mesma_direcao <- sign(common_table$SUI_log2FC_Chen2006) == sign(common_table$POP_logFC_GSE53868)
print(common_table)
write.csv(common_table, "results/POP_x_SUI_genes_comuns_FINAL.csv", row.names = FALSE)

common_genes <- final$common
cat("\nTotal de genes em comum:", length(common_genes), "\n")

## =========================================================================
## PASSO 4: GO e KEGG dos genes em comum
## =========================================================================

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

all_symbols <- keys(org.Hs.eg.db, keytype = "SYMBOL")

cat("\n=== GO Biological Process (genes em comum) ===\n")
go_common <- run_enrichment(common_genes, all_symbols, "GO", ont_filter = "BP")
if (!is.null(go_common)) {
  terms <- suppressMessages(select(GO.db, keys = go_common$TERM_ID, keytype = "GOID", columns = "TERM"))
  go_common <- merge(go_common, terms, by.x = "TERM_ID", by.y = "GOID")
  go_common <- go_common[order(go_common$pvalue), ]
  write.csv(go_common, "results/GO_BP_genes_comuns_FINAL.csv", row.names = FALSE)
  cat("Termos com padj<0.05:", sum(go_common$p.adjust < 0.05, na.rm = TRUE), "de", nrow(go_common), "\n")
  print(head(go_common[, c("TERM", "Count", "pvalue", "p.adjust", "geneID")], 15))
} else cat("Nenhum termo GO encontrado (lista de genes em comum pode ser pequena demais).\n")

cat("\n=== KEGG (genes em comum) ===\n")
kegg_common <- run_enrichment(common_genes, all_symbols, "PATH")
if (!is.null(kegg_common)) {
  write.csv(kegg_common, "results/KEGG_genes_comuns_FINAL.csv", row.names = FALSE)
  cat("Vias com padj<0.05:", sum(kegg_common$p.adjust < 0.05, na.rm = TRUE), "de", nrow(kegg_common), "\n")
  print(head(kegg_common[, c("TERM_ID", "Count", "pvalue", "p.adjust", "geneID")], 15))
} else cat("Nenhuma via KEGG encontrada.\n")

## =========================================================================
## PASSO 5: PPI e hub genes dos genes em comum
##
## STRING ao vivo esta bloqueado neste sandbox (confirmado em varios testes
## ao longo deste projeto). Exportando a lista de genes para a usuaria
## submeter no STRING (string-db.org -> Multiple proteins -> Homo sapiens),
## e tentando tambem montar uma rede a partir de interacoes ja conhecidas
## nas redes STRING que a usuaria ja exportou antes (resultado parcial).
## =========================================================================

writeLines(common_genes, "results/STRING_input_genes_comuns_FINAL.txt")
cat("\nLista de genes em comum salva para submissao no STRING:",
    "results/STRING_input_genes_comuns_FINAL.txt (", length(common_genes), "genes )\n")

## Tenta reaproveitar arestas ja conhecidas das redes STRING exportadas
## anteriormente pela usuaria (POP e SUI 36h/72h) - cobertura parcial, so
## para os genes que ja apareciam nessas listas anteriores.
existing_networks <- c("results/STRING_network_POP.tsv",
                        "results/STRING_network_SUI_36h.tsv",
                        "results/STRING_network_SUI_72h.tsv")
edges_found <- list()
for (net in existing_networks) {
  if (file.exists(net)) {
    d <- read.delim(net)
    d <- d[d[[1]] %in% common_genes & d[[2]] %in% common_genes, ]
    if (nrow(d) > 0) edges_found[[net]] <- d[, 1:2]
  }
}
edges_all <- do.call(rbind, edges_found)

if (!is.null(edges_all) && nrow(edges_all) > 0) {
  colnames(edges_all)[1:2] <- c("from", "to")
  g <- graph_from_data_frame(unique(edges_all), directed = FALSE)
  deg <- degree(g); btw <- betweenness(g)
  hub_table <- data.frame(gene = names(deg), degree = as.integer(deg), betweenness = round(btw, 2))
  hub_table <- hub_table[order(-hub_table$degree), ]
  write.csv(hub_table, "results/hub_genes_comuns_FINAL.csv", row.names = FALSE)
  cat("\n=== Hub genes (a partir de redes STRING ja exportadas antes - cobertura PARCIAL) ===\n")
  cat("Nos:", vcount(g), "| Arestas:", ecount(g), "\n")
  print(hub_table)
  cat("\nATENCAO: essa rede so cobre genes que ja apareciam nas listas ANTERIORES enviadas ao\n",
      "STRING - nao e uma rede completa e nova para ESTA lista especifica de genes em comum.\n",
      "Para hub genes de verdade e completos, exporte results/STRING_input_genes_comuns_FINAL.txt\n",
      "no STRING (Multiple proteins -> Homo sapiens -> Search -> Exports -> rede em TSV) e me\n",
      "envie o arquivo - calculo os hub genes de verdade com igraph a partir dai.\n")
} else {
  cat("\nNenhuma aresta encontrada nas redes STRING ja exportadas para os genes em comum -\n",
      "precisa de uma exportacao NOVA do STRING com a lista atual (arquivo ja gerado acima)\n",
      "para calcular hub genes de verdade.\n")
}
