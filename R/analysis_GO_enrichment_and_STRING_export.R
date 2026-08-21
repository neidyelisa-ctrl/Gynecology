# =============================================================================
# Enriquecimento GO (offline) para SUI e POP, e exportacao das listas de genes
# para o STRING (necessarias para a etapa de hub genes / rede PPI e KEGG,
# que dependem do site string-db.org - nao acessivel a partir deste sandbox).
#
# ESTE SCRIPT REPRODUZ A ANALISE REAL rodada em 20/08/2026: ver
# results/GO_BP_SUI_36h.csv e results/GO_BP_POP.csv para os resultados, e
# results/STRING_input_SUI_36h_genes.txt / STRING_input_POP_genes.txt para
# as listas que devem ser coladas em https://string-db.org (Multiple
# proteins -> cola a lista -> Homo sapiens -> Search).
#
# O enriquecimento GO abaixo NAO usa clusterProfiler - e feito manualmente
# com org.Hs.eg.db + GO.db + teste hipergeometrico (phyper) + correcao BH,
# que e equivalente ao que enrichGO faz por baixo dos panos. Isso evita a
# dependencia de instalar clusterProfiler via Bioconductor.
# =============================================================================

if (!requireNamespace("org.Hs.eg.db", quietly = TRUE) || !requireNamespace("GO.db", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install(c("org.Hs.eg.db", "GO.db"), update = FALSE, ask = FALSE)
}
if (!requireNamespace("homologene", quietly = TRUE)) install.packages("homologene")

suppressMessages({
  library(org.Hs.eg.db)
  library(GO.db)
  library(AnnotationDbi)
})

dir.create("results", showWarnings = FALSE)

## -----------------------------------------------------------------------
## Funcao de enriquecimento GO via teste hipergeometrico + BH
## -----------------------------------------------------------------------
run_go_enrichment <- function(hit_genes, universe_genes, ont = "BP", min_gs = 5, max_gs = 500) {
  hit_genes <- unique(intersect(hit_genes, universe_genes))
  universe_genes <- unique(universe_genes)

  ann <- suppressWarnings(select(org.Hs.eg.db, keys = universe_genes, keytype = "SYMBOL", columns = "GO"))
  ann <- ann[!is.na(ann$GO) & ann$ONTOLOGY == ont, c("SYMBOL", "GO")]
  ann <- unique(ann)

  gs_size <- table(ann$GO)
  valid_go <- names(gs_size)[gs_size >= min_gs & gs_size <= max_gs]
  ann <- ann[ann$GO %in% valid_go, ]

  N <- length(universe_genes)
  n <- length(hit_genes)

  res <- lapply(valid_go, function(go_id) {
    genes_in_term <- unique(ann$SYMBOL[ann$GO == go_id])
    M <- length(genes_in_term)
    hits_in_term <- intersect(genes_in_term, hit_genes)
    k <- length(hits_in_term)
    if (k == 0) return(NULL)
    pval <- phyper(k - 1, M, N - M, n, lower.tail = FALSE)
    data.frame(GO_ID = go_id, GeneRatio = paste0(k, "/", n), BgRatio = paste0(M, "/", N),
               Count = k, pvalue = pval, geneID = paste(hits_in_term, collapse = "/"),
               stringsAsFactors = FALSE)
  })
  res <- do.call(rbind, res)
  if (is.null(res) || nrow(res) == 0) return(res)

  terms <- suppressMessages(select(GO.db, keys = res$GO_ID, keytype = "GOID", columns = "TERM"))
  res <- merge(res, terms, by.x = "GO_ID", by.y = "GOID")
  res$p.adjust <- p.adjust(res$pvalue, method = "BH")
  res <- res[order(res$pvalue), ]
  res[, c("GO_ID", "TERM", "GeneRatio", "BgRatio", "Count", "pvalue", "p.adjust", "geneID")]
}

## -----------------------------------------------------------------------
## SUI 36h
## -----------------------------------------------------------------------
sui_full <- read.csv("results/SUI_36hr_DESeq2_completo.csv")
sui72_full <- read.csv("results/SUI_72hr_DESeq2_completo.csv")
all_rat_genes <- unique(c(sui_full$Rat_Gene_Symbol, sui72_full$Rat_Gene_Symbol))
orth <- homologene::homologene(all_rat_genes, inTax = 10116, outTax = 9606)
colnames(orth)[1:2] <- c("Rat_Gene_Symbol", "Human_Ortholog_Symbol")
sui_full_orth <- merge(sui_full, orth[, c("Rat_Gene_Symbol", "Human_Ortholog_Symbol")],
                        by = "Rat_Gene_Symbol", all.x = TRUE)

sui_universe <- unique(na.omit(sui_full_orth$Human_Ortholog_Symbol))
sui_hits <- unique(na.omit(sui_full_orth$Human_Ortholog_Symbol[
  !is.na(sui_full_orth$padj) & sui_full_orth$padj < 0.05 & abs(sui_full_orth$log2FoldChange) > 0.5]))

cat("SUI 36h - universo (com ortologo):", length(sui_universe), "| DEGs:", length(sui_hits), "\n")
go_sui <- run_go_enrichment(sui_hits, sui_universe, ont = "BP")
write.csv(go_sui, "results/GO_BP_SUI_36h.csv", row.names = FALSE)
writeLines(sui_hits, "results/STRING_input_SUI_36h_genes.txt")

## -----------------------------------------------------------------------
## SUI 72h
## -----------------------------------------------------------------------
sui72_full_orth <- merge(sui72_full, orth[, c("Rat_Gene_Symbol", "Human_Ortholog_Symbol")],
                          by = "Rat_Gene_Symbol", all.x = TRUE)

sui72_universe <- unique(na.omit(sui72_full_orth$Human_Ortholog_Symbol))
sui72_hits <- unique(na.omit(sui72_full_orth$Human_Ortholog_Symbol[
  !is.na(sui72_full_orth$padj) & sui72_full_orth$padj < 0.05 & abs(sui72_full_orth$log2FoldChange) > 0.5]))

cat("SUI 72h - universo (com ortologo):", length(sui72_universe), "| DEGs:", length(sui72_hits), "\n")
go_sui72 <- run_go_enrichment(sui72_hits, sui72_universe, ont = "BP")
write.csv(go_sui72, "results/GO_BP_SUI_72h.csv", row.names = FALSE)
writeLines(sui72_hits, "results/STRING_input_SUI_72h_genes.txt")

cat("\nTop 15 termos GO BP - SUI 72h (por p-valor bruto):\n")
print(head(go_sui72[, c("TERM", "Count", "pvalue", "p.adjust")], 15))

## -----------------------------------------------------------------------
## POP
## -----------------------------------------------------------------------
pop_full <- read.csv("results/POP_GSE208261_DESeq2_completo.csv")
pop_full$Human_Entrez_ID <- as.character(pop_full$Human_Entrez_ID)
pop_map_all <- select(org.Hs.eg.db, keys = pop_full$Human_Entrez_ID, keytype = "ENTREZID", columns = "SYMBOL")
pop_full_sym <- merge(pop_full, pop_map_all, by.x = "Human_Entrez_ID", by.y = "ENTREZID", all.x = TRUE)

pop_universe <- unique(na.omit(pop_full_sym$SYMBOL))
pop_hits <- unique(na.omit(pop_full_sym$SYMBOL[
  !is.na(pop_full_sym$padj) & pop_full_sym$padj < 0.05 & abs(pop_full_sym$log2FoldChange) > 0.5]))

cat("POP - universo:", length(pop_universe), "| DEGs:", length(pop_hits), "\n")
go_pop <- run_go_enrichment(pop_hits, pop_universe, ont = "BP")
write.csv(go_pop, "results/GO_BP_POP.csv", row.names = FALSE)
writeLines(pop_hits, "results/STRING_input_POP_genes.txt")

## -----------------------------------------------------------------------
## Comparacao SUI x POP
## -----------------------------------------------------------------------
sig_sui <- go_sui[!is.na(go_sui$p.adjust) & go_sui$p.adjust < 0.05, ]
sig_pop <- go_pop[!is.na(go_pop$p.adjust) & go_pop$p.adjust < 0.05, ]
common_terms <- intersect(sig_sui$TERM, sig_pop$TERM)
cat("\nTermos GO BP significativos (p.adjust<0.05) em comum SUI x POP:", length(common_terms), "\n")
if (length(common_terms) > 0) print(common_terms)

cat("\nTop 15 termos GO BP - SUI 36h (por p-valor bruto):\n")
print(head(go_sui[, c("TERM", "Count", "pvalue", "p.adjust")], 15))
cat("\nTop 15 termos GO BP - POP (por p-valor bruto):\n")
print(head(go_pop[, c("TERM", "Count", "pvalue", "p.adjust")], 15))

## -----------------------------------------------------------------------
## PROXIMO PASSO (fora deste script): hub genes + KEGG
## -----------------------------------------------------------------------
## 1. Abra https://string-db.org -> "Multiple proteins" -> cole o conteudo
##    de results/STRING_input_SUI_36h_genes.txt (organismo: Homo sapiens)
##    -> Search. Repita para STRING_input_POP_genes.txt.
## 2. Na pagina de resultado, va em "Exports" -> baixe o TSV da rede
##    (network as text, simple tabular text output) para cada lista.
## 3. Na aba "Analysis" da mesma pagina, o STRING ja mostra enriquecimento
##    KEGG/GO/Reactome para a lista - pode exportar essa tabela tambem.
## 4. Com o TSV da rede em maos, use o pacote igraph para calcular hub
##    genes (grau, closeness, betweenness) - script a ser adicionado apos
##    os arquivos serem fornecidos.
