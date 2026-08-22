# =============================================================================
# Enriquecimento GO (offline) para o cruzamento SUI(36h) x POP SEM filtro de
# baixa contagem (ver R/analysis_SUI_POP_nofilter.R) + exportacao das listas
# de genes para o STRING (necessario para hub genes / KEGG, dependem de
# string-db.org - nao acessivel deste sandbox).
#
# Reutiliza a mesma funcao de enriquecimento GO (hipergeometrico + BH,
# offline) de R/analysis_GO_enrichment_and_STRING_export.R.
# =============================================================================

suppressMessages({
  library(org.Hs.eg.db)
  library(GO.db)
  library(AnnotationDbi)
})

dir.create("results", showWarnings = FALSE)

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
## SUI 36h (sem filtro) - universo = todos os genes com ortologo humano
## -----------------------------------------------------------------------
sui_full <- read.csv("results/SUI_36hr_nofilter_DESeq2_completo.csv")
load("data/homologeneData2.rda")
homologene_lookup <- function(rat_symbols) {
  rat_rows <- homologeneData2[homologeneData2$Taxonomy == 10116 &
                                 homologeneData2$Gene.Symbol %in% rat_symbols, c("HID", "Gene.Symbol")]
  colnames(rat_rows) <- c("HID", "Rat_Gene_Symbol")
  human_rows <- homologeneData2[homologeneData2$Taxonomy == 9606, c("HID", "Gene.Symbol", "Gene.ID")]
  colnames(human_rows) <- c("HID", "Human_Ortholog_Symbol", "Human_Entrez_ID")
  merge(rat_rows, human_rows, by = "HID")[, c("Rat_Gene_Symbol", "Human_Ortholog_Symbol", "Human_Entrez_ID")]
}
orth_all <- homologene_lookup(unique(sui_full$Rat_Gene_Symbol))
sui_full_orth <- merge(sui_full, orth_all, by = "Rat_Gene_Symbol", all.x = TRUE)

sui_universe <- unique(na.omit(sui_full_orth$Human_Ortholog_Symbol))
sui_hits <- unique(na.omit(sui_full_orth$Human_Ortholog_Symbol[
  !is.na(sui_full_orth$padj) & sui_full_orth$padj < 0.05 & abs(sui_full_orth$log2FoldChange) > 0.5]))

cat("SUI 36h (sem filtro) - universo (com ortologo):", length(sui_universe), "| DEGs:", length(sui_hits), "\n")
go_sui <- run_go_enrichment(sui_hits, sui_universe, ont = "BP")
write.csv(go_sui, "results/GO_BP_SUI_36h_nofilter.csv", row.names = FALSE)
writeLines(sui_hits, "results/STRING_input_SUI_36h_nofilter_genes.txt")

cat("\nTop 15 termos GO BP - SUI 36h sem filtro (por p-valor bruto):\n")
print(head(go_sui[, c("TERM", "Count", "pvalue", "p.adjust")], 15))

## -----------------------------------------------------------------------
## POP (sem filtro)
## -----------------------------------------------------------------------
pop_full <- read.csv("results/POP_GSE208261_nofilter_DESeq2_completo.csv")
pop_full$Human_Entrez_ID <- as.character(pop_full$Human_Entrez_ID)
pop_map_all <- select(org.Hs.eg.db, keys = pop_full$Human_Entrez_ID, keytype = "ENTREZID", columns = "SYMBOL")
pop_full_sym <- merge(pop_full, pop_map_all, by.x = "Human_Entrez_ID", by.y = "ENTREZID", all.x = TRUE)

pop_universe <- unique(na.omit(pop_full_sym$SYMBOL))
pop_hits <- unique(na.omit(pop_full_sym$SYMBOL[
  !is.na(pop_full_sym$padj) & pop_full_sym$padj < 0.05 & abs(pop_full_sym$log2FoldChange) > 0.5]))

cat("\nPOP (sem filtro) - universo:", length(pop_universe), "| DEGs:", length(pop_hits), "\n")
go_pop <- run_go_enrichment(pop_hits, pop_universe, ont = "BP")
write.csv(go_pop, "results/GO_BP_POP_nofilter.csv", row.names = FALSE)
writeLines(pop_hits, "results/STRING_input_POP_nofilter_genes.txt")

cat("\nTop 15 termos GO BP - POP sem filtro (por p-valor bruto):\n")
print(head(go_pop[, c("TERM", "Count", "pvalue", "p.adjust")], 15))

sig_sui <- go_sui[!is.na(go_sui$p.adjust) & go_sui$p.adjust < 0.05, ]
sig_pop <- go_pop[!is.na(go_pop$p.adjust) & go_pop$p.adjust < 0.05, ]
common_terms <- intersect(sig_sui$TERM, sig_pop$TERM)
cat("\nTermos GO BP significativos (p.adjust<0.05) em comum SUI x POP (sem filtro):", length(common_terms), "\n")
if (length(common_terms) > 0) print(common_terms)
