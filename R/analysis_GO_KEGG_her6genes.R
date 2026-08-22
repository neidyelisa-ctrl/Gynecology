# =============================================================================
# Verificacao independente: GO e KEGG para a lista de 6 genes em comum da
# usuaria (KRT10, SERPINB2, CALML5, CWH43, INPP4B, DMKN), extraida do
# thesis_proposal_draft.pdf dela (pipeline via Perplexity).
#
# Objetivo: rodar enriquecimento GO/KEGG de forma independente (nao copiando
# o resultado dela) para checar se os mesmos temas (fosfatidilinositol,
# glicerofosfolipideo, cornificacao, sinalizacao de estrogenio) aparecem.
#
# GO: teste hipergeometrico (phyper) + BH, via org.Hs.eg.db + GO.db (mesmo
# metodo usado no resto do repositorio).
# KEGG: org.Hs.eg.db tem uma coluna PATH (mapeamento antigo gene->KEGG
# pathway ID, ainda mantido no pacote) - usada aqui do mesmo jeito, ja que
# nao ha acesso a internet para KEGGREST/enrichKEGG neste sandbox. Os nomes
# dos pathway IDs (numeros KEGG sao estaveis) foram preenchidos manualmente
# a partir do conhecimento do meu treinamento, nao de uma consulta ao vivo.
#
# Testado com DOIS universos (denominador), para mostrar o efeito da escolha:
#   (a) genoma inteiro anotado (o que enrichGO usa por padrao se voce nao
#       passa "universe" - provavelmente o que o pipeline dela fez)
#   (b) universo restrito aos genes realmente testados no meu DESeq2 (POP +
#       ortologos do SUI) - mais conservador, e o que eu uso no resto do
#       repositorio.
# =============================================================================

suppressMessages({
  library(org.Hs.eg.db)
  library(GO.db)
  library(AnnotationDbi)
})

her_genes <- c("CALML5","KRT10","DMKN","CWH43","SERPINB2","INPP4B")

## -----------------------------------------------------------------------
## Universo (a): genoma inteiro anotado
## -----------------------------------------------------------------------
all_symbols <- keys(org.Hs.eg.db, keytype = "SYMBOL")

## -----------------------------------------------------------------------
## Universo (b): genes realmente testados (POP + SUI ortologos)
## -----------------------------------------------------------------------
pop_full <- read.csv("results/POP_GSE208261_DESeq2_completo.csv")
pop_full$Human_Entrez_ID <- as.character(pop_full$Human_Entrez_ID)
pop_map <- select(org.Hs.eg.db, keys = pop_full$Human_Entrez_ID, keytype = "ENTREZID", columns = "SYMBOL")
pop_universe <- unique(na.omit(pop_map$SYMBOL))

load("data/homologeneData2.rda")
sui_full <- read.csv("results/SUI_36hr_DESeq2_completo.csv")
homologene_lookup <- function(rat_symbols) {
  rat_rows <- homologeneData2[homologeneData2$Taxonomy == 10116 &
                                 homologeneData2$Gene.Symbol %in% rat_symbols, c("HID", "Gene.Symbol")]
  colnames(rat_rows) <- c("HID", "Rat_Gene_Symbol")
  human_rows <- homologeneData2[homologeneData2$Taxonomy == 9606, c("HID", "Gene.Symbol")]
  colnames(human_rows) <- c("HID", "Human_Ortholog_Symbol")
  merge(rat_rows, human_rows, by = "HID")[, c("Rat_Gene_Symbol", "Human_Ortholog_Symbol")]
}
orth <- homologene_lookup(unique(sui_full$Rat_Gene_Symbol))
sui_universe <- unique(na.omit(orth$Human_Ortholog_Symbol))

tested_universe <- unique(c(pop_universe, sui_universe))
cat("Universo (a) genoma inteiro:", length(all_symbols), "genes\n")
cat("Universo (b) genes testados (POP + SUI ortologos):", length(tested_universe), "genes\n")
cat("Dos 6 genes dela, quantos estao no universo (b)?", sum(her_genes %in% tested_universe), "de 6\n\n")

## -----------------------------------------------------------------------
## Funcao generica de enriquecimento (GO ou KEGG/PATH) via hipergeometrico
## -----------------------------------------------------------------------
run_enrichment <- function(hit_genes, universe_genes, keytype_col, ont_filter = NULL, min_gs = 2, max_gs = 2000) {
  hit_genes <- unique(intersect(hit_genes, universe_genes))
  universe_genes <- unique(universe_genes)

  cols_needed <- keytype_col
  ann <- suppressWarnings(select(org.Hs.eg.db, keys = universe_genes, keytype = "SYMBOL", columns = cols_needed))
  ann <- ann[!is.na(ann[[keytype_col]]), c("SYMBOL", keytype_col)]
  if (!is.null(ont_filter)) ann <- ann[ann$ONTOLOGY == ont_filter, ]
  colnames(ann) <- c("SYMBOL", "TERM_ID")
  ann <- unique(ann)

  gs_size <- table(ann$TERM_ID)
  valid <- names(gs_size)[gs_size >= min_gs & gs_size <= max_gs]
  ann <- ann[ann$TERM_ID %in% valid, ]

  N <- length(universe_genes)
  n <- length(hit_genes)

  res <- lapply(valid, function(term_id) {
    genes_in_term <- unique(ann$SYMBOL[ann$TERM_ID == term_id])
    M <- length(genes_in_term)
    hits_in_term <- intersect(genes_in_term, hit_genes)
    k <- length(hits_in_term)
    if (k == 0) return(NULL)
    pval <- phyper(k - 1, M, N - M, n, lower.tail = FALSE)
    data.frame(TERM_ID = term_id, Count = k, M = M, N = N, n = n,
               pvalue = pval, geneID = paste(hits_in_term, collapse = "/"), stringsAsFactors = FALSE)
  })
  res <- do.call(rbind, res)
  if (is.null(res) || nrow(res) == 0) return(NULL)
  res$p.adjust <- p.adjust(res$pvalue, method = "BH")
  res[order(res$pvalue), ]
}

## -----------------------------------------------------------------------
## GO Biological Process
## -----------------------------------------------------------------------
for (uni_name in c("genoma_inteiro", "genes_testados")) {
  uni <- if (uni_name == "genoma_inteiro") all_symbols else tested_universe
  cat("\n\n########## GO Biological Process - universo:", uni_name, "##########\n")
  go_res <- run_enrichment(her_genes, uni, "GO", ont_filter = NULL)
  if (!is.null(go_res)) {
    ann_bp <- suppressWarnings(select(org.Hs.eg.db, keys = uni, keytype = "SYMBOL", columns = c("GO","ONTOLOGY")))
    bp_ids <- unique(ann_bp$GO[ann_bp$ONTOLOGY == "BP" & !is.na(ann_bp$GO)])
    go_res <- go_res[go_res$TERM_ID %in% bp_ids, ]
    terms <- suppressMessages(select(GO.db, keys = go_res$TERM_ID, keytype = "GOID", columns = "TERM"))
    go_res <- merge(go_res, terms, by.x = "TERM_ID", by.y = "GOID")
    go_res <- go_res[order(go_res$pvalue), ]
    print(head(go_res[, c("TERM","Count","M","N","pvalue","p.adjust","geneID")], 20))
    write.csv(go_res, paste0("results/GO_BP_her6genes_", uni_name, ".csv"), row.names = FALSE)
  } else {
    cat("Nenhum termo GO encontrado.\n")
  }
}

## -----------------------------------------------------------------------
## KEGG (via PATH do org.Hs.eg.db)
## -----------------------------------------------------------------------
kegg_names <- c(
  "00562" = "Inositol phosphate metabolism",
  "01100" = "Metabolic pathways",
  "04020" = "Calcium signaling pathway",
  "04070" = "Phosphatidylinositol signaling system",
  "04114" = "Oocyte meiosis",
  "04270" = "Vascular smooth muscle contraction",
  "04720" = "Long-term potentiation",
  "04722" = "Neurotrophin signaling pathway",
  "04740" = "Olfactory transduction",
  "04744" = "Phototransduction",
  "04910" = "Insulin signaling pathway",
  "04912" = "GnRH signaling pathway",
  "04916" = "Melanogenesis",
  "04970" = "Salivary secretion",
  "04971" = "Gastric acid secretion",
  "05010" = "Alzheimer disease",
  "05214" = "Glioma",
  "05150" = "Staphylococcus aureus infection",
  "05146" = "Amoebiasis"
)

for (uni_name in c("genoma_inteiro", "genes_testados")) {
  uni <- if (uni_name == "genoma_inteiro") all_symbols else tested_universe
  cat("\n\n########## KEGG (PATH) - universo:", uni_name, "##########\n")
  kegg_res <- run_enrichment(her_genes, uni, "PATH")
  if (!is.null(kegg_res)) {
    kegg_res$PATH_NAME <- kegg_names[kegg_res$TERM_ID]
    kegg_res <- kegg_res[order(kegg_res$pvalue), ]
    print(kegg_res[, c("TERM_ID","PATH_NAME","Count","M","N","pvalue","p.adjust","geneID")])
    write.csv(kegg_res, paste0("results/KEGG_her6genes_", uni_name, ".csv"), row.names = FALSE)
  } else {
    cat("Nenhuma via KEGG encontrada.\n")
  }
}
