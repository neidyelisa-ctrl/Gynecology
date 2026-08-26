# =============================================================================
# Pipeline completo (COM filtro de baixa contagem) - POP e SUI 72h (prioridade),
# depois SUI 36h: DEGs -> ortologos humanos (via babelgene, inclui Ensembl) ->
# GO e KEGG dos DEGs -> hub genes -> GO e KEGG dos hub genes.
#
# Usa os resultados DESeq2 JA CALCULADOS com filtro de baixa contagem
# (rowSums(contagem>=10)>=N amostras), covariavel de idade no POP, shrinkage
# tipo "normal" - o pipeline principal do repositorio (results/*_DESeq2_
# completo.csv), NAO a versao sem filtro pedida em outro momento.
#
# Ortologos: pacote babelgene (offline, MIT, ver R/analysis_SUI_POP_ensembl_
# orthologs.R para o histórico de por que HomoloGene sozinho nao bastava),
# com filtro de ortologia 1-para-1.
# =============================================================================

suppressMessages({
  library(org.Hs.eg.db)
  library(GO.db)
  library(AnnotationDbi)
})

FDR_CUT <- 0.05
LFC_CUT <- 0.5
dir.create("results", showWarnings = FALSE)

## -----------------------------------------------------------------------
## Ortologos (babelgene, 1-para-1)
## -----------------------------------------------------------------------
load("data/babelgene_orthologs.rda")
rat_orth_all <- unique(orthologs_df[orthologs_df$taxon_id == 10116,
                                     c("symbol", "human_symbol", "human_entrez")])
colnames(rat_orth_all) <- c("Rat_Gene_Symbol", "Human_Ortholog_Symbol", "Human_Entrez_ID")
rat_n <- table(rat_orth_all$Rat_Gene_Symbol)
hum_n <- table(rat_orth_all$Human_Ortholog_Symbol)
rat_orth <- rat_orth_all[rat_n[rat_orth_all$Rat_Gene_Symbol] == 1 &
                            hum_n[rat_orth_all$Human_Ortholog_Symbol] == 1, ]

## -----------------------------------------------------------------------
## Funcoes de enriquecimento (reaproveitadas do resto do repositorio)
## -----------------------------------------------------------------------
run_enrichment <- function(hit_genes, universe_genes, keytype_col, ont_filter = NULL, min_gs = 3, max_gs = 2000) {
  hit_genes <- unique(intersect(hit_genes, universe_genes))
  universe_genes <- unique(universe_genes)
  if (length(hit_genes) == 0) return(NULL)

  ann <- suppressWarnings(select(org.Hs.eg.db, keys = universe_genes, keytype = "SYMBOL", columns = keytype_col))
  ann <- ann[!is.na(ann[[keytype_col]]), c("SYMBOL", keytype_col)]
  if (!is.null(ont_filter)) {
    ann2 <- suppressWarnings(select(org.Hs.eg.db, keys = universe_genes, keytype = "SYMBOL", columns = c("GO","ONTOLOGY")))
    keep_go <- unique(ann2$GO[ann2$ONTOLOGY == ont_filter & !is.na(ann2$GO)])
    ann <- ann[ann[[keytype_col]] %in% keep_go, ]
  }
  colnames(ann) <- c("SYMBOL", "TERM_ID")
  ann <- unique(ann)

  gs_size <- table(ann$TERM_ID)
  valid <- names(gs_size)[gs_size >= min_gs & gs_size <= max_gs]
  ann <- ann[ann$TERM_ID %in% valid, ]

  N <- length(universe_genes); n <- length(hit_genes)
  res <- lapply(valid, function(term_id) {
    genes_in_term <- unique(ann$SYMBOL[ann$TERM_ID == term_id])
    M <- length(genes_in_term)
    hits_in_term <- intersect(genes_in_term, hit_genes)
    k <- length(hits_in_term)
    if (k == 0) return(NULL)
    pval <- phyper(k - 1, M, N - M, n, lower.tail = FALSE)
    data.frame(TERM_ID = term_id, Count = k, M = M, N = N, n = n, pvalue = pval,
               geneID = paste(hits_in_term, collapse = "/"), stringsAsFactors = FALSE)
  })
  res <- do.call(rbind, res)
  if (is.null(res) || nrow(res) == 0) return(NULL)
  res$p.adjust <- p.adjust(res$pvalue, method = "BH")
  res[order(res$pvalue), ]
}

kegg_names <- c(
  "00562"="Inositol phosphate metabolism","01100"="Metabolic pathways",
  "04020"="Calcium signaling pathway","04070"="Phosphatidylinositol signaling system",
  "04151"="PI3K-Akt signaling pathway","04350"="TGF-beta signaling pathway",
  "04510"="Focal adhesion","04512"="ECM-receptor interaction",
  "04514"="Cell adhesion molecules","04670"="Leukocyte transendothelial migration",
  "04810"="Regulation of actin cytoskeleton","05205"="Proteoglycans in cancer",
  "04330"="Notch signaling pathway","04390"="Hippo signaling pathway",
  "04115"="p53 signaling pathway","04110"="Cell cycle","04114"="Oocyte meiosis",
  "04150"="mTOR signaling pathway","04910"="Insulin signaling pathway",
  "04150"="mTOR signaling pathway","04022"="cGMP-PKG signaling pathway",
  "04270"="Vascular smooth muscle contraction","04261"="Adrenergic signaling in cardiomyocytes",
  "04217"="Necroptosis","04371"="Apelin signaling pathway",
  "04970"="Salivary secretion","04742"="Taste transduction",
  "05414"="Dilated cardiomyopathy","00230"="Purine metabolism",
  "05410"="Hypertrophic cardiomyopathy","04260"="Cardiac muscle contraction",
  "04540"="Gap junction","04724"="Glutamatergic synapse","04360"="Axon guidance",
  "04024"="cAMP signaling pathway","04972"="Pancreatic secretion",
  "04973"="Carbohydrate digestion and absorption","04974"="Protein digestion and absorption",
  "04923"="Regulation of lipolysis in adipocytes","04928"="Parathyroid hormone synthesis, secretion and action",
  "05412"="Arrhythmogenic right ventricular cardiomyopathy","04725"="Cholinergic synapse",
  "04713"="Circadian entrainment","04730"="Long-term depression"
)

run_go_kegg <- function(hit_genes, universe_genes, label) {
  go <- run_enrichment(hit_genes, universe_genes, "GO", ont_filter = "BP")
  if (!is.null(go)) {
    terms <- suppressMessages(select(GO.db, keys = go$TERM_ID, keytype = "GOID", columns = "TERM"))
    go <- merge(go, terms, by.x = "TERM_ID", by.y = "GOID")
    go <- go[order(go$pvalue), ]
    write.csv(go, paste0("results/GO_BP_", label, ".csv"), row.names = FALSE)
    cat(label, "- GO BP termos com padj<0.05:", sum(go$p.adjust < 0.05, na.rm=TRUE), "de", nrow(go), "testados\n")
  } else cat(label, "- GO BP: nenhum termo\n")

  kegg <- run_enrichment(hit_genes, universe_genes, "PATH")
  if (!is.null(kegg)) {
    kegg$PATH_NAME <- ifelse(kegg$TERM_ID %in% names(kegg_names), kegg_names[kegg$TERM_ID], NA)
    kegg <- kegg[order(kegg$pvalue), ]
    write.csv(kegg, paste0("results/KEGG_", label, ".csv"), row.names = FALSE)
    cat(label, "- KEGG vias com padj<0.05:", sum(kegg$p.adjust < 0.05, na.rm=TRUE), "de", nrow(kegg), "testadas\n")
  } else cat(label, "- KEGG: nenhuma via\n")
  list(go = go, kegg = kegg)
}

## -----------------------------------------------------------------------
## POP - ja em simbolo humano, nao precisa de ortologo
## -----------------------------------------------------------------------
pop_full <- read.csv("results/POP_GSE208261_DESeq2_completo.csv")
pop_full$Human_Entrez_ID <- as.character(pop_full$Human_Entrez_ID)
pop_map <- suppressWarnings(select(org.Hs.eg.db, keys = pop_full$Human_Entrez_ID, keytype = "ENTREZID", columns = "SYMBOL"))
pop_full_sym <- merge(pop_full, pop_map, by.x = "Human_Entrez_ID", by.y = "ENTREZID", all.x = TRUE)
pop_universe <- unique(na.omit(pop_full_sym$SYMBOL))
pop_deg <- unique(na.omit(pop_full_sym$SYMBOL[!is.na(pop_full_sym$padj) &
             pop_full_sym$padj < FDR_CUT & abs(pop_full_sym$log2FoldChange) > LFC_CUT]))
cat("\n=== POP (com filtro) - DEGs:", length(pop_deg), "===\n")
writeLines(pop_deg, "results/STRING_input_POP_filtrado_genes.txt")
res_pop_deg <- run_go_kegg(pop_deg, pop_universe, "POP_DEG_filtrado")

## -----------------------------------------------------------------------
## SUI 72h (PRIORIDADE) - converter para ortologos humanos primeiro
## -----------------------------------------------------------------------
sui72_full <- read.csv("results/SUI_72hr_DESeq2_completo.csv")
sui72_orth <- merge(sui72_full, rat_orth, by = "Rat_Gene_Symbol", all.x = TRUE)
sui72_universe <- unique(na.omit(sui72_orth$Human_Ortholog_Symbol))
sui72_deg <- unique(na.omit(sui72_orth$Human_Ortholog_Symbol[!is.na(sui72_orth$padj) &
               sui72_orth$padj < FDR_CUT & abs(sui72_orth$log2FoldChange) > LFC_CUT]))
cat("\n=== SUI 72h (com filtro) - DEGs (ratos):", sum(!is.na(sui72_full$padj) & sui72_full$padj<FDR_CUT & abs(sui72_full$log2FoldChange)>LFC_CUT),
    "| com ortologo humano 1-para-1:", length(sui72_deg), "===\n")
write.csv(sui72_orth[!is.na(sui72_orth$Human_Ortholog_Symbol), ], "results/SUI_72hr_DESeq2_com_ortologos_babelgene.csv", row.names = FALSE)
writeLines(sui72_deg, "results/STRING_input_SUI_72h_filtrado_genes.txt")
res_sui72_deg <- run_go_kegg(sui72_deg, sui72_universe, "SUI_72h_DEG_filtrado")

## -----------------------------------------------------------------------
## SUI 36h (prioridade menor) - mesma logica
## -----------------------------------------------------------------------
sui36_full <- read.csv("results/SUI_36hr_DESeq2_completo.csv")
sui36_orth <- merge(sui36_full, rat_orth, by = "Rat_Gene_Symbol", all.x = TRUE)
sui36_universe <- unique(na.omit(sui36_orth$Human_Ortholog_Symbol))
sui36_deg <- unique(na.omit(sui36_orth$Human_Ortholog_Symbol[!is.na(sui36_orth$padj) &
               sui36_orth$padj < FDR_CUT & abs(sui36_orth$log2FoldChange) > LFC_CUT]))
cat("\n=== SUI 36h (com filtro) - DEGs (ratos):", sum(!is.na(sui36_full$padj) & sui36_full$padj<FDR_CUT & abs(sui36_full$log2FoldChange)>LFC_CUT),
    "| com ortologo humano 1-para-1:", length(sui36_deg), "===\n")
write.csv(sui36_orth[!is.na(sui36_orth$Human_Ortholog_Symbol), ], "results/SUI_36hr_DESeq2_com_ortologos_babelgene.csv", row.names = FALSE)
writeLines(sui36_deg, "results/STRING_input_SUI_36h_filtrado_genes.txt")
res_sui36_deg <- run_go_kegg(sui36_deg, sui36_universe, "SUI_36h_DEG_filtrado")

cat("\n=== RESUMO DEGs (com filtro, ortologos babelgene) ===\n")
cat("POP:", length(pop_deg), "| SUI 72h:", length(sui72_deg), "| SUI 36h:", length(sui36_deg), "\n")
