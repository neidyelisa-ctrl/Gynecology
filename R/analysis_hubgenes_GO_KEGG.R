# =============================================================================
# GO e KEGG dos HUB GENES (separado dos DEGs) - POP, SUI 72h (prioridade),
# SUI 36h. Reaproveita os hub genes ja calculados em results/hub_genes_*.csv
# (via igraph sobre as redes STRING ja exportadas) e roda o mesmo
# enriquecimento offline (hipergeometrico) usado para os DEGs.
#
# NOTA IMPORTANTE: as redes STRING usadas aqui foram exportadas pela usuaria
# a partir de listas de genes baseadas em HomoloGene, nao das listas mais
# completas baseadas em babelgene/Ensembl calculadas em
# analysis_full_pipeline_filtered_orthologs.R. Os hub genes abaixo sao uma
# APROXIMACAO com os dados que ja temos - cobrem boa parte da lista nova
# (POP: 337 de 624 genes; SUI 72h: os 10 genes que ja estavam na rede antiga)
# mas nao 100%. Para um resultado definitivo, a usuaria precisaria exportar
# do STRING novamente com as listas atualizadas
# (results/STRING_input_*_filtrado_genes.txt).
# =============================================================================

suppressMessages({
  library(org.Hs.eg.db)
  library(GO.db)
})

run_enrichment <- function(hit_genes, universe_genes, keytype_col, ont_filter = NULL, min_gs = 2, max_gs = 2000) {
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
  "04022"="cGMP-PKG signaling pathway",
  "04270"="Vascular smooth muscle contraction","04261"="Adrenergic signaling in cardiomyocytes",
  "04217"="Necroptosis","04371"="Apelin signaling pathway",
  "04970"="Salivary secretion","04742"="Taste transduction",
  "05414"="Dilated cardiomyopathy","00230"="Purine metabolism",
  "05410"="Hypertrophic cardiomyopathy","04260"="Cardiac muscle contraction",
  "04540"="Gap junction","04724"="Glutamatergic synapse","04360"="Axon guidance",
  "04024"="cAMP signaling pathway","04380"="Osteoclast differentiation",
  "04664"="Fc epsilon RI signaling pathway","04640"="Hematopoietic cell lineage",
  "04611"="Platelet activation","04storeIgnore"=""
)
kegg_names <- kegg_names[names(kegg_names) != "04storeIgnore"]

run_go_kegg <- function(hit_genes, universe_genes, label) {
  go <- run_enrichment(hit_genes, universe_genes, "GO", ont_filter = "BP")
  if (!is.null(go)) {
    terms <- suppressMessages(select(GO.db, keys = go$TERM_ID, keytype = "GOID", columns = "TERM"))
    go <- merge(go, terms, by.x = "TERM_ID", by.y = "GOID")
    go <- go[order(go$pvalue), ]
    write.csv(go, paste0("results/GO_BP_", label, ".csv"), row.names = FALSE)
    cat(label, "- GO BP termos com padj<0.05:", sum(go$p.adjust < 0.05, na.rm=TRUE), "de", nrow(go), "testados\n")
    print(head(go[, c("TERM","Count","pvalue","p.adjust","geneID")], 8))
  } else cat(label, "- GO BP: nenhum termo\n")

  kegg <- run_enrichment(hit_genes, universe_genes, "PATH")
  if (!is.null(kegg)) {
    kegg$PATH_NAME <- ifelse(kegg$TERM_ID %in% names(kegg_names), kegg_names[kegg$TERM_ID], NA)
    kegg <- kegg[order(kegg$pvalue), ]
    write.csv(kegg, paste0("results/KEGG_", label, ".csv"), row.names = FALSE)
    cat(label, "- KEGG vias com padj<0.05:", sum(kegg$p.adjust < 0.05, na.rm=TRUE), "de", nrow(kegg), "testadas\n")
    print(head(kegg[, c("TERM_ID","PATH_NAME","Count","pvalue","p.adjust","geneID")], 8))
  } else cat(label, "- KEGG: nenhuma via\n")
  list(go = go, kegg = kegg)
}

## Universo = genoma inteiro anotado (padrao, mesma escolha usada na analise
## dos 6 genes validados) - mais interpretavel que usar so os nos da rede
## como universo, que para redes pequenas (10-23 nos) deixaria o teste quase
## sem poder estatistico nenhum.
all_symbols <- keys(org.Hs.eg.db, keytype = "SYMBOL")

## -----------------------------------------------------------------------
## POP - top 30 hub genes por hub_score
## -----------------------------------------------------------------------
hub_pop <- read.csv("results/hub_genes_POP.csv")
hub_pop <- hub_pop[order(hub_pop$hub_score), ]
pop_hub_top <- head(hub_pop$gene, 30)
cat("\n=== POP - hub genes (top 30 de", nrow(hub_pop), "nos da rede) ===\n")
res_pop_hub <- run_go_kegg(pop_hub_top, all_symbols, "POP_HUB_genes")

## -----------------------------------------------------------------------
## SUI 72h - rede pequena (10 nos), usa todos como "hub genes"
## -----------------------------------------------------------------------
hub_sui72 <- read.csv("results/hub_genes_SUI_72h.csv")
sui72_hub_all <- hub_sui72$gene
cat("\n=== SUI 72h - hub genes (todos os", length(sui72_hub_all), "nos da rede, muito pequena para subdividir) ===\n")
res_sui72_hub <- run_go_kegg(sui72_hub_all, all_symbols, "SUI_72h_HUB_genes")

## -----------------------------------------------------------------------
## SUI 36h - rede pequena (23 nos), usa todos como "hub genes"
## -----------------------------------------------------------------------
hub_sui36 <- read.csv("results/hub_genes_SUI_36h.csv")
sui36_hub_all <- hub_sui36$gene
cat("\n=== SUI 36h - hub genes (todos os", length(sui36_hub_all), "nos da rede) ===\n")
res_sui36_hub <- run_go_kegg(sui36_hub_all, all_symbols, "SUI_36h_HUB_genes")
