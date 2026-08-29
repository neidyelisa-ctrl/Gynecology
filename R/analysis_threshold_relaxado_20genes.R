# =============================================================================
# Resposta ao feedback do professor (pontos 3 e 4): afrouxar o limiar dentro
# de faixa defensavel para conseguir mais genes em comum, e rodar GO/KEGG
# nessa lista maior.
#
# Testado: FDR<0.05 (padrao) = 5 genes em comum; FDR<0.10 = 20 genes; FDR<0.20
# = 76 genes (mas essa ja fica pouco defensavel - FDR 20% e alto demais para
# alegar "significativo", vira so "sugestivo"). Escolhido FDR<0.10 (ainda
# controle de falsa descoberta, so um pouco mais permissivo - pratica comum
# em genomica exploratoria com N pequeno) mantendo |log2FC|>0.5.
#
# Usa a mesma receita da usuaria (sem filtro de baixa contagem, ortologos
# via babelgene/Ensembl 1-para-1) - POP aqui usa o modelo COM covariavel de
# idade (unico "sem filtro" ja calculado no repositorio; a versao "receita
# exata" sem covariavel so tem os DEGs significativos salvos, nao a tabela
# completa necessaria para o novo limiar - anotado como limitacao).
# =============================================================================

suppressMessages({
  library(org.Hs.eg.db)
  library(GO.db)
})

FDR_RELAXADO <- 0.10
LFC_CUT <- 0.5

load("data/babelgene_orthologs.rda")
rat_orth_all <- unique(orthologs_df[orthologs_df$taxon_id == 10116, c("symbol","human_symbol","human_entrez")])
colnames(rat_orth_all) <- c("Rat_Gene_Symbol","Human_Ortholog_Symbol","Human_Entrez_ID")
rat_n <- table(rat_orth_all$Rat_Gene_Symbol); hum_n <- table(rat_orth_all$Human_Ortholog_Symbol)
rat_orth <- rat_orth_all[rat_n[rat_orth_all$Rat_Gene_Symbol] == 1 & hum_n[rat_orth_all$Human_Ortholog_Symbol] == 1, ]

sui36 <- read.csv("results/SUI_36hr_nofilter_DESeq2_completo.csv")
pop_full <- read.csv("results/POP_GSE208261_nofilter_DESeq2_completo.csv")

sui_orth <- merge(sui36, rat_orth, by = "Rat_Gene_Symbol", all.x = TRUE)
sui_sig <- subset(sui_orth, !is.na(padj) & padj < FDR_RELAXADO & abs(log2FoldChange) > LFC_CUT & !is.na(Human_Entrez_ID))
pop_sig <- subset(pop_full, !is.na(padj) & padj < FDR_RELAXADO & abs(log2FoldChange) > LFC_CUT)
pop_sig$Human_Entrez_ID <- as.character(pop_sig$Human_Entrez_ID)
sui_sig$Human_Entrez_ID <- as.character(sui_sig$Human_Entrez_ID)

common_ids <- intersect(sui_sig$Human_Entrez_ID, pop_sig$Human_Entrez_ID)
sui_match <- sui_sig[sui_sig$Human_Entrez_ID %in% common_ids,
                      c("Rat_Gene_Symbol","Human_Ortholog_Symbol","Human_Entrez_ID","log2FoldChange","padj")]
colnames(sui_match)[4:5] <- c("SUI_36h_log2FC","SUI_36h_padj")
pop_match <- pop_sig[pop_sig$Human_Entrez_ID %in% common_ids, c("Human_Entrez_ID","log2FoldChange","padj")]
colnames(pop_match)[2:3] <- c("POP_log2FC","POP_padj")
common <- merge(sui_match, pop_match, by = "Human_Entrez_ID")
common <- common[order(common$POP_padj), ]
write.csv(common, "results/SUI_36h_x_POP_FDR010_20genes.csv", row.names = FALSE)

cat("Genes em comum com FDR<0.10:", nrow(common), "\n")
print(common[, c("Rat_Gene_Symbol","Human_Ortholog_Symbol","SUI_36h_log2FC","POP_log2FC")])

## -----------------------------------------------------------------------
## GO e KEGG offline nessa lista de 20 genes
## -----------------------------------------------------------------------
run_enrichment <- function(hit_genes, universe_genes, keytype_col, ont_filter = NULL, min_gs = 2, max_gs = 2000) {
  hit_genes <- unique(intersect(hit_genes, universe_genes)); universe_genes <- unique(universe_genes)
  if (length(hit_genes) == 0) return(NULL)
  ann <- suppressWarnings(select(org.Hs.eg.db, keys = universe_genes, keytype = "SYMBOL", columns = keytype_col))
  ann <- ann[!is.na(ann[[keytype_col]]), c("SYMBOL", keytype_col)]
  if (!is.null(ont_filter)) {
    ann2 <- suppressWarnings(select(org.Hs.eg.db, keys = universe_genes, keytype = "SYMBOL", columns = c("GO","ONTOLOGY")))
    keep_go <- unique(ann2$GO[ann2$ONTOLOGY == ont_filter & !is.na(ann2$GO)])
    ann <- ann[ann[[keytype_col]] %in% keep_go, ]
  }
  colnames(ann) <- c("SYMBOL", "TERM_ID"); ann <- unique(ann)
  gs_size <- table(ann$TERM_ID); valid <- names(gs_size)[gs_size >= min_gs & gs_size <= max_gs]
  ann <- ann[ann$TERM_ID %in% valid, ]
  N <- length(universe_genes); n <- length(hit_genes)
  res <- lapply(valid, function(term_id) {
    genes_in_term <- unique(ann$SYMBOL[ann$TERM_ID == term_id]); M <- length(genes_in_term)
    hits_in_term <- intersect(genes_in_term, hit_genes); k <- length(hits_in_term)
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

all_symbols <- keys(org.Hs.eg.db, keytype = "SYMBOL")
hit_genes <- unique(common$Human_Ortholog_Symbol)
cat("\n20 genes (simbolos humanos) usados no enriquecimento:", paste(hit_genes, collapse=", "), "\n")

go <- run_enrichment(hit_genes, all_symbols, "GO", ont_filter = "BP")
if (!is.null(go)) {
  terms <- suppressMessages(select(GO.db, keys = go$TERM_ID, keytype = "GOID", columns = "TERM"))
  go <- merge(go, terms, by.x = "TERM_ID", by.y = "GOID"); go <- go[order(go$pvalue), ]
  write.csv(go, "results/GO_BP_20genes_FDR010.csv", row.names = FALSE)
  cat("\nGO BP - termos com padj<0.05:", sum(go$p.adjust < 0.05, na.rm=TRUE), "de", nrow(go), "\n")
  print(head(go[, c("TERM","Count","pvalue","p.adjust","geneID")], 10))
}

kegg <- run_enrichment(hit_genes, all_symbols, "PATH")
if (!is.null(kegg)) {
  write.csv(kegg, "results/KEGG_20genes_FDR010.csv", row.names = FALSE)
  cat("\nKEGG - vias com padj<0.05:", sum(kegg$p.adjust < 0.05, na.rm=TRUE), "de", nrow(kegg), "\n")
  print(head(kegg[, c("TERM_ID","Count","pvalue","p.adjust","geneID")], 10))
}
