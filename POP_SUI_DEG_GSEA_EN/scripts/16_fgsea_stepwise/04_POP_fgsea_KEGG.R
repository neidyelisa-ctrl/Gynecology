# =============================================================================
# STEP 4 of 13 - KEGG pathway enrichment (fgsea) in POP
# =============================================================================
# RUN ORDER: needs 01-03 already sourced in this same R session (uses:
# fit_pop, kegg_label from steps 1 and 3).
#
# WHAT THIS STEP DOES: runs fgsea - the ONLY GSEA method used anywhere in
# this pipeline - on POP's full ranked gene list (all tested genes, ranked
# by moderated t-statistic, not just the 163 DEG). This is what lets a
# pathway show up as a real signal even when no single one of its genes
# individually reached significance.
#
# OUTPUT: results/POP_fgsea_KEGG_full.csv - every KEGG pathway tested.
# =============================================================================

setwd("E:/POP+SUI FGSEA")  # safe to repeat - keeps this correct even in a fresh R session

cat("================ STEP 4: GSEA (fgsea) in POP ================\n\n")
cat("What GSEA asks: even with no individual gene reaching significance, is a\n")
cat("WHOLE GROUP of genes from the same biological pathway shifted in the same\n")
cat("direction as a block, across the FULL ranked gene list (not just the\n")
cat("163 DEG)? This detects weak, distributed signal a gene-by-gene test\n")
cat("alone would miss.\n\n")

ranked_full <- fit_pop$t[, "group_fullPOP"]
ranked_genes_pop <- names(ranked_full)[order(-ranked_full)]

# For each ranked gene, find out which KEGG pathways it belongs to.
ann_pop <- suppressWarnings(select(org.Hs.eg.db, keys = ranked_genes_pop, keytype = "SYMBOL", columns = "PATH"))
ann_pop <- ann_pop[!is.na(ann_pop$PATH), ]
gs_sizes_pop <- table(ann_pop$PATH)
valid_paths_pop <- names(gs_sizes_pop)[gs_sizes_pop >= 5 & gs_sizes_pop <= 200]  # neither too small nor too large
gene_sets_pop <- split(ann_pop$SYMBOL[ann_pop$PATH %in% valid_paths_pop], ann_pop$PATH[ann_pop$PATH %in% valid_paths_pop])
cat("KEGG pathways tested (5-200 genes):", length(gene_sets_pop), "\n\n")

cat("Running fgsea on POP...\n")
set.seed(208261)
fgsea_pop <- as.data.frame(fgsea(pathways = gene_sets_pop, stats = ranked_full, minSize = 5, maxSize = 200))
fgsea_pop$leadingEdge <- sapply(fgsea_pop$leadingEdge, paste, collapse = "/")
fgsea_pop$PathwayName <- kegg_label(fgsea_pop$pathway)
fgsea_pop <- fgsea_pop[order(fgsea_pop$pval), ]
write.csv(fgsea_pop, "results/POP_fgsea_KEGG_full.csv", row.names = FALSE)

cat("=== RESULT: KEGG pathway enrichment in POP ===\n")
cat("Pathways tested:", nrow(fgsea_pop), "\n")
cat("Significant at FDR<0.25:", sum(fgsea_pop$padj < 0.25, na.rm = TRUE), "\n")
cat("Significant at FDR<0.05:", sum(fgsea_pop$padj < 0.05, na.rm = TRUE), "\n\n")

cat("=== STEP 4 DONE === Now Source 05_SUI_data_DEG_fgsea_KEGG.R\n")
