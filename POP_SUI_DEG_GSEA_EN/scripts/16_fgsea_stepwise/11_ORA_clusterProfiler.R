# =============================================================================
# STEP 11 of 13 - ORA (Over-Representation Analysis) via clusterProfiler
# (optional - skips itself if clusterProfiler is not installed)
# =============================================================================
# RUN ORDER: needs 01, 03, 04, 05 and 09 already sourced in this same R
# session (uses: has_clusterProfiler/kegg_tab from step 1, pop_deg/pop_full
# from step 3, gene_sets_pop from step 4, sui_full from step 5,
# gene_sets_go from step 9).
#
# WHAT THIS STEP ADDS: fgsea (steps 4, 5, 9) and GSVA (step 10) both use the
# FULL ranked gene list. ORA is a third, simpler kind of test: it takes just
# a FIXED LIST of significant genes (POP's 163 DEG, SUI's Wei 2020 DEG
# list) and asks "is this pathway over-represented in that list, compared to
# what would be expected by chance from the full universe of tested genes".
# This is the same style of test the STRING website's "Analysis" tab runs.
#
# DESIGN CHOICE: clusterProfiler's own convenience functions (enrichKEGG(),
# enrichGO()) contact an online server live, every time they run. Using
# clusterProfiler's generic enricher() instead, with the SAME KEGG and GO
# gene sets already built from org.Hs.eg.db earlier in this pipeline, keeps
# this fully offline.
#
# OUTPUT: results/POP_ORA_KEGG.csv, results/POP_ORA_GO_BP.csv,
# results/SUI_ORA_KEGG.csv, results/SUI_ORA_GO_BP.csv (only if
# clusterProfiler is installed).
# =============================================================================

setwd("E:/POP+SUI FGSEA")  # safe to repeat - keeps this correct even in a fresh R session

cat("================ STEP 11: ORA via clusterProfiler ================\n\n")

if (!has_clusterProfiler) {
  cat("clusterProfiler is not installed/available in this R session -\n")
  cat("SKIPPING step 11 entirely. This does not affect any other step.\n\n")
} else {

gene_sets_to_term2gene <- function(gene_sets) {
  do.call(rbind, lapply(names(gene_sets), function(nm) data.frame(term = nm, gene = gene_sets[[nm]])))
}
kegg_term2gene <- gene_sets_to_term2gene(gene_sets_pop)
kegg_term2name <- data.frame(term = kegg_tab$PATH5, name = kegg_tab$Name)
go_term2gene <- gene_sets_to_term2gene(gene_sets_go)
go_term2name <- data.frame(term = names(gene_sets_go), name = go_term_name(names(gene_sets_go)))

## --- ORA on POP (proper universe available: all 22,253 tested genes) -----
cat("Running ORA on POP's 163 DEG (universe = all", nrow(pop_full), "genes DESeq2 tested)...\n")
ora_kegg_pop <- as.data.frame(clusterProfiler::enricher(
  gene = pop_deg$Gene, universe = pop_full$Gene,
  TERM2GENE = kegg_term2gene, TERM2NAME = kegg_term2name,
  pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1))
write.csv(ora_kegg_pop, "results/POP_ORA_KEGG.csv", row.names = FALSE)
cat("KEGG pathways significant at FDR<0.25:", sum(ora_kegg_pop$p.adjust < 0.25, na.rm = TRUE),
    "of", nrow(ora_kegg_pop), "\n")

ora_go_pop <- as.data.frame(clusterProfiler::enricher(
  gene = pop_deg$Gene, universe = pop_full$Gene,
  TERM2GENE = go_term2gene, TERM2NAME = go_term2name,
  pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1))
write.csv(ora_go_pop, "results/POP_ORA_GO_BP.csv", row.names = FALSE)
cat("GO BP terms significant at FDR<0.25:", sum(ora_go_pop$p.adjust < 0.25, na.rm = TRUE),
    "of", nrow(ora_go_pop), "\n\n")

## --- ORA on SUI (NO proper universe available - see note below) ----------
cat("Running ORA on SUI's DEG list (", nrow(sui_full), "genes)...\n")
cat("CAVEAT: Wei 2020 never published the full array background, only their\n")
cat("own already-significant DEG list. Without a real universe, enricher()\n")
cat("falls back to using every gene annotated in org.Hs.eg.db as the\n")
cat("background, which is NOT the actual array tested - this makes the SUI\n")
cat("ORA p-values anti-conservative (biased toward looking MORE significant\n")
cat("than they really are). Treat results/SUI_ORA_*.csv as qualitative /\n")
cat("hypothesis-generating only - a limitation of the published SUI data,\n")
cat("not of this code.\n\n")
ora_kegg_sui <- as.data.frame(clusterProfiler::enricher(
  gene = sui_full$GeneSymbol,
  TERM2GENE = kegg_term2gene, TERM2NAME = kegg_term2name,
  pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1))
write.csv(ora_kegg_sui, "results/SUI_ORA_KEGG.csv", row.names = FALSE)
cat("KEGG pathways significant at FDR<0.25:", sum(ora_kegg_sui$p.adjust < 0.25, na.rm = TRUE),
    "of", nrow(ora_kegg_sui), "(caveat above applies)\n")

ora_go_sui <- as.data.frame(clusterProfiler::enricher(
  gene = sui_full$GeneSymbol,
  TERM2GENE = go_term2gene, TERM2NAME = go_term2name,
  pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1))
write.csv(ora_go_sui, "results/SUI_ORA_GO_BP.csv", row.names = FALSE)
cat("GO BP terms significant at FDR<0.25:", sum(ora_go_sui$p.adjust < 0.25, na.rm = TRUE),
    "of", nrow(ora_go_sui), "(caveat above applies)\n\n")

} # end if (has_clusterProfiler)

cat("=== STEP 11 DONE === Now Source 12_ECM_junction_panel.R\n")
