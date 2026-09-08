# =============================================================================
# STEP 6 of 13 - pathways common to POP and SUI
# =============================================================================
# RUN ORDER: needs 01, 04 and 05 already sourced in this same R session
# (uses: fgsea_pop from step 4, fgsea_sui from step 5, kegg_label from step 1).
#
# WHAT THIS STEP DOES: merges POP's and SUI's fgsea results by KEGG pathway
# ID, flags direction (same/opposite NES sign) and significance in each
# disease, and produces the "significant in both diseases" table that steps
# 8, 12 and 13 (Pathview) all read from directly - the automatic,
# non-hardcoded pathway selection for every downstream step comes from here.
#
# OUTPUT:
#   results/KEGG_all_common_pathways_fgsea.csv        - every pathway tested
#                                                        in both diseases
#   results/KEGG_significant_in_BOTH_diseases_fgsea.csv - only padj<0.25 in
#                                                          BOTH
# =============================================================================

setwd("E:/POP+SUI ONE")  # safe to repeat - keeps this correct even in a fresh R session

cat("================ STEP 6: pathways common to POP and SUI ================\n\n")

kegg_common <- merge(
  data.frame(PATH = fgsea_pop$pathway, PathwayName = kegg_label(fgsea_pop$pathway),
             NES_POP = fgsea_pop$NES, padj_POP = fgsea_pop$padj),
  data.frame(PATH = fgsea_sui$pathway, NES_SUI = fgsea_sui$NES, padj_SUI = fgsea_sui$padj),
  by = "PATH"
)

# The full table: every pathway tested in BOTH diseases (not filtered to
# significant only), flagged with direction and significance columns, so
# borderline/near-significant pathways stay visible - relevant given SUI's
# small sample size (n=3 vs 3), where a pathway just above the FDR<0.25 line
# today could easily fall below it with a few more samples. Sorted by
# whichever disease's own padj is smaller, so the strongest and the
# closest-to-significant pathways float to the top together.
kegg_common_full <- kegg_common
kegg_common_full$Same_direction <- sign(kegg_common_full$NES_POP) == sign(kegg_common_full$NES_SUI)
kegg_common_full$Significant_POP_FDR025 <- kegg_common_full$padj_POP < 0.25
kegg_common_full$Significant_SUI_FDR025 <- kegg_common_full$padj_SUI < 0.25
kegg_common_full$Significant_both_FDR025 <- kegg_common_full$Significant_POP_FDR025 & kegg_common_full$Significant_SUI_FDR025
kegg_common_full <- kegg_common_full[order(pmin(kegg_common_full$padj_POP, kegg_common_full$padj_SUI)), ]
write.csv(kegg_common_full, "results/KEGG_all_common_pathways_fgsea.csv", row.names = FALSE)
cat("Saved: results/KEGG_all_common_pathways_fgsea.csv (", nrow(kegg_common_full), "pathways total )\n")
cat(" ", sum(kegg_common_full$Same_direction), "same-direction,",
    sum(!kegg_common_full$Same_direction), "opposite-direction;\n")
cat(" ", sum(kegg_common_full$Significant_POP_FDR025), "significant in POP,",
    sum(kegg_common_full$Significant_SUI_FDR025), "significant in SUI,",
    sum(kegg_common_full$Significant_both_FDR025), "significant in BOTH )\n\n")

# The small table: only the pathways significant in BOTH diseases (the
# strict double bar - padj<0.25 in POP AND in SUI).
kegg_both_significant <- subset(kegg_common_full, Significant_both_FDR025)
kegg_both_significant <- kegg_both_significant[order(pmin(kegg_both_significant$padj_POP, kegg_both_significant$padj_SUI)), ]
write.csv(kegg_both_significant, "results/KEGG_significant_in_BOTH_diseases_fgsea.csv", row.names = FALSE)
cat("Saved: results/KEGG_significant_in_BOTH_diseases_fgsea.csv (", nrow(kegg_both_significant),
    "pathways significant in both diseases, out of", nrow(kegg_common_full), "tested in both -",
    sum(kegg_both_significant$Same_direction), "of", nrow(kegg_both_significant), "are direction-concordant )\n\n")

cat("=== STEP 6 DONE === Now Source 07_figures_volcano_heatmaps.R\n")
