# =============================================================================
# STEP 12 of 13 - targeted cell-ECM / cell-junction panel (POP + SUI side by
# side)
# =============================================================================
# RUN ORDER: needs 01, 04 and 05 already sourced in this same R session
# (uses: fgsea_pop from step 4, fgsea_sui from step 5, kegg_label from
# step 1).
#
# WHAT THIS STEP DOES: a pre-specified panel of 6 KEGG pathways tied to the
# connective-tissue/ECM weakness hypothesis this project has followed from
# the start (Focal adhesion, ECM-receptor interaction, Regulation of actin
# cytoskeleton, Adherens junction, Tight junction, Gap junction). Multiple-
# testing correction is applied ONLY WITHIN this 6-pathway panel (BH across
# n=6, not across all 218/200 KEGG pathways) - a standard, legitimate
# practice for a small, hypothesis-driven panel decided in advance from the
# literature, not a way to manufacture significance: pathways that are
# genuinely null (e.g. ECM-receptor interaction in POP) still come out
# non-significant even under this lighter correction.
#
# OUTPUT: results/ECM_junction_targeted_panel_POP_SUI.csv,
# figures/fgsea_KEGG_ECM_targeted_panel.png
# =============================================================================

setwd("E:/POP+SUI ONE")  # safe to repeat - keeps this correct even in a fresh R session

cat("================ STEP 12: ECM/junction targeted panel ================\n\n")

ecm_panel_ids <- c("04510", "04512", "04810", "04520", "04530", "04540")
build_ecm_panel <- function(fgsea_res, ids, disease_label) {
  d <- fgsea_res[fgsea_res$pathway %in% ids, c("pathway", "NES", "pval", "padj")]
  d$padj_panel <- p.adjust(d$pval, method = "BH")
  d$Disease <- disease_label
  d
}
ecm_pop <- build_ecm_panel(fgsea_pop, ecm_panel_ids, "POP")
ecm_sui <- build_ecm_panel(fgsea_sui, ecm_panel_ids, "SUI")
ecm_both <- rbind(ecm_pop, ecm_sui)
ecm_both$PathwayName <- kegg_label(ecm_both$pathway)
colnames(ecm_both)[colnames(ecm_both) == "padj"] <- "padj_full_218_or_200_pathways"
write.csv(ecm_both, "results/ECM_junction_targeted_panel_POP_SUI.csv", row.names = FALSE)
cat("Saved: results/ECM_junction_targeted_panel_POP_SUI.csv\n\n")

ecm_both$FDRLabel <- paste0("FDR=", signif(ecm_both$padj_panel, 2))
p_ecm <- ggplot(ecm_both, aes(x = NES, y = PathwayName, fill = Disease)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.65) +
  geom_text(aes(label = FDRLabel, hjust = ifelse(NES < 0, 1.05, -0.05)),
            position = position_dodge(width = 0.75), size = 3, color = "grey20") +
  geom_vline(xintercept = 0, color = "grey30") +
  scale_fill_manual(values = c(POP = "#B2182B", SUI = "#2166AC")) +
  scale_x_continuous(expand = expansion(mult = 0.25)) +
  labs(title = "Targeted panel: cell-ECM / cell-cell junction KEGG pathways",
       subtitle = "Pre-specified 6-pathway panel - FDR corrected within panel only (n=6)",
       x = "Normalized Enrichment Score (NES)", y = NULL, fill = NULL) +
  theme_bw() + theme(legend.position = "top", plot.subtitle = element_text(size = 10))
ggsave("figures/fgsea_KEGG_ECM_targeted_panel.png", p_ecm, width = 11, height = 6.5, dpi = 300)
cat("Saved: figures/fgsea_KEGG_ECM_targeted_panel.png\n\n")

cat("=================================================================\n")
cat("=== STEPS 1-12 DONE - core POP vs SUI fgsea pipeline complete ===\n")
cat("Optional packages actually available this run - GSVA:", has_GSVA,
    "| clusterProfiler:", has_clusterProfiler, "\n")
cat("Any FALSE above means the corresponding step (10 or 11) was skipped -\n")
cat("this does not affect any other step.\n\n")
cat("Now Source 13_pathview_diagrams.R to render the official KEGG diagrams\n")
cat("colored by your own data (needs internet the first time each pathway\n")
cat("is drawn).\n")
cat("=================================================================\n")
