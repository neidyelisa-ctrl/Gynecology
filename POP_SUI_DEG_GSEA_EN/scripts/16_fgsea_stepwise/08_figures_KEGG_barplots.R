# =============================================================================
# STEP 8 of 13 - KEGG pathway barplots (POP + SUI) and shared-pathway heatmap
# =============================================================================
# RUN ORDER: needs 01, 04, 05 and 06 already sourced in this same R session
# (uses: fgsea_pop from step 4, fgsea_sui from step 5, kegg_both_significant
# from step 6, kegg_label from step 1).
#
# WHAT THIS STEP DOES: four related KEGG barplots/heatmap, each answering a
# slightly different question - see the comment above each block below for
# exactly which pathways can and cannot appear in that specific figure.
#
# OUTPUT:
#   figures/fgsea_KEGG_barplot_POP.png             - top 15 in POP by p-value
#   figures/fgsea_KEGG_barplot_SUI.png              - ALL significant in SUI
#   figures/fgsea_KEGG_barplot_POP_FDR025.png       - POP, strong+moderate bands
#   figures/fgsea_KEGG_significant_both_barplot.png - ALL significant in BOTH
#   figures/heatmap_shared_KEGG_NES_fgsea.png       - NES heatmap, significant
#                                                      in BOTH
# =============================================================================

setwd("E:/POP+SUI ONE")  # safe to repeat - keeps this correct even in a fresh R session

cat("================ STEP 8: KEGG pathway barplots ================\n\n")

## --- quick-glance top-15-by-p-value barplot (POP) ---------------------------
make_fgsea_barplot <- function(fgsea_res, title, n_top = 15) {
  fgsea_res <- fgsea_res[!is.na(fgsea_res$NES), ]
  d <- head(fgsea_res[order(fgsea_res$pval), ], n_top)
  d$Sig <- ifelse(d$padj < 0.05, "FDR<0.05", ifelse(d$padj < 0.25, "FDR<0.25", "NS"))
  d$Label <- factor(kegg_label(d$pathway), levels = rev(kegg_label(d$pathway)))
  ggplot(d, aes(x = NES, y = Label, fill = Sig)) +
    geom_col() +
    scale_fill_manual(values = c("FDR<0.05" = "#B2182B", "FDR<0.25" = "#F4A582", "NS" = "grey70")) +
    geom_vline(xintercept = 0, color = "grey30") +
    labs(title = title, x = "Normalized Enrichment Score (NES)", y = NULL, fill = "Significance") +
    theme_bw() + theme(axis.text.y = element_text(size = 8), plot.title = element_text(size = 12))
}
ggsave("figures/fgsea_KEGG_barplot_POP.png",
       make_fgsea_barplot(fgsea_pop, "GSEA via fgsea - top 15 KEGG pathways in POP (by p-value)"),
       width = 12, height = 6.5, dpi = 300)
cat("Saved: figures/fgsea_KEGG_barplot_POP.png (top 15 by p-value only - POP has\n")
cat(sum(fgsea_pop$padj < 0.25, na.rm = TRUE), "pathways at FDR<0.25 total, too many for one\n")
cat("readable chart; see fgsea_KEGG_barplot_POP_FDR025.png below for the complete picture)\n")

# SUI has far fewer significant pathways than POP, so unlike POP's chart
## above, filtering to padj<0.25 FIRST and only then taking the top N
## comfortably fits all of them - no pathway is cut off by an arbitrary
## top-15-by-p-value-across-everything-tested limit.
n_sig_sui <- sum(fgsea_sui$padj < 0.25, na.rm = TRUE)
ggsave("figures/fgsea_KEGG_barplot_SUI.png",
       make_fgsea_barplot(subset(fgsea_sui, padj < 0.25),
                          paste0("GSEA via fgsea - all ", n_sig_sui, " KEGG pathways significant in SUI (FDR<0.25)"),
                          n_top = n_sig_sui),
       width = 12, height = max(6.5, 0.4 * n_sig_sui + 1.5), dpi = 300)
cat("Saved: figures/fgsea_KEGG_barplot_SUI.png (all", n_sig_sui, "significant pathways shown)\n\n")

## --- KEGG barplot for POP, strong (FDR<0.05) and moderate (FDR 0.05-0.25)
## bands ranked and shown separately, so a moderate-but-real hit like Focal
## adhesion is not crowded out by the strongest hits ------------------------
make_fgsea_band_barplot <- function(fgsea_res, title, n_each = 15) {
  fgsea_res <- fgsea_res[!is.na(fgsea_res$NES), ]
  strong <- fgsea_res[fgsea_res$padj < 0.05, ]
  strong <- strong[order(strong$pval), ][seq_len(min(n_each, nrow(strong))), ]
  moderate <- fgsea_res[fgsea_res$padj >= 0.05 & fgsea_res$padj < 0.25, ]
  moderate <- moderate[order(moderate$pval), ][seq_len(min(n_each, nrow(moderate))), ]
  d <- rbind(strong, moderate)
  d$Band <- ifelse(d$padj < 0.05, "FDR<0.05 (strong)", "FDR 0.05-0.25 (moderate)")
  d$Band <- factor(d$Band, levels = c("FDR<0.05 (strong)", "FDR 0.05-0.25 (moderate)"))
  d$Label <- factor(kegg_label(d$pathway), levels = rev(unique(kegg_label(d$pathway))))
  ggplot(d, aes(x = NES, y = Label, fill = Band)) +
    geom_col() +
    scale_fill_manual(values = c("FDR<0.05 (strong)" = "#B2182B", "FDR 0.05-0.25 (moderate)" = "#F4A582")) +
    geom_vline(xintercept = 0, color = "grey30") +
    facet_wrap(~Band, scales = "free_y", ncol = 1) +
    labs(title = title,
         subtitle = "Top 15 pathways in each band, ranked separately, so moderate hits are not\ncrowded out by the strongest ones",
         x = "Normalized Enrichment Score (NES)", y = NULL) +
    theme_bw() + theme(axis.text.y = element_text(size = 8), legend.position = "none",
                        strip.text = element_text(size = 10, face = "bold"),
                        plot.subtitle = element_text(size = 9))
}
n_fdr025_pop <- sum(fgsea_pop$padj < 0.25, na.rm = TRUE)
p_fdr025_pop <- make_fgsea_band_barplot(fgsea_pop, paste0("GSEA via fgsea - KEGG pathways in POP (", n_fdr025_pop, " total at FDR<0.25)"))
ggsave("figures/fgsea_KEGG_barplot_POP_FDR025.png", p_fdr025_pop, width = 12, height = 11, dpi = 300)
cat("Saved: figures/fgsea_KEGG_barplot_POP_FDR025.png (", n_fdr025_pop,
    "pathways at FDR<0.25 in POP total, top 15 of each band shown )\n\n")

## --- dedicated barplot of EVERY pathway significant in BOTH diseases -------
## Built directly from kegg_both_significant (step 6) - NOT a top-N-by-
## p-value cut of anything, so unlike the two charts above, a pathway here
## can never be crowded out by other, more extreme pathways elsewhere in the
## genome. The two barplots above rank ALL of POP's (or SUI's) significant
## pathways by p-value and show the top 15 per band - correct for "what
## stands out most in POP/SUI on its own", but a pathway that matters for
## the CROSS-DISEASE comparison specifically (i.e. it made the strict
## double-FDR<0.25 bar) can still rank outside the top 15 within its own
## band if enough OTHER, unrelated pathways (e.g. broad housekeeping/
## neurodegenerative pathways sharing ribosomal genes) happen to have even
## smaller p-values. This chart exists specifically so that never happens:
## every pathway significant in both diseases is guaranteed a bar here,
## always.
if (nrow(kegg_both_significant) > 0) {
  both_plot_data <- kegg_both_significant[, c("PATH", "PathwayName", "NES_POP", "padj_POP", "NES_SUI", "padj_SUI")]
  both_long <- rbind(
    data.frame(PathwayName = both_plot_data$PathwayName, NES = both_plot_data$NES_POP,
               padj = both_plot_data$padj_POP, Disease = "POP"),
    data.frame(PathwayName = both_plot_data$PathwayName, NES = both_plot_data$NES_SUI,
               padj = both_plot_data$padj_SUI, Disease = "SUI")
  )
  both_long$FDRLabel <- paste0("FDR=", signif(both_long$padj, 2))
  p_both <- ggplot(both_long, aes(x = NES, y = PathwayName, fill = Disease)) +
    geom_col(position = position_dodge(width = 0.75), width = 0.65) +
    geom_text(aes(label = FDRLabel, hjust = ifelse(NES < 0, 1.05, -0.05)),
              position = position_dodge(width = 0.75), size = 3, color = "grey20") +
    geom_vline(xintercept = 0, color = "grey30") +
    scale_fill_manual(values = c(POP = "#B2182B", SUI = "#2166AC")) +
    scale_x_continuous(expand = expansion(mult = 0.3)) +
    labs(title = paste0("KEGG pathways significant in BOTH POP and SUI (n=", nrow(kegg_both_significant), ")"),
         subtitle = "Every pathway that passed FDR<0.25 in both diseases - none cut off by a top-N limit",
         x = "Normalized Enrichment Score (NES)", y = NULL, fill = NULL) +
    theme_bw() + theme(legend.position = "top", plot.subtitle = element_text(size = 10))
  ggsave("figures/fgsea_KEGG_significant_both_barplot.png", p_both,
         width = 11, height = max(4, 0.6 * nrow(kegg_both_significant) + 1.5), dpi = 300)
  cat("Saved: figures/fgsea_KEGG_significant_both_barplot.png (all", nrow(kegg_both_significant),
      "pathways significant in both diseases, guaranteed complete)\n\n")
} else {
  cat("SKIPPED fgsea_KEGG_significant_both_barplot.png (no pathway significant in both diseases)\n\n")
}

## --- shared-pathway NES heatmap (pathways significant in BOTH) ------------
if (nrow(kegg_both_significant) > 0) {
  heatmap_mat <- as.matrix(kegg_both_significant[, c("NES_POP", "NES_SUI")])
  rownames(heatmap_mat) <- kegg_both_significant$PathwayName
  colnames(heatmap_mat) <- c("POP", "SUI")
  png("figures/heatmap_shared_KEGG_NES_fgsea.png", width = 2600, height = 2600, res = 300)
  pheatmap(heatmap_mat, cluster_cols = FALSE,
           main = "KEGG pathways significant in BOTH diseases - NES in POP vs SUI",
           fontsize_row = 9, color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100))
  dev.off()
  cat("Saved: figures/heatmap_shared_KEGG_NES_fgsea.png\n\n")
} else {
  cat("SKIPPED heatmap_shared_KEGG_NES_fgsea.png (no pathway significant in both diseases)\n\n")
}

cat("=== STEP 8 DONE === Now Source 09_GO_biological_process.R\n")
