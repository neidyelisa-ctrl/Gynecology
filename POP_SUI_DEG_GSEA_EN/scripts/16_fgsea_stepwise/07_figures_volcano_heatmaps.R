# =============================================================================
# STEP 7 of 13 - figures: volcano plots, expression heatmaps, shared-pathway
# NES scatter plot
# =============================================================================
# RUN ORDER: needs 01, 02, 03, 05 and 06 already sourced in this same R
# session (uses: pop_full/pop_deg/voom_fit/group_full from steps 2-3,
# sui_full/sui_mat/group_sui from step 5, kegg_common_full from step 6).
#
# WHAT THIS STEP DOES: builds the gene-level figures (volcano plots for POP
# and SUI, expression heatmaps of the top DEG in each) plus a pathway-level
# scatter plot comparing NES between POP and SUI for every pathway
# significant in at least one disease.
#
# OUTPUT:
#   figures/volcano_POP.png
#   figures/volcano_SUI.png
#   figures/heatmap_POP_topDEG.png
#   figures/heatmap_SUI_topDEG.png
#   figures/shared_pathways_NES_comparison.png
# =============================================================================

setwd("E:/POP+SUI ONE")  # safe to repeat - keeps this correct even in a fresh R session

cat("================ STEP 7: figures ================\n\n")
cat("All figures are saved into the 'figures' subfolder of your working\n")
cat("directory (", getwd(), "\\figures).\n\n", sep = "")

## --- volcano plot - POP -----------------------------------------------------
pop_full$sig <- "NS"
pop_full$sig[pop_full$logFC > 1 & pop_full$adj.P.Val < 0.05] <- "Up"
pop_full$sig[pop_full$logFC < -1 & pop_full$adj.P.Val < 0.05] <- "Down"
pop_full$sig <- factor(pop_full$sig, levels = c("Down", "NS", "Up"))

n_label_each <- 15
top_up_pop <- pop_full[pop_full$sig == "Up", ]
top_up_pop <- top_up_pop[order(top_up_pop$adj.P.Val), ][seq_len(min(n_label_each, nrow(top_up_pop))), ]
top_down_pop <- pop_full[pop_full$sig == "Down", ]
top_down_pop <- top_down_pop[order(top_down_pop$adj.P.Val), ][seq_len(min(n_label_each, nrow(top_down_pop))), ]
labels_pop <- rbind(top_up_pop, top_down_pop)

p_volcano_pop <- ggplot(pop_full, aes(x = logFC, y = -log10(P.Value), color = sig)) +
  geom_point(alpha = 0.6, size = 1.2) +
  scale_color_manual(values = c(Down = "#2166AC", NS = "grey75", Up = "#B2182B")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey40") +
  ggrepel::geom_text_repel(data = labels_pop, aes(label = Gene), size = 2.8,
                            max.overlaps = 100, segment.size = 0.2, show.legend = FALSE) +
  labs(title = paste0("Volcano - POP (GSE208261, 12x12), ", nrow(pop_deg), " DEG at FDR<0.05"),
       subtitle = paste0("Top ", nrow(top_up_pop), " Up and top ", nrow(top_down_pop), " Down genes labeled (by FDR)"),
       x = "log2(Fold Change)", y = "-log10(p-value)", color = NULL) +
  theme_bw() + theme(legend.position = "top")
ggsave("figures/volcano_POP.png", p_volcano_pop, width = 9, height = 7, dpi = 300)
cat("Saved: figures/volcano_POP.png\n\n")

## --- volcano plot - SUI ------------------------------------------------------
sui_full$sig <- "NS"
sui_full$sig[sui_full$logFC > 1 & sui_full$FDR < 0.05] <- "Up"
sui_full$sig[sui_full$logFC < -1 & sui_full$FDR < 0.05] <- "Down"
sui_full$sig <- factor(sui_full$sig, levels = c("Down", "NS", "Up"))

top_up_sui <- sui_full[sui_full$sig == "Up", ]
top_up_sui <- top_up_sui[order(top_up_sui$FDR), ][seq_len(min(n_label_each, nrow(top_up_sui))), ]
top_down_sui <- sui_full[sui_full$sig == "Down", ]
top_down_sui <- top_down_sui[order(top_down_sui$FDR), ][seq_len(min(n_label_each, nrow(top_down_sui))), ]
labels_sui <- rbind(top_up_sui, top_down_sui)

p_volcano_sui <- ggplot(sui_full, aes(x = logFC, y = -log10(PValue), color = sig)) +
  geom_point(alpha = 0.6, size = 1.2) +
  scale_color_manual(values = c(Down = "#2166AC", NS = "grey75", Up = "#B2182B")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey40") +
  ggrepel::geom_text_repel(data = labels_sui, aes(label = GeneSymbol), size = 2.8,
                            max.overlaps = 100, segment.size = 0.2, show.legend = FALSE) +
  labs(title = paste0("Volcano - SUI (Wei 2020, 3x3), ",
                       sum(sui_full$sig != "NS"), " genes at |log2FC|>1 & FDR<0.05"),
       subtitle = paste0("Top ", nrow(top_up_sui), " Up and top ", nrow(top_down_sui), " Down genes labeled (by FDR)"),
       x = "log2(Fold Change)", y = "-log10(p-value)", color = NULL) +
  theme_bw() + theme(legend.position = "top")
ggsave("figures/volcano_SUI.png", p_volcano_sui, width = 9, height = 7, dpi = 300)
cat("Saved: figures/volcano_SUI.png\n\n")

## --- expression heatmap - top DEG, POP --------------------------------------
top_n_heatmap <- min(40, nrow(pop_deg))
top_genes_pop <- pop_deg$Gene[order(pop_deg$adj.P.Val)][seq_len(top_n_heatmap)]
mat_pop <- voom_fit$E[top_genes_pop, , drop = FALSE]
mat_pop_z <- t(scale(t(mat_pop)))
ann_col_pop <- data.frame(Group = group_full, row.names = colnames(mat_pop_z))
png("figures/heatmap_POP_topDEG.png", width = 2400, height = 3000, res = 300)
pheatmap(mat_pop_z, annotation_col = ann_col_pop, show_colnames = FALSE,
         main = paste0("Top ", top_n_heatmap, " DEG (by FDR) - POP vs Control (z-score)"),
         fontsize_row = 7, color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100))
dev.off()
cat("Saved: figures/heatmap_POP_topDEG.png\n\n")

## --- expression heatmap - top DEG, SUI ---------------------------------------
top_genes_sui <- sui_full$GeneSymbol[order(sui_full$FDR)][seq_len(min(40, nrow(sui_full)))]
top_genes_sui <- intersect(top_genes_sui, rownames(sui_mat))
mat_sui <- sui_mat[top_genes_sui, , drop = FALSE]
mat_sui_z <- t(scale(t(mat_sui)))
ann_col_sui <- data.frame(Group = group_sui, row.names = colnames(mat_sui_z))
png("figures/heatmap_SUI_topDEG.png", width = 2400, height = 3000, res = 300)
pheatmap(mat_sui_z, annotation_col = ann_col_sui, show_colnames = TRUE,
         main = paste0("Top ", length(top_genes_sui), " DEG (by FDR) - SUI vs Ctrl (z-score)"),
         fontsize_row = 7, color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100))
dev.off()
cat("Saved: figures/heatmap_SUI_topDEG.png\n\n")

## --- shared-pathway NES scatter (all pathways significant in at least one
## disease, colored by whether they are significant in BOTH) ----------------
scatter_data <- subset(kegg_common_full, Significant_POP_FDR025 | Significant_SUI_FDR025)
p_shared <- ggplot(scatter_data, aes(x = NES_POP, y = NES_SUI, label = PathwayName,
                                       color = Significant_both_FDR025)) +
  geom_hline(yintercept = 0, color = "grey70") + geom_vline(xintercept = 0, color = "grey70") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(size = 3) +
  scale_color_manual(values = c("TRUE" = "#B2182B", "FALSE" = "#2166AC"),
                      labels = c("TRUE" = "Significant in both (FDR<0.25)", "FALSE" = "Significant in one only")) +
  geom_text_repel(data = subset(scatter_data, Significant_both_FDR025), size = 3, max.overlaps = 30) +
  labs(title = "KEGG pathways significant in POP and/or SUI: NES side by side",
       x = "NES - POP", y = "NES - SUI", color = NULL) +
  theme_bw() + theme(legend.position = "top")
ggsave("figures/shared_pathways_NES_comparison.png", p_shared, width = 9, height = 7, dpi = 300)
cat("Saved: figures/shared_pathways_NES_comparison.png\n\n")

cat("=== STEP 7 DONE === Now Source 08_figures_KEGG_barplots.R\n")
