# =============================================================================
# Regenerate KEGG pathway names in saved results + figures, using the full,
# audited lookup table (data/kegg_pathway_names.csv, 223 IDs), WITHOUT
# recomputing any statistic.
# =============================================================================
# Every number in every results/*.csv file (NES, p-value, FDR, DEG calls,
# concordance) is untouched by this script - it only re-derives the
# human-readable PathwayName column from the numeric KEGG PATH ID (which was
# always correct) and redraws the plots that display that name, using the
# exact same plotting code as the script that first produced each figure.
# Run this any time data/kegg_pathway_names.csv is edited, instead of
# re-running the ~10-20 minute permutation loop in scripts 01-06.
# =============================================================================

suppressMessages({
  library(ggplot2)
  library(ggrepel)
})

kegg_tab <- read.csv("data/kegg_pathway_names.csv", colClasses = c("character", "character"))
kegg_names <- setNames(kegg_tab$Name, kegg_tab$PATH5)
kegg_label <- function(id) {
  nm <- kegg_names[id]
  ifelse(is.na(nm), paste0("KEGG ", id), paste0("KEGG ", id, " - ", nm))
}

make_gsea_barplot <- function(gsea_res, title, n_top = 15) {
  gsea_res <- gsea_res[!is.na(gsea_res$NES), ]
  d <- head(gsea_res[order(gsea_res$pvalue), ], n_top)
  d$Sig <- ifelse(d$p.adjust < 0.05, "FDR<0.05", ifelse(d$p.adjust < 0.25, "FDR<0.25", "NS"))
  d$Label <- factor(d$PathwayName, levels = rev(d$PathwayName))
  ggplot(d, aes(x = NES, y = Label, fill = Sig)) +
    geom_col() +
    scale_fill_manual(values = c("FDR<0.05" = "#B2182B", "FDR<0.25" = "#F4A582", "NS" = "grey70")) +
    geom_vline(xintercept = 0, color = "grey30") +
    labs(title = title, x = "Normalized Enrichment Score (NES)", y = NULL, fill = "Significance") +
    theme_bw() + theme(axis.text.y = element_text(size = 8))
}

make_shared_plot <- function(shared_df, title, x_lab, y_lab) {
  ggplot(shared_df, aes(x = NES_POP, y = NES_SUI, label = PathwayName)) +
    geom_hline(yintercept = 0, color = "grey70") + geom_vline(xintercept = 0, color = "grey70") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
    geom_point(aes(color = Same_direction), size = 3) +
    scale_color_manual(values = c("TRUE" = "#2166AC", "FALSE" = "#B2182B"),
                        labels = c("TRUE" = "Same direction", "FALSE" = "Opposite direction")) +
    geom_text_repel(size = 3, max.overlaps = 30) +
    labs(title = title, x = x_lab, y = y_lab, color = NULL) +
    theme_bw() + theme(legend.position = "top")
}

# --- Re-derive PathwayName in every saved GSEA table and re-save + re-plot -
gsea_jobs <- list(
  list(csv = "results/03_GSEA_classic_POP_KEGG.csv",
       fig = "figures/05_GSEA_barplot_POP.png",
       title = "GSEA (classic) - top KEGG pathways in POP"),
  list(csv = "results/04_GSEA_preranked_SUI_KEGG.csv",
       fig = "figures/06_GSEA_barplot_SUI.png",
       title = "GSEA (preranked) - top KEGG pathways in SUI"),
  list(csv = "results/07_GSEA_classic_POP_GSE208261_KEGG.csv",
       fig = "figures/12_GSEA_barplot_POP_GSE208261.png",
       title = "GSEA (classic) - top KEGG pathways in POP (GSE208261, RNA-seq)"),
  list(csv = "results/10_GSEA_preranked_Chen2006_KEGG.csv",
       fig = "figures/16_GSEA_barplot_SUI_Chen2006.png",
       title = "GSEA (preranked) - KEGG pathways in SUI (Chen 2006)\n(exploratory - only 19 pathways testable in a 59-gene panel)"),
  list(csv = "results/14_GSEA_classic_POP_GSE267852VFB_KEGG.csv",
       fig = "figures/21_GSEA_barplot_POP_GSE267852VFB.png",
       title = "GSEA (classic) - top KEGG pathways in POP (GSE267852, VFB)"),
  list(csv = "results/17_GSEA_classic_POP_GSE208261_FULL24_KEGG.csv",
       fig = "figures/26_GSEA_barplot_POP_GSE208261_FULL24.png",
       title = "GSEA (classic) - top KEGG pathways in POP (GSE208261, FULL 12v12)"),
  list(csv = "results/21_GSEA_classic_POP_GSE208261_6v6_POPY_KEGG.csv",
       fig = "figures/28_GSEA_barplot_POP_GSE208261_6v6_POPY.png",
       title = "GSEA (classic) - POP (GSE208261, 6v6, POP_Y half)"),
  list(csv = "results/21_GSEA_classic_POP_GSE208261_6v6_POPD_KEGG.csv",
       fig = "figures/29_GSEA_barplot_POP_GSE208261_6v6_POPD.png",
       title = "GSEA (classic) - POP (GSE208261, 6v6, POP_D half)")
)

for (job in gsea_jobs) {
  d <- read.csv(job$csv, colClasses = c(PATH = "character"))
  d$PathwayName <- kegg_label(d$PATH)
  write.csv(d, job$csv, row.names = FALSE)
  p <- make_gsea_barplot(d, job$title) + theme(plot.title = element_text(size = 12))
  ggsave(job$fig, p, width = 12, height = 6.5, dpi = 300)
  cat("Updated:", job$csv, "+", job$fig, "\n")
}

# --- Same for shared-pathway comparison tables + scatter plots -------------
shared_jobs <- list(
  list(csv = "results/05_shared_pathways_FDR025.csv",
       fig = "figures/07_shared_pathways_NES_comparison.png",
       title = "Shared pathways (FDR<0.25): NES in POP vs NES in SUI",
       x = "NES - POP (classic GSEA)", y = "NES - SUI (preranked GSEA)"),
  list(csv = "results/08_shared_pathways_GSE208261xSUI_FDR025.csv",
       fig = "figures/13_shared_pathways_NES_comparison_GSE208261.png",
       title = "Shared pathways (FDR<0.25): NES in POP (GSE208261, RNA-seq) vs NES in SUI",
       x = "NES - POP GSE208261 (classic GSEA)", y = "NES - SUI (preranked GSEA)"),
  list(csv = "results/15_shared_pathways_GSE267852VFBxSUI_FDR025.csv",
       fig = "figures/22_shared_pathways_NES_comparison_GSE267852VFB.png",
       title = "Shared pathways (FDR<0.25): NES in POP (GSE267852 VFB) vs NES in SUI",
       x = "NES - POP GSE267852 VFB (classic GSEA)", y = "NES - SUI (preranked GSEA)"),
  list(csv = "results/19_shared_pathways_GSE208261FULL24xSUI_FDR025.csv",
       fig = "figures/27_shared_pathways_NES_comparison_GSE208261_FULL24.png",
       title = "Shared pathways (FDR<0.25): NES in POP (GSE208261 FULL24) vs NES in SUI",
       x = "NES - POP GSE208261 FULL24 (classic GSEA)", y = "NES - SUI (preranked GSEA)"),
  list(csv = "results/22_shared_pathways_6v6_POPY_FDR025.csv",
       fig = "figures/30_shared_pathways_6v6_POPY.png",
       title = "Shared pathways (FDR<0.25): NES in POP (GSE208261 6v6, POP_Y) vs NES in SUI",
       x = "NES - POP (classic GSEA, 6v6 POP_Y)", y = "NES - SUI (preranked GSEA)"),
  list(csv = "results/22_shared_pathways_6v6_POPD_FDR025.csv",
       fig = "figures/31_shared_pathways_6v6_POPD.png",
       title = "Shared pathways (FDR<0.25): NES in POP (GSE208261 6v6, POP_D) vs NES in SUI",
       x = "NES - POP (classic GSEA, 6v6 POP_D)", y = "NES - SUI (preranked GSEA)")
)

for (job in shared_jobs) {
  if (!file.exists(job$csv)) { cat("SKIPPED (no file):", job$csv, "\n"); next }
  d <- read.csv(job$csv, colClasses = c(PATH = "character"))
  if (nrow(d) == 0) { cat("SKIPPED (empty):", job$csv, "\n"); next }
  d$PathwayName <- kegg_label(d$PATH)
  write.csv(d, job$csv, row.names = FALSE)
  d_plot <- d[!is.na(d$NES_POP) & !is.na(d$NES_SUI), ]
  if (nrow(d_plot) == 0) { cat("SKIPPED (no comparable NES):", job$fig, "\n"); next }
  p <- make_shared_plot(d_plot, job$title, job$x, job$y)
  ggsave(job$fig, p, width = 9, height = 7, dpi = 300)
  cat("Updated:", job$csv, "+", job$fig, "\n")
}

# results/10_GSEA_preranked_Chen2006_KEGG.csv has no shared-pathway partner
# file (script 3 found 0 shared pathways at FDR<0.25 - figure 18 was never
# generated in the first place, nothing to relabel).

# The gene-direction heatmaps (09, 14, 23) and gene-concordance scatter (17)
# use gene symbols, not KEGG IDs, as labels - unaffected by this table and
# not touched here.

cat("\n=== KEGG label regeneration complete - no statistic was recomputed ===\n")
