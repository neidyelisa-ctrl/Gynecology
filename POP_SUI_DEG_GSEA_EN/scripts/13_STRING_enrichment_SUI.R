# =============================================================================
# STRING functional enrichment - SUI PPI network
# Reads the enrichment tables exported directly from string-db.org (Multiple
# proteins search on the SUI DEG, high-confidence 0.700 interactions) and
# builds publication-ready lollipop charts from them. No numbers are computed
# here - every value plotted (signal, strength, FDR, gene count) is exactly
# what STRING itself reported; this script only visualizes it.
#
# Companion to 12_STRING_enrichment_POP.R - same method, same chart style,
# so the two sides are visually and statistically comparable side by side.
# =============================================================================
#
# INPUT FILES (place in data/string_SUI/, exactly as exported from STRING's
# "Exports" -> TSV button on the Enrichment tab):
#   enrichment.KEGG.tsv       - KEGG pathways
#   enrichment.Process.tsv    - GO Biological Process
#   enrichment.Function.tsv   - GO Molecular Function
#
# OUTPUT: figures/STRING_SUI_enrichment_KEGG.png,
#         figures/STRING_SUI_enrichment_GO_Process.png,
#         figures/STRING_SUI_enrichment_GO_Function.png
#
# Same project root as scripts 11-12 (setwd below) - "results/" and
# "figures/" are direct children of it, same convention throughout.

setwd("E:/POP+SUI 63")  # SAME root as scripts 11-12 - change this one line only

suppressMessages({
  library(ggplot2)
})

dir.create("figures", showWarnings = FALSE)

read_string_tsv <- function(path) {
  d <- read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  colnames(d) <- trimws(sub("^#", "", colnames(d)))
  d
}

make_lollipop <- function(d, title, subtitle, point_color) {
  d$Label <- paste0(d$`term description`, "  (", d$`term ID`, ")")
  d <- d[order(d$signal), ]
  d$Label <- factor(d$Label, levels = d$Label)

  ggplot(d, aes(x = signal, y = Label)) +
    geom_segment(aes(x = 0, xend = signal, y = Label, yend = Label),
                 color = point_color, linewidth = 1) +
    geom_point(aes(size = `observed gene count`), color = point_color) +
    geom_text(aes(label = paste0("FDR=", formatC(`false discovery rate`, format = "e", digits = 1))),
              hjust = -0.15, size = 2.8, color = "grey30") +
    scale_size_continuous(name = "Genes\nmatched", range = c(2, 7)) +
    scale_x_continuous(expand = expansion(mult = c(0.02, 0.42))) +
    labs(title = title, subtitle = paste(strwrap(subtitle, width = 95), collapse = "\n"),
         x = "Signal (STRING enrichment score)", y = NULL) +
    theme_bw() +
    theme(axis.text.y = element_text(size = 8), plot.title = element_text(size = 12),
          plot.subtitle = element_text(size = 9, color = "grey30"),
          panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
          legend.position = "right")
}

## --- KEGG pathways ----------------------------------------------------------
kegg_tab <- read_string_tsv("data/string_SUI/enrichment.KEGG.tsv")
cat("Loaded enrichment.KEGG.tsv:", nrow(kegg_tab), "pathway(s)\n")

p_kegg <- make_lollipop(
  kegg_tab,
  "STRING KEGG pathway enrichment - SUI PPI network",
  paste0(nrow(kegg_tab), " enriched KEGG pathway(s) at FDR<0.25"),
  point_color = "#1baf7a"
)
ggsave("figures/STRING_SUI_enrichment_KEGG.png", p_kegg,
       width = 12, height = max(3, 0.5 * nrow(kegg_tab) + 1.8), dpi = 300)
cat("Saved: figures/STRING_SUI_enrichment_KEGG.png (", nrow(kegg_tab), "term(s) )\n\n")

## --- GO Process --------------------------------------------------------------
proc_tab <- read_string_tsv("data/string_SUI/enrichment.Process.tsv")
cat("Loaded enrichment.Process.tsv:", nrow(proc_tab), "term(s)\n")

p_proc <- make_lollipop(
  proc_tab,
  "STRING GO Biological Process enrichment - SUI PPI network",
  paste0(nrow(proc_tab), " enriched GO Process term(s) at FDR<0.25"),
  point_color = "#2a78d6"
)
ggsave("figures/STRING_SUI_enrichment_GO_Process.png", p_proc,
       width = 12, height = max(3, 0.5 * nrow(proc_tab) + 1.8), dpi = 300)
cat("Saved: figures/STRING_SUI_enrichment_GO_Process.png (", nrow(proc_tab), "term(s) )\n\n")

## --- GO Molecular Function ----------------------------------------------------
func_tab <- read_string_tsv("data/string_SUI/enrichment.Function.tsv")
cat("Loaded enrichment.Function.tsv:", nrow(func_tab), "term(s)\n")

p_func <- make_lollipop(
  func_tab,
  "STRING GO Molecular Function enrichment - SUI PPI network",
  paste0(nrow(func_tab), " enriched GO Function term(s) at FDR<0.25"),
  point_color = "#eda100"
)
ggsave("figures/STRING_SUI_enrichment_GO_Function.png", p_func,
       width = 12, height = max(3, 0.5 * nrow(func_tab) + 1.8), dpi = 300)
cat("Saved: figures/STRING_SUI_enrichment_GO_Function.png (", nrow(func_tab), "term(s) )\n\n")

## --- Note on overlap with POP ------------------------------------------------
cat("=== Note ===\n")
cat("SUI's top hits (Ribosome, Cytoplasmic translation, Translation,\n")
cat("Structural constituent of ribosome) are all ribosomal/translation-\n")
cat("machinery terms - a DIFFERENT biological axis from POP's cornified-\n")
cat("envelope/skin-barrier signal (see 12_STRING_enrichment_POP.R output).\n")
cat("This matches the fgsea KEGG finding already established for this\n")
cat("project: 'Ribosome' (hsa03010) is one of the pathways significant in\n")
cat("BOTH POP and SUI - the STRING PPI network confirms it independently,\n")
cat("from a completely different method (protein-protein interaction\n")
cat("enrichment instead of gene-set enrichment on the ranked gene list).\n\n")

cat("=== DONE - 3 figures saved to figures/ ===\n")
