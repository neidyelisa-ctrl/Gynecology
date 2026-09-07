# =============================================================================
# STRING functional enrichment - POP+SUI common-gene PPI network ("conjunta")
# Reads the enrichment table exported directly from string-db.org for the 34
# genes common to BOTH the POP and SUI DEG lists (see
# POP_SUI_common_34_genes_STRING.txt) and builds the same lollipop chart used
# for the POP-only (script 12) and SUI-only (script 13) networks, so all
# three are visually and statistically comparable side by side. No numbers
# are computed here - every value plotted is exactly what STRING reported.
# =============================================================================
#
# INPUT FILE (place in data/string_common/, exactly as exported from
# STRING's "Exports" -> TSV button on the Enrichment tab, "all" categories):
#   enrichment.all.tsv
#
# OUTPUT: figures/STRING_common_enrichment_all.png
#
# Same project root as scripts 11-13 (setwd below).

setwd("E:/POP+SUI 63")  # SAME root as scripts 11-13 - change this one line only

suppressMessages({
  library(ggplot2)
})

dir.create("figures", showWarnings = FALSE)

read_string_tsv <- function(path) {
  d <- read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  colnames(d) <- trimws(sub("^#", "", colnames(d)))
  d
}

# Same categorical palette as script 12 (POP), for consistency across all
# three STRING figures in this project.
category_colors <- c(
  "GO Process"        = "#2a78d6",  # blue
  "GO Component"      = "#eb6834",  # orange
  "GO Function"        = "#eda100",  # yellow (kept distinct from Monarch below)
  "KEGG"               = "#1baf7a",  # aqua
  "STRING clusters"   = "#1baf7a",  # aqua
  "Monarch"           = "#eda100",  # yellow
  "DISEASES"          = "#e87ba4",  # magenta
  "UniProt Keywords"  = "#4a3aa7"   # violet
)

make_lollipop <- function(d, title, subtitle, color_by_category = TRUE) {
  d$Label <- paste0(d$`term description`, "  (", d$`term ID`, ")")
  d <- d[order(d$signal), ]
  d$Label <- factor(d$Label, levels = d$Label)

  p <- ggplot(d, aes(x = signal, y = Label))
  if (color_by_category) {
    p <- p + geom_segment(aes(x = 0, xend = signal, y = Label, yend = Label, color = category),
                           linewidth = 1) +
      geom_point(aes(color = category, size = `observed gene count`)) +
      scale_color_manual(values = category_colors, name = "Category")
  } else {
    p <- p + geom_segment(aes(x = 0, xend = signal, y = Label, yend = Label),
                           color = "#2a78d6", linewidth = 1) +
      geom_point(aes(size = `observed gene count`), color = "#2a78d6")
  }
  p +
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

## --- Common (POP+SUI) network, all categories ------------------------------
common_tab <- read_string_tsv("data/string_common/enrichment.all.tsv")
cat("Loaded enrichment.all.tsv:", nrow(common_tab), "term(s) across",
    length(unique(common_tab$category)), "categor(ies)\n")

p_common <- make_lollipop(
  common_tab,
  "STRING functional enrichment - POP+SUI common-gene network",
  paste0(nrow(common_tab), " enriched term(s) at FDR<0.25 (34 genes shared between the POP and SUI DEG lists)"),
  color_by_category = TRUE
)
ggsave("figures/STRING_common_enrichment_all.png", p_common,
       width = 12, height = max(3, 0.5 * nrow(common_tab) + 1.8), dpi = 300)
cat("Saved: figures/STRING_common_enrichment_all.png (", nrow(common_tab), "term(s) )\n\n")

## --- Note -------------------------------------------------------------------
cat("=== Note ===\n")
cat("Only 1 term (Ichthyosis, UniProt KW-0977, FDR=0.048) reaches FDR<0.25\n")
cat("here - expected, not a failure: with only 34 genes in the network,\n")
cat("STRING has very little statistical power left after correcting for\n")
cat("multiple testing across thousands of possible terms.\n\n")
cat("This is the EXACT SAME UniProt term (KW-0977) already significant in\n")
cat("the POP-only network (FDR=2.18e-05, 7 genes: KRT1, KRT10, SDR9C7,\n")
cat("SLC27A4, ALOXE3, ELOVL4, ALOX12B - see 12_STRING_enrichment_POP.R). The\n")
cat("3 genes driving it here (ALOXE3, ELOVL4, ALOX12B) are a direct subset\n")
cat("of those same 7 - not just a similar theme, the identical STRING term\n")
cat("and 3 of the same genes. That is a real, verifiable, targeted finding:\n")
cat("genes dysregulated in BOTH POP and SUI that independently converge on\n")
cat("the same lipid-metabolism/skin-barrier term POP's much larger network\n")
cat("already flagged as significant on its own.\n\n")

cat("=== DONE - 1 figure saved to figures/ ===\n")
