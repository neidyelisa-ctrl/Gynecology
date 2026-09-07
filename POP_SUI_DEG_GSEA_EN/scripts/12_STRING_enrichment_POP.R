# =============================================================================
# STRING functional enrichment - POP PPI network
# Reads the enrichment tables exported directly from string-db.org (Multiple
# proteins search on the 163 POP DEG, high-confidence 0.700 interactions) and
# builds publication-ready lollipop charts from them. No numbers are computed
# here - every value plotted (signal, strength, FDR, gene count) is exactly
# what STRING itself reported; this script only visualizes it.
# =============================================================================
#
# INPUT FILES (place in data/string_POP/, exactly as exported from STRING's
# "Exports" -> TSV button on the Enrichment tab):
#   enrichment.all.tsv        - every category together (GO Process, GO
#                                Component, STRING clusters, Monarch,
#                                DISEASES, UniProt Keywords)
#   enrichment.Component.tsv  - GO Component only
#   enrichment.Process.tsv    - GO Process only
#
# OUTPUT: figures/STRING_POP_enrichment_all.png,
#         figures/STRING_POP_enrichment_GO_Component.png,
#         figures/STRING_POP_enrichment_GO_Process.png
#
# Same project root as script 11 (setwd below) - "results/" and "figures/"
# are direct children of it, same convention throughout this project.

setwd("E:/POP+SUI 63")  # SAME root as script 11 - change this one line only

suppressMessages({
  library(ggplot2)
})

dir.create("figures", showWarnings = FALSE)

read_string_tsv <- function(path) {
  d <- read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  colnames(d) <- trimws(sub("^#", "", colnames(d)))
  d
}

# Categorical palette (fixed hue order, validated for colorblind-safe
# adjacent contrast) - one slot per STRING category.
category_colors <- c(
  "GO Process"        = "#2a78d6",  # blue
  "GO Component"      = "#eb6834",  # orange
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

## --- All categories together -----------------------------------------------
all_tab <- read_string_tsv("data/string_POP/enrichment.all.tsv")
cat("Loaded enrichment.all.tsv:", nrow(all_tab), "terms across",
    length(unique(all_tab$category)), "categories\n")

p_all <- make_lollipop(
  all_tab,
  "STRING functional enrichment - POP PPI network",
  paste0(nrow(all_tab), " enriched terms across ", length(unique(all_tab$category)),
         " categories (163 DEG, high-confidence interactions)"),
  color_by_category = TRUE
)
ggsave("figures/STRING_POP_enrichment_all.png", p_all, width = 14,
       height = max(4, 0.42 * nrow(all_tab) + 1.5), dpi = 300)
cat("Saved: figures/STRING_POP_enrichment_all.png (", nrow(all_tab), "terms )\n\n")

## --- GO Component only -------------------------------------------------------
comp_tab <- read_string_tsv("data/string_POP/enrichment.Component.tsv")
comp_tab$category <- "GO Component"
cat("Loaded enrichment.Component.tsv:", nrow(comp_tab), "terms\n")

p_comp <- make_lollipop(
  comp_tab,
  "STRING GO Component enrichment - POP PPI network",
  paste0(nrow(comp_tab), " enriched GO Component terms (163 DEG)"),
  color_by_category = FALSE
)
ggsave("figures/STRING_POP_enrichment_GO_Component.png", p_comp, width = 12,
       height = max(3, 0.5 * nrow(comp_tab) + 1.5), dpi = 300)
cat("Saved: figures/STRING_POP_enrichment_GO_Component.png (", nrow(comp_tab), "terms )\n\n")

## --- GO Process only ----------------------------------------------------------
proc_tab <- read_string_tsv("data/string_POP/enrichment.Process.tsv")
proc_tab$category <- "GO Process"
cat("Loaded enrichment.Process.tsv:", nrow(proc_tab), "term(s)\n")

p_proc <- make_lollipop(
  proc_tab,
  "STRING GO Biological Process enrichment - POP PPI network",
  paste0(nrow(proc_tab), " enriched GO Process term(s) at FDR<0.25 (163 DEG)"),
  color_by_category = FALSE
)
ggsave("figures/STRING_POP_enrichment_GO_Process.png", p_proc, width = 12, height = 3.5, dpi = 300)
cat("Saved: figures/STRING_POP_enrichment_GO_Process.png (", nrow(proc_tab), "term(s) )\n")
cat("Note: only", nrow(proc_tab), "broad GO Process term passed FDR<0.25 - this is real, not a\n")
cat("rendering issue. A 163-gene network from one tissue commonly surfaces few BROAD Process terms\n")
cat("even when specific Component/cluster terms are strong (see the other two figures).\n\n")

cat("=== DONE - 3 figures saved to figures/ ===\n")
