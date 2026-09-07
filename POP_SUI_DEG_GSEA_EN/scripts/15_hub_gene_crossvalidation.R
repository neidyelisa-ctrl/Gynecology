# =============================================================================
# Hub-gene cross-validation - do the STRING/Cytoscape hub genes agree with
# the functional enrichment already established for POP and SUI?
# =============================================================================
# INPUT: the hub-gene node tables exported from Cytoscape (stringApp "Export
# table"), and the enrichment TSVs already used in scripts 12-13. IMPORTANT:
# these hub-gene files carry NO degree/MCC/centrality SCORE column - they are
# whichever nodes were selected/exported in Cytoscape, not a formal cytoHubba
# ranking. If a numeric hub ranking is needed for the methods section, re-
# export from Cytoscape's Network Analyzer or cytoHubba WITH those score
# columns included; this script only checks membership, it invents no score.
#
# WHAT THIS DOES: for each hub gene, lists which of the already-established
# enriched terms (script 12/13's STRING output) it belongs to. This is a
# sanity/consistency check, not a new statistical test - no p-value or FDR is
# computed here.
#
# OUTPUT: results/hub_gene_crossvalidation.csv,
#         figures/hub_gene_crossvalidation.png

setwd("E:/POP+SUI 63")  # SAME root as scripts 11-14 - change this one line only

suppressMessages({
  library(ggplot2)
})

dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

read_hub_genes <- function(path) {
  d <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  unique(d[["display name"]])
}

read_string_tsv <- function(path) {
  d <- read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  colnames(d) <- trimws(sub("^#", "", colnames(d)))
  d
}

match_terms <- function(gene, enrich_tab) {
  hit <- grepl(paste0("(^|,)", gene, "(,|$)"), enrich_tab$`matching proteins in your network (labels)`)
  enrich_tab$`term description`[hit]
}

## --- POP -----------------------------------------------------------------
pop_hubs <- read_hub_genes("data/string_hubs/POP_hub_genes.csv")
pop_enrich <- read_string_tsv("data/string_POP/enrichment.all.tsv")
cat("POP hub genes (", length(pop_hubs), "):", paste(pop_hubs, collapse = ", "), "\n\n")

pop_rows <- lapply(pop_hubs, function(g) {
  terms <- match_terms(g, pop_enrich)
  data.frame(Group = "POP", Gene = g, N_enriched_terms_matched = length(terms),
             Matched_terms = paste(terms, collapse = "; "), stringsAsFactors = FALSE)
})
pop_df <- do.call(rbind, pop_rows)

## --- SUI -------------------------------------------------------------------
sui_hubs <- read_hub_genes("data/string_hubs/SUI_hub_genes.csv")
sui_enrich <- rbind(
  read_string_tsv("data/string_SUI/enrichment.KEGG.tsv")[, c("term description", "matching proteins in your network (labels)")],
  read_string_tsv("data/string_SUI/enrichment.Process.tsv")[, c("term description", "matching proteins in your network (labels)")],
  read_string_tsv("data/string_SUI/enrichment.Function.tsv")[, c("term description", "matching proteins in your network (labels)")]
)
cat("SUI hub genes (", length(sui_hubs), "):", paste(sui_hubs, collapse = ", "), "\n\n")

sui_rows <- lapply(sui_hubs, function(g) {
  terms <- match_terms(g, sui_enrich)
  data.frame(Group = "SUI", Gene = g, N_enriched_terms_matched = length(terms),
             Matched_terms = paste(terms, collapse = "; "), stringsAsFactors = FALSE)
})
sui_df <- do.call(rbind, sui_rows)

## --- combine, save, summarize -----------------------------------------------
all_df <- rbind(pop_df, sui_df)
all_df$Confirmed_by_enrichment <- all_df$N_enriched_terms_matched > 0
write.csv(all_df, "results/hub_gene_crossvalidation.csv", row.names = FALSE)

cat("=== Summary ===\n")
cat("POP:", sum(pop_df$N_enriched_terms_matched > 0), "of", nrow(pop_df),
    "hub genes fall inside at least one already-enriched term\n")
cat("SUI:", sum(sui_df$N_enriched_terms_matched > 0), "of", nrow(sui_df),
    "hub genes fall inside at least one already-enriched term\n\n")
cat("Genes NOT matching any enriched term (real, not a failure - these are\n")
cat("PPI-network hubs by CONNECTIVITY, which is a different criterion than\n")
cat("enrichment membership; worth checking individually in GeneCards/STRING):\n")
print(all_df[!all_df$Confirmed_by_enrichment, c("Group", "Gene")], row.names = FALSE)
cat("\n")

## --- figure: how many enriched terms does each hub gene match? -------------
all_df$Gene <- factor(all_df$Gene, levels = all_df$Gene[order(all_df$Group, all_df$N_enriched_terms_matched)])
p <- ggplot(all_df, aes(x = N_enriched_terms_matched, y = Gene, fill = Group)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = ifelse(N_enriched_terms_matched == 0, "no enriched term", N_enriched_terms_matched)),
            hjust = -0.1, size = 3, color = "grey30") +
  scale_fill_manual(values = c(POP = "#eb6834", SUI = "#2a78d6")) +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.25))) +
  facet_wrap(~Group, scales = "free_y", ncol = 1) +
  labs(title = "Hub genes vs. already-established functional enrichment",
       subtitle = "Number of enriched STRING terms each hub gene belongs to (POP: all categories; SUI: KEGG+Process+Function)",
       x = "Enriched terms matched", y = NULL) +
  theme_bw() +
  theme(legend.position = "none", strip.text = element_text(size = 10, face = "bold"),
        plot.title = element_text(size = 12), plot.subtitle = element_text(size = 9, color = "grey30"))
ggsave("figures/hub_gene_crossvalidation.png", p, width = 10, height = 7, dpi = 300)
cat("Saved: figures/hub_gene_crossvalidation.png\n")
cat("Saved: results/hub_gene_crossvalidation.csv\n")
