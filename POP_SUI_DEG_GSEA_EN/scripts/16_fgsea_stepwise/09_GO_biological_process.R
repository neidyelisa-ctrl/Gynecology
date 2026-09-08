# =============================================================================
# STEP 9 of 13 - GO Biological Process enrichment via fgsea
# =============================================================================
# RUN ORDER: needs 01, 04 and 05 already sourced in this same R session
# (uses: ranked_full/ranked_genes_pop from step 4, t_obs_sui from step 5,
# kegg_label is not used here but select/org.Hs.eg.db from step 1 are).
# gene_sets_go built here is also needed later by step 11 (ORA).
#
# WHY THIS STEP EXISTS: KEGG (steps 4-8) covers ~218 broad pathways. GO
# Biological Process has far more specific terms - including some that map
# directly onto a connective-tissue/extracellular-matrix weakness hypothesis
# (e.g. "collagen fibril organization", "extracellular matrix organization")
# that KEGG has no equivalent for.
#
# OUTPUT:
#   results/POP_fgsea_GO_BP_full.csv
#   results/SUI_fgsea_GO_BP_full.csv
#   results/shared_GO_BP_FDR025.csv (only if any term is shared)
#   figures/GO_BP_barplot_POP.png
#   figures/GO_BP_barplot_SUI.png
# =============================================================================

setwd("E:/POP+SUI ONE")  # safe to repeat - keeps this correct even in a fresh R session

cat("================ STEP 9: GO Biological Process enrichment ================\n\n")

## --- build GO Biological Process gene sets from org.Hs.eg.db ---------------
cat("Building GO Biological Process gene sets from org.Hs.eg.db...\n")
ann_go_pop <- suppressWarnings(select(org.Hs.eg.db, keys = ranked_genes_pop, keytype = "SYMBOL", columns = c("GO", "ONTOLOGY")))
ann_go_pop <- ann_go_pop[!is.na(ann_go_pop$GO) & ann_go_pop$ONTOLOGY == "BP", ]
gs_sizes_go <- table(ann_go_pop$GO)
valid_go <- names(gs_sizes_go)[gs_sizes_go >= 5 & gs_sizes_go <= 200]
gene_sets_go <- split(ann_go_pop$SYMBOL[ann_go_pop$GO %in% valid_go], ann_go_pop$GO[ann_go_pop$GO %in% valid_go])
cat("GO Biological Process terms tested (5-200 genes):", length(gene_sets_go), "\n\n")

## go_term_name(): translates a GO ID into its readable name - GO.db ships
## with org.Hs.eg.db, no internet needed.
go_term_name <- function(ids) {
  nm <- suppressMessages(tryCatch(AnnotationDbi::Term(ids), error = function(e) rep(NA_character_, length(ids))))
  unname(ifelse(is.na(nm), ids, paste0(ids, " - ", nm)))
}

## --- run fgsea (GO Biological Process) on POP ------------------------------
cat("Running fgsea (GO Biological Process) on POP...\n")
set.seed(208261)
go_pop <- as.data.frame(fgsea(pathways = gene_sets_go, stats = ranked_full, minSize = 5, maxSize = 200))
go_pop$GOName <- go_term_name(go_pop$pathway)
go_pop$leadingEdge <- sapply(go_pop$leadingEdge, paste, collapse = "/")
go_pop <- go_pop[order(go_pop$pval), ]
write.csv(go_pop, "results/POP_fgsea_GO_BP_full.csv", row.names = FALSE)
cat("Terms significant at FDR<0.25:", sum(go_pop$padj < 0.25, na.rm = TRUE), "of", nrow(go_pop), "\n")
cat("Terms significant at FDR<0.05:", sum(go_pop$padj < 0.05, na.rm = TRUE), "\n\n")

## --- run fgsea (GO Biological Process) on SUI ------------------------------
cat("Running fgsea (GO Biological Process) on SUI...\n")
set.seed(2020)
go_sui <- as.data.frame(fgsea(pathways = gene_sets_go, stats = t_obs_sui, minSize = 5, maxSize = 200))
go_sui$GOName <- go_term_name(go_sui$pathway)
go_sui$leadingEdge <- sapply(go_sui$leadingEdge, paste, collapse = "/")
go_sui <- go_sui[order(go_sui$pval), ]
write.csv(go_sui, "results/SUI_fgsea_GO_BP_full.csv", row.names = FALSE)
cat("Terms significant at FDR<0.25:", sum(go_sui$padj < 0.25, na.rm = TRUE), "of", nrow(go_sui), "\n")
cat("Terms significant at FDR<0.05:", sum(go_sui$padj < 0.05, na.rm = TRUE), "\n\n")

## --- shared GO BP terms between POP and SUI --------------------------------
report_shared_go <- function(fdr_cut) {
  pop_sig <- subset(go_pop, padj < fdr_cut)
  sui_sig <- subset(go_sui, padj < fdr_cut)
  shared_ids <- intersect(pop_sig$pathway, sui_sig$pathway)
  cat("--- FDR <", fdr_cut, ": POP significant =", nrow(pop_sig),
      "| SUI significant =", nrow(sui_sig), "| SHARED =", length(shared_ids), "---\n")
  if (length(shared_ids) == 0) return(data.frame())
  out <- merge(pop_sig[pop_sig$pathway %in% shared_ids, c("pathway", "GOName", "NES", "padj")],
               sui_sig[sui_sig$pathway %in% shared_ids, c("pathway", "NES", "padj")],
               by = "pathway", suffixes = c("_POP", "_SUI"))
  out$Same_direction <- sign(out$NES_POP) == sign(out$NES_SUI)
  out <- out[order(out$padj_POP), ]
  cat("Same direction in both diseases:", sum(out$Same_direction, na.rm = TRUE), "of", nrow(out), "shared terms\n\n")
  out
}
shared_go_025 <- report_shared_go(0.25)
if (nrow(shared_go_025) > 0) {
  write.csv(shared_go_025, "results/shared_GO_BP_FDR025.csv", row.names = FALSE)
}

## --- direct check on connective-tissue/ECM hypothesis terms ---------------
cat("Checking specific GO terms tied to a connective-tissue/ECM weakness\n")
cat("hypothesis (present here if they had 5-200 genes in the ranked list):\n\n")
ecm_go_terms <- c(
  "GO:0030198",  # extracellular matrix organization
  "GO:0030199",  # collagen fibril organization
  "GO:0030574",  # collagen catabolic process
  "GO:0007044",  # cell-substrate junction assembly
  "GO:0031012"   # extracellular matrix (a CC term, kept in case also annotated BP)
)
for (gid in ecm_go_terms) {
  row_pop <- go_pop[go_pop$pathway == gid, ]
  row_sui <- go_sui[go_sui$pathway == gid, ]
  if (nrow(row_pop) == 1 || nrow(row_sui) == 1) {
    cat(go_term_name(gid), ":\n")
    if (nrow(row_pop) == 1) {
      cat("  POP: NES =", round(row_pop$NES, 2), ", FDR =", signif(row_pop$padj, 3), "\n")
    } else {
      cat("  POP: not tested here (fewer than 5 or more than 200 genes in the ranked list)\n")
    }
    if (nrow(row_sui) == 1) {
      cat("  SUI: NES =", round(row_sui$NES, 2), ", FDR =", signif(row_sui$padj, 3), "\n")
    } else {
      cat("  SUI: not tested here (fewer than 5 or more than 200 genes in the ranked list)\n")
    }
    cat("\n")
  }
}

## --- figure: top GO Biological Process terms in POP and SUI ---------------
make_go_barplot <- function(go_res, title, n_top = 15) {
  go_res <- go_res[!is.na(go_res$NES), ]
  d <- head(go_res[order(go_res$pval), ], n_top)
  d$Sig <- ifelse(d$padj < 0.05, "FDR<0.05", ifelse(d$padj < 0.25, "FDR<0.25", "NS"))
  d$Label <- factor(d$GOName, levels = rev(d$GOName))
  ggplot(d, aes(x = NES, y = Label, fill = Sig)) +
    geom_col() +
    scale_fill_manual(values = c("FDR<0.05" = "#B2182B", "FDR<0.25" = "#F4A582", "NS" = "grey70")) +
    geom_vline(xintercept = 0, color = "grey30") +
    labs(title = title, x = "Normalized Enrichment Score (NES)", y = NULL, fill = "Significance") +
    theme_bw() + theme(axis.text.y = element_text(size = 7), plot.title = element_text(size = 12))
}
ggsave("figures/GO_BP_barplot_POP.png", make_go_barplot(go_pop, "GO Biological Process (fgsea) - top terms in POP"), width = 12, height = 6.5, dpi = 300)
cat("Saved: figures/GO_BP_barplot_POP.png\n")
ggsave("figures/GO_BP_barplot_SUI.png", make_go_barplot(go_sui, "GO Biological Process (fgsea) - top terms in SUI"), width = 12, height = 6.5, dpi = 300)
cat("Saved: figures/GO_BP_barplot_SUI.png\n\n")

cat("With thousands of GO terms tested, expect the 'shared between POP and\n")
cat("SUI' list to look thin even when both diseases show real, biologically\n")
cat("coherent signal individually - requiring FDR<0.25 in BOTH datasets\n")
cat("independently is a much harder bar with ~2,000-4,000 tests than with\n")
cat("KEGG's 218. Check each term of interest in the POP-only and SUI-only\n")
cat("files too, not only the 'shared' file, before concluding there is no\n")
cat("overlap.\n\n")

cat("=== STEP 9 DONE === Now Source 10_GSVA.R\n")
