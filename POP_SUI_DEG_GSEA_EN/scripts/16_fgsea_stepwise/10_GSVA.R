# =============================================================================
# STEP 10 of 13 - GSVA on KEGG pathways (optional - skips itself if GSVA is
# not installed)
# =============================================================================
# RUN ORDER: needs 01, 02, 03, 04 and 05 already sourced in this same R
# session (uses: has_GSVA/kegg_label from step 1, group_full from step 2,
# voom_fit from step 3, gene_sets_pop from step 4, sui_mat/group_sui/
# gene_sets_sui from step 5).
#
# WHAT THIS STEP ADDS: fgsea (steps 4, 5, 9) all start from a group
# comparison (POP vs Control, SUI vs Ctrl) and ask "is this pathway shifted
# between the two groups". GSVA works the other way around: for EACH
# INDIVIDUAL SAMPLE, it scores how active each pathway is, using only that
# sample's own gene ranking (no group label used in the scoring step). Only
# afterwards are those per-sample scores compared between groups (a plain
# t-test here) - a genuinely different statistical approach, not a re-run
# of GSEA.
#
# OUTPUT: results/POP_GSVA_KEGG.csv, results/SUI_GSVA_KEGG.csv (only if
# GSVA is installed - this step is not part of the fgsea method itself, it
# is an independent cross-check).
# =============================================================================

setwd("E:/POP+SUI FGSEA")  # safe to repeat - keeps this correct even in a fresh R session

cat("================ STEP 10: GSVA on KEGG pathways ================\n\n")

if (!has_GSVA) {
  cat("GSVA is not installed/available in this R session - SKIPPING step 10\n")
  cat("entirely (usually a Windows compilation issue - GSVA's dependencies\n")
  cat("need Rtools to build from source). This does not affect any other step.\n\n")
} else {

## GSVA's function interface changed between versions (older: gsva(expr,
## gene_sets, method="gsva", ...) directly; newer: build a gsvaParam() first,
## then call gsva(param)) - checking which interface this version has.
run_gsva <- function(expr, gene_sets) {
  if (exists("gsvaParam", where = asNamespace("GSVA"), inherits = FALSE)) {
    cat("(using the newer GSVA interface: gsvaParam() + gsva(param))\n")
    param <- GSVA::gsvaParam(expr, gene_sets, kcdf = "Gaussian")
    as.matrix(GSVA::gsva(param, verbose = FALSE))
  } else {
    cat("(using the older GSVA interface: gsva(expr, gene_sets, method='gsva'))\n")
    as.matrix(GSVA::gsva(expr, gene_sets, method = "gsva", kcdf = "Gaussian", verbose = FALSE))
  }
}

gsva_group_test <- function(gsva_mat, group_labels) {
  lv <- levels(group_labels)
  pvals <- apply(gsva_mat, 1, function(x) {
    tryCatch(t.test(x[group_labels == lv[2]], x[group_labels == lv[1]])$p.value, error = function(e) NA)
  })
  data.frame(PATH = rownames(gsva_mat), PathwayName = kegg_label(rownames(gsva_mat)),
             mean_diff = rowMeans(gsva_mat[, group_labels == lv[2], drop = FALSE]) -
                         rowMeans(gsva_mat[, group_labels == lv[1], drop = FALSE]),
             pvalue = pvals, padj = p.adjust(pvals, "BH"))
}

cat("Running GSVA on POP (24 samples x", length(gene_sets_pop), "KEGG pathways)...\n")
gsva_pop <- run_gsva(voom_fit$E, gene_sets_pop)
gsva_pop_res <- gsva_group_test(gsva_pop, group_full)
gsva_pop_res <- gsva_pop_res[order(gsva_pop_res$pvalue), ]
write.csv(gsva_pop_res, "results/POP_GSVA_KEGG.csv", row.names = FALSE)
cat("Pathways significant at FDR<0.25:", sum(gsva_pop_res$padj < 0.25, na.rm = TRUE),
    "of", nrow(gsva_pop_res), "\n")
cat("Pathways significant at FDR<0.05:", sum(gsva_pop_res$padj < 0.05, na.rm = TRUE), "\n\n")

cat("Running GSVA on SUI (6 samples x", length(gene_sets_sui), "KEGG pathways)...\n")
gsva_sui <- run_gsva(sui_mat, gene_sets_sui)
gsva_sui_res <- gsva_group_test(gsva_sui, group_sui)
gsva_sui_res <- gsva_sui_res[order(gsva_sui_res$pvalue), ]
write.csv(gsva_sui_res, "results/SUI_GSVA_KEGG.csv", row.names = FALSE)
cat("Pathways significant at FDR<0.25:", sum(gsva_sui_res$padj < 0.25, na.rm = TRUE),
    "of", nrow(gsva_sui_res), "\n")
cat("Pathways significant at FDR<0.05:", sum(gsva_sui_res$padj < 0.05, na.rm = TRUE), "\n\n")
cat("NOTE: with only 3 vs 3 samples, a per-pathway t-test on SUI has very\n")
cat("little statistical power - treat the SUI GSVA p-values as suggestive,\n")
cat("not conclusive (same caveat as everywhere else in this project on the\n")
cat("SUI side of any comparison).\n\n")

} # end if (has_GSVA)

cat("=== STEP 10 DONE === Now Source 11_ORA_clusterProfiler.R\n")
