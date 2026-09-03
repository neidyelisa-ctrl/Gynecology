# =============================================================================
# ROC / AUC screen of candidate biomarker genes - POP (GSE208261, 12v12)
# =============================================================================
# Standalone script. Re-derives the exact same 12-Control-vs-12-POP DESeq2
# model as script 05 (same filtering, same design ~group) so the candidate
# gene list is read directly off real, already-reported DEG results, not a
# separately re-fit or re-tuned model.
#
# WHAT THIS DOES: for the top individually significant genes from the
# 12x12 DEG list (|log2FC|>1, FDR<0.05), computes how well each gene alone
# discriminates POP vs Control by ROC/AUC, and gives one honest read on
# whether a small combined panel could work as a diagnostic candidate.
#
# WHAT THIS DELIBERATELY DOES NOT CLAIM:
#   - AUC computed on the same 24 samples used to pick and fit the genes is
#     an OPTIMISTIC, in-sample estimate, not a validated biomarker
#     performance. The panel score below is evaluated with leave-one-out
#     cross-validation (LOOCV) specifically to be honest about this, but
#     even LOOCV here is not full nested CV (the gene LIST itself was
#     chosen using all 24 samples, not re-selected inside each fold) - so
#     even the LOOCV AUC is still an upper bound, not a validated estimate.
#     True validation requires an independent cohort. This caveat is
#     printed with every number below and belongs next to any of these
#     numbers if quoted in the thesis.
#   - n=24 is far too small for a clinical AUC confidence interval to be
#     meaningful; only point estimates are reported.
# =============================================================================

suppressMessages({
  library(DESeq2)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(ggplot2)
})

dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

## --- STEP 1: rebuild the exact 12x12 naive DESeq2 model from script 05 ----
meta <- read.csv("data/GSE208261_sample_metadata.csv")
counts_raw <- read.delim("data/GSE208261_raw_counts.tsv", row.names = 1, check.names = FALSE)
counts_all <- as.matrix(counts_raw[, meta$GSM])
stopifnot(identical(colnames(counts_all), meta$GSM))

ann_id <- suppressWarnings(select(org.Hs.eg.db, keys = rownames(counts_all), keytype = "ENTREZID", columns = "SYMBOL"))
ann_id <- ann_id[!is.na(ann_id$SYMBOL) & !duplicated(ann_id$ENTREZID), ]
counts_all <- counts_all[rownames(counts_all) %in% ann_id$ENTREZID, ]
rownames(counts_all) <- ann_id$SYMBOL[match(rownames(counts_all), ann_id$ENTREZID)]
counts_all <- rowsum(counts_all, group = rownames(counts_all))

group_full <- factor(meta$Treatment, levels = c("Control", "POP"))
counts_int <- counts_all; storage.mode(counts_int) <- "integer"
coldata_naive <- data.frame(row.names = meta$GSM, group = group_full)
dds_naive <- DESeqDataSetFromMatrix(countData = counts_int, colData = coldata_naive, design = ~group)
dds_naive <- dds_naive[rowSums(counts(dds_naive) >= 10) >= 12, ]
dds_naive <- DESeq(dds_naive, quiet = TRUE)
res_naive <- as.data.frame(results(dds_naive, contrast = c("group", "POP", "Control"), alpha = 0.05))
res_naive$Gene <- rownames(res_naive)
cat("Re-derived 12x12 DESeq2 model: ", nrow(res_naive), "genes tested (must match script 05's",
    "results/16_POP_GSE208261_FULL24_naive_DESeq2_full.csv - checked below)\n")

check <- read.csv("results/16_POP_GSE208261_FULL24_naive_DESeq2_full.csv")
stopifnot(nrow(check) == nrow(res_naive))
cat("Consistency check OK: same number of tested genes as script 05's saved result.\n\n")

## --- STEP 2: candidate gene list - top significant, protein-coding DEGs ---
sig <- subset(res_naive, !is.na(padj) & padj < 0.05 & abs(log2FoldChange) > 1)
sig <- sig[order(sig$padj), ]
is_probably_ncrna <- grepl("^(RNU|MIR|LINC|LOC|SNOR|RNA5|RNVU)", sig$Gene)
sig_coding <- sig[!is_probably_ncrna, ]
n_candidates <- min(10, nrow(sig_coding))
candidates <- head(sig_coding$Gene, n_candidates)
cat("Candidate genes (top", n_candidates, "significant, named protein-coding DEGs, FDR<0.05, |log2FC|>1):\n")
print(sig_coding[seq_len(n_candidates), c("Gene", "log2FoldChange", "padj")])
cat("\n(", sum(is_probably_ncrna), "of", nrow(sig), "significant DEGs were excluded as non-coding",
    "RNA symbols (RNU/MIR/LINC/LOC/SNOR prefixes) - not biologically interpretable as a protein",
    "biomarker candidate, though they remain valid DEGs in the main results.)\n\n")

## --- STEP 3: variance-stabilized expression for ROC (VST, not raw counts) -
vst_mat <- assay(vst(dds_naive, blind = FALSE))
labels <- ifelse(group_full == "POP", 1L, 0L)  # POP = positive class

## --- STEP 4: per-gene ROC/AUC (rank-based, exact, no external package) ----
## AUC for a continuous score vs a binary label is exactly the
## Mann-Whitney U statistic scaled to [0,1] (Hanley & McNeil 1982) -
## this avoids depending on an internet-installed ROC package while being
## numerically identical to what pROC::roc()$auc would return.
compute_roc <- function(score, label) {
  ord <- order(score, decreasing = TRUE)
  score <- score[ord]; label <- label[ord]
  P <- sum(label == 1); N <- sum(label == 0)
  tpr <- cumsum(label == 1) / P
  fpr <- cumsum(label == 0) / N
  tpr <- c(0, tpr); fpr <- c(0, fpr)
  auc <- sum(diff(fpr) * (head(tpr, -1) + tail(tpr, -1)) / 2)
  # if score is inversely related to label, AUC<0.5; report the
  # direction-agnostic discriminative power (max of the two orientations)
  if (auc < 0.5) {
    auc <- 1 - auc
    tpr <- 1 - tpr; fpr <- 1 - fpr
    o2 <- order(fpr)
    fpr <- fpr[o2]; tpr <- tpr[o2]
  }
  list(auc = auc, roc = data.frame(FPR = fpr, TPR = tpr))
}

roc_table <- data.frame(Gene = character(), logFC = numeric(), padj = numeric(), AUC = numeric())
roc_curves <- list()
for (g in candidates) {
  sc <- vst_mat[g, ]
  r <- compute_roc(sc, labels)
  roc_table <- rbind(roc_table, data.frame(
    Gene = g,
    logFC = sig_coding$log2FoldChange[sig_coding$Gene == g],
    padj = sig_coding$padj[sig_coding$Gene == g],
    AUC = round(r$auc, 3)
  ))
  roc_curves[[g]] <- cbind(r$roc, Gene = g)
}
roc_table <- roc_table[order(-roc_table$AUC), ]
write.csv(roc_table, "results/23_ROC_AUC_candidate_genes.csv", row.names = FALSE)
cat("=== Single-gene AUC (POP vs Control, n=24, IN-SAMPLE - see caveats in header) ===\n")
print(roc_table)
cat("\n")

## --- Figure: ROC curves, top 6 genes by AUC --------------------------------
top6 <- head(roc_table$Gene, 6)
curves_df <- do.call(rbind, roc_curves[top6])
auc_labels <- setNames(paste0(top6, " (AUC=", roc_table$AUC[match(top6, roc_table$Gene)], ")"), top6)
curves_df$GeneLabel <- factor(auc_labels[curves_df$Gene], levels = auc_labels[top6])
p_roc <- ggplot(curves_df, aes(x = FPR, y = TPR, color = GeneLabel)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey60") +
  geom_step(linewidth = 0.9) +
  coord_equal() +
  labs(title = "ROC curves - top single-gene candidates, POP vs Control (GSE208261, n=24)",
       subtitle = "In-sample AUC - exploratory screen, not a validated biomarker estimate (see script header)",
       x = "1 - Specificity (FPR)", y = "Sensitivity (TPR)", color = NULL) +
  theme_bw() + theme(legend.position = "right")
ggsave("figures/32_ROC_curves_candidate_genes.png", p_roc, width = 9, height = 7, dpi = 300)
cat("Saved: figures/32_ROC_curves_candidate_genes.png\n\n")

## --- STEP 5: combined panel score, evaluated with LOOCV (honest estimate) -
## Score = mean of sign-adjusted z-scores across the top 5 candidate genes
## (signed so that "higher score" always means "more POP-like", per each
## gene's own DEG direction). Evaluated two ways for direct comparison:
##   (a) in-sample AUC (optimistic, includes the same samples used to
##       standardize/select genes)
##   (b) leave-one-out cross-validated AUC (each sample's score computed
##       using only the OTHER 23 samples' mean/SD for standardization -
##       the gene LIST itself is still fixed from the full-data DEG call,
##       so this is a partial, not full, cross-validation - see header)
panel_genes <- head(roc_table$Gene, 5)
panel_sign <- sign(sig_coding$log2FoldChange[match(panel_genes, sig_coding$Gene)])

score_in_sample <- rowMeans(sapply(seq_along(panel_genes), function(i) {
  x <- vst_mat[panel_genes[i], ]
  panel_sign[i] * scale(x)[, 1]
}))
auc_in_sample <- compute_roc(score_in_sample, labels)$auc

score_loocv <- numeric(length(labels))
for (i in seq_along(labels)) {
  s_i <- sapply(seq_along(panel_genes), function(j) {
    x <- vst_mat[panel_genes[j], ]
    mu <- mean(x[-i]); sdv <- sd(x[-i])
    panel_sign[j] * (x[i] - mu) / sdv
  })
  score_loocv[i] <- mean(s_i)
}
auc_loocv <- compute_roc(score_loocv, labels)$auc

cat("=== Combined 5-gene panel (", paste(panel_genes, collapse = ", "), ") ===\n")
cat("In-sample AUC (optimistic):", round(auc_in_sample, 3), "\n")
cat("Leave-one-out cross-validated AUC (more honest, still not fully nested):",
    round(auc_loocv, 3), "\n")
cat("Gap between the two =", round(auc_in_sample - auc_loocv, 3),
    "- this gap IS the overfitting this script is trying to be transparent about.\n\n")

panel_roc_in <- compute_roc(score_in_sample, labels)$roc; panel_roc_in$Type <- paste0("In-sample (AUC=", round(auc_in_sample,3), ")")
panel_roc_loo <- compute_roc(score_loocv, labels)$roc; panel_roc_loo$Type <- paste0("LOOCV (AUC=", round(auc_loocv,3), ")")
panel_df <- rbind(panel_roc_in, panel_roc_loo)
p_panel <- ggplot(panel_df, aes(x = FPR, y = TPR, color = Type)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey60") +
  geom_step(linewidth = 1) +
  coord_equal() +
  labs(title = "5-gene combined panel - in-sample vs leave-one-out cross-validated AUC",
       subtitle = paste0("Genes: ", paste(panel_genes, collapse = ", "),
                          " | gap between curves = overfitting, not signal"),
       x = "1 - Specificity (FPR)", y = "Sensitivity (TPR)", color = NULL) +
  theme_bw() + theme(legend.position = "top", plot.title = element_text(size = 13))
ggsave("figures/33_ROC_panel_score.png", p_panel, width = 9.5, height = 8, dpi = 300)
cat("Saved: figures/33_ROC_panel_score.png\n\n")

panel_summary <- data.frame(
  Genes = paste(panel_genes, collapse = ";"),
  AUC_in_sample = round(auc_in_sample, 3),
  AUC_LOOCV = round(auc_loocv, 3)
)
write.csv(panel_summary, "results/24_ROC_panel_summary.csv", row.names = FALSE)

cat("=== DONE. Files written: ===\n")
cat("results/23_ROC_AUC_candidate_genes.csv\n")
cat("results/24_ROC_panel_summary.csv\n")
cat("figures/32_ROC_curves_candidate_genes.png\n")
cat("figures/33_ROC_panel_score.png\n")
