# =============================================================================
# STEP 3 of 13 - differentially expressed genes (DEG) in POP
# =============================================================================
# RUN ORDER: needs 01 and 02 already sourced in this same R session (uses:
# counts_all, meta, group_full from step 2).
#
# WHAT THIS STEP DOES: runs DESeq2 (primary method) to find genes
# differentially expressed between 12 Control and 12 POP samples, then runs
# edgeR+voom+limma as an independent cross-check - the voom/limma fit is
# also what step 4 (fgsea) ranks genes by.
#
# OUTPUT:
#   results/POP_DESeq2_full_table.csv       - every gene tested, DESeq2
#   results/POP_DEG_logFC1_FDR05.csv        - DEG only (|log2FC|>1, FDR<0.05)
#   results/POP_voom_limma_full_table.csv   - every gene tested, voom/limma
# =============================================================================

setwd("E:/POP+SUI FGSEA")  # safe to repeat - keeps this correct even in a fresh R session

cat("================ STEP 3: differentially expressed genes (DEG) in POP ================\n\n")

counts_int <- counts_all
storage.mode(counts_int) <- "integer"  # DESeq2 requires integer counts

coldata_naive <- data.frame(row.names = meta$GSM, group = group_full)
dds_naive <- DESeqDataSetFromMatrix(countData = counts_int, colData = coldata_naive, design = ~group)

# Low-count filter: keep genes with at least 10 reads in at least 12 of the
# 24 samples (the smaller group's size) - DESeq2's own recommended practice.
dds_naive <- dds_naive[rowSums(counts(dds_naive) >= 10) >= 12, ]
cat("Genes kept after the low-count filter:", nrow(dds_naive), "\n")

dds_naive <- DESeq(dds_naive, quiet = TRUE)
res_naive <- as.data.frame(results(dds_naive, contrast = c("group", "POP", "Control"), alpha = 0.05))
res_naive$Gene <- rownames(res_naive)
res_naive <- res_naive[order(res_naive$pvalue), ]
pop_full <- res_naive[, c("Gene", "log2FoldChange", "baseMean", "stat", "pvalue", "padj")]
colnames(pop_full) <- c("Gene", "logFC", "baseMean", "stat", "P.Value", "adj.P.Val")
write.csv(pop_full, "results/POP_DESeq2_full_table.csv", row.names = FALSE)

# DEG = |log2FC|>1 AND FDR<0.05.
pop_deg <- subset(pop_full, !is.na(adj.P.Val) & adj.P.Val < 0.05 & abs(logFC) > 1)
pop_deg <- pop_deg[order(pop_deg$adj.P.Val), ]
write.csv(pop_deg, "results/POP_DEG_logFC1_FDR05.csv", row.names = FALSE)

cat("\n=== RESULT: DEG in POP (12 Control vs 12 POP) ===\n")
cat("Genes tested:", nrow(pop_full), "\n")
cat("DEG (|log2FC|>1, FDR<0.05):", nrow(pop_deg), "(",
    sum(pop_deg$logFC > 0), "up /", sum(pop_deg$logFC < 0), "down )\n\n")

# edgeR+voom+limma: a second, independent statistical method, run as a
# cross-check and to produce the moderated-t ranking fgsea uses in step 4.
dge <- DGEList(counts = counts_all, group = group_full)
design_naive_lm <- model.matrix(~group_full)
keep_expr <- filterByExpr(dge, design_naive_lm)
dge <- dge[keep_expr, , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge, method = "TMM")
voom_fit <- voom(dge, design_naive_lm)
fit_pop <- eBayes(lmFit(voom_fit, design_naive_lm))
pop_full_limma <- topTable(fit_pop, coef = "group_fullPOP", number = Inf, sort.by = "P")
pop_full_limma$Gene <- rownames(pop_full_limma)
write.csv(pop_full_limma[, c("Gene","logFC","AveExpr","t","P.Value","adj.P.Val")],
          "results/POP_voom_limma_full_table.csv", row.names = FALSE)
pop_deg_limma <- subset(pop_full_limma, adj.P.Val < 0.05 & abs(logFC) > 1)
cat("Cross-check with edgeR+voom+limma (secondary method):", nrow(pop_deg_limma), "DEG\n\n")

cat("=== STEP 3 DONE === Now Source 04_POP_fgsea_KEGG.R\n")
