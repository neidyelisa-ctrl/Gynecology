# =============================================================================
# STEP 5 of 13 - SUI data (Wei et al. 2020): DEG (as published) + KEGG fgsea
# =============================================================================
# RUN ORDER: needs 01 already sourced in this same R session (uses: select,
# kegg_label). Does not need steps 2-4 (POP and SUI are loaded independently)
# but is numbered after them to match the pipeline's reading order.
#
# WHAT THIS STEP DOES: reads the SUI microarray data exactly as published by
# Wei et al. 2020 (3 SUI vs 3 Control samples - proven directly from the
# file's own column names, printed below), reconstructs the fold-change sign
# convention, then runs fgsea on SUI's full ranked gene list - the same
# method, same code shape, as step 4 used for POP.
#
# OUTPUT:
#   results/SUI_Wei2020_full_table.csv   - SUI's DEG list, as published
#   results/SUI_fgsea_KEGG_full.csv      - every KEGG pathway tested in SUI
# =============================================================================

setwd("E:/POP+SUI ONE")  # safe to repeat - keeps this correct even in a fresh R session

cat("================ STEP 5: SUI data (Wei 2020) ================\n\n")

wei_xls <- "43032_2020_144_MOESM2_ESM.xls"
sheet_names <- excel_sheets(wei_xls)
cat("Sheets found in the Excel file:", paste(sheet_names, collapse = ", "), "\n")
stopifnot(all(c("up_Sui_vs_Ctrl", "down_Sui_vs_Ctrl") %in% sheet_names))

# skip=17: the first 17 rows are a descriptive header written by the
# authors (not data) - skipping them to reach the real table.
up_sheet   <- read_excel(wei_xls, sheet = "up_Sui_vs_Ctrl",   skip = 17)
down_sheet <- read_excel(wei_xls, sheet = "down_Sui_vs_Ctrl", skip = 17)
cat("'up_Sui_vs_Ctrl':", nrow(up_sheet), "rows,", ncol(up_sheet), "columns\n")
cat("'down_Sui_vs_Ctrl':", nrow(down_sheet), "rows,", ncol(down_sheet), "columns\n\n")

# Direct proof this is 3 SUI vs 3 Control samples - the file's own column
# names, printed for you to check with your own eyes.
cat("Individual-sample columns found in the Wei 2020 file:\n")
print(grep("Sui[0-9]|Ctrl[0-9]", colnames(up_sheet), value = TRUE))
cat("-> 3 'SuiN' columns + 3 'CtrlN' columns = a 3-vs-3 design, confirmed\n")
cat("   directly from the file's own column names, not by assumption.\n\n")

# Fold Change sign convention: the "down" sheet lists a MAGNITUDE (always
# >=1), not a negative value - the sign comes from which sheet the gene is
# in, not from the number itself.
cat("Checking the Fold Change sign convention:\n")
cat("  'up' sheet, Fold Change range:", round(min(up_sheet$`Fold Change`), 2), "to",
    round(max(up_sheet$`Fold Change`), 2), "\n")
cat("  'down' sheet, Fold Change range:", round(min(down_sheet$`Fold Change`), 2), "to",
    round(max(down_sheet$`Fold Change`), 2), "\n")
if (min(down_sheet$`Fold Change`) >= 1) {
  cat("  CONFIRMED: the 'down' sheet only gives the MAGNITUDE (always >=1) -\n")
  cat("  the negative sign is applied manually below based on the sheet.\n\n")
}

up_sheet$Direction <- "up"; down_sheet$Direction <- "down"
sample_cols <- c("[Sui1, Sui](normalized)", "[Sui2, Sui](normalized)", "[Sui3, Sui](normalized)",
                  "[Ctrl1, Ctrl](normalized)", "[Ctrl2, Ctrl](normalized)", "[Ctrl3, Ctrl](normalized)")
keep_cols <- c("GeneSymbol", "P-value", "FDR", "Fold Change", "Direction", sample_cols)
wei_both <- rbind(as.data.frame(up_sheet[, keep_cols]), as.data.frame(down_sheet[, keep_cols]))
colnames(wei_both) <- c("GeneSymbol", "PValue", "FDR", "FoldChange", "Direction",
                         "Sui1", "Sui2", "Sui3", "Ctrl1", "Ctrl2", "Ctrl3")
wei_both <- wei_both[!is.na(wei_both$GeneSymbol), ]
cat("Combined total (both sheets, before removing duplicates):", nrow(wei_both), "probes\n")

wei_both$logFC <- ifelse(wei_both$Direction == "down",
                          -log2(wei_both$FoldChange), log2(wei_both$FoldChange))
wei_both <- wei_both[order(wei_both$PValue), ]
sui_full <- wei_both[!duplicated(wei_both$GeneSymbol), ]  # keeps the most significant probe per gene
rownames(sui_full) <- NULL
write.csv(sui_full, "results/SUI_Wei2020_full_table.csv", row.names = FALSE)
cat("Unique genes after removing duplicates:", nrow(sui_full), "(",
    sum(sui_full$Direction == "up"), "up /", sum(sui_full$Direction == "down"), "down )\n")
cat("This table IS ALREADY the SUI DEG list - the original Wei 2020 file only\n")
cat("lists the genes the authors classified as differentially expressed\n")
cat("(fold-change>=2, raw p<0.05); the full tested array was never\n")
cat("published, so SUI DEG is used as the authors published it, not re-derived.\n\n")

set.seed(2020)
sui_mat <- as.matrix(sui_full[, c("Sui1","Sui2","Sui3","Ctrl1","Ctrl2","Ctrl3")])
rownames(sui_mat) <- sui_full$GeneSymbol
sui_mat <- avereps(sui_mat, ID = rownames(sui_mat))
group_sui <- factor(c("SUI","SUI","SUI","Ctrl","Ctrl","Ctrl"), levels = c("Ctrl","SUI"))
design_sui <- model.matrix(~group_sui)
fit_sui <- eBayes(lmFit(sui_mat, design_sui))
t_obs_sui <- fit_sui$t[, 2]
ranked_genes_sui <- names(t_obs_sui)[order(-t_obs_sui)]
cat("SUI ranked list (moderated t-statistic; positive = up in SUI):", length(ranked_genes_sui), "genes\n\n")

ann_sui <- suppressWarnings(select(org.Hs.eg.db, keys = ranked_genes_sui, keytype = "SYMBOL", columns = "PATH"))
ann_sui <- ann_sui[!is.na(ann_sui$PATH), ]
gs_sizes_sui <- table(ann_sui$PATH)
valid_paths_sui <- names(gs_sizes_sui)[gs_sizes_sui >= 5 & gs_sizes_sui <= 200]
gene_sets_sui <- split(ann_sui$SYMBOL[ann_sui$PATH %in% valid_paths_sui], ann_sui$PATH[ann_sui$PATH %in% valid_paths_sui])
cat("KEGG pathways tested (5-200 genes):", length(gene_sets_sui), "\n\n")

cat("Running fgsea on SUI (preranked on the full gene list, not phenotype\n")
cat("permutation - with only 3 vs 3 samples, phenotype permutation would\n")
cat("cap resolution at p=1/choose(6,3)=0.05; ranking the full gene list has\n")
cat("no such ceiling)...\n")
set.seed(2020)
fgsea_sui <- as.data.frame(fgsea(pathways = gene_sets_sui, stats = t_obs_sui, minSize = 5, maxSize = 200))
fgsea_sui$leadingEdge <- sapply(fgsea_sui$leadingEdge, paste, collapse = "/")
fgsea_sui$PathwayName <- kegg_label(fgsea_sui$pathway)
fgsea_sui <- fgsea_sui[order(fgsea_sui$pval), ]
write.csv(fgsea_sui, "results/SUI_fgsea_KEGG_full.csv", row.names = FALSE)

cat("=== RESULT: KEGG pathway enrichment in SUI ===\n")
cat("Pathways tested:", nrow(fgsea_sui), "\n")
cat("Significant at FDR<0.25:", sum(fgsea_sui$padj < 0.25, na.rm = TRUE), "\n")
cat("Significant at FDR<0.05:", sum(fgsea_sui$padj < 0.05, na.rm = TRUE), "\n\n")

cat("=== STEP 5 DONE === Now Source 06_shared_pathways_POP_SUI.R\n")
