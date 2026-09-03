# =============================================================================
# GSE208261 (POP, RNA-seq, ALL 24 samples) vs SUI (Wei 2020): DEG (DESeq2),
# GSEA on both, shared pathways - written FROM SCRATCH per explicit request,
# NOT reusing script 02's code, results, or cached files.
# =============================================================================
# THE REQUEST THIS SCRIPT ANSWERS: use all 12 "Control" samples (6 uterosacral
# ligament + 6 anterior vaginal wall combined) vs all 12 "POP" samples (all
# anterior vaginal wall) - i.e. do NOT drop the 6 ligament controls as script
# 02 did. The hypothesis behind this request: the original authors would not
# have included ligament samples "for nothing," so they must be intended for
# use, and since there is no POP-side ligament arm, the only way to use them
# is folded into the general Control group against all of POP.
#
# THIS SCRIPT DOES EXACTLY THAT, AS REQUESTED - but also runs a second,
# tissue-adjusted model and a direct diagnostic of the tissue effect, so the
# question "is this academically valid?" is answered with evidence from this
# actual dataset, not just a general statistical argument. See the STEP 1c
# diagnostic and the conclusion in this header's final paragraph (verified
# against the STEP 1 output below, not asserted a priori).
#
# ACADEMIC ASSESSMENT (write-up; read the STEP 1c numbers before trusting this
# paragraph blindly - it is written to be checked against them, not the other
# way around): combining ligament and vaginal-wall samples into one "Control"
# group and comparing it to an all-vaginal-wall "POP" group is a CONFOUNDED
# design for the specific 6 ligament samples: tissue type (ligament vs
# vaginal wall) is a real, likely LARGE source of expression variance on its
# own (different cell type composition, different anatomical structure)
# unrelated to POP status, and in this design tissue type is only variable
# among controls (all 12 POP samples are the single vaginal-wall tissue) -
# so any "Control vs POP" difference driven by the 6 ligament samples cannot
# be distinguished from a pure tissue-of-origin difference. This is different
# from, and worse than, a simple unbalanced design (like script 02's 6 vs 12,
# same tissue both sides) - there, only sample SIZE differs; here, the two
# groups differ in COMPOSITION in a way that is entangled with the biological
# question being asked. Adding tissue as a covariate (~tissue + group,
# included below) is the standard partial fix, but it is only a partial fix
# here: because no POP-side ligament samples exist, the "tissue effect" is
# estimated using controls only and then extrapolated to correct the POP
# comparison - a real assumption (that the ligament-vs-vaginal-wall shift is
# the same magnitude and direction in POP as in controls), not a directly
# tested fact.
#
# WHY THE AUTHORS LIKELY INCLUDED THE LIGAMENT SAMPLES ANYWAY (a plausible,
# non-speculative explanation, though their actual stated rationale is not in
# the metadata available here): POP surgical repair is typically performed
# on, and biopsies taken from, the vaginal wall - the site of the visible
# prolapse and the tissue actually operated on. Uterosacral ligament tissue
# is more readily available from CONTROL patients undergoing hysterectomy
# for unrelated benign indications (where the ligament is more routinely
# accessible/excised) than from POP patients, whose surgery is targeted at
# the prolapsed vaginal segment and may not always include ligament excision
# (procedure-dependent). Under this explanation, the ligament controls were
# not "wasted" - they likely serve a DIFFERENT comparison the original study
# may have made (e.g. characterizing baseline differences between pelvic
# support tissues in unaffected women), not a POP-vs-Control test, which is
# consistent with GSE208261's own metadata (`tissue: uterosacral ligaments`
# only ever appears in the Control arm - see data/GSE208261_sample_metadata.csv).
# This is this script's best inference, not a claim to certainty about the
# authors' intent.
# =============================================================================

suppressMessages({
  library(limma)
  library(edgeR)
  library(DESeq2)
  library(org.Hs.eg.db)
  library(ggplot2)
  library(pheatmap)
})
if (!requireNamespace("readxl", quietly = TRUE)) install.packages("readxl")
if (!requireNamespace("ggrepel", quietly = TRUE)) install.packages("ggrepel")
library(readxl)
library(ggrepel)

dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

# Full KEGG pathway ID -> name lookup, loaded from the shared, audited table
# data/kegg_pathway_names.csv (223 IDs covering every pathway used anywhere
# in this project) instead of a per-script partial dict - one source of
# truth so every figure across every script names the same ID the same way.
kegg_tab <- read.csv("data/kegg_pathway_names.csv", colClasses = c("character", "character"))
kegg_names <- setNames(kegg_tab$Name, kegg_tab$PATH5)
kegg_label <- function(id) {
  nm <- kegg_names[id]
  ifelse(is.na(nm), paste0("KEGG ", id), paste0("KEGG ", id, " - ", nm))
}
calc_es <- function(hit_idx, scores_abs, N) {
  Nh <- length(hit_idx); Nm <- N - Nh
  if (Nm <= 0 || Nh == 0) return(NA)
  sum_hit <- sum(scores_abs[hit_idx])
  step <- rep(-1 / Nm, N); step[hit_idx] <- scores_abs[hit_idx] / sum_hit
  running <- cumsum(step); running[which.max(abs(running))]
}


## =============================================================================
## STEP 1: Load data, verify sample metadata from scratch, build the 12-vs-12
##         design as requested (all 24 samples, ligament controls INCLUDED)
## =============================================================================
cat("\n================ STEP 1: Data and design (ALL 24 GSE208261 samples) ================\n\n")

meta <- read.csv("data/GSE208261_sample_metadata.csv")
cat("Full sample table (verify this against the GEO metadata before trusting it):\n")
print(meta)
cat("\nTissue x Group cross-tab:\n")
print(table(meta$Tissue, meta$Treatment))
cat("\n")

## --- STEP 1a: load counts, map Entrez -> symbol, collapse duplicates -------
counts_raw <- read.delim("data/GSE208261_raw_counts.tsv", row.names = 1, check.names = FALSE)
counts_all <- as.matrix(counts_raw[, meta$GSM])
stopifnot(identical(colnames(counts_all), meta$GSM))
cat("Raw count matrix: all 24 samples,", nrow(counts_all), "Entrez genes\n")

ann_id <- suppressWarnings(select(org.Hs.eg.db, keys = rownames(counts_all), keytype = "ENTREZID", columns = "SYMBOL"))
ann_id <- ann_id[!is.na(ann_id$SYMBOL) & !duplicated(ann_id$ENTREZID), ]
counts_all <- counts_all[rownames(counts_all) %in% ann_id$ENTREZID, ]
rownames(counts_all) <- ann_id$SYMBOL[match(rownames(counts_all), ann_id$ENTREZID)]
counts_all <- rowsum(counts_all, group = rownames(counts_all))
cat("After Entrez->symbol mapping and collapsing duplicates:", nrow(counts_all), "genes\n\n")

## --- STEP 1b: the 12-vs-12 design, exactly as requested --------------------
group_full <- factor(meta$Treatment, levels = c("Control", "POP"))  # 12 vs 12, tissue mixed within Control
tissue_full <- factor(meta$Tissue)
cat("Design used below: 12 Control (6 uterosacral ligament + 6 anterior\n")
cat("vaginal wall COMBINED) vs 12 POP (12 anterior vaginal wall). Control:",
    sum(group_full == "Control"), "| POP:", sum(group_full == "POP"), "\n\n")

## --- STEP 1c: DIAGNOSTIC - how large is the tissue effect on its own? -----
## Direct empirical check of the confound, using ONLY the 12 Control samples
## (where both tissues exist) - ligament vs vaginal wall, holding disease
## status (Control) constant. This measures the size of the very effect that
## would be folded into "Control vs POP" below for the 6 ligament samples.
cat("--- STEP 1c: DIAGNOSTIC - tissue effect size (ligament vs vaginal wall,\n")
cat("    CONTROLS ONLY, disease status held constant) ---\n\n")
ctrl_idx <- meta$Treatment == "Control"
dge_diag <- DGEList(counts = counts_all[, ctrl_idx], group = tissue_full[ctrl_idx])
keep_diag <- filterByExpr(dge_diag, model.matrix(~tissue_full[ctrl_idx]))
dge_diag <- dge_diag[keep_diag, , keep.lib.sizes = FALSE]
dge_diag <- calcNormFactors(dge_diag, method = "TMM")
design_diag <- model.matrix(~tissue_full[ctrl_idx])
voom_diag <- voom(dge_diag, design_diag)
fit_diag <- eBayes(lmFit(voom_diag, design_diag))
tt_diag <- topTable(fit_diag, coef = 2, number = Inf)
n_tissue_deg <- sum(tt_diag$adj.P.Val < 0.05 & abs(tt_diag$logFC) > 1)
cat("Genes kept:", nrow(dge_diag), "\n")
cat("Genes differing by TISSUE ALONE within controls (|log2FC|>1, FDR<0.05):",
    n_tissue_deg, "of", nrow(tt_diag), "(",
    round(100 * n_tissue_deg / nrow(tt_diag), 1), "% )\n")
cat("For comparison, script 02's tissue-matched POP-vs-Control DEG count\n")
cat("(6 vs 12, both vaginal wall) was 0-2 genes out of ~24,000-25,000 tested.\n")
cat("CONCLUSION OF THIS DIAGNOSTIC (numeric, not asserted): if the count above\n")
cat("is much larger than 0-2, tissue-of-origin alone moves far more genes than\n")
cat("POP status does in the tissue-matched comparison - meaning a naive\n")
cat("12-vs-12 comparison below risks being dominated by the tissue confound,\n")
cat("not by POP biology. Read the STEP 2 result in that light.\n\n")


## =============================================================================
## STEP 2: DEG, POP (all vaginal wall) vs Control (ligament+vaginal wall) -
##          the comparison exactly as requested, DESeq2 primary
## =============================================================================
cat("================ STEP 2: DEG for POP vs Control (12 vs 12, AS REQUESTED) ================\n\n")

counts_int <- counts_all; storage.mode(counts_int) <- "integer"

## --- Model A: naive, ~group only (exactly what was requested) -------------
coldata_naive <- data.frame(row.names = meta$GSM, group = group_full)
dds_naive <- DESeqDataSetFromMatrix(countData = counts_int, colData = coldata_naive, design = ~group)
dds_naive <- dds_naive[rowSums(counts(dds_naive) >= 10) >= 12, ]
dds_naive <- DESeq(dds_naive, quiet = TRUE)
res_naive <- as.data.frame(results(dds_naive, contrast = c("group", "POP", "Control"), alpha = 0.05))
res_naive$Gene <- rownames(res_naive)
res_naive <- res_naive[order(res_naive$pvalue), ]
pop_full <- res_naive[, c("Gene", "log2FoldChange", "baseMean", "stat", "pvalue", "padj")]
colnames(pop_full) <- c("Gene", "logFC", "baseMean", "stat", "P.Value", "adj.P.Val")
write.csv(pop_full, "results/16_POP_GSE208261_FULL24_naive_DESeq2_full.csv", row.names = FALSE)
pop_deg_naive <- subset(pop_full, !is.na(adj.P.Val) & adj.P.Val < 0.05 & abs(logFC) > 1)
pop_deg_naive <- pop_deg_naive[order(pop_deg_naive$adj.P.Val), ]
write.csv(pop_deg_naive, "results/16_POP_GSE208261_FULL24_naive_DEG_logFC1_FDR05.csv", row.names = FALSE)

cat("=== Model A: naive ~group (12 mixed-tissue Control vs 12 vaginal-wall POP) ===\n")
cat("Genes tested:", nrow(pop_full), "\n")
cat("DEG (|log2FC|>1, FDR<0.05):", nrow(pop_deg_naive), "(",
    sum(pop_deg_naive$logFC > 0), "up /", sum(pop_deg_naive$logFC < 0), "down )\n")
cat("Minimum FDR reached:", signif(min(pop_full$adj.P.Val, na.rm = TRUE), 3), "\n\n")

## --- Model B: tissue-adjusted, ~tissue + group (partial correction) -------
coldata_adj <- data.frame(row.names = meta$GSM, group = group_full, tissue = tissue_full)
dds_adj <- DESeqDataSetFromMatrix(countData = counts_int, colData = coldata_adj, design = ~tissue + group)
dds_adj <- dds_adj[rowSums(counts(dds_adj) >= 10) >= 12, ]
dds_adj <- DESeq(dds_adj, quiet = TRUE)
res_adj <- as.data.frame(results(dds_adj, name = "group_POP_vs_Control", alpha = 0.05))
res_adj$Gene <- rownames(res_adj)
res_adj <- res_adj[order(res_adj$pvalue), ]
pop_full_adj <- res_adj[, c("Gene", "log2FoldChange", "baseMean", "stat", "pvalue", "padj")]
colnames(pop_full_adj) <- c("Gene", "logFC", "baseMean", "stat", "P.Value", "adj.P.Val")
write.csv(pop_full_adj, "results/16_POP_GSE208261_FULL24_tissueAdj_DESeq2_full.csv", row.names = FALSE)
pop_deg_adj <- subset(pop_full_adj, !is.na(adj.P.Val) & adj.P.Val < 0.05 & abs(logFC) > 1)
write.csv(pop_deg_adj, "results/16_POP_GSE208261_FULL24_tissueAdj_DEG_logFC1_FDR05.csv", row.names = FALSE)

cat("=== Model B: tissue-adjusted ~tissue + group (partial correction, see header caveat) ===\n")
cat("DEG (|log2FC|>1, FDR<0.05):", nrow(pop_deg_adj), "\n")
cat("Minimum FDR reached:", signif(min(pop_full_adj$adj.P.Val, na.rm = TRUE), 3), "\n\n")

cat("=== COMPARISON: naive vs tissue-adjusted vs script 02's tissue-matched-only result ===\n")
cat("Naive (~group, 12v12, mixed tissue):        ", nrow(pop_deg_naive), "DEG\n")
cat("Tissue-adjusted (~tissue+group, 12v12):     ", nrow(pop_deg_adj), "DEG\n")
cat("Tissue-matched only (script 02, 6v12, same tissue both sides): 2 DEG (LOC105375520, COMP)\n")
cat("If naive >> tissue-matched-only, that is direct evidence the naive 12v12\n")
cat("comparison is picking up tissue composition, not (only) POP biology.\n\n")

## --- edgeR + voom + limma (secondary DEG check + GSEA ranking source) -----
dge <- DGEList(counts = counts_all, group = group_full)
design_naive_lm <- model.matrix(~group_full)
keep_expr <- filterByExpr(dge, design_naive_lm)
dge <- dge[keep_expr, , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge, method = "TMM")
voom_fit <- voom(dge, design_naive_lm)
fit_pop <- eBayes(lmFit(voom_fit, design_naive_lm))
pop_full_limma <- topTable(fit_pop, coef = "group_fullPOP", number = Inf, sort.by = "P")
pop_full_limma$Gene <- rownames(pop_full_limma)
pop_full_limma <- pop_full_limma[, c("Gene", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val")]
write.csv(pop_full_limma, "results/16_POP_GSE208261_FULL24_naive_voom_limma_full.csv", row.names = FALSE)
pop_deg_limma <- subset(pop_full_limma, adj.P.Val < 0.05 & abs(logFC) > 1)
cat("edgeR+voom+limma (naive ~group, secondary check): ", nrow(pop_deg_limma), "DEG\n\n")


## =============================================================================
## STEP 3: GSEA for POP (GSE208261, naive 12-vs-12 as requested) - CLASSIC
## =============================================================================
cat("================ STEP 3: GSEA for POP (GSE208261, 12v12) - CLASSIC method ================\n\n")
cat("Method: classic GSEA, phenotype permutation (shuffling the 24 group\n")
cat("labels, unpaired, balanced 12 vs 12). choose(24,12)=2,704,156 possible\n")
cat("relabelings - excellent resolution.\n")
cat("IMPORTANT INTERPRETIVE CAVEAT (see header): the permutation TEST is valid\n")
cat("for its literal null (no difference between these labels vs random\n")
cat("relabelings of the same 24 samples) - but because tissue composition\n")
cat("differs between the labeled groups (see STEP 1c), a significant result\n")
cat("cannot be attributed to POP status alone with confidence. Ranking =\n")
cat("moderated t-statistic from the naive edgeR+voom+limma fit above.\n\n")

set.seed(208261)
ranked_full <- fit_pop$t[, "group_fullPOP"]
ord2 <- order(-ranked_full)
ranked_genes_pop <- names(ranked_full)[ord2]
ranked_scores_pop <- ranked_full[ord2]
N_pop <- length(ranked_genes_pop)

ann_pop <- suppressWarnings(select(org.Hs.eg.db, keys = ranked_genes_pop, keytype = "SYMBOL", columns = "PATH"))
ann_pop <- ann_pop[!is.na(ann_pop$PATH), ]
gs_sizes_pop <- table(ann_pop$PATH)
valid_paths_pop <- names(gs_sizes_pop)[gs_sizes_pop >= 5 & gs_sizes_pop <= 200]
gene_sets_pop <- split(ann_pop$SYMBOL[ann_pop$PATH %in% valid_paths_pop], ann_pop$PATH[ann_pop$PATH %in% valid_paths_pop])
cat("KEGG pathways tested (5-200 members):", length(gene_sets_pop), "\n")

hit_idx_pop <- lapply(gene_sets_pop, function(g) which(ranked_genes_pop %in% g))
hit_idx_pop <- hit_idx_pop[sapply(hit_idx_pop, length) >= 3]
cat("Pathways with >=3 genes in the ranked list:", length(hit_idx_pop), "\n")
es_obs_pop <- sapply(hit_idx_pop, calc_es, scores_abs = abs(ranked_scores_pop), N = N_pop)

n_perm <- 500
cat("Running", n_perm, "phenotype (group-label) permutations...\n")
t0 <- Sys.time()
gene_sets_syms_pop <- lapply(hit_idx_pop, function(idx) ranked_genes_pop[idx])
perm_es_pop <- matrix(NA_real_, nrow = n_perm, ncol = length(hit_idx_pop))
for (i in seq_len(n_perm)) {
  perm_group <- sample(group_full)
  perm_design <- model.matrix(~perm_group)
  perm_fit <- eBayes(lmFit(voom_fit, perm_design))
  t_perm <- perm_fit$t[, 2]
  rank_of_gene <- rank(-t_perm, ties.method = "first")
  scores_abs_sorted <- sort(abs(t_perm), decreasing = TRUE)
  for (j in seq_along(gene_sets_syms_pop)) {
    hidx <- rank_of_gene[gene_sets_syms_pop[[j]]]
    if (length(hidx) >= 3) perm_es_pop[i, j] <- calc_es(hidx, scores_abs_sorted, N_pop)
  }
}
cat("Done in", round(difftime(Sys.time(), t0, units = "secs"), 1), "seconds\n\n")

pval_pop <- numeric(length(hit_idx_pop)); nes_pop <- numeric(length(hit_idx_pop))
for (j in seq_along(hit_idx_pop)) {
  pe <- perm_es_pop[, j]; pe <- pe[!is.na(pe)]
  if (es_obs_pop[j] >= 0) {
    pval_pop[j] <- (sum(pe >= es_obs_pop[j]) + 1) / (length(pe) + 1)
    base <- mean(pe[pe >= 0]); if (is.nan(base) || base == 0) base <- NA
  } else {
    pval_pop[j] <- (sum(pe <= es_obs_pop[j]) + 1) / (length(pe) + 1)
    base <- mean(abs(pe[pe < 0])); if (is.nan(base) || base == 0) base <- NA
  }
  nes_pop[j] <- es_obs_pop[j] / base
}

leading_edge_pop <- sapply(hit_idx_pop, function(idx) paste(ranked_genes_pop[idx], collapse = "/"))
gsea_pop <- data.frame(PATH = names(hit_idx_pop), Nh = sapply(hit_idx_pop, length),
                        ES = es_obs_pop, NES = nes_pop, pvalue = pval_pop, leadingEdge = leading_edge_pop)
gsea_pop$p.adjust <- p.adjust(gsea_pop$pvalue, "BH")
gsea_pop$PathwayName <- kegg_label(gsea_pop$PATH)
gsea_pop <- gsea_pop[order(gsea_pop$pvalue), ]
write.csv(gsea_pop, "results/17_GSEA_classic_POP_GSE208261_FULL24_KEGG.csv", row.names = FALSE)

cat("=== GSEA classic, POP GSE208261 FULL 12v12 (KEGG) ===\n")
cat("Pathways significant at FDR<0.25:", sum(gsea_pop$p.adjust < 0.25, na.rm = TRUE), "of", nrow(gsea_pop), "\n")
cat("Pathways significant at FDR<0.05:", sum(gsea_pop$p.adjust < 0.05, na.rm = TRUE), "\n\n")
print(head(gsea_pop[, c("PathwayName","Nh","NES","pvalue","p.adjust")], 10))
cat("\n")


## =============================================================================
## STEP 4: SUI (Wei 2020) - freshly re-read and re-verified from scratch
## =============================================================================
cat("================ STEP 4: DEG and GSEA for SUI (Wei 2020) - fresh re-derivation ================\n\n")
cat("Per explicit request, this is written fresh, not copy-pasted from any\n")
cat("earlier script - re-verifying each step rather than assuming it is\n")
cat("still correct.\n\n")

wei_xls <- "data/Wei2020_TableS2_mRNA.xls"
sheet_names <- excel_sheets(wei_xls)
cat("Sheets found in the Excel file:", paste(sheet_names, collapse = ", "), "\n")
stopifnot(all(c("up_Sui_vs_Ctrl", "down_Sui_vs_Ctrl") %in% sheet_names))

up_sheet   <- read_excel(wei_xls, sheet = "up_Sui_vs_Ctrl",   skip = 17)
down_sheet <- read_excel(wei_xls, sheet = "down_Sui_vs_Ctrl", skip = 17)
cat("'up_Sui_vs_Ctrl':", nrow(up_sheet), "rows,", ncol(up_sheet), "columns\n")
cat("'down_Sui_vs_Ctrl':", nrow(down_sheet), "rows,", ncol(down_sheet), "columns\n")

# Re-verify the Fold Change sign convention independently, rather than
# assuming the earlier finding still holds.
cat("\nRe-checking Fold Change sign convention in each sheet:\n")
cat("  up_Sui_vs_Ctrl Fold Change range:", round(min(up_sheet$`Fold Change`), 2), "to",
    round(max(up_sheet$`Fold Change`), 2), "\n")
cat("  down_Sui_vs_Ctrl Fold Change range:", round(min(down_sheet$`Fold Change`), 2), "to",
    round(max(down_sheet$`Fold Change`), 2), "\n")
if (min(down_sheet$`Fold Change`) >= 1) {
  cat("  CONFIRMED: down-sheet Fold Change is an unsigned MAGNITUDE (always >=1)\n")
  cat("  - sign must come from the sheet/Direction, not the value. Same finding\n")
  cat("  as before, independently re-verified here.\n\n")
}

up_sheet$Direction <- "up"; down_sheet$Direction <- "down"
sample_cols <- c("[Sui1, Sui](normalized)", "[Sui2, Sui](normalized)", "[Sui3, Sui](normalized)",
                  "[Ctrl1, Ctrl](normalized)", "[Ctrl2, Ctrl](normalized)", "[Ctrl3, Ctrl](normalized)")
keep_cols <- c("GeneSymbol", "P-value", "FDR", "Fold Change", "Direction", sample_cols)
wei_both <- rbind(as.data.frame(up_sheet[, keep_cols]), as.data.frame(down_sheet[, keep_cols]))
colnames(wei_both) <- c("GeneSymbol", "PValue", "FDR", "FoldChange", "Direction",
                         "Sui1", "Sui2", "Sui3", "Ctrl1", "Ctrl2", "Ctrl3")
wei_both <- wei_both[!is.na(wei_both$GeneSymbol), ]
cat("Combined (both sheets, before dedup):", nrow(wei_both), "probes -",
    "expect 4615+2487=7102:", nrow(wei_both) == 7102, "\n")

wei_both$logFC <- ifelse(wei_both$Direction == "down",
                          -log2(wei_both$FoldChange), log2(wei_both$FoldChange))
wei_both <- wei_both[order(wei_both$PValue), ]
sui_full <- wei_both[!duplicated(wei_both$GeneSymbol), ]
rownames(sui_full) <- NULL
sui_deg <- sui_full
write.csv(sui_deg, "results/18_SUI_Wei2020_full_FRESH.csv", row.names = FALSE)
cat("Unique genes after dedup:", nrow(sui_deg), "(",
    sum(sui_deg$Direction == "up"), "up /", sum(sui_deg$Direction == "down"), "down )\n\n")

set.seed(2020)
sui_mat <- as.matrix(sui_full[, c("Sui1","Sui2","Sui3","Ctrl1","Ctrl2","Ctrl3")])
rownames(sui_mat) <- sui_full$GeneSymbol
sui_mat <- avereps(sui_mat, ID = rownames(sui_mat))
group_sui <- factor(c("SUI","SUI","SUI","Ctrl","Ctrl","Ctrl"), levels = c("Ctrl","SUI"))
design_sui <- model.matrix(~group_sui)
fit_sui <- eBayes(lmFit(sui_mat, design_sui))
t_obs_sui <- fit_sui$t[, 2]
ord_s <- order(-t_obs_sui)
ranked_genes_sui <- names(t_obs_sui)[ord_s]
ranked_scores_sui <- t_obs_sui[ord_s]
N_sui <- length(ranked_genes_sui)
cat("SUI ranked list (moderated t, positive = up in SUI):", N_sui, "genes\n\n")

cat("Method: preranked GSEA for SUI (gene-set-label permutation) - classic\n")
cat("phenotype permutation is not viable with only 3 vs 3 samples\n")
cat("(choose(6,3)=20 relabelings caps resolution at p=1/20=0.05).\n\n")

ann_sui <- suppressWarnings(select(org.Hs.eg.db, keys = ranked_genes_sui, keytype = "SYMBOL", columns = "PATH"))
ann_sui <- ann_sui[!is.na(ann_sui$PATH), ]
gs_sizes_sui <- table(ann_sui$PATH)
valid_paths_sui <- names(gs_sizes_sui)[gs_sizes_sui >= 5 & gs_sizes_sui <= 200]
gene_sets_sui <- split(ann_sui$SYMBOL[ann_sui$PATH %in% valid_paths_sui], ann_sui$PATH[ann_sui$PATH %in% valid_paths_sui])
hit_idx_sui <- lapply(gene_sets_sui, function(g) which(ranked_genes_sui %in% g))
hit_idx_sui <- hit_idx_sui[sapply(hit_idx_sui, length) >= 3]
es_obs_sui <- sapply(hit_idx_sui, calc_es, scores_abs = abs(ranked_scores_sui), N = N_sui)

n_perm2 <- 1000
cat("Running", n_perm2, "gene-set-label permutations for SUI...\n")
t0 <- Sys.time()
scores_abs_sui <- abs(ranked_scores_sui)
perm_es_sui <- matrix(NA_real_, nrow = n_perm2, ncol = length(hit_idx_sui))
for (i in seq_len(n_perm2)) {
  for (j in seq_along(hit_idx_sui)) {
    hidx <- sample.int(N_sui, length(hit_idx_sui[[j]]))
    perm_es_sui[i, j] <- calc_es(hidx, scores_abs_sui, N_sui)
  }
}
cat("Done in", round(difftime(Sys.time(), t0, units = "secs"), 1), "seconds\n\n")

pval_sui <- numeric(length(hit_idx_sui)); nes_sui <- numeric(length(hit_idx_sui))
for (j in seq_along(hit_idx_sui)) {
  pe <- perm_es_sui[, j]; pe <- pe[!is.na(pe)]
  if (es_obs_sui[j] >= 0) {
    pval_sui[j] <- (sum(pe >= es_obs_sui[j]) + 1) / (length(pe) + 1)
    base <- mean(pe[pe >= 0]); if (is.nan(base) || base == 0) base <- NA
  } else {
    pval_sui[j] <- (sum(pe <= es_obs_sui[j]) + 1) / (length(pe) + 1)
    base <- mean(abs(pe[pe < 0])); if (is.nan(base) || base == 0) base <- NA
  }
  nes_sui[j] <- es_obs_sui[j] / base
}

leading_edge_sui <- sapply(hit_idx_sui, function(idx) paste(ranked_genes_sui[idx], collapse = "/"))
gsea_sui <- data.frame(PATH = names(hit_idx_sui), Nh = sapply(hit_idx_sui, length),
                        ES = es_obs_sui, NES = nes_sui, pvalue = pval_sui, leadingEdge = leading_edge_sui)
gsea_sui$p.adjust <- p.adjust(gsea_sui$pvalue, "BH")
gsea_sui$PathwayName <- kegg_label(gsea_sui$PATH)
gsea_sui <- gsea_sui[order(gsea_sui$pvalue), ]
write.csv(gsea_sui, "results/18_GSEA_preranked_SUI_KEGG_FRESH.csv", row.names = FALSE)

cat("=== GSEA preranked, SUI (KEGG) ===\n")
cat("Pathways significant at FDR<0.25:", sum(gsea_sui$p.adjust < 0.25, na.rm = TRUE), "of", nrow(gsea_sui), "\n")
cat("Pathways significant at FDR<0.05:", sum(gsea_sui$p.adjust < 0.05, na.rm = TRUE), "\n\n")


## =============================================================================
## STEP 5: Shared pathways (POP GSE208261 FULL24 x SUI) + gene direction table
## =============================================================================
cat("================ STEP 5: Shared pathways - the professor's original request ================\n\n")
cat("分别对POP数据集和SUI数据集运行GSEA, 找出两者共同显著富集的通路\n")
cat("(Run GSEA separately on the POP and SUI datasets, find the pathways\n")
cat("significantly enriched in common between them.)\n\n")

report_shared5 <- function(fdr_cut) {
  pop_sig <- subset(gsea_pop, p.adjust < fdr_cut)
  sui_sig <- subset(gsea_sui, p.adjust < fdr_cut)
  shared_ids <- intersect(pop_sig$PATH, sui_sig$PATH)
  cat("--- FDR <", fdr_cut, ": POP(GSE208261 FULL24) significant =", nrow(pop_sig),
      "| SUI significant =", nrow(sui_sig), "| SHARED =", length(shared_ids), "---\n")
  if (length(shared_ids) == 0) return(data.frame())
  out <- merge(pop_sig[pop_sig$PATH %in% shared_ids, c("PATH","PathwayName","Nh","NES","p.adjust")],
               sui_sig[sui_sig$PATH %in% shared_ids, c("PATH","Nh","NES","p.adjust")],
               by = "PATH", suffixes = c("_POP", "_SUI"))
  out$Same_direction <- sign(out$NES_POP) == sign(out$NES_SUI)
  out <- out[order(out$p.adjust_POP), ]
  print(out[, c("PathwayName","Nh_POP","NES_POP","p.adjust_POP","Nh_SUI","NES_SUI","p.adjust_SUI","Same_direction")])
  n_na <- sum(is.na(out$Same_direction))
  cat("Same direction in both diseases:", sum(out$Same_direction, na.rm = TRUE), "of",
      sum(!is.na(out$Same_direction)), "pathways with a comparable NES on both sides")
  if (n_na > 0) cat(" (", n_na, "excluded - NES undefined on at least one side)")
  cat("\n\n")
  out
}

shared5_025 <- report_shared5(0.25)
shared5_005 <- report_shared5(0.05)
if (nrow(shared5_025) > 0) write.csv(shared5_025, "results/19_shared_pathways_GSE208261FULL24xSUI_FDR025.csv", row.names = FALSE)
if (nrow(shared5_005) > 0) write.csv(shared5_005, "results/19_shared_pathways_GSE208261FULL24xSUI_FDR005.csv", row.names = FALSE)

shared5_for_table <- if (nrow(shared5_005) > 0) shared5_005 else shared5_025
cat("\n--- Gene-level direction table for the shared pathways ---\n")
if (nrow(shared5_for_table) == 0) {
  cat("No shared pathway at either threshold.\n\n")
  gene_dir_table5 <- data.frame()
} else {
  pop_lookup5 <- setNames(pop_full$logFC, pop_full$Gene)
  sui_lookup5 <- setNames(sui_full$logFC, sui_full$GeneSymbol)
  pop5_tested <- ranked_genes_pop
  sui_tested5 <- ranked_genes_sui

  rows5 <- lapply(seq_len(nrow(shared5_for_table)), function(i) {
    pid <- shared5_for_table$PATH[i]
    pname <- shared5_for_table$PathwayName[i]
    full_members <- suppressWarnings(select(org.Hs.eg.db, keys = pid, keytype = "PATH", columns = "SYMBOL")$SYMBOL)
    members <- intersect(unique(full_members), union(pop5_tested, sui_tested5))
    data.frame(
      Pathway = pname, KEGG_ID = pid, Gene = members,
      logFC_POP = unname(pop_lookup5[members]), logFC_SUI = unname(sui_lookup5[members]),
      stringsAsFactors = FALSE)
  })
  gene_dir_table5 <- do.call(rbind, rows5)
  gene_dir_table5$Direction_POP <- ifelse(is.na(gene_dir_table5$logFC_POP), "not tested",
                                           ifelse(gene_dir_table5$logFC_POP > 0, "up", "down"))
  gene_dir_table5$Direction_SUI <- ifelse(is.na(gene_dir_table5$logFC_SUI), "not tested",
                                           ifelse(gene_dir_table5$logFC_SUI > 0, "up", "down"))
  gene_dir_table5$Concordant <- with(gene_dir_table5,
                                      Direction_POP != "not tested" & Direction_SUI != "not tested" &
                                        Direction_POP == Direction_SUI)
  gene_dir_table5 <- gene_dir_table5[order(gene_dir_table5$Pathway, -gene_dir_table5$Concordant), ]
  write.csv(gene_dir_table5, "results/19_shared_pathways_gene_direction_table.csv", row.names = FALSE)

  both_tested5 <- subset(gene_dir_table5, Direction_POP != "not tested" & Direction_SUI != "not tested")
  cat("Genes in shared pathways tested in both diseases:", nrow(both_tested5),
      "| concordant direction:", sum(both_tested5$Concordant),
      "(", round(100 * mean(both_tested5$Concordant), 1), "% )\n\n")
  print(head(gene_dir_table5, 20))
}
cat("\n")


## =============================================================================
## STEP 6: GRAPHICS
## =============================================================================
cat("================ STEP 6: Graphics ================\n\n")

make_volcano <- function(df, logfc_col, p_col, label_col, title, fc_cut = 1, p_cut = 0.05) {
  df <- df[!is.na(df[[logfc_col]]) & !is.na(df[[p_col]]), ]
  df$negLog10P <- -log10(df[[p_col]])
  df$sig <- "NS"
  df$sig[df[[logfc_col]] > fc_cut & df[[p_col]] < p_cut] <- "Up"
  df$sig[df[[logfc_col]] < -fc_cut & df[[p_col]] < p_cut] <- "Down"
  df$sig <- factor(df$sig, levels = c("Down", "NS", "Up"))
  top_lab <- df[order(df[[p_col]]), ][1:15, ]
  ggplot(df, aes(x = .data[[logfc_col]], y = negLog10P, color = sig)) +
    geom_point(alpha = 0.6, size = 1.3) +
    scale_color_manual(values = c(Down = "#2166AC", NS = "grey75", Up = "#B2182B")) +
    geom_vline(xintercept = c(-fc_cut, fc_cut), linetype = "dashed", color = "grey40") +
    geom_hline(yintercept = -log10(p_cut), linetype = "dashed", color = "grey40") +
    geom_text_repel(data = top_lab, aes(label = .data[[label_col]]), size = 3, color = "black",
                     max.overlaps = 20, segment.size = 0.2) +
    labs(title = title, x = "log2(Fold Change)",
         y = expression(-log[10](italic(p)~value)), color = NULL) +
    theme_bw() + theme(legend.position = "top", plot.title = element_text(face = "bold"))
}

pop_full$Gene_label <- pop_full$Gene
deg_note <- if (nrow(pop_deg_naive) == 0) "none reach FDR<0.05" else paste0(nrow(pop_deg_naive), " genes reach FDR<0.05")
v_pop <- make_volcano(pop_full, "logFC", "P.Value", "Gene_label",
                       paste0("Volcano plot - POP (GSE208261, FULL 12v12, DESeq2, naive)\n",
                              "12 Control (ligament+vaginal wall) vs 12 POP (vaginal wall)\n",
                              "(", deg_note, " - tissue-confound caveat, see README)"))
ggsave("figures/24_volcano_POP_GSE208261_FULL24.png", v_pop, width = 9.5, height = 6.5, dpi = 300)
cat("Saved: figures/24_volcano_POP_GSE208261_FULL24.png\n\n")

# PCA showing tissue vs disease structure directly (the core diagnostic,
# visualized) - all 24 samples, colored by tissue, shaped by group.
logcpm_all <- voom_fit$E
pca <- prcomp(t(logcpm_all), scale. = TRUE)
pca_df <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2],
                      Tissue = meta$Tissue, Group = meta$Treatment)
var_exp <- round(100 * summary(pca)$importance[2, 1:2], 1)
p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Tissue, shape = Group)) +
  geom_point(size = 4, alpha = 0.85) +
  labs(title = "PCA, all 24 GSE208261 samples - tissue vs disease structure",
       x = paste0("PC1 (", var_exp[1], "%)"), y = paste0("PC2 (", var_exp[2], "%)")) +
  theme_bw() + theme(legend.position = "right")
ggsave("figures/25_PCA_tissue_vs_disease_GSE208261.png", p_pca, width = 8.5, height = 6, dpi = 300)
cat("Saved: figures/25_PCA_tissue_vs_disease_GSE208261.png (the key diagnostic figure)\n\n")

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
p_bar_pop <- make_gsea_barplot(gsea_pop, "GSEA (classic) - top KEGG pathways in POP (GSE208261, FULL 12v12)")
ggsave("figures/26_GSEA_barplot_POP_GSE208261_FULL24.png", p_bar_pop, width = 10, height = 6, dpi = 300)
cat("Saved: figures/26_GSEA_barplot_POP_GSE208261_FULL24.png\n\n")

compare_pool5 <- if (nrow(shared5_025) > 0) shared5_025 else data.frame()
if (nrow(compare_pool5) > 0) {
  p_shared5 <- ggplot(compare_pool5, aes(x = NES_POP, y = NES_SUI, label = PathwayName)) +
    geom_hline(yintercept = 0, color = "grey70") + geom_vline(xintercept = 0, color = "grey70") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
    geom_point(aes(color = Same_direction), size = 3) +
    scale_color_manual(values = c("TRUE" = "#2166AC", "FALSE" = "#B2182B"),
                        labels = c("TRUE" = "Same direction", "FALSE" = "Opposite direction")) +
    geom_text_repel(size = 3, max.overlaps = 30) +
    labs(title = "Shared pathways (FDR<0.25): NES in POP (GSE208261 FULL24) vs NES in SUI",
         x = "NES - POP GSE208261 FULL24 (classic GSEA)", y = "NES - SUI (preranked GSEA)", color = NULL) +
    theme_bw() + theme(legend.position = "top")
  ggsave("figures/27_shared_pathways_NES_comparison_GSE208261_FULL24.png", p_shared5, width = 9, height = 7, dpi = 300)
  cat("Saved: figures/27_shared_pathways_NES_comparison_GSE208261_FULL24.png\n\n")
} else {
  cat("SKIPPED: figures/27_shared_pathways_NES_comparison_GSE208261_FULL24.png (no shared pathway at FDR<0.25)\n\n")
}

cat("\n=== ALL GSE208261 FULL24 FIGURES SAVED TO figures/ ===\n")
cat("24_volcano_POP_GSE208261_FULL24.png\n25_PCA_tissue_vs_disease_GSE208261.png (key diagnostic)\n")
cat("26_GSEA_barplot_POP_GSE208261_FULL24.png\n27_shared_pathways_NES_comparison_GSE208261_FULL24.png (if any shared)\n")
cat("\n=== DONE (script 05) ===\n")
