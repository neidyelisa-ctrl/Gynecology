# =============================================================================
# GSE208261 (POP, RNA-seq), BALANCED 6-vs-6, vaginal wall only, vs SUI (Wei 2020)
# =============================================================================
# Follow-up to script 05's tissue-confound finding: script 02 already used
# only vaginal-wall samples (no tissue confound), but was UNBALANCED (6
# Control vs 12 POP). This script tests whether balancing to 6-vs-6 changes
# anything - a different, secondary concern from the tissue confound (group
# SIZE/variance-estimate stability, not tissue composition).
#
# WHICH 6 OF THE 12 POP SAMPLES? The "_D"/"_Y" suffix does not correspond to
# tissue for POP (both are anterior vaginal wall - see script 02's header) -
# so there is no principled reason to prefer POP_D or POP_Y over the other.
# Rather than pick one (which would be an arbitrary, unreported choice this
# project's stated values reject), THIS SCRIPT RUNS BOTH as independent
# checks: 6 Control_Y vs 6 POP_Y, and 6 Control_Y vs 6 POP_D. If conclusions
# agree regardless of which half of POP is used, that is itself useful
# evidence the choice doesn't matter; if they disagree, that is worth
# knowing too.
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

meta <- read.csv("data/GSE208261_sample_metadata.csv")
counts_raw <- read.delim("data/GSE208261_raw_counts.tsv", row.names = 1, check.names = FALSE)


## =============================================================================
## STEP 1: DEG, both balanced 6v6 subsets (DESeq2)
## =============================================================================
cat("\n================ STEP 1: DEG, balanced 6v6 (vaginal wall only) ================\n\n")

run_deg <- function(pop_titles, label) {
  ctrl <- meta[meta$Title %in% paste0("Control_Y", 1:6), ]
  pop  <- meta[meta$Title %in% pop_titles, ]
  sub_meta <- rbind(ctrl, pop)
  counts <- as.matrix(counts_raw[, sub_meta$GSM])
  ann_id <- suppressWarnings(select(org.Hs.eg.db, keys = rownames(counts), keytype = "ENTREZID", columns = "SYMBOL"))
  ann_id <- ann_id[!is.na(ann_id$SYMBOL) & !duplicated(ann_id$ENTREZID), ]
  counts <- counts[rownames(counts) %in% ann_id$ENTREZID, ]
  rownames(counts) <- ann_id$SYMBOL[match(rownames(counts), ann_id$ENTREZID)]
  counts <- rowsum(counts, group = rownames(counts))
  counts_int <- counts; storage.mode(counts_int) <- "integer"

  group <- factor(sub_meta$Treatment, levels = c("Control", "POP"))
  coldata <- data.frame(row.names = sub_meta$GSM, group = group)
  dds <- DESeqDataSetFromMatrix(countData = counts_int, colData = coldata, design = ~group)
  dds <- dds[rowSums(counts(dds) >= 10) >= 6, ]
  dds <- DESeq(dds, quiet = TRUE)
  res <- as.data.frame(results(dds, contrast = c("group", "POP", "Control"), alpha = 0.05))
  res$Gene <- rownames(res)
  res <- res[order(res$pvalue), ]
  full <- res[, c("Gene", "log2FoldChange", "baseMean", "stat", "pvalue", "padj")]
  colnames(full) <- c("Gene", "logFC", "baseMean", "stat", "P.Value", "adj.P.Val")
  deg <- subset(full, !is.na(adj.P.Val) & adj.P.Val < 0.05 & abs(logFC) > 1)
  cat("---", label, "---\n")
  cat("Genes tested:", nrow(full), "| DEG (|log2FC|>1, FDR<0.05):", nrow(deg),
      "(", sum(deg$logFC > 0), "up /", sum(deg$logFC < 0), "down ) | min FDR:",
      signif(min(full$adj.P.Val, na.rm = TRUE), 3), "\n")
  if (nrow(deg) > 0) print(head(deg[, c("Gene","logFC","adj.P.Val")], 15))
  cat("\n")
  list(full = full, deg = deg, counts = counts, meta = sub_meta, dds = dds)
}

res_Y <- run_deg(paste0("POP_Y", 1:6), "6 Control_Y vs 6 POP_Y")
res_D <- run_deg(paste0("POP_D", 1:6), "6 Control_Y vs 6 POP_D")

write.csv(res_Y$full, "results/20_POP_GSE208261_6v6_POPY_DESeq2_full.csv", row.names = FALSE)
write.csv(res_Y$deg,  "results/20_POP_GSE208261_6v6_POPY_DEG.csv", row.names = FALSE)
write.csv(res_D$full, "results/20_POP_GSE208261_6v6_POPD_DESeq2_full.csv", row.names = FALSE)
write.csv(res_D$deg,  "results/20_POP_GSE208261_6v6_POPD_DEG.csv", row.names = FALSE)

overlap <- intersect(res_Y$deg$Gene, res_D$deg$Gene)
cat("=== Consistency check: DEG overlap between the two independent 6v6 subsets ===\n")
cat("POP_Y DEG:", nrow(res_Y$deg), "| POP_D DEG:", nrow(res_D$deg),
    "| overlap:", length(overlap), "genes\n")
if (length(overlap) > 0) cat("Shared genes:", paste(overlap, collapse = ", "), "\n")
cat("\n")


## =============================================================================
## STEP 2: GSEA, classic, both 6v6 subsets
## =============================================================================
cat("================ STEP 2: GSEA classic, both 6v6 subsets ================\n\n")
cat("choose(12,6)=924 possible relabelings per subset - decent resolution,\n")
cat("similar to script 04's VFB design.\n\n")

run_gsea <- function(res_obj, label, seed) {
  set.seed(seed)
  group <- factor(res_obj$meta$Treatment, levels = c("Control", "POP"))
  dge <- DGEList(counts = res_obj$counts, group = group)
  design <- model.matrix(~group)
  keep_expr <- filterByExpr(dge, design)
  dge <- dge[keep_expr, , keep.lib.sizes = FALSE]
  dge <- calcNormFactors(dge, method = "TMM")
  voom_fit <- voom(dge, design)
  fit <- eBayes(lmFit(voom_fit, design))
  ranked_full <- fit$t[, 2]
  ord <- order(-ranked_full)
  ranked_genes <- names(ranked_full)[ord]
  ranked_scores <- ranked_full[ord]
  N <- length(ranked_genes)

  ann <- suppressWarnings(select(org.Hs.eg.db, keys = ranked_genes, keytype = "SYMBOL", columns = "PATH"))
  ann <- ann[!is.na(ann$PATH), ]
  gs_sizes <- table(ann$PATH)
  valid_paths <- names(gs_sizes)[gs_sizes >= 5 & gs_sizes <= 200]
  gene_sets <- split(ann$SYMBOL[ann$PATH %in% valid_paths], ann$PATH[ann$PATH %in% valid_paths])
  hit_idx <- lapply(gene_sets, function(g) which(ranked_genes %in% g))
  hit_idx <- hit_idx[sapply(hit_idx, length) >= 3]
  cat(label, "- KEGG pathways tested:", length(hit_idx), "\n")
  es_obs <- sapply(hit_idx, calc_es, scores_abs = abs(ranked_scores), N = N)

  n_perm <- 500
  t0 <- Sys.time()
  gene_sets_syms <- lapply(hit_idx, function(idx) ranked_genes[idx])
  perm_es <- matrix(NA_real_, nrow = n_perm, ncol = length(hit_idx))
  for (i in seq_len(n_perm)) {
    perm_group <- sample(group)
    perm_design <- model.matrix(~perm_group)
    perm_fit <- eBayes(lmFit(voom_fit, perm_design))
    t_perm <- perm_fit$t[, 2]
    rank_of_gene <- rank(-t_perm, ties.method = "first")
    scores_abs_sorted <- sort(abs(t_perm), decreasing = TRUE)
    for (j in seq_along(gene_sets_syms)) {
      hidx <- rank_of_gene[gene_sets_syms[[j]]]
      if (length(hidx) >= 3) perm_es[i, j] <- calc_es(hidx, scores_abs_sorted, N)
    }
  }
  cat(label, "- permutations done in", round(difftime(Sys.time(), t0, units = "secs"), 1), "seconds\n")

  pval <- numeric(length(hit_idx)); nes <- numeric(length(hit_idx))
  for (j in seq_along(hit_idx)) {
    pe <- perm_es[, j]; pe <- pe[!is.na(pe)]
    if (es_obs[j] >= 0) {
      pval[j] <- (sum(pe >= es_obs[j]) + 1) / (length(pe) + 1)
      base <- mean(pe[pe >= 0]); if (is.nan(base) || base == 0) base <- NA
    } else {
      pval[j] <- (sum(pe <= es_obs[j]) + 1) / (length(pe) + 1)
      base <- mean(abs(pe[pe < 0])); if (is.nan(base) || base == 0) base <- NA
    }
    nes[j] <- es_obs[j] / base
  }
  leading_edge <- sapply(hit_idx, function(idx) paste(ranked_genes[idx], collapse = "/"))
  gsea_res <- data.frame(PATH = names(hit_idx), Nh = sapply(hit_idx, length),
                          ES = es_obs, NES = nes, pvalue = pval, leadingEdge = leading_edge)
  gsea_res$p.adjust <- p.adjust(gsea_res$pvalue, "BH")
  gsea_res$PathwayName <- kegg_label(gsea_res$PATH)
  gsea_res <- gsea_res[order(gsea_res$pvalue), ]
  cat(label, "- significant at FDR<0.25:", sum(gsea_res$p.adjust < 0.25, na.rm = TRUE),
      "| FDR<0.05:", sum(gsea_res$p.adjust < 0.05, na.rm = TRUE), "\n\n")
  list(gsea = gsea_res, ranked_genes = ranked_genes, fit = fit)
}

gsea_Y <- run_gsea(res_Y, "POP_Y", 20261)
gsea_D <- run_gsea(res_D, "POP_D", 20262)

write.csv(gsea_Y$gsea, "results/21_GSEA_classic_POP_GSE208261_6v6_POPY_KEGG.csv", row.names = FALSE)
write.csv(gsea_D$gsea, "results/21_GSEA_classic_POP_GSE208261_6v6_POPD_KEGG.csv", row.names = FALSE)


## =============================================================================
## STEP 3: SUI (Wei 2020), fresh
## =============================================================================
cat("================ STEP 3: SUI (Wei 2020) ================\n\n")

wei_xls <- "data/Wei2020_TableS2_mRNA.xls"
up_sheet   <- read_excel(wei_xls, sheet = "up_Sui_vs_Ctrl",   skip = 17)
down_sheet <- read_excel(wei_xls, sheet = "down_Sui_vs_Ctrl", skip = 17)
up_sheet$Direction <- "up"; down_sheet$Direction <- "down"
sample_cols <- c("[Sui1, Sui](normalized)", "[Sui2, Sui](normalized)", "[Sui3, Sui](normalized)",
                  "[Ctrl1, Ctrl](normalized)", "[Ctrl2, Ctrl](normalized)", "[Ctrl3, Ctrl](normalized)")
keep_cols <- c("GeneSymbol", "P-value", "FDR", "Fold Change", "Direction", sample_cols)
wei_both <- rbind(as.data.frame(up_sheet[, keep_cols]), as.data.frame(down_sheet[, keep_cols]))
colnames(wei_both) <- c("GeneSymbol", "PValue", "FDR", "FoldChange", "Direction",
                         "Sui1", "Sui2", "Sui3", "Ctrl1", "Ctrl2", "Ctrl3")
wei_both <- wei_both[!is.na(wei_both$GeneSymbol), ]
wei_both$logFC <- ifelse(wei_both$Direction == "down",
                          -log2(wei_both$FoldChange), log2(wei_both$FoldChange))
wei_both <- wei_both[order(wei_both$PValue), ]
sui_full <- wei_both[!duplicated(wei_both$GeneSymbol), ]
rownames(sui_full) <- NULL

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

ann_sui <- suppressWarnings(select(org.Hs.eg.db, keys = ranked_genes_sui, keytype = "SYMBOL", columns = "PATH"))
ann_sui <- ann_sui[!is.na(ann_sui$PATH), ]
gs_sizes_sui <- table(ann_sui$PATH)
valid_paths_sui <- names(gs_sizes_sui)[gs_sizes_sui >= 5 & gs_sizes_sui <= 200]
gene_sets_sui <- split(ann_sui$SYMBOL[ann_sui$PATH %in% valid_paths_sui], ann_sui$PATH[ann_sui$PATH %in% valid_paths_sui])
hit_idx_sui <- lapply(gene_sets_sui, function(g) which(ranked_genes_sui %in% g))
hit_idx_sui <- hit_idx_sui[sapply(hit_idx_sui, length) >= 3]
es_obs_sui <- sapply(hit_idx_sui, calc_es, scores_abs = abs(ranked_scores_sui), N = N_sui)

n_perm2 <- 1000
scores_abs_sui <- abs(ranked_scores_sui)
perm_es_sui <- matrix(NA_real_, nrow = n_perm2, ncol = length(hit_idx_sui))
for (i in seq_len(n_perm2)) {
  for (j in seq_along(hit_idx_sui)) {
    hidx <- sample.int(N_sui, length(hit_idx_sui[[j]]))
    perm_es_sui[i, j] <- calc_es(hidx, scores_abs_sui, N_sui)
  }
}
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
cat("SUI GSEA significant at FDR<0.25:", sum(gsea_sui$p.adjust < 0.25, na.rm = TRUE),
    "| FDR<0.05:", sum(gsea_sui$p.adjust < 0.05, na.rm = TRUE), "\n\n")


## =============================================================================
## STEP 4: Shared pathways, both subsets
## =============================================================================
cat("================ STEP 4: Shared pathways (both 6v6 subsets vs SUI) ================\n\n")

report_shared <- function(gsea_pop, label, fdr_cut = 0.25) {
  pop_sig <- subset(gsea_pop, p.adjust < fdr_cut)
  sui_sig <- subset(gsea_sui, p.adjust < fdr_cut)
  shared_ids <- intersect(pop_sig$PATH, sui_sig$PATH)
  cat("---", label, ", FDR<", fdr_cut, ": POP sig =", nrow(pop_sig),
      "| SUI sig =", nrow(sui_sig), "| SHARED =", length(shared_ids), "---\n")
  if (length(shared_ids) == 0) return(data.frame())
  out <- merge(pop_sig[pop_sig$PATH %in% shared_ids, c("PATH","PathwayName","Nh","NES","p.adjust")],
               sui_sig[sui_sig$PATH %in% shared_ids, c("PATH","Nh","NES","p.adjust")],
               by = "PATH", suffixes = c("_POP", "_SUI"))
  out$Same_direction <- sign(out$NES_POP) == sign(out$NES_SUI)
  out <- out[order(out$p.adjust_POP), ]
  n_na <- sum(is.na(out$Same_direction))
  cat("Same direction:", sum(out$Same_direction, na.rm = TRUE), "of",
      sum(!is.na(out$Same_direction)), "comparable (", n_na, "NA )\n\n")
  out
}

shared_Y <- report_shared(gsea_Y$gsea, "POP_Y")
shared_D <- report_shared(gsea_D$gsea, "POP_D")
if (nrow(shared_Y) > 0) write.csv(shared_Y, "results/22_shared_pathways_6v6_POPY_FDR025.csv", row.names = FALSE)
if (nrow(shared_D) > 0) write.csv(shared_D, "results/22_shared_pathways_6v6_POPD_FDR025.csv", row.names = FALSE)

cat("=== SUMMARY: comparing to script 02 (6 vs 12, unbalanced) ===\n")
cat("Script 02 (6v12): 25 shared at FDR<0.25, 15/18 (83%) concordant\n")
cat("This script, POP_Y (6v6):", nrow(shared_Y), "shared\n")
cat("This script, POP_D (6v6):", nrow(shared_D), "shared\n\n")


## =============================================================================
## STEP 5: GRAPHICS
## =============================================================================
cat("================ STEP 5: Graphics ================\n\n")

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
ggsave("figures/28_GSEA_barplot_POP_GSE208261_6v6_POPY.png",
       make_gsea_barplot(gsea_Y$gsea, "GSEA (classic) - POP (GSE208261, 6v6, POP_Y half)"),
       width = 10, height = 6, dpi = 300)
ggsave("figures/29_GSEA_barplot_POP_GSE208261_6v6_POPD.png",
       make_gsea_barplot(gsea_D$gsea, "GSEA (classic) - POP (GSE208261, 6v6, POP_D half)"),
       width = 10, height = 6, dpi = 300)
cat("Saved: figures/28_*.png, figures/29_*.png\n\n")

plot_shared <- function(shared_df, label, fname) {
  if (nrow(shared_df) == 0) { cat("SKIPPED", fname, "(no shared pathways)\n"); return(invisible(NULL)) }
  p <- ggplot(shared_df, aes(x = NES_POP, y = NES_SUI, label = PathwayName)) +
    geom_hline(yintercept = 0, color = "grey70") + geom_vline(xintercept = 0, color = "grey70") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
    geom_point(aes(color = Same_direction), size = 3) +
    scale_color_manual(values = c("TRUE" = "#2166AC", "FALSE" = "#B2182B"),
                        labels = c("TRUE" = "Same direction", "FALSE" = "Opposite direction")) +
    geom_text_repel(size = 3, max.overlaps = 30) +
    labs(title = paste0("Shared pathways (FDR<0.25): NES in POP (", label, ") vs NES in SUI"),
         x = "NES - POP (classic GSEA)", y = "NES - SUI (preranked GSEA)", color = NULL) +
    theme_bw() + theme(legend.position = "top")
  ggsave(fname, p, width = 9, height = 7, dpi = 300)
  cat("Saved:", fname, "\n")
}
plot_shared(shared_Y, "GSE208261 6v6, POP_Y", "figures/30_shared_pathways_6v6_POPY.png")
plot_shared(shared_D, "GSE208261 6v6, POP_D", "figures/31_shared_pathways_6v6_POPD.png")

cat("\n=== DONE (script 06) ===\n")
