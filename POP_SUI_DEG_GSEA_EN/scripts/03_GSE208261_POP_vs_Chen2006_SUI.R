# =============================================================================
# GSE208261 (POP, RNA-seq) vs Chen 2006 (SUI, literature DEG panel):
# GSEA and shared pathways
# =============================================================================
# A third comparison in this project: same POP dataset as script 02
# (GSE208261, RNA-seq), now against a DIFFERENT SUI source - Chen B, Wen Y,
# Zhang Z, Guo Y, Warrington JA, Polan ML, "Microarray analysis of
# differentially expressed genes in vaginal tissues from women with stress
# urinary incontinence compared with asymptomatic women." Hum Reprod.
# 2006;21(1):22-29. 5 SUI vs 5 continent pairs, periurethral vaginal wall,
# Affymetrix U133A. The PDF only reports the article's own FINAL DEG list
# (Tables II/III: 79 genes total, common to their MAS5.0 and RMA pipelines,
# p<0.05 after multiple-testing correction) - not the full array, so, like
# Wei2020 in script 01/02, this already IS the DEG list, not raw data we
# can re-threshold. `data/Chen2006_SUI_panel.csv` curates 60 of the 79
# (34 up / 26 down) with a current, confidently-mapped HGNC symbol; the
# other 19 use 2005-era ambiguous names ("hypothetical protein FLJxxxxx",
# generic "Zinc finger protein", etc.) and are excluded rather than guessed.
#
# METHOD - why Chen2006 can ONLY use preranked GSEA, never classic:
# classic (phenotype) GSEA needs a full per-sample expression matrix to
# permute - the original study's raw array data was never published, only
# this final 60-gene table (p-value + fold-change per gene, no per-sample
# values). Preranked (gene-set-label permutation) is the only option here,
# same as the reasoning already used for the small-panel Chen sources in
# the wider project. Ranking score = sign(direction) x -log10(p-value)
# (RMA p-value column; the 2 pipelines the paper reports agree on which
# genes make the final list, RMA is used here for the continuous p-value).
# A KNOWN LIMITATION of using preranked GSEA on a list this short: a KEGG
# pathway is only testable if >=2 of its members happen to be among these
# 60 genes - most pathways will have 0 or 1 by chance, so few are testable,
# and this is a non-standard, lower-powered use of GSEA (designed for
# whole-transcriptome ranked lists, not a pre-filtered 60-gene panel).
# Read these numbers as an exploratory, low-power signal - not on the same
# footing as the GSE208261 classic GSEA (whole transcriptome, phenotype
# permutation) it is compared against below.
#
# METHOD - POP (GSE208261): reuses the DEG and classic-GSEA results ALREADY
# COMPUTED by script 02 (`results/06_POP_GSE208261_DESeq2_full.csv`,
# `results/07_GSEA_classic_POP_GSE208261_KEGG.csv`) rather than recomputing
# the ~500-permutation voom+limma loop again (identical computation,
# ~18 minutes) - if you have not run script 02 yet, run it first. Only the
# Chen2006 side and the cross-comparison are computed fresh here.
# =============================================================================

suppressMessages({
  library(org.Hs.eg.db)
  library(ggplot2)
  library(pheatmap)
})
if (!requireNamespace("ggrepel", quietly = TRUE)) install.packages("ggrepel")
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
## STEP 1: Load and describe the Chen 2006 SUI panel
## =============================================================================
cat("\n================ STEP 1: SUI panel (Chen 2006) ================\n\n")

chen <- read.csv("data/Chen2006_SUI_panel.csv")
chen <- chen[!duplicated(chen$GeneSymbol), ]
chen$logFC <- log2(chen$RMA_FoldChange)  # RMA_FoldChange is a signed ratio here (unlike Wei2020's GeneSpring export) - <1 for down genes, so log2() alone gives the correct sign, no fix needed.
write.csv(chen, "results/09_SUI_Chen2006_panel.csv", row.names = FALSE)

cat("Chen 2006 (Hum Reprod): 5 SUI vs 5 continent pairs, periurethral vaginal\n")
cat("wall, Affymetrix U133A. This CSV is the article's own final DEG list\n")
cat("(Tables II/III, common to their MAS5.0 and RMA pipelines, p<0.05) -\n")
cat("60 of 79 reported genes with a confidently-mapped current HGNC symbol.\n")
cat("DEG SUI (Chen 2006, authors' criterion):", nrow(chen), "(",
    sum(chen$Direction == "up"), "up /", sum(chen$Direction == "down"), "down )\n\n")
print(head(chen[, c("GeneSymbol","Direction","RMA_pvalue","RMA_FoldChange","logFC")], 10))
cat("\n")


## =============================================================================
## STEP 2: GSEA for SUI (Chen 2006) - PRERANKED ONLY (see header note)
## =============================================================================
cat("================ STEP 2: GSEA for SUI (Chen 2006) - PRERANKED only ================\n\n")
cat("Classic (phenotype) GSEA is not possible here - the original study's raw\n")
cat("array data was never published, only this final 60-gene table. Ranking\n")
cat("score = sign(direction) x -log10(RMA p-value); p=0 rows floored to half\n")
cat("the smallest non-zero p-value in the panel to avoid an infinite score.\n\n")

build_score <- function(direction, pvalue) {
  p <- pvalue
  floor_p <- min(p[p > 0], na.rm = TRUE) / 2
  p[p <= 0] <- floor_p
  sign <- ifelse(direction == "up", 1, -1)
  sign * -log10(p)
}
score_chen <- build_score(chen$Direction, chen$RMA_pvalue)
ord_c <- order(-score_chen)
ranked_genes_chen <- chen$GeneSymbol[ord_c]
ranked_scores_chen <- score_chen[ord_c]
N_chen <- length(ranked_genes_chen)

all_sym <- keys(org.Hs.eg.db, keytype = "SYMBOL")
ann_path_all <- suppressWarnings(select(org.Hs.eg.db, keys = all_sym, keytype = "SYMBOL", columns = "PATH"))
ann_path_all <- ann_path_all[!is.na(ann_path_all$PATH), ]
kegg_sets_all <- split(ann_path_all$SYMBOL, ann_path_all$PATH)
cat("KEGG pathways loaded (full universe, org.Hs.eg.db):", length(kegg_sets_all), "\n")

min_gs_chen <- 2
testable_chen <- lapply(kegg_sets_all, function(g) which(ranked_genes_chen %in% g))
testable_chen <- testable_chen[sapply(testable_chen, length) >= min_gs_chen]
cat("KEGG pathways testable (>=", min_gs_chen, "members present in this", N_chen, "-gene list):",
    length(testable_chen), "\n")

scores_abs_chen <- abs(ranked_scores_chen)
es_obs_chen <- sapply(testable_chen, calc_es, scores_abs = scores_abs_chen, N = N_chen)

set.seed(123)
n_perm5 <- 1000
cat("Running", n_perm5, "gene-set-label permutations...\n")
perm_mat_chen <- matrix(NA_real_, nrow = n_perm5, ncol = length(testable_chen))
for (i in seq_len(n_perm5)) {
  for (j in seq_along(testable_chen)) {
    hidx <- sample.int(N_chen, length(testable_chen[[j]]))
    perm_mat_chen[i, j] <- calc_es(hidx, scores_abs_chen, N_chen)
  }
}

pval_chen <- numeric(length(testable_chen)); nes_chen <- numeric(length(testable_chen))
for (j in seq_along(testable_chen)) {
  pe <- perm_mat_chen[, j]; pe <- pe[!is.na(pe)]
  if (is.na(es_obs_chen[j])) { pval_chen[j] <- NA; nes_chen[j] <- NA; next }
  if (es_obs_chen[j] >= 0) {
    pval_chen[j] <- (sum(pe >= es_obs_chen[j]) + 1) / (length(pe) + 1)
    base <- mean(pe[pe >= 0]); if (is.nan(base) || base == 0) base <- NA
  } else {
    pval_chen[j] <- (sum(pe <= es_obs_chen[j]) + 1) / (length(pe) + 1)
    base <- mean(abs(pe[pe < 0])); if (is.nan(base) || base == 0) base <- NA
  }
  nes_chen[j] <- es_obs_chen[j] / base
}

leading_edge_chen <- sapply(testable_chen, function(idx) paste(ranked_genes_chen[idx], collapse = "/"))
gsea_chen <- data.frame(PATH = names(testable_chen), Nh = sapply(testable_chen, length),
                         ES = es_obs_chen, NES = nes_chen, pvalue = pval_chen, leadingEdge = leading_edge_chen)
gsea_chen$p.adjust <- p.adjust(gsea_chen$pvalue, "BH")
gsea_chen$PathwayName <- kegg_label(gsea_chen$PATH)
gsea_chen <- gsea_chen[order(gsea_chen$pvalue), ]
write.csv(gsea_chen, "results/10_GSEA_preranked_Chen2006_KEGG.csv", row.names = FALSE)

cat("\n=== GSEA preranked, SUI Chen 2006 (KEGG) ===\n")
cat("Pathways significant at FDR<0.25:", sum(gsea_chen$p.adjust < 0.25, na.rm = TRUE), "of", nrow(gsea_chen), "\n")
cat("Pathways significant at FDR<0.05:", sum(gsea_chen$p.adjust < 0.05, na.rm = TRUE), "\n\n")
print(head(gsea_chen[, c("PathwayName","Nh","NES","pvalue","p.adjust")], 10))
cat("\n")


## =============================================================================
## STEP 3: Load POP (GSE208261) results from script 02
## =============================================================================
cat("================ STEP 3: POP (GSE208261) - loaded from script 02 ================\n\n")

pop2_file <- "results/06_POP_GSE208261_DESeq2_full.csv"
gsea_pop2_file <- "results/07_GSEA_classic_POP_GSE208261_KEGG.csv"
if (!file.exists(pop2_file) || !file.exists(gsea_pop2_file)) {
  stop("Run scripts/02_GSE208261_POP_RNAseq_vs_SUI.R first - this script reuses ",
       "its POP DEG and classic-GSEA results instead of recomputing the ",
       "~18-minute permutation loop again.")
}
pop2_full <- read.csv(pop2_file)
gsea_pop2 <- read.csv(gsea_pop2_file)
ranked_genes_pop2 <- pop2_full$Gene
cat("Loaded POP (GSE208261) DEG table:", nrow(pop2_full), "genes tested\n")
cat("Loaded POP (GSE208261) classic GSEA:", nrow(gsea_pop2), "KEGG pathways tested,",
    sum(gsea_pop2$p.adjust < 0.05, na.rm = TRUE), "significant at FDR<0.05\n\n")


## =============================================================================
## STEP 4: Shared pathways (POP GSE208261 x SUI Chen2006) + gene direction table
## =============================================================================
cat("================ STEP 4: Shared pathways GSE208261(POP) x Chen2006(SUI) ================\n\n")

report_shared3 <- function(fdr_cut) {
  pop_sig <- subset(gsea_pop2, p.adjust < fdr_cut)
  sui_sig <- subset(gsea_chen, p.adjust < fdr_cut)
  shared_ids <- intersect(pop_sig$PATH, sui_sig$PATH)
  cat("--- FDR <", fdr_cut, ": POP(GSE208261) significant =", nrow(pop_sig),
      "| SUI(Chen2006) significant =", nrow(sui_sig), "| SHARED =", length(shared_ids), "---\n")
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

shared3_025 <- report_shared3(0.25)
shared3_005 <- report_shared3(0.05)
if (nrow(shared3_025) > 0) write.csv(shared3_025, "results/11_shared_pathways_GSE208261xChen2006_FDR025.csv", row.names = FALSE)
if (nrow(shared3_005) > 0) write.csv(shared3_005, "results/11_shared_pathways_GSE208261xChen2006_FDR005.csv", row.names = FALSE)

shared3_for_table <- if (nrow(shared3_005) > 0) shared3_005 else shared3_025
cat("\n--- Gene-level direction table for the shared pathways ---\n")
if (nrow(shared3_for_table) == 0) {
  cat("No shared pathway at either threshold - no gene-level table to build.\n\n")
  gene_dir_table3 <- data.frame()
} else {
  pop_lookup3 <- setNames(pop2_full$logFC, pop2_full$Gene)
  sui_lookup3 <- setNames(chen$logFC, chen$GeneSymbol)
  pop3_tested <- ranked_genes_pop2
  sui_tested3 <- ranked_genes_chen

  rows3 <- lapply(seq_len(nrow(shared3_for_table)), function(i) {
    pid <- shared3_for_table$PATH[i]
    pname <- shared3_for_table$PathwayName[i]
    full_members <- suppressWarnings(select(org.Hs.eg.db, keys = pid, keytype = "PATH", columns = "SYMBOL")$SYMBOL)
    members <- intersect(unique(full_members), union(pop3_tested, sui_tested3))
    data.frame(
      Pathway = pname, KEGG_ID = pid, Gene = members,
      logFC_POP = unname(pop_lookup3[members]), logFC_SUI = unname(sui_lookup3[members]),
      stringsAsFactors = FALSE)
  })
  gene_dir_table3 <- do.call(rbind, rows3)
  gene_dir_table3$Direction_POP <- ifelse(is.na(gene_dir_table3$logFC_POP), "not tested",
                                           ifelse(gene_dir_table3$logFC_POP > 0, "up", "down"))
  gene_dir_table3$Direction_SUI <- ifelse(is.na(gene_dir_table3$logFC_SUI), "not tested",
                                           ifelse(gene_dir_table3$logFC_SUI > 0, "up", "down"))
  gene_dir_table3$Concordant <- with(gene_dir_table3,
                                      Direction_POP != "not tested" & Direction_SUI != "not tested" &
                                        Direction_POP == Direction_SUI)
  gene_dir_table3 <- gene_dir_table3[order(gene_dir_table3$Pathway, -gene_dir_table3$Concordant), ]
  write.csv(gene_dir_table3, "results/11_shared_pathways_gene_direction_table.csv", row.names = FALSE)

  both_tested3 <- subset(gene_dir_table3, Direction_POP != "not tested" & Direction_SUI != "not tested")
  cat("Genes in shared pathways tested in both diseases:", nrow(both_tested3),
      "| concordant direction:", sum(both_tested3$Concordant),
      "(", round(100 * mean(both_tested3$Concordant), 1), "% )\n\n")
  print(head(gene_dir_table3, 20))
}
cat("\n")


## =============================================================================
## STEP 4b: Direct gene-level concordance (whole 59-gene panel, not limited
##           to shared-pathway members) - the more sensitive test for a
##           panel this short, and the same check already used elsewhere in
##           this project for Chen2006 against other POP datasets.
## =============================================================================
cat("================ STEP 4b: Direct gene-level concordance (all 59 Chen2006 genes) ================\n\n")

pop_lookup_all <- setNames(pop2_full$logFC, pop2_full$Gene)
chen_cross <- chen[chen$GeneSymbol %in% names(pop_lookup_all), ]
chen_cross$logFC_POP <- unname(pop_lookup_all[chen_cross$GeneSymbol])
chen_cross$Direction_POP <- ifelse(chen_cross$logFC_POP > 0, "up", "down")
chen_cross$Concordant <- chen_cross$Direction == chen_cross$Direction_POP
write.csv(chen_cross, "results/12_Chen2006_x_GSE208261_gene_concordance.csv", row.names = FALSE)

n_concordant <- sum(chen_cross$Concordant)
n_tested <- nrow(chen_cross)
bt <- binom.test(n_concordant, n_tested, p = 0.5)
cat("Chen2006 genes tested in GSE208261:", n_tested, "of", nrow(chen), "\n")
cat("Concordant direction:", n_concordant, "of", n_tested, "(",
    round(100 * n_concordant / n_tested, 1), "% )\n")
cat("Binomial sign test vs 50%: p =", signif(bt$p.value, 3), "\n\n")
print(chen_cross[, c("GeneSymbol","Direction","Direction_POP","Concordant","logFC","logFC_POP")])
cat("\n")


## =============================================================================
## STEP 5: GRAPHICS
## =============================================================================
cat("================ STEP 5: Graphics ================\n\n")

# Volcano-style plot of the Chen2006 panel (all 59 genes are already "DEG"
# by the authors' own criterion - this shows their fold-change/p-value
# spread and labels every gene, since the panel is short enough to do so).
p_chen_volcano <- ggplot(chen, aes(x = logFC, y = -log10(RMA_pvalue), color = Direction)) +
  geom_point(size = 2.2, alpha = 0.8) +
  scale_color_manual(values = c(down = "#2166AC", up = "#B2182B")) +
  geom_vline(xintercept = 0, color = "grey40", linetype = "dashed") +
  geom_text_repel(aes(label = GeneSymbol), size = 2.8, color = "black", max.overlaps = 60) +
  labs(title = "SUI panel - Chen 2006 (59 genes, authors' own final DEG list)",
       x = "log2(Fold Change)", y = expression(-log[10](italic(p)~value)), color = NULL) +
  theme_bw() + theme(legend.position = "top", plot.title = element_text(face = "bold"))
ggsave("figures/15_volcano_SUI_Chen2006.png", p_chen_volcano, width = 9, height = 7, dpi = 300)
cat("Saved: figures/15_volcano_SUI_Chen2006.png\n\n")

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
p_bar_chen <- make_gsea_barplot(gsea_chen, "GSEA (preranked) - KEGG pathways in SUI (Chen 2006)\n(exploratory - only 19 pathways testable in a 59-gene panel)")
ggsave("figures/16_GSEA_barplot_SUI_Chen2006.png", p_bar_chen, width = 10, height = 6, dpi = 300)
cat("Saved: figures/16_GSEA_barplot_SUI_Chen2006.png\n\n")

# Gene-level concordance plot (STEP 4b): every Chen2006 gene's logFC in
# Chen2006 vs its logFC in GSE208261, colored by direction agreement.
p_concord <- ggplot(chen_cross, aes(x = logFC, y = logFC_POP, label = GeneSymbol)) +
  geom_hline(yintercept = 0, color = "grey70") + geom_vline(xintercept = 0, color = "grey70") +
  geom_point(aes(color = Concordant), size = 2.5) +
  scale_color_manual(values = c("TRUE" = "#2166AC", "FALSE" = "#B2182B"),
                      labels = c("TRUE" = "Same direction", "FALSE" = "Opposite direction")) +
  geom_text_repel(size = 2.6, max.overlaps = 40) +
  labs(title = paste0("Gene-level direction: Chen2006 (SUI) vs GSE208261 (POP)\n",
                       n_concordant, " of ", n_tested, " genes concordant (",
                       round(100 * n_concordant / n_tested, 1), "%), binomial p=", signif(bt$p.value, 3)),
       x = "log2(Fold Change) - Chen 2006 (SUI)", y = "log2(Fold Change) - GSE208261 (POP)", color = NULL) +
  theme_bw() + theme(legend.position = "top")
ggsave("figures/17_gene_concordance_Chen2006_GSE208261.png", p_concord, width = 9, height = 7.5, dpi = 300)
cat("Saved: figures/17_gene_concordance_Chen2006_GSE208261.png\n\n")

if (nrow(shared3_025) > 0) {
  p_shared3 <- ggplot(shared3_025, aes(x = NES_POP, y = NES_SUI, label = PathwayName)) +
    geom_hline(yintercept = 0, color = "grey70") + geom_vline(xintercept = 0, color = "grey70") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
    geom_point(aes(color = Same_direction), size = 3) +
    scale_color_manual(values = c("TRUE" = "#2166AC", "FALSE" = "#B2182B"),
                        labels = c("TRUE" = "Same direction", "FALSE" = "Opposite direction")) +
    geom_text_repel(size = 3, max.overlaps = 30) +
    labs(title = "Shared pathways (FDR<0.25): NES in POP (GSE208261) vs NES in SUI (Chen2006)",
         x = "NES - POP GSE208261 (classic GSEA)", y = "NES - SUI Chen2006 (preranked GSEA)", color = NULL) +
    theme_bw() + theme(legend.position = "top")
  ggsave("figures/18_shared_pathways_NES_comparison_Chen2006.png", p_shared3, width = 9, height = 7, dpi = 300)
  cat("Saved: figures/18_shared_pathways_NES_comparison_Chen2006.png\n\n")
} else {
  cat("SKIPPED: figures/18_shared_pathways_NES_comparison_Chen2006.png (no shared pathway at FDR<0.25)\n\n")
}

cat("\n=== ALL GSE208261 x Chen2006 FIGURES SAVED TO figures/ ===\n")
cat("15_volcano_SUI_Chen2006.png\n16_GSEA_barplot_SUI_Chen2006.png\n")
cat("17_gene_concordance_Chen2006_GSE208261.png\n18_shared_pathways_NES_comparison_Chen2006.png (if any shared)\n")
cat("\n=== DONE (script 03) ===\n")
