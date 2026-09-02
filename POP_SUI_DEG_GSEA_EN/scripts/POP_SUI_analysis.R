# =============================================================================
# POP vs SUI: Differential Expression, GSEA, and Shared-Pathway Analysis
# =============================================================================
# Self-contained script. Run it top to bottom from this folder
# (POP_SUI_DEG_GSEA_EN) in RStudio (Session > Set Working Directory >
# To Source File Location, or just open the .Rproj/this folder).
#
# Required packages (Bioconductor + CRAN):
#   limma, org.Hs.eg.db, readxl, ggplot2, ggrepel, pheatmap
# Install once with:
#   install.packages(c("readxl","ggplot2","ggrepel","pheatmap","BiocManager"))
#   BiocManager::install(c("limma","org.Hs.eg.db"))
#
# -----------------------------------------------------------------------------
# DATASETS
# -----------------------------------------------------------------------------
# POP: GSE53868 (Kerkhof et al.) - 12 women with pelvic organ prolapse,
#      PAIRED biopsy per patient (prolapse-site tissue vs non-prolapse-site
#      tissue from the same woman), Agilent 4x44K array, already
#      log2-normalized. File: data/GSE53868_series_matrix.txt
# SUI: Wei et al. 2020 (Reprod Sci) Supplementary Table S2 - 3 women with
#      stress urinary incontinence vs 3 continent controls, periurethral
#      vaginal wall tissue, Arraystar Human LncRNA+mRNA V4.0 array. This
#      file contains TWO sheets ("up_Sui_vs_Ctrl" and "down_Sui_vs_Ctrl")
#      and already only lists the genes the authors called differentially
#      expressed (fold-change>=2, raw p<0.05) - NOT the full tested array.
#      File: data/Wei2020_TableS2_mRNA.xls
#
# -----------------------------------------------------------------------------
# METHODS (and why)
# -----------------------------------------------------------------------------
# DEG:
#   - POP: standard limma moderated t-test, paired design (~individual +
#     tissue), |log2FC|>1 and FDR(BH)<0.05.
#   - SUI: the file we have IS already the authors' DEG list (fold-change
#     >=2, raw p<0.05) - the full untested array was not published, so we
#     cannot re-derive a POP-style FDR-based DEG list from scratch (doing so
#     from an already-pre-filtered list would be statistically circular).
#     We report their list as the primary SUI DEG set, and additionally
#     apply FDR<0.05 to their own precomputed FDR column (real, calculated
#     on the full array by the authors' software before their p<0.05 filter
#     was applied - not a circular re-derivation) as a stricter,
#     POP-comparable subset.
#
# GSEA:
#   - POP: CLASSIC GSEA (phenotype permutation - paired sign-flip across the
#     12 patients). This is the statistically preferred method whenever a
#     full sample-level expression matrix is available (Subramanian et al.
#     2005), and 12 paired patients give reasonable permutation resolution.
#   - SUI: PRERANKED GSEA (gene-set-label permutation), ranked by a REAL
#     moderated t-statistic computed from the per-sample intensities in the
#     supplementary table (3 SUI vs 3 Ctrl), not just the authors' summary
#     p-value. Classic phenotype permutation is NOT used for SUI because
#     with only 3 vs 3 samples there are only choose(6,3) = 20 possible
#     relabelings, capping resolution at p = 1/20 = 0.05 and making
#     FDR < 0.05 structurally unreachable no matter how real the signal is.
#     Preranked (gene-set) permutation does not have this ceiling because it
#     permutes over the ~200 pathways being tested, not over the 6 samples.
#   Both use the same weighted running-sum enrichment score (ES, weight=1)
#   from Subramanian et al. 2005, PNAS, reimplemented here because this
#   environment has no live access to CRAN/Bioconductor to install
#   fgsea/msigdbr; it uses the KEGG pathway annotation already bundled
#   offline in org.Hs.eg.db (works with no internet, including in RStudio).
#   Gene sets are labelled by their KEGG numeric ID with a small built-in
#   name lookup for readability in the figures/tables.
# =============================================================================

suppressMessages({
  library(limma)
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

# A small lookup of common KEGG pathway names, for readable figures/tables
# (org.Hs.eg.db only stores pathway IDs, not names, and this environment has
# no live KEGG access to fetch the full name table).
kegg_names <- c(
  "00190"="Oxidative phosphorylation","01100"="Metabolic pathways",
  "03010"="Ribosome","03013"="RNA transport","03015"="mRNA surveillance pathway",
  "03018"="RNA degradation","03020"="RNA polymerase","03022"="Basal transcription factors",
  "03030"="DNA replication","03040"="Spliceosome","03050"="Proteasome",
  "03060"="Protein export","03320"="PPAR signaling pathway",
  "04010"="MAPK signaling pathway","04012"="ErbB signaling pathway",
  "04020"="Calcium signaling pathway","04060"="Cytokine-cytokine receptor interaction",
  "04062"="Chemokine signaling pathway","04066"="HIF-1 signaling pathway",
  "04068"="FoxO signaling pathway","04110"="Cell cycle",
  "04115"="p53 signaling pathway","04120"="Ubiquitin mediated proteolysis",
  "04140"="Autophagy","04141"="Protein processing in endoplasmic reticulum",
  "04142"="Lysosome","04144"="Endocytosis","04145"="Phagosome",
  "04150"="mTOR signaling pathway","04151"="PI3K-Akt signaling pathway",
  "04210"="Apoptosis","04310"="Wnt signaling pathway",
  "04330"="Notch signaling pathway","04340"="Hedgehog signaling pathway",
  "04350"="TGF-beta signaling pathway","04370"="VEGF signaling pathway",
  "04510"="Focal adhesion","04512"="ECM-receptor interaction",
  "04520"="Adherens junction","04530"="Tight junction","04540"="Gap junction",
  "04610"="Complement and coagulation cascades","04620"="Toll-like receptor signaling pathway",
  "04621"="NOD-like receptor signaling pathway","04630"="JAK-STAT signaling pathway",
  "04660"="T cell receptor signaling pathway","04662"="B cell receptor signaling pathway",
  "04664"="Fc epsilon RI signaling pathway","04670"="Leukocyte transendothelial migration",
  "04810"="Regulation of actin cytoskeleton","04910"="Insulin signaling pathway",
  "04914"="Progesterone-mediated oocyte maturation","05010"="Alzheimer disease",
  "05012"="Parkinson disease","05014"="Amyotrophic lateral sclerosis",
  "05016"="Huntington disease","05020"="Prion diseases","05140"="Leishmaniasis",
  "05142"="Chagas disease (American trypanosomiasis)","05144"="Malaria",
  "05145"="Toxoplasmosis","05146"="Amoebiasis","05160"="Hepatitis C",
  "05161"="Hepatitis B","05164"="Influenza A","05166"="HTLV-I infection",
  "05169"="Epstein-Barr virus infection","05200"="Pathways in cancer",
  "05203"="Viral carcinogenesis","05205"="Proteoglycans in cancer",
  "05219"="Bladder cancer","05222"="Small cell lung cancer",
  "05223"="Non-small cell lung cancer")
kegg_label <- function(id) {
  nm <- kegg_names[id]
  ifelse(is.na(nm), paste0("KEGG ", id), paste0("KEGG ", id, " - ", nm))
}


## =============================================================================
## STEP 1: DEG for POP (GSE53868)
## =============================================================================
cat("\n================ STEP 1: DEG for POP (GSE53868) ================\n\n")

raw_lines <- readLines("data/GSE53868_series_matrix.txt")
start_row <- grep("^!series_matrix_table_begin", raw_lines) + 1
end_row   <- grep("^!series_matrix_table_end", raw_lines) - 1
expr_pop <- read.delim("data/GSE53868_series_matrix.txt", skip = start_row - 1,
                        nrows = end_row - start_row, header = TRUE,
                        row.names = 1, check.names = FALSE, quote = "\"")

sample_title_line <- raw_lines[grep("^!Sample_title", raw_lines)]
sample_titles <- gsub('"', "", strsplit(sample_title_line, "\t")[[1]][-1])
individual_line <- raw_lines[grep("^!Sample_characteristics_ch1.*individual:", raw_lines)][1]
individuals <- gsub("individual: ", "", gsub('"', "", strsplit(individual_line, "\t")[[1]][-1]))
tissue <- ifelse(grepl("\\(POP site\\)", sample_titles), "POP_site", "NonPOP_site")

coldata_pop <- data.frame(row.names = colnames(expr_pop),
                           tissue = factor(tissue, levels = c("NonPOP_site", "POP_site")),
                           individual = factor(individuals))

design_pop <- model.matrix(~ individual + tissue, data = coldata_pop)
fit_pop <- eBayes(lmFit(as.matrix(expr_pop), design_pop))
pop_full <- topTable(fit_pop, coef = "tissuePOP_site", number = Inf, sort.by = "P")
pop_full$Gene <- rownames(pop_full)
pop_full <- pop_full[, c("Gene", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val")]
write.csv(pop_full, "results/01_POP_limma_full.csv", row.names = FALSE)

pop_deg <- subset(pop_full, adj.P.Val < 0.05 & abs(logFC) > 1)
pop_deg <- pop_deg[order(pop_deg$adj.P.Val), ]
write.csv(pop_deg, "results/01_POP_DEG_logFC1_FDR05.csv", row.names = FALSE)

cat("Design: 12 women with POP, paired biopsy per patient (prolapse site vs\n")
cat("non-prolapse site). Model: ~individual + tissue (paired moderated t-test).\n")
cat("Genes tested:", nrow(pop_full), "\n")
cat("DEG (|log2FC|>1, FDR<0.05):", nrow(pop_deg), "(", sum(pop_deg$logFC > 0), "up /",
    sum(pop_deg$logFC < 0), "down )\n\n")
print(head(pop_deg, 10))
cat("\n")


## =============================================================================
## STEP 2: DEG for SUI (Wei et al. 2020) - reading BOTH sheets
## =============================================================================
cat("================ STEP 2: DEG for SUI (Wei 2020) ================\n\n")

wei_xls <- "data/Wei2020_TableS2_mRNA.xls"
up_sheet   <- read_excel(wei_xls, sheet = "up_Sui_vs_Ctrl",   skip = 17)
down_sheet <- read_excel(wei_xls, sheet = "down_Sui_vs_Ctrl", skip = 17)
cat("Sheet 'up_Sui_vs_Ctrl'  :", nrow(up_sheet), "probes\n")
cat("Sheet 'down_Sui_vs_Ctrl':", nrow(down_sheet), "probes\n")

up_sheet$Direction <- "up"; down_sheet$Direction <- "down"
sample_cols <- c("[Sui1, Sui](normalized)", "[Sui2, Sui](normalized)", "[Sui3, Sui](normalized)",
                  "[Ctrl1, Ctrl](normalized)", "[Ctrl2, Ctrl](normalized)", "[Ctrl3, Ctrl](normalized)")
keep_cols <- c("GeneSymbol", "P-value", "FDR", "Fold Change", "Direction", sample_cols)
wei_both <- rbind(as.data.frame(up_sheet[, keep_cols]), as.data.frame(down_sheet[, keep_cols]))
colnames(wei_both) <- c("GeneSymbol", "PValue", "FDR", "FoldChange", "Direction",
                         "Sui1", "Sui2", "Sui3", "Ctrl1", "Ctrl2", "Ctrl3")
wei_both <- wei_both[!is.na(wei_both$GeneSymbol), ]
cat("Combined (both sheets):", nrow(wei_both), "probes\n")

# IMPORTANT: the authors' "Fold Change" column is a MAGNITUDE only (always
# >=2 in BOTH sheets) - the direction/sign comes from the sheet itself
# (Direction column), not from the numeric value. We restore a correctly
# signed log2FC here (negative for down-regulated genes).
wei_both$logFC <- ifelse(wei_both$Direction == "down",
                          -log2(wei_both$FoldChange), log2(wei_both$FoldChange))

# Collapse duplicate probes per gene (keep the one with the lowest p-value)
wei_both <- wei_both[order(wei_both$PValue), ]
sui_full <- wei_both[!duplicated(wei_both$GeneSymbol), ]
rownames(sui_full) <- NULL
write.csv(sui_full, "results/02_SUI_all_genes_both_sheets.csv", row.names = FALSE)
cat("Unique genes after collapsing duplicate probes:", nrow(sui_full), "\n\n")

cat("NOTE on the SUI DEG definition: this supplementary table already ONLY\n")
cat("contains the genes the authors classified as differentially expressed\n")
cat("(fold-change>=2, i.e. |log2FC|>1, raw p<0.05); the full tested array was\n")
cat("not published, so this entire table already IS the DEG list, by the\n")
cat("authors' own original criterion.\n\n")

sui_deg <- sui_full
write.csv(sui_deg, "results/02_SUI_DEG_primary.csv", row.names = FALSE)
cat("DEG SUI (authors' original criterion, fold-change>=2 & raw p<0.05):",
    nrow(sui_deg), "(", sum(sui_deg$Direction == "up"), "up /",
    sum(sui_deg$Direction == "down"), "down )\n\n")

cat("For an apples-to-apples check against POP's FDR<0.05 cutoff: applying\n")
cat("FDR<0.05 to the authors' own precomputed FDR column (calculated by their\n")
cat("software - GeneSpring GX - on the full array BEFORE their p<0.05 filter\n")
cat("was applied; it has real variation, from", round(min(sui_deg$FDR), 4), "to",
    round(max(sui_deg$FDR), 4), "- confirming it is not a circular re-derivation)\n")
sui_deg_fdr <- subset(sui_deg, FDR < 0.05)
write.csv(sui_deg_fdr, "results/02_SUI_DEG_FDR005_subset.csv", row.names = FALSE)
cat("gives a stricter subset of", nrow(sui_deg_fdr), "genes (",
    sum(sui_deg_fdr$Direction == "up"), "up /", sum(sui_deg_fdr$Direction == "down"),
    "down ) - results/02_SUI_DEG_FDR005_subset.csv.\n")
cat("The primary SUI DEG list used below is the full", nrow(sui_deg), "-gene\n")
cat("authors' list (results/02_SUI_DEG_primary.csv), matching how the original\n")
cat("paper itself defines its DEGs.\n\n")
print(head(sui_deg[, c("GeneSymbol","logFC","PValue","FDR","Direction")], 10))
cat("\n")


## =============================================================================
## STEP 3: GSEA for POP - CLASSIC (phenotype permutation)
## =============================================================================
cat("================ STEP 3: GSEA for POP - CLASSIC method ================\n\n")
cat("Method: classic GSEA, phenotype permutation (paired sign-flip across the\n")
cat("12 patients) - the statistically preferred method here because we have\n")
cat("the full per-sample expression matrix and enough paired replicates for\n")
cat("real permutation resolution. Ranking = paired moderated t-statistic.\n\n")

set.seed(42)
expr_mat_pop <- as.matrix(expr_pop)
pop_cols    <- colnames(expr_mat_pop)[tissue == "POP_site"]
nonpop_cols <- colnames(expr_mat_pop)[tissue == "NonPOP_site"]
pop_ind     <- individuals[tissue == "POP_site"]
nonpop_ind  <- individuals[tissue == "NonPOP_site"]
stopifnot(setequal(pop_ind, nonpop_ind))
pop_cols <- pop_cols[match(nonpop_ind, pop_ind)]  # align by patient

D_pop <- expr_mat_pop[, pop_cols] - expr_mat_pop[, nonpop_cols]  # + = up at prolapse site
valid_rows <- stats::complete.cases(D_pop) & (apply(D_pop, 1, sd) > 0)
D_pop <- D_pop[valid_rows, , drop = FALSE]
n_pat <- ncol(D_pop)
cat("Paired difference matrix:", nrow(D_pop), "genes x", n_pat, "patients\n")

paired_t <- function(D) { m <- rowMeans(D); s <- apply(D, 1, sd); m / (s / sqrt(ncol(D))) }
t_obs_pop <- paired_t(D_pop)
ord <- order(-t_obs_pop)
ranked_genes_pop <- rownames(D_pop)[ord]
ranked_scores_pop <- t_obs_pop[ord]
N_pop <- length(ranked_genes_pop)

ann_pop <- suppressWarnings(select(org.Hs.eg.db, keys = ranked_genes_pop, keytype = "SYMBOL", columns = "PATH"))
ann_pop <- ann_pop[!is.na(ann_pop$PATH), ]
gs_sizes_pop <- table(ann_pop$PATH)
valid_paths_pop <- names(gs_sizes_pop)[gs_sizes_pop >= 5 & gs_sizes_pop <= 200]
gene_sets_pop <- split(ann_pop$SYMBOL[ann_pop$PATH %in% valid_paths_pop], ann_pop$PATH[ann_pop$PATH %in% valid_paths_pop])
cat("KEGG pathways tested (5-200 members, offline org.Hs.eg.db):", length(gene_sets_pop), "\n")

calc_es <- function(hit_idx, scores_abs, N) {
  Nh <- length(hit_idx); Nm <- N - Nh
  if (Nm <= 0 || Nh == 0) return(NA)
  sum_hit <- sum(scores_abs[hit_idx])
  step <- rep(-1 / Nm, N); step[hit_idx] <- scores_abs[hit_idx] / sum_hit
  running <- cumsum(step); running[which.max(abs(running))]
}

hit_idx_pop <- lapply(gene_sets_pop, function(g) which(ranked_genes_pop %in% g))
hit_idx_pop <- hit_idx_pop[sapply(hit_idx_pop, length) >= 3]
cat("Pathways with >=3 genes in the ranked list:", length(hit_idx_pop), "\n")
es_obs_pop <- sapply(hit_idx_pop, calc_es, scores_abs = abs(ranked_scores_pop), N = N_pop)

n_perm <- 500
cat("Running", n_perm, "phenotype (sign-flip) permutations...\n")
t0 <- Sys.time()
gene_sets_syms_pop <- lapply(hit_idx_pop, function(idx) ranked_genes_pop[idx])
perm_es_pop <- matrix(NA_real_, nrow = n_perm, ncol = length(hit_idx_pop))
for (i in seq_len(n_perm)) {
  signs <- sample(c(-1, 1), n_pat, replace = TRUE)
  t_perm <- paired_t(sweep(D_pop, 2, signs, `*`))
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
write.csv(gsea_pop, "results/03_GSEA_classic_POP_KEGG.csv", row.names = FALSE)

cat("=== GSEA classic, POP (KEGG) ===\n")
cat("Pathways significant at FDR<0.25:", sum(gsea_pop$p.adjust < 0.25, na.rm = TRUE), "of", nrow(gsea_pop), "\n")
cat("Pathways significant at FDR<0.05:", sum(gsea_pop$p.adjust < 0.05, na.rm = TRUE), "\n\n")
print(head(gsea_pop[, c("PathwayName","Nh","NES","pvalue","p.adjust")], 10))
cat("\n")


## =============================================================================
## STEP 4: GSEA for SUI - PRERANKED (gene-set-label permutation)
## =============================================================================
cat("================ STEP 4: GSEA for SUI - PRERANKED method ================\n\n")
cat("Method: preranked GSEA, gene-set-label permutation. Ranking = REAL\n")
cat("moderated t-statistic from limma on the per-sample intensities (3 SUI vs\n")
cat("3 Ctrl, both sheets combined), not just the authors' summary p-value.\n")
cat("Classic phenotype permutation is NOT used here: with only 3 vs 3 samples\n")
cat("there are only choose(6,3)=20 possible relabelings, capping resolution at\n")
cat("p=1/20=0.05 and making FDR<0.05 structurally unreachable regardless of the\n")
cat("real signal. Preranked (gene-set) permutation avoids that ceiling because\n")
cat("it permutes over pathways, not over the 6 samples.\n\n")

set.seed(2020)
sui_mat <- as.matrix(sui_full[, c("Sui1","Sui2","Sui3","Ctrl1","Ctrl2","Ctrl3")])
rownames(sui_mat) <- sui_full$GeneSymbol
sui_mat <- avereps(sui_mat, ID = rownames(sui_mat))
cat("Per-sample matrix (probes collapsed by gene):", nrow(sui_mat), "genes x", ncol(sui_mat), "samples\n")

group_sui <- factor(c("SUI","SUI","SUI","Ctrl","Ctrl","Ctrl"), levels = c("Ctrl","SUI"))
design_sui <- model.matrix(~group_sui)
fit_sui <- eBayes(lmFit(sui_mat, design_sui))
t_obs_sui <- fit_sui$t[, 2]  # positive = up in SUI
ord_s <- order(-t_obs_sui)
ranked_genes_sui <- names(t_obs_sui)[ord_s]
ranked_scores_sui <- t_obs_sui[ord_s]
N_sui <- length(ranked_genes_sui)
cat("Ranked list (moderated t, positive = up in SUI):", N_sui, "genes\n\n")

ann_sui <- suppressWarnings(select(org.Hs.eg.db, keys = ranked_genes_sui, keytype = "SYMBOL", columns = "PATH"))
ann_sui <- ann_sui[!is.na(ann_sui$PATH), ]
gs_sizes_sui <- table(ann_sui$PATH)
valid_paths_sui <- names(gs_sizes_sui)[gs_sizes_sui >= 5 & gs_sizes_sui <= 200]
gene_sets_sui <- split(ann_sui$SYMBOL[ann_sui$PATH %in% valid_paths_sui], ann_sui$PATH[ann_sui$PATH %in% valid_paths_sui])
cat("KEGG pathways tested (5-200 members, offline org.Hs.eg.db):", length(gene_sets_sui), "\n")

hit_idx_sui <- lapply(gene_sets_sui, function(g) which(ranked_genes_sui %in% g))
hit_idx_sui <- hit_idx_sui[sapply(hit_idx_sui, length) >= 3]
cat("Pathways with >=3 genes in the ranked list:", length(hit_idx_sui), "\n")
es_obs_sui <- sapply(hit_idx_sui, calc_es, scores_abs = abs(ranked_scores_sui), N = N_sui)

n_perm2 <- 1000
cat("Running", n_perm2, "gene-set-label permutations...\n")
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
write.csv(gsea_sui, "results/04_GSEA_preranked_SUI_KEGG.csv", row.names = FALSE)

cat("=== GSEA preranked, SUI (KEGG) ===\n")
cat("Pathways significant at FDR<0.25:", sum(gsea_sui$p.adjust < 0.25, na.rm = TRUE), "of", nrow(gsea_sui), "\n")
cat("Pathways significant at FDR<0.05:", sum(gsea_sui$p.adjust < 0.05, na.rm = TRUE), "\n\n")
print(head(gsea_sui[, c("PathwayName","Nh","NES","pvalue","p.adjust")], 10))
cat("\n")


## =============================================================================
## STEP 5: Shared pathways (POP GSEA x SUI GSEA) + gene-level direction table
## =============================================================================
cat("================ STEP 5: Shared pathways POP x SUI ================\n\n")

report_shared <- function(fdr_cut) {
  pop_sig <- subset(gsea_pop, p.adjust < fdr_cut)
  sui_sig <- subset(gsea_sui, p.adjust < fdr_cut)
  shared_ids <- intersect(pop_sig$PATH, sui_sig$PATH)
  cat("--- FDR <", fdr_cut, ": POP significant =", nrow(pop_sig),
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
  if (n_na > 0) cat(" (", n_na, "excluded - NES undefined on at least one side; this is a",
                     "per-pathway normalization edge case (no permuted ES fell on the same",
                     "side as the observed one) and is distinct from the SUI sample-size",
                     "resolution limit discussed above)")
  cat("\n\n")
  out
}

shared_025 <- report_shared(0.25)
shared_005 <- report_shared(0.05)
if (nrow(shared_025) > 0) write.csv(shared_025, "results/05_shared_pathways_FDR025.csv", row.names = FALSE)
if (nrow(shared_005) > 0) write.csv(shared_005, "results/05_shared_pathways_FDR005.csv", row.names = FALSE)

# Pick the set of shared pathways to build the gene-level direction table for:
# prefer the FDR<0.05 list; fall back to FDR<0.25 if that one is empty.
shared_for_table <- if (nrow(shared_005) > 0) shared_005 else shared_025

cat("\n--- Gene-level direction table for the shared pathways ---\n")
if (nrow(shared_for_table) == 0) {
  cat("No shared pathway at either threshold - no gene-level table to build.\n\n")
  gene_dir_table <- data.frame()
} else {
  # POP and SUI logFC/direction per gene, keyed by symbol
  pop_lookup <- setNames(pop_full$logFC, pop_full$Gene)
  sui_lookup <- setNames(sui_full$logFC, sui_full$GeneSymbol)
  # "Tested" means the gene was present in that GSEA's ranked list (i.e. on
  # the array/table for that disease), not that it was individually
  # significant - GSEA looks at whole pathways, not single-gene cutoffs.
  pop_tested <- ranked_genes_pop
  sui_tested <- ranked_genes_sui

  rows <- lapply(seq_len(nrow(shared_for_table)), function(i) {
    pid <- shared_for_table$PATH[i]
    pname <- shared_for_table$PathwayName[i]
    # Full KEGG pathway membership (not restricted to genes tested in only
    # one of the two ranked lists), then keep genes tested in EITHER disease.
    full_members <- suppressWarnings(select(org.Hs.eg.db, keys = pid, keytype = "PATH", columns = "SYMBOL")$SYMBOL)
    members <- intersect(unique(full_members), union(pop_tested, sui_tested))
    data.frame(
      Pathway = pname, KEGG_ID = pid, Gene = members,
      logFC_POP = unname(pop_lookup[members]), logFC_SUI = unname(sui_lookup[members]),
      stringsAsFactors = FALSE)
  })
  gene_dir_table <- do.call(rbind, rows)
  gene_dir_table$Direction_POP <- ifelse(is.na(gene_dir_table$logFC_POP), "not tested",
                                          ifelse(gene_dir_table$logFC_POP > 0, "up", "down"))
  gene_dir_table$Direction_SUI <- ifelse(is.na(gene_dir_table$logFC_SUI), "not tested",
                                          ifelse(gene_dir_table$logFC_SUI > 0, "up", "down"))
  gene_dir_table$Concordant <- with(gene_dir_table,
                                     Direction_POP != "not tested" & Direction_SUI != "not tested" &
                                       Direction_POP == Direction_SUI)
  gene_dir_table <- gene_dir_table[order(gene_dir_table$Pathway, -gene_dir_table$Concordant), ]
  write.csv(gene_dir_table, "results/05_shared_pathways_gene_direction_table.csv", row.names = FALSE)

  both_tested <- subset(gene_dir_table, Direction_POP != "not tested" & Direction_SUI != "not tested")
  cat("Genes in shared pathways tested in both diseases:", nrow(both_tested),
      "| concordant direction:", sum(both_tested$Concordant),
      "(", round(100 * mean(both_tested$Concordant), 1), "% )\n\n")
  print(head(gene_dir_table, 20))
}
cat("\n")


## =============================================================================
## STEP 6: GRAPHICS - all figures saved to figures/
## =============================================================================
cat("================ STEP 6: Graphics ================\n\n")

# --- 6a. Volcano plots (POP and SUI) ----------------------------------------
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

v_pop <- make_volcano(pop_full, "logFC", "adj.P.Val", "Gene",
                       "Volcano plot - POP (GSE53868)\nDEG: prolapse site vs non-prolapse site (FDR<0.05, |log2FC|>1)")
ggsave("figures/01_volcano_POP.png", v_pop, width = 9, height = 6.5, dpi = 300)
cat("Saved: figures/01_volcano_POP.png\n")

v_sui <- make_volcano(sui_full, "logFC", "PValue", "GeneSymbol",
                       "Volcano plot - SUI (Wei 2020)\nDEG: SUI vs continent controls (P<0.05, |log2FC|>1, authors' criterion)")
ggsave("figures/02_volcano_SUI.png", v_sui, width = 9.5, height = 6.5, dpi = 300)
cat("Saved: figures/02_volcano_SUI.png\n\n")

# --- 6b. Heatmaps of top DEGs --------------------------------------------
ord_cols_pop <- order(coldata_pop$individual, coldata_pop$tissue)
ann_col_pop <- data.frame(Site = coldata_pop$tissue, Patient = coldata_pop$individual,
                           row.names = colnames(expr_pop))
top_pop_genes <- head(pop_deg$Gene, 40)
mat_pop_top <- expr_mat_pop[rownames(expr_mat_pop) %in% top_pop_genes, , drop = FALSE]
png("figures/03_heatmap_DEG_POP.png", width = 2600, height = 3200, res = 300)
pheatmap(mat_pop_top[, ord_cols_pop], scale = "row",
         annotation_col = ann_col_pop[ord_cols_pop, , drop = FALSE],
         cluster_cols = FALSE, cluster_rows = TRUE,
         main = "Top DEG genes - POP (GSE53868)\n(row z-score, paired patients)",
         fontsize = 7, show_colnames = FALSE)
dev.off()
cat("Saved: figures/03_heatmap_DEG_POP.png\n")

top_sui_genes <- head(sui_deg$GeneSymbol, 40)
mat_sui_top <- sui_mat[rownames(sui_mat) %in% top_sui_genes, , drop = FALSE]
ann_col_sui <- data.frame(Group = c("SUI","SUI","SUI","Ctrl","Ctrl","Ctrl"),
                           row.names = colnames(sui_mat))
png("figures/04_heatmap_DEG_SUI.png", width = 2200, height = 3200, res = 300)
pheatmap(mat_sui_top, scale = "row", annotation_col = ann_col_sui,
         cluster_cols = FALSE, cluster_rows = TRUE,
         main = "Top DEG genes - SUI (Wei 2020)\n(row z-score, 3 SUI vs 3 controls)",
         fontsize = 7)
dev.off()
cat("Saved: figures/04_heatmap_DEG_SUI.png\n\n")

# --- 6c. GSEA barplots (POP and SUI) ----------------------------------------
make_gsea_barplot <- function(gsea_res, title, n_top = 15) {
  # NES can be NA for a handful of pathways where all permuted ES values fell
  # on one side of zero (permutation-floor edge case) - drop those before
  # picking the top N, or the barplot shows empty rows with a label but no bar.
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

p_bar_pop <- make_gsea_barplot(gsea_pop, "GSEA (classic) - top KEGG pathways in POP")
ggsave("figures/05_GSEA_barplot_POP.png", p_bar_pop, width = 10, height = 6, dpi = 300)
cat("Saved: figures/05_GSEA_barplot_POP.png\n")

p_bar_sui <- make_gsea_barplot(gsea_sui, "GSEA (preranked) - top KEGG pathways in SUI")
ggsave("figures/06_GSEA_barplot_SUI.png", p_bar_sui, width = 10, height = 6, dpi = 300)
cat("Saved: figures/06_GSEA_barplot_SUI.png\n\n")

# --- 6d. Shared-pathway NES comparison plot ---------------------------------
compare_pool <- if (nrow(shared_025) > 0) shared_025 else data.frame()
if (nrow(compare_pool) > 0) {
  compare_pool$Label <- reorder(compare_pool$PathwayName, compare_pool$NES_POP)
  p_shared <- ggplot(compare_pool, aes(x = NES_POP, y = NES_SUI, label = PathwayName)) +
    geom_hline(yintercept = 0, color = "grey70") + geom_vline(xintercept = 0, color = "grey70") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
    geom_point(aes(color = Same_direction), size = 3) +
    scale_color_manual(values = c("TRUE" = "#2166AC", "FALSE" = "#B2182B"),
                        labels = c("TRUE" = "Same direction", "FALSE" = "Opposite direction")) +
    geom_text_repel(size = 3, max.overlaps = 30) +
    labs(title = "Shared pathways (FDR<0.25): NES in POP vs NES in SUI",
         x = "NES - POP (classic GSEA)", y = "NES - SUI (preranked GSEA)", color = NULL) +
    theme_bw() + theme(legend.position = "top")
  ggsave("figures/07_shared_pathways_NES_comparison.png", p_shared, width = 9, height = 7, dpi = 300)
  cat("Saved: figures/07_shared_pathways_NES_comparison.png\n\n")
} else {
  cat("SKIPPED: figures/07_shared_pathways_NES_comparison.png (no shared pathway at FDR<0.25)\n\n")
}

# --- 6e. Classic GSEA running-score enrichment plot (top pathway, POP) -----
gsea_running_score <- function(hit_idx, scores_abs, N) {
  Nh <- length(hit_idx); Nm <- N - Nh
  sum_hit <- sum(scores_abs[hit_idx])
  step <- rep(-1 / Nm, N); step[hit_idx] <- scores_abs[hit_idx] / sum_hit
  cumsum(step)
}

top_path_pop <- gsea_pop$PATH[1]
running_pop <- gsea_running_score(hit_idx_pop[[top_path_pop]], abs(ranked_scores_pop), N_pop)
df_run <- data.frame(rank = seq_len(N_pop), es = running_pop)
df_hits <- data.frame(rank = hit_idx_pop[[top_path_pop]])
p_run <- ggplot(df_run, aes(x = rank, y = es)) +
  geom_line(color = "#2166AC", linewidth = 1) +
  geom_hline(yintercept = 0, color = "grey40") +
  geom_rug(data = df_hits, aes(x = rank), inherit.aes = FALSE, sides = "b", color = "black") +
  labs(title = paste0("GSEA enrichment plot (classic, POP) - top pathway:\n", gsea_pop$PathwayName[1],
                       "\nNES=", round(gsea_pop$NES[1], 2), ", p=", signif(gsea_pop$pvalue[1], 3),
                       ", FDR=", signif(gsea_pop$p.adjust[1], 3)),
       x = "Rank in ranked gene list (POP, high to low t-statistic)",
       y = "Running enrichment score") +
  theme_bw()
ggsave("figures/08_GSEA_enrichment_plot_top_POP.png", p_run, width = 9, height = 5.5, dpi = 300)
cat("Saved: figures/08_GSEA_enrichment_plot_top_POP.png (top pathway:", gsea_pop$PathwayName[1], ")\n\n")

# --- 6f. Gene-direction heatmap for shared-pathway genes --------------------
# Fill color = direction (up/down) so the concordance pattern is visible
# regardless of magnitude - POP fold-changes are much smaller than SUI's, so
# a single shared continuous color scale (fill = logFC) makes the POP column
# look uniformly pale and hides the pattern. The actual logFC values are
# still shown as text in each cell.
if (nrow(gene_dir_table) > 0) {
  both_tested <- subset(gene_dir_table, Direction_POP != "not tested" & Direction_SUI != "not tested")
  if (nrow(both_tested) >= 2) {
    fc_mat <- as.matrix(both_tested[, c("logFC_POP", "logFC_SUI")])
    dir_mat <- ifelse(fc_mat > 0, 1, -1)
    rownames(dir_mat) <- rownames(fc_mat) <- make.unique(both_tested$Gene)
    ord <- order(both_tested$Pathway, -both_tested$Concordant, -abs(rowSums(fc_mat, na.rm = TRUE)))
    dir_mat <- dir_mat[ord, , drop = FALSE]; fc_mat <- fc_mat[ord, , drop = FALSE]
    dir_mat <- head(dir_mat, 50); fc_mat <- head(fc_mat, 50)
    ann_row <- data.frame(Pathway = both_tested$Pathway[ord][seq_len(nrow(dir_mat))],
                           row.names = rownames(dir_mat))
    colnames(dir_mat) <- colnames(fc_mat) <- c("POP", "SUI")
    png("figures/09_gene_direction_heatmap_shared_pathways.png",
        width = 2600, height = max(1800, 45 * nrow(dir_mat)), res = 300)
    pheatmap(dir_mat, cluster_cols = FALSE, cluster_rows = FALSE,
             annotation_row = ann_row,
             color = c("#2166AC", "#B2182B"), breaks = c(-2, 0, 2),
             legend = FALSE,
             display_numbers = matrix(sprintf("%.2f", fc_mat), nrow(fc_mat)),
             number_color = "white", fontsize_number = 7,
             main = "Gene direction, shared pathways (POP vs SUI)\n(fill = up/down; number = log2 Fold Change)",
             fontsize = 7, fontsize_row = 6)
    dev.off()
    cat("Saved: figures/09_gene_direction_heatmap_shared_pathways.png\n\n")
  } else {
    cat("SKIPPED: figures/09_gene_direction_heatmap_shared_pathways.png (fewer than 2 genes tested in both)\n\n")
  }
} else {
  cat("SKIPPED: figures/09_gene_direction_heatmap_shared_pathways.png (no shared-pathway gene table)\n\n")
}

cat("\n=== ALL FIGURES SAVED TO figures/ ===\n")
cat("01_volcano_POP.png, 02_volcano_SUI.png\n")
cat("03_heatmap_DEG_POP.png, 04_heatmap_DEG_SUI.png\n")
cat("05_GSEA_barplot_POP.png, 06_GSEA_barplot_SUI.png\n")
cat("07_shared_pathways_NES_comparison.png\n")
cat("08_GSEA_enrichment_plot_top_POP.png\n")
cat("09_gene_direction_heatmap_shared_pathways.png\n")
cat("\n=== DONE ===\n")
