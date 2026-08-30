# =============================================================================
# Preranked GSEA on POP (GSE53868 microarray) vs. pathway-level analysis of
# SUI (Chen et al. 2006, 58 curated genes) - looking for CONVERGENT pathways.
#
# WHY TWO DIFFERENT METHODS ON THE TWO SIDES (data-availability constraint,
# not an inconsistency):
#   - POP (GSE53868): the FULL ranked transcriptome is available (all
#     ~31,000 genes, paired design, 12 patients) -> a real preranked GSEA
#     is run below.
#   - SUI (Chen et al. 2006): the original paper only reports a curated
#     table of 58 candidate genes ALREADY selected as significant by the
#     original authors, not the full ~33,000-probe ranked output of their
#     array. GSEA requires a full ranked list to build its null
#     distribution; it cannot be run on a 58-gene pre-filtered list. This
#     was explained before (R/GSEA_POP_exemplo.R) and remains true. The
#     closest defensible pathway-level analysis for a fixed candidate list
#     is over-representation analysis (ORA, hypergeometric test), already
#     computed in results/KEGG_Chen2006_58genes.csv.
#   We compare the POP GSEA results against the Chen 2006 ORA results to
#   find CONVERGENT pathways significant in both.
#
# GSEA IMPLEMENTATION NOTE: no internet access is available here to install
# fgsea/clusterProfiler/msigdbr (CRAN and Bioconductor both unreachable,
# confirmed again below) and MSigDB (Hallmark) cannot be downloaded. We
# implement preranked GSEA from scratch (Subramanian et al. 2005, PNAS):
# weighted running-sum enrichment score (ES) against KEGG pathway gene sets
# from org.Hs.eg.db (offline), and, critically, PHENOTYPE PERMUTATION
# (not gene-set-label permutation) for the null distribution.
#
# A FIRST VERSION of this script used gene-set-label permutation (shuffling
# which genes count as "hits" on a FIXED ranking) - the simpler fallback
# GSEA-Preranked uses when no per-sample matrix is available. That first
# run flagged 140 of 218 pathways (64%) as significant at FDR<0.25, an
# implausibly high hit rate. This is a KNOWN, documented weakness of
# gene-set-label permutation: real genes within a biological pathway are
# co-expressed/correlated in actual tissue data, so shuffling gene labels
# underestimates the true null variance and inflates false positives
# (Subramanian et al. 2005 explicitly warn about this and recommend
# phenotype permutation whenever a per-sample matrix is available - which
# we have here, since GSE53868 is a paired design with raw expression
# values, not just a pre-computed ranked list). This corrected version re-
# ranks the ENTIRE gene list under each permutation of the sample labels
# (sign-flipping the paired POP-site/non-POP-site difference per patient,
# the standard permutation scheme for a paired/matched design, equivalent
# to phenotype permutation for a paired t-test) and is the statistically
# valid approach - flagged as a self-correction, not swept under the rug.
# =============================================================================

suppressMessages(library(org.Hs.eg.db))
set.seed(42)

## -----------------------------------------------------------------------
## 1) Rebuild the PAIRED per-patient difference matrix for GSE53868 (same
##    12 patients, POP-site vs non-POP-site, as in
##    R/analysis_GSE53868_limma_x_Chen2006.R) - needed for phenotype
##    (sign-flip) permutation, not just the final limma table.
## -----------------------------------------------------------------------
raw_lines <- readLines("data/GSE53868_series_matrix.txt")
start_row <- grep("^!series_matrix_table_begin", raw_lines) + 1
end_row   <- grep("^!series_matrix_table_end", raw_lines) - 1
expr <- read.delim("data/GSE53868_series_matrix.txt", skip = start_row - 1,
                    nrows = end_row - start_row, header = TRUE,
                    row.names = 1, check.names = FALSE, quote = "\"")

sample_title_line <- raw_lines[grep("^!Sample_title", raw_lines)]
sample_titles <- gsub('"', "", strsplit(sample_title_line, "\t")[[1]][-1])
individual_line <- raw_lines[grep("^!Sample_characteristics_ch1.*individual:", raw_lines)][1]
individuals <- gsub("individual: ", "", gsub('"', "", strsplit(individual_line, "\t")[[1]][-1]))
tissue <- ifelse(grepl("\\(POP site\\)", sample_titles), "POP_site", "NonPOP_site")

expr <- as.matrix(expr)
pop_cols <- colnames(expr)[tissue == "POP_site"]
nonpop_cols <- colnames(expr)[tissue == "NonPOP_site"]
pop_ind <- individuals[tissue == "POP_site"]
nonpop_ind <- individuals[tissue == "NonPOP_site"]
stopifnot(setequal(pop_ind, nonpop_ind))
pop_cols <- pop_cols[match(nonpop_ind, pop_ind)]  # align by individual

D <- expr[, pop_cols] - expr[, nonpop_cols]  # genes x 12 patients, + = up in POP site

# Drop rows with any missing value or zero variance across patients (a
# paired t-statistic is undefined for these - 59 such probes found here,
# e.g. constant/flagged probes) so every permutation operates on the exact
# same, fully-defined gene universe (avoids NA/length-mismatch bugs from
# ranking a vector that contains undefined scores).
valid_rows <- stats::complete.cases(D) & (apply(D, 1, sd) > 0)
cat("Dropping", sum(!valid_rows), "probes with missing values or zero variance across patients\n")
D <- D[valid_rows, , drop = FALSE]
n_pat <- ncol(D)
cat("Paired difference matrix D:", nrow(D), "genes x", n_pat, "patients\n")

paired_t <- function(D) {
  m <- rowMeans(D)
  s <- apply(D, 1, sd)
  n <- ncol(D)
  m / (s / sqrt(n))
}

t_obs <- paired_t(D)
ord <- order(-t_obs)
ranked_genes <- rownames(D)[ord]
ranked_scores <- t_obs[ord]
N <- length(ranked_genes)
cat("Ranked list (paired t-statistic, positive = up in POP site):", N, "genes\n")
cat("Top 5:", paste(head(ranked_genes, 5), collapse = ", "), "\n")
cat("Bottom 5:", paste(tail(ranked_genes, 5), collapse = ", "), "\n\n")

## -----------------------------------------------------------------------
## 2) KEGG pathway gene sets (offline, org.Hs.eg.db), size 5-200.
## -----------------------------------------------------------------------
ann <- suppressWarnings(select(org.Hs.eg.db, keys = ranked_genes, keytype = "SYMBOL", columns = "PATH"))
ann <- ann[!is.na(ann$PATH), ]
gs_sizes <- table(ann$PATH)
valid_paths <- names(gs_sizes)[gs_sizes >= 5 & gs_sizes <= 200]
gene_sets <- split(ann$SYMBOL[ann$PATH %in% valid_paths], ann$PATH[ann$PATH %in% valid_paths])
cat("KEGG pathways tested (5-200 genes, offline org.Hs.eg.db):", length(gene_sets), "\n\n")

## -----------------------------------------------------------------------
## 3) Enrichment score (weighted, p=1) on a given ranking.
## -----------------------------------------------------------------------
calc_es <- function(hit_idx, scores_abs, N) {
  Nh <- length(hit_idx); Nm <- N - Nh
  sum_hit <- sum(scores_abs[hit_idx])
  step <- rep(-1 / Nm, N)
  step[hit_idx] <- scores_abs[hit_idx] / sum_hit
  running <- cumsum(step)
  running[which.max(abs(running))]
}

gene_pos <- match(ranked_genes, ranked_genes)  # identity, kept for clarity
hit_idx_list <- lapply(gene_sets, function(g) which(ranked_genes %in% g))
hit_idx_list <- hit_idx_list[sapply(hit_idx_list, length) >= 3]
cat("Pathways with >=3 genes in ranked list:", length(hit_idx_list), "\n")

es_obs <- sapply(hit_idx_list, calc_es, scores_abs = abs(ranked_scores), N = N)

## -----------------------------------------------------------------------
## 4) PHENOTYPE permutation: sign-flip each patient's paired difference
##    (standard permutation null for a paired design; equivalent to
##    permuting which sample is "POP site" vs "non-POP site" within each
##    matched pair), recompute the paired t-statistic for ALL genes, and
##    RE-RANK the whole gene list, for every permutation. This preserves
##    the real gene-gene correlation structure under the null, unlike
##    gene-set-label permutation.
## -----------------------------------------------------------------------
n_perm <- 500
cat("Running", n_perm, "phenotype (sign-flip) permutations across", length(hit_idx_list), "pathways...\n")
t0 <- Sys.time()

# pathway gene symbols (fixed across permutations) for fast name-based lookup
gene_sets_syms <- lapply(hit_idx_list, function(idx) ranked_genes[idx])

perm_es_mat <- matrix(NA_real_, nrow = n_perm, ncol = length(hit_idx_list))
colnames(perm_es_mat) <- names(hit_idx_list)

for (i in seq_len(n_perm)) {
  signs <- sample(c(-1, 1), n_pat, replace = TRUE)
  D_perm <- sweep(D, 2, signs, `*`)
  t_perm <- paired_t(D_perm)                       # named vector, names = rownames(D)
  rank_of_gene <- rank(-t_perm, ties.method = "first")  # 1 = most up under this permutation
  scores_abs_sorted <- sort(abs(t_perm), decreasing = TRUE)  # scores in rank order
  for (j in seq_along(gene_sets_syms)) {
    hidx <- rank_of_gene[gene_sets_syms[[j]]]       # O(Nh) name lookup, no O(N) scan
    if (length(hidx) >= 3) perm_es_mat[i, j] <- calc_es(hidx, scores_abs_sorted, N)
  }
}
cat("Done in", round(difftime(Sys.time(), t0, units = "secs"), 1), "seconds\n\n")

pval <- numeric(length(hit_idx_list))
nes <- numeric(length(hit_idx_list))
for (j in seq_along(hit_idx_list)) {
  pe <- perm_es_mat[, j]
  pe <- pe[!is.na(pe)]
  if (es_obs[j] >= 0) {
    pval[j] <- (sum(pe >= es_obs[j]) + 1) / (length(pe) + 1)
    base <- mean(pe[pe >= 0]); if (is.nan(base) || base == 0) base <- NA
  } else {
    pval[j] <- (sum(pe <= es_obs[j]) + 1) / (length(pe) + 1)
    base <- mean(abs(pe[pe < 0])); if (is.nan(base) || base == 0) base <- NA
  }
  nes[j] <- es_obs[j] / base
}

leading_edge <- sapply(hit_idx_list, function(idx) paste(ranked_genes[idx], collapse = "/"))
gsea_res <- data.frame(PATH = names(hit_idx_list), Nh = sapply(hit_idx_list, length),
                        ES = es_obs, NES = nes, pvalue = pval, leadingEdge = leading_edge)
gsea_res$p.adjust <- p.adjust(gsea_res$pvalue, "BH")
gsea_res <- gsea_res[order(gsea_res$pvalue), ]

write.csv(gsea_res, "results/GSEA_POP_GSE53868_KEGG.csv", row.names = FALSE)
cat("POP GSEA (phenotype permutation) - pathways FDR<0.25 (standard GSEA screening threshold):",
    sum(gsea_res$p.adjust < 0.25, na.rm = TRUE), "of", nrow(gsea_res), "tested\n")
cat("POP GSEA - pathways FDR<0.05:", sum(gsea_res$p.adjust < 0.05, na.rm = TRUE), "\n")
print(head(gsea_res[, c("PATH", "Nh", "ES", "NES", "pvalue", "p.adjust")], 15))

## -----------------------------------------------------------------------
## 5) Chen 2006 (SUI) side: KEGG ORA on the 58-gene panel (real GSEA is not
##    possible - see header note). RE-COMPUTED HERE rather than reusing
##    results/KEGG_Chen2006_58genes.csv, because building this script
##    surfaced a real universe-size bug in that earlier file (and in every
##    other "standalone panel" ORA run this session): the background
##    universe was built as keys(org.Hs.eg.db, keytype="SYMBOL"), which in
##    this package snapshot returns 191,076 unique symbols - NCBI's full
##    Entrez Gene catalogue for human (pseudogenes, pending loci, etc
##    included), not the ~5,900 genes that actually carry a KEGG pathway
##    annotation. Using the bloated, mostly-unannotated 191,076 as N in the
##    hypergeometric test inflates significance (the earlier file flagged
##    74 of 76 KEGG pathways as significant at FDR<0.05, which is
##    implausibly high and was the tell that something was wrong). The fix,
##    standard ORA practice, is to restrict the universe to genes that
##    could in principle be annotated to the collection being tested (here:
##    genes with >=1 KEGG PATH entry, ~5,868 of them) - done below.
##    NOTE: this bug likely also affected other "standalone panel" ORA
##    results delivered earlier this session (Chen2003, Poelmans, Tong2010
##    panels tested alone) - flagged to the user, not silently corrected
##    everywhere, since re-running all of them is a separate task.
## -----------------------------------------------------------------------
chen58 <- read.csv("data/chen2006_79genes.csv")
chen58_genes <- unique(chen58$GeneSymbol)

kegg_universe <- {
  ann_full <- suppressWarnings(select(org.Hs.eg.db, keys = keys(org.Hs.eg.db, keytype = "SYMBOL"),
                                       keytype = "SYMBOL", columns = "PATH"))
  unique(ann_full$SYMBOL[!is.na(ann_full$PATH)])
}
cat("\nCorrected KEGG-annotatable universe:", length(kegg_universe), "genes",
    "(vs. 191,076 used in the earlier, inflated run)\n")

run_ora_kegg <- function(hit_genes, universe_genes) {
  hit_genes <- unique(intersect(hit_genes, universe_genes))
  ann <- suppressWarnings(select(org.Hs.eg.db, keys = universe_genes, keytype = "SYMBOL", columns = "PATH"))
  ann <- ann[!is.na(ann$PATH), ]
  gs <- table(ann$PATH); valid <- names(gs)[gs >= 2 & gs <= 2000]
  ann <- ann[ann$PATH %in% valid, ]
  N <- length(universe_genes); n <- length(hit_genes)
  res <- lapply(valid, function(t) {
    g <- unique(ann$SYMBOL[ann$PATH == t]); M <- length(g)
    h <- intersect(g, hit_genes); k <- length(h)
    if (k == 0) return(NULL)
    p <- phyper(k - 1, M, N - M, n, lower.tail = FALSE)
    data.frame(TERM_ID = t, Count = k, M = M, N = N, n = n, pvalue = p, geneID = paste(h, collapse = "/"))
  })
  res <- do.call(rbind, res)
  res$p.adjust <- p.adjust(res$pvalue, "BH")
  res[order(res$pvalue), ]
}

chen_kegg <- run_ora_kegg(chen58_genes, kegg_universe)
write.csv(chen_kegg, "results/KEGG_Chen2006_58genes_CORRECTED_universe.csv", row.names = FALSE)
cat("SUI (Chen 2006) - KEGG ORA on the 58-gene panel, CORRECTED universe:",
    nrow(chen_kegg), "pathways tested,",
    sum(chen_kegg$p.adjust < 0.05, na.rm = TRUE), "significant at FDR<0.05",
    "(vs. 74 of 76 with the old, inflated universe)\n")

## -----------------------------------------------------------------------
## 6) Convergent pathways.
## -----------------------------------------------------------------------
pop_sig <- subset(gsea_res, p.adjust < 0.25)
chen_sig <- subset(chen_kegg, p.adjust < 0.05)
common_ids <- intersect(pop_sig$PATH, chen_sig$TERM_ID)
cat("\n=== Convergent KEGG pathways (POP GSEA FDR<0.25 AND Chen2006 ORA FDR<0.05) ===\n")
cat("Common pathway IDs:", length(common_ids), "\n")
if (length(common_ids) > 0) {
  convergent <- merge(pop_sig[pop_sig$PATH %in% common_ids, c("PATH","Nh","NES","pvalue","p.adjust","leadingEdge")],
                       chen_sig[chen_sig$TERM_ID %in% common_ids, c("TERM_ID","Count","pvalue","p.adjust","geneID")],
                       by.x = "PATH", by.y = "TERM_ID", suffixes = c("_POP_GSEA", "_Chen_ORA"))
  print(convergent)
  write.csv(convergent, "results/GSEA_convergent_pathways_POP_Chen2006.csv", row.names = FALSE)
} else {
  cat("None at these thresholds.\n")
}

cat("\n=== END OF SCRIPT ===\n")
