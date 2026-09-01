# =============================================================================
# Lower POP's threshold to |log2FC|>0.5 (from the stricter 1.0 used in the
# last GSEA-preranked run) and check KEGG pathway convergence via a
# DIFFERENT, better-powered route than small-panel preranked GSEA: run
# GO/KEGG ORA directly on the genes that are individually direction-
# concordant between each paper and POP (the same approach already used
# for Chen 2006, now completed for Chen 2003 and Tong 2010 too), using the
# REAL POP microarray gene list (31,072 genes) as background universe -
# not the tiny 60-90 gene panel used in the small-panel preranked GSEA,
# which has almost no statistical power for KEGG (see the previous
# script/README section for why).
#
# NOTE ON THE 0.5 THRESHOLD: this does NOT change POP's real preranked
# GSEA at all (GSEA always used the FULL ranked transcriptome, no
# threshold, in every version run so far) - it also does NOT change the
# concordant-gene list itself for any of the 3 papers (concordance is
# defined by matching DIRECTION only, with no significance threshold on
# either side - always was, see analysis_Chen2003_Poelmans_x_GSE53868.R,
# analysis_Tong2010_Chen2006_x_GSE53868_genes_e_vias.R). The 0.5 threshold
# only affects: (a) how many POP genes count as "individually significant"
# for reporting, and (b) which POP-standalone DEG list feeds the
# already-existing results/KEGG_GSE53868.csv (534 genes, |log2FC|>0.5,
# already computed earlier this session using the correct/real universe -
# reused here unchanged).
# =============================================================================

suppressMessages({
  library(org.Hs.eg.db)
  library(GO.db)
})

pop <- read.csv("results/GSE53868_limma_completo.csv")
universe_53868 <- pop$Gene

run_enrichment <- function(hit_genes, universe_genes, keytype_col, ont_filter = NULL, min_gs = 2, max_gs = 2000) {
  hit_genes <- unique(intersect(hit_genes, universe_genes)); universe_genes <- unique(universe_genes)
  if (length(hit_genes) == 0) return(NULL)
  ann <- suppressWarnings(select(org.Hs.eg.db, keys = universe_genes, keytype = "SYMBOL", columns = keytype_col))
  ann <- ann[!is.na(ann[[keytype_col]]), c("SYMBOL", keytype_col)]
  if (!is.null(ont_filter)) {
    ann2 <- suppressWarnings(select(org.Hs.eg.db, keys = universe_genes, keytype = "SYMBOL", columns = c("GO","ONTOLOGY")))
    keep <- unique(ann2$GO[ann2$ONTOLOGY == ont_filter & !is.na(ann2$GO)])
    ann <- ann[ann[[keytype_col]] %in% keep, ]
  }
  colnames(ann) <- c("SYMBOL", "TERM_ID"); ann <- unique(ann)
  gs <- table(ann$TERM_ID); valid <- names(gs)[gs >= min_gs & gs <= max_gs]; ann <- ann[ann$TERM_ID %in% valid, ]
  N <- length(universe_genes); n <- length(hit_genes)
  res <- lapply(valid, function(t) {
    g <- unique(ann$SYMBOL[ann$TERM_ID == t]); M <- length(g)
    h <- intersect(g, hit_genes); k <- length(h)
    if (k == 0) return(NULL)
    p <- phyper(k - 1, M, N - M, n, lower.tail = FALSE)
    data.frame(TERM_ID = t, Count = k, M = M, N = N, n = n, pvalue = p, geneID = paste(h, collapse = "/"))
  })
  res <- do.call(rbind, res); if (is.null(res)) return(NULL)
  res$p.adjust <- p.adjust(res$pvalue, "BH"); res[order(res$pvalue), ]
}

run_go_kegg_concordant <- function(csv_path, tested_col, concordant_col, tag) {
  cat("===", tag, "===\n")
  df <- read.csv(csv_path)
  testv <- df[df[[tested_col]] == TRUE, ]
  concordant <- unique(testv$GeneSymbol[testv[[concordant_col]] == TRUE])
  cat("Concordant genes:", length(concordant), "\n")

  go <- run_enrichment(concordant, universe_53868, "GO", ont_filter = "BP")
  if (!is.null(go)) {
    terms <- suppressMessages(select(GO.db, keys = go$TERM_ID, keytype = "GOID", columns = "TERM"))
    go <- merge(go, terms, by.x = "TERM_ID", by.y = "GOID"); go <- go[order(go$pvalue), ]
    cat("GO BP sig:", sum(go$p.adjust < 0.05, na.rm = TRUE), "of", nrow(go), "\n")
  }
  kegg <- run_enrichment(concordant, universe_53868, "PATH")
  if (!is.null(kegg)) cat("KEGG sig:", sum(kegg$p.adjust < 0.05, na.rm = TRUE), "of", nrow(kegg), "\n")
  cat("\n")
  list(go = go, kegg = kegg)
}

chen03 <- run_go_kegg_concordant("results/Chen2003_x_GSE53868.csv", "Testado_no_GSE53868", "Concordante", "Chen 2003 (concordant genes)")
tong10 <- run_go_kegg_concordant("results/Tong2010_x_GSE53868.csv", "Testado", "Concordante", "Tong 2010 (concordant genes)")
chen06_kegg <- read.csv("results/KEGG_Chen_x_GSE53868_concordantes_CORRECTED.csv")  # already done, reused

if (!is.null(chen03$go))   write.csv(chen03$go,   "results/GO_BP_Chen2003_concordant_GSE53868_CORRECTED.csv", row.names = FALSE)
if (!is.null(chen03$kegg)) write.csv(chen03$kegg, "results/KEGG_Chen2003_concordant_GSE53868_CORRECTED.csv", row.names = FALSE)
if (!is.null(tong10$go))   write.csv(tong10$go,   "results/GO_BP_Tong2010_concordant_GSE53868_CORRECTED.csv", row.names = FALSE)
if (!is.null(tong10$kegg)) write.csv(tong10$kegg, "results/KEGG_Tong2010_concordant_GSE53868_CORRECTED.csv", row.names = FALSE)

## -----------------------------------------------------------------------
## Check overlap with POP's own significant KEGG pathways - two POP
## reference lists, both already computed and reused unchanged: the
## ORA-based list (534 DEGs, |log2FC|>0.5, real universe) and the real
## preranked GSEA list (full transcriptome, phenotype permutation).
## -----------------------------------------------------------------------
# NORMALIZE all KEGG IDs to plain integers before comparing - org.Hs.eg.db
# returns zero-padded 5-char IDs (e.g. "04350") kept as-is in in-memory
# objects, but a CSV round-trip (write.csv then read.csv) silently strips
# the leading zero (R/pandas both auto-parse "04350" as the number 4350
# when reading a CSV). Comparing without normalizing first is a real bug
# that was caught here: an in-memory "04350" vs a CSV-read "4350" look
# different as strings even though they are the same pathway - always
# convert to integer (which is format-independent) before any ID matching.
pop_ora <- read.csv("results/KEGG_GSE53868.csv"); pop_ora$TERM_ID <- as.integer(pop_ora$TERM_ID)
pop_ora_sig <- subset(pop_ora, p.adjust < 0.05)
pop_gsea <- read.csv("results/GSEA_POP_GSE53868_KEGG.csv"); pop_gsea$PATH <- as.integer(pop_gsea$PATH)
pop_gsea_sig <- subset(pop_gsea, p.adjust < 0.25)

check_overlap <- function(kegg_tab, tag) {
  if (is.null(kegg_tab)) { cat(tag, ": no KEGG pathways\n"); return(invisible(NULL)) }
  kegg_tab$TERM_ID <- as.integer(kegg_tab$TERM_ID)
  sig <- subset(kegg_tab, p.adjust < 0.05)
  ov_ora <- intersect(sig$TERM_ID, pop_ora_sig$TERM_ID)
  ov_gsea <- intersect(sig$TERM_ID, pop_gsea_sig$PATH)
  cat(tag, "- sig KEGG:", nrow(sig), "| overlap w/ POP ORA (8):", length(ov_ora),
      paste(ov_ora, collapse=","), "| overlap w/ POP GSEA (94):", length(ov_gsea),
      paste(ov_gsea, collapse=","), "\n")
}

cat("\n=== Overlap with POP's own significant KEGG pathways ===\n")
check_overlap(chen03$kegg, "Chen 2003")
check_overlap(chen06_kegg, "Chen 2006")
check_overlap(tong10$kegg, "Tong 2010")

cat("\n=== END ===\n")
