# =============================================================================
# CORRECTION: re-run standalone-panel GO/KEGG ORA (Chen 2006, Chen 2003,
# Poelmans, Tong 2010) with a properly-sized background universe.
#
# BUG FOUND while building R/analysis_GSEA_POP_vs_Chen2006_convergent_pathways.R:
# every "panel alone" ORA call this session used
# keys(org.Hs.eg.db, keytype = "SYMBOL") as the background universe. In this
# package snapshot that returns 191,076 unique symbols - NCBI's full Entrez
# Gene catalogue for human (pseudogenes, uncharacterized loci, etc
# included), NOT the much smaller set of genes that can actually be
# annotated to a GO term or KEGG pathway. Using this inflated, mostly
# unannotated number as N in the hypergeometric test makes the population
# background proportion (M/N) look artificially small, so any real hit
# count looks artificially enriched - deflating p-values project-wide for
# every "panel alone" run (NOT the DEG-crossing runs, which used the real
# POP dataset gene list as universe_pop and are unaffected).
#
# FIX (standard ORA practice): restrict the universe to genes that could in
# principle receive an annotation in the collection being tested - here,
# genes with >=1 GO BP term (~20,700) or >=1 KEGG PATH entry (~5,868).
#
# This script re-runs GO BP + KEGG ORA for all four standalone gene panels
# used this session, and re-checks the two "pathways in common with POP"
# comparisons that depended on them (Tong 2010 vs POP, Chen 2006 vs POP).
# =============================================================================

suppressMessages({
  library(org.Hs.eg.db)
  library(GO.db)
})

## -----------------------------------------------------------------------
## Corrected universes (computed once)
## -----------------------------------------------------------------------
all_sym <- keys(org.Hs.eg.db, keytype = "SYMBOL")
ann_path_all <- suppressWarnings(select(org.Hs.eg.db, keys = all_sym, keytype = "SYMBOL", columns = "PATH"))
kegg_universe <- unique(ann_path_all$SYMBOL[!is.na(ann_path_all$PATH)])
ann_go_all <- suppressWarnings(select(org.Hs.eg.db, keys = all_sym, keytype = "SYMBOL", columns = c("GO","ONTOLOGY")))
go_bp_universe <- unique(ann_go_all$SYMBOL[!is.na(ann_go_all$GO) & ann_go_all$ONTOLOGY == "BP"])
cat("Corrected KEGG universe:", length(kegg_universe), "genes (old, wrong: 191,076)\n")
cat("Corrected GO-BP universe:", length(go_bp_universe), "genes (old, wrong: 191,076)\n\n")

run_ora <- function(hit_genes, universe_genes, keytype_col, ont_filter = NULL, min_gs = 2, max_gs = 2000) {
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

redo_panel <- function(genes, tag) {
  cat("---", tag, "(n=", length(genes), "genes) ---\n")
  go <- run_ora(genes, go_bp_universe, "GO", ont_filter = "BP")
  if (!is.null(go)) {
    terms <- suppressMessages(select(GO.db, keys = go$TERM_ID, keytype = "GOID", columns = "TERM"))
    go <- merge(go, terms, by.x = "TERM_ID", by.y = "GOID"); go <- go[order(go$pvalue), ]
    cat("GO BP: ", sum(go$p.adjust < 0.05, na.rm = TRUE), " sig (FDR<0.05) of ", nrow(go), " tested\n", sep = "")
  } else { cat("GO: no terms\n") }
  kegg <- run_ora(genes, kegg_universe, "PATH")
  if (!is.null(kegg)) cat("KEGG: ", sum(kegg$p.adjust < 0.05, na.rm = TRUE), " sig (FDR<0.05) of ", nrow(kegg), " tested\n\n", sep = "")
  else cat("KEGG: no pathways\n\n")
  list(go = go, kegg = kegg)
}

## -----------------------------------------------------------------------
## Re-run the four standalone panels
## -----------------------------------------------------------------------
chen06 <- unique(read.csv("data/chen2006_79genes.csv")$GeneSymbol)
chen03 <- unique(read.csv("data/chen2003_90genes.csv")$GeneSymbol)
poelmans <- unique(read.csv("data/poelmans_2023_SUI_GWAS_188genes.csv")$Gene)
tong10 <- unique(read.csv("data/tong2010_75genes.csv")$GeneSymbol)

r_chen06 <- redo_panel(chen06, "Chen 2006 (58 genes)")
r_chen03 <- redo_panel(chen03, "Chen 2003 (69 genes)")
r_poelmans <- redo_panel(poelmans, "Poelmans GWAS (183 genes)")
r_tong10 <- redo_panel(tong10, "Tong 2010 (66 genes)")

for (nm in c("r_chen06","r_chen03","r_poelmans","r_tong10")) {
  obj <- get(nm)
  tag <- sub("^r_", "", nm)
  if (!is.null(obj$go)) write.csv(obj$go, paste0("results/GO_BP_", tag, "_CORRECTED_universe.csv"), row.names = FALSE)
  if (!is.null(obj$kegg)) write.csv(obj$kegg, paste0("results/KEGG_", tag, "_CORRECTED_universe.csv"), row.names = FALSE)
}

## -----------------------------------------------------------------------
## Re-check the two "pathways in common with POP" comparisons that used
## these (now corrected) standalone panels. POP-alone GO/KEGG
## (results/GO_BP_GSE53868.csv, KEGG_GSE53868.csv) already used the real
## POP dataset gene list as universe - unaffected, reused as-is.
## -----------------------------------------------------------------------
compare_pathways <- function(tabA, tabB, id_col = "TERM_ID") {
  if (is.null(tabA) || is.null(tabB)) return(data.frame())
  sigA <- subset(tabA, p.adjust < 0.05); sigB <- subset(tabB, p.adjust < 0.05)
  common_ids <- intersect(sigA[[id_col]], sigB[[id_col]])
  if (length(common_ids) == 0) return(data.frame())
  merge(sigA[sigA[[id_col]] %in% common_ids, ], sigB[sigB[[id_col]] %in% common_ids, ],
        by = id_col, suffixes = c("_A", "_B"))
}

go_pop <- read.csv("results/GO_BP_GSE53868.csv")
kegg_pop <- read.csv("results/KEGG_GSE53868.csv")

cat("\n=== RE-CHECKED: Tong2010 (corrected) vs POP pathway overlap ===\n")
go_common_tong_fixed <- compare_pathways(r_tong10$go, go_pop)
kegg_common_tong_fixed <- compare_pathways(r_tong10$kegg, kegg_pop)
cat("GO BP common:", nrow(go_common_tong_fixed), "(previously reported: 19, with the inflated universe)\n")
cat("KEGG common:", nrow(kegg_common_tong_fixed), "(previously reported: 0)\n")
if (nrow(go_common_tong_fixed) > 0) write.csv(go_common_tong_fixed, "results/GO_BP_vias_comuns_Tong2010_GSE53868_CORRECTED.csv", row.names = FALSE)

cat("\n=== RE-CHECKED: Chen2006 (corrected) vs POP pathway overlap ===\n")
go_common_chen_fixed <- compare_pathways(r_chen06$go, go_pop)
kegg_common_chen_fixed <- compare_pathways(r_chen06$kegg, kegg_pop)
cat("GO BP common:", nrow(go_common_chen_fixed), "(previously reported: 15, with the inflated universe)\n")
cat("KEGG common:", nrow(kegg_common_chen_fixed), "(previously reported: 4, with the inflated universe)\n")
if (nrow(go_common_chen_fixed) > 0) write.csv(go_common_chen_fixed, "results/GO_BP_vias_comuns_Chen2006_GSE53868_CORRECTED.csv", row.names = FALSE)
if (nrow(kegg_common_chen_fixed) > 0) write.csv(kegg_common_chen_fixed, "results/KEGG_vias_comuns_Chen2006_GSE53868_CORRECTED.csv", row.names = FALSE)

cat("\n=== END ===\n")
