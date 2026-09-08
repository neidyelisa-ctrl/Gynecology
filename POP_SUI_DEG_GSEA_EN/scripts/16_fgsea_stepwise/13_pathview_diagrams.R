# =============================================================================
# STEP 13 of 13 - Pathview diagrams: official KEGG pathway maps, colored by
# YOUR own POP and SUI fold-change data - one diagram per pathway, POP and
# SUI shown side by side inside each gene box.
# =============================================================================
# RUN ORDER: needs step 6 (results/KEGG_significant_in_BOTH_diseases_
# fgsea.csv) and step 3 (results/POP_DESeq2_full_table.csv) and step 5
# (results/SUI_Wei2020_full_table.csv) already run at least once - but
# unlike every other step in this pipeline, this one reads everything back
# from the CSV files on disk rather than from R objects still in memory, so
# it can also be run on its own, in a brand-new R session, any time after
# steps 1-6 have been run once.
#
# WHICH PATHWAYS GET A DIAGRAM: decided automatically from the fgsea results
# step 6 already produced - results/KEGG_significant_in_BOTH_diseases_
# fgsea.csv, i.e. every KEGG pathway with padj<0.25 in BOTH POP and SUI.
# This is deliberately NOT a hardcoded pathway list: whatever comes out of
# your fgsea run this time (5 pathways, 8, whatever it is) is what gets
# rendered - if you re-run steps 1-6 later and the significant set changes
# slightly, re-running this step picks that up automatically, no manual
# editing.
#
# REQUIRES: steps 1-6 already run at least once in this same folder (needs
# their output files in results/). Also requires internet the first time
# each pathway ID is rendered - pathview downloads the official KEGG diagram
# (XML+PNG) for that pathway and caches it locally; re-running later for the
# same pathway ID reuses the cached file, no internet needed then.
#
# OUTPUT: figures/pathview/hsa<PATHID>.POP_vs_SUI.png - one per significant
# pathway, each gene box split into two colors (left = POP, right = SUI;
# red = upregulated, blue = downregulated, matching the color scheme used
# throughout this project's other figures).
# =============================================================================

setwd("E:/POP+SUI ONE")  # SAME folder as the rest of the pipeline - change this one line only
cat("Working directory set to:", getwd(), "\n\n")

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
for (pkg in c("pathview", "org.Hs.eg.db", "AnnotationDbi")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("Installing (Bioconductor):", pkg, "...\n")
    BiocManager::install(pkg, update = FALSE, ask = FALSE)
  }
}
suppressMessages({
  library(pathview)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
})

required_files <- c(
  "results/KEGG_significant_in_BOTH_diseases_fgsea.csv",
  "results/POP_DESeq2_full_table.csv",
  "results/SUI_Wei2020_full_table.csv"
)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop(
    "\n\nMissing file(s) from earlier steps' output:\n  - ", paste(missing_files, collapse = "\n  - "),
    "\n\nRun steps 01-06 in this same folder first (01_setup_packages_and_data_check.R\n",
    "through 06_shared_pathways_POP_SUI.R), then run this step again.\n"
  )
}

## --- which pathways get a diagram - decided from THIS run's fgsea results,
## not a fixed list -----------------------------------------------------------
kegg_both <- read.csv("results/KEGG_significant_in_BOTH_diseases_fgsea.csv", colClasses = c(PATH = "character"))
cat("Pathways significant in BOTH POP and SUI (FDR<0.25) this run:", nrow(kegg_both), "\n")
if (nrow(kegg_both) > 0) {
  print(kegg_both[, c("PATH", "PathwayName", "NES_POP", "padj_POP", "NES_SUI", "padj_SUI")])
}
cat("\n")

if (nrow(kegg_both) == 0) {
  stop(
    "\n\nNo pathway is significant in both POP and SUI in this run's fgsea\n",
    "results - nothing to draw. This is a real result, not a script error:\n",
    "check results/KEGG_all_common_pathways_fgsea.csv for near-significant\n",
    "pathways if you want to pick a few by hand instead.\n"
  )
}
pathway_ids <- kegg_both$PATH

## --- per-gene fold-change data for POP and SUI, keyed by Entrez ID (KEGG
## diagrams are keyed by Entrez ID, not gene symbol) -------------------------
pop_full <- read.csv("results/POP_DESeq2_full_table.csv")
sui_full <- read.csv("results/SUI_Wei2020_full_table.csv")

sym_to_entrez <- function(symbols) {
  ann <- suppressWarnings(select(org.Hs.eg.db, keys = unique(symbols), keytype = "SYMBOL", columns = "ENTREZID"))
  ann[!is.na(ann$ENTREZID) & !duplicated(ann$SYMBOL), ]
}

pop_ann <- sym_to_entrez(pop_full$Gene)
pop_lfc <- setNames(pop_full$logFC[match(pop_ann$SYMBOL, pop_full$Gene)], pop_ann$ENTREZID)

sui_ann <- sym_to_entrez(sui_full$GeneSymbol)
sui_lfc <- setNames(sui_full$logFC[match(sui_ann$SYMBOL, sui_full$GeneSymbol)], sui_ann$ENTREZID)

all_entrez <- union(names(pop_lfc), names(sui_lfc))
gene_mat <- cbind(POP = pop_lfc[all_entrez], SUI = sui_lfc[all_entrez])
rownames(gene_mat) <- all_entrez
cat("Combined POP+SUI logFC matrix ready:", nrow(gene_mat), "genes (Entrez ID) x 2 diseases\n")
cat("(NA for a gene in one column just means that gene wasn't in that\n")
cat("disease's own tested/DEG list - pathview draws it grey, not an error)\n\n")

## --- render one diagram per significant pathway ----------------------------
dir.create("figures/pathview", showWarnings = FALSE, recursive = TRUE)
orig_wd <- getwd()
setwd("figures/pathview")  # pathview writes its output PNGs to the working directory

cat("Rendering", length(pathway_ids), "pathway diagram(s) - needs internet the first\n")
cat("time each pathway ID is drawn (downloads the official KEGG map)...\n\n")

## pathview's own KEGG download step can fail with just a WARNING (not an R
## error) - e.g. a transient server hiccup - in which case it silently skips
## that pathway and returns without producing an output file, but WITHOUT
## raising anything a plain tryCatch(error=...) would catch. So success is
## verified here by checking the actual output PNG exists, not by the
## absence of an R error, and a failed download gets a couple of automatic
## retries (transient KEGG server issues are common and usually clear up in
## a few seconds).
n_retries <- 2
rendered <- character(0)
failed <- character(0)
for (pid in pathway_ids) {
  cat("Pathway", pid, "(", kegg_both$PathwayName[kegg_both$PATH == pid], ") ... ")
  expected_file <- paste0("hsa", pid, ".POP_vs_SUI.png")
  ok <- FALSE
  for (attempt in seq_len(n_retries + 1)) {
    tryCatch({
      suppressWarnings(pathview(gene.data = gene_mat, pathway.id = pid, species = "hsa",
               gene.idtype = "ENTREZID", out.suffix = "POP_vs_SUI",
               kegg.native = TRUE, same.layer = TRUE,
               low = list(gene = "#2166AC"), mid = list(gene = "#F7F7F7"), high = list(gene = "#B2182B"),
               na.col = "grey85"))
    }, error = function(e) NULL)
    if (file.exists(expected_file)) { ok <- TRUE; break }
    if (attempt <= n_retries) { cat("(retry", attempt, ") "); Sys.sleep(2) }
  }
  if (ok) {
    cat("done\n")
    rendered <- c(rendered, pid)
  } else {
    cat("FAILED after", n_retries + 1, "attempt(s) - likely a transient KEGG server issue\n")
    failed <- c(failed, pid)
  }
}
setwd(orig_wd)

cat("\n=== Summary ===\n")
cat("Rendered:", length(rendered), "of", length(pathway_ids), "pathways -> figures/pathview/\n")
if (length(failed) > 0) {
  cat("Could not render:", paste(failed, collapse = ", "), "even after retrying -\n")
  cat("usually a transient KEGG server issue, occasionally means no internet\n")
  cat("reachable or that ID has no KEGG map (check the pathway name against\n")
  cat("kegg.jp manually if that's the case). Re-running this script later will\n")
  cat("retry only these - the rest are already cached on disk and pathview will\n")
  cat("reuse them, not re-download.\n")
}
cat("\nEach output file has POP's fold-change on the left half of every gene\n")
cat("box and SUI's on the right half - red = up, blue = down, grey = gene\n")
cat("not in that disease's tested list. Compare directly against the\n")
cat("Focal adhesion / Ribosome / etc. discussion already written up for the\n")
cat("thesis - this is the same numbers, drawn onto the real KEGG diagram.\n")

cat("\n=== STEP 13 DONE - all 13 steps of the pipeline complete ===\n")
