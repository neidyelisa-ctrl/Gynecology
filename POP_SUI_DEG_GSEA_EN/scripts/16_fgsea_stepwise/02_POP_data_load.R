# =============================================================================
# STEP 2 of 13 - load POP data (GSE208261) and build the 12x12 design
# =============================================================================
# RUN ORDER: needs 01_setup_packages_and_data_check.R already sourced in this
# same R session (uses: select, dir "results"/"figures" already created).
#
# WHAT THIS STEP DOES: loads the sample metadata and the raw RNA-seq count
# matrix for GSE208261, converts gene IDs from Entrez to symbol, and builds
# the 12 Control vs 12 POP group design used by every POP step after this one.
#
# OUTPUT: no files written yet (that starts at step 3) - this step only
# builds R objects in memory: meta, counts_all, group_full, tissue_full.
# =============================================================================

setwd("E:/POP+SUI ONE")  # safe to repeat - keeps this correct even in a fresh R session

cat("\n================ STEP 2: POP data (GSE208261) ================\n\n")

meta <- read.csv("GSE208261_sample_metadata.csv")
cat("Sample table (check this against the GEO metadata before trusting it blindly):\n")
print(meta)
cat("\nTissue x Group cross-tab:\n")
print(table(meta$Tissue, meta$Treatment))
cat("\n")

counts_raw <- read.delim("GSE208261_raw_counts_GRCh38.p13_NCBI.tsv", row.names = 1, check.names = FALSE)
counts_all <- as.matrix(counts_raw[, meta$GSM])
stopifnot(identical(colnames(counts_all), meta$GSM))
cat("Count matrix: 24 samples,", nrow(counts_all), "genes (by Entrez ID)\n")

# Genes come identified by number (Entrez ID) - convert to symbol for
# readability and comparability with the SUI data (which already uses names).
ann_id <- suppressWarnings(select(org.Hs.eg.db, keys = rownames(counts_all), keytype = "ENTREZID", columns = "SYMBOL"))
ann_id <- ann_id[!is.na(ann_id$SYMBOL) & !duplicated(ann_id$ENTREZID), ]
counts_all <- counts_all[rownames(counts_all) %in% ann_id$ENTREZID, ]
rownames(counts_all) <- ann_id$SYMBOL[match(rownames(counts_all), ann_id$ENTREZID)]
counts_all <- rowsum(counts_all, group = rownames(counts_all))  # sums duplicate genes
cat("After ID->symbol conversion and summing duplicates:", nrow(counts_all), "genes\n\n")

# The 12x12 design: 12 Control (6 ligament + 6 vaginal wall, COMBINED) vs
# 12 POP (vaginal wall).
group_full <- factor(meta$Treatment, levels = c("Control", "POP"))
tissue_full <- factor(meta$Tissue)
cat("Design: 12 Control (6 ligament + 6 vaginal wall) vs 12 POP (vaginal wall). Control:",
    sum(group_full == "Control"), "| POP:", sum(group_full == "POP"), "\n\n")

cat("=== STEP 2 DONE === Now Source 03_POP_DEG_DESeq2.R\n")
