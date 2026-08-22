library(DESeq2)
pop_counts <- read.delim("data/GSE208261_raw_counts_GRCh38.p13_NCBI.tsv", row.names = 1, check.names = FALSE)
pop_counts <- round(as.matrix(pop_counts)); mode(pop_counts) <- "integer"

sample_labels <- c(
  GSM6339911="Control_D1", GSM6339912="Control_D2", GSM6339913="Control_D3",
  GSM6339914="Control_D4", GSM6339915="Control_D5", GSM6339916="Control_D6",
  GSM6339917="Control_Y1", GSM6339918="Control_Y2", GSM6339919="Control_Y3",
  GSM6339920="Control_Y4", GSM6339921="Control_Y5", GSM6339922="Control_Y6",
  GSM6339923="POP_D1", GSM6339924="POP_D2", GSM6339925="POP_D3",
  GSM6339926="POP_D4", GSM6339927="POP_D5", GSM6339928="POP_D6",
  GSM6339929="POP_Y1", GSM6339930="POP_Y2", GSM6339931="POP_Y3",
  GSM6339932="POP_Y4", GSM6339933="POP_Y5", GSM6339934="POP_Y6"
)
label <- sample_labels[colnames(pop_counts)]
condition <- factor(ifelse(grepl("^POP", label), "POP", "Control"), levels = c("Control", "POP"))
age_group <- factor(ifelse(grepl("_D", label), "D", "Y"))

genes_of_interest <- c(CWH43="80157", INPP4B="8821", CALML5="51806", KRT10="3858", SERPINB2="5055", DMKN="93099", ZCCHC12="170261")

run_combo <- function(use_age, use_shrink) {
  coldata <- if (use_age) data.frame(row.names=colnames(pop_counts), condition=condition, age_group=age_group) else data.frame(row.names=colnames(pop_counts), condition=condition)
  design <- if (use_age) ~ age_group + condition else ~ condition
  dds <- DESeqDataSetFromMatrix(countData=pop_counts, colData=coldata, design=design)
  dds <- DESeq(dds)
  res <- results(dds, contrast=c("condition","POP","Control"), alpha=0.05)
  if (use_shrink) res <- lfcShrink(dds, contrast=c("condition","POP","Control"), res=res, type="normal")
  df <- as.data.frame(res)
  df$id <- rownames(df)
  df[df$id %in% genes_of_interest, c("id","log2FoldChange","padj")]
}

combos <- list(
  "age+shrink"    = run_combo(TRUE, TRUE),
  "age+noshrink"  = run_combo(TRUE, FALSE),
  "noage+shrink"  = run_combo(FALSE, TRUE),
  "noage+noshrink"= run_combo(FALSE, FALSE)
)

for (nm in names(combos)) {
  cat("\n===", nm, "===\n")
  d <- combos[[nm]]
  d$gene <- names(genes_of_interest)[match(d$id, genes_of_interest)]
  d$sig <- !is.na(d$padj) & d$padj < 0.05 & abs(d$log2FoldChange) > 0.5
  print(d[, c("gene","log2FoldChange","padj","sig")])
}
