# =============================================================================
# Reproducao do "recipe" exato da usuaria (visto em ANALYSIS_SUI.docx) para
# explicar por que ela encontrou CALML5 (e outros) e a analise anterior
# (analysis_SUI_POP_nofilter.R) nao encontrou.
#
# Diferencas identificadas no documento dela, ALEM do filtro de baixa
# contagem (ja resolvido em analysis_SUI_POP_nofilter.R):
#   1. POP: design ~ condition (SEM covariavel de grupo etario D/Y)
#   2. POP e SUI: usa log2FoldChange BRUTO de results() - NUNCA chama
#      lfcShrink(). O filtro |log2FoldChange|>0.5 dela e sobre a estimativa
#      NAO encolhida (MLE), que para genes esparsos/baixa contagem pode ser
#      bem maior (ou, mais raramente, menor) que a versao shrunk que eu uso.
#   3. Ortologos via BioMart (nao reproduzido aqui - eu mantenho HomoloGene,
#      por instrucao explicita da usuaria de nao copiar os ortologos dela).
#
# Este script isola o efeito de (1) e (2), mantendo (3) = HomoloGene, para
# ver quanto disso explica a diferenca de genes em comum.
# =============================================================================

library(DESeq2)

FDR_CUT <- 0.05
LFC_CUT <- 0.5

load("data/homologeneData2.rda")
homologene_lookup <- function(rat_symbols) {
  rat_rows <- homologeneData2[homologeneData2$Taxonomy == 10116 &
                                 homologeneData2$Gene.Symbol %in% rat_symbols, c("HID", "Gene.Symbol")]
  colnames(rat_rows) <- c("HID", "Rat_Gene_Symbol")
  human_rows <- homologeneData2[homologeneData2$Taxonomy == 9606, c("HID", "Gene.Symbol", "Gene.ID")]
  colnames(human_rows) <- c("HID", "Human_Ortholog_Symbol", "Human_Entrez_ID")
  merge(rat_rows, human_rows, by = "HID")[, c("Rat_Gene_Symbol", "Human_Ortholog_Symbol", "Human_Entrez_ID")]
}

## -----------------------------------------------------------------------
## SUI 36h - SEM filtro, SEM lfcShrink (log2FC bruto de results())
## -----------------------------------------------------------------------
sui_counts <- read.csv("data/GSE149072_rawCounts.csv", row.names = 1, check.names = FALSE)
sui_counts <- round(as.matrix(sui_counts)); mode(sui_counts) <- "integer"

untreated <- grep("Rat_Urethra_Untreated_36hr", colnames(sui_counts), value = TRUE)
treated   <- grep("Rat_Urethra_Treated_36hr", colnames(sui_counts), value = TRUE)
sel <- c(untreated, treated)
sui_counts_sub <- sui_counts[, sel]

group <- factor(c(rep("Untreated", length(untreated)), rep("Treated", length(treated))),
                 levels = c("Untreated", "Treated"))
coldata_sui <- data.frame(row.names = sel, group = group)

dds_sui <- DESeqDataSetFromMatrix(countData = sui_counts_sub, colData = coldata_sui, design = ~ group)
dds_sui <- DESeq(dds_sui)
res_sui <- results(dds_sui, contrast = c("group", "Treated", "Untreated"), alpha = FDR_CUT)  # SEM shrink

sui_df <- as.data.frame(res_sui)
sui_df$Rat_Gene_Symbol <- rownames(sui_df)
sui_df <- sui_df[order(sui_df$padj), ]
rownames(sui_df) <- NULL

sui_sig <- subset(sui_df, !is.na(padj) & padj < FDR_CUT & abs(log2FoldChange) > LFC_CUT)
cat("SUI 36h (sem filtro, SEM shrink) - DEGs:", nrow(sui_sig), "\n")

orth <- homologene_lookup(sui_sig$Rat_Gene_Symbol)
sui_sig_orth <- merge(sui_sig, orth, by = "Rat_Gene_Symbol", all.x = TRUE)
sui_sig_orth <- sui_sig_orth[order(sui_sig_orth$padj), ]

cat("Genes de interesse da usuaria estao no DEG list do SUI (sem shrink)?\n")
for (g in c("Cwh43","Inpp4b","Calml5","Krt10","Serpinb2","Dmkn")) {
  in_list <- g %in% sui_sig$Rat_Gene_Symbol
  cat(" ", g, ":", in_list, "\n")
}

## -----------------------------------------------------------------------
## POP - SEM filtro, SEM age covariate, SEM lfcShrink
## -----------------------------------------------------------------------
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
coldata_pop <- data.frame(row.names = colnames(pop_counts), condition = condition)  # SEM age_group

dds_pop <- DESeqDataSetFromMatrix(countData = pop_counts, colData = coldata_pop, design = ~ condition)
dds_pop <- DESeq(dds_pop)
res_pop <- results(dds_pop, contrast = c("condition", "POP", "Control"), alpha = FDR_CUT)  # SEM shrink

pop_df <- as.data.frame(res_pop)
pop_df$Human_Entrez_ID <- rownames(pop_df)
pop_df <- pop_df[order(pop_df$padj), ]
rownames(pop_df) <- NULL

pop_sig <- subset(pop_df, !is.na(padj) & padj < FDR_CUT & abs(log2FoldChange) > LFC_CUT)
cat("\nPOP (sem filtro, sem age covariate, SEM shrink) - DEGs:", nrow(pop_sig), "\n")

## -----------------------------------------------------------------------
## Cruzamento
## -----------------------------------------------------------------------
pop_sig$Human_Entrez_ID <- as.character(pop_sig$Human_Entrez_ID)
sui_valid <- subset(sui_sig_orth, !is.na(Human_Entrez_ID))
sui_valid$Human_Entrez_ID <- as.character(sui_valid$Human_Entrez_ID)

common_ids <- intersect(sui_valid$Human_Entrez_ID, pop_sig$Human_Entrez_ID)
common_genes <- sui_valid[sui_valid$Human_Entrez_ID %in% common_ids, c("Rat_Gene_Symbol","Human_Ortholog_Symbol","Human_Entrez_ID")]
cat("\nGenes em comum (recipe da usuaria, ortologos HomoloGene):", length(common_ids), "\n")
print(common_genes)

cat("\n=== Comparacao com a lista dela (CWH43, INPP4B, CALML5, KRT10, SERPINB2, DMKN) ===\n")
dela <- c("CWH43","INPP4B","CALML5","KRT10","SERPINB2","DMKN")
minha <- toupper(common_genes$Human_Ortholog_Symbol)
cat("Bateram:", paste(intersect(dela, minha), collapse=", "), "\n")
cat("So na dela (nao achei):", paste(setdiff(dela, minha), collapse=", "), "\n")
cat("So na minha (ela nao tem):", paste(setdiff(minha, dela), collapse=", "), "\n")

## Verificar se os "so na dela" tem ortologo no HomoloGene (rato->humano)
cat("\nEsses genes tem ortologo de rato no HomoloGene?\n")
for (hg in setdiff(dela, minha)) {
  rat_candidates <- homologeneData2[homologeneData2$Taxonomy == 9606 &
                                       toupper(homologeneData2$Gene.Symbol) == hg, "HID"]
  if (length(rat_candidates) == 0) { cat(" ", hg, ": nao encontrado no HomoloGene em humano\n"); next }
  hid <- rat_candidates[1]
  taxa_no_hid <- homologeneData2$Taxonomy[homologeneData2$HID == hid]
  cat(" ", hg, "- HID", hid, "- taxons no grupo:", paste(sort(unique(taxa_no_hid)), collapse=","),
      "- tem rato (10116)?", 10116 %in% taxa_no_hid, "\n")
}
