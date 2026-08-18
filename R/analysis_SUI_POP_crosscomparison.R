# =============================================================================
# Cruzamento SUI (GSE149072, rato) x POP (GSE208261, humano)
#
# ESTE SCRIPT REPRODUZ A ANALISE REAL rodada sobre:
#   - data/GSE149072_rawCounts.csv          (SUI, contagens brutas de rato)
#   - data/GSE208261_raw_counts_GRCh38.p13_NCBI.tsv  (POP, contagens brutas humanas)
# cujo resultado final esta em results/SUI_x_POP_36h_genes_comuns.xlsx
#
# Passos:
#   1) SUI: DEGs em uretra de rata lesionada SEM tratamento (Untreated) vs.
#      COM tratamento com hMSC (Treated), separadamente em 36h e 72h.
#   2) POP: DEGs em POP (n=12) vs. Controle (n=12), ajustando por grupo
#      etario (D=idosa, Y=jovem) como covariavel.
#   3) Converter os DEGs de rato (SUI) para ortologos humanos via HomoloGene.
#   4) Cruzar (interseccao) os ortologos humanos do SUI com os DEGs de POP,
#      em 36h e em 72h, e usar o tempo que efetivamente tiver correspondencias.
#
# Criterios: DESeq2, FDR (padj) < 0.05, |log2FoldChange (shrunk)| > 0.5
#
# IMPORTANTE: "gene em comum" significa que o MESMO gene aparece como
# significativo nas duas analises - isso NAO garante que ele muda na MESMA
# direcao nas duas doencas. O script classifica cada gene em comum como
# "Concordante" (mesma direcao de desregulacao) ou "Discordante" (direcoes
# opostas) - confira essa coluna antes de tirar qualquer conclusao biologica.
# =============================================================================

## -----------------------------------------------------------------------
## 1. Pacotes
## -----------------------------------------------------------------------

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("DESeq2", quietly = TRUE)) BiocManager::install("DESeq2", update = FALSE, ask = FALSE)
if (!requireNamespace("homologene", quietly = TRUE)) install.packages("homologene")
if (!requireNamespace("openxlsx", quietly = TRUE)) install.packages("openxlsx")

library(DESeq2)
library(homologene)
library(openxlsx)

dir.create("results", showWarnings = FALSE, recursive = TRUE)

FDR_CUT <- 0.05
LFC_CUT <- 0.5

## -----------------------------------------------------------------------
## 2. SUI (GSE149072) - DEGs em 36h e 72h
## -----------------------------------------------------------------------

sui_counts <- read.csv("data/GSE149072_rawCounts.csv", row.names = 1, check.names = FALSE)
sui_counts <- round(as.matrix(sui_counts)); mode(sui_counts) <- "integer"

run_sui_timepoint <- function(hour_label) {
  untreated <- grep(paste0("Rat_Urethra_Untreated_", hour_label), colnames(sui_counts), value = TRUE)
  treated   <- grep(paste0("Rat_Urethra_Treated_", hour_label), colnames(sui_counts), value = TRUE)
  sel <- c(untreated, treated)
  counts_sub <- sui_counts[, sel]

  group <- factor(c(rep("Untreated", length(untreated)), rep("Treated", length(treated))),
                   levels = c("Untreated", "Treated"))
  coldata <- data.frame(row.names = sel, group = group)

  keep <- rowSums(counts_sub >= 10) >= 3
  counts_sub <- counts_sub[keep, ]

  dds <- DESeqDataSetFromMatrix(countData = counts_sub, colData = coldata, design = ~ group)
  dds <- DESeq(dds)
  res <- results(dds, contrast = c("group", "Treated", "Untreated"), alpha = FDR_CUT)
  res_shrunk <- lfcShrink(dds, contrast = c("group", "Treated", "Untreated"), res = res, type = "normal")

  res_df <- as.data.frame(res_shrunk)
  res_df$Rat_Gene_Symbol <- rownames(res_df)
  res_df <- res_df[order(res_df$padj), ]
  rownames(res_df) <- NULL

  sig <- subset(res_df, !is.na(padj) & padj < FDR_CUT & abs(log2FoldChange) > LFC_CUT)
  cat("SUI", hour_label, "- genes testados:", nrow(res_df), "| DEGs:", nrow(sig), "\n")

  write.csv(res_df, paste0("results/SUI_", hour_label, "_DESeq2_completo.csv"), row.names = FALSE)
  sig
}

sui_sig_36h <- run_sui_timepoint("36hr")
sui_sig_72h <- run_sui_timepoint("72hr")

## -----------------------------------------------------------------------
## 3. Ortologos humanos para os DEGs de rato (HomoloGene)
## -----------------------------------------------------------------------

map_orthologs <- function(sig_df) {
  orth <- homologene(sig_df$Rat_Gene_Symbol, inTax = 10116, outTax = 9606)
  colnames(orth)[1:2] <- c("Rat_Gene_Symbol", "Human_Ortholog_Symbol")
  human_id_col <- grep("^9606_ID$|Gene.ID.*9606|human.*id", colnames(orth), ignore.case = TRUE, value = TRUE)[1]
  orth$Human_Entrez_ID <- orth[[human_id_col]]
  merged <- merge(sig_df, orth[, c("Rat_Gene_Symbol", "Human_Ortholog_Symbol", "Human_Entrez_ID")],
                   by = "Rat_Gene_Symbol", all.x = TRUE)
  merged[order(merged$padj), ]
}

sui_36h_orth <- map_orthologs(sui_sig_36h)
sui_72h_orth <- map_orthologs(sui_sig_72h)

write.csv(sui_36h_orth, "results/SUI_36hr_DEG_com_ortologos.csv", row.names = FALSE)
write.csv(sui_72h_orth, "results/SUI_72hr_DEG_com_ortologos.csv", row.names = FALSE)

cat("SUI 36h com ortologo humano:", sum(!is.na(sui_36h_orth$Human_Entrez_ID)), "\n")
cat("SUI 72h com ortologo humano:", sum(!is.na(sui_72h_orth$Human_Entrez_ID)), "\n")

## -----------------------------------------------------------------------
## 4. POP (GSE208261) - DEGs, ajustando por grupo etario
## -----------------------------------------------------------------------

pop_counts <- read.delim("data/GSE208261_raw_counts_GRCh38.p13_NCBI.tsv", row.names = 1, check.names = FALSE)
pop_counts <- round(as.matrix(pop_counts)); mode(pop_counts) <- "integer"

## Rotulos das amostras (GSM) conforme a pagina do GEO do GSE208261
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
coldata_pop <- data.frame(row.names = colnames(pop_counts), condition = condition, age_group = age_group)

keep_pop <- rowSums(pop_counts >= 10) >= 6
pop_counts_sub <- pop_counts[keep_pop, ]

dds_pop <- DESeqDataSetFromMatrix(countData = pop_counts_sub, colData = coldata_pop,
                                   design = ~ age_group + condition)
dds_pop <- DESeq(dds_pop)
res_pop <- results(dds_pop, contrast = c("condition", "POP", "Control"), alpha = FDR_CUT)
res_pop_shrunk <- lfcShrink(dds_pop, contrast = c("condition", "POP", "Control"), res = res_pop, type = "normal")

pop_df <- as.data.frame(res_pop_shrunk)
pop_df$Human_Entrez_ID <- rownames(pop_df)
pop_df <- pop_df[order(pop_df$padj), ]
rownames(pop_df) <- NULL

pop_sig <- subset(pop_df, !is.na(padj) & padj < FDR_CUT & abs(log2FoldChange) > LFC_CUT)
cat("POP - genes testados:", nrow(pop_df), "| DEGs:", nrow(pop_sig), "\n")

write.csv(pop_df, "results/POP_GSE208261_DESeq2_completo.csv", row.names = FALSE)
write.csv(pop_sig, "results/POP_GSE208261_DEG_sig.csv", row.names = FALSE)

## -----------------------------------------------------------------------
## 5. Cruzamento SUI x POP (36h e 72h) - usar o tempo com correspondencias
## -----------------------------------------------------------------------

pop_sig$Human_Entrez_ID <- as.character(pop_sig$Human_Entrez_ID)

cross_timepoint <- function(sui_orth_df, hour_label) {
  sui_valid <- subset(sui_orth_df, !is.na(Human_Entrez_ID))
  sui_valid$Human_Entrez_ID <- as.character(sui_valid$Human_Entrez_ID)
  common_ids <- intersect(sui_valid$Human_Entrez_ID, pop_sig$Human_Entrez_ID)
  cat(hour_label, "- genes em comum com POP:", length(common_ids), "\n")

  if (length(common_ids) == 0) return(NULL)

  sui_match <- sui_valid[sui_valid$Human_Entrez_ID %in% common_ids,
                          c("Rat_Gene_Symbol","Human_Ortholog_Symbol","Human_Entrez_ID","log2FoldChange","padj")]
  colnames(sui_match)[4:5] <- c(paste0("SUI_",hour_label,"_log2FC"), paste0("SUI_",hour_label,"_padj"))

  pop_match <- pop_sig[pop_sig$Human_Entrez_ID %in% common_ids, c("Human_Entrez_ID","log2FoldChange","padj")]
  colnames(pop_match)[2:3] <- c("POP_log2FC","POP_padj")

  merged <- merge(sui_match, pop_match, by = "Human_Entrez_ID")

  ## SUI e um contraste Treated vs Untreated: log2FC negativo = maior no
  ## grupo Untreated (lesao/doenca). POP e POP vs Controle: log2FC positivo
  ## = maior no POP (doenca). "Concordante" = mesma direcao de desregulacao
  ## nas duas doencas.
  sui_col <- paste0("SUI_", hour_label, "_log2FC")
  merged$Interpretacao_direcao <- ifelse(
    (merged[[sui_col]] < 0) == (merged$POP_log2FC > 0),
    "Concordante (mesma direcao nas 2 doencas)",
    "Discordante (direcoes opostas)"
  )
  merged <- merged[order(merged$POP_padj), ]
  write.csv(merged, paste0("results/SUI_", hour_label, "_x_POP_common_genes.csv"), row.names = FALSE)
  merged
}

cross_36h <- cross_timepoint(sui_36h_orth, "36hr")
cross_72h <- cross_timepoint(sui_72h_orth, "72hr")

cat("\n=== RESUMO FINAL ===\n")
cat("SUI 36h DEGs:", nrow(sui_sig_36h), "| SUI 72h DEGs:", nrow(sui_sig_72h), "\n")
cat("POP DEGs:", nrow(pop_sig), "\n")
cat("Genes em comum (36h):", if (is.null(cross_36h)) 0 else nrow(cross_36h), "\n")
cat("Genes em comum (72h):", if (is.null(cross_72h)) 0 else nrow(cross_72h), "\n")
cat("\nUse o resultado de 36h se 72h nao tiver correspondencias (era o caso nos dados originais).\n")

## -----------------------------------------------------------------------
## 6. Excel final
## -----------------------------------------------------------------------

if (!is.null(cross_36h)) {
  wb <- createWorkbook()
  addWorksheet(wb, "Genes_em_comum_36h")
  writeData(wb, "Genes_em_comum_36h", cross_36h)
  addWorksheet(wb, "SUI_36h_DEGs_completo")
  writeData(wb, "SUI_36h_DEGs_completo", sui_36h_orth)
  addWorksheet(wb, "SUI_72h_DEGs_completo")
  writeData(wb, "SUI_72h_DEGs_completo", sui_72h_orth)
  addWorksheet(wb, "POP_DEGs_completo")
  writeData(wb, "POP_DEGs_completo", pop_sig)
  saveWorkbook(wb, "results/SUI_x_POP_36h_genes_comuns.xlsx", overwrite = TRUE)
  cat("\nExcel final salvo em results/SUI_x_POP_36h_genes_comuns.xlsx\n")
}
