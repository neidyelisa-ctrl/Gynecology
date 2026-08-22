# =============================================================================
# SUI (36h) x POP - SEM filtro de baixa contagem, a pedido da usuaria
#
# Motivacao: a usuaria identificou que Cwh43 e Calml5 (genes com contagem
# quase zero, so 1-2 amostras com sinal) foram removidos pelo filtro de
# baixa expressao (rowSums(contagem>=10) >= N amostras) usado no pipeline
# principal (analysis_SUI_POP_crosscomparison.R) antes mesmo de entrarem no
# DESeq2. Ela quer reproduzir sem esse filtro para ver o efeito.
#
# ATENCAO METODOLOGICA: sem o filtro de baixa contagem, genes esparsos
# (quase todos zeros, com 1-2 outliers altos) tendem a gerar fold-change
# artificialmente alto e p-valor artificialmente baixo - um problema
# classico de RNA-seq (ver Love, Huber & Anders 2014, a documentacao do
# DESeq2 recomenda EXPLICITAMENTE pre-filtrar genes de baixa contagem por
# esse motivo). Os resultados deste script devem ser lidos com essa ressalva:
# qualquer DEG que so aparece aqui (e nao no pipeline com filtro) e candidato
# a falso positivo tecnico, nao biologico - merece confirmacao por qPCR
# antes de qualquer conclusao.
#
# Mantido igual ao pipeline principal: FDR (padj) < 0.05, |log2FC| > 0.5,
# lfcShrink tipo "normal", POP com covariavel de grupo etario (D/Y), SUI 36h
# = Treated vs Untreated, ortologos via NCBI HomoloGene (metodologia propria,
# nao a do BioMart usada pela usuaria).
# =============================================================================

library(DESeq2)
library(openxlsx)

dir.create("results", showWarnings = FALSE, recursive = TRUE)

FDR_CUT <- 0.05
LFC_CUT <- 0.5

## Mapeamento de ortologos via NCBI HomoloGene (dados brutos, mesma fonte
## usada no pacote 'homologene' - agrupamento por HID entre Taxonomy IDs
## 10116=rato e 9606=humano). Reimplementado aqui em vez de depender do
## pacote 'homologene' (indisponivel neste sandbox sem acesso a CRAN).
load("data/homologeneData2.rda")  # -> homologeneData2 (HID, Gene.Symbol, Taxonomy, Gene.ID)

homologene_lookup <- function(rat_symbols) {
  rat_rows <- homologeneData2[homologeneData2$Taxonomy == 10116 &
                                 homologeneData2$Gene.Symbol %in% rat_symbols, c("HID", "Gene.Symbol")]
  colnames(rat_rows) <- c("HID", "Rat_Gene_Symbol")
  human_rows <- homologeneData2[homologeneData2$Taxonomy == 9606, c("HID", "Gene.Symbol", "Gene.ID")]
  colnames(human_rows) <- c("HID", "Human_Ortholog_Symbol", "Human_Entrez_ID")
  merge(rat_rows, human_rows, by = "HID")[, c("Rat_Gene_Symbol", "Human_Ortholog_Symbol", "Human_Entrez_ID")]
}

## -----------------------------------------------------------------------
## 1. SUI (GSE149072) - 36h, SEM filtro de baixa contagem
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

## SEM filtro de baixa contagem (a pedido da usuaria) - todos os genes do
## arquivo de contagem entram direto no DESeq2.
dds_sui <- DESeqDataSetFromMatrix(countData = sui_counts_sub, colData = coldata_sui, design = ~ group)
dds_sui <- DESeq(dds_sui)
res_sui <- results(dds_sui, contrast = c("group", "Treated", "Untreated"), alpha = FDR_CUT)
res_sui_shrunk <- lfcShrink(dds_sui, contrast = c("group", "Treated", "Untreated"), res = res_sui, type = "normal")

sui_df <- as.data.frame(res_sui_shrunk)
sui_df$Rat_Gene_Symbol <- rownames(sui_df)
sui_df <- sui_df[order(sui_df$padj), ]
rownames(sui_df) <- NULL

sui_sig <- subset(sui_df, !is.na(padj) & padj < FDR_CUT & abs(log2FoldChange) > LFC_CUT)
cat("SUI 36h (SEM filtro) - genes testados:", nrow(sui_df), "| DEGs:", nrow(sui_sig), "\n")

write.csv(sui_df, "results/SUI_36hr_nofilter_DESeq2_completo.csv", row.names = FALSE)

## -----------------------------------------------------------------------
## 2. Ortologos humanos para os DEGs de rato (HomoloGene) - ANTES de cruzar
## -----------------------------------------------------------------------

orth <- homologene_lookup(sui_sig$Rat_Gene_Symbol)

sui_sig_orth <- merge(sui_sig, orth, by = "Rat_Gene_Symbol", all.x = TRUE)
sui_sig_orth <- sui_sig_orth[order(sui_sig_orth$padj), ]
write.csv(sui_sig_orth, "results/SUI_36hr_nofilter_DEG_com_ortologos.csv", row.names = FALSE)

cat("SUI 36h (sem filtro) com ortologo humano encontrado:",
    sum(!is.na(sui_sig_orth$Human_Entrez_ID)), "de", nrow(sui_sig_orth), "DEGs\n")

## -----------------------------------------------------------------------
## 3. POP (GSE208261) - SEM filtro de baixa contagem
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
age_group <- factor(ifelse(grepl("_D", label), "D", "Y"))
coldata_pop <- data.frame(row.names = colnames(pop_counts), condition = condition, age_group = age_group)

## SEM filtro de baixa contagem (a pedido da usuaria).
dds_pop <- DESeqDataSetFromMatrix(countData = pop_counts, colData = coldata_pop,
                                   design = ~ age_group + condition)
dds_pop <- DESeq(dds_pop)
res_pop <- results(dds_pop, contrast = c("condition", "POP", "Control"), alpha = FDR_CUT)
res_pop_shrunk <- lfcShrink(dds_pop, contrast = c("condition", "POP", "Control"), res = res_pop, type = "normal")

pop_df <- as.data.frame(res_pop_shrunk)
pop_df$Human_Entrez_ID <- rownames(pop_df)
pop_df <- pop_df[order(pop_df$padj), ]
rownames(pop_df) <- NULL

pop_sig <- subset(pop_df, !is.na(padj) & padj < FDR_CUT & abs(log2FoldChange) > LFC_CUT)
cat("POP (SEM filtro) - genes testados:", nrow(pop_df), "| DEGs:", nrow(pop_sig), "\n")

write.csv(pop_df, "results/POP_GSE208261_nofilter_DESeq2_completo.csv", row.names = FALSE)
write.csv(pop_sig, "results/POP_GSE208261_nofilter_DEG_sig.csv", row.names = FALSE)

## -----------------------------------------------------------------------
## 4. Cruzamento SUI(36h, ortologo humano) x POP
## -----------------------------------------------------------------------

pop_sig$Human_Entrez_ID <- as.character(pop_sig$Human_Entrez_ID)
sui_valid <- subset(sui_sig_orth, !is.na(Human_Entrez_ID))
sui_valid$Human_Entrez_ID <- as.character(sui_valid$Human_Entrez_ID)

common_ids <- intersect(sui_valid$Human_Entrez_ID, pop_sig$Human_Entrez_ID)
cat("\nGenes em comum SUI(36h, sem filtro) x POP(sem filtro):", length(common_ids), "\n")

common_genes <- NULL
if (length(common_ids) > 0) {
  sui_match <- sui_valid[sui_valid$Human_Entrez_ID %in% common_ids,
                          c("Rat_Gene_Symbol","Human_Ortholog_Symbol","Human_Entrez_ID","log2FoldChange","padj")]
  colnames(sui_match)[4:5] <- c("SUI_36h_log2FC","SUI_36h_padj")

  pop_match <- pop_sig[pop_sig$Human_Entrez_ID %in% common_ids, c("Human_Entrez_ID","log2FoldChange","padj")]
  colnames(pop_match)[2:3] <- c("POP_log2FC","POP_padj")

  common_genes <- merge(sui_match, pop_match, by = "Human_Entrez_ID")
  common_genes$Interpretacao_direcao <- ifelse(
    (common_genes$SUI_36h_log2FC < 0) == (common_genes$POP_log2FC > 0),
    "Concordante (mesma direcao nas 2 doencas)",
    "Discordante (direcoes opostas)"
  )
  common_genes <- common_genes[order(common_genes$POP_padj), ]
  print(common_genes[, c("Rat_Gene_Symbol","Human_Ortholog_Symbol","SUI_36h_log2FC","POP_log2FC","Interpretacao_direcao")])
  write.csv(common_genes, "results/SUI_36hr_nofilter_x_POP_common_genes.csv", row.names = FALSE)
}

cat("\n=== RESUMO ===\n")
cat("SUI 36h DEGs (sem filtro):", nrow(sui_sig), "| POP DEGs (sem filtro):", nrow(pop_sig), "\n")
cat("Genes em comum:", length(common_ids), "\n")
