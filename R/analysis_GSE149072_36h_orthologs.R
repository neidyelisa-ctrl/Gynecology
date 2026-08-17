# =============================================================================
# GSE149072 - Gene expression profiling of tissue and hMSC xenografts in a
# rat postpartum urinary injury model (Sadeghi et al. 2020, Tissue Eng Part A,
# DOI 10.1089/ten.tea.2020.0033)
#
# ESTE SCRIPT REPRODUZ A ANALISE REAL RODADA SOBRE data/GSE149072_rawCounts.csv
# (contagens brutas de RNA-seq baixadas pelo usuario da pagina do GSE149072
# no GEO), cujo resultado esta em results/GSE149072_36h_DEGs_ortologos_humanos.xlsx.
#
# Comparacao: amostras de uretra de rata em 36h pos-lesao por distensao
# vaginal, SEM tratamento com hMSC (Rat_Urethra_Untreated_36hr, n=3) vs. COM
# tratamento com hMSC (Rat_Urethra_Treated_36hr, n=3). Esta e a unica
# comparacao possivel no arquivo de contagens brutas para o tempo de 36h -
# nao ha amostras de uretra normal/nao-lesionada (controle saudavel) neste
# arquivo, apenas lesao-sem-tratamento, lesao-com-tratamento e as celulas
# hMSC humanas isoladas (colunas "MSC_*").
#
# log2FoldChange positivo = maior expressao no grupo Treated (efeito do hMSC)
# log2FoldChange negativo = maior expressao no grupo Untreated (doente)
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

## -----------------------------------------------------------------------
## 2. Ler as contagens brutas
## -----------------------------------------------------------------------

counts_file <- "data/GSE149072_rawCounts.csv"
df <- read.csv(counts_file, row.names = 1, check.names = FALSE)
df <- round(as.matrix(df)); mode(df) <- "integer"

cat("Total de amostras no arquivo:", ncol(df), "| genes:", nrow(df), "\n")

## -----------------------------------------------------------------------
## 3. Selecionar as amostras de uretra de rato em 36h
## -----------------------------------------------------------------------

untreated_36h <- grep("Rat_Urethra_Untreated_36hr", colnames(df), value = TRUE)
treated_36h   <- grep("Rat_Urethra_Treated_36hr", colnames(df), value = TRUE)

cat("\nAmostras Untreated_36h (lesao, sem tratamento):\n"); print(untreated_36h)
cat("\nAmostras Treated_36h (lesao + hMSC):\n"); print(treated_36h)

sel <- c(untreated_36h, treated_36h)
counts_sub <- df[, sel]

group <- factor(c(rep("Untreated", length(untreated_36h)),
                   rep("Treated", length(treated_36h))),
                levels = c("Untreated", "Treated"))
coldata <- data.frame(row.names = sel, group = group)

## Aviso de QC: confira a profundidade de sequenciamento das amostras usadas
cat("\nProfundidade de sequenciamento (soma de contagens) por amostra:\n")
print(colSums(counts_sub))

## -----------------------------------------------------------------------
## 4. Filtro de baixa expressao + DESeq2
## -----------------------------------------------------------------------

keep <- rowSums(counts_sub >= 10) >= 3
cat("\nGenes antes do filtro:", nrow(counts_sub), "| depois do filtro:", sum(keep), "\n")
counts_sub <- counts_sub[keep, ]

dds <- DESeqDataSetFromMatrix(countData = counts_sub, colData = coldata, design = ~ group)
dds <- DESeq(dds)

res <- results(dds, contrast = c("group", "Treated", "Untreated"), alpha = 0.05)
## Shrinkage do log2FoldChange (reduz inflacao em genes de baixa contagem)
res_shrunk <- lfcShrink(dds, contrast = c("group", "Treated", "Untreated"), res = res, type = "normal")

res_df <- as.data.frame(res_shrunk)
res_df$Rat_Gene_Symbol <- rownames(res_df)
res_df <- res_df[order(res_df$padj), ]
rownames(res_df) <- NULL

sig <- subset(res_df, !is.na(padj) & padj < 0.05 & abs(log2FoldChange) > 1)
cat("\nGenes testados:", nrow(res_df), "| DEGs significativos (padj<0.05, |log2FC|>1):", nrow(sig), "\n")

## -----------------------------------------------------------------------
## 5. Ortologos humanos (pacote CRAN 'homologene', base NCBI HomoloGene)
## -----------------------------------------------------------------------

orth <- homologene(sig$Rat_Gene_Symbol, inTax = 10116, outTax = 9606)
## colunas tipicas: <symbol_rato>, <symbol_humano>, <geneID_rato>, <geneID_humano>
colnames(orth)[1:2] <- c("Rat_Gene_Symbol", "Human_Ortholog_Symbol")

deg_com_ortologos <- merge(sig, orth[, c("Rat_Gene_Symbol", "Human_Ortholog_Symbol")],
                            by = "Rat_Gene_Symbol", all.x = TRUE)
deg_com_ortologos <- deg_com_ortologos[order(deg_com_ortologos$padj), ]

n_com <- sum(!is.na(deg_com_ortologos$Human_Ortholog_Symbol))
cat("Genes com ortologo humano identificado:", n_com, "de", nrow(sig), "\n")

## -----------------------------------------------------------------------
## 6. Exportar
## -----------------------------------------------------------------------

write.csv(res_df, "results/GSE149072_36h_DESeq2_completo.csv", row.names = FALSE)
write.csv(deg_com_ortologos, "results/GSE149072_36h_DEGs_com_ortologos.csv", row.names = FALSE)

wb <- createWorkbook()
addWorksheet(wb, "DEGs_significativos_ortologos")
writeData(wb, "DEGs_significativos_ortologos", deg_com_ortologos)
addWorksheet(wb, "DEGs_completo_36h")
writeData(wb, "DEGs_completo_36h", res_df)
saveWorkbook(wb, "results/GSE149072_36h_DEGs_ortologos_humanos.xlsx", overwrite = TRUE)

cat("\nAnalise concluida. Resultado final em results/GSE149072_36h_DEGs_ortologos_humanos.xlsx\n")
