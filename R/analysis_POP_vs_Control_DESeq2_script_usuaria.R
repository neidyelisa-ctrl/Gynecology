# =============================================================================
# DEG: POP (n=12) vs. Controle (n=12) - GSE208261
# Script pronto para rodar no RStudio da usuaria, na pasta "POP+SUI FG".
#
# AJUSTE SO ESTA LINHA antes de rodar, para o caminho real da sua pasta:
setwd("E:/POP+SUI FG")
# =============================================================================

# Se DESeq2 nao estiver instalado, rode uma vez (pode demorar alguns minutos):
# if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install("DESeq2")

library(DESeq2)

# --- 1. Carregar as contagens brutas -----------------------------------------
counts <- read.delim("GSE208261_raw_counts_GRCh38.p13_NCBI.tsv", row.names = 1, check.names = FALSE)
counts <- round(as.matrix(counts))
mode(counts) <- "integer"

# --- 2. Definir os grupos (POP vs Controle), por amostra ---------------------
# Rotulos oficiais do GEO (GSE208261_family_1.soft) - 12 Controle + 12 POP.
sample_labels <- c(
  GSM6339911 = "Control_D1", GSM6339912 = "Control_D2", GSM6339913 = "Control_D3",
  GSM6339914 = "Control_D4", GSM6339915 = "Control_D5", GSM6339916 = "Control_D6",
  GSM6339917 = "Control_Y1", GSM6339918 = "Control_Y2", GSM6339919 = "Control_Y3",
  GSM6339920 = "Control_Y4", GSM6339921 = "Control_Y5", GSM6339922 = "Control_Y6",
  GSM6339923 = "POP_D1", GSM6339924 = "POP_D2", GSM6339925 = "POP_D3",
  GSM6339926 = "POP_D4", GSM6339927 = "POP_D5", GSM6339928 = "POP_D6",
  GSM6339929 = "POP_Y1", GSM6339930 = "POP_Y2", GSM6339931 = "POP_Y3",
  GSM6339932 = "POP_Y4", GSM6339933 = "POP_Y5", GSM6339934 = "POP_Y6"
)
label <- sample_labels[colnames(counts)]
condition <- factor(ifelse(grepl("^POP", label), "POP", "Control"), levels = c("Control", "POP"))

# --- 3. Filtro de baixa contagem ----------------------------------------------
# Mantem genes com contagem >=10 em pelo menos 12 das 24 amostras.
keep <- rowSums(counts >= 10) >= 12
counts <- counts[keep, ]
cat("Genes retidos apos filtro:", nrow(counts), "de", length(keep), "\n")

# --- 4. DESeq2 -----------------------------------------------------------------
coldata <- data.frame(row.names = colnames(counts), condition = condition)
dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~condition)
dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "POP", "Control"), alpha = 0.05)

res_df <- as.data.frame(res)
res_df$GeneID <- rownames(res_df)
res_df <- res_df[order(res_df$padj), ]

# --- 5. Salvar resultados -------------------------------------------------------
write.csv(res_df, "POP_GSE208261_24samples_DESeq2_completo.csv", row.names = FALSE)

sig <- subset(res_df, padj < 0.05 & abs(log2FoldChange) > 1)
write.csv(sig, "POP_GSE208261_24samples_DEG_sig.csv", row.names = FALSE)

cat("\nGenes testados:", nrow(res_df), "\n")
cat("DEGs significativos (padj<0.05, |log2FC|>1):", nrow(sig), "\n")
cat("  up em POP:", sum(sig$log2FoldChange > 0), "  down em POP:", sum(sig$log2FoldChange < 0), "\n")
