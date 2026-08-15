# =============================================================================
# Analise bioinformatica: Prolapso de Orgaos Pelvicos (POP) e
# Incontinencia Urinaria de Esforco (SUI) - datasets publicos do GEO
#
# Fonte da lista de datasets: data/POP_SUI_GEO_datasets.xlsx
# Como usar: abra este arquivo no RStudio (File > Open File...) e execute
# por blocos (Ctrl+Enter). Ajuste GSE_ID na secao 1 para trocar de dataset.
# =============================================================================

## -----------------------------------------------------------------------
## 1. Pacotes e configuracao
## -----------------------------------------------------------------------

required_bioc <- c("GEOquery", "limma", "Biobase", "clusterProfiler",
                    "org.Hs.eg.db", "AnnotationDbi", "enrichplot")
required_cran <- c("ggplot2", "pheatmap", "dplyr", "readxl", "openxlsx",
                    "RColorBrewer", "tibble")

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
missing_bioc <- required_bioc[!sapply(required_bioc, requireNamespace, quietly = TRUE)]
if (length(missing_bioc) > 0) BiocManager::install(missing_bioc, update = FALSE, ask = FALSE)

missing_cran <- required_cran[!sapply(required_cran, requireNamespace, quietly = TRUE)]
if (length(missing_cran) > 0) install.packages(missing_cran)

library(GEOquery)
library(limma)
library(Biobase)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(pheatmap)
library(dplyr)
library(readxl)
library(openxlsx)

## Diretorios de saida
dir.create("data/geo_cache", showWarnings = FALSE, recursive = TRUE)
dir.create("results", showWarnings = FALSE, recursive = TRUE)
dir.create("results/figures", showWarnings = FALSE, recursive = TRUE)

## -----------------------------------------------------------------------
## 2. Ler a lista de datasets (gerada pela busca de bioinformatica)
## -----------------------------------------------------------------------

datasets <- read_excel("data/POP_SUI_GEO_datasets.xlsx", sheet = "Datasets")
print(datasets[, c("Condição", "Accession", "Título")])

## Dataset alvo desta analise - troque para qualquer GSE valido da lista acima
## (ex.: "GSE53868", "GSE12852", "GSE28660")
GSE_ID <- "GSE53868"

## -----------------------------------------------------------------------
## 3. Download do dataset a partir do GEO
## -----------------------------------------------------------------------
## NOTA: esta etapa requer acesso a internet a partir do RStudio (o download
## nao foi executado neste ambiente porque o acesso a ncbi.nlm.nih.gov esta
## bloqueado pelo proxy da sessao). Rode este script localmente no RStudio.

gse_list <- getGEO(GSE_ID, GSEMatrix = TRUE, AnnotGPL = TRUE,
                    destdir = "data/geo_cache")
eset <- gse_list[[1]]

cat("Amostras:", ncol(eset), "| Sondas/genes:", nrow(eset), "\n")
print(head(pData(eset)))

## -----------------------------------------------------------------------
## 4. Definir os grupos (POP vs. controle)
## -----------------------------------------------------------------------
## Ajuste esta linha conforme a coluna de metadados real do dataset
## (varia entre series - inspecione pData(eset) para achar a coluna certa,
## ex.: "disease state:ch1", "tissue:ch1", "group:ch1")

pd <- pData(eset)
print(colnames(pd))

## Exemplo generico - EDITAR conforme a coluna real encontrada acima:
group_col <- grep("disease|group|status|state", colnames(pd),
                   ignore.case = TRUE, value = TRUE)[1]
stopifnot(!is.na(group_col))
cat("Usando coluna de grupo:", group_col, "\n")

group <- factor(ifelse(grepl("prolapse|POP|case", pd[[group_col]],
                              ignore.case = TRUE), "POP", "Control"))
table(group)

## -----------------------------------------------------------------------
## 5. Controle de qualidade
## -----------------------------------------------------------------------

exprs_mat <- exprs(eset)
if (max(exprs_mat, na.rm = TRUE) > 100) {
  exprs_mat <- log2(exprs_mat + 1)
}

pca <- prcomp(t(exprs_mat), scale. = TRUE)
pca_df <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], Grupo = group)

p_pca <- ggplot(pca_df, aes(PC1, PC2, color = Grupo)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(title = paste("PCA -", GSE_ID))
ggsave("results/figures/pca_plot.png", p_pca, width = 6, height = 5, dpi = 300)

## -----------------------------------------------------------------------
## 6. Expressao diferencial (limma)
## -----------------------------------------------------------------------

design <- model.matrix(~ group)
fit <- lmFit(exprs_mat, design)
fit <- eBayes(fit)

deg_table <- topTable(fit, coef = 2, number = Inf, adjust.method = "BH") %>%
  tibble::rownames_to_column("Probe_ID") %>%
  arrange(adj.P.Val)

deg_sig <- deg_table %>% filter(adj.P.Val < 0.05, abs(logFC) > 1)
cat("Genes diferencialmente expressos (padj<0.05, |logFC|>1):", nrow(deg_sig), "\n")

write.csv(deg_table, file.path("results", paste0(GSE_ID, "_DEG_completo.csv")),
          row.names = FALSE)
write.csv(deg_sig, file.path("results", paste0(GSE_ID, "_DEG_significativos.csv")),
          row.names = FALSE)

## Volcano plot
deg_table$sig <- ifelse(deg_table$adj.P.Val < 0.05 & abs(deg_table$logFC) > 1,
                         "Significativo", "Nao significativo")
p_volcano <- ggplot(deg_table, aes(logFC, -log10(adj.P.Val), color = sig)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("grey70", "firebrick")) +
  theme_minimal() +
  labs(title = paste("Volcano plot -", GSE_ID), color = NULL)
ggsave("results/figures/volcano_plot.png", p_volcano, width = 6, height = 5, dpi = 300)

## Heatmap dos top 50 genes
top50 <- head(deg_sig$Probe_ID, 50)
if (length(top50) > 1) {
  pheatmap(exprs_mat[top50, ], scale = "row",
           annotation_col = data.frame(Grupo = group, row.names = colnames(exprs_mat)),
           show_rownames = FALSE,
           filename = "results/figures/heatmap_top50.png",
           width = 7, height = 8)
}

## -----------------------------------------------------------------------
## 7. Enriquecimento funcional (GO / KEGG)
## -----------------------------------------------------------------------

gene_symbols <- deg_sig$Gene.symbol
if (is.null(gene_symbols)) gene_symbols <- deg_sig$Probe_ID
gene_symbols <- unique(na.omit(gene_symbols))

if (length(gene_symbols) > 5) {
  ego <- enrichGO(gene = gene_symbols, OrgDb = org.Hs.eg.db, keyType = "SYMBOL",
                   ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05)
  write.csv(as.data.frame(ego), file.path("results", paste0(GSE_ID, "_GO_BP.csv")),
            row.names = FALSE)

  p_go <- dotplot(ego, showCategory = 15) +
    labs(title = paste("GO Biological Process -", GSE_ID))
  ggsave("results/figures/go_dotplot.png", p_go, width = 8, height = 6, dpi = 300)
}

## -----------------------------------------------------------------------
## 8. Exportar resumo final em Excel
## -----------------------------------------------------------------------

wb <- createWorkbook()
addWorksheet(wb, "DEG_significativos")
writeData(wb, "DEG_significativos", deg_sig)
addWorksheet(wb, "DEG_completo")
writeData(wb, "DEG_completo", deg_table)
if (exists("ego")) {
  addWorksheet(wb, "GO_enrichment")
  writeData(wb, "GO_enrichment", as.data.frame(ego))
}
saveWorkbook(wb, file.path("results", paste0(GSE_ID, "_resultado_analise.xlsx")),
             overwrite = TRUE)

cat("\nAnalise concluida. Resultados em: results/\n")

## -----------------------------------------------------------------------
## 9. Para rodar outro dataset da lista
## -----------------------------------------------------------------------
## Basta trocar GSE_ID na secao 3 (ex.: "GSE12852", "GSE28660", "GSE151202")
## e reexecutar o script do inicio. Para GSE151202 (single-cell RNA-seq),
## use os pacotes Seurat/SingleCellExperiment em vez de limma - o desenho
## acima e voltado para dados de microarray/RNA-seq bulk.
