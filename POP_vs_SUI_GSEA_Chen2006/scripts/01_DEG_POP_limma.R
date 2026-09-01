# =============================================================================
# Script 1/4 — DEG do POP (GSE53868) com limma, desenho pareado
#
# Dataset: GSE53868 — "Micro-array analysis of the anterior vaginal wall from
# premenopausal patients with Pelvic Organ Prolapse (POP)" (Kerkhof et al.,
# VU University Medical Centre). 12 mulheres com POP, biópsia PAREADA por
# paciente: tecido do sítio do prolapso ("POP site") vs tecido da mesma
# paciente fora do prolapso ("non-POP site", precervical). Plataforma Agilent
# 4x44K (GPL18142). Os dados já vêm como matrix normalizada (log2), com
# símbolos de gene como ID_REF — não precisa de anotação de plataforma
# separada (conferido diretamente no arquivo).
#
# Critério de DEG pedido: |log2FC| > 1 e FDR (Benjamini-Hochberg) < 0.05.
#
# Pacote: limma (Bioconductor) — o padrão para microarray (diferente de
# DESeq2/edgeR, que são para contagens de RNA-seq). moderated t-test com
# eBayes (shrinkage de variância), bloco por paciente ("individual") para
# respeitar o desenho pareado.
# =============================================================================

suppressMessages(library(limma))

## --- 1) Ler o series matrix (formato texto padrão do GEO) -------------------
raw_lines <- readLines("data/GSE53868_series_matrix.txt")
start_row <- grep("^!series_matrix_table_begin", raw_lines) + 1
end_row   <- grep("^!series_matrix_table_end", raw_lines) - 1

expr <- read.delim("data/GSE53868_series_matrix.txt", skip = start_row - 1,
                    nrows = end_row - start_row, header = TRUE,
                    row.names = 1, check.names = FALSE, quote = "\"")

cat("Matriz de expressão:", nrow(expr), "genes x", ncol(expr), "amostras\n")

## --- 2) Metadados das amostras (do cabeçalho do series matrix) -------------
sample_title_line <- raw_lines[grep("^!Sample_title", raw_lines)]
sample_titles <- gsub('"', "", strsplit(sample_title_line, "\t")[[1]][-1])

individual_line <- raw_lines[grep("^!Sample_characteristics_ch1.*individual:", raw_lines)][1]
individuals <- gsub("individual: ", "", gsub('"', "", strsplit(individual_line, "\t")[[1]][-1]))

tissue <- ifelse(grepl("\\(POP site\\)", sample_titles), "POP_site", "NonPOP_site")
coldata <- data.frame(row.names = colnames(expr),
                       tissue = factor(tissue, levels = c("NonPOP_site", "POP_site")),
                       individual = factor(individuals))
stopifnot(all(rownames(coldata) == colnames(expr)))
print(table(coldata$tissue))

## --- 3) limma pareado: design ~ individual + tissue -------------------------
## Positivo em "tissue" = para cima no sítio do prolapso (lado afetado);
## negativo = para cima no sítio sem prolapso (lado não afetado, mesma
## paciente). Convenção usada em todo este pipeline: positivo = "para cima
## no lado mais afetado/doente".
design <- model.matrix(~ individual + tissue, data = coldata)
fit <- eBayes(lmFit(as.matrix(expr), design))

res <- topTable(fit, coef = "tissuePOP_site", number = Inf, sort.by = "P")
res$Gene <- rownames(res)
res <- res[, c("Gene", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val")]

dir.create("results", showWarnings = FALSE)
write.csv(res, "results/GSE53868_limma_completo.csv", row.names = FALSE)
cat("\nTabela completa (todos os genes testados) salva em results/GSE53868_limma_completo.csv\n")

## --- 4) DEGs: |log2FC| > 1 e FDR < 0.05 -------------------------------------
deg <- subset(res, adj.P.Val < 0.05 & abs(logFC) > 1)
deg <- deg[order(deg$adj.P.Val), ]
write.csv(deg, "results/GSE53868_DEG_logFC1_FDR05.csv", row.names = FALSE)

cat("\n=== DEGs do POP (GSE53868), |log2FC| > 1 e FDR < 0.05 ===\n")
cat(nrow(deg), "de", nrow(res), "genes testados\n")
cat(" - para cima no sítio do prolapso:", sum(deg$logFC > 0), "\n")
cat(" - para baixo no sítio do prolapso:", sum(deg$logFC < 0), "\n")
cat("\nEsta tabela (results/GSE53868_DEG_logFC1_FDR05.csv) é o ranked list completo\n")
cat("(results/GSE53868_limma_completo.csv) que alimenta o script 02 (GSEA clássico).\n")
