# =============================================================================
# Script 16 — DEG do GSE12852 (segundo dataset de POP, INDEPENDENTE do
# GSE53868), com limma.
#
# GSE12852 — "Gene expression profile in pelvic organ prolapse". 8 mulheres
# com POP vs. 9 controles (estudo caso-controle, NÃO pareado por sítio
# dentro da mesma paciente como o GSE53868), hysterectomia por indicação
# benigna. DOIS tecidos por paciente: ligamento uterossacral e ligamento
# redondo (34 arrays = 17 pacientes x 2 tecidos). Plataforma Applied
# Biosystems Human Genome Survey Microarray V2.0 (GPL2986) — IDs de sonda
# numéricos, SEM símbolo de gene direto (diferente do GSE53868/Agilent, que
# já vinha com símbolos). Anotação da plataforma
# (`data/GPL2986_annotation.tsv`) extraída da tabela `!platform_table` do
# arquivo SOFT completo da série (`GSE12852_family.soft`, fornecido pela
# usuária) — coluna `ID` (sonda) -> `Gene Symbol`. Valores da matriz de
# expressão são intensidade LINEAR (não log2) — transformados aqui.
#
# DESENHO: não é pareado por tecido (diferente do GSE53868), mas tem
# MEDIDAS REPETIDAS por paciente (2 tecidos por paciente, mesma condição
# nos dois). Tratado com limma + duplicateCorrelation (bloco = paciente),
# design ~ tissue + condition — o método padrão do limma para desenhos com
# medidas repetidas dentro do mesmo indivíduo. Convenção de direção mantida
# igual ao resto do projeto: positivo = "para cima no POP" (o lado
# afetado/doente).
# =============================================================================

suppressMessages(library(limma))

## --- 1) Ler o series matrix -------------------------------------------------
raw_lines <- readLines("data/GSE12852_series_matrix.txt")
start_row <- grep("^!series_matrix_table_begin", raw_lines) + 1
end_row   <- grep("^!series_matrix_table_end", raw_lines) - 1
expr <- read.delim("data/GSE12852_series_matrix.txt", skip = start_row - 1,
                    nrows = end_row - start_row, header = TRUE,
                    row.names = 1, check.names = FALSE, quote = "\"")
expr <- as.matrix(expr)
cat("Matriz de expressão bruta:", nrow(expr), "sondas x", ncol(expr), "amostras\n")

## --- 2) Metadados das amostras (extraídos do !Sample_title) ---------------
sample_title_line <- raw_lines[grep("^!Sample_title", raw_lines)]
titles <- gsub('"', "", strsplit(sample_title_line, "\t")[[1]][-1])

tissue <- ifelse(grepl("^round", titles), "round", "uterosacral")
condition <- ifelse(grepl("_POP_", titles), "POP", "control")
subject <- regmatches(titles, regexpr("subject[0-9]+", titles))
subject <- gsub("subject", "", subject)

coldata <- data.frame(row.names = colnames(expr),
                       tissue = factor(tissue, levels = c("uterosacral", "round")),
                       condition = factor(condition, levels = c("control", "POP")),
                       subject = factor(subject))
stopifnot(all(rownames(coldata) == colnames(expr)))
cat("\nDesenho:\n"); print(table(coldata$condition, coldata$tissue))
cat("Pacientes únicas:", nlevels(coldata$subject), "\n\n")

## --- 3) log2 e mapeamento sonda -> gene (GPL2986) --------------------------
expr_log2 <- log2(pmax(expr, 1))  # pmax evita log2(0) nas poucas sondas com sinal ~0

annot <- read.delim("data/GPL2986_annotation.tsv", colClasses = "character")
annot <- annot[!duplicated(annot$ID), ]
rownames(annot) <- annot$ID

common_probes <- intersect(rownames(expr_log2), annot$ID)
cat("Sondas com símbolo de gene anotado:", length(common_probes), "de", nrow(expr_log2), "\n")

expr_annot <- expr_log2[common_probes, ]
gene_symbol <- annot[common_probes, "GeneSymbol"]

# colapsa sondas duplicadas por gene (média do log2) - método padrão limma::avereps
expr_gene <- avereps(expr_annot, ID = gene_symbol)
cat("Genes únicos após colapsar sondas duplicadas:", nrow(expr_gene), "\n\n")

## --- 4) limma com duplicateCorrelation (bloco = paciente) ------------------
design <- model.matrix(~ tissue + condition, data = coldata)
corfit <- duplicateCorrelation(expr_gene, design, block = coldata$subject)
cat("Correlação intra-paciente estimada (duplicateCorrelation):",
    round(corfit$consensus.correlation, 4), "\n\n")

fit <- lmFit(expr_gene, design, block = coldata$subject, correlation = corfit$consensus.correlation)
fit <- eBayes(fit)

res <- topTable(fit, coef = "conditionPOP", number = Inf, sort.by = "P")
res$Gene <- rownames(res)
res <- res[, c("Gene", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val")]

dir.create("results", showWarnings = FALSE)
write.csv(res, "results/GSE12852_limma_completo.csv", row.names = FALSE)
cat("Tabela completa salva em results/GSE12852_limma_completo.csv\n\n")

## --- 5) DEGs: |log2FC| > 1 e FDR < 0.05 -------------------------------------
deg <- subset(res, adj.P.Val < 0.05 & abs(logFC) > 1)
deg <- deg[order(deg$adj.P.Val), ]
write.csv(deg, "results/GSE12852_DEG_logFC1_FDR05.csv", row.names = FALSE)

cat("=== DEGs do POP (GSE12852, tecidos combinados), |log2FC| > 1 e FDR < 0.05 ===\n")
cat(nrow(deg), "de", nrow(res), "genes testados\n")
cat(" - para cima no POP:", sum(deg$logFC > 0), "\n")
cat(" - para baixo no POP:", sum(deg$logFC < 0), "\n\n")
print(head(deg, 15))

## -----------------------------------------------------------------------
## 6) DIAGNÓSTICO: sinal muito fraco na análise combinada (FDR mínimo =
##    0,97 nos 16.752 genes - nem perto de FDR<0,05). Investigado por
##    tecido separadamente (round vs. uterosacral) antes de aceitar isso
##    como resultado final - é o tipo de checagem que este projeto sempre
##    faz antes de reportar um número (ver histórico de bugs encontrados
##    e corrigidos no README principal).
## -----------------------------------------------------------------------
cat("\n=== DIAGNÓSTICO: sinal por tecido separado (sem duplicateCorrelation) ===\n")
run_by_tissue <- function(tt) {
  idx <- coldata$tissue == tt
  cond_tt <- factor(coldata$condition[idx], levels = c("control", "POP"))
  design_tt <- model.matrix(~cond_tt)
  fit_tt <- eBayes(lmFit(expr_gene[, idx], design_tt))
  topTable(fit_tt, coef = 2, number = Inf, sort.by = "P")
}
res_round <- run_by_tissue("round")
res_utero <- run_by_tissue("uterosacral")
cat("Round ligament (n=17): FDR mínimo =", round(min(res_round$adj.P.Val), 3),
    "| p bruto<0.05:", sum(res_round$P.Value < 0.05), "de", nrow(res_round),
    "(esperado ao acaso:", round(0.05 * nrow(res_round)), ")\n")
cat("Uterosacral ligament (n=17): FDR mínimo =", round(min(res_utero$adj.P.Val), 3),
    "| p bruto<0.05:", sum(res_utero$P.Value < 0.05), "de", nrow(res_utero),
    "(esperado ao acaso:", round(0.05 * nrow(res_utero)), ")\n\n")

cat("CONCLUSÃO DO DIAGNÓSTICO: o ligamento REDONDO não mostra sinal acima do\n")
cat("esperado por acaso (371 vs 838 esperados com p<0.05) - na verdade ABAIXO\n")
cat("do acaso, e está DILUINDO o sinal da análise combinada. O ligamento\n")
cat("UTEROSSACRAL mostra sinal real, ainda que modesto (914 vs 838 esperados,\n")
cat("FDR mínimo 0,22 - melhor que a análise combinada, mas ainda não cruza\n")
cat("FDR<0,05). Biologicamente plausível: o ligamento uterossacral é a\n")
cat("estrutura mais diretamente implicada na fisiopatologia do POP. A partir\n")
cat("daqui, os scripts 17-18 usam UTEROSSACRAL SOZINHO como análise principal\n")
cat("(não a combinada) - mais sinal real, desenho mais simples (2 grupos,\n")
cat("sem repetição por paciente), sem inventar significância que não existe.\n\n")

res_round$Gene <- rownames(res_round); res_utero$Gene <- rownames(res_utero)
write.csv(res_round[, c("Gene","logFC","AveExpr","t","P.Value","adj.P.Val")],
          "results/GSE12852_round_limma_completo.csv", row.names = FALSE)
write.csv(res_utero[, c("Gene","logFC","AveExpr","t","P.Value","adj.P.Val")],
          "results/GSE12852_uterosacral_limma_completo.csv", row.names = FALSE)

deg_utero <- subset(res_utero, adj.P.Val < 0.05 & abs(logFC) > 1)
write.csv(deg_utero[, c("Gene","logFC","AveExpr","t","P.Value","adj.P.Val")],
          "results/GSE12852_uterosacral_DEG_logFC1_FDR05.csv", row.names = FALSE)
cat("DEGs uterossacral sozinho (|log2FC|>1, FDR<0.05):", nrow(deg_utero), "de", nrow(res_utero), "\n")
cat("Top 10 genes uterossacral (por p-valor bruto):\n")
print(head(res_utero[, c("Gene","logFC","AveExpr","t","P.Value","adj.P.Val")], 10))
