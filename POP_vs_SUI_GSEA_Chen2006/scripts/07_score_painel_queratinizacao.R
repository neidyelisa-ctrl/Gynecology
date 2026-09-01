# =============================================================================
# Script 7 — Escore de painel (module score) para os 6 genes do eixo de
# queratinização, testado diretamente nas 12 pacientes do POP.
#
# POR QUÊ ESTE TESTE (proposta de "resultado satisfatório sem forçar"):
# GSEA contra milhares de vias KEGG/GO precisa de muitos genes para ter
# poder estatístico - por isso o GSEA preranked do Chen 2006 (script 03) deu
# um resultado fraco (painel de 59 genes é pequeno demais). Mas o eixo de
# queratinização já não é mais uma hipótese "descoberta nos próprios dados"
# (o que seria caça-resultado/forçar) - ele é uma hipótese EXTERNA,
# levantada de forma independente pela combinação de: (1) genes concordantes
# Chen2006 x POP (script 06), (2) GO enrichment nesses genes concordantes
# (script 06), e (3) um estudo publicado independente (Zhang et al. 2024,
# Exp Cell Res, scRNA-seq humano de SUI) que achou o mesmo eixo por um
# método totalmente diferente. Testar ESSE painel específico, pré-definido,
# diretamente nas 12 pacientes do POP é uma confirmação dirigida por
# hipótese - estatisticamente muito mais forte para um painel pequeno do que
# testar contra o universo inteiro do GO/KEGG.
#
# OS 6 GENES (todos concordantes em direção Chen2006 x POP; escolhidos por
# serem os genes que aparecem nos 6 termos GO mais significativos do
# cruzamento do script 06 - não escolhidos à mão):
#   KRT14, KRT16, KRT17, PKP1, S100A7, COL17A1
# (TP63 aparece no Chen 2006 mas não foi medido no array GSE53868 - excluído,
# não forçado.)
#
# MÉTODO: para cada paciente, calcula-se a diferença pareada (sítio do
# prolapso menos sítio sem prolapso) de cada um dos 6 genes; cada gene é
# ESCALADO (dividido pelo seu próprio desvio-padrão entre as 12 pacientes,
# SEM subtrair a média) para que nenhum gene de alta variância (ex.: S100A7,
# que varia muito mais que os outros) domine sozinho o escore, mas SEM
# apagar o deslocamento real do grupo. O escore do painel de cada paciente é
# a MÉDIA desses 6 valores escalados.
#
# ATENÇÃO - ARMADILHA ESTATÍSTICA JÁ ENCONTRADA E CORRIGIDA: a primeira
# versão deste script usava `scale()` padrão do R, que SUBTRAI a média de
# cada gene antes de dividir pelo desvio-padrão (z-score completo, não só
# escala). Como cada gene já é, por construção, centrado em zero entre as
# mesmas 12 pacientes, a média do escore do painel dava EXATAMENTE zero em
# qualquer painel, para qualquer conjunto de genes - um resultado nulo
# garantido pela matemática do teste, não um achado biológico. Corrigido
# aqui dividindo pelo desvio-padrão SEM centralizar.
#
# Testa-se se esse escore é
# consistentemente diferente de zero nas 12 pacientes com três testes
# (redundantes de propósito, para robustez):
#   - teste t pareado de uma amostra (paramétrico)
#   - teste de Wilcoxon signed-rank de uma amostra (não-paramétrico)
#   - permutação de fenótipo (sign-flip, 10.000 permutações - mesmo
#     princípio do script 02, sem depender de suposição de normalidade)
# =============================================================================

set.seed(2024)

## --- 1) Reconstruir a matriz de diferença pareada (mesmo método do script 02)
raw_lines <- readLines("data/GSE53868_series_matrix.txt")
start_row <- grep("^!series_matrix_table_begin", raw_lines) + 1
end_row   <- grep("^!series_matrix_table_end", raw_lines) - 1
expr <- read.delim("data/GSE53868_series_matrix.txt", skip = start_row - 1,
                    nrows = end_row - start_row, header = TRUE,
                    row.names = 1, check.names = FALSE, quote = "\"")

sample_title_line <- raw_lines[grep("^!Sample_title", raw_lines)]
sample_titles <- gsub('"', "", strsplit(sample_title_line, "\t")[[1]][-1])
individual_line <- raw_lines[grep("^!Sample_characteristics_ch1.*individual:", raw_lines)][1]
individuals <- gsub("individual: ", "", gsub('"', "", strsplit(individual_line, "\t")[[1]][-1]))
tissue <- ifelse(grepl("\\(POP site\\)", sample_titles), "POP_site", "NonPOP_site")

expr <- as.matrix(expr)
pop_cols <- colnames(expr)[tissue == "POP_site"]
nonpop_cols <- colnames(expr)[tissue == "NonPOP_site"]
pop_ind <- individuals[tissue == "POP_site"]
nonpop_ind <- individuals[tissue == "NonPOP_site"]
pop_cols <- pop_cols[match(nonpop_ind, pop_ind)]

D_all <- expr[, pop_cols] - expr[, nonpop_cols]   # genes x 12 pacientes, + = para cima no sitio do prolapso
colnames(D_all) <- nonpop_ind                      # nomear colunas pelo ID da paciente

## --- 2) Painel de 6 genes ---------------------------------------------------
panel <- c("KRT14", "KRT16", "KRT17", "PKP1", "S100A7", "COL17A1")
panel_present <- intersect(panel, rownames(D_all))
cat("Genes do painel presentes no array:", paste(panel_present, collapse = ", "),
    "(", length(panel_present), "de", length(panel), ")\n\n")
stopifnot(length(panel_present) == length(panel))  # os 6 devem estar presentes (já conferido)

D_panel <- D_all[panel_present, , drop = FALSE]
print(round(D_panel, 2))

## --- 3) Escore por paciente (média dos 6 genes, escalados por DP, sem centralizar) ---
gene_sd <- apply(D_panel, 1, sd)
z_panel <- D_panel / gene_sd  # escala por gene (divide por DP), preserva o deslocamento real
score <- colMeans(z_panel)
cat("\nEscore do painel por paciente (positivo = para cima no sítio do prolapso):\n")
print(round(sort(score, decreasing = TRUE), 3))

## --- 4) Testes estatísticos --------------------------------------------------
t_test <- t.test(score, mu = 0)
w_test <- wilcox.test(score, mu = 0)

n_perm <- 10000
n_pat <- length(score)
perm_scores <- numeric(n_perm)
for (i in seq_len(n_perm)) {
  signs <- sample(c(-1, 1), n_pat, replace = TRUE)
  D_perm <- sweep(D_panel, 2, signs, `*`)
  sd_perm <- apply(D_perm, 1, sd)
  z_perm <- D_perm / sd_perm
  perm_scores[i] <- mean(colMeans(z_perm))
}
obs_mean <- mean(score)
p_perm <- (sum(abs(perm_scores) >= abs(obs_mean)) + 1) / (n_perm + 1)

cat("\n=== Testes: escore do painel é diferente de zero nas 12 pacientes? ===\n")
cat("Média do escore:", round(obs_mean, 3), "\n")
cat("Teste t pareado (1 amostra): t =", round(t_test$statistic, 2),
    ", p =", format(t_test$p.value, digits = 4), "\n")
cat("Wilcoxon signed-rank: V =", w_test$statistic,
    ", p =", format(w_test$p.value, digits = 4), "\n")
cat("Permutação de fenótipo (sign-flip, 10.000 perms): p =", format(p_perm, digits = 4), "\n")

## --- 5) Salvar resultados + gráfico -----------------------------------------
out <- data.frame(paciente = names(score), escore_painel = round(score, 4))
out <- out[order(-out$escore_painel), ]
dir.create("results", showWarnings = FALSE)
write.csv(out, "results/07_escore_painel_queratinizacao_por_paciente.csv", row.names = FALSE)

resumo <- data.frame(
  metodo = c("t pareado (1 amostra)", "Wilcoxon signed-rank", "Permutacao de fenotipo (10000x)"),
  estatistica = c(round(t_test$statistic, 3), w_test$statistic, NA),
  pvalue = c(t_test$p.value, w_test$p.value, p_perm)
)
write.csv(resumo, "results/07_escore_painel_testes.csv", row.names = FALSE)

png("results/07_escore_painel_boxplot.png", width = 900, height = 700, res = 130)
par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
boxplot(score, ylab = "Escore do painel (z-score médio, 6 genes)",
        main = "Escore por paciente\n(POP site - non-POP site)", col = "lightblue")
stripchart(score, vertical = TRUE, add = TRUE, pch = 16, col = "darkblue")
abline(h = 0, lty = 2, col = "red")

matplot(t(D_panel), type = "b", pch = 1, lty = 1, col = 1:6,
        xlab = "Paciente", ylab = "Diferença pareada (log2, POP - nonPOP)",
        main = "Os 6 genes individualmente\n(1 linha por gene)", xaxt = "n")
axis(1, at = seq_len(ncol(D_panel)), labels = colnames(D_panel), cex.axis = 0.7)
abline(h = 0, lty = 2, col = "grey40")
legend("topright", legend = panel_present, col = 1:6, lty = 1, pch = 1, cex = 0.6, bty = "n")
dev.off()

cat("\nResultados salvos em results/07_escore_painel_*.csv e .png\n")
cat("=== FIM ===\n")
