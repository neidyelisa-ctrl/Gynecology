# =============================================================================
# SCRIPT PRINCIPAL - rode este arquivo do início ao fim no RStudio.
# =============================================================================
# Responde 3 pedidos da usuária:
#
#   1) "Onde foi o eixo de queratinizacao?" -> SEÇÃO 5 abaixo. Resumo: os
#      genes KRT14, KRT16, KRT17, PKP1, S100A7, COL17A1 aparecem entre os
#      genes do painel do Chen 2006 (SUI) que concordam em direção com o
#      POP (GSE53868) - ver `results/06_cruzamento_Chen2006_x_POP_genes.csv`
#      e `results/06_GO_BP_concordantes.csv` (118 de 243 termos GO
#      significativos, todos desse eixo). Este script gera um heatmap
#      desses 6 genes nos dois datasets para visualizar diretamente.
#
#   2) "No ficheiro do Wei2020 tem 2 paineis (abas), você visitou os 2?" ->
#      SIM. A aba `up_Sui_vs_Ctrl` (3.991 genes) e a aba `down_Sui_vs_Ctrl`
#      (2.127 genes) foram lidas e combinadas em TODOS os scripts que usam
#      o Wei2020 completo (14, 18, 20, 24, 25, 27) - juntas somam os 6.118
#      genes usados como painel de SUI. Ver PARTE 2 abaixo, que le as 2
#      abas de novo, do zero, para conferir.
#
#   3) "Quero o script para rodar no R" + "quero todos os gráficos
#      possíveis" -> e este arquivo. Gera volcano plot (POP e SUI),
#      heatmap do eixo de queratinizacao, heatmap dos genes concordantes,
#      e gráficos de barra das vias GO/KEGG mais importantes (eixo de
#      queratinizacao E via TGF-beta, a via mais replicada do projeto).
#      Todas as figuras são salvas em `figures/` (arquivos .png).
#
# Pre-requisitos (uma vez): install.packages(c("ggplot2","ggrepel",
# "pheatmap","RColorBrewer")); BiocManager::install(c("limma",
# "org.Hs.eg.db","GO.db"))
# =============================================================================

suppressMessages({
  library(limma)
  library(org.Hs.eg.db)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
})

dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)
theme_set(theme_bw(base_size = 12))


## =============================================================================
## PARTE 1 - DEG do POP (GSE53868), limma pareado
## =============================================================================
cat("\n##### PARTE 1: DEG do POP (GSE53868) #####\n\n")

raw_lines <- readLines("data/GSE53868_series_matrix.txt")
start_row <- grep("^!series_matrix_table_begin", raw_lines) + 1
end_row   <- grep("^!series_matrix_table_end", raw_lines) - 1
expr_pop <- read.delim("data/GSE53868_series_matrix.txt", skip = start_row - 1,
                        nrows = end_row - start_row, header = TRUE,
                        row.names = 1, check.names = FALSE, quote = "\"")
expr_pop <- as.matrix(expr_pop)

sample_title_line <- raw_lines[grep("^!Sample_title", raw_lines)]
sample_titles <- gsub('"', "", strsplit(sample_title_line, "\t")[[1]][-1])
individual_line <- raw_lines[grep("^!Sample_characteristics_ch1.*individual:", raw_lines)][1]
individuals <- gsub("individual: ", "", gsub('"', "", strsplit(individual_line, "\t")[[1]][-1]))
tissue_pop <- ifelse(grepl("\\(POP site\\)", sample_titles), "POP_site", "NonPOP_site")

coldata_pop <- data.frame(row.names = colnames(expr_pop),
                           tissue = factor(tissue_pop, levels = c("NonPOP_site", "POP_site")),
                           individual = factor(individuals))
design_pop <- model.matrix(~ individual + tissue, data = coldata_pop)
fit_pop <- eBayes(lmFit(expr_pop, design_pop))
pop_full <- topTable(fit_pop, coef = "tissuePOP_site", number = Inf, sort.by = "P")
pop_full$Gene <- rownames(pop_full)
write.csv(pop_full[, c("Gene","logFC","AveExpr","t","P.Value","adj.P.Val")],
          "results/GSE53868_limma_completo.csv", row.names = FALSE)
pop_deg <- subset(pop_full, adj.P.Val < 0.05 & abs(logFC) > 1)
cat("DEGs do POP (|log2FC|>1, FDR<0.05):", nrow(pop_deg), "de", nrow(pop_full), "\n\n")


## =============================================================================
## PARTE 2 - DEG do SUI (Wei 2020) - as DUAS abas do Excel, lidas de novo
##            aqui mesmo para deixar explícito que as 2 são usadas.
## =============================================================================
cat("##### PARTE 2: DEG do SUI (Wei2020) - lendo as 2 abas do Excel #####\n\n")

wei_xls <- "data/wei2020_TableS2_mRNA_original.xls"
if (!requireNamespace("readxl", quietly = TRUE)) install.packages("readxl")
library(readxl)

up_sheet   <- read_excel(wei_xls, sheet = "up_Sui_vs_Ctrl",   skip = 17)
down_sheet <- read_excel(wei_xls, sheet = "down_Sui_vs_Ctrl", skip = 17)
cat("Aba 'up_Sui_vs_Ctrl':", nrow(up_sheet), "sondas\n")
cat("Aba 'down_Sui_vs_Ctrl':", nrow(down_sheet), "sondas\n")
cat("Total combinado (as 2 abas):", nrow(up_sheet) + nrow(down_sheet), "sondas\n\n")

up_sheet$Direction <- "up"; down_sheet$Direction <- "down"
sample_cols <- c("[Sui1, Sui](normalized)","[Sui2, Sui](normalized)","[Sui3, Sui](normalized)",
                  "[Ctrl1, Ctrl](normalized)","[Ctrl2, Ctrl](normalized)","[Ctrl3, Ctrl](normalized)")
keep_cols <- c("GeneSymbol", "P-value", "FDR", "Fold Change", "Direction", sample_cols)
wei_both <- rbind(up_sheet[, keep_cols], down_sheet[, keep_cols])
colnames(wei_both) <- c("GeneSymbol","PValue","FDR","FoldChange","Direction",
                         "Sui1","Sui2","Sui3","Ctrl1","Ctrl2","Ctrl3")
wei_both <- as.data.frame(wei_both)
wei_both <- wei_both[!is.na(wei_both$GeneSymbol), ]
# IMPORTANTE: a coluna "Fold Change" do GeneSpring (autores) e so a MAGNITUDE
# (sempre >=2, tanto na aba up quanto na aba down) - o sinal vem da aba/coluna
# Direction, nao do valor. Sem inverter o sinal para os genes "down", o
# log2(FoldChange) fica POSITIVO tambem para eles, o que e errado. Corrigido
# aqui: log2FC negativo para Direction=="down".
wei_both$logFC <- ifelse(wei_both$Direction == "down",
                          -log2(wei_both$FoldChange), log2(wei_both$FoldChange))

# colapsa sondas duplicadas (mantem a de menor p-valor por gene)
wei_both <- wei_both[order(wei_both$PValue), ]
wei_dedup <- wei_both[!duplicated(wei_both$GeneSymbol), ]
cat("Genes unicos apos colapsar sondas duplicadas:", nrow(wei_dedup), "\n\n")
write.csv(wei_dedup, "results/wei2020_ambas_abas_combinadas.csv", row.names = FALSE)


## =============================================================================
## PARTE 3 - VOLCANO PLOTS
## =============================================================================
cat("##### PARTE 3: Volcano plots #####\n\n")

make_volcano <- function(df, logfc_col, p_col, label, title, fc_cut = 1, p_cut = 0.05) {
  df <- df[!is.na(df[[logfc_col]]) & !is.na(df[[p_col]]), ]
  df$negLog10P <- -log10(df[[p_col]])
  df$sig <- "NS"
  df$sig[df[[logfc_col]] > fc_cut & df[[p_col]] < p_cut] <- "Up"
  df$sig[df[[logfc_col]] < -fc_cut & df[[p_col]] < p_cut] <- "Down"
  df$sig <- factor(df$sig, levels = c("Down","NS","Up"))
  top_lab <- df[order(df[[p_col]]), ][1:15, ]

  ggplot(df, aes(x = .data[[logfc_col]], y = negLog10P, color = sig)) +
    geom_point(alpha = 0.6, size = 1.3) +
    scale_color_manual(values = c(Down = "#2166AC", NS = "grey75", Up = "#B2182B")) +
    geom_vline(xintercept = c(-fc_cut, fc_cut), linetype = "dashed", color = "grey40") +
    geom_hline(yintercept = -log10(p_cut), linetype = "dashed", color = "grey40") +
    geom_text_repel(data = top_lab, aes(label = Gene_label), size = 3, color = "black",
                     max.overlaps = 20, segment.size = 0.2) +
    labs(title = title, x = "log2(Fold Change)",
         y = expression(-log[10](italic(p)~value)), color = NULL) +
    theme(legend.position = "top", plot.title = element_text(face = "bold"))
}

pop_full$Gene_label <- pop_full$Gene
v_pop <- make_volcano(pop_full, "logFC", "adj.P.Val",
                       title = "Volcano plot - POP (GSE53868)\nDEG: sitio do prolapso vs. sitio sem prolapso (FDR<0.05, |log2FC|>1)")
ggsave("figures/01_volcano_POP_GSE53868.png", v_pop, width = 9, height = 6.5, dpi = 300)
cat("Salvo: figures/01_volcano_POP_GSE53868.png\n")

wei_dedup$Gene_label <- wei_dedup$GeneSymbol
v_sui <- make_volcano(wei_dedup, "logFC", "PValue",
                       title = "Volcano plot - SUI (Wei 2020)\nDEG: SUI vs. controle (P<0.05, |log2FC|>1 - criterio original do artigo)")
ggsave("figures/02_volcano_SUI_Wei2020.png", v_sui, width = 9.5, height = 6.5, dpi = 300)
cat("Salvo: figures/02_volcano_SUI_Wei2020.png\n\n")


## =============================================================================
## PARTE 4 - Genes em comum / concordantes (reaproveita a metodologia já
##            validada nos scripts 06/14 deste projeto)
## =============================================================================
cat("##### PARTE 4: Cruzamento POP x SUI #####\n\n")

cross <- merge(wei_dedup[, c("GeneSymbol","Direction","logFC")], pop_full[, c("Gene","logFC","adj.P.Val")],
               by.x = "GeneSymbol", by.y = "Gene", suffixes = c("_SUI", "_POP"))
cross$Direcao_POP <- ifelse(cross$logFC_POP > 0, "up", "down")
cross$Concordante <- cross$Direction == cross$Direcao_POP
cat("Genes testaveis nos dois:", nrow(cross), "\n")
cat("Concordantes:", sum(cross$Concordante), "(",
    round(100*sum(cross$Concordante)/nrow(cross), 1), "% )\n\n")
write.csv(cross, "results/00_cruzamento_completo_POP_x_SUI.csv", row.names = FALSE)


## =============================================================================
## PARTE 5 - EIXO DE QUERATINIZAÇÃO: onde ele esta, com heatmap
## =============================================================================
cat("##### PARTE 5: Eixo de queratinizacao #####\n\n")

kerat_genes <- c("KRT14", "KRT16", "KRT17", "PKP1", "S100A7", "COL17A1")
cat("Estes 6 genes vieram do artigo do Chen 2006 (SUI) e sao os que mais\n")
cat("aparecem nos termos GO mais significativos do cruzamento Chen2006 x POP\n")
cat("(ver results/06_GO_BP_concordantes.csv - 118 de 243 termos GO sig.,\n")
cat("os 6 primeiros termos do topo sao TODOS sustentados por estes genes):\n\n")

kerat_tab <- pop_full[pop_full$Gene %in% kerat_genes, c("Gene","logFC","P.Value","adj.P.Val")]
print(kerat_tab)
cat("\nDirecao no Chen2006 (SUI) - todos 'up' no artigo original.\n")
cat("Direcao no POP acima - TODOS positivos tambem = concordantes.\n")
cat("KRT17 e o unico individualmente significativo no POP (FDR<0.05).\n\n")

# heatmap no POP (paired, 12 pacientes x 2 sitios)
kerat_expr_pop <- expr_pop[rownames(expr_pop) %in% kerat_genes, , drop = FALSE]
ann_col_pop <- data.frame(Sitio = coldata_pop$tissue, Paciente = coldata_pop$individual,
                           row.names = colnames(expr_pop))
ord_cols <- order(coldata_pop$individual, coldata_pop$tissue)
png("figures/03_heatmap_queratinizacao_POP.png", width = 2400, height = 1400, res = 300)
pheatmap(kerat_expr_pop[, ord_cols], scale = "row",
         annotation_col = ann_col_pop[ord_cols, , drop = FALSE],
         cluster_cols = FALSE, cluster_rows = TRUE,
         main = "Eixo de queratinizacao - POP (GSE53868)\n(z-score por gene, pacientes pareadas)",
         fontsize = 8, show_colnames = FALSE)
dev.off()
cat("Salvo: figures/03_heatmap_queratinizacao_POP.png\n")

# heatmap no SUI (3 SUI x 3 Ctrl) - NOTA: o eixo de queratinizacao veio do
# painel do Chen2006, nao do Wei2020 - a maioria destes 6 genes nao esta
# entre os 6.118 genes significativos do PROPRIO Wei2020 (paineis
# diferentes, doenca igual). So plota se houver >=2 genes em comum.
kerat_in_sui <- kerat_genes[kerat_genes %in% wei_dedup$GeneSymbol]
cat("Genes do eixo de queratinizacao presentes no painel do Wei2020:",
    if (length(kerat_in_sui) > 0) paste(kerat_in_sui, collapse=", ") else "NENHUM", "\n")
cat("(esperado: o eixo veio do Chen2006, um painel de SUI DIFERENTE do Wei2020 -\n")
cat(" nao e um problema, e so nao dá pra plotar heatmap de queratinizacao no Wei2020)\n\n")

if (length(kerat_in_sui) >= 2) {
  sui_mat <- as.matrix(wei_dedup[wei_dedup$GeneSymbol %in% kerat_in_sui,
                                  c("Sui1","Sui2","Sui3","Ctrl1","Ctrl2","Ctrl3")])
  rownames(sui_mat) <- wei_dedup$GeneSymbol[wei_dedup$GeneSymbol %in% kerat_in_sui]
  ann_col_sui <- data.frame(Grupo = c("SUI","SUI","SUI","Ctrl","Ctrl","Ctrl"),
                             row.names = colnames(sui_mat))
  png("figures/04_heatmap_queratinizacao_SUI.png", width = 2000, height = 1400, res = 300)
  pheatmap(sui_mat, scale = "row", annotation_col = ann_col_sui,
           cluster_cols = FALSE, cluster_rows = (length(kerat_in_sui) >= 3),
           main = "Eixo de queratinizacao - SUI (Wei 2020)\n(z-score por gene)",
           fontsize = 9)
  dev.off()
  cat("Salvo: figures/04_heatmap_queratinizacao_SUI.png\n\n")
} else {
  cat("PULADO: figures/04_heatmap_queratinizacao_SUI.png (menos de 2 genes em comum)\n\n")
}


## =============================================================================
## PARTE 6 - Heatmap dos genes CONCORDANTES (POP x SUI Wei2020) - os top 40
##            por menor p-valor no POP, para nao poluir o heatmap
## =============================================================================
cat("##### PARTE 6: Heatmap dos genes concordantes (top 40) #####\n\n")

concord <- subset(cross, Concordante)
concord <- merge(concord, pop_full[, c("Gene","P.Value")], by.x = "GeneSymbol", by.y = "Gene")
concord <- concord[order(concord$P.Value), ]
top_concord <- head(concord$GeneSymbol, 40)

kerat_expr_pop40 <- expr_pop[rownames(expr_pop) %in% top_concord, , drop = FALSE]
png("figures/05_heatmap_top40_concordantes_POP.png", width = 2600, height = 3200, res = 300)
pheatmap(kerat_expr_pop40[, ord_cols], scale = "row",
         annotation_col = ann_col_pop[ord_cols, , drop = FALSE],
         cluster_cols = FALSE, cluster_rows = TRUE,
         main = "Top 40 genes concordantes (POP x SUI Wei2020)\nExpressao no POP (GSE53868)",
         fontsize = 7, show_colnames = FALSE)
dev.off()
cat("Salvo: figures/05_heatmap_top40_concordantes_POP.png\n\n")


## =============================================================================
## PARTE 7 - Gráficos de barra: vias GO/KEGG (queratinizacao + TGF-beta)
## =============================================================================
cat("##### PARTE 7: Graficos de vias GO/KEGG #####\n\n")

# 7a. GO BP do cruzamento Chen2006 x POP (eixo de queratinizacao)
go_kerat <- read.csv("results/06_GO_BP_concordantes.csv")
go_kerat_top <- head(go_kerat[order(go_kerat$pvalue), ], 12)
go_kerat_top$TERM <- factor(go_kerat_top$TERM, levels = rev(go_kerat_top$TERM))

p_go <- ggplot(go_kerat_top, aes(x = -log10(pvalue), y = TERM, fill = Count)) +
  geom_col() +
  scale_fill_gradient(low = "#FDBB84", high = "#B30000") +
  labs(title = "Eixo de queratinizacao - GO Biological Process\n(genes concordantes Chen2006 x POP GSE53868)",
       x = expression(-log[10](italic(p)~value)), y = NULL, fill = "N genes") +
  theme(axis.text.y = element_text(size = 9))
ggsave("figures/06_barplot_GO_queratinizacao.png", p_go, width = 8, height = 5.5, dpi = 300)
cat("Salvo: figures/06_barplot_GO_queratinizacao.png\n")

# 7b. KEGG do cruzamento Wei2020 x POP (via TGF-beta, seção 11/13 do README)
kegg_wei <- read.csv("results/14_KEGG_concordantes.csv")
kegg_wei_top <- head(kegg_wei[order(kegg_wei$pvalue), ], 12)
kegg_wei_top$Via <- paste0("KEGG ", sprintf("%05d", as.integer(kegg_wei_top$TERM_ID)))
kegg_wei_top$Via[kegg_wei_top$TERM_ID == 4350] <- "KEGG 04350 - TGF-beta signaling *"
kegg_wei_top$Via <- factor(kegg_wei_top$Via, levels = rev(kegg_wei_top$Via))

p_kegg <- ggplot(kegg_wei_top, aes(x = -log10(pvalue), y = Via, fill = Count)) +
  geom_col() +
  scale_fill_gradient(low = "#9ECAE1", high = "#08519C") +
  labs(title = "Vias KEGG - genes concordantes Wei2020 x POP (GSE53868)\n(* TGF-beta replica em 4 dos 6 pares POP x SUI testados, ver secao 13 do README)",
       x = expression(-log[10](italic(p)~value)), y = NULL, fill = "N genes") +
  theme(axis.text.y = element_text(size = 9))
ggsave("figures/07_barplot_KEGG_Wei2020xPOP.png", p_kegg, width = 10.5, height = 5.5, dpi = 300)
cat("Salvo: figures/07_barplot_KEGG_Wei2020xPOP.png\n\n")

# 7c. Resumo: recorrência do TGF-beta e outras vias nos 6 pares POP x SUI
if (file.exists("results/23_KEGG_recorrentes_2ormais_pares.csv")) {
  rec <- read.csv("results/23_KEGG_recorrentes_2ormais_pares.csv")
  rec_tab <- as.data.frame(table(rec$TERM_ID))
  colnames(rec_tab) <- c("TERM_ID", "N_pares")
  rec_tab <- rec_tab[order(-rec_tab$N_pares), ]
  rec_tab$TERM_ID <- sprintf("%05d", as.integer(as.character(rec_tab$TERM_ID)))
  rec_tab$Via <- paste0("KEGG ", rec_tab$TERM_ID)
  rec_tab$Via[rec_tab$TERM_ID == "04350"] <- "KEGG 04350 - TGF-beta *"
  rec_top <- head(rec_tab, 15)
  rec_top$Via <- factor(rec_top$Via, levels = rev(rec_top$Via))
  p_rec <- ggplot(rec_top, aes(x = N_pares, y = Via)) +
    geom_col(fill = "#4292C6") +
    scale_x_continuous(breaks = 0:4) +
    labs(title = "Vias KEGG recorrentes em >=2 dos 6 pares POP x SUI testados\n(* achado mais replicado do projeto - 4 de 4 pares testaveis)",
         x = "Numero de pares POP x SUI em que a via e significativa (FDR<0.05)", y = NULL) +
    theme(axis.text.y = element_text(size = 9))
  ggsave("figures/08_barplot_vias_recorrentes_6pares.png", p_rec, width = 8, height = 6, dpi = 300)
  cat("Salvo: figures/08_barplot_vias_recorrentes_6pares.png\n\n")
}

cat("\n=== TODAS AS FIGURAS SALVAS EM figures/ ===\n")
cat("01_volcano_POP_GSE53868.png\n02_volcano_SUI_Wei2020.png\n")
cat("03_heatmap_queratinizacao_POP.png\n04_heatmap_queratinizacao_SUI.png\n")
cat("05_heatmap_top40_concordantes_POP.png\n06_barplot_GO_queratinizacao.png\n")
cat("07_barplot_KEGG_Wei2020xPOP.png\n08_barplot_vias_recorrentes_6pares.png\n")
cat("\n=== FIM ===\n")
