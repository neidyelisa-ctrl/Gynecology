# =============================================================================
# Script 5/4 (OPCIONAL) — versão "oficial", para rodar no SEU computador com
# internet normal (não neste sandbox, que bloqueia CRAN/Bioconductor/MSigDB).
#
# Os scripts 01-04 reimplementam GSEA do zero usando só as vias KEGG que já
# vêm dentro do pacote offline org.Hs.eg.db (uma anotação congelada, sem
# nomes de via legíveis, e sem Hallmark/Reactome/GO). Este script usa as
# ferramentas padrão da comunidade (fgsea + msigdbr) para o MESMO par de
# análises, mas com:
#   - a base COMPLETA e atualizada do MSigDB (Hallmark, KEGG, Reactome, GO,
#     etc. — não só KEGG), com nomes de via legíveis;
#   - fgsea, a implementação de referência do algoritmo (mais rápida e mais
#     testada que a reimplementação manual dos scripts 01-04);
#   - o mesmo desenho pareado (limma) para o POP e o mesmo painel de 60
#     genes do Chen 2006 para a lista pré-ranqueada do SUI.
#
# COMO RODAR:
#   1. Abra este arquivo no RStudio, no seu computador (com internet normal).
#   2. Rode a Seção 0 uma vez para instalar os pacotes.
#   3. Rode o script inteiro, ou por blocos, a partir da pasta deste projeto
#      (setwd() para a pasta POP_vs_SUI_GSEA_Chen2006/, ou abra-a como
#      projeto do RStudio).
# =============================================================================

## --- Seção 0: instalar pacotes (uma vez) -----------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("limma", "fgsea"), update = FALSE, ask = FALSE)
install.packages(c("msigdbr", "dplyr"))

suppressMessages({
  library(limma)
  library(fgsea)
  library(msigdbr)
  library(dplyr)
})

## --- Seção 1: DEG do POP (idêntico ao script 01) ---------------------------
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
coldata <- data.frame(row.names = colnames(expr),
                       tissue = factor(tissue, levels = c("NonPOP_site", "POP_site")),
                       individual = factor(individuals))

design <- model.matrix(~ individual + tissue, data = coldata)
fit <- eBayes(lmFit(as.matrix(expr), design))
res <- topTable(fit, coef = "tissuePOP_site", number = Inf, sort.by = "P")
res$Gene <- rownames(res)

deg <- subset(res, adj.P.Val < 0.05 & abs(logFC) > 1)
write.csv(deg, "results/GSE53868_DEG_logFC1_FDR05_fgsea.csv", row.names = FALSE)
cat("POP DEGs (|log2FC|>1, FDR<0.05):", nrow(deg), "\n")

## --- Seção 2: gene sets do MSigDB (Hallmark + KEGG + Reactome, humano) ----
## collection = "H" (Hallmark), "C2"/subcollection "CP:KEGG" e "CP:REACTOME"
## (msigdbr >= 7.5.1 usa collection/subcollection; versões mais antigas usam
## category/subcategory — ajuste se necessário conforme a versão instalada).
gs_hallmark <- msigdbr(species = "Homo sapiens", collection = "H")
gs_kegg     <- msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:KEGG")
gs_reactome <- msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:REACTOME")
gs_all <- bind_rows(gs_hallmark, gs_kegg, gs_reactome)
pathways <- split(gs_all$gene_symbol, gs_all$gs_name)
cat("Gene sets carregados (Hallmark+KEGG+Reactome):", length(pathways), "\n")

## --- Seção 3: GSEA "clássico" do POP (transcriptoma inteiro, ranking = t) -
ranks_pop <- setNames(res$t, res$Gene)
ranks_pop <- ranks_pop[!is.na(ranks_pop)]
ranks_pop <- sort(ranks_pop, decreasing = TRUE)

set.seed(42)
fgsea_pop <- fgsea(pathways = pathways, stats = ranks_pop, minSize = 5, maxSize = 500, eps = 0)
fgsea_pop <- fgsea_pop[order(fgsea_pop$pval), ]
fwrite_or_write_csv <- function(x, path) {
  x$leadingEdge <- vapply(x$leadingEdge, paste, collapse = "/", FUN.VALUE = character(1))
  write.csv(x, path, row.names = FALSE)
}
fwrite_or_write_csv(as.data.frame(fgsea_pop), "results/fgsea_POP_full_MSigDB.csv")
cat("POP (fgsea, MSigDB completo) — significativos FDR<0.25:",
    sum(fgsea_pop$padj < 0.25), "| FDR<0.05:", sum(fgsea_pop$padj < 0.05), "\n")

## --- Seção 4: GSEA preranked do Chen 2006 (painel de 60 genes) -----------
chen <- read.csv("data/chen2006_79genes.csv")
chen <- chen[!duplicated(chen$GeneSymbol), ]
p <- chen$RMA_pvalue
p[p <= 0] <- min(p[p > 0], na.rm = TRUE) / 2
score_chen <- ifelse(chen$Direction == "up", 1, -1) * -log10(p)
ranks_chen <- setNames(score_chen, chen$GeneSymbol)
ranks_chen <- sort(ranks_chen, decreasing = TRUE)

set.seed(123)
fgsea_chen <- fgsea(pathways = pathways, stats = ranks_chen, minSize = 2, maxSize = 500, eps = 0)
fgsea_chen <- fgsea_chen[order(fgsea_chen$pval), ]
fwrite_or_write_csv(as.data.frame(fgsea_chen), "results/fgsea_Chen2006_full_MSigDB.csv")
cat("Chen2006 (fgsea, MSigDB completo) — significativos FDR<0.25:",
    sum(fgsea_chen$padj < 0.25), "| FDR<0.05:", sum(fgsea_chen$padj < 0.05), "\n")

## --- Seção 5: vias compartilhadas (por gs_name, legível) -------------------
pop_sig_025  <- fgsea_pop$pathway[fgsea_pop$padj < 0.25]
chen_sig_025 <- fgsea_chen$pathway[fgsea_chen$padj < 0.25]
shared_025 <- intersect(pop_sig_025, chen_sig_025)
cat("\nVias compartilhadas (FDR<0.25, MSigDB completo):", length(shared_025), "\n")
if (length(shared_025) > 0) print(shared_025)

pop_sig_005  <- fgsea_pop$pathway[fgsea_pop$padj < 0.05]
chen_sig_005 <- fgsea_chen$pathway[fgsea_chen$padj < 0.05]
shared_005 <- intersect(pop_sig_005, chen_sig_005)
cat("Vias compartilhadas (FDR<0.05, MSigDB completo):", length(shared_005), "\n")
if (length(shared_005) > 0) print(shared_005)

shared_tbl <- data.frame(pathway = union(shared_025, shared_005))
if (nrow(shared_tbl) > 0) {
  shared_tbl <- merge(shared_tbl, fgsea_pop[, c("pathway","NES","pval","padj")], by = "pathway")
  shared_tbl <- merge(shared_tbl, fgsea_chen[, c("pathway","NES","pval","padj")], by = "pathway",
                       suffixes = c("_POP", "_Chen2006"))
  write.csv(shared_tbl, "results/fgsea_shared_pathways_full_MSigDB.csv", row.names = FALSE)
}

cat("\n=== FIM (versão fgsea + MSigDB completo) ===\n")
cat("Compare estes resultados com results/GSEA_shared_KEGG_FDR025.csv (versão\n")
cat("offline, só KEGG, gerada nos scripts 01-04) — o padrão qualitativo deve\n")
cat("ser semelhante (poucas ou nenhuma via nomeada compartilhada, dado o\n")
cat("painel pequeno do Chen 2006), mas os números exatos podem diferir por\n")
cat("causa da base de vias maior/mais atual e da implementação mais precisa\n")
cat("do fgsea.\n")
