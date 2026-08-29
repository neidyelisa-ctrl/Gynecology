# =============================================================================
# GSEA (Gene Set Enrichment Analysis) no dataset REAL de POP (GSE208261).
#
# Este script e para a USUARIA RODAR NO PROPRIO AMBIENTE (com internet), nao
# neste sandbox - GSEA precisa baixar os gene sets do MSigDB (Hallmark/KEGG/
# GO/Reactome), o que exige acesso ao vivo (via msigdbr ou clusterProfiler),
# ja confirmado bloqueado neste ambiente em tentativa anterior deste projeto.
#
# IMPORTANTE sobre por que isso NAO resolve o lado do SUI: GSEA exige uma
# lista RANQUEADA de TODOS os genes testados (com um escore continuo por
# gene - aqui, a estatistica do DESeq2), nao so os genes "destaque" que
# aparecem em texto de artigo. Os artigos de SUI humano que temos (Chen 2006,
# Tong 2010) so tem 13 genes citados na literatura de revisao - nao a tabela
# completa dos ~79 genes nem o ranking do transcriptoma inteiro. Ou seja, dá
# para rodar GSEA no POP (dado real completo), mas NAO dá para rodar GSEA
# comparavel no SUI sem conseguir a tabela suplementar completa do artigo
# original (verifique se o Chen 2006/Tong 2010 tem arquivo suplementar com
# a lista inteira - às vezes esses artigos antigos tem Excel supplementary
# que a busca não indexa bem).
# =============================================================================

if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
  BiocManager::install("clusterProfiler")
}
if (!requireNamespace("msigdbr", quietly = TRUE)) {
  install.packages("msigdbr")  # baixa os gene sets do MSigDB - precisa internet
}
library(clusterProfiler)
library(msigdbr)
library(org.Hs.eg.db)

## -----------------------------------------------------------------------
## 1) Lista RANQUEADA - TODOS os genes testados no POP, nao so os DEGs
## -----------------------------------------------------------------------
pop_full <- read.csv("results/POP_GSE208261_DESeq2_completo.csv")
pop_full$Human_Entrez_ID <- as.character(pop_full$Human_Entrez_ID)
pop_full <- pop_full[!is.na(pop_full$stat), ]  # remove genes sem estatistica (baixa contagem)

## Ranking recomendado: a estatistica de Wald do DESeq2 ("stat") - já combina
## magnitude e confiança da mudança, é o que clusterProfiler/fgsea esperam.
## Alternativa comum: sign(log2FoldChange) * -log10(pvalue).
ranked <- setNames(pop_full$stat, pop_full$Human_Entrez_ID)
ranked <- sort(ranked, decreasing = TRUE)

## -----------------------------------------------------------------------
## 2) Gene sets do MSigDB - Hallmark (50 vias amplas, bom ponto de partida)
##    e KEGG. category="H" = Hallmark; category="C2", subcategory="CP:KEGG"
##    = KEGG.
## -----------------------------------------------------------------------
hallmark <- msigdbr(species = "Homo sapiens", category = "H")
hallmark_t2g <- hallmark[, c("gs_name", "entrez_gene")]

kegg_sets <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG")
kegg_t2g <- kegg_sets[, c("gs_name", "entrez_gene")]

## -----------------------------------------------------------------------
## 3) Rodar GSEA
## -----------------------------------------------------------------------
gsea_hallmark <- GSEA(geneList = ranked, TERM2GENE = hallmark_t2g,
                       pvalueCutoff = 0.05, pAdjustMethod = "BH", seed = TRUE)
gsea_kegg <- GSEA(geneList = ranked, TERM2GENE = kegg_t2g,
                   pvalueCutoff = 0.05, pAdjustMethod = "BH", seed = TRUE)

print(as.data.frame(gsea_hallmark)[, c("ID","NES","pvalue","p.adjust")])
print(as.data.frame(gsea_kegg)[, c("ID","NES","pvalue","p.adjust")])

write.csv(as.data.frame(gsea_hallmark), "GSEA_POP_Hallmark.csv", row.names = FALSE)
write.csv(as.data.frame(gsea_kegg), "GSEA_POP_KEGG.csv", row.names = FALSE)

## Visualizar uma via especifica (troque pelo nome exato que aparecer na
## coluna ID do resultado, ex. "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"):
## library(enrichplot)
## gseaplot2(gsea_hallmark, geneSetID = "HALLMARK_...")
