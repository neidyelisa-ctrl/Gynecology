# =============================================================================
# Genes candidatos da via neurotrofica/neuropeptidica (literatura) testados
# contra os DEGs reais de SUI (GSE149072, 36h/72h) e POP (GSE208261)
#
# Lista curada a partir de:
#   - Shen et al. 2024, "Research progress on pelvic nerve in pelvic organ
#     prolapse and stress urinary incontinence: Systematic review"
#     (preprint, DOI 10.22541/au.171220777.79433562/v1)
#   - Masyhuroh et al. 2024, "Navigating the Complexities of Pelvic Organ
#     Prolapse and Stress Urinary Incontinence: A Review of Current
#     Literature" (Int J Scientific Advances, DOI 10.51542/ijscia.v5i6.40)
#
# Resultado real (rodado em 19/08/2026): ver results/candidate_neuro_genes_table.csv
# e a aba "Genes_candidatos" de results/Proposta_tema_neuro_POP_SUI.xlsx
#
# Pre-requisito: rodar analysis_SUI_POP_crosscomparison.R antes (gera os
# arquivos SUI_36hr_DESeq2_completo.csv, SUI_72hr_DESeq2_completo.csv e
# POP_GSE208261_DESeq2_completo.csv usados aqui).
# =============================================================================

if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install("org.Hs.eg.db", update = FALSE, ask = FALSE)
}
if (!requireNamespace("homologene", quietly = TRUE)) install.packages("homologene")

library(org.Hs.eg.db)

sui36 <- read.csv("results/SUI_36hr_DESeq2_completo.csv")
sui72 <- read.csv("results/SUI_72hr_DESeq2_completo.csv")
pop   <- read.csv("results/POP_GSE208261_DESeq2_completo.csv")

## Lista de genes candidatos: vias NGF/Trk/p75NTR, VIP/PACAP, NPY e receptores,
## familia GDNF/GFRA, endotelina, marcadores gliais/neuronais, Wnt/beta-catenina
## e inflamacao (IL-1/NLRP3) - todos citados nos dois artigos acima como
## associados a POP e/ou SUI via mecanismo de nervo pelvico ou inflamacao.
neuro_genes <- c("NGF","NTRK1","NGFR","VIP","ADCYAP1","NPY","NPY1R","NPY2R","NPY5R",
                  "BDNF","NTF3","GFRA1","GFRA2","GFRA3","GFRA4","EDNRA","EDNRB",
                  "S100B","UCHL1","GFAP","CTNNB1","IL1B","IL1A","IL1R1","NLRP3")

## Ortologos humano->rato via NCBI HomoloGene
orth <- homologene::homologene(neuro_genes, inTax = 9606, outTax = 10116)
colnames(orth)[1:2] <- c("Human_Symbol", "Rat_Symbol")

entrez_map <- select(org.Hs.eg.db, keys = neuro_genes, keytype = "SYMBOL", columns = "ENTREZID")

rows <- list()
for (g in neuro_genes) {
  rat_sym <- orth$Rat_Symbol[match(g, orth$Human_Symbol)]
  row_sui36 <- if (!is.na(rat_sym)) sui36[sui36$Rat_Gene_Symbol == rat_sym, ] else data.frame()
  row_sui72 <- if (!is.na(rat_sym)) sui72[sui72$Rat_Gene_Symbol == rat_sym, ] else data.frame()
  eid <- entrez_map$ENTREZID[entrez_map$SYMBOL == g][1]
  row_pop <- if (!is.na(eid)) pop[as.character(pop$Human_Entrez_ID) == eid, ] else data.frame()

  rows[[g]] <- data.frame(
    Gene_Humano = g,
    Ortologo_rato = ifelse(is.na(rat_sym), "N/D", rat_sym),
    SUI_36h_log2FC = if (nrow(row_sui36) > 0) round(row_sui36$log2FoldChange[1], 3) else NA,
    SUI_36h_padj   = if (nrow(row_sui36) > 0) row_sui36$padj[1] else NA,
    SUI_72h_log2FC = if (nrow(row_sui72) > 0) round(row_sui72$log2FoldChange[1], 3) else NA,
    SUI_72h_padj   = if (nrow(row_sui72) > 0) row_sui72$padj[1] else NA,
    POP_log2FC     = if (nrow(row_pop) > 0) round(row_pop$log2FoldChange[1], 3) else NA,
    POP_padj       = if (nrow(row_pop) > 0) row_pop$padj[1] else NA,
    stringsAsFactors = FALSE
  )
}

final <- do.call(rbind, rows)
final$Significativo_SUI36h <- ifelse(!is.na(final$SUI_36h_padj) & final$SUI_36h_padj < 0.05, "SIM", "nao")
final$Significativo_SUI72h <- ifelse(!is.na(final$SUI_72h_padj) & final$SUI_72h_padj < 0.05, "SIM", "nao")
final$Significativo_POP    <- ifelse(!is.na(final$POP_padj) & final$POP_padj < 0.05, "SIM", "nao")

write.csv(final, "results/candidate_neuro_genes_table.csv", row.names = FALSE)
cat("Genes candidatos testados:", nrow(final), "\n")
cat("Significativos em SUI 36h:", sum(final$Significativo_SUI36h == "SIM"), "\n")
cat("Significativos em SUI 72h:", sum(final$Significativo_SUI72h == "SIM"), "\n")
cat("Significativos em POP:", sum(final$Significativo_POP == "SIM"), "\n")
