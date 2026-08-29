# =============================================================================
# Cruzamento SUI HUMANO (literatura, ja que nao ha dataset bruto publico) x POP
# HUMANO real (GSE208261, DESeq2 ja calculado).
#
# Fonte da lista de SUI: Chen B, Wen Y, Zhang Z, Guo Y, Warrington JA, Polan ML
# (2006) "Microarray analysis of differentially expressed genes in vaginal
# tissues from women with stress urinary incontinence compared with
# asymptomatic women", Human Reproduction 21(1):22-29. Microarray Affymetrix
# HG-U133, n=5 pares SUI x continentes, parede vaginal periuretral,
# pre-menopausa. 79 DEGs totais no artigo original - só consegui confirmar
# via busca os 13 genes nomeados explicitamente na literatura de revisão
# (Miličić et al. 2023) - NAO É a lista completa de 79, é um subconjunto
# citado. Mais: Tong X, Lang J, Zhu L (2010) Int Urogynecol J - GBA,
# microarray de tecido pelvico degenerado.
#
# CORRECAO: "SKALP/elafin" (nome usado na revisao) e o gene PI3 (peptidase
# inhibitor 3), NAO SLPI (secretory leukocyte protease inhibitor) - sao genes
# homologos mas DIFERENTES. Usei PI3 aqui, meu uso anterior de "SLPI" estava
# errado.
#
# Sem dataset bruto disponivel para SUI humano, esta lista e o retangulo de
# genes candidatos "conhecidos verdadeiros" contra os quais cruzo o POP -
# metodologia = candidate gene lookup, nao DEG-vs-DEG cego.
# =============================================================================

suppressMessages(library(org.Hs.eg.db))

sui_lit <- data.frame(
  Gene = c("PI3","COL17A1","PKP1","KRT16","DCN","BGN","BICD2","GRB2","STAT3",
           "APOE","GOSR1","FMOD","GBA"),
  Direcao_SUI_literatura = c("up","up","up","up","up","up","up","up","up",
                              "up","up","down","down"),
  Fonte = c(rep("Chen et al. 2006, Hum Reprod 21:22-29", 12), "Tong et al. 2010, Int Urogynecol J"),
  stringsAsFactors = FALSE
)

pop_full <- read.csv("results/POP_GSE208261_DESeq2_completo.csv")
pop_full$Human_Entrez_ID <- as.character(pop_full$Human_Entrez_ID)
pop_map <- suppressWarnings(select(org.Hs.eg.db, keys = pop_full$Human_Entrez_ID, keytype = "ENTREZID", columns = "SYMBOL"))
pop_sym <- merge(pop_full, pop_map, by.x = "Human_Entrez_ID", by.y = "ENTREZID", all.x = TRUE)

result <- merge(sui_lit, pop_sym[, c("SYMBOL","baseMean","log2FoldChange","padj")],
                 by.x = "Gene", by.y = "SYMBOL", all.x = TRUE)
result$Testado_no_POP <- !is.na(result$log2FoldChange)
result$Significativo_no_POP_FDR05 <- !is.na(result$padj) & result$padj < 0.05 & abs(result$log2FoldChange) > 0.5
result$Direcao_POP <- ifelse(is.na(result$log2FoldChange), NA, ifelse(result$log2FoldChange > 0, "up", "down"))
result$Concordante_com_literatura_SUI <- result$Direcao_SUI_literatura == result$Direcao_POP

write.csv(result, "results/SUI_humano_literatura_x_POP_real.csv", row.names = FALSE)
cat("=== Cruzamento SUI (literatura humana) x POP (dado real, GSE208261) ===\n")
print(result[, c("Gene","Direcao_SUI_literatura","Testado_no_POP","log2FoldChange","padj",
                  "Significativo_no_POP_FDR05","Direcao_POP","Concordante_com_literatura_SUI")])

cat("\nGenes SIGNIFICATIVOS no POP entre os candidatos de SUI da literatura:",
    sum(result$Significativo_no_POP_FDR05, na.rm=TRUE), "de", nrow(result), "\n")
