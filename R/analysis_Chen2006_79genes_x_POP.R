# =============================================================================
# Cruzamento COMPLETO: Chen et al. 2006 (79 DEGs de SUI humano, Tables II/III
# do PDF que a usuaria enviou) x POP humano real (GSE208261, DESeq2 ja
# calculado).
#
# Chen B, Wen Y, Zhang Z, Guo Y, Warrington JA, Polan ML (2006) "Microarray
# analysis of differentially expressed genes in vaginal tissues from women
# with stress urinary incontinence compared with asymptomatic women",
# Human Reproduction 21(1):22-29. Microarray Affymetrix HG-U133A, n=5 pares
# SUI x continentes, parede vaginal periuretral, fase secretora do ciclo.
# Tabelas II (39 up) e III (40 down) = 79 genes totais.
#
# CURADORIA: extrai do PDF, mapeei para simbolo HGNC atual. De 79 linhas nas
# tabelas, consegui mapear com confianca 55 para simbolo de gene atual -
# excluidas ~24 entradas ambiguas ("hypothetical protein FLJxxxxx", "KIAAxxxx
# protein" sem nome atualizado, "Zinc finger protein" generico, "Minor
# histocompatibility antigen HA-8", probe duplicado de KCNN3) porque nao da
# pra mapear com seguranca para um symbol atual sem o accession do probeset
# verificado ao vivo (blast/annotation, bloqueado neste sandbox).
# =============================================================================

suppressMessages(library(org.Hs.eg.db))

chen <- read.csv("data/chen2006_79genes.csv", stringsAsFactors = FALSE)
cat("Genes de Chen 2006 curados com simbolo atual:", nrow(chen),
    "(", sum(chen$Direction=="up"), "up,", sum(chen$Direction=="down"), "down )\n")

pop_full <- read.csv("results/POP_GSE208261_DESeq2_completo.csv")
pop_full$Human_Entrez_ID <- as.character(pop_full$Human_Entrez_ID)
pop_map <- suppressWarnings(select(org.Hs.eg.db, keys = pop_full$Human_Entrez_ID, keytype = "ENTREZID", columns = "SYMBOL"))
pop_sym <- merge(pop_full, pop_map, by.x = "Human_Entrez_ID", by.y = "ENTREZID", all.x = TRUE)

## Um gene do Chen pode ter mais de uma probeset (ex. PI3 aparece 2x) - dedup
chen_genes <- unique(chen[, c("GeneSymbol","Direction")])
## Se um symbol tiver linhas conflitantes de direcao (nao deveria acontecer
## aqui, mas por seguranca) mantem a primeira
chen_genes <- chen_genes[!duplicated(chen_genes$GeneSymbol), ]

result <- merge(chen_genes, pop_sym[, c("SYMBOL","baseMean","log2FoldChange","padj")],
                 by.x = "GeneSymbol", by.y = "SYMBOL", all.x = TRUE)
result$Testado_no_POP <- !is.na(result$log2FoldChange)
result$Significativo_no_POP_FDR05 <- !is.na(result$padj) & result$padj < 0.05 & abs(result$log2FoldChange) > 0.5
result$Direcao_POP <- ifelse(is.na(result$log2FoldChange), NA, ifelse(result$log2FoldChange > 0, "up", "down"))
result$Concordante <- result$Direction == result$Direcao_POP
result <- result[order(result$padj), ]

write.csv(result, "results/Chen2006_79genes_x_POP_completo.csv", row.names = FALSE)

cat("\nGenes testados no POP:", sum(result$Testado_no_POP), "de", nrow(result), "\n")
cat("Genes SIGNIFICATIVOS no POP (padj<0.05, |log2FC|>0.5):", sum(result$Significativo_no_POP_FDR05, na.rm=TRUE), "\n")

sig <- subset(result, Significativo_no_POP_FDR05)
cat("\n=== Genes significativos no POP entre os 55 candidatos de Chen 2006 ===\n")
print(sig[, c("GeneSymbol","Direction","log2FoldChange","padj","Direcao_POP","Concordante")])

testados <- subset(result, Testado_no_POP & !is.na(Direcao_POP))
n_concordante <- sum(testados$Concordante, na.rm = TRUE)
n_total <- nrow(testados)
cat("\nConcordancia de direcao (sinal bruto) em TODOS os", n_total, "genes testaveis:",
    n_concordante, "de", n_total, "\n")
bt <- binom.test(n_concordante, n_total, p = 0.5, alternative = "two.sided")
cat("Teste binomial de sinal: p =", format(bt$p.value, digits=4), "\n")
