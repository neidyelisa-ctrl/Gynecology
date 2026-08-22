# =============================================================================
# SUI(36h) x POP - reproduzindo a receita exata da usuaria, mas com ortologos
# via Ensembl (nao HomoloGene), para tentar bater exatamente com os 6 genes
# dela: KRT10, SERPINB2, CALML5, CWH43, INPP4B, DMKN.
#
# POR QUE HomoloGene nao batia 100%: KRT10 nao tem grupo de ortologo no
# HomoloGene em nenhuma especie de rato; DMKN tem grupo no HomoloGene mas
# sem rato (so humano/chimpanze/macaco/camundongo). Ambos SO aparecem via
# Ensembl.
#
# COMO CONSEGUI ORTOLOGOS DO ENSEMBL OFFLINE: o BioMart/Ensembl ao vivo esta
# bloqueado neste sandbox (mesmo teste de sempre, reconfirmado - CONNECT
# tunnel 403 para ensembl.org, rest.ensembl.org, biomart). Em vez disso, usei
# o pacote R "babelgene" (CRAN/GitHub, MIT license) que empacota, de forma
# OFFLINE (nao acessa nenhum servidor em tempo de execucao), uma tabela de
# ortologos JA compilada a partir de 9 bases, entre elas Ensembl, HomoloGene,
# NCBI, OMA, OrthoDB, Panther, Treefam, EggNOG e Inparanoid (fonte: HCOP -
# HGNC Comparison of Orthology Predictions). Baixei o arquivo de dados
# interno do pacote (R/sysdata.rda, MIT) direto do repositorio GitHub
# (raw.githubusercontent.com, que este sandbox consegue acessar) e salvei em
# data/babelgene_orthologs.rda. Conferi manualmente: KRT10->Krt10 (Ensembl
# entre as fontes) e DMKN->Dmkn (Ensembl entre as fontes) - os dois genes que
# faltavam no HomoloGene aparecem aqui.
#
# Resto da receita = exatamente a da usuaria (ANALYSIS_SUI.docx):
#   - SEM filtro de baixa contagem
#   - SEM covariavel de idade no POP (design = ~ condition)
#   - SEM lfcShrink (log2FoldChange bruto de results())
#   - FDR (padj) < 0.05, |log2FoldChange| > 0.5
# =============================================================================

library(DESeq2)

FDR_CUT <- 0.05
LFC_CUT <- 0.5

load("data/babelgene_orthologs.rda")  # -> orthologs_df (human_symbol, taxon_id, symbol, entrez, support, ...)
rat_orth_all <- orthologs_df[orthologs_df$taxon_id == 10116, c("symbol", "human_symbol", "human_entrez")]
colnames(rat_orth_all) <- c("Rat_Gene_Symbol", "Human_Ortholog_Symbol", "Human_Entrez_ID")
rat_orth_all <- unique(rat_orth_all)

## Filtro de ortologia 1-para-1 (a usuaria relata no thesis_proposal_draft.pdf: "197 (DEG)
## rat -> human orthologs. - 192 human genes with one-to-one orthology" - ela usou SO os
## pares 1-para-1 do BioMart). Reproduzido aqui: mantem so pares onde o gene de rato mapeia
## para exatamente 1 gene humano E o gene humano mapeia para exatamente 1 gene de rato.
rat_n <- table(rat_orth_all$Rat_Gene_Symbol)
hum_n <- table(rat_orth_all$Human_Ortholog_Symbol)
rat_orth <- rat_orth_all[rat_n[rat_orth_all$Rat_Gene_Symbol] == 1 &
                            hum_n[rat_orth_all$Human_Ortholog_Symbol] == 1, ]
cat("Tabela de ortologos rato->humano (Ensembl+8 outras bases, via babelgene):",
    nrow(rat_orth_all), "pares totais |", nrow(rat_orth), "pares 1-para-1\n")

## -----------------------------------------------------------------------
## SUI 36h - SEM filtro, SEM shrink
## -----------------------------------------------------------------------
sui_counts <- read.csv("data/GSE149072_rawCounts.csv", row.names = 1, check.names = FALSE)
sui_counts <- round(as.matrix(sui_counts)); mode(sui_counts) <- "integer"

untreated <- grep("Rat_Urethra_Untreated_36hr", colnames(sui_counts), value = TRUE)
treated   <- grep("Rat_Urethra_Treated_36hr", colnames(sui_counts), value = TRUE)
sel <- c(untreated, treated)
sui_counts_sub <- sui_counts[, sel]

group <- factor(c(rep("Untreated", length(untreated)), rep("Treated", length(treated))),
                 levels = c("Untreated", "Treated"))
coldata_sui <- data.frame(row.names = sel, group = group)

dds_sui <- DESeqDataSetFromMatrix(countData = sui_counts_sub, colData = coldata_sui, design = ~ group)
dds_sui <- DESeq(dds_sui)
res_sui <- results(dds_sui, contrast = c("group", "Treated", "Untreated"), alpha = FDR_CUT)  # SEM shrink

sui_df <- as.data.frame(res_sui)
sui_df$Rat_Gene_Symbol <- rownames(sui_df)
sui_df <- sui_df[order(sui_df$padj), ]
rownames(sui_df) <- NULL

sui_sig <- subset(sui_df, !is.na(padj) & padj < FDR_CUT & abs(log2FoldChange) > LFC_CUT)
cat("SUI 36h (sem filtro, sem shrink) - DEGs:", nrow(sui_sig), "\n")

sui_sig_orth <- merge(sui_sig, rat_orth, by = "Rat_Gene_Symbol", all.x = TRUE)
sui_sig_orth <- sui_sig_orth[order(sui_sig_orth$padj), ]
write.csv(sui_sig_orth, "results/SUI_36hr_ensembl_DEG_com_ortologos.csv", row.names = FALSE)
cat("SUI 36h DEGs com ortologo humano (via Ensembl+outras bases):",
    sum(!is.na(sui_sig_orth$Human_Entrez_ID)), "de", nrow(sui_sig), "\n")

for (g in c("Cwh43","Inpp4b","Calml5","Krt10","Serpinb2","Dmkn")) {
  cat(" ", g, "no DEG list do SUI?", g %in% sui_sig$Rat_Gene_Symbol, "\n")
}

## -----------------------------------------------------------------------
## POP - SEM filtro, SEM covariavel de idade, SEM shrink
## -----------------------------------------------------------------------
pop_counts <- read.delim("data/GSE208261_raw_counts_GRCh38.p13_NCBI.tsv", row.names = 1, check.names = FALSE)
pop_counts <- round(as.matrix(pop_counts)); mode(pop_counts) <- "integer"

sample_labels <- c(
  GSM6339911="Control_D1", GSM6339912="Control_D2", GSM6339913="Control_D3",
  GSM6339914="Control_D4", GSM6339915="Control_D5", GSM6339916="Control_D6",
  GSM6339917="Control_Y1", GSM6339918="Control_Y2", GSM6339919="Control_Y3",
  GSM6339920="Control_Y4", GSM6339921="Control_Y5", GSM6339922="Control_Y6",
  GSM6339923="POP_D1", GSM6339924="POP_D2", GSM6339925="POP_D3",
  GSM6339926="POP_D4", GSM6339927="POP_D5", GSM6339928="POP_D6",
  GSM6339929="POP_Y1", GSM6339930="POP_Y2", GSM6339931="POP_Y3",
  GSM6339932="POP_Y4", GSM6339933="POP_Y5", GSM6339934="POP_Y6"
)
label <- sample_labels[colnames(pop_counts)]
condition <- factor(ifelse(grepl("^POP", label), "POP", "Control"), levels = c("Control", "POP"))
coldata_pop <- data.frame(row.names = colnames(pop_counts), condition = condition)

dds_pop <- DESeqDataSetFromMatrix(countData = pop_counts, colData = coldata_pop, design = ~ condition)
dds_pop <- DESeq(dds_pop)
res_pop <- results(dds_pop, contrast = c("condition", "POP", "Control"), alpha = FDR_CUT)  # SEM shrink

pop_df <- as.data.frame(res_pop)
pop_df$Human_Entrez_ID <- rownames(pop_df)
pop_df <- pop_df[order(pop_df$padj), ]
rownames(pop_df) <- NULL

pop_sig <- subset(pop_df, !is.na(padj) & padj < FDR_CUT & abs(log2FoldChange) > LFC_CUT)
cat("\nPOP (sem filtro, sem covariavel de idade, sem shrink) - DEGs:", nrow(pop_sig), "\n")
write.csv(pop_sig, "results/POP_GSE208261_recipe_usuaria_DEG_sig.csv", row.names = FALSE)

## -----------------------------------------------------------------------
## Cruzamento
## -----------------------------------------------------------------------
pop_sig$Human_Entrez_ID <- as.character(pop_sig$Human_Entrez_ID)
sui_valid <- subset(sui_sig_orth, !is.na(Human_Entrez_ID))
sui_valid$Human_Entrez_ID <- as.character(sui_valid$Human_Entrez_ID)

common_ids <- intersect(sui_valid$Human_Entrez_ID, pop_sig$Human_Entrez_ID)
common_genes <- sui_valid[sui_valid$Human_Entrez_ID %in% common_ids,
                           c("Rat_Gene_Symbol","Human_Ortholog_Symbol","Human_Entrez_ID","log2FoldChange","padj")]
colnames(common_genes)[4:5] <- c("SUI_36h_log2FC","SUI_36h_padj")
pop_match <- pop_sig[pop_sig$Human_Entrez_ID %in% common_ids, c("Human_Entrez_ID","log2FoldChange","padj")]
colnames(pop_match)[2:3] <- c("POP_log2FC","POP_padj")
common_genes <- merge(common_genes, pop_match, by = "Human_Entrez_ID")
common_genes$SUI_dir <- ifelse(common_genes$SUI_36h_log2FC > 0, "up", "down")
common_genes$POP_dir <- ifelse(common_genes$POP_log2FC > 0, "up", "down")
common_genes$same_direction <- common_genes$SUI_dir == common_genes$POP_dir
common_genes <- common_genes[order(common_genes$POP_padj), ]

cat("\n=== Genes em comum (receita da usuaria + ortologos via Ensembl/babelgene) ===\n")
print(common_genes[, c("Rat_Gene_Symbol","Human_Ortholog_Symbol","SUI_36h_log2FC","SUI_36h_padj",
                        "POP_log2FC","POP_padj","same_direction")])
write.csv(common_genes, "results/SUI_36hr_ensembl_x_POP_common_genes.csv", row.names = FALSE)

cat("\n=== Comparacao com a lista dela ===\n")
dela <- c("CWH43","INPP4B","CALML5","KRT10","SERPINB2","DMKN")
minha <- toupper(common_genes$Human_Ortholog_Symbol)
cat("Bateram:", paste(intersect(dela, minha), collapse=", "), "(", length(intersect(dela,minha)), "de", length(dela), ")\n")
cat("So na dela:", paste(setdiff(dela, minha), collapse=", "), "\n")
cat("So na minha:", paste(setdiff(minha, dela), collapse=", "), "\n")
