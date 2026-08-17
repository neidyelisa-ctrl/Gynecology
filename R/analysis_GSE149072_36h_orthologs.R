# =============================================================================
# GSE149072 - Gene expression profiling of tissue and hMSC xenografts in a
# rat postpartum urinary injury model (Sadeghi et al. 2020, Tissue Eng Part A,
# DOI 10.1089/ten.tea.2020.0033)
#
# Objetivo deste script:
#   1) Baixar o GSE149072 do GEO.
#   2) Isolar as amostras de tecido de RATO (Rattus norvegicus) - o xenograft
#      de hMSC (Homo sapiens) fica em outra plataforma/subset do mesmo Series
#      e nao entra nesta comparacao.
#   3) Filtrar as amostras do tempo de 36h, grupo "lesao sem tratamento"
#      (vaginal distension / injury, SEM injecao de hMSC).
#   4) Comparar esse grupo contra o grupo controle/normal para achar os genes
#      diferencialmente expressos (DEGs).
#   5) Converter os DEGs de rato para ortologos humanos (biomaRt + orthogene).
#   6) Exportar tudo para Excel.
#
# IMPORTANTE - LEIA ANTES DE RODAR:
#   Este script NAO foi executado neste momento porque o ambiente onde ele
#   foi gerado bloqueia acesso a ncbi.nlm.nih.gov. Rode-o no seu RStudio,
#   com internet liberada.
#
#   Os nomes exatos das colunas/valores de metadados (grupo, tempo,
#   tratamento) do GSE149072 variam conforme como os autores rotularam as
#   amostras no GEO. A Secao 2 IMPRIME esses metadados brutos (pData) antes
#   de qualquer filtro - confira a saida no console e ajuste os padroes de
#   regex da Secao 3 (marcados com "AJUSTAR AQUI") para bater exatamente com
#   os rotulos reais, caso os padroes automaticos abaixo nao encontrem as
#   amostras certas ou encontrem amostras demais/de menos.
# =============================================================================

## -----------------------------------------------------------------------
## 1. Pacotes
## -----------------------------------------------------------------------

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

required_bioc <- c("GEOquery", "limma", "Biobase", "biomaRt")
missing_bioc <- required_bioc[!sapply(required_bioc, requireNamespace, quietly = TRUE)]
if (length(missing_bioc) > 0) BiocManager::install(missing_bioc, update = FALSE, ask = FALSE)

if (!requireNamespace("orthogene", quietly = TRUE)) {
  # orthogene facilita a conversao rato -> humano; se a instalacao falhar,
  # o script cai automaticamente no metodo biomaRt (Secao 6B).
  try(BiocManager::install("orthogene", update = FALSE, ask = FALSE), silent = TRUE)
}

required_cran <- c("dplyr", "tibble", "openxlsx", "stringr")
missing_cran <- required_cran[!sapply(required_cran, requireNamespace, quietly = TRUE)]
if (length(missing_cran) > 0) install.packages(missing_cran)

library(GEOquery)
library(limma)
library(Biobase)
library(dplyr)
library(tibble)
library(stringr)
library(openxlsx)

dir.create("data/geo_cache", showWarnings = FALSE, recursive = TRUE)
dir.create("results", showWarnings = FALSE, recursive = TRUE)

## -----------------------------------------------------------------------
## 2. Download do GSE149072 e identificacao do subset de rato
## -----------------------------------------------------------------------

GSE_ID <- "GSE149072"

gse_list <- getGEO(GSE_ID, GSEMatrix = TRUE, AnnotGPL = TRUE,
                    destdir = "data/geo_cache")

cat("Numero de subsets/plataformas encontrados:", length(gse_list), "\n")
for (i in seq_along(gse_list)) {
  cat("\n--- Subset", i, "| Plataforma:", annotation(gse_list[[i]]), "---\n")
  print(table(pData(gse_list[[i]])$organism_ch1))
}

## Seleciona o(s) subset(s) cujas amostras sao de Rattus norvegicus
is_rat <- sapply(gse_list, function(x) {
  any(grepl("Rattus norvegicus", pData(x)$organism_ch1, ignore.case = TRUE))
})
stopifnot("Nenhum subset de Rattus norvegicus encontrado - confira gse_list manualmente" = any(is_rat))

eset_rat <- gse_list[[which(is_rat)[1]]]
cat("\nUsando subset de rato com", ncol(eset_rat), "amostras e",
    nrow(eset_rat), "genes/sondas.\n")

## -----------------------------------------------------------------------
## 3. Inspecionar metadados reais das amostras (OBRIGATORIO antes de filtrar)
## -----------------------------------------------------------------------

pd <- pData(eset_rat)

## Imprime todas as colunas de caracteristicas/titulo/fonte para voce
## conferir os rotulos EXATOS usados neste dataset
meta_cols <- grep("characteristics|title|source_name|description",
                   colnames(pd), ignore.case = TRUE, value = TRUE)
print(pd[, meta_cols])

## AJUSTAR AQUI se necessario: nomes das colunas de tempo e de grupo/tratamento
time_col  <- grep("time", colnames(pd), ignore.case = TRUE, value = TRUE)[1]
group_col <- grep("treat|agent|group|injury|condition",
                   colnames(pd), ignore.case = TRUE, value = TRUE)[1]

cat("\nColuna de tempo detectada:", time_col, "\n")
cat("Coluna de grupo/tratamento detectada:", group_col, "\n")
cat("\nValores unicos de tempo:\n"); print(unique(pd[[time_col]]))
cat("\nValores unicos de grupo:\n"); print(unique(pd[[group_col]]))

## -----------------------------------------------------------------------
## 4. Definir os dois grupos de comparacao
##    Grupo A = lesao (vaginal distension), SEM tratamento com hMSC, 36h
##    Grupo B = controle / tecido normal (nao lesionado)
## -----------------------------------------------------------------------
## AJUSTAR AQUI conforme os valores reais impressos acima.

is_36h <- str_detect(as.character(pd[[time_col]]), "36")

is_injury_untreated <- str_detect(as.character(pd[[group_col]]),
                                   regex("injur|distension|VD|lesao|trauma", ignore_case = TRUE)) &
                        !str_detect(as.character(pd[[group_col]]),
                                    regex("hMSC|MSC|stem cell|treated|therap", ignore_case = TRUE))

is_control <- str_detect(as.character(pd[[group_col]]),
                          regex("control|normal|sham|uninjured|naive", ignore_case = TRUE))

grupo_lesao_36h <- rownames(pd)[is_36h & is_injury_untreated]
grupo_controle  <- rownames(pd)[is_control]

cat("\nAmostras no grupo LESAO 36h SEM tratamento:", length(grupo_lesao_36h), "\n")
print(grupo_lesao_36h)
cat("\nAmostras no grupo CONTROLE:", length(grupo_controle), "\n")
print(grupo_controle)

## Se o dataset nao tiver um grupo controle "puro" separado por tempo,
## uma alternativa e comparar contra o grupo lesao+hMSC no mesmo tempo (36h)
## para isolar o efeito do tratamento, em vez do efeito da lesao:
# grupo_alt_tratado_36h <- rownames(pd)[is_36h & str_detect(as.character(pd[[group_col]]),
#                            regex("hMSC|MSC|treated", ignore_case = TRUE))]

stopifnot("Grupo de lesao 36h vazio - ajuste os regex da Secao 3/4" = length(grupo_lesao_36h) > 0)
stopifnot("Grupo controle vazio - ajuste os regex da Secao 3/4" = length(grupo_controle) > 0)

amostras_usadas <- c(grupo_controle, grupo_lesao_36h)
grupo <- factor(c(rep("Controle", length(grupo_controle)),
                   rep("Lesao_36h_SemTratamento", length(grupo_lesao_36h))),
                levels = c("Controle", "Lesao_36h_SemTratamento"))

eset_comp <- eset_rat[, amostras_usadas]

## -----------------------------------------------------------------------
## 5. Expressao diferencial (limma)
## -----------------------------------------------------------------------

exprs_mat <- exprs(eset_comp)
if (max(exprs_mat, na.rm = TRUE) > 100) {
  exprs_mat <- log2(exprs_mat + 1)
}

design <- model.matrix(~ grupo)
fit <- lmFit(exprs_mat, design)
fit <- eBayes(fit)

deg_table <- topTable(fit, coef = 2, number = Inf, adjust.method = "BH") %>%
  tibble::rownames_to_column("Probe_ID") %>%
  arrange(adj.P.Val)

## Anexa o simbolo do gene de rato, se disponivel na anotacao da plataforma
fdata <- fData(eset_rat)
symbol_col <- grep("gene.symbol|symbol", colnames(fdata), ignore.case = TRUE, value = TRUE)[1]
if (!is.na(symbol_col)) {
  deg_table$Rat_Gene_Symbol <- fdata[deg_table$Probe_ID, symbol_col]
} else {
  deg_table$Rat_Gene_Symbol <- deg_table$Probe_ID
}

deg_sig <- deg_table %>% filter(adj.P.Val < 0.05, abs(logFC) > 1)
cat("\nGenes diferencialmente expressos (Lesao_36h_SemTratamento vs Controle,",
    "padj<0.05, |logFC|>1):", nrow(deg_sig), "\n")

write.csv(deg_table, file.path("results", "GSE149072_36h_lesao_vs_controle_completo.csv"),
          row.names = FALSE)
write.csv(deg_sig, file.path("results", "GSE149072_36h_lesao_vs_controle_significativos.csv"),
          row.names = FALSE)

## -----------------------------------------------------------------------
## 6A. Ortologos humanos - metodo orthogene (recomendado, mais simples)
## -----------------------------------------------------------------------

rat_genes <- unique(na.omit(deg_sig$Rat_Gene_Symbol))
orth_table <- NULL

if (requireNamespace("orthogene", quietly = TRUE) && length(rat_genes) > 0) {
  orth_result <- tryCatch({
    orthogene::convert_orthologs(
      gene_df = rat_genes,
      input_species = "rat",
      output_species = "human",
      non121_strategy = "drop_both_species",
      method = "homologene"
    )
  }, error = function(e) { message("orthogene falhou: ", e$message); NULL })

  if (!is.null(orth_result)) {
    orth_table <- orth_result %>%
      tibble::rownames_to_column("Human_Ortholog_Symbol") %>%
      rename(Rat_Gene_Symbol = input_gene)
  }
}

## -----------------------------------------------------------------------
## 6B. Ortologos humanos - metodo biomaRt (fallback caso orthogene falhe)
## -----------------------------------------------------------------------

if (is.null(orth_table) && length(rat_genes) > 0) {
  rat_mart   <- biomaRt::useMart("ensembl", dataset = "rnorvegicus_gene_ensembl")
  human_mart <- biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")

  orth_biomart <- biomaRt::getLDS(
    attributes  = "external_gene_name",
    filters     = "external_gene_name",
    values      = rat_genes,
    mart        = rat_mart,
    attributesL = "hgnc_symbol",
    martL       = human_mart,
    uniqueRows  = TRUE
  )
  colnames(orth_biomart) <- c("Rat_Gene_Symbol", "Human_Ortholog_Symbol")
  orth_table <- orth_biomart %>% filter(Human_Ortholog_Symbol != "")
}

## -----------------------------------------------------------------------
## 7. Juntar DEGs + ortologos humanos e exportar
## -----------------------------------------------------------------------

if (!is.null(orth_table)) {
  deg_com_ortologos <- deg_sig %>%
    inner_join(orth_table, by = "Rat_Gene_Symbol") %>%
    arrange(adj.P.Val)

  cat("\nGenes com ortologo humano identificado:", nrow(deg_com_ortologos), "de",
      nrow(deg_sig), "DEGs significativos.\n")

  write.csv(deg_com_ortologos,
            file.path("results", "GSE149072_36h_DEGs_ortologos_humanos.csv"),
            row.names = FALSE)

  wb <- createWorkbook()
  addWorksheet(wb, "DEG_rato_36h_vs_controle")
  writeData(wb, "DEG_rato_36h_vs_controle", deg_sig)
  addWorksheet(wb, "DEG_com_ortologos_humanos")
  writeData(wb, "DEG_com_ortologos_humanos", deg_com_ortologos)
  saveWorkbook(wb, file.path("results", "GSE149072_36h_resultado_final.xlsx"),
               overwrite = TRUE)

  cat("\nResultado final salvo em results/GSE149072_36h_resultado_final.xlsx\n")
} else {
  warning("Nao foi possivel mapear ortologos humanos (orthogene e biomaRt falharam ",
          "ou nao ha genes significativos). Confira conexao com Ensembl/Homologene.")
}

cat("\nAnalise concluida.\n")
