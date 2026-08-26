# =============================================================================
# Mendelian randomization (MR) dos 6 genes validados (CWH43, INPP4B, CALML5,
# KRT10, SERPINB2, DMKN) contra POP - e, se achar GWAS de SUI, contra SUI.
#
# IMPORTANTE - POR QUE ESTE SCRIPT NAO RODA NESTE SANDBOX:
# MR de verdade exige baixar dados de GWAS/eQTL em tempo real (via
# TwoSampleMR/ieugwasr, que consultam a API do IEU OpenGWAS
# https://gwas.mrcieu.ac.uk). Testei e confirmei bloqueado (erro 403, mesma
# politica de rede de todo o projeto - o mesmo motivo pelo qual STRING/
# Ensembl/KEGG ao vivo tambem nao funcionam aqui). Este script foi escrito
# para voce RODAR NO SEU PROPRIO AMBIENTE (RStudio com internet normal), nao
# aqui. Ele documenta exatamente os IDs de GWAS reais que encontrei via busca
# na literatura - nao inventei nenhum accession.
#
# DESENHO DO MR (two-sample MR classico, expression -> disease):
#   - Exposicao: eQTLs (variantes que afetam a expressao) de cada um dos 6
#     genes - fonte recomendada: eQTLGen Consortium (sangue, sample size
#     grande, ~31.000 pessoas) via OpenGWAS, ou GTEx v8 se quiser eQTL de
#     tecido especifico (pele para KRT10/DMKN/CALML5, por exemplo).
#   - Desfecho (POP): FinnGen "finn-b-N14_FEMGENPROL" (prolapso genital
#     feminino) - 9.092 casos / 68.969 controles. Accession real, usado no
#     estudo "Unraveling the Causal Linkages of RBP7 and SCGB3A1 on Pelvic
#     Organ Prolapse" (PMC12765987) com a MESMA metodologia que voce quer
#     aplicar aos seus 6 genes.
#   - Desfecho (SUI/incontinencia urinaria): ainda NAO existe um GWAS
#     dedicado so a SUI com accession publico confirmado que eu tenha
#     achado. O mais promissor e um meta-analise de 2026 (medRxiv, ainda
#     preprint) com 54 loci para incontinencia urinaria e subtipos,
#     >1 milhao de individuos (HUNT+UK Biobank+FinnGen+MGI) - confira se ja
#     tem GWAS Catalog/OpenGWAS ID publicado antes de usar.
# =============================================================================

if (!requireNamespace("TwoSampleMR", quietly = TRUE)) {
  remotes::install_github("MRCIEU/TwoSampleMR")  # requer internet
}
library(TwoSampleMR)

genes <- c("CWH43", "INPP4B", "CALML5", "KRT10", "SERPINB2", "DMKN")

## 1) Exposicao: eQTLs de cada gene (eQTLGen via OpenGWAS, ou substitua pelo
##    ID do GTEx do tecido que fizer mais sentido - pele para KRT10/DMKN/
##    CALML5, ou "whole blood" do eQTLGen para todos)
##    Ex.: ao_agree() e authentication podem ser necessarios - ver
##    documentacao do pacote ieugwasr/TwoSampleMR.
exposure_list <- lapply(genes, function(g) {
  extract_instruments(outcomes = paste0("eqtl-a-", g))  # ID ilustrativo -
  # confirme o ID exato de cada gene no catalogo eQTLGen dentro do
  # OpenGWAS (plataforma "eQTLGen") antes de rodar de verdade.
})
names(exposure_list) <- genes

## 2) Desfecho: FinnGen POP (accession real, confirmado por busca)
outcome_pop <- extract_outcome_data(
  snps = unlist(lapply(exposure_list, function(x) x$SNP)),
  outcomes = "finn-b-N14_FEMGENPROL"
)

## 3) Harmonizar e rodar MR gene a gene
resultados <- list()
for (g in genes) {
  exp_g <- exposure_list[[g]]
  if (nrow(exp_g) == 0) { cat(g, "- sem instrumentos (SNPs) suficientes, pular\n"); next }
  dat <- harmonise_data(exposure_dat = exp_g, outcome_dat = outcome_pop)
  res <- mr(dat, method_list = c("mr_ivw", "mr_egger_regression", "mr_weighted_median"))
  res$gene <- g
  resultados[[g]] <- res
}

resultado_final <- do.call(rbind, resultados)
print(resultado_final)
write.csv(resultado_final, "MR_6genes_POP_resultado.csv", row.names = FALSE)

## Sensibilidade (rodar para cada gene com instrumentos suficientes,
## sobretudo se >=3 SNPs):
## mr_heterogeneity(dat); mr_pleiotropy_test(dat); mr_scatter_plot(res, dat)
