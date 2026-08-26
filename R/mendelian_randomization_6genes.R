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
library(ieugwasr)

## -----------------------------------------------------------------------
## AUTENTICACAO (obrigatoria desde 1/maio/2024) - PASSO NOVO que faltava
## -----------------------------------------------------------------------
## 1) No navegador, va em https://api.opengwas.io/profile/ e faca login
##    (Microsoft, GitHub ou e-mail).
## 2) Nessa mesma pagina de perfil, clique em "Generate new token" (ou
##    "Manage tokens") e copie o token gerado - e uma string longa, tipo
##    senha. NAO compartilhe esse token com ninguem.
## 3) Salve o token no seu arquivo .Renviron (NAO no script, para nao
##    vazar sem querer se voce compartilhar o codigo). No RStudio, rode:
##      usethis::edit_r_environ()
##    (instale o pacote "usethis" se nao tiver: install.packages("usethis"))
##    Isso abre o arquivo .Renviron. Adicione esta linha (troque pelo seu
##    token de verdade) e salve:
##      OPENGWAS_JWT=coloque_aqui_o_token_que_voce_copiou
##    Feche e REINICIE a sessao do R (Session > Restart R no RStudio) para
##    o .Renviron ser recarregado.
## 4) Confirme que funcionou - deve imprimir uma string longa (o token) e
##    NAO ficar vazio/NULL:
print(ieugwasr::get_opengwas_jwt())
## Se aparecer vazio, o .Renviron nao foi lido corretamente - confira se
## salvou no arquivo certo (rode usethis::edit_r_environ() de novo para
## abrir o mesmo arquivo) e se reiniciou a sessao do R depois de salvar.
##
## 5) Teste rapido de conexao (deve retornar informacoes da sua conta, sem
##    erro 401):
print(ieugwasr::user())
## -----------------------------------------------------------------------

## CORRECAO: o dataset "eqtl-a-*" do OpenGWAS usa ID Ensembl do GENE
## (ENSG...), NAO o simbolo (era esse o erro do "argument is of length
## zero" - a query com "eqtl-a-CWH43" nao existe, a API devolveu vazio).
## IDs confirmados localmente via org.Hs.eg.db (Bioconductor, nao inventados):
gene_ensembl <- c(
  CWH43    = "ENSG00000109182",
  INPP4B   = "ENSG00000109452",
  CALML5   = "ENSG00000178372",
  KRT10    = "ENSG00000186395",
  SERPINB2 = "ENSG00000197632",
  DMKN     = "ENSG00000161249"
)
genes <- names(gene_ensembl)

## 1) Exposicao: eQTLs de cada gene (eQTLGen via OpenGWAS - catalogo "eqtl-a-",
##    baseado em sangue). Se algum gene nao tiver instrumento nesse catalogo
##    (comum para genes de expressao restrita a pele/epiderme, como CALML5/
##    DMKN/CWH43/KRT10 - o catalogo eqtl-a e majoritariamente de sangue),
##    considere GTEx v8 pele (nao disponivel via OpenGWAS "eqtl-a", precisa
##    de outra fonte, ex. GTEx portal diretamente) como alternativa.
## tryCatch por gene para um ID sem dado no catalogo nao travar o loop todo.
exposure_list <- lapply(genes, function(g) {
  ensg <- gene_ensembl[[g]]
  out <- tryCatch(
    extract_instruments(outcomes = paste0("eqtl-a-", ensg)),
    error = function(e) {
      cat(g, "(", ensg, ") - sem instrumentos nesse catalogo eQTL:", conditionMessage(e), "\n")
      NULL
    }
  )
  if (is.null(out) || nrow(out) == 0) {
    cat(g, "(", ensg, ") - NENHUM instrumento encontrado no catalogo eqtl-a (sangue).",
        "Provavel gene de expressao tecido-especifica (pele) nao coberto por esse painel.\n")
    return(data.frame())
  }
  out
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
