# Datasets de Bioinformática — Prolapso de Órgãos Pélvicos (POP) e Incontinência Urinária de Esforço (SUI)

Levantamento de datasets públicos de expressão gênica (GEO/ArrayExpress) relacionados a
Prolapso de Órgãos Pélvicos (POP) e Incontinência Urinária de Esforço (SUI), com script de
análise em R.

## Conteúdo

- `data/POP_SUI_GEO_datasets.xlsx` — lista consolidada dos datasets encontrados (accession,
  organismo, tecido, plataforma, desenho experimental, referência e link), com abas de resumo
  e notas metodológicas.
- `R/analysis_geo_pop_sui.R` — script para RStudio que baixa um dataset do GEO (`GEOquery`),
  roda controle de qualidade (PCA), expressão diferencial (`limma`), volcano plot, heatmap e
  enriquecimento funcional GO/KEGG (`clusterProfiler`).

## Como usar

1. Abra `R/analysis_geo_pop_sui.R` no RStudio.
2. Rode a Seção 1 para instalar os pacotes necessários (Bioconductor + CRAN).
3. Ajuste `GSE_ID` na Seção 3 para o accession desejado (ex.: `GSE53868`, `GSE12852`,
   `GSE28660`) — a lista completa está em `data/POP_SUI_GEO_datasets.xlsx`.
4. Execute o script por blocos. Os resultados (CSV, gráficos e Excel final) são salvos em
   `results/`.

## Principal achado da busca

Foram encontrados vários datasets públicos no GEO dedicados a **POP** (parede vaginal
anterior e ligamentos uterossacrais/redondos, incluindo um dataset de single-cell RNA-seq).
Para **SUI** isoladamente, não foi localizado nenhum dataset de expressão gênica com
accession pública confirmada — a literatura de referência usa microarray/proteômica sem
depósito público identificado nesta busca. Detalhes e recomendações para expandir a busca
estão na aba "Notas_metodologicas" da planilha.
