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
- `R/analysis_GSE149072_36h_orthologs.R` — pipeline **real e já executado** sobre o
  **GSE149072** ("Gene expression profiling of tissue and hMSC xenografts in a rat
  postpartum urinary injury model", Sadeghi et al. 2020, DOI 10.1089/ten.tea.2020.0033).
  Roda DESeq2 sobre `data/GSE149072_rawCounts.csv` comparando uretra de rata em 36h
  pós-lesão **sem tratamento** (Untreated) vs. **com hMSC** (Treated), e converte os DEGs
  de rato para ortólogos humanos via HomoloGene (pacote `homologene`).
- `data/GSE149072_rawCounts.csv` — contagens brutas de RNA-seq do GSE149072 (36 amostras:
  uretra de rato tratada/não tratada em 4 tempos × 3 réplicas, mais células hMSC humanas
  isoladas), fornecidas pelo usuário a partir da página do GEO.
- `data/GSE208261_raw_counts_GRCh38.p13_NCBI.tsv` — contagens brutas de RNA-seq do GSE208261
  (24 amostras humanas: POP vs. Controle, 12 idosas + 12 jovens), fornecidas pelo usuário.
- `results/GSE149072_36h_DEGs_ortologos_humanos.xlsx` — DEGs de SUI em 36h (padj<0.05,
  |log2FC|>1) com ortólogos humanos.
- `R/analysis_SUI_POP_crosscomparison.R` — pipeline **real e já executado** que roda os DEGs
  de SUI (GSE149072, 36h e 72h) e de POP (GSE208261), converte os DEGs de rato para
  ortólogos humanos via HomoloGene, e cruza os dois conjuntos para achar genes em comum.
- `results/SUI_x_POP_36h_genes_comuns.xlsx` — resultado real do cruzamento: os genes em
  comum entre SUI e POP em 36h (única janela de tempo com correspondências), já classificados
  como concordantes ou discordantes na direção da desregulação.
- `R/analysis_neuro_candidate_genes.R` — testa uma lista de 25 genes candidatos da via
  neurotrófica/neuropeptídica (curada de Shen et al. 2024 e Masyhuroh et al. 2024) contra os
  DEGs reais de SUI e POP já calculados — a base da proposta de tema de pesquisa abaixo.
- `results/Proposta_tema_neuro_POP_SUI.xlsx` — proposta de tema de pesquisa: sinalização
  neurotrófica/neuropeptídica como eixo convergente entre POP e SUI, com fundamentação na
  literatura, tabela de genes candidatos testados contra os dados reais, e roteiro
  metodológico para expandir o trabalho.

## Como usar

1. Abra o script desejado no RStudio.
2. Rode a Seção 1 para instalar os pacotes necessários (Bioconductor + CRAN).
3. Em `analysis_geo_pop_sui.R`, ajuste `GSE_ID` na Seção 3 para o accession desejado (ex.:
   `GSE53868`, `GSE12852`, `GSE28660`) — a lista completa está em
   `data/POP_SUI_GEO_datasets.xlsx`.
4. `analysis_GSE149072_36h_orthologs.R` e `analysis_SUI_POP_crosscomparison.R` já estão
   prontos para rodar como estão — leem os arquivos de `data/` e reproduzem os resultados
   salvos em `results/`.
5. Execute o script por blocos. Os resultados (CSV e Excel final) são salvos em `results/`.

## Principal achado da busca de datasets

Foram encontrados vários datasets públicos no GEO dedicados a **POP** (parede vaginal
anterior e ligamentos uterossacrais/redondos, incluindo um dataset de single-cell RNA-seq).
Para **SUI** isoladamente, não foi localizado nenhum dataset de expressão gênica com
accession pública confirmada por busca automática — mas o usuário identificou e forneceu o
**GSE149072** (modelo de SUI/lesão pós-parto em ratas, com xenograft de hMSC humano), que foi
efetivamente baixado e analisado (ver acima). Detalhes e recomendações para expandir a busca
estão na aba "Notas_metodologicas" da planilha `POP_SUI_GEO_datasets.xlsx`.

## Resultado da análise do GSE149072 (36h, lesão sem tratamento vs. tratada com hMSC)

- Comparação: `Rat_Urethra_Untreated_36hr` (n=3) vs. `Rat_Urethra_Treated_36hr` (n=3) — única
  comparação possível no arquivo de contagens para esse tempo, já que não há amostras de
  uretra normal/não lesionada no arquivo fornecido.
- Método: DESeq2, filtro de baixa expressão, shrinkage de log2FoldChange tipo "normal".
- 12.806 genes testados após filtro; **110 DEGs significativos** (padj<0.05, |log2FC|>1);
  **88 com ortólogo humano** identificado via NCBI HomoloGene.
- Aviso de qualidade: a amostra `Rat_Urethra_Treated_12hr_M2` tem profundidade de
  sequenciamento muito abaixo das demais — não afeta esta comparação (é do tempo de 12h),
  mas vale conferir antes de usar essa amostra em outra análise.

## Cruzamento SUI x POP (FDR<0.05, |log2FC|>0.5)

- **SUI (GSE149072)**: Untreated vs. Treated (hMSC) em uretra de rata — 36h: **110 DEGs**
  (88 com ortólogo humano); 72h: **58 DEGs** (47 com ortólogo humano).
- **POP (GSE208261)**: POP (n=12) vs. Controle (n=12), ajustando por grupo etário (idosa/jovem)
  — **642 DEGs**.
- Cruzamento dos ortólogos humanos do SUI com os DEGs de POP: **4 genes em comum em 36h,
  0 em 72h** — por isso 36h é a única janela de tempo viável para este cruzamento.
- Os 4 genes em comum: **INPP4B**, **ECM1**, **BEND3** (direção discordante entre as duas
  doenças) e **KREMEN1** (única com direção concordante — regulado para cima tanto na lesão
  não tratada quanto no POP). Ver a aba `Genes_em_comum_36h` da planilha para os números
  completos (log2FC e padj de cada análise) e a lógica de classificação de direção.

## Proposta de tema: convergência via sinalização neural (não gene-a-gene)

Cruzar o transcriptoma inteiro gene-a-gene deu uma amostra pequena demais (4-6 genes) para
testar enriquecimento estatístico de sobreposição. Alternativa mais robusta, ancorada em
literatura publicada (Shen et al. 2024 — revisão sistemática sobre nervo pélvico em POP/SUI;
Masyhuroh et al. 2024): testar se o **mesmo sistema biológico** (sinalização de
neurotrofinas/neuropeptídeos ligada à inervação pélvica) está desregulado nas duas condições,
mesmo que por genes diferentes dentro desse sistema.

Achado real (25 genes candidatos testados, ver `results/candidate_neuro_genes_table.csv`):
- **EDNRA** e **GFRA3** (receptor da família GDNF) — significativos no **SUI em 36h**
  (padj = 0.0003 e 0.013), não no POP.
- **NPY1R**, **NPY5R** e **NTF3** — significativos no **POP** (padj = 0.002–0.034), não
  testados/detectados no SUI.
- **ADCYAP1** (gene do PACAP) — tendência forte no POP (padj = 0.079), consistente com PACAP
  reduzido relatado na literatura em tecido humano de SUI/POP.
- Nenhum gene bateu significativo nos dois lados ao mesmo tempo, mas o padrão inteiro aponta
  para o mesmo sistema (neurotrofinas/neuropeptídeos/receptores GDNF) desregulado
  independentemente nos dois lados — uma conclusão honesta e ancorada em literatura, em vez
  de forçar uma sobreposição gene-a-gene estatisticamente frágil.

Ver `results/Proposta_tema_neuro_POP_SUI.xlsx` para a fundamentação completa, evidência da
literatura citada e roteiro metodológico proposto para expandir este trabalho.
