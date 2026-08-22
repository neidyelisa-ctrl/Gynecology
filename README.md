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

## Reanálise SEM filtro de baixa contagem (a pedido da usuária)

- `R/analysis_SUI_POP_nofilter.R` e `R/analysis_GO_nofilter.R` — reproduzem o mesmo
  cruzamento SUI(36h) x POP acima, mas **sem** o filtro de baixa expressão
  (`rowSums(contagem≥10) ≥ N amostras`) que o pipeline principal aplica antes do DESeq2.
  Motivação: a usuária identificou, pelo próprio histórico de análise, que **Cwh43** e
  **Calml5** (genes com contagem quase zero, só 1-2 réplicas com sinal) eram removidos por
  esse filtro antes mesmo de entrarem no teste estatístico.
- `results/SUI_x_POP_36h_SEM_FILTRO.xlsx` — resultado real e completo (6 abas).

**Resultado**: sem o filtro, **245 DEGs** no SUI 36h (181 com ortólogo humano) e **631 DEGs**
no POP — bem mais que com filtro (110 e 642 → note que o total de POP já era alto mesmo com
filtro; o salto grande é no SUI). **6 genes em comum** (ZCCHC12, INPP4B, ECM1, NTF3, BEND3,
KREMEN1) — o mesmo número que a usuária tinha encontrado na análise dela (sem filtro, sem
covariável de idade, ortólogos via BioMart), confirmando que o filtro de baixa contagem é
de fato o que explica a diferença de 4 → 6 genes. A composição exata difere um pouco porque
este pipeline mantém a metodologia própria (HomoloGene + covariável de idade no POP).

**Ressalva importante**: 4 dos 6 genes em comum (**ZCCHC12, INPP4B, NTF3, KREMEN1**) têm ≥2
de 6 amostras com contagem bruta ZERO no SUI 36h — **KREMEN1** inclusive com separação
perfeita 3 amostras positivas / 3 zeradas — exatamente o padrão que gera fold-change
inflado e p-valor artificialmente baixo (o problema técnico que a usuária identificou). Só
**ECM1** e **BEND3** são bem expressos nas 6 amostras sem nenhum zero. Recomendação: tratar
os 4 genes esparsos como candidatos a falso positivo técnico até confirmação por qPCR.
GO Biological Process (offline): nenhum termo significativo no SUI mesmo sem filtro; no POP
(lista maior) aparecem vários termos significativos de cálcio/músculo liso, mas zero em
comum com o SUI. KEGG e hub genes desta reanálise dependem de novas exportações do STRING
(listas em `results/STRING_input_SUI_36h_nofilter_genes.txt` e
`results/STRING_input_POP_nofilter_genes.txt`) — ainda não fornecidas.

## Causa raiz completa da diferença de genes (ex.: CALML5) entre a análise da usuária e a minha

A usuária forneceu o histórico real do console dela (`ANALYSIS_SUI.docx`). Comparando linha a
linha com o meu pipeline, além do filtro de baixa contagem (já resolvido acima), havia mais
duas diferenças em POP:
1. O modelo dela é `design = ~ condition` — **sem** a covariável de grupo etário (D/Y) que eu
   uso (`~ age_group + condition`).
2. Ela nunca chama `lfcShrink()` — usa o `log2FoldChange` bruto (MLE) de `results()`.

`R/analysis_SUI_POP_reproduce_userrecipe.R` testa isso isoladamente. Resultado decisivo
(`R/pop_age_shrink_4combos_diagnostic.R`, testado nos 7 genes de interesse — CWH43, INPP4B,
CALML5, KRT10, SERPINB2, DMKN, ZCCHC12): **`lfcShrink` não muda nenhum `padj`** (nunca muda — encolhimento
só afeta a estimativa de log2FC, não o teste de Wald) e não muda o veredito de significância
nesses genes. **A covariável de idade é a causa real**: com `~ age_group + condition`, só
INPP4B e ZCCHC12 desses 7 passam de padj<0,05; sem a covariável (`~ condition`, como no
script dela), **todos os 7 passam** — incluindo **CALML5** (padj vai de 0,17 → 0,033).

Rodando o pipeline completo com a receita exata dela (sem filtro, sem covariável de idade,
sem shrink, ortólogos via HomoloGene — não BioMart, por instrução da própria usuária): **5
dos 6 genes da lista dela batem exatamente** (CWH43, INPP4B, CALML5, SERPINB2, e mais
ZCCHC12 que não está na lista dela). Os únicos 2 que não bati são **KRT10** (não existe no
HomoloGene em nenhum grupo de ortólogo — falha real de cobertura da base) e **DMKN** (existe
no HomoloGene, mas o grupo de ortólogos dele não inclui rato — só humano/chimpanzé/macaco/
camundongo) — ambos só aparecem via BioMart, que tem uma base de mapeamento diferente. Não é
um bug em nenhum dos dois pipelines — é a diferença esperada entre duas fontes de ortologia
(HomoloGene vs. BioMart/Ensembl) e entre incluir ou não a covariável de idade no desenho do
POP. Ver commit com `R/analysis_SUI_POP_reproduce_userrecipe.R` para o script completo.

## Verificação independente de GO/KEGG para os 6 genes em comum da usuária

- `R/analysis_GO_KEGG_her6genes.R` — roda GO Biological Process e KEGG nos 6 genes da
  usuária (KRT10, SERPINB2, CALML5, CWH43, INPP4B, DMKN, extraídos de
  `thesis_proposal_draft.pdf`), do zero, sem copiar o resultado do Perplexity dela: teste
  hipergeométrico manual (org.Hs.eg.db + GO.db para GO; coluna `PATH` do org.Hs.eg.db para
  KEGG, já que `enrichKEGG`/KEGGREST exigem internet, indisponível aqui). Testado com dois
  universos (genoma inteiro anotado e só os genes realmente testados no meu DESeq2).
- `results/Verificacao_GO_KEGG_seus6genes.xlsx` — resultado completo (5 abas).

**Confirmado de forma independente** (sobrevive à correção BH nos dois universos, com apenas
6 genes de entrada — sinal de que não é ruído): KEGG **Phosphatidylinositol signaling
system** (hsa04070, via CALML5+INPP4B, p.adjust 0,00005–0,0025) e GO **epidermis
development** (via KRT10+CALML5, p.adjust 0,00005–0,0022). Termos relacionados aparecem logo
atrás: cornified envelope assembly (DMKN), keratinocyte differentiation (KRT10),
phosphatidylinositol/inositol phosphate metabolism (INPP4B).

**Não confirmado** com as ferramentas offline disponíveis: via de sinalização de estrogênio
(hsa04915) e PI3K/AKT (hsa04151) — nenhum dos 6 genes aparece ligado a essas vias no
mapeamento KEGG do `org.Hs.eg.db` (que é uma tabela antiga, congelada, não atualizada pelo
Bioconductor há anos por restrição de redistribuição da KEGG). Pode ser uma limitação da
minha fonte de dados, ou pode ser que esse resultado específico tenha vindo da rede
"expandida" (genes vizinhos adicionados via STRING para robustecer os hub genes) em vez dos
6 genes originais — recomendado à usuária conferir no próprio código qual lista de genes
alimentou esse `enrichKEGG()` específico.

**Confirmação externa real**: nenhum dos 6 símbolos aparece nominalmente no artigo de scRNA-seq
que a usuária enviou (Zhang et al. 2024, *Exp Cell Res* 442:114280 — scRNA-seq de parede
vaginal anterior em mulheres com SUI), mas o artigo relata, de forma totalmente independente
(dataset humano real, scRNA-seq, não bulk RNA-seq de rato), que os genes UP-regulados nas
células epiteliais do grupo SUI se concentram exatamente no mesmo eixo temático: *epidermis
development, epidermal cell differentiation, skin development, keratinocyte differentiation*
(seção 3.4, Fig. 5A do artigo), com queratinização aumentada confirmada por
imuno-histoquímica (KRT8). Três métodos independentes (pipeline da usuária, meu recálculo do
zero, e um estudo publicado peer-reviewed em tecido humano) convergem no mesmo eixo biológico
— diferenciação/cornificação epitelial alterada — o que reduz bastante a chance de os
achados da usuária serem um artefato do código. Ressalva que continua valendo: 6 genes é uma
amostra pequena para enriquecimento estatístico robusto; tratar como hipótese forte para
validação (qPCR/imuno-histoquímica de CALML5/KRT10/DMKN em tecido de POP humano), não como
achado definitivo.

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

## GO enrichment (offline) + próximo passo: hub genes e KEGG via STRING

- `R/analysis_GO_enrichment_and_STRING_export.R` — enriquecimento GO Biological Process
  rodado **sem depender de clusterProfiler/internet** (teste hipergeométrico manual via
  `org.Hs.eg.db` + `GO.db`, equivalente ao `enrichGO`), para **SUI 36h, SUI 72h e POP**
  separadamente. Resultado real em `results/GO_BP_SUI_36h.csv`, `results/GO_BP_SUI_72h.csv`
  e `results/GO_BP_POP.csv`; DEGs+ortólogos e GO de SUI consolidados em
  `results/SUI_36h_72h_DEG_ortologos_GO.xlsx`.
- Nenhum termo passou de padj<0.05 após correção BH em nenhum dos três (esperado com listas
  de DEG desse tamanho contra milhares de termos GO), mas os termos de menor p-valor bruto
  contam uma história temporal real no SUI:
  - **36h**: *sympathetic nervous system development*, *cellular response to BDNF
    stimulus*, *dendrite extension* — tema dominante: **sistema nervoso** (resposta aguda
    à lesão).
  - **72h**: *mitotic cell cycle phase transition*, *collagen catabolic process*,
    *extracellular matrix disassembly*, *mitotic spindle organization* — tema dominante:
    **ciclo celular / remodelação de matriz extracelular** (fase de reparo tecidual).
  - Essa mudança de tema 36h→72h é biologicamente coerente com a cronologia
    lesão→reparo descrita em Shen et al. 2024.
- `results/STRING_input_SUI_36h_genes.txt` (88 genes), `results/STRING_input_SUI_72h_genes.txt`
  (47 genes) e `results/STRING_input_POP_genes.txt` (624 genes) — listas prontas para colar
  em [string-db.org](https://string-db.org) ("Multiple proteins" → Homo sapiens → Search)
  para obter a rede PPI (hub genes) e o enriquecimento KEGG — string-db.org não é acessível
  diretamente deste ambiente. O KEGG do SUI 72h também foi verificado pela usuária no STRING:
  sem via significativa (mesmo padrão do SUI 36h).

## Hub genes (rede PPI) e comparação de vias SUI 36h x SUI 72h x POP

- `R/analysis_hub_genes_igraph.R` — calcula hub genes com `igraph` (grau, betweenness,
  força ponderada pelo `combined_score` do STRING — equivalente ao cytoHubba) a partir das
  redes exportadas do STRING pelo usuário (`results/STRING_network_SUI_36h.tsv`,
  `results/STRING_network_SUI_72h.tsv`, `results/STRING_network_POP.tsv`) e compara com o
  enriquecimento KEGG exportado (`results/KEGG_enrichment_POP.tsv`).
- `results/HubGenes_GO_KEGG_SUI_x_POP.xlsx` — resultado real e completo (5 abas: Resumo,
  Hub_genes_SUI_36h, Hub_genes_SUI_72h, Hub_genes_POP_top50, KEGG_POP).

**Achados reais:**
- Rede do SUI 36h: **23 nós / 18 interações** — esparsa (só 23 dos 88 genes de entrada
  tinham interação de alta confiança no STRING).
- Rede do SUI 72h: **10 nós / 10 interações** — ainda mais esparsa que a de 36h (10 dos 47
  genes de entrada). Hub gene principal: **CCNB1** (grau 4), seguido por **DLGAP5, UBE2C,
  KIF18A** (formam um cluster denso de ciclo celular/mitose), depois **L1CAM** (adesão
  neural/guia axonal) e **MMP11/MMP7** (remodelação de matriz).
- Rede do POP: **342 nós / 787 interações** — bem mais densa que as duas do SUI.
- **Nenhuma via GO ou KEGG significativa para o SUI em nenhum dos dois tempos** (confirmado
  tanto pelo STRING quanto pelo teste hipergeométrico offline deste repositório) — esperado
  dado o tamanho pequeno das listas e as redes esparsas, não uma falha da análise.
- **POP**: 9 vias KEGG significativas (FDR<0,05), dominadas por sinalização de cálcio/
  músculo liso/cardiomiócito (Calcium signaling, cGMP-PKG, Vascular smooth muscle
  contraction, Adrenergic signaling in cardiomyocytes) e metabolismo de purina.
- **Zero genes em comum entre qualquer par das três redes** (SUI 36h x SUI 72h, SUI 36h x
  POP, SUI 72h x POP) — nenhuma convergência de hub gene ou de via GO/KEGG entre SUI e POP
  com os dados atuais, em nenhum dos tempos. Isso reforça que a via de convergência mais
  defensável continua sendo a curadoria por literatura (aba de proposta de tema), não a
  rede PPI bruta.
- Curiosidade: 2 dos hub genes do SUI 36h (**AGRN**, **GPC4**) são genes de sinapse/junção
  neuromuscular, e o SUI 72h também tem um hub neural independente (**L1CAM**, adesão/guia
  axonal) — reforça o eixo neural mesmo dentro desta análise independente, mesmo com o
  tema dominante do 72h migrando para ciclo celular/ECM.
