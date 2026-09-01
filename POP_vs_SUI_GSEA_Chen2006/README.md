# DEG (POP, GSE53868) + GSEA preranked (SUI, Chen 2006) x GSEA clássico (POP) — vias compartilhadas

Pasta autocontida com a análise pedida: DEGs do POP (logFC>1, FDR<0,05), GSEA
preranked nos DEGs do artigo de Chen et al. (SUI) e GSEA clássico (transcriptoma
inteiro) nos genes do POP, e a comparação de vias compartilhadas entre os dois.
Todo o código roda do zero a partir dos dados em `data/` — não depende do
resto do repositório.

## Os dois conjuntos de dados

- **POP**: `data/GSE53868_series_matrix.txt` — GSE53868, "Micro-array analysis
  of the anterior vaginal wall from premenopausal patients with Pelvic Organ
  Prolapse (POP)" (Kerkhof et al., VU University Medical Centre). 12 mulheres
  com POP, biópsia **pareada** por paciente: tecido do sítio do prolapso vs.
  tecido da mesma paciente fora do prolapso. Agilent 4x44K (GPL18142), já
  normalizado (log2), símbolos de gene como ID_REF.
- **SUI**: `data/chen2006_79genes.csv` — Chen B, Wen Y, Zhang Z, Guo Y,
  Warrington JA, Polan ML. *"Microarray analysis of differentially expressed
  genes in vaginal tissues from women with stress urinary incontinence
  compared with asymptomatic women."* Hum Reprod. 2006;21(1):22-29 (publicado
  online em ago/2005 — é o PDF enviado como "chen2005_DEG.pdf"). 5 pares SUI
  x continentes, parede vaginal periuretral, array Affymetrix U133A. O artigo
  reporta 79 DEGs finais (Tabelas II/III, comuns aos algoritmos MAS 5.0 e
  RMA); **60 linhas foram mapeadas com confiança para símbolo HGNC atual**
  (59 símbolos únicos após remover 1 duplicata de sonda — PI3/elafin aparece
  em 2 sondas do array original). As outras 19 permanecem genuinamente
  ambíguas na nomenclatura de 2005 ("hypothetical protein FLJxxxxx",
  "KIAAxxxx protein", "Zinc finger protein" genérico) e foram excluídas —
  não inventamos símbolo para elas.

## O que cada script faz (rode em ordem, a partir desta pasta)

| Script | O que faz | Saída |
|---|---|---|
| `scripts/01_DEG_POP_limma.R` | limma pareado (design `~individual+tissue`) no GSE53868 | `results/GSE53868_limma_completo.csv` (todos os genes), `results/GSE53868_DEG_logFC1_FDR05.csv` (DEGs) |
| `scripts/02_GSEA_classic_POP_KEGG.R` | GSEA clássico (transcriptoma inteiro, ranking = t pareado), vias KEGG, **permutação de fenótipo** (sign-flip pareado, 500 permutações) | `results/GSEA_classic_POP_KEGG.csv` |
| `scripts/03_GSEA_preranked_Chen2006.R` | GSEA preranked nos 59 genes do Chen 2006 (ranking = sinal×-log10(p)), vias KEGG, permutação de rótulo de gene set (1000 permutações) | `results/GSEA_preranked_Chen2006_KEGG.csv` |
| `scripts/04_shared_pathways.R` | Compara os dois resultados de GSEA, reporta vias com o mesmo ID de via KEGG significativo nos dois lados | `results/GSEA_shared_KEGG_FDR025.csv` / `FDR005.csv` (ou `_VAZIO.txt` se nenhuma) |
| `scripts/05_OPCIONAL_upgrade_msigdbr_fgsea.R` | **Rode no seu computador**, com internet normal — versão com `fgsea`+`msigdbr` (MSigDB completo: Hallmark+KEGG+Reactome, nomes de via legíveis) | `results/fgsea_*` |

Pré-requisitos (scripts 1-4, já resolvidos neste ambiente via `apt`, mas
seguem os nomes padrão do Bioconductor para quem for rodar em outro lugar):
`limma`, `org.Hs.eg.db`. Nenhum dos dois precisa de internet em tempo de
execução.

## Resultados (rodados e conferidos nesta sessão)

### 1) DEGs do POP — |log2FC| > 1 e FDR < 0,05

**117 de 31.072 genes testados** (73 para cima no sítio do prolapso, 44 para
baixo) — `results/GSE53868_DEG_logFC1_FDR05.csv`.

### 2) GSEA clássico do POP (KEGG, transcriptoma inteiro, permutação de fenótipo)

- 218 vias KEGG testadas (5-200 genes cada).
- **75 de 218 significativas a FDR<0,05** (94 a FDR<0,25).
- Ressalva: com 500 permutações, muitas vias empatam no p-valor mínimo
  possível (1/501 ≈ 0,002) — não dá para refinar esse ranking sem rodar mais
  permutações (o script 05, com `fgsea`, não tem esse teto).
- Plausível dado o desenho do GSE53868: compara tecido do sítio de POP vs.
  sítio sem POP **na mesma paciente** — uma diferença ampla de remodelação
  tecidual pode legitimamente afetar muitas vias ao mesmo tempo, não é sinal
  de erro.

### 3) GSEA preranked do Chen 2006 (SUI, 59 genes, KEGG)

- Só 19 das 229 vias KEGG têm ≥2 dos 59 genes presentes (a maioria das vias
  não é testável com uma lista tão curta).
- **0 de 19 significativas a FDR<0,05; 3 a FDR<0,25** (05223, 05140, 04062 —
  todas movidas por só 2 genes cada, `PRKCB` + `FOXO3`/`HLA-DRB4`).
- **Leia com cautela**: isto é um uso não-padrão do GSEA (painel curto,
  pré-filtrado pelos próprios autores do artigo, não o transcriptoma
  original do estudo deles — que não foi disponibilizado). A permutação
  usada (rótulo de gene set) é mais liberal que a permutação de fenótipo
  usada no POP acima. Ver nota metodológica completa no cabeçalho do script
  03.

### 4) Vias compartilhadas — a resposta à pergunta principal

| Limiar | POP significativo | Chen2006 significativo | **Compartilhadas** |
|---|---|---|---|
| FDR < 0,25 | 94 | 3 | **1** (via KEGG `05140`) |
| FDR < 0,05 | 75 | 0 | **0** |

A única via que bate nos dois lados a FDR<0,25 é a **05140 (Leishmaniasis,
categoria "doença infecciosa" do KEGG)** — no POP ela é movida por 63 genes,
majoritariamente imunes/inflamatórios (TLR2/TLR4, complemento, HLA classe
II, TNF/IL1B/IL4/IFNG, MAPKs); no Chen 2006 só por 2 genes (`PRKCB`,
`HLA-DRB4`). **Não é um achado mecanístico defensável para citar
isoladamente** — é uma via "guarda-chuva" de resposta imune/inflamatória
genérica do KEGG, capturada por qualquer lista com 2 genes imunes
quaisquer, o mesmo padrão de ruído de via-única-gene já documentado no
projeto principal (ver `README.md` do repositório). A FDR<0,05 (o limiar
mais rigoroso, e o mesmo usado para os DEGs), **não há nenhuma via KEGG
nomeada compartilhada** entre o GSEA do POP e o GSEA preranked do Chen 2006.

**Isto é uma resposta negativa real, não um erro de código** — e é
consistente com o resto deste projeto de pesquisa (ver contexto abaixo):
convergência a nível de VIA NOMEADA entre um transcriptoma inteiro (POP,
~31 mil genes) e um painel curto pré-filtrado (Chen 2006, 59 genes) é um
critério estatisticamente muito mais difícil de bater do que convergência a
nível de GENE INDIVIDUAL ou de termo GO — porque uma via KEGG só "conta"
como testável no lado do Chen 2006 se ≥2 dos seus membros, por acaso,
estiverem entre os 59 genes daquela lista curta.

## Contexto complementar (achados reais, já validados antes neste repositório)

Embora a comparação GSEA-vs-GSEA acima não encontre vias KEGG nomeadas em
comum, dois métodos complementares — mais adequados para uma lista curta
como a do Chen 2006 — **encontram convergência real e robusta**:

- **Nível de gene individual** (`results/CONTEXTO_genes_individuais_Chen2006_x_GSE53868.csv`):
  dos 59 genes do Chen 2006 testáveis no GSE53868, **2 são individualmente
  significativos** no POP (FDR<0,05, |log2FC|>0,5) — `SERPINB8` e `KRT17`,
  **ambos na mesma direção** relatada no Chen 2006 (para cima). Olhando o
  painel inteiro (não só os 2 significativos): **36 de 53 genes testáveis
  (68%) concordam em direção** entre SUI (Chen 2006) e POP — teste binomial
  de sinal **p = 0,0127**, estatisticamente improvável por acaso.
- **Nível de via GO** (`results/CONTEXTO_GO_BP_genes_concordantes.csv`, GO
  Biological Process, teste hipergeométrico com universo correto = os
  31.072 genes do array): rodado nos genes que são ao mesmo tempo do painel
  Chen 2006 E concordantes em direção com o GSE53868 — **118 de 243 termos
  GO significativos (FDR<0,05)**, dominados por um eixo temático único e
  consistente: **queratinização/diferenciação epitelial** — *intermediate
  filament bundle assembly* (KRT14/PKP1), *intermediate filament
  organization* (KRT17/KRT14/KRT16), *keratinocyte differentiation*,
  *epidermis development*, *epithelial cell differentiation*.

Esse eixo de queratinização (via KRT14/KRT16/KRT17/PKP1/S100A7/COL17A1/TP63)
é o achado central deste projeto de pesquisa mais amplo — replicado de forma
independente em pelo menos 7 análises diferentes ao longo do trabalho (ver
`../README.md` do repositório para o histórico completo), incluindo uma
confirmação externa real e totalmente independente: Zhang et al. 2024 (*Exp
Cell Res*, scRNA-seq humano de parede vaginal em SUI) relata exatamente o
mesmo eixo (*epidermis development*, *keratinocyte differentiation*) nas
células epiteliais do grupo SUI, com queratinização confirmada por
imuno-histoquímica (KRT8).

**Leitura honesta e resumida**: GSEA "via nomeada contra via nomeada" (a
pergunta feita) dá uma resposta negativa/fraca (0 vias a FDR<0,05) — mas
isso é esperado dado o tamanho do painel do Chen 2006, não uma ausência
real de convergência biológica. A nível de gene individual e de GO
(critérios mais sensíveis para uma lista curta), a convergência SUI-POP é
real, robusta e tematicamente específica (queratinização), não genérica.

## Limitações deste ambiente (leia antes de levar os números para a tese)

- **Sem acesso a CRAN, Bioconductor, MSigDB, KEGG REST, Reactome ou Enrichr
  ao vivo** neste sandbox (confirmado por teste direto: `CONNECT` bloqueado
  para `cloud.r-project.org`, `bioconductor.org`, `zenodo.org` — inclusive a
  fonte de dados oficial que o próprio pacote `msigdbr` usa hoje —,
  `rest.kegg.jp`, `reactome.org`, `maayanlab.cloud`). Por isso os scripts
  1-4 usam **apenas KEGG**, e a versão **congelada** dentro do pacote
  `org.Hs.eg.db` (não a mais atual do KEGG.jp).
- **Sem nomes de via legíveis** nos scripts 1-4 (só o código numérico KEGG,
  ex. `05140`) — `org.Hs.eg.db` não inclui uma tabela de nomes de via. Para
  ver o nome, consulte `https://www.kegg.jp/pathway/hsa<ID>` no seu
  navegador, ou rode o script 05 (`fgsea`+`msigdbr`), que já traz nomes.
- **`scripts/05_OPCIONAL_upgrade_msigdbr_fgsea.R` não foi executado aqui**
  (precisa de internet normal) — está pronto para você rodar no seu
  computador/RStudio e obter a versão definitiva, com Hallmark+KEGG+Reactome
  completos e atualizados, no lugar do KEGG congelado usado acima. É
  esperado que o padrão qualitativo se mantenha (poucas/nenhuma via nomeada
  compartilhada, dado o painel pequeno do Chen 2006), mas os números exatos
  podem mudar com a base mais completa.
- **GSEA no painel do Chen 2006 é um uso não-padrão do método** (ver nota no
  cabeçalho do script 03) — trate os números da seção 3 como um sinal
  exploratório de baixo poder, não no mesmo patamar do GSEA do POP (que usa
  o transcriptoma completo com permutação de fenótipo, a forma
  metodologicamente correta).
- **Variabilidade de permutação**: os p-valores dos scripts 03-04 usam 500-
  1000 permutações; rodar de novo com uma semente (`set.seed`) diferente
  pode mover 1-2 vias para dentro/fora do limiar FDR<0,25 nas bordas (ex.:
  a contagem de "vias significativas do Chen2006 a FDR<0,25" variou entre 3
  e 5 em duas rodadas desta mesma análise, com pequenas diferenças na lista
  de genes de entrada) — não mude a conclusão principal (nenhuma via
  temática de ECM/queratinização aparece do lado do KEGG-GSEA do Chen 2006,
  e nada sobrevive a FDR<0,05).
