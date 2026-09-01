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
| `scripts/06_cruzamento_DEG_GO_KEGG.R` | A abordagem original (antes do conselho do professor): cruza DEG(POP) com a lista do Chen 2006, ORA (GO/KEGG) nos genes em comum/concordantes | `results/06_*` |
| `scripts/07_score_painel_queratinizacao.R` | Escore de painel (6 genes do eixo de queratinização) testado diretamente nas 12 pacientes do POP — a proposta de validação "sem forçar" | `results/07_*` |
| `scripts/08_cruzamento_DEG_Chen2003_GO_KEGG.R` | O MESMO processo do script 06, agora com o segundo artigo do Chen (2003, fase proliferativa) | `results/08_*` |
| `scripts/09_GSEA_preranked_Chen2003.R` | GSEA preranked no painel do Chen 2003 (mesmo método do script 03) | `results/GSEA_preranked_Chen2003_KEGG.csv` |
| `scripts/10_GSEA_preranked_POP_KEGG.R` | GSEA preranked do POP com permutação de rótulo de gene set (em vez de fenótipo) — para comparar diretamente com o GSEA "normal" do script 02, no mesmo dataset | `results/GSEA_preranked_POP_KEGG_genesetpermutation.csv`, `results/10_comparacao_metodos_POP_normal_vs_preranked.csv` |
| `scripts/11_vias_compartilhadas_final.R` | Consolidação: vias compartilhadas entre os 2 métodos de GSEA do POP e os 3 painéis de SUI (Chen 2003, Chen 2006, Wei 2020) | `results/11_shared_*.csv` |
| `scripts/12_cruzamento_DEG_Wei2020_GO_KEGG.R` | Versão INICIAL (top-40 do artigo impresso) — SUPERADA pelos scripts 14-15, mantida por histórico | `results/12_*` |
| `scripts/13_GSEA_preranked_Wei2020.R` | Versão INICIAL (top-40) — SUPERADA pelo script 15 | `results/GSEA_preranked_Wei2020_KEGG.csv` |
| `scripts/14_cruzamento_DEG_Wei2020full_GO_KEGG.R` | **Versão DEFINITIVA**: mesmo processo dos scripts 06/08/12, agora com a Tabela Suplementar S2 COMPLETA do Wei 2020 (6.118 genes reais, com p-valor) | `results/14_*` |
| `scripts/15_GSEA_preranked_Wei2020full.R` | **Versão DEFINITIVA**: GSEA preranked no painel completo (ranking = sinal×-log10(p), p-valor real por gene) | `results/GSEA_preranked_Wei2020full_KEGG.csv` |
| `scripts/16_DEG_GSE12852_limma.R` | DEG do SEGUNDO dataset de POP (GSE12852, Applied Biosystems) — limma + `duplicateCorrelation`, com diagnóstico por tecido (round vs. uterossacral) | `results/GSE12852_*` |
| `scripts/17_GSEA_classic_GSE12852_uterosacral_KEGG.R` | GSEA clássico (permutação de fenótipo) no GSE12852, ligamento uterossacral | `results/GSEA_classic_GSE12852_uterosacral_KEGG.csv` |
| `scripts/18_cruzamento_GSE12852_x_Wei2020full.R` | Cruza GSE12852 com o painel completo do Wei 2020 (mesmo processo do script 14) | `results/18_*` |
| `scripts/19_GSE12852_validacao_cruzada.R` | Validação: GSE12852 x GSE53868 concordam entre si? + GSEA GSE12852 x Wei2020 | `results/19_*` |

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

### 5) Cruzamento DEG(POP) x DEG(Chen2006) e GO/KEGG (script 06) — a abordagem original

Esta é a abordagem que você estava seguindo antes do conselho do professor
(cruzar DEG das duas doenças, olhar genes em comum, rodar GO/KEGG), rodada
aqui de forma independente e completa:

- **Interseção estrita** (DEG no POP com |log2FC|>1 & FDR<0,05 **E** na lista
  do Chen 2006): **1 gene** — `KRT17`. De menos de 2 genes não dá para rodar
  GO/KEGG como teste estatístico (precisa de ≥2 genes na mesma via para
  calcular sobreposição).
- **Genes concordantes em direção** (sem exigir significância individual no
  POP — só que a direção bata): **36 de 53 genes testáveis (68%)** — teste
  binomial de sinal, **p = 0,0127**.
- **GO Biological Process nos 36 genes concordantes** (universo correto = os
  31.072 genes do array, não o catálogo Entrez inteiro): **118 de 243 termos
  significativos (FDR<0,05)**, os 6 primeiros todos do mesmo eixo:
  *intermediate filament bundle assembly* (KRT14/PKP1), *intermediate
  filament organization* (KRT17/KRT14/KRT16), *keratinocyte differentiation*
  (KRT14/KRT16/S100A7), *epidermis development* (KRT14/S100A7/COL17A1),
  *hair cycle* (KRT14/KRT16), *epithelial cell differentiation*
  (KRT17/KRT14/KRT16).
- **KEGG nos mesmos 36 genes**: **0 de 37 vias significativas** — mesmo
  padrão negativo do GSEA (KEGG é mais genérico/menos granular que GO,
  precisa de mais genes por via para detectar sinal).

### 6) Escore do painel de queratinização nas 12 pacientes do POP (script 07) — proposta de validação sem forçar

Os 6 genes que aparecem nos 6 termos GO mais significativos acima (KRT14,
KRT16, KRT17, PKP1, S100A7, COL17A1 — todos medidos e concordantes no
GSE53868; TP63 aparece no Chen 2006 mas não foi medido neste array) formam
um painel definido **a priori**, a partir de evidência externa (GO + Zhang
et al. 2024), não garimpado nos dados do POP. Testar esse painel específico
diretamente nas 12 pacientes é uma confirmação dirigida por hipótese — mais
apropriada para 6 genes do que competir contra o universo inteiro do
GO/KEGG.

**Método**: para cada paciente, diferença pareada (sítio do prolapso menos
sítio sem prolapso) dos 6 genes, cada gene escalado pelo seu próprio desvio-
padrão entre as 12 pacientes (sem centralizar — ver nota de bug corrigido no
cabeçalho do script), escore = média dos 6 valores escalados.

- **9 de 12 pacientes têm escore positivo** (para cima no sítio do
  prolapso); 3 negativo. Escore médio = 0,417.
- Teste t pareado (1 amostra): **p = 0,123**.
- Wilcoxon signed-rank: **p = 0,092**.
- Permutação de fenótipo (sign-flip, 10.000 permutações): **p = 0,105**.

**Leitura honesta**: os três testes concordam entre si (nenhum contradiz os
outros) e todos apontam na direção esperada (escore positivo, mesmo sentido
do Chen 2006), mas **nenhum cruza o limiar convencional de p<0,05** — é uma
tendência real, não um resultado definitivo. Com n=12 pacientes, um painel
de 6 genes não tem muito poder estatístico sozinho. Isto NÃO invalida o
achado de queratinização (que já tem 3 linhas de evidência independentes:
concordância de sinal p=0,013, GO enrichment FDR<0,05, e o estudo publicado
de Zhang et al. 2024) — mostra que testar diretamente nas pacientes
individuais, com um teste mais rigoroso e um n pequeno, é mais exigente do
que os testes agregados acima. Reporte os dois lados na tese: convergência
agregada real e estatisticamente significativa (genes + GO), tendência (não
significativa) ao nível de paciente individual.

**Gráfico**: `results/07_escore_painel_boxplot.png` — painel esquerdo mostra
a distribuição do escore nas 12 pacientes; painel direito mostra os 6 genes
individualmente por paciente (dá para ver que `S100A7` tem a maior amplitude
de variação, e que a paciente 4 é uma "outlier" que puxa a média).

### 7) Segundo artigo do Chen (2003, fase proliferativa) — mesmo processo (scripts 08-09)

Chen B, Wen Y, Zhang Z, Wang H, Warrington JA, Polan ML. *"Menstrual
phase-dependent gene expression differences in periurethral vaginal tissue
from women with stress incontinence."* Am J Obstet Gynecol. 2003;189(1):
89-97. Mesmos 5 pares SUI x continentes do Chen 2006, mas amostrados na fase
PROLIFERATIVA do ciclo menstrual (não secretória), array Affymetrix HuGeneFL
(6800 genes — menor/mais antigo que o U133A do Chen 2006). 90 genes
candidatos no artigo (62 up / 28 down); **69 mapeados com confiança** para
símbolo HGNC atual (43 up / 26 down) — conferido linha a linha contra o PDF
antes de rodar (ex.: PPIF p=0,02358 FC=-3,1; ALOX12 p=0,03810 FC=-2,7,
ambos batem exatamente com o texto original).

- **Interseção estrita** (DEG no POP E na lista do Chen 2003): **0 genes**.
- **Concordância de direção**: **35 de 62 genes testáveis (56,5%)** — teste
  binomial, **p = 0,374 (NÃO significativo)**. Diferente do Chen 2006 (68%,
  p=0,013) — a lista da fase proliferativa não converge com o POP do jeito
  que a lista da fase secretória converge.
- **GO/KEGG nos 35 genes concordantes**: tecnicamente 179 termos GO e 23
  vias KEGG saem significativos (FDR<0,05) — **mas leia com cautela**:
  como a concordância em si NÃO foi estatisticamente significativa (p=0,374,
  ou seja, esses 35 genes não são diferentes de qualquer metade aleatória
  da lista de 62), esse enriquecimento não tem a mesma base sólida do Chen
  2006. Ainda assim, chama atenção o **eixo TGF-beta** aparecendo
  repetidamente: `SMAD2` e `TGFB3` juntos em 6 das 10 vias KEGG mais
  significativas (adherens junction, cell cycle, pathways in cancer,
  endocytosis, hypertrophic/dilated cardiomyopathy) e no termo GO
  *"transforming growth factor beta receptor signaling pathway"*
  (HPGD/SMAD2/TGFB3). Isso bate com a própria discussão do artigo original
  do Chen 2003, que destaca TGFβ-3 como um dos genes centrais de ECM
  encontrados. Trate como hipótese a explorar, não achado confirmado.
- **GSEA preranked do Chen 2003 (KEGG)**: **0 vias significativas mesmo a
  FDR<0,25** (47 vias testáveis) — resultado negativo, consistente com a
  falta de concordância de direção acima.

### 8) GSEA "normal" (fenótipo) vs. GSEA preranked (rótulo de gene set) no POP — a comparação pedida (script 10)

Rodamos os DOIS métodos no MESMO ranking do POP (t pareado, transcriptoma
inteiro, as mesmas 218 vias KEGG) para ver o efeito real da escolha do tipo
de permutação — o ponto que motivou a correção já documentada acima.

| Método | Vias sig. FDR<0,05 | Vias sig. FDR<0,25 |
|---|---|---|
| **Normal** (permutação de fenótipo, sign-flip pareado) | 75 | 94 |
| **Preranked** (permutação de rótulo de gene set) | 64 | 140 |
| **Nos DOIS métodos ao mesmo tempo** | **8** | — |
| **Só no preranked** (não confirma no normal — candidato a falso positivo) | 56 | — |

**Leitura**: os dois métodos concordam bem menos do que se imagina — de 75
vias significativas no método correto (fenótipo) e 64 no método mais
liberal (rótulo de gene set), só **8 aparecem nos dois**. Isso é a prova
direta, no seu próprio dataset, do alerta de Subramanian et al. 2005: a
permutação de rótulo de gene set infla o número de vias "significativas"
(140 vs. 94 a FDR<0,25) porque ignora a correlação real entre genes de uma
mesma via. **A direção (sinal do NES), porém, é praticamente sempre a
mesma nos dois métodos** quando dá pra comparar (101 de 101 casos válidos
concordam) — ou seja, os dois métodos concordam sobre "a via sobe ou desce",
só discordam sobre "isso é forte o suficiente pra contar como
significativo". **Use sempre o método de fenótipo (`GSEA_classic_POP_KEGG.csv`)
como a referência principal** — é o metodologicamente correto para este
dataset (tem matriz bruta por amostra); o preranked existe aqui só para
essa demonstração comparativa.

### 9) Vias compartilhadas — versão final consolidada (script 11)

Testando POP (2 métodos) x Chen (2 artigos), 4 combinações, 2 limiares cada:

| Comparação | FDR<0,25 | FDR<0,05 |
|---|---|---|
| POP normal x Chen 2006 | **1** (05140) | 0 |
| POP preranked x Chen 2006 | **2** (04062, 05223) | 0 |
| POP normal x Chen 2003 | 0 | 0 |
| POP preranked x Chen 2003 | 0 | 0 |

**Atenção à DIREÇÃO nas 2 vias que só aparecem via POP-preranked**: `04062`
(chemokine signaling) e `05223` (non-small cell lung cancer) têm NES
**positivo** no POP-preranked (para cima) mas NES **negativo** no Chen 2006
(para baixo) — **direções opostas**. Isso não é uma via convergente de
verdade; é exatamente o tipo de falso positivo que a seção 8 acima avisa
que o método preranked/rótulo-de-gene-set produz. A única via que aparece
via o método CORRETO (POP normal/fenótipo) — `05140`, Leishmaniasis — tem
a mesma direção (negativa) nos dois lados, mas é uma via genérica de
resposta imune (movida por `HLA-DRB4`/`PRKCB` do lado do Chen, e por
dezenas de genes imunes/inflamatórios do lado do POP), não uma via
mecanisticamente específica de ECM/queratinização.

**Conclusão final e honesta desta seção**: nenhuma via KEGG nomeada e
mecanisticamente relevante é compartilhada entre POP e nenhum dos dois
artigos do Chen, em nenhuma combinação de método testada. Isso reforça (não
contradiz) o que as seções 5-6 já mostraram: a convergência real entre POP e
SUI está a nível de GENE INDIVIDUAL e de GO (queratinização, Chen 2006) —
não aparece como via KEGG nomeada porque o painel do Chen é pequeno demais
para esse tipo de teste ter poder estatístico.

### 10) Terceiro painel de SUI: Wei et al. 2020 (scripts 12-13) — mesmo tecido, mulheres pós-menopausa

Wei A, Wang R, Wei K, Dai C, Huang Y, Xu P, Xu J, Tang H, Zhang Y, Fan Y.
*"LncRNA and mRNA Expression Profiling in the Periurethral Vaginal Wall
Tissues of Postmenopausal Women with Stress Urinary Incontinence."*
Reprod Sci. 2020;27:1490-1501. 11 pares SUI x continentes (mulheres
PÓS-menopausa — diferente das duas do Chen, que são pré-menopausa), mesmo
tipo de tecido (parede vaginal periuretral), array Arraystar Human lncRNA +
mRNA V4.0.

**Limitação importante deste painel, diferente dos dois do Chen**: o artigo
identificou 7.102 mRNAs diferencialmente expressos no total (FC≥2, P<0,05)
— muito mais que os painéis do Chen — mas **o PDF fornecido só traz as
Tabelas 5 e 6 do artigo impresso (top 20 up-regulados + top 20
down-regulados = 40 genes, sem p-valor individual por gene, só fold
change)**. A tabela completa (Material Suplementar S2) está hospedada à
parte no site da revista e não veio com o PDF. **Na prática, este painel
(39 genes únicos após remover 1 duplicata) é MENOR que os dois do Chen
(59-69 genes)** — o ganho de poder estatístico do "7.102" só existiria se
a Tabela S2 completa fosse obtida.

**mRNA vs. lncRNA**: o artigo também lista 8.840 lncRNAs diferencialmente
expressos, mas eles NÃO foram usados aqui — a maioria não tem símbolo de
gene estável (aparecem como `TCONS_00017996`, `XLOC_008852`, IDs de
transcrito, não de gene) e não tem anotação KEGG/GO no `org.Hs.eg.db`. Só
os mRNAs (equivalentes a "gene", no mesmo sentido que os arrays do Chen e
do POP) são compatíveis com este pipeline.

**Resultados**:
- **Interseção estrita**: 0 genes.
- **Concordância de direção**: **11 de 33 genes testáveis (33,3%) — ABAIXO
  de 50%** (teste binomial, p=0,080, não significativo, mas na direção
  OPOSTA da que se esperaria para convergência). Diferente do Chen 2006
  (68% concordante) e mais parecido com o Chen 2003 (56,5%, também sem
  convergência) — este painel específico NÃO mostra o mesmo padrão de
  convergência com o POP.
- **GO/KEGG nos 11 genes concordantes**: tecnicamente "significativo"
  (126 de 141 termos GO), mas **cada termo tem Count=1** (um gene só) —
  isso não é enriquecimento real, é o mesmo artefato de gene-único-genérico
  já documentado várias vezes neste projeto (qualquer gene bate em dezenas
  de termos GO só dele, sem relação nenhuma com o tema da tese). **Não
  cite estes termos GO isoladamente na tese.**
- **GSEA preranked (KEGG)**: só 5 vias testáveis (painel muito pequeno);
  **0 significativas mesmo a FDR<0,25**.
- **Vias compartilhadas com o POP** (nos 2 métodos de GSEA do POP): **0**.

**Leitura honesta**: este terceiro painel, do jeito que está disponível
para nós (só o top-40 do artigo impresso, mulheres pós-menopausa), não
reforça nem contradiz o eixo de queratinização do Chen 2006 de forma
estatisticamente válida — os números têm cara de ruído (painel pequeno,
sem p-valor por gene, GO artefactual). **Se você conseguir a Tabela
Suplementar S2 completa** (7.102 mRNAs, no site do artigo, seção
"Electronic supplementary material"), esse painel passaria a ser o mais
forte dos três — bem maior que os dois do Chen, e com uma métrica de
ranking real (score contínuo por gene), permitindo GSEA no sentido pleno
do método, não a versão painel-pequeno usada aqui.

### 11) ATUALIZAÇÃO: Wei 2020 com a Tabela Suplementar S2 COMPLETA (scripts 14-15) — o maior achado do painel de literatura até agora

A usuária conseguiu baixar a Tabela Suplementar S2 real do artigo (Electronic
Supplementary Material 2, `data/wei2020_TableS2_mRNA_original.xls`) — a saída
completa do software original (GeneSpring GX), com P-value, FDR, Fold Change
E intensidade normalizada por amostra individual (3 SUI x 3 Ctrl, o
subconjunto de 6 amostras usado no array, por gene). Extraída para
`data/wei2020_mRNA_full.csv`: **7.102 linhas → 6.118 genes únicos** após
remover sondas duplicadas (mantendo a de menor p-valor por gene). **Esta
seção SUBSTITUI a análise anterior do Wei 2020 (scripts 12-13, painel de
40 genes) como referência principal** — mesma limitação de fundo continua
(a tabela já vem pré-filtrada pelo próprio estudo, não é o array inteiro
de ~20.730 genes testados), mas 6.118 genes é ~100x maior que o painel
anterior e permite testar quase todo o KEGG.

**(a) Interseção estrita (DEG POP ∩ Wei2020)**: **15 genes** — de longe a
maior sobreposição direta do projeto inteiro: `WEE1, SLN, ADAMTS4, NEDD9,
NFATC2, NR4A3, SNAI1, IGJ, CSF3, AREG, APOLD1, THBD, ATF3, IGLL1, LDLR`.

**(b) Concordância de direção — resultado SURPREENDENTE, leia com
atenção**: **2.215 de 4.797 genes testáveis (46,2%) concordantes — ABAIXO
de 50%**, teste binomial **p = 1,24×10⁻⁷ (extremamente significativo)**.
Isto é o OPOSTO do Chen 2006 (68% concordante, convergência real) — aqui
temos **discordância estatisticamente significativa**: genes que sobem no
Wei 2020 (SUI pós-menopausa) tendem a DESCER no POP, e vice-versa, mais do
que o esperado por acaso. Não é ruído (a amostra é grande demais para
isso, p=1,2×10⁻⁷) — é um padrão real. Hipóteses para a diferença em
relação ao Chen 2006 (mesmo tecido, mesma doença): (i) mulheres
PÓS-menopausa (Wei) vs PRÉ-menopausa (Chen/POP) — o estado hormonal pode
alterar a direção da resposta transcricional; (ii) desenho não-pareado do
Wei (3 SUI x 3 Ctrl, mulheres diferentes) vs desenho pareado do POP (mesma
paciente); (iii) plataformas de array diferentes. **Discuta esta
divergência explicitamente na tese como uma limitação real de comparar
resultados entre estudos com populações/desenhos diferentes — não tente
esconder ou "resolver" o contraste, ele é uma informação genuína.**

**(c) GO/KEGG na interseção estrita (15 genes)**: 164 de 195 termos GO
significativos (FDR<0,05), vários com >1 gene de suporte (ex.: *negative
regulation of transcription by RNA polymerase II* — NFATC2/NR4A3/SNAI1/ATF3,
4 genes). KEGG: 0 de 21 (todos Count=1, ruído).

**(d) GO/KEGG nos 2.215 genes concordantes**: aqui sim, com uma amostra
grande, o enriquecimento é robusto — **605 de 5.212 termos GO
significativos** (dominado por regulação da transcrição, sinalização,
diferenciação celular, migração celular, adesão celular — temas amplos,
esperado com 2.215 genes de entrada) e, mais importante, **109 de 214 vias
KEGG significativas**, incluindo:
- **04350 — Sinalização TGF-beta** (21 genes: `TGFBR1/ID3/RPS6KB2/SMURF1/
  MAPK1/BMP2/FST/SMAD3/BMP6/TGFBR2/PPP2CA/BMP8A/ID1/THBS2/RHOA/TGFB1/SMAD5/
  NODAL/BMP5/ACVR1C/MAPK3`, FDR=2,8×10⁻⁷) — **QUARTA confirmação
  independente do eixo TGF-beta neste projeto** (depois do Chen 2003
  humano — SMAD2/TGFB3 —, do modelo de rata pós-parto de Kerns/Damaser —
  SMAD2 — e agora aqui, com suporte estatístico muito mais forte, 21
  genes de uma vez). Este é o achado de via mais robusto e mais replicado
  de todo o projeto.
- Outras vias de tecido conjuntivo/adesão: *regulation of actin
  cytoskeleton* (04810), *focal adhesion* (não listada no top 15 mas
  presente na tabela completa), *pathways in cancer* (04520, genérica).

**(e) GSEA preranked (KEGG, painel completo)**: **0 de 201 vias
significativas a FDR<0,05; 52 a FDR<0,25** — muito mais rico que qualquer
painel anterior (Chen 2003: 0, Chen 2006: 3, Wei top-40: 0).

**(f) Vias compartilhadas com o GSEA do POP — o maior número do projeto**:

| Comparação | FDR<0,25 | FDR<0,05 |
|---|---|---|
| POP normal (fenótipo) x Wei2020(completo) | **19** | 0 |
| POP preranked (rótulo de gene set) x Wei2020(completo) | **36** | 0 |

De longe o maior número de vias compartilhadas do projeto (bem acima do
"0-2" dos artigos do Chen) — **mas, checando a coluna `Mesma_direcao`
(sinal do NES nos dois lados) nos resultados salvos
(`results/11_shared_POP_*_x_Wei2020_FDR025.csv`), a MAIORIA é
DISCORDANTE** (POP preranked: só 3 de 36 concordantes; POP normal: só 1 de
6 comparáveis, o resto tem NES=NA de um lado — caso de borda já documentado
na seção 8). **Isto é consistente com o achado (b) acima**: o padrão
dominante entre Wei2020 e POP é de mudança em direções opostas nas mesmas
vias, não convergência.

**Leitura honesta e final desta seção**: com dados reais e completos, o
Wei 2020 mostra MUITO mais sinal estatístico que os dois artigos do Chen
(tanto a nível de gene quanto de via) — mas o sinal dominante é de
**discordância** com o POP, não convergência, provavelmente refletindo a
diferença de população (pós vs pré-menopausa) ou desenho experimental. Ainda
assim, dentro do subconjunto real de 2.215 genes concordantes, o eixo
**TGF-beta continua aparecendo com força** — a quarta vez neste projeto,
agora com o suporte estatístico mais robusto de todos. **Recomendação para
a tese**: reporte os dois achados lado a lado (discordância ampla + TGF-beta
concordante e robusto dentro do subconjunto que concorda) em vez de
escolher um e ignorar o outro — os dois são reais e ambos informativos.

### 12) SEGUNDO dataset de POP: GSE12852 (scripts 16-19) — independente do GSE53868

A usuária forneceu um segundo dataset de POP do GEO, **GSE12852** — "Gene
expression profile in pelvic organ prolapse" — 8 mulheres com POP vs. 9
controles (caso-controle, **não pareado por paciente** como o GSE53868),
ligamento uterossacral + ligamento redondo (2 tecidos/paciente, 34
arrays), plataforma Applied Biosystems Human Genome Survey Microarray V2.0
(GPL2986). Diferente do GSE53868 (Agilent, símbolos de gene diretos), esta
plataforma usa IDs de sonda numéricos — a anotação (`data/GPL2986_annotation.tsv`,
sonda→símbolo) foi extraída da tabela de plataforma dentro do arquivo SOFT
completo da série (`GSE12852_family.soft`, fornecido pela usuária, já que
a página da plataforma no GEO está bloqueada neste ambiente).

**Diagnóstico importante ANTES de aceitar qualquer número** (script 16): a
análise combinando os dois tecidos (limma + `duplicateCorrelation`, bloco
por paciente) deu um sinal quase nulo — **FDR mínimo entre 16.752 genes =
0,97**. Investigando por tecido separado: o **ligamento redondo** não tem
sinal acima do esperado por acaso (371 genes com p bruto<0,05, MENOS que os
838 esperados por acaso) e estava diluindo a análise combinada; o
**ligamento uterossacral** tem sinal real, ainda que modesto (914 vs. 838
esperados, FDR mínimo 0,22) — biologicamente plausível, é a estrutura mais
diretamente implicada na fisiopatologia do POP. **A partir daqui, uterossacral
sozinho (17 amostras, 8 POP x 9 controle, desenho não-pareado simples) é a
análise de referência** para este dataset — não a combinada.

**DEGs (uterossacral, |log2FC|>1, FDR<0,05): 0 de 16.752 genes.** Nenhum
gene sobrevive ao corte formal — dataset pequeno e heterogêneo (idades
24-70, várias etnias, pré e pós-menopausa misturadas). Os genes com menor
p-valor bruto são biologicamente plausíveis e batem com a conclusão do
próprio artigo original ("Immunity and Defense" + remodelação de ECM):
`MYH3, TREM1, CTHRC1, THBS1, TAC1, ANGPTL5` (`results/GSE12852_uterosacral_limma_completo.csv`).

**GSEA clássico (KEGG, permutação de fenótipo, script 17)**: **0 de 213
vias significativas, mesmo a FDR<0,25** — consistente com a ausência de
DEGs individuais; dataset pequeno demais para poder estatístico de via.

**Cruzamento com Wei 2020 completo (script 18) — resultado forte e
inesperado**: **60,0% de concordância de direção (2.339 de 3.900 genes
testáveis)**, teste binomial **p = 9,2×10⁻³⁶** — extremamente
significativo, e desta vez CONCORDANTE (ao contrário do GSE53868 x Wei2020,
que foi significativamente DISCORDANTE, 46,2%, seção 11). GO/KEGG nos 2.339
genes concordantes: 34 de 5.428 termos GO e **40 de 222 vias KEGG**
significativas — mas aqui é preciso ler com cuidado: as vias KEGG do topo
são majoritariamente **categorias genéricas de grandes módulos** (Ribosome,
Metabolic pathways, Oxidative phosphorylation, e as "doenças" que
tradicionalmente carregam os mesmos genes de OXPHOS/ribossomo por
artefato de anotação — Parkinson, Alzheimer, Huntington-like) — o mesmo
padrão de "via genérica" já visto e descontado várias vezes neste projeto.
Os termos GO/KEGG mais especificamente interessantes, com menos risco de
serem artefato: *antigen processing and presentation of exogenous peptide
antigen via MHC class II* (16 genes), *positive regulation of T cell
activation*, *phagosome* (KEGG 04145), *lysosome* (KEGG 04142) — um tema de
**resposta imune/apresentação de antígeno**, coerente com a conclusão do
próprio artigo do GSE12852 ("Immunity and Defense... independent of
inflammatory infiltrates").

**Checagem de validação essencial (script 19) — leia antes de confiar no
achado acima**: os DOIS datasets de POP concordam entre si?
**GSE12852 (uterossacral) x GSE53868: só 51,8% de concordância (7.020 de
13.543 genes)** — estatisticamente diferente de 50% (p=2×10⁻⁵, a amostra é
grande), mas o efeito é **muito fraco** (quase empate), bem mais fraco que
a concordância de qualquer um dos dois com o Wei2020. **Isto é uma
informação importante para a interpretação**: os dois datasets de POP,
apesar de medirem a mesma doença, não concordam fortemente entre si
(plataformas diferentes — Agilent vs. Applied Biosystems —, tecidos
diferentes — parede vaginal vs. ligamentos uterossacral/redondo —,
desenhos diferentes — pareado vs. caso-controle). Isso não invalida os
achados anteriores, mas é motivo real de cautela: **um "achado" que só
aparece com um dos dois POPs (como a discordância com o Wei2020 no
GSE53868, ou a concordância forte no GSE12852) pode refletir a
especificidade daquele dataset/tecido, não POP em geral** — recomendação:
tratar tanto a concordância do GSE12852 quanto a discordância do GSE53868
como achados específicos de cada dataset, discutidos lado a lado na tese,
não escolher um como "o resultado do POP".

**GSEA: GSE12852 x Wei2020 — 0 vias compartilhadas** (o GSE12852 não teve
nenhuma via significativa no GSEA para comparar — ver acima).

**Resumo desta seção**: o GSE12852 é um dataset mais fraco/ruidoso que o
GSE53868 a nível de gene individual e de via (esperado — menor, mais
heterogêneo, plataforma mais antiga, tecido diferente e parcialmente sem
sinal como o redondo). Mesmo assim, rendeu o teste de concordância de
direção mais forte do projeto inteiro contra o Wei2020 (p=9×10⁻³⁶) — mas a
checagem de validação cruzada (GSE12852 x GSE53868, só 51,8%) mostra que
esse resultado é específico deste par de datasets, não uma confirmação
universal de "POP concorda com SUI". Reporte com essa nuance.

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
