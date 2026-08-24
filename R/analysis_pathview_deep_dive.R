# =============================================================================
# Tentativa de usar pathview (Bioconductor) para aprofundar as vias GO/KEGG dos
# 6 genes validados (CWH43, INPP4B, CALML5, KRT10, SERPINB2, DMKN), e alternativa
# quando pathview nao funciona neste ambiente.
#
# COMO PATHVIEW FUNCIONA: ele baixa o diagrama KGML/PNG da via diretamente do
# KEGG (rest.kegg.jp) a cada chamada - nao ha modo 100% offline sem esse
# arquivo. Testado e CONFIRMADO bloqueado (mesma politica de rede de todo o
# projeto): CONNECT tunnel 403 para rest.kegg.jp/get/hsa04070/kgml. Ou seja,
# mesmo instalando o pacote (que exigiria compilar do zero via GitHub, ja que
# nao esta disponivel via apt nem CRAN/Bioconductor neste sandbox), a chamada
# pathview() falharia no mesmo ponto de sempre: download do KEGG.
#
# ALTERNATIVA USADA (a pedido explicito da usuaria, "se esse pacote nao
# funcionar, use outras maneiras"): em vez de colorir o diagrama, fui direto
# na bioquimica primaria de cada gene/enzima (a mesma fonte de dados que o
# KEGG usa para desenhar o mapa) e cruzei com a DIRECAO real de expressao
# (log2FC) de cada gene nos seus dados validados (results/SUI_36hr_ensembl_x_
# POP_common_genes.csv) para inferir a direcao provavel de cada metabolito/
# lipideo river abaixo. Isso da uma resposta MAIS especifica do que o
# pathview daria (que so pinta o gene de vermelho/verde, sem inferir o
# metabolito), mas e uma INFERENCIA mecanistica a partir de expressao genica -
# nao e dado de lipidomica medido diretamente (esse existe separadamente, ver
# results/Datasets_e_lipidios_vias.xlsx, aba Lipidomica_real_Zhang2025).
# =============================================================================

cat("Testando se o pathview consegue baixar o KGML da via KEGG (mesmo teste de\n",
    "rede feito para todo o projeto - ja confirmado bloqueado por policy 403):\n")
test <- tryCatch({
  con <- url("https://rest.kegg.jp/get/hsa04070/kgml", open = "rt")
  readLines(con, n = 1)
  close(con)
  "OK - conseguiu baixar"
}, error = function(e) paste("FALHOU:", conditionMessage(e)))
cat(test, "\n\n")

if (grepl("FALHOU", test)) {
  cat("Confirmado: pathview nao pode funcionar neste ambiente (download do KEGG\n",
      "bloqueado). Prosseguindo com a analise mecanistica manual abaixo, que usa\n",
      "a direcao real (log2FC) de cada gene validado.\n\n")
}

## -----------------------------------------------------------------------
## Direcao real de expressao dos 6 genes validados (SUI 36h e POP)
## -----------------------------------------------------------------------
common <- read.csv("results/SUI_36hr_ensembl_x_POP_common_genes.csv")
print(common[, c("Human_Ortholog_Symbol", "SUI_dir", "POP_dir", "same_direction")])

## -----------------------------------------------------------------------
## Anotacao mecanistica manual: para cada via GO/KEGG relevante, o gene que a
## dirige, a reacao/funcao bioquimica exata, e a inferencia de direcao do
## metabolito/lipideo dado o sentido real da expressao genica.
## Fontes bioquimicas primarias por gene (nao inferidas por mim, ver
## comentarios de cada bloco):
##  - INPP4B: PMC3248162, PMC7136497 (revisoes sobre INPP4B e PI3K/AKT)
##  - CALML5: cusabio.com/pathway/Phosphatidylinositol-signaling-system.html
##  - CWH43: PubMed 17761529 (Cwh43p remodela ancora GPI)
##  - SERPINB2: literatura classica de PAI-2/keratinocito (ver
##    results/Literatura_recente_apoio_tese.xlsx e conversas anteriores)
## -----------------------------------------------------------------------
pathway_deep_dive <- data.frame(
  Via = c(
    "KEGG hsa04070 - Phosphatidylinositol signaling system",
    "GO/KEGG - Inositol phosphate metabolism / phosphatidylinositol biosynthesis/dephosphorylation",
    "GO - GPI anchor biosynthetic process",
    "GO - Epidermis development / cornified envelope assembly / keratinocyte differentiation / intermediate filament organization / epithelial cell differentiation / peptide cross-linking / protein heterotetramerization",
    "GO - Fibrinolysis / negative regulation of endopeptidase activity / negative regulation of apoptotic process",
    "~14 vias KEGG genericas (Phototransduction, Glioma, Long-term potentiation, Gastric/Salivary secretion, GnRH signaling, Melanogenesis, Oocyte meiosis, Vascular smooth muscle contraction, Neurotrophin signaling, Insulin signaling, Alzheimer disease, Calcium signaling, Olfactory transduction) + Staphylococcus aureus infection + Amoebiasis"
  ),
  Gene_motor = c(
    "INPP4B (down/down) + CALML5 (up/up)",
    "INPP4B (down/down)",
    "CWH43 (up/up)",
    "KRT10 (up/up) + DMKN (up/up) + CALML5 (up/up)",
    "SERPINB2 (up/up)",
    "CALML5 (12 vias) / KRT10 (1) / SERPINB2 (1) - Count=1 em cada, gene unico"
  ),
  Alteracao_especifica_inferida = c(
    "INPP4B para de degradar PI(3,4)P2 -> PI(3)P na mesma proporcao (enzima reduzida): PI(3,4)P2 tende a ACUMULAR (subir); PI(3)P (produto) tende a CAIR. Como INPP4B tambem limita a sinalizacao PI3K/AKT (e classificado como supressor tumoral via esse mecanismo em varios canceres), a queda de INPP4B aponta para sinalizacao PI3K/AKT POTENCIALMENTE MAIS ATIVA, nao necessariamente deficiente - vale refinar essa nuance no rascunho. CALML5 nao mexe em lipideo, e sensor de Ca2+ 2 passos depois (via PLC/IP3).",
    "Mesma logica do INPP4B acima: menos INPP4B = menos producao de PI(3)P e de fosfatos de inositol soluveis rio abaixo (Ins(1,3,4)P3 etc.) - via de sinalizacao de fosfoinositideo reduzida nesse braco especifico.",
    "CWH43 aumentado = mais troca da porcao lipidica da ancora GPI de diacilglicerol para CERAMIDA. Previsao: mais proteinas ancoradas por GPI do tipo ceramida na membrana (mais ceramida localmente nesse pool lipidico especifico, nao no lipidoma total).",
    "Todos os genes SOBEM: sinal coerente de programa de diferenciacao terminal/cornificacao ATIVADO - bate com o achado real de queratinizacao aumentada em epitelio de SUI por scRNA-seq (Zhang et al. 2024) e com a hipotese ja presente no rascunho da usuaria (CALML5 como resposta compensatoria). Mecanismo efetor provavel: mais reticulacao por transglutaminase (KRT10 cross-linking + DMKN no envelope cornificado + CALML5 fisicamente associado a TGM3, ja citado no rascunho).",
    "SERPINB2 (PAI-2) aumentado = MAIS inibicao do ativador de plasminogenio tipo uroquinase (uPA) -> MENOS geracao de plasmina -> fibrinolise/proteolise de matriz REDUZIDA. Tambem tem papel anti-apoptotico documentado. Sugere um eixo adicional (nao destacado no rascunho atual): possivel REDUCAO do turnover de matriz extracelular, complementar a discussao de MMP/TIMP ja presente na literatura de POP.",
    "NAO recomendado citar individualmente - sao vias muito amplas e nao especificas de assoalho pelvico, cada uma \"significativa\" so por causa de 1 gene generico (CALML5 pertence a familia calmodulina, anotada em dezenas de vias de calcio de tecidos nao relacionados; mesma logica para KRT10/SERPINB2). Citar isso individualmente seria o tipo de 'forcar achado' que a usuaria quer evitar."
  ),
  Tipo_de_evidencia = c(
    "Inferencia mecanistica (bioquimica da enzima + direcao real do log2FC)",
    "Inferencia mecanistica (idem)",
    "Inferencia mecanistica (idem)",
    "Inferencia mecanistica + confirmacao empirica externa (scRNA-seq Zhang 2024)",
    "Inferencia mecanistica (bioquimica classica do PAI-2)",
    "Artefato estatistico de lista pequena - nao e achado biologico especifico"
  ),
  stringsAsFactors = FALSE
)

write.csv(pathway_deep_dive, "results/pathway_deep_dive_6genes.csv", row.names = FALSE)
cat("\nTabela salva em results/pathway_deep_dive_6genes.csv\n")
print(pathway_deep_dive[, c("Via", "Gene_motor")])
