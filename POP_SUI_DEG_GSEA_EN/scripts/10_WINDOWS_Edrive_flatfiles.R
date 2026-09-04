# =============================================================================
# ANÁLISE 12x12: POP (GSE208261, RNA-seq) vs SUI (Wei 2020, microarranjo)
# Script autocontido para rodar no SEU computador (Windows), do zero.
# =============================================================================
#
# COMO USAR (leia isto antes de rodar):
#
#   1. A pasta E:\POP+SUI 63 já existe e já tem 3 dos 4 arquivos necessários
#      (você confirmou com list.files()). Ficam TODOS soltos direto na
#      pasta, sem subpasta "data" desta vez - versão deste script ajustada
#      exatamente para isso.
#   2. Confirme que estes 4 arquivos estão em E:\POP+SUI 63 (nomes EXATOS):
#        - GSE208261_raw_counts_GRCh38.p13_NCBI.tsv  (você já tem)
#        - 43032_2020_144_MOESM2_ESM.xls              (você já tem - o
#          correto é o MOESM2, não o MOESM1, que é uma tabela diferente)
#        - GSE208261_sample_metadata.csv              (reenviado - falta
#          adicionar, sem ele não dá pra saber quem é Controle/POP)
#        - kegg_pathway_names.csv                     (reenviado - falta
#          adicionar, só usado para nomear as vias nos resultados/gráficos)
#   3. Abra o RStudio, abra este arquivo (.R), e rode do início ao fim
#      (Ctrl+Alt+R roda o script inteiro, ou clique em "Source").
#   4. Na PRIMEIRA vez que rodar, a Parte 0 vai instalar os pacotes que
#      faltarem no seu computador - isso pode demorar 5-15 minutos e só
#      acontece uma vez. Da segunda vez em diante já roda direto.
#   5. Resultados vão aparecer em E:\POP+SUI 63\results\ e as figuras em
#      E:\POP+SUI 63\figures\ (o script cria essas pastas sozinho).
#
# O QUE ESTE SCRIPT FAZ (resumo): pega os 24 dados de RNA-seq de GSE208261
# (12 mulheres com POP, 12 controles) e compara com os dados de microarranjo
# de Wei et al. 2020 (SUI) - a MESMA análise 12x12 já documentada no projeto
# (DEG via DESeq2, GSEA nas duas doenças, comparação de vias compartilhadas).
# Se os números que aparecerem no final baterem com os já reportados
# (163 DEG, 117/218 vias significativas em POP, 22 vias compartilhadas,
# 83,3% de concordância direcional), isso VALIDA de forma independente que
# a análise está correta - é exatamente para isso que este script existe.
#
# POR QUE UM SCRIPT NOVO EM VEZ DO scripts/05_...R QUE JÁ EXISTIA: aquele
# script foi escrito e testado num ambiente Linux onde os pacotes R já
# estavam pré-instalados via apt - ele nunca verificava/instalava pacotes
# sozinho. Rodado num Windows limpo, a primeira linha (library(DESeq2))
# já dava erro "there is no package called 'DESeq2'" e o script parava ALI,
# antes até de carregar qualquer dado - por isso parecia que "não carregava
# os datasets": na verdade nunca chegava a esse ponto. Este script corrige
# isso com a Parte 0 abaixo, e explica o "porquê" de cada parte, não só
# o "o quê".
# =============================================================================


## =============================================================================
## PARTE 0: pasta de trabalho + instalação automática de pacotes
## =============================================================================
# setwd() diz ao R "a partir de agora, procure os arquivos dentro desta
# pasta". Ajuste o caminho abaixo SE a sua pasta não for exatamente esta.
setwd("E:/POP+SUI 63")
cat("Pasta de trabalho definida como:", getwd(), "\n\n")

# Cada pacote abaixo faz uma coisa específica na análise:
#   - limma, edgeR: normalização e modelos estatísticos para RNA-seq/microarranjo
#   - DESeq2: o método principal para achar genes diferencialmente expressos
#     (DEG) no RNA-seq de POP
#   - org.Hs.eg.db + AnnotationDbi: dicionário que traduz IDs de genes
#     (Entrez <-> símbolo do gene, ex. "7157" <-> "TP53") e dá a lista de
#     genes de cada via do KEGG
#   - readxl: para ler o arquivo .xls do Wei 2020
#   - ggplot2, ggrepel, pheatmap: só para desenhar os gráficos no final
#
# BiocManager é o "instalador" usado para pacotes de Bioconductor (os 4
# primeiros da lista); install.packages() é o instalador padrão do R, usado
# para os demais. O bloco abaixo verifica um por um: se já está instalado,
# pula; se não está, instala - e avisa o que está fazendo, para você
# acompanhar o progresso (isso pode demorar bastante na primeira vez).

cat("=== Verificando pacotes necessários (só instala o que faltar) ===\n\n")

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  cat("Instalando o BiocManager (instalador de pacotes de Bioconductor)...\n")
  install.packages("BiocManager")
}

pacotes_bioc <- c("limma", "edgeR", "DESeq2", "org.Hs.eg.db", "AnnotationDbi")
for (pkg in pacotes_bioc) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("Instalando (Bioconductor):", pkg, "... (pode demorar alguns minutos)\n")
    BiocManager::install(pkg, update = FALSE, ask = FALSE)
  } else {
    cat("OK, já instalado:", pkg, "\n")
  }
}

pacotes_cran <- c("readxl", "ggplot2", "ggrepel", "pheatmap")
for (pkg in pacotes_cran) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("Instalando (CRAN):", pkg, "...\n")
    install.packages(pkg)
  } else {
    cat("OK, já instalado:", pkg, "\n")
  }
}

cat("\n=== Carregando os pacotes na sessão ===\n\n")
suppressMessages({
  library(limma)
  library(edgeR)
  library(DESeq2)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(ggplot2)
  library(pheatmap)
  library(readxl)
  library(ggrepel)
})
cat("Todos os pacotes carregados com sucesso.\n\n")

# select() existe em mais de um pacote carregado (AnnotationDbi e, se você
# tiver o tidyverse instalado, dplyr também) - forçamos explicitamente qual
# usar, para nunca dar o erro silencioso "objeto de tipo errado" por causa
# de um select() errado ter sido chamado.
select <- AnnotationDbi::select

dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)
cat("Pastas 'results' e 'figures' prontas dentro de", getwd(), "\n\n")


## =============================================================================
## PARTE 0b: conferir se os 4 arquivos de dados estão no lugar certo
## =============================================================================
# Em vez de deixar o script quebrar mais adiante com um erro genérico tipo
# "cannot open file" (que não diz QUAL arquivo nem ONDE deveria estar), a
# gente confere os 4 arquivos JÁ AQUI, no início, e para com uma mensagem
# clara em português se algum estiver faltando.

arquivos_necessarios <- c(
  "GSE208261_raw_counts_GRCh38.p13_NCBI.tsv",
  "GSE208261_sample_metadata.csv",
  "43032_2020_144_MOESM2_ESM.xls",
  "kegg_pathway_names.csv"
)
faltando <- arquivos_necessarios[!file.exists(arquivos_necessarios)]
if (length(faltando) > 0) {
  stop(
    "\n\nFALTAM ARQUIVOS DE DADOS. O script parou aqui de propósito, ANTES\n",
    "de tentar rodar qualquer análise, para não gerar erros confusos depois.\n\n",
    "Pasta de trabalho atual: ", getwd(), "\n",
    "Arquivo(s) que não foram encontrados:\n  - ", paste(faltando, collapse = "\n  - "), "\n\n",
    "Solução: copie esse(s) arquivo(s) para dentro de ", getwd(), "\n",
    "com exatamente esse nome, e rode o script de novo.\n"
  )
}
cat("=== Os 4 arquivos de dados foram encontrados. Prosseguindo. ===\n\n")


## Funções auxiliares usadas mais abaixo -------------------------------------
# kegg_label(): transforma um número de via do KEGG (ex. "04510") no nome
# legível ("KEGG 04510 - Focal adhesion"), usando a tabela que vem no arquivo
# kegg_pathway_names.csv (sem isso, os gráficos mostrariam só números).
kegg_tab <- read.csv("kegg_pathway_names.csv", colClasses = c("character", "character"))
kegg_names <- setNames(kegg_tab$Name, kegg_tab$PATH5)
kegg_label <- function(id) {
  nm <- kegg_names[id]
  ifelse(is.na(nm), paste0("KEGG ", id), paste0("KEGG ", id, " - ", nm))
}

# calc_es(): calcula o "Enrichment Score" (ES) do GSEA - a estatística que
# mede se os genes de uma via aparecem mais no topo (ou na base) da lista
# ranqueada de genes do que se estivessem espalhados ao acaso. Implementado
# manualmente aqui (fórmula de Subramanian et al. 2005) porque o pacote
# oficial (fgsea) não está disponível sem acesso à internet no ambiente
# onde este projeto foi originalmente construído - rodando no seu PC com
# internet, o resultado é matematicamente equivalente.
calc_es <- function(hit_idx, scores_abs, N) {
  Nh <- length(hit_idx); Nm <- N - Nh
  if (Nm <= 0 || Nh == 0) return(NA)
  sum_hit <- sum(scores_abs[hit_idx])
  step <- rep(-1 / Nm, N); step[hit_idx] <- scores_abs[hit_idx] / sum_hit
  running <- cumsum(step); running[which.max(abs(running))]
}


## =============================================================================
## PARTE 1: carregar os dados de POP (GSE208261) e montar o desenho 12x12
## =============================================================================
cat("\n================ PARTE 1: dados de POP (GSE208261) ================\n\n")

# read.csv: lê a tabela de metadados - qual amostra (GSM) é Controle ou POP,
# e de qual tecido (parede vaginal ou ligamento uterossacral) ela veio.
meta <- read.csv("GSE208261_sample_metadata.csv")
cat("Tabela de amostras (confira contra os metadados do GEO antes de confiar cegamente):\n")
print(meta)
cat("\nCruzamento Tecido x Grupo:\n")
print(table(meta$Tissue, meta$Treatment))
cat("\n")

# read.delim: lê a matriz de contagens brutas (quantas "leituras" de
# sequenciamento caíram em cada gene, em cada amostra) - a matéria-prima
# de qualquer análise de RNA-seq.
counts_raw <- read.delim("GSE208261_raw_counts_GRCh38.p13_NCBI.tsv", row.names = 1, check.names = FALSE)
counts_all <- as.matrix(counts_raw[, meta$GSM])
stopifnot(identical(colnames(counts_all), meta$GSM))
cat("Matriz de contagens: 24 amostras,", nrow(counts_all), "genes (por ID Entrez)\n")

# Os genes vêm identificados por número (ID Entrez, ex. "7157"), não por
# nome (símbolo, ex. "TP53") - convertendo aqui para os nomes ficarem
# legíveis no resto do script e comparáveis com o Wei2020 (que já usa nomes).
ann_id <- suppressWarnings(select(org.Hs.eg.db, keys = rownames(counts_all), keytype = "ENTREZID", columns = "SYMBOL"))
ann_id <- ann_id[!is.na(ann_id$SYMBOL) & !duplicated(ann_id$ENTREZID), ]
counts_all <- counts_all[rownames(counts_all) %in% ann_id$ENTREZID, ]
rownames(counts_all) <- ann_id$SYMBOL[match(rownames(counts_all), ann_id$ENTREZID)]
counts_all <- rowsum(counts_all, group = rownames(counts_all))  # soma genes duplicados
cat("Depois de converter ID->símbolo e somar duplicados:", nrow(counts_all), "genes\n\n")

# O desenho 12x12: 12 Controle (6 ligamento + 6 parede vaginal, JUNTOS) vs
# 12 POP (parede vaginal). Ver o relatório para a discussão completa sobre
# por que os controles têm dois tecidos e o que isso implica.
group_full <- factor(meta$Treatment, levels = c("Control", "POP"))
tissue_full <- factor(meta$Tissue)
cat("Desenho: 12 Controle (6 ligamento + 6 parede vaginal) vs 12 POP (parede vaginal). Controle:",
    sum(group_full == "Control"), "| POP:", sum(group_full == "POP"), "\n\n")


## =============================================================================
## PARTE 2: DEG em POP via DESeq2
## =============================================================================
cat("================ PARTE 2: genes diferencialmente expressos (DEG) em POP ================\n\n")

counts_int <- counts_all
storage.mode(counts_int) <- "integer"  # DESeq2 exige contagens inteiras

coldata_naive <- data.frame(row.names = meta$GSM, group = group_full)
dds_naive <- DESeqDataSetFromMatrix(countData = counts_int, colData = coldata_naive, design = ~group)

# FILTRO DE BAIXA CONTAGEM - este é o passo que você perguntou se tinha sido
# esquecido em outra análise. Aqui ele está: mantemos só genes com pelo
# menos 10 leituras em pelo menos 12 das 24 amostras (o tamanho do menor
# grupo). Genes com contagem muito baixa não têm poder estatístico e, se
# deixados na tabela, só atrapalham a correção de múltiplos testes (FDR) -
# é a prática recomendada pela própria documentação do DESeq2.
dds_naive <- dds_naive[rowSums(counts(dds_naive) >= 10) >= 12, ]
cat("Genes mantidos após o filtro de baixa contagem:", nrow(dds_naive), "\n")

dds_naive <- DESeq(dds_naive, quiet = TRUE)
res_naive <- as.data.frame(results(dds_naive, contrast = c("group", "POP", "Control"), alpha = 0.05))
res_naive$Gene <- rownames(res_naive)
res_naive <- res_naive[order(res_naive$pvalue), ]
pop_full <- res_naive[, c("Gene", "log2FoldChange", "baseMean", "stat", "pvalue", "padj")]
colnames(pop_full) <- c("Gene", "logFC", "baseMean", "stat", "P.Value", "adj.P.Val")
write.csv(pop_full, "results/POP_DESeq2_full_table.csv", row.names = FALSE)

# DEG = genes com mudança de pelo menos 2x (|log2FC|>1) E estatisticamente
# significativos depois de corrigir para múltiplos testes (FDR<0,05).
pop_deg <- subset(pop_full, !is.na(adj.P.Val) & adj.P.Val < 0.05 & abs(logFC) > 1)
pop_deg <- pop_deg[order(pop_deg$adj.P.Val), ]
write.csv(pop_deg, "results/POP_DEG_logFC1_FDR05.csv", row.names = FALSE)

cat("\n=== RESULTADO: DEG em POP (12 Controle vs 12 POP) ===\n")
cat("Genes testados:", nrow(pop_full), "\n")
cat("DEG (|log2FC|>1, FDR<0,05):", nrow(pop_deg), "(",
    sum(pop_deg$logFC > 0), "para cima /", sum(pop_deg$logFC < 0), "para baixo )\n")
cat("Valor esperado, já documentado no projeto: 163 DEG (123 up / 40 down).\n")
cat("Se o número acima bater com 163, esta parte está validada.\n\n")

# edgeR+voom+limma: um SEGUNDO método estatístico, independente do DESeq2,
# rodado aqui só para conferência e para gerar o ranking de genes usado no
# GSEA (Parte 3). Não é o método principal de DEG deste script - é um
# cross-check, e serve para mostrar que os dois métodos concordam ou
# discordam de forma esperada (ver discussão no relatório sobre DESeq2 vs
# edgeR).
dge <- DGEList(counts = counts_all, group = group_full)
design_naive_lm <- model.matrix(~group_full)
keep_expr <- filterByExpr(dge, design_naive_lm)  # o filtro equivalente do pacote edgeR
dge <- dge[keep_expr, , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge, method = "TMM")
voom_fit <- voom(dge, design_naive_lm)
fit_pop <- eBayes(lmFit(voom_fit, design_naive_lm))
pop_full_limma <- topTable(fit_pop, coef = "group_fullPOP", number = Inf, sort.by = "P")
pop_full_limma$Gene <- rownames(pop_full_limma)
write.csv(pop_full_limma[, c("Gene","logFC","AveExpr","t","P.Value","adj.P.Val")],
          "results/POP_voom_limma_full_table.csv", row.names = FALSE)
pop_deg_limma <- subset(pop_full_limma, adj.P.Val < 0.05 & abs(logFC) > 1)
cat("Conferência com edgeR+voom+limma (método secundário):", nrow(pop_deg_limma), "DEG\n\n")


## =============================================================================
## PARTE 3: GSEA em POP (método clássico)
## =============================================================================
cat("================ PARTE 3: GSEA em POP (método clássico) ================\n\n")
cat("GSEA pergunta: mesmo sem nenhum gene individual 'dar significativo', será\n")
cat("que um GRUPO inteiro de genes de uma mesma via biológica está, em conjunto,\n")
cat("deslocado na mesma direção? Isso detecta sinal fraco e distribuído que um\n")
cat("teste gene-a-gene sozinho perderia.\n\n")
cat("Método CLÁSSICO (permutação de fenótipo): válido aqui porque temos as 24\n")
cat("amostras completas e um bom número de reordenações possíveis\n")
cat("(escolher 12 de 24 = 2.704.156 combinações) - resolução estatística\n")
cat("excelente.\n\n")

set.seed(208261)  # fixa a semente aleatória: rodar de novo dá o MESMO resultado
ranked_full <- fit_pop$t[, "group_fullPOP"]
ord2 <- order(-ranked_full)
ranked_genes_pop <- names(ranked_full)[ord2]
ranked_scores_pop <- ranked_full[ord2]
N_pop <- length(ranked_genes_pop)

# Para cada gene ranqueado, descobrir a quais vias do KEGG ele pertence.
ann_pop <- suppressWarnings(select(org.Hs.eg.db, keys = ranked_genes_pop, keytype = "SYMBOL", columns = "PATH"))
ann_pop <- ann_pop[!is.na(ann_pop$PATH), ]
gs_sizes_pop <- table(ann_pop$PATH)
valid_paths_pop <- names(gs_sizes_pop)[gs_sizes_pop >= 5 & gs_sizes_pop <= 200]  # vias nem pequenas nem enormes demais
gene_sets_pop <- split(ann_pop$SYMBOL[ann_pop$PATH %in% valid_paths_pop], ann_pop$PATH[ann_pop$PATH %in% valid_paths_pop])
cat("Vias do KEGG testadas (entre 5 e 200 genes):", length(gene_sets_pop), "\n")

hit_idx_pop <- lapply(gene_sets_pop, function(g) which(ranked_genes_pop %in% g))
hit_idx_pop <- hit_idx_pop[sapply(hit_idx_pop, length) >= 3]
cat("Vias com pelo menos 3 genes na lista ranqueada:", length(hit_idx_pop), "\n")
es_obs_pop <- sapply(hit_idx_pop, calc_es, scores_abs = abs(ranked_scores_pop), N = N_pop)

# A permutação: embaralhar os rótulos Controle/POP 500 vezes, recalcular o
# ES em cada embaralhamento, e comparar o ES real contra essa distribuição
# "ao acaso" - isso é o que dá o p-valor de cada via.
n_perm <- 500
cat("Rodando", n_perm, "permutações (embaralhamentos) do rótulo Controle/POP...\n")
cat("(isso é a parte mais demorada do script - pode levar alguns minutos)\n")
t0 <- Sys.time()
gene_sets_syms_pop <- lapply(hit_idx_pop, function(idx) ranked_genes_pop[idx])
perm_es_pop <- matrix(NA_real_, nrow = n_perm, ncol = length(hit_idx_pop))
for (i in seq_len(n_perm)) {
  perm_group <- sample(group_full)
  perm_design <- model.matrix(~perm_group)
  perm_fit <- eBayes(lmFit(voom_fit, perm_design))
  t_perm <- perm_fit$t[, 2]
  rank_of_gene <- rank(-t_perm, ties.method = "first")
  scores_abs_sorted <- sort(abs(t_perm), decreasing = TRUE)
  for (j in seq_along(gene_sets_syms_pop)) {
    hidx <- rank_of_gene[gene_sets_syms_pop[[j]]]
    if (length(hidx) >= 3) perm_es_pop[i, j] <- calc_es(hidx, scores_abs_sorted, N_pop)
  }
}
cat("Concluído em", round(difftime(Sys.time(), t0, units = "secs"), 1), "segundos\n\n")

pval_pop <- numeric(length(hit_idx_pop)); nes_pop <- numeric(length(hit_idx_pop))
for (j in seq_along(hit_idx_pop)) {
  pe <- perm_es_pop[, j]; pe <- pe[!is.na(pe)]
  if (es_obs_pop[j] >= 0) {
    pval_pop[j] <- (sum(pe >= es_obs_pop[j]) + 1) / (length(pe) + 1)
    base <- mean(pe[pe >= 0]); if (is.nan(base) || base == 0) base <- NA
  } else {
    pval_pop[j] <- (sum(pe <= es_obs_pop[j]) + 1) / (length(pe) + 1)
    base <- mean(abs(pe[pe < 0])); if (is.nan(base) || base == 0) base <- NA
  }
  nes_pop[j] <- es_obs_pop[j] / base
}

leading_edge_pop <- sapply(hit_idx_pop, function(idx) paste(ranked_genes_pop[idx], collapse = "/"))
gsea_pop <- data.frame(PATH = names(hit_idx_pop), Nh = sapply(hit_idx_pop, length),
                        ES = es_obs_pop, NES = nes_pop, pvalue = pval_pop, leadingEdge = leading_edge_pop)
gsea_pop$p.adjust <- p.adjust(gsea_pop$pvalue, "BH")
gsea_pop$PathwayName <- kegg_label(gsea_pop$PATH)
gsea_pop <- gsea_pop[order(gsea_pop$pvalue), ]
write.csv(gsea_pop, "results/POP_GSEA_classic_KEGG.csv", row.names = FALSE)

cat("=== RESULTADO: GSEA clássico em POP ===\n")
cat("Vias significativas a FDR<0,25:", sum(gsea_pop$p.adjust < 0.25, na.rm = TRUE), "de", nrow(gsea_pop), "\n")
cat("Vias significativas a FDR<0,05:", sum(gsea_pop$p.adjust < 0.05, na.rm = TRUE), "\n")
cat("Valor esperado, já documentado: 117 a FDR<0,25 e 97 a FDR<0,05, de 218 testadas.\n\n")


## =============================================================================
## PARTE 4: dados de SUI (Wei 2020) - DEG (já vem pronto) + GSEA pré-ranqueado
## =============================================================================
cat("================ PARTE 4: dados de SUI (Wei 2020) ================\n\n")

wei_xls <- "43032_2020_144_MOESM2_ESM.xls"
sheet_names <- excel_sheets(wei_xls)
cat("Abas encontradas no Excel:", paste(sheet_names, collapse = ", "), "\n")
stopifnot(all(c("up_Sui_vs_Ctrl", "down_Sui_vs_Ctrl") %in% sheet_names))

# skip=17: as primeiras 17 linhas do Excel são um cabeçalho descritivo dos
# autores (não são dados) - pulamos elas para chegar na tabela de verdade.
up_sheet   <- read_excel(wei_xls, sheet = "up_Sui_vs_Ctrl",   skip = 17)
down_sheet <- read_excel(wei_xls, sheet = "down_Sui_vs_Ctrl", skip = 17)
cat("'up_Sui_vs_Ctrl':", nrow(up_sheet), "linhas,", ncol(up_sheet), "colunas\n")
cat("'down_Sui_vs_Ctrl':", nrow(down_sheet), "linhas,", ncol(down_sheet), "colunas\n\n")

# PROVA DIRETA DE QUE SÃO 3 AMOSTRAS SUI E 3 CONTROLE (não é uma suposição -
# são os nomes reais das colunas do arquivo, impressos abaixo para você
# conferir com os próprios olhos):
cat("Colunas de amostra individual encontradas no arquivo do Wei 2020:\n")
print(grep("Sui[0-9]|Ctrl[0-9]", colnames(up_sheet), value = TRUE))
cat("-> 3 colunas 'SuiN' + 3 colunas 'CtrlN' = desenho 3 vs 3, confirmado\n")
cat("   diretamente pelos nomes de coluna do arquivo, não por suposição.\n\n")

# Conferindo de novo (é a segunda vez que esse teste é feito no projeto,
# de propósito) a convenção de sinal do "Fold Change": a aba "down" lista
# a MAGNITUDE (sempre >=1), não um valor negativo - o sinal vem da aba em
# que o gene está, não do número em si.
cat("Conferindo a convenção de sinal do Fold Change:\n")
cat("  aba 'up', faixa de Fold Change:", round(min(up_sheet$`Fold Change`), 2), "a",
    round(max(up_sheet$`Fold Change`), 2), "\n")
cat("  aba 'down', faixa de Fold Change:", round(min(down_sheet$`Fold Change`), 2), "a",
    round(max(down_sheet$`Fold Change`), 2), "\n")
if (min(down_sheet$`Fold Change`) >= 1) {
  cat("  CONFIRMADO: a aba 'down' traz só a MAGNITUDE (sempre >=1) - o sinal\n")
  cat("  negativo precisa ser aplicado manualmente com base na aba/Direction,\n")
  cat("  não vem pronto no número. O código abaixo já faz essa correção.\n\n")
}

up_sheet$Direction <- "up"; down_sheet$Direction <- "down"
sample_cols <- c("[Sui1, Sui](normalized)", "[Sui2, Sui](normalized)", "[Sui3, Sui](normalized)",
                  "[Ctrl1, Ctrl](normalized)", "[Ctrl2, Ctrl](normalized)", "[Ctrl3, Ctrl](normalized)")
keep_cols <- c("GeneSymbol", "P-value", "FDR", "Fold Change", "Direction", sample_cols)
wei_both <- rbind(as.data.frame(up_sheet[, keep_cols]), as.data.frame(down_sheet[, keep_cols]))
colnames(wei_both) <- c("GeneSymbol", "PValue", "FDR", "FoldChange", "Direction",
                         "Sui1", "Sui2", "Sui3", "Ctrl1", "Ctrl2", "Ctrl3")
wei_both <- wei_both[!is.na(wei_both$GeneSymbol), ]
cat("Total combinado (as duas abas, antes de remover duplicatas):", nrow(wei_both), "sondas\n")

# Aqui é onde a correção de sinal acontece de fato: genes 'down' recebem
# log2(FoldChange) NEGATIVO; genes 'up' recebem log2(FoldChange) positivo.
wei_both$logFC <- ifelse(wei_both$Direction == "down",
                          -log2(wei_both$FoldChange), log2(wei_both$FoldChange))
wei_both <- wei_both[order(wei_both$PValue), ]
sui_full <- wei_both[!duplicated(wei_both$GeneSymbol), ]  # mantém o probe mais significativo por gene
rownames(sui_full) <- NULL
write.csv(sui_full, "results/SUI_Wei2020_full_table.csv", row.names = FALSE)
cat("Genes únicos após remover duplicatas:", nrow(sui_full), "(",
    sum(sui_full$Direction == "up"), "para cima /", sum(sui_full$Direction == "down"), "para baixo )\n")
cat("Esta tabela JÁ É a lista de DEG de SUI - o arquivo original do Wei 2020\n")
cat("só lista os genes que os autores classificaram como diferencialmente\n")
cat("expressos (fold-change>=2, p bruto<0,05); o array completo testado não\n")
cat("foi publicado. Por isso não recalculamos DEG de SUI do zero - usamos a\n")
cat("lista dos próprios autores.\n\n")

set.seed(2020)
sui_mat <- as.matrix(sui_full[, c("Sui1","Sui2","Sui3","Ctrl1","Ctrl2","Ctrl3")])
rownames(sui_mat) <- sui_full$GeneSymbol
sui_mat <- avereps(sui_mat, ID = rownames(sui_mat))
group_sui <- factor(c("SUI","SUI","SUI","Ctrl","Ctrl","Ctrl"), levels = c("Ctrl","SUI"))
design_sui <- model.matrix(~group_sui)
fit_sui <- eBayes(lmFit(sui_mat, design_sui))
t_obs_sui <- fit_sui$t[, 2]
ord_s <- order(-t_obs_sui)
ranked_genes_sui <- names(t_obs_sui)[ord_s]
ranked_scores_sui <- t_obs_sui[ord_s]
N_sui <- length(ranked_genes_sui)
cat("Lista ranqueada de SUI (estatística t moderada; positivo = para cima em SUI):", N_sui, "genes\n\n")

cat("Método PRÉ-RANQUEADO para o GSEA de SUI (diferente do método usado em\n")
cat("POP): com só 3 vs 3 amostras, a permutação de FENÓTIPO teria apenas\n")
cat("choose(6,3)=20 reordenações possíveis, limitando a resolução a\n")
cat("p=1/20=0,05 - não dá para chegar em FDR<0,05 nunca, não importa quão\n")
cat("real seja o sinal. A permutação PRÉ-RANQUEADA (embaralha as ~200 vias\n")
cat("testadas, não as 6 amostras) não tem esse teto.\n\n")

ann_sui <- suppressWarnings(select(org.Hs.eg.db, keys = ranked_genes_sui, keytype = "SYMBOL", columns = "PATH"))
ann_sui <- ann_sui[!is.na(ann_sui$PATH), ]
gs_sizes_sui <- table(ann_sui$PATH)
valid_paths_sui <- names(gs_sizes_sui)[gs_sizes_sui >= 5 & gs_sizes_sui <= 200]
gene_sets_sui <- split(ann_sui$SYMBOL[ann_sui$PATH %in% valid_paths_sui], ann_sui$PATH[ann_sui$PATH %in% valid_paths_sui])
hit_idx_sui <- lapply(gene_sets_sui, function(g) which(ranked_genes_sui %in% g))
hit_idx_sui <- hit_idx_sui[sapply(hit_idx_sui, length) >= 3]
es_obs_sui <- sapply(hit_idx_sui, calc_es, scores_abs = abs(ranked_scores_sui), N = N_sui)

n_perm2 <- 1000
cat("Rodando", n_perm2, "permutações (embaralhamento das vias) para SUI...\n")
t0 <- Sys.time()
scores_abs_sui <- abs(ranked_scores_sui)
perm_es_sui <- matrix(NA_real_, nrow = n_perm2, ncol = length(hit_idx_sui))
for (i in seq_len(n_perm2)) {
  for (j in seq_along(hit_idx_sui)) {
    hidx <- sample.int(N_sui, length(hit_idx_sui[[j]]))
    perm_es_sui[i, j] <- calc_es(hidx, scores_abs_sui, N_sui)
  }
}
cat("Concluído em", round(difftime(Sys.time(), t0, units = "secs"), 1), "segundos\n\n")

pval_sui <- numeric(length(hit_idx_sui)); nes_sui <- numeric(length(hit_idx_sui))
for (j in seq_along(hit_idx_sui)) {
  pe <- perm_es_sui[, j]; pe <- pe[!is.na(pe)]
  if (es_obs_sui[j] >= 0) {
    pval_sui[j] <- (sum(pe >= es_obs_sui[j]) + 1) / (length(pe) + 1)
    base <- mean(pe[pe >= 0]); if (is.nan(base) || base == 0) base <- NA
  } else {
    pval_sui[j] <- (sum(pe <= es_obs_sui[j]) + 1) / (length(pe) + 1)
    base <- mean(abs(pe[pe < 0])); if (is.nan(base) || base == 0) base <- NA
  }
  nes_sui[j] <- es_obs_sui[j] / base
}

leading_edge_sui <- sapply(hit_idx_sui, function(idx) paste(ranked_genes_sui[idx], collapse = "/"))
gsea_sui <- data.frame(PATH = names(hit_idx_sui), Nh = sapply(hit_idx_sui, length),
                        ES = es_obs_sui, NES = nes_sui, pvalue = pval_sui, leadingEdge = leading_edge_sui)
gsea_sui$p.adjust <- p.adjust(gsea_sui$pvalue, "BH")
gsea_sui$PathwayName <- kegg_label(gsea_sui$PATH)
gsea_sui <- gsea_sui[order(gsea_sui$pvalue), ]
write.csv(gsea_sui, "results/SUI_GSEA_preranked_KEGG.csv", row.names = FALSE)

cat("=== RESULTADO: GSEA pré-ranqueado em SUI ===\n")
cat("Vias significativas a FDR<0,25:", sum(gsea_sui$p.adjust < 0.25, na.rm = TRUE), "de", nrow(gsea_sui), "\n")
cat("Vias significativas a FDR<0,05:", sum(gsea_sui$p.adjust < 0.05, na.rm = TRUE), "\n\n")


## =============================================================================
## PARTE 5: vias compartilhadas entre POP e SUI + tabela de genes
## =============================================================================
cat("================ PARTE 5: vias compartilhadas entre POP e SUI ================\n\n")

reportar_compartilhadas <- function(corte_fdr) {
  pop_sig <- subset(gsea_pop, p.adjust < corte_fdr)
  sui_sig <- subset(gsea_sui, p.adjust < corte_fdr)
  ids_compartilhados <- intersect(pop_sig$PATH, sui_sig$PATH)
  cat("--- FDR <", corte_fdr, ": POP significativas =", nrow(pop_sig),
      "| SUI significativas =", nrow(sui_sig), "| COMPARTILHADAS =", length(ids_compartilhados), "---\n")
  if (length(ids_compartilhados) == 0) return(data.frame())
  out <- merge(pop_sig[pop_sig$PATH %in% ids_compartilhados, c("PATH","PathwayName","Nh","NES","p.adjust")],
               sui_sig[sui_sig$PATH %in% ids_compartilhados, c("PATH","Nh","NES","p.adjust")],
               by = "PATH", suffixes = c("_POP", "_SUI"))
  out$Same_direction <- sign(out$NES_POP) == sign(out$NES_SUI)
  out <- out[order(out$p.adjust_POP), ]
  n_na <- sum(is.na(out$Same_direction))
  cat("Mesma direção nas duas doenças:", sum(out$Same_direction, na.rm = TRUE), "de",
      sum(!is.na(out$Same_direction)), "vias com NES comparável")
  if (n_na > 0) cat(" (", n_na, "excluídas - NES indefinido em algum lado)")
  cat("\n\n")
  out
}

compartilhadas_025 <- reportar_compartilhadas(0.25)
if (nrow(compartilhadas_025) > 0) {
  write.csv(compartilhadas_025, "results/shared_pathways_FDR025.csv", row.names = FALSE)
}
cat("Valor esperado, já documentado: 22 vias compartilhadas a FDR<0,25,\n")
cat("das quais 10 de 12 comparáveis (83,3%) concordantes em direção.\n\n")


## =============================================================================
## PARTE 6: gráfico rápido de conferência (volcano plot de POP)
## =============================================================================
cat("================ PARTE 6: gráfico de conferência ================\n\n")

pop_full$sig <- "NS"
pop_full$sig[pop_full$logFC > 1 & pop_full$adj.P.Val < 0.05] <- "Up"
pop_full$sig[pop_full$logFC < -1 & pop_full$adj.P.Val < 0.05] <- "Down"
pop_full$sig <- factor(pop_full$sig, levels = c("Down", "NS", "Up"))
p_volcano <- ggplot(pop_full, aes(x = logFC, y = -log10(P.Value), color = sig)) +
  geom_point(alpha = 0.6, size = 1.2) +
  scale_color_manual(values = c(Down = "#2166AC", NS = "grey75", Up = "#B2182B")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey40") +
  labs(title = paste0("Volcano - POP (GSE208261, 12x12), ", nrow(pop_deg), " DEG a FDR<0,05"),
       x = "log2(Fold Change)", y = "-log10(p-valor)", color = NULL) +
  theme_bw() + theme(legend.position = "top")
ggsave("figures/volcano_POP.png", p_volcano, width = 8, height = 6, dpi = 300)
cat("Salvo: figures/volcano_POP.png\n\n")

cat("=================================================================\n")
cat("=== FIM DO SCRIPT ===\n")
cat("Confira os números marcados 'Valor esperado, já documentado' acima\n")
cat("contra o que saiu aqui no seu computador. Se baterem, a análise está\n")
cat("validada de forma independente. Resultados completos em results/,\n")
cat("figura de conferência em figures/volcano_POP.png\n")
cat("=================================================================\n")
