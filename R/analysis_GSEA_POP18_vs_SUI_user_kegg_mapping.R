# =============================================================================
# GSEA (KEGG) de POP e de SUI, usando:
#   - POP (GSE208261): SO as 18 amostras do MESMO TECIDO (parede vaginal
#     anterior) - 12 POP + 6 Controle_Y, excluindo os 6 Controle_D
#     (ligamento uterossacral, tecido diferente - ver secao de correcao no
#     README sobre D/Y = tecido, nao idade). GSEA REAL (permutacao de
#     fenotipo), identica ao metodo ja estabelecido em
#     R/analysis_GSEA_POP_vs_SUI36h_real_RNAseq.R, so que agora com 18
#     amostras tissue-matched em vez das 24 originais (que misturavam
#     tecido no lado Controle).
#   - SUI (tabela 43032_2020_144, fornecida pela usuaria): GSEA pequeno-
#     painel, mesmo metodo ja usado antes (permutacao de rotulo de gene
#     set - essa tabela so tem os genes ja filtrados pelos autores, nao
#     tem matriz bruta completa).
#   - Vias (KEGG gene sets) desta vez vem do arquivo da propria usuaria
#     (data/KEGG_gene_pathway_mapping_usuaria.csv, ENTREZID/SYMBOL/PATH),
#     nao do org.Hs.eg.db, para os dois lados usarem exatamente a mesma
#     fonte de vias.
# =============================================================================

suppressMessages(library(DESeq2))
set.seed(2024)

## -----------------------------------------------------------------------
## 0) KEGG gene sets da usuaria
## -----------------------------------------------------------------------
kegg_map <- read.csv("data/KEGG_gene_pathway_mapping_usuaria.csv")
kegg_map$ENTREZID <- as.character(kegg_map$ENTREZID)
kegg_map$PATH <- sprintf("%05d", as.integer(kegg_map$PATH))
kegg_sets <- split(kegg_map$ENTREZID, kegg_map$PATH)
cat("KEGG pathway definitions (arquivo da usuaria):", length(kegg_sets), "vias\n\n")

## -----------------------------------------------------------------------
## Funcoes reaproveitadas sem alteracao dos scripts ja estabelecidos
## -----------------------------------------------------------------------
calc_es <- function(hit_idx, scores_abs, N) {
  Nh <- length(hit_idx); Nm <- N - Nh
  if (Nm <= 0 || Nh == 0) return(NA)
  sum_hit <- sum(scores_abs[hit_idx])
  step <- rep(-1 / Nm, N)
  step[hit_idx] <- scores_abs[hit_idx] / sum_hit
  running <- cumsum(step)
  running[which.max(abs(running))]
}

fast_welch_t <- function(logmat, group) {
  lv <- levels(group)
  g1 <- logmat[, group == lv[1], drop = FALSE]
  g2 <- logmat[, group == lv[2], drop = FALSE]
  n1 <- ncol(g1); n2 <- ncol(g2)
  m1 <- rowMeans(g1); m2 <- rowMeans(g2)
  v1 <- rowSums((g1 - m1)^2) / (n1 - 1)
  v2 <- rowSums((g2 - m2)^2) / (n2 - 1)
  se <- sqrt(v1 / n1 + v2 / n2)
  (m2 - m1) / se
}

run_gsea_real <- function(gene_ids, observed_stat, logmat, group, gene_sets,
                           min_gs = 5, max_gs = 200, n_perm = 500) {
  N <- length(gene_ids)
  ord <- order(-observed_stat)
  ranked_ids <- gene_ids[ord]
  ranked_scores <- observed_stat[ord]
  scores_abs <- abs(ranked_scores)

  hit_idx_list <- lapply(gene_sets, function(g) which(ranked_ids %in% g))
  keep <- sapply(hit_idx_list, length) >= min_gs & sapply(hit_idx_list, length) <= max_gs
  hit_idx_list <- hit_idx_list[keep]
  cat("  Vias testaveis (", min_gs, "-", max_gs, " membros presentes):", length(hit_idx_list), "\n", sep = "")
  if (length(hit_idx_list) == 0) return(NULL)

  es_obs <- sapply(hit_idx_list, calc_es, scores_abs = scores_abs, N = N)
  gene_sets_ids <- lapply(hit_idx_list, function(idx) ranked_ids[idx])

  perm_mat <- matrix(NA_real_, nrow = n_perm, ncol = length(hit_idx_list))
  rownames_logmat <- rownames(logmat)
  for (i in seq_len(n_perm)) {
    perm_group <- sample(group)
    t_perm <- fast_welch_t(logmat, perm_group)
    names(t_perm) <- rownames_logmat
    t_perm <- t_perm[ranked_ids]
    rank_of_gene <- rank(-t_perm, ties.method = "first")
    scores_abs_sorted <- sort(abs(t_perm), decreasing = TRUE)
    for (j in seq_along(gene_sets_ids)) {
      hidx <- rank_of_gene[as.character(gene_sets_ids[[j]])]
      if (length(hidx) > 0 && !any(is.na(hidx))) {
        perm_mat[i, j] <- calc_es(hidx, scores_abs_sorted, N)
      }
    }
  }

  pval <- numeric(length(hit_idx_list)); nes <- numeric(length(hit_idx_list))
  for (j in seq_along(hit_idx_list)) {
    pe <- perm_mat[, j]; pe <- pe[!is.na(pe)]
    if (length(pe) == 0) { pval[j] <- NA; nes[j] <- NA; next }
    if (es_obs[j] >= 0) {
      pval[j] <- (sum(pe >= es_obs[j]) + 1) / (length(pe) + 1)
      base <- mean(pe[pe >= 0]); if (is.nan(base) || base == 0) base <- NA
    } else {
      pval[j] <- (sum(pe <= es_obs[j]) + 1) / (length(pe) + 1)
      base <- mean(abs(pe[pe < 0])); if (is.nan(base) || base == 0) base <- NA
    }
    nes[j] <- es_obs[j] / base
  }

  leading_edge <- sapply(hit_idx_list, function(idx) paste(ranked_ids[idx], collapse = "/"))
  res <- data.frame(PATH = names(hit_idx_list), Nh = sapply(hit_idx_list, length),
                     ES = es_obs, NES = nes, pvalue = pval, leadingEdge = leading_edge)
  res$p.adjust <- p.adjust(res$pvalue, "BH")
  res[order(res$pvalue), ]
}

run_gsea_small_panel <- function(ranked_genes, ranked_scores, gene_sets, min_gs = 2, n_perm = 1000) {
  N <- length(ranked_genes)
  scores_abs <- abs(ranked_scores)
  testable <- lapply(gene_sets, function(g) which(ranked_genes %in% g))
  testable <- testable[sapply(testable, length) >= min_gs]
  cat("  Vias testaveis (>=", min_gs, "membros presentes nesta lista de", N, "genes):", length(testable), "\n")
  if (length(testable) == 0) return(NULL)
  es_obs <- sapply(testable, calc_es, scores_abs = scores_abs, N = N)
  perm_mat <- matrix(NA_real_, nrow = n_perm, ncol = length(testable))
  for (i in seq_len(n_perm)) {
    for (j in seq_along(testable)) {
      hidx <- sample.int(N, length(testable[[j]]))
      perm_mat[i, j] <- calc_es(hidx, scores_abs, N)
    }
  }
  pval <- numeric(length(testable)); nes <- numeric(length(testable))
  for (j in seq_along(testable)) {
    pe <- perm_mat[, j]; pe <- pe[!is.na(pe)]
    if (is.na(es_obs[j])) { pval[j] <- NA; nes[j] <- NA; next }
    if (es_obs[j] >= 0) {
      pval[j] <- (sum(pe >= es_obs[j]) + 1) / (length(pe) + 1)
      base <- mean(pe[pe >= 0]); if (is.nan(base) || base == 0) base <- NA
    } else {
      pval[j] <- (sum(pe <= es_obs[j]) + 1) / (length(pe) + 1)
      base <- mean(abs(pe[pe < 0])); if (is.nan(base) || base == 0) base <- NA
    }
    nes[j] <- es_obs[j] / base
  }
  leading_edge <- sapply(testable, function(idx) paste(ranked_genes[idx], collapse = "/"))
  res <- data.frame(PATH = names(testable), Nh = sapply(testable, length),
                     ES = es_obs, NES = nes, pvalue = pval, leadingEdge = leading_edge)
  res$p.adjust <- p.adjust(res$pvalue, "BH")
  res[order(res$pvalue), ]
}

build_score <- function(direction, pvalue) {
  p <- pvalue
  floor_p <- min(p[p > 0], na.rm = TRUE) / 2
  p[p <= 0] <- floor_p
  sign <- ifelse(direction == "up", 1, -1)
  sign * -log10(p)
}

## =============================================================================
## PART 1 - POP (GSE208261), 18 amostras do mesmo tecido (parede vaginal anterior)
## =============================================================================
cat("============================================================\n")
cat("PART 1: POP (GSE208261), 18 amostras tissue-matched\n")
cat("============================================================\n\n")

pop_counts_all <- read.delim("data/GSE208261_raw_counts_GRCh38.p13_NCBI.tsv", row.names = 1, check.names = FALSE)
pop_counts_all <- round(as.matrix(pop_counts_all)); mode(pop_counts_all) <- "integer"

sample_labels <- c(
  GSM6339917 = "Control_Y1", GSM6339918 = "Control_Y2", GSM6339919 = "Control_Y3",
  GSM6339920 = "Control_Y4", GSM6339921 = "Control_Y5", GSM6339922 = "Control_Y6",
  GSM6339923 = "POP_D1", GSM6339924 = "POP_D2", GSM6339925 = "POP_D3",
  GSM6339926 = "POP_D4", GSM6339927 = "POP_D5", GSM6339928 = "POP_D6",
  GSM6339929 = "POP_Y1", GSM6339930 = "POP_Y2", GSM6339931 = "POP_Y3",
  GSM6339932 = "POP_Y4", GSM6339933 = "POP_Y5", GSM6339934 = "POP_Y6"
)
pop_counts <- pop_counts_all[, names(sample_labels)]
label <- sample_labels[colnames(pop_counts)]
condition <- factor(ifelse(grepl("^POP", label), "POP", "Control"), levels = c("Control", "POP"))
cat("Amostras (18):\n"); print(table(condition))

keep_pop <- rowSums(pop_counts >= 10) >= 9  # metade do menor grupo nao se aplica (6); usa metade das 18
counts_f <- pop_counts[keep_pop, ]
cat("Genes retidos apos filtro (contagem>=10 em >=9 de 18):", nrow(counts_f), "de", nrow(pop_counts), "\n")

coldata <- data.frame(row.names = colnames(counts_f), condition = condition)
dds <- DESeqDataSetFromMatrix(countData = counts_f, colData = coldata, design = ~condition)
dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "POP", "Control"), alpha = 0.05)
res_df <- as.data.frame(res)
res_df$GeneID <- rownames(res_df)
res_df <- res_df[!is.na(res_df$stat), ]
write.csv(res_df, "results/POP_GSE208261_18samples_tissuematched_DESeq2_completo.csv", row.names = FALSE)
cat("DEGs (padj<0.05, |log2FC|>1):", sum(res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1, na.rm = TRUE), "\n\n")

sf_pop <- estimateSizeFactorsForMatrix(counts_f)
pop_logmat <- log2(sweep(counts_f, 2, sf_pop, "/") + 1)
pop_logmat <- pop_logmat[res_df$GeneID, ]
stopifnot(identical(rownames(pop_logmat), res_df$GeneID))

cat("POP ranked gene list:", nrow(res_df), "genes\n")
pop_gsea <- run_gsea_real(
  gene_ids = res_df$GeneID, observed_stat = res_df$stat,
  logmat = pop_logmat, group = condition, gene_sets = kegg_sets
)
write.csv(pop_gsea, "results/GSEA_POP_GSE208261_18samples_tissuematched_KEGG.csv", row.names = FALSE)
cat("POP GSEA (18 amostras) - significativas FDR<0.05:", sum(pop_gsea$p.adjust < 0.05, na.rm = TRUE),
    "| FDR<0.25:", sum(pop_gsea$p.adjust < 0.25, na.rm = TRUE), "de", nrow(pop_gsea), "\n\n")

## =============================================================================
## PART 2 - SUI (tabela 43032_2020_144), mesma fonte de vias KEGG agora
## =============================================================================
cat("============================================================\n")
cat("PART 2: SUI (43032_2020_144)\n")
cat("============================================================\n\n")

sui_tab <- read.csv("results/SUI_43032_2020_144_combined_deduplicated.csv")
sui_tab$EntrezID <- as.character(sui_tab$EntrezID)
score <- build_score(sui_tab$Regulation, sui_tab$pvalue)
ord <- order(-score)
ranked_genes <- sui_tab$EntrezID[ord]
ranked_scores <- score[ord]

sui_gsea <- run_gsea_small_panel(ranked_genes, ranked_scores, kegg_sets)
cat("SUI GSEA - significativas FDR<0.05:", sum(sui_gsea$p.adjust < 0.05, na.rm = TRUE),
    "| FDR<0.25:", sum(sui_gsea$p.adjust < 0.25, na.rm = TRUE), "de", nrow(sui_gsea), "\n\n")
write.csv(sui_gsea, "results/GSEA_SUI_43032_2020_144_usermapping_KEGG.csv", row.names = FALSE)

## =============================================================================
## PART 3 - Vias compartilhadas (FDR<0.25 dos dois lados, sinal de ES)
## =============================================================================
cat("============================================================\n")
cat("PART 3: Vias compartilhadas POP (18 amostras) vs SUI\n")
cat("============================================================\n\n")

pop_sig <- subset(pop_gsea, p.adjust < 0.25)
sui_sig <- subset(sui_gsea, p.adjust < 0.25)
common_ids <- intersect(pop_sig$PATH, sui_sig$PATH)
cat("POP significativas:", nrow(pop_sig), "| SUI significativas:", nrow(sui_sig),
    "| COMPARTILHADAS:", length(common_ids), "\n\n")

if (length(common_ids) > 0) {
  shared <- merge(pop_sig[pop_sig$PATH %in% common_ids, c("PATH","Nh","ES","NES","pvalue","p.adjust")],
                   sui_sig[sui_sig$PATH %in% common_ids, c("PATH","Nh","ES","NES","pvalue","p.adjust")],
                   by = "PATH", suffixes = c("_POP","_SUI"))
  shared$Same_direction <- sign(shared$ES_POP) == sign(shared$ES_SUI)
  print(shared[, c("PATH","Nh_POP","ES_POP","Nh_SUI","ES_SUI","Same_direction","p.adjust_POP","p.adjust_SUI")])
  write.csv(shared, "results/GSEA_shared_KEGG_POP18samples_vs_SUI_usermapping.csv", row.names = FALSE)
  cat("\nConcordantes em direcao:", sum(shared$Same_direction), "de", nrow(shared), "\n")
}

cat("\n=== FIM DO SCRIPT ===\n")
