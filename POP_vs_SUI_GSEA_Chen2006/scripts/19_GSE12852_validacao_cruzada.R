# =============================================================================
# Script 19 — Duas checagens de validação cruzada envolvendo o GSE12852:
#
#   (a) GSE12852 (uterossacral) x GSE53868 — os DOIS datasets de POP
#       independentes do projeto concordam entre si? (checagem essencial
#       antes de confiar em qualquer resultado que dependa de qual dataset
#       de POP foi usado - especialmente dado o achado surpreendente do
#       script 18: GSE12852 CONCORDA fortemente com o Wei2020 (60%,
#       p=9x10^-36), enquanto o GSE53868 tinha DISCORDADO do mesmo Wei2020
#       (46%, p=10^-7) - script 14. Se os dois POPs nem concordam entre si,
#       isso explica a diferença; se concordam, o puzzle fica mais
#       interessante ainda.)
#
#   (b) GSE12852-GSEA (script 17) x Wei2020-GSEA (script 15) - vias
#       compartilhadas, com checagem de direção (NES), no mesmo formato do
#       script 11.
# =============================================================================

pop53868 <- read.csv("results/GSE53868_limma_completo.csv")
pop12852 <- read.csv("results/GSE12852_uterosacral_limma_completo.csv")

## --- (a) GSE12852 x GSE53868: os dois POPs concordam entre si? -----------
cross_pop <- merge(pop12852[, c("Gene","logFC","P.Value","adj.P.Val")],
                    pop53868[, c("Gene","logFC","P.Value","adj.P.Val")],
                    by = "Gene", suffixes = c("_GSE12852", "_GSE53868"))
cross_pop$Dir_12852 <- ifelse(cross_pop$logFC_GSE12852 > 0, "up", "down")
cross_pop$Dir_53868 <- ifelse(cross_pop$logFC_GSE53868 > 0, "up", "down")
cross_pop$Concordante <- cross_pop$Dir_12852 == cross_pop$Dir_53868

n_concord <- sum(cross_pop$Concordante)
n_total <- nrow(cross_pop)
cat("=== (a) GSE12852 (uterossacral) x GSE53868 - concordância de direção ===\n")
cat("Genes em comum nos dois arrays:", n_total, "\n")
cat("Concordantes:", n_concord, "(", round(100 * n_concord / n_total, 1), "% )\n")
bt <- binom.test(n_concord, n_total, p = 0.5)
cat("H0: concordância = 50% (acaso). p-valor =", format(bt$p.value, digits = 6), "\n\n")

dir.create("results", showWarnings = FALSE)
write.csv(cross_pop, "results/19_GSE12852_x_GSE53868_genes.csv", row.names = FALSE)

## --- (b) GSEA: GSE12852(uterossacral) x Wei2020(completo) ------------------
gse12852_gsea <- tryCatch(read.csv("results/GSEA_classic_GSE12852_uterosacral_KEGG.csv"), error = function(e) NULL)
wei_gsea <- tryCatch(read.csv("results/GSEA_preranked_Wei2020full_KEGG.csv"), error = function(e) NULL)

if (!is.null(gse12852_gsea) && !is.null(wei_gsea)) {
  gse12852_gsea$PATH <- sprintf("%05d", as.integer(gse12852_gsea$PATH))
  wei_gsea$PATH <- sprintf("%05d", as.integer(wei_gsea$PATH))

  compare <- function(a, b, name_a, name_b, fdr) {
    a_sig <- subset(a, p.adjust < fdr)
    b_sig <- subset(b, p.adjust < fdr)
    common <- intersect(a_sig$PATH, b_sig$PATH)
    cat(sprintf("%s (n=%d sig) x %s (n=%d sig), FDR<%.2f -> COMPARTILHADAS: %d\n",
                name_a, nrow(a_sig), name_b, nrow(b_sig), fdr, length(common)))
    if (length(common) > 0) {
      out <- merge(a_sig[a_sig$PATH %in% common, c("PATH","Nh","NES","p.adjust")],
                   b_sig[b_sig$PATH %in% common, c("PATH","Nh","NES","p.adjust")],
                   by = "PATH", suffixes = c(paste0("_", name_a), paste0("_", name_b)))
      nes_cols <- grep("^NES_", names(out), value = TRUE)
      out$Mesma_direcao <- sign(out[[nes_cols[1]]]) == sign(out[[nes_cols[2]]])
      print(out)
      return(out)
    }
    data.frame()
  }

  cat("\n=== (b) GSEA: GSE12852 (uterossacral) x Wei2020 (completo) ===\n")
  r1 <- compare(gse12852_gsea, wei_gsea, "GSE12852", "Wei2020", 0.25)
  r2 <- compare(gse12852_gsea, wei_gsea, "GSE12852", "Wei2020", 0.05)
  if (nrow(r1) > 0) write.csv(r1, "results/19_shared_GSE12852_x_Wei2020_FDR025.csv", row.names = FALSE)
  if (nrow(r2) > 0) write.csv(r2, "results/19_shared_GSE12852_x_Wei2020_FDR005.csv", row.names = FALSE)
} else {
  cat("Arquivos de GSEA não encontrados - rode os scripts 15 e 17 primeiro.\n")
}

cat("\n=== FIM ===\n")
