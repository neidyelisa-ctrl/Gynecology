# =============================================================================
# Script 26 — A comparação pedida pelo professor, agora com GSEA CLÁSSICO
# (permutação de fenótipo) dos DOIS lados — POP (GSE53868, script 02/20) e
# SUI (Wei2020, script 25) — nenhum dos dois usa mais preranked/permutação
# de rótulo de gene set. Vias compartilhadas, com checagem de direção.
# =============================================================================

pop <- read.csv("results/GSEA_classic_POP_KEGG.csv")
sui <- read.csv("results/GSEA_classic_Wei2020_KEGG.csv")

pop$PATH <- sprintf("%05d", as.integer(pop$PATH))
sui$PATH <- sprintf("%05d", as.integer(sui$PATH))

cat("=== GSEA clássico do POP (GSE53868, KEGG) ===\n")
cat("Vias testadas:", nrow(pop), "| sig FDR<0.05:", sum(pop$p.adjust < 0.05, na.rm=TRUE),
    "| sig FDR<0.25:", sum(pop$p.adjust < 0.25, na.rm=TRUE), "\n\n")

cat("=== GSEA clássico do SUI (Wei2020, KEGG) ===\n")
cat("Vias testadas:", nrow(sui), "| sig FDR<0.05:", sum(sui$p.adjust < 0.05, na.rm=TRUE),
    "| sig FDR<0.25:", sum(sui$p.adjust < 0.25, na.rm=TRUE),
    "(teto de resolucao: p-valor minimo possivel = 0.05, ver script 25)\n\n")

compare <- function(fdr) {
  pop_sig <- subset(pop, p.adjust < fdr)
  sui_sig <- subset(sui, p.adjust < fdr)
  common <- intersect(pop_sig$PATH, sui_sig$PATH)
  cat("--- FDR<", fdr, ": POP sig=", nrow(pop_sig), " | SUI sig=", nrow(sui_sig),
      " | COMPARTILHADAS=", length(common), " ---\n", sep = "")
  if (length(common) == 0) { cat("\n"); return(data.frame()) }
  out <- merge(pop_sig[pop_sig$PATH %in% common, c("PATH","Nh","NES","p.adjust")],
               sui_sig[sui_sig$PATH %in% common, c("PATH","Nh","NES","p.adjust")],
               by = "PATH", suffixes = c("_POP", "_SUI"))
  out$Mesma_direcao <- sign(out$NES_POP) == sign(out$NES_SUI)
  print(out)
  cat("Concordantes:", sum(out$Mesma_direcao, na.rm=TRUE), "de", sum(!is.na(out$Mesma_direcao)), "\n\n")
  out
}

r25 <- compare(0.25)
r05 <- compare(0.05)

dir.create("results", showWarnings = FALSE)
if (nrow(r25) > 0) write.csv(r25, "results/26_shared_GSEA_classico_FDR025.csv", row.names = FALSE)
if (nrow(r05) > 0) write.csv(r05, "results/26_shared_GSEA_classico_FDR005.csv", row.names = FALSE)

cat("=== FIM ===\n")
