# =============================================================================
# Script 4/4 — Vias compartilhadas entre o GSEA clássico do POP (script 02)
# e o GSEA preranked do Chen 2006/SUI (script 03).
#
# Critério: mesmo ID de via KEGG significativo nos DOIS lados, ao limiar
# padrão de rastreio do GSEA (FDR<0.25) e, separadamente, ao limiar mais
# rígido (FDR<0.05) — reportados os dois para deixar claro o quão frágil (ou
# não) é qualquer sobreposição encontrada.
# =============================================================================

pop <- read.csv("results/GSEA_classic_POP_KEGG.csv")
chen <- read.csv("results/GSEA_preranked_Chen2006_KEGG.csv")

pop$PATH <- sprintf("%05d", as.integer(pop$PATH))
chen$PATH <- sprintf("%05d", as.integer(chen$PATH))

cat("=== GSEA clássico do POP (KEGG) ===\n")
cat("Vias testadas:", nrow(pop), "\n")
cat("Significativas FDR<0.25:", sum(pop$p.adjust < 0.25, na.rm = TRUE),
    " | FDR<0.05:", sum(pop$p.adjust < 0.05, na.rm = TRUE), "\n\n")

cat("=== GSEA preranked do Chen 2006/SUI (KEGG) ===\n")
cat("Vias testáveis:", nrow(chen), "\n")
cat("Significativas FDR<0.25:", sum(chen$p.adjust < 0.25, na.rm = TRUE),
    " | FDR<0.05:", sum(chen$p.adjust < 0.05, na.rm = TRUE), "\n\n")

report_shared <- function(fdr_cut) {
  pop_sig <- subset(pop, p.adjust < fdr_cut)
  chen_sig <- subset(chen, p.adjust < fdr_cut)
  common <- intersect(pop_sig$PATH, chen_sig$PATH)
  cat("--- Limiar FDR<", fdr_cut, " ---\n", sep = "")
  cat("POP significativo:", nrow(pop_sig), " | Chen2006 significativo:", nrow(chen_sig),
      " | COMPARTILHADAS:", length(common), "\n")
  if (length(common) > 0) {
    out <- merge(pop_sig[pop_sig$PATH %in% common, c("PATH","Nh","NES","pvalue","p.adjust","leadingEdge")],
                 chen_sig[chen_sig$PATH %in% common, c("PATH","Nh","NES","pvalue","p.adjust","leadingEdge")],
                 by = "PATH", suffixes = c("_POP", "_Chen2006"))
    print(out[, 1:5])
    return(out)
  }
  cat("Nenhuma.\n")
  data.frame()
}

shared_025 <- report_shared(0.25)
cat("\n")
shared_005 <- report_shared(0.05)

dir.create("results", showWarnings = FALSE)
if (nrow(shared_025) > 0) {
  write.csv(shared_025, "results/GSEA_shared_KEGG_FDR025.csv", row.names = FALSE)
} else {
  writeLines("Nenhuma via KEGG compartilhada entre o GSEA classico do POP e o GSEA preranked do Chen 2006 a FDR<0.25.",
             "results/GSEA_shared_KEGG_FDR025_VAZIO.txt")
}
if (nrow(shared_005) > 0) {
  write.csv(shared_005, "results/GSEA_shared_KEGG_FDR005.csv", row.names = FALSE)
} else {
  writeLines("Nenhuma via KEGG compartilhada entre o GSEA classico do POP e o GSEA preranked do Chen 2006 a FDR<0.05.",
             "results/GSEA_shared_KEGG_FDR005_VAZIO.txt")
}

cat("\n=== FIM ===\n")
cat("Resultado escrito em results/GSEA_shared_KEGG_FDR025* e results/GSEA_shared_KEGG_FDR005*\n")
