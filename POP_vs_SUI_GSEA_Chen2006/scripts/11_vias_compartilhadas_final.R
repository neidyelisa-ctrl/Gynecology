# =============================================================================
# Script 11 — Consolidação final: vias compartilhadas entre POP (2 métodos
# de GSEA: normal/fenótipo e preranked/rótulo de gene set) e os DOIS artigos
# do Chen (2003 e 2006), todos testados par a par.
# =============================================================================

pop_normal    <- read.csv("results/GSEA_classic_POP_KEGG.csv")
pop_preranked <- read.csv("results/GSEA_preranked_POP_KEGG_genesetpermutation.csv")
chen06 <- read.csv("results/GSEA_preranked_Chen2006_KEGG.csv")
chen03 <- read.csv("results/GSEA_preranked_Chen2003_KEGG.csv")
wei20  <- if (file.exists("results/GSEA_preranked_Wei2020full_KEGG.csv")) read.csv("results/GSEA_preranked_Wei2020full_KEGG.csv") else NULL

for (df in c("pop_normal", "pop_preranked", "chen06", "chen03")) {
  d <- get(df); d$PATH <- sprintf("%05d", as.integer(d$PATH)); assign(df, d)
}
if (!is.null(wei20)) wei20$PATH <- sprintf("%05d", as.integer(wei20$PATH))

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

cat("=========================================================\n")
cat("POP normal (fenótipo) vs Chen 2006 (preranked)\n")
cat("=========================================================\n")
r1a <- compare(pop_normal, chen06, "POP_normal", "Chen2006", 0.25)
r1b <- compare(pop_normal, chen06, "POP_normal", "Chen2006", 0.05)

cat("\n=========================================================\n")
cat("POP preranked (rótulo de gene set) vs Chen 2006 (preranked) - MESMO tipo de permutação nos dois lados\n")
cat("=========================================================\n")
r2a <- compare(pop_preranked, chen06, "POP_preranked", "Chen2006", 0.25)
r2b <- compare(pop_preranked, chen06, "POP_preranked", "Chen2006", 0.05)

cat("\n=========================================================\n")
cat("POP normal (fenótipo) vs Chen 2003 (preranked)\n")
cat("=========================================================\n")
r3a <- compare(pop_normal, chen03, "POP_normal", "Chen2003", 0.25)
r3b <- compare(pop_normal, chen03, "POP_normal", "Chen2003", 0.05)

cat("\n=========================================================\n")
cat("POP preranked (rótulo de gene set) vs Chen 2003 (preranked) - MESMO tipo de permutação nos dois lados\n")
cat("=========================================================\n")
r4a <- compare(pop_preranked, chen03, "POP_preranked", "Chen2003", 0.25)
r4b <- compare(pop_preranked, chen03, "POP_preranked", "Chen2003", 0.05)

all_results <- list(
  POP_normal_x_Chen2006_FDR025 = r1a, POP_normal_x_Chen2006_FDR005 = r1b,
  POP_preranked_x_Chen2006_FDR025 = r2a, POP_preranked_x_Chen2006_FDR005 = r2b,
  POP_normal_x_Chen2003_FDR025 = r3a, POP_normal_x_Chen2003_FDR005 = r3b,
  POP_preranked_x_Chen2003_FDR025 = r4a, POP_preranked_x_Chen2003_FDR005 = r4b
)

if (!is.null(wei20)) {
  cat("\n=========================================================\n")
  cat("POP normal (fenótipo) vs Wei 2020 (preranked)\n")
  cat("=========================================================\n")
  r5a <- compare(pop_normal, wei20, "POP_normal", "Wei2020", 0.25)
  r5b <- compare(pop_normal, wei20, "POP_normal", "Wei2020", 0.05)
  cat("\n=========================================================\n")
  cat("POP preranked (rótulo de gene set) vs Wei 2020 (preranked)\n")
  cat("=========================================================\n")
  r6a <- compare(pop_preranked, wei20, "POP_preranked", "Wei2020", 0.25)
  r6b <- compare(pop_preranked, wei20, "POP_preranked", "Wei2020", 0.05)
  all_results$POP_normal_x_Wei2020_FDR025 <- r5a
  all_results$POP_normal_x_Wei2020_FDR005 <- r5b
  all_results$POP_preranked_x_Wei2020_FDR025 <- r6a
  all_results$POP_preranked_x_Wei2020_FDR005 <- r6b
}

dir.create("results", showWarnings = FALSE)
for (nm in names(all_results)) {
  if (nrow(all_results[[nm]]) > 0) {
    write.csv(all_results[[nm]], paste0("results/11_shared_", nm, ".csv"), row.names = FALSE)
  }
}

cat("\n=== RESUMO ===\n")
for (nm in names(all_results)) cat(nm, ":", nrow(all_results[[nm]]), "vias compartilhadas\n")
cat("\n=== FIM ===\n")
