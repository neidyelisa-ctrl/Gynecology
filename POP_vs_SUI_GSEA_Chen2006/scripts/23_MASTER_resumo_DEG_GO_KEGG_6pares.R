# =============================================================================
# Script 23 (MASTER) — Consolida os 6 cruzamentos DEG->GO/KEGG (2 datasets
# de POP x 3 painéis de SUI) e procura vias/termos que se repetem em MAIS
# DE UM par POP x SUI - o critério de robustez mais rigoroso deste
# projeto: um achado que sobrevive à troca de dataset de POP e de fonte de
# SUI é muito mais defensável do que um achado de um único par.
#
# Os 6 pares (script que gerou cada um):
#   1. Chen2006  x GSE53868  (script 06)
#   2. Chen2003  x GSE53868  (script 08)
#   3. Wei2020   x GSE53868  (script 14)
#   4. Chen2006  x GSE12852  (script 21)
#   5. Chen2003  x GSE12852  (script 22)
#   6. Wei2020   x GSE12852  (script 18)
# =============================================================================

pares_kegg <- list(
  "Chen2006_x_GSE53868" = "results/06_KEGG_concordantes.csv",
  "Chen2003_x_GSE53868" = "results/08_KEGG_concordantes.csv",
  "Wei2020_x_GSE53868"  = "results/14_KEGG_concordantes.csv",
  "Chen2006_x_GSE12852" = "results/21_KEGG_concordantes.csv",
  "Chen2003_x_GSE12852" = "results/22_KEGG_concordantes.csv",
  "Wei2020_x_GSE12852"  = "results/18_KEGG_concordantes.csv"
)
pares_go <- list(
  "Chen2006_x_GSE53868" = "results/06_GO_BP_concordantes.csv",
  "Chen2003_x_GSE53868" = "results/08_GO_BP_concordantes.csv",
  "Wei2020_x_GSE53868"  = "results/14_GO_BP_concordantes.csv",
  "Chen2006_x_GSE12852" = "results/21_GO_BP_concordantes.csv",
  "Chen2003_x_GSE12852" = "results/22_GO_BP_concordantes.csv",
  "Wei2020_x_GSE12852"  = "results/18_GO_BP_concordantes.csv"
)

read_sig <- function(path, id_col, fdr = 0.05) {
  if (!file.exists(path)) return(NULL)
  d <- read.csv(path)
  if (!id_col %in% names(d)) return(NULL)
  subset(d, p.adjust < fdr)[, c(id_col, "Count", "geneID", "p.adjust")]
}

cat("=== KEGG: vias significativas (FDR<0.05) por par ===\n")
kegg_sig <- lapply(names(pares_kegg), function(nm) {
  d <- read_sig(pares_kegg[[nm]], "TERM_ID")
  cat(sprintf("%-22s: %d vias sig.\n", nm, if (is.null(d)) 0 else nrow(d)))
  if (!is.null(d) && nrow(d) > 0) d$par <- nm
  d
})
names(kegg_sig) <- names(pares_kegg)
kegg_sig <- kegg_sig[sapply(kegg_sig, function(d) !is.null(d) && nrow(d) > 0)]
kegg_all <- if (length(kegg_sig) > 0) do.call(rbind, kegg_sig) else NULL

cat("\n=== KEGG: vias que aparecem em MAIS DE 1 par (o achado mais robusto) ===\n")
if (!is.null(kegg_all) && nrow(kegg_all) > 0) {
  tab <- table(kegg_all$TERM_ID)
  recorrentes <- names(tab)[tab > 1]
  if (length(recorrentes) > 0) {
    out <- kegg_all[kegg_all$TERM_ID %in% recorrentes, ]
    out <- out[order(out$TERM_ID), ]
    print(out[, c("TERM_ID", "par", "Count", "geneID", "p.adjust")])
    write.csv(out, "results/23_KEGG_recorrentes_2ormais_pares.csv", row.names = FALSE)
  } else {
    cat("Nenhuma via KEGG significativa em mais de 1 par.\n")
  }
}

cat("\n\n=== GO BP: termos significativos (FDR<0.05) por par ===\n")
go_sig <- lapply(names(pares_go), function(nm) {
  d <- read_sig(pares_go[[nm]], "TERM_ID")
  cat(sprintf("%-22s: %d termos sig.\n", nm, if (is.null(d)) 0 else nrow(d)))
  if (!is.null(d) && nrow(d) > 0) d$par <- nm
  d
})
names(go_sig) <- names(pares_go)
go_sig <- go_sig[sapply(go_sig, function(d) !is.null(d) && nrow(d) > 0)]
go_all <- if (length(go_sig) > 0) do.call(rbind, go_sig) else NULL

cat("\n=== GO BP: termos que aparecem em MAIS DE 1 par ===\n")
if (!is.null(go_all) && nrow(go_all) > 0) {
  tab <- table(go_all$TERM_ID)
  recorrentes <- names(tab)[tab > 1]
  if (length(recorrentes) > 0) {
    out <- go_all[go_all$TERM_ID %in% recorrentes, ]
    out <- out[order(out$TERM_ID), ]
    print(out[, c("TERM_ID", "par", "Count", "geneID", "p.adjust")])
    write.csv(out, "results/23_GO_recorrentes_2ormais_pares.csv", row.names = FALSE)
  } else {
    cat("Nenhum termo GO significativo em mais de 1 par.\n")
  }
}

## Caso especial: KEGG 04350 (TGF-beta) - ver se está em algum dos 6 pares,
## mesmo que só uma vez (já apareceu em varios lugares neste projeto por
## fora deste teste específico - checagem direta)
cat("\n\n=== Checagem específica: KEGG 04350 (TGF-beta) em cada par ===\n")
for (nm in names(pares_kegg)) {
  path <- pares_kegg[[nm]]
  if (!file.exists(path)) next
  d <- read.csv(path)
  d$TERM_ID <- sprintf("%05d", suppressWarnings(as.integer(d$TERM_ID)))
  hit <- d[d$TERM_ID == "04350", ]
  if (nrow(hit) > 0) {
    cat(sprintf("%-22s: presente (Count=%d, p.adjust=%.4g, genes=%s)\n",
                nm, hit$Count[1], hit$p.adjust[1], hit$geneID[1]))
  } else {
    cat(sprintf("%-22s: não testável ou não presente na tabela KEGG deste par\n", nm))
  }
}

cat("\n=== FIM ===\n")
