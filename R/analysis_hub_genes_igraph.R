# =============================================================================
# Hub genes via rede PPI (STRING) + comparacao de vias GO/KEGG entre SUI e POP
#
# Pre-requisito: exportar do STRING (string-db.org -> Multiple proteins ->
# cola o conteudo de results/STRING_input_SUI_36h_genes.txt ou
# STRING_input_POP_genes.txt -> Homo sapiens -> Search):
#   - aba "Exports" -> network (TSV) -> salvar como
#     data/STRING_network_SUI_36h.tsv e data/STRING_network_POP.tsv
#   - aba "Analysis" -> KEGG Pathways -> exportar como
#     data/KEGG_enrichment_POP.tsv (e o de SUI, se houver)
#
# ESTE SCRIPT REPRODUZ A ANALISE REAL rodada em 20/08/2026 com os arquivos
# fornecidos pelo usuario (ver results/STRING_network_*.tsv e
# results/hub_genes_*.csv para os dados/resultados usados).
#
# Hub genes = media dos rankings de grau (degree), betweenness e forca
# ponderada pelo combined_score do STRING - metodo equivalente aos algoritmos
# do cytoHubba (Degree, MCC aproximado via subgrafos densos), sem precisar
# do Cytoscape.
# =============================================================================

if (!requireNamespace("igraph", quietly = TRUE)) install.packages("igraph")
suppressMessages(library(igraph))

compute_hubs <- function(tsv_path, label) {
  edges <- read.delim(tsv_path, check.names = FALSE)
  colnames(edges)[1:2] <- c("node1", "node2")
  edges <- edges[!is.na(edges$node1) & !is.na(edges$node2) & edges$node1 != "", ]

  g <- graph_from_data_frame(edges[, c("node1", "node2", "combined_score")], directed = FALSE)
  g <- simplify(g, remove.multiple = TRUE, remove.loops = TRUE,
                edge.attr.comb = list(combined_score = "max"))

  deg <- degree(g)
  btw <- betweenness(g, weights = NA)
  clo <- closeness(g, weights = NA)
  str_wt <- strength(g, weights = E(g)$combined_score)

  tab <- data.frame(gene = names(deg), degree = deg, betweenness = round(btw, 2),
                     closeness = round(clo, 5), weighted_strength = round(str_wt, 2))
  tab$rank_degree <- rank(-tab$degree, ties.method = "min")
  tab$rank_betweenness <- rank(-tab$betweenness, ties.method = "min")
  tab$rank_strength <- rank(-tab$weighted_strength, ties.method = "min")
  tab$hub_score <- (tab$rank_degree + tab$rank_betweenness + tab$rank_strength) / 3
  tab <- tab[order(tab$hub_score), ]
  rownames(tab) <- NULL

  cat("\n===", label, "===\nNos:", vcount(g), "| Arestas:", ecount(g), "\n")
  print(head(tab, 10))
  write.csv(tab, file.path("results", paste0("hub_genes_", label, ".csv")), row.names = FALSE)
  tab
}

hubs_sui36 <- compute_hubs("results/STRING_network_SUI_36h.tsv", "SUI_36h")
hubs_sui72 <- compute_hubs("results/STRING_network_SUI_72h.tsv", "SUI_72h")
hubs_pop <- compute_hubs("results/STRING_network_POP.tsv", "POP")

## Genes em comum entre cada par de redes inteiras (nao so os hubs)
get_nodes <- function(tsv_path) {
  d <- read.delim(tsv_path)
  unique(c(d[[1]], d[[2]]))
}
sui36_nodes <- get_nodes("results/STRING_network_SUI_36h.tsv")
sui72_nodes <- get_nodes("results/STRING_network_SUI_72h.tsv")
pop_nodes <- get_nodes("results/STRING_network_POP.tsv")

cat("\nGenes em comum SUI_36h x SUI_72h:", length(intersect(sui36_nodes, sui72_nodes)), "\n")
cat("Genes em comum SUI_36h x POP:", length(intersect(sui36_nodes, pop_nodes)), "\n")
cat("Genes em comum SUI_72h x POP:", length(intersect(sui72_nodes, pop_nodes)), "\n")

## KEGG do POP (exportado do STRING)
kegg_pop <- read.delim("results/KEGG_enrichment_POP.tsv")
cat("\nVias KEGG significativas no POP (FDR<0.05):", nrow(kegg_pop), "\n")
print(kegg_pop[, c("term.description", "false.discovery.rate")])
