# =============================================================================
# PPI / hub genes para os genes CONCORDANTES do Chen 2006 (58 genes curados)
# contra os DOIS datasets de POP, feitos SEPARADAMENTE como pedido:
#   - RNA-seq (GSE208261): 45 genes concordantes em direcao (de 55 testaveis)
#   - Microarray (GSE53868): 35 genes concordantes em direcao (de 52 testaveis)
#
# Fonte da rede: reaproveita o STRING_network_POP.tsv ja baixado nesta sessao
# (rede STRING real, baixada quando a usuaria tinha acesso, para os DEGs do
# POP GSE208261) - e a UNICA rede STRING real disponivel localmente. Nao
# fabrica nenhuma aresta - so verifica quais pares do nosso gene set JA
# aparecem conectados nessa rede real (direto ou via 1 vizinho em comum).
# =============================================================================

net <- read.delim("results/STRING_network_POP.tsv", check.names = FALSE)
colnames(net) <- sub("^X\\.", "", colnames(net))
colnames(net)[1] <- "node1"

rnaseq_df <- read.csv("results/Chen2006_79genes_x_POP_completo.csv")
rnaseq_concordant <- unique(subset(rnaseq_df, Testado_no_POP & Concordante)$GeneSymbol)
cat("RNA-seq (GSE208261) - genes concordantes:", length(rnaseq_concordant), "\n")

micro_df <- read.csv("results/Chen2006_x_GSE53868.csv")
micro_concordant <- unique(subset(micro_df, Testado_no_GSE53868 & Concordante)$GeneSymbol)
cat("Microarray (GSE53868) - genes concordantes:", length(micro_concordant), "\n\n")

writeLines(rnaseq_concordant, "results/STRING_input_Chen2006_concordantes_RNAseq_GSE208261.txt")
writeLines(micro_concordant, "results/STRING_input_Chen2006_concordantes_Microarray_GSE53868.txt")

find_direct_edges <- function(genes, net) {
  sub <- net[net$node1 %in% genes & net$node2 %in% genes, ]
  sub[, c("node1", "node2", "combined_score")]
}

find_shared_neighbors <- function(genes, net) {
  genes_in_net <- intersect(genes, union(net$node1, net$node2))
  neigh <- lapply(genes_in_net, function(g) union(net$node2[net$node1 == g], net$node1[net$node2 == g]))
  names(neigh) <- genes_in_net
  out <- list()
  if (length(genes_in_net) >= 2) {
    combs <- combn(genes_in_net, 2, simplify = FALSE)
    for (pair in combs) {
      shared <- intersect(neigh[[pair[1]]], neigh[[pair[2]]])
      if (length(shared) > 0) out[[paste(pair, collapse = "-")]] <- shared
    }
  }
  list(genes_in_net = genes_in_net, neighbors = neigh, shared = out)
}

cat("=== RNA-seq (45 genes concordantes) ===\n")
cat("Genes presentes na rede STRING_network_POP.tsv (real, ja baixada):",
    paste(intersect(rnaseq_concordant, union(net$node1, net$node2)), collapse = ", "), "\n")
edges_rnaseq <- find_direct_edges(rnaseq_concordant, net)
cat("Arestas DIRETAS entre genes concordantes:", nrow(edges_rnaseq), "\n")
sn_rnaseq <- find_shared_neighbors(rnaseq_concordant, net)
cat("Pares com vizinho em comum (conexao indireta, 1 salto):\n")
print(sn_rnaseq$shared)
cat("\n")

cat("=== Microarray (35 genes concordantes) ===\n")
cat("Genes presentes na rede STRING_network_POP.tsv (real, ja baixada):",
    paste(intersect(micro_concordant, union(net$node1, net$node2)), collapse = ", "), "\n")
edges_micro <- find_direct_edges(micro_concordant, net)
cat("Arestas DIRETAS entre genes concordantes:", nrow(edges_micro), "\n")
sn_micro <- find_shared_neighbors(micro_concordant, net)
cat("Pares com vizinho em comum (conexao indireta, 1 salto):\n")
print(sn_micro$shared)

cat("\n=== FIM ===\n")
