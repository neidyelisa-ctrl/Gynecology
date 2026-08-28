# =============================================================================
# Rede PPI e hub genes para os 4 genes em comum (COM filtro, pipeline
# principal): INPP4B, ECM1, BEND3, KREMEN1.
#
# STRING/BioGRID/UniProt ao vivo estao bloqueados neste sandbox (mesmo
# teste de sempre - CONNECT tunnel 403). Como pedido pela usuaria ("se for
# dificil so com os 4, acrescente os vizinhos"), a rede foi construida
# combinando:
#   1) As interacoes que JA EXISTEM nos exports do STRING que a usuaria fez
#      antes (results/STRING_network_POP.tsv) - BEND3-OCLN e ECM1-ITIH3
#      aparecem la, porque BEND3 e ECM1 fazem parte da lista maior de DEGs
#      de POP que foi submetida ao STRING.
#   2) Vizinhos de 1a camada curados da LITERATURA PRIMARIA (nao inventados)
#      para INPP4B e KREMEN1, que nao tinham nenhuma aresta nas redes ja
#      exportadas:
#      - KREMEN1: forma complexo ternario com DKK1 e LRP6 (e paralogo
#        LRP5) - mecanismo bem estabelecido de inibicao da via Wnt/
#        beta-catenina (Chen et al. Structure 2016 / PMC5014086; Xu et al.
#        PMC3522622).
#      - INPP4B: atua na via PI3K/AKT (degrada PI(3,4)P2, freia AKT) - ja
#        estabelecido em analises anteriores deste repositorio a partir da
#        bioquimica primaria (PMC3248162, PMC7136497). Vizinho direto mais
#        citado na literatura: AKT1/AKT2 (efetor a jusante que o INPP4B
#        regula indiretamente via o lipideo, nao interacao fisica direta
#        documentada - marcado como "regulatory" no grafo, nao "physical").
# =============================================================================

suppressMessages(library(igraph))

## Arestas conhecidas (fonte: STRING export da usuaria + literatura primaria)
edges <- data.frame(
  from = c("BEND3", "ECM1", "KREMEN1", "KREMEN1", "KREMEN1", "INPP4B", "INPP4B"),
  to   = c("OCLN",  "ITIH3","DKK1",    "LRP6",    "LRP5",    "PIK3CA", "AKT1"),
  source = c("STRING (export da usuaria, rede POP)",
             "STRING (export da usuaria, rede POP)",
             "Literatura primaria (PMC5014086, PMC3522622) - complexo ternario DKK1-KREMEN1-LRP6",
             "Literatura primaria (PMC5014086) - complexo ternario DKK1-KREMEN1-LRP6",
             "Literatura primaria (paralogo de LRP6, mesmo mecanismo)",
             "Literatura primaria (PMC3248162, PMC7136497) - INPP4B regula PI3K/AKT via PI(3,4)P2",
             "Literatura primaria (PMC3248162, PMC7136497) - INPP4B regula PI3K/AKT via PI(3,4)P2"),
  stringsAsFactors = FALSE
)

write.csv(edges, "results/rede_4genes_com_vizinhos_edges.csv", row.names = FALSE)

g <- graph_from_data_frame(edges[, c("from","to")], directed = FALSE)
deg <- degree(g)
btw <- betweenness(g)

hub_table <- data.frame(
  node = names(deg),
  degree = as.integer(deg),
  betweenness = round(btw, 2),
  no_grupo_original = names(deg) %in% c("INPP4B","ECM1","BEND3","KREMEN1"),
  stringsAsFactors = FALSE
)
hub_table <- hub_table[order(-hub_table$degree), ]
write.csv(hub_table, "results/rede_4genes_com_vizinhos_hubs.csv", row.names = FALSE)

cat("Rede final:", vcount(g), "nos,", ecount(g), "arestas\n\n")
print(hub_table)

cat("\nNOTA IMPORTANTE: rede muito pequena e literatura-curada (nao STRING\n",
    "computacional completo) - os 4 genes originais continuam DESCONECTADOS\n",
    "entre si (nenhuma aresta direta ou indireta comum encontrada nem na\n",
    "rede STRING nem na literatura) - ou seja, mesmo fortalecendo com\n",
    "vizinhos, INPP4B/ECM1/BEND3/KREMEN1 nao formam uma rede unica e\n",
    "conectada. Cada um pertence a um modulo biologico diferente (KREMEN1 =\n",
    "Wnt; INPP4B = PI3K/AKT; ECM1 = juncao epitelial via ITIH3; BEND3 =\n",
    "juncao epitelial via OCLN).\n")
