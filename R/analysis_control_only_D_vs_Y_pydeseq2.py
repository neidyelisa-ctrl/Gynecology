"""
A pedido explicito da usuaria: DEG test SO entre as 12 amostras de Controle do
GSE208261, comparando as duas subdivisoes de 6 vs 6 ja usadas no projeto como
covariavel "age_group" (D = idosa, Y = jovem) nas analises POP-vs-Controle
anteriores (ver R/pop_age_shrink_4combos_diagnostic.R, README linha 80). A
usuaria descreveu essa divisao como "tipo de tecido" - o rotulo real nos
metadados do GEO/no pipeline e grupo etario (idosa vs jovem), nao um sitio
anatomico diferente (nao ha anotacao de tecido diferente por amostra no
GSE208261). Mantive a comparacao exatamente como pedida (D vs Y dentro dos
12 Controles), so documentando a diferenca de rotulo.

Usa pydeseq2 (reimplementacao Python do DESeq2, mesmo metodo estatistico -
binomial negativa + teste de Wald) porque este ambiente nao tem R instalado.
"""
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

counts_raw = pd.read_csv("data/GSE208261_raw_counts_GRCh38.p13_NCBI.tsv", sep="\t", index_col=0)

sample_labels = {
    "GSM6339911": "Control_D1", "GSM6339912": "Control_D2", "GSM6339913": "Control_D3",
    "GSM6339914": "Control_D4", "GSM6339915": "Control_D5", "GSM6339916": "Control_D6",
    "GSM6339917": "Control_Y1", "GSM6339918": "Control_Y2", "GSM6339919": "Control_Y3",
    "GSM6339920": "Control_Y4", "GSM6339921": "Control_Y5", "GSM6339922": "Control_Y6",
}
control_samples = list(sample_labels.keys())

counts = counts_raw[control_samples].T  # pydeseq2 quer amostras nas linhas
counts.index.name = "sample"

metadata = pd.DataFrame(
    {"age_group": ["D" if "_D" in sample_labels[s] else "Y" for s in control_samples]},
    index=control_samples,
)

# Filtro de baixa contagem: padrao DESeq2 - manter genes com contagem >=10 em
# pelo menos no numero de amostras do menor grupo (6, ja que e 6 vs 6 aqui).
keep = (counts_raw[control_samples] >= 10).sum(axis=1) >= 6
counts_filtered = counts.loc[:, keep.values]
print(f"Genes antes do filtro: {counts_raw.shape[0]}, depois: {counts_filtered.shape[1]}")

dds = DeseqDataSet(
    counts=counts_filtered,
    metadata=metadata,
    design="~age_group",
    refit_cooks=True,
)
dds.deseq2()

stat = DeseqStats(dds, contrast=["age_group", "D", "Y"], alpha=0.05)
stat.summary()
res = stat.results_df.copy()
res["GeneID"] = res.index
res = res.sort_values("padj")

res.to_csv("results/POP_GSE208261_Control_D_vs_Y_DESeq2_completo.csv", index=False)

sig = res[(res["padj"] < 0.05) & (res["log2FoldChange"].abs() > 1)]
sig.to_csv("results/POP_GSE208261_Control_D_vs_Y_DEG_sig.csv", index=False)

print(f"\nGenes testados: {res.shape[0]}")
print(f"DEGs significativos (padj<0.05, |log2FC|>1): {sig.shape[0]}")
print("\nTop 20 por padj:")
print(res.head(20)[["GeneID", "log2FoldChange", "pvalue", "padj"]].to_string(index=False))
