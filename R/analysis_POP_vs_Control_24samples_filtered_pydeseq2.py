"""
DEG POP (n=12) vs Controle (n=12), GSE208261, usando a matriz ja filtrada
(contagem>=10 em >=12 das 24 amostras, data/GSE208261_raw_counts_filtered_ge10in12of24.tsv).
Design simples ~condition (sem covariavel D/Y), a pedido explicito da usuaria.
Limiar de significancia: FDR<0,05, |log2FC|>1.
"""
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

counts_raw = pd.read_csv("data/GSE208261_raw_counts_filtered_ge10in12of24.tsv", sep="\t", index_col=0)

sample_labels = {
    "GSM6339911": "Control_D1", "GSM6339912": "Control_D2", "GSM6339913": "Control_D3",
    "GSM6339914": "Control_D4", "GSM6339915": "Control_D5", "GSM6339916": "Control_D6",
    "GSM6339917": "Control_Y1", "GSM6339918": "Control_Y2", "GSM6339919": "Control_Y3",
    "GSM6339920": "Control_Y4", "GSM6339921": "Control_Y5", "GSM6339922": "Control_Y6",
    "GSM6339923": "POP_D1", "GSM6339924": "POP_D2", "GSM6339925": "POP_D3",
    "GSM6339926": "POP_D4", "GSM6339927": "POP_D5", "GSM6339928": "POP_D6",
    "GSM6339929": "POP_Y1", "GSM6339930": "POP_Y2", "GSM6339931": "POP_Y3",
    "GSM6339932": "POP_Y4", "GSM6339933": "POP_Y5", "GSM6339934": "POP_Y6",
}

counts = counts_raw.T
counts.index.name = "sample"
condition = ["POP" if sample_labels[s].startswith("POP") else "Control" for s in counts.index]
metadata = pd.DataFrame({"condition": condition}, index=counts.index)

print(f"Genes de entrada (ja filtrados): {counts.shape[1]}")
print(metadata["condition"].value_counts().to_string())

dds = DeseqDataSet(counts=counts, metadata=metadata, design="~condition", refit_cooks=True)
dds.deseq2()

stat = DeseqStats(dds, contrast=["condition", "POP", "Control"], alpha=0.05)
stat.summary()
res = stat.results_df.copy()
res["GeneID"] = res.index
res = res.sort_values("padj")

res.to_csv("results/POP_GSE208261_24samples_nocovariate_DESeq2_completo.csv", index=False)

sig = res[(res["padj"] < 0.05) & (res["log2FoldChange"].abs() > 1)]
sig.to_csv("results/POP_GSE208261_24samples_nocovariate_DEG_sig.csv", index=False)

print(f"\nGenes testados: {res.shape[0]}")
print(f"DEGs significativos (padj<0.05, |log2FC|>1): {sig.shape[0]}")
print(f"  up em POP: {(sig['log2FoldChange']>0).sum()}  down em POP: {(sig['log2FoldChange']<0).sum()}")
