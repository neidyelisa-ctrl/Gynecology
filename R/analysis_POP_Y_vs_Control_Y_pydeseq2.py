"""
DEG POP_Y (n=6) vs Control_Y (n=6), GSE208261 - comparacao tissue-matched
(ambos "parede vaginal anterior" pelos metadados reais do GEO).
"""
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

counts_raw = pd.read_csv("data/GSE208261_raw_counts_GRCh38.p13_NCBI.tsv", sep="\t", index_col=0)

sample_labels = {
    "GSM6339917": "Control_Y1", "GSM6339918": "Control_Y2", "GSM6339919": "Control_Y3",
    "GSM6339920": "Control_Y4", "GSM6339921": "Control_Y5", "GSM6339922": "Control_Y6",
    "GSM6339929": "POP_Y1", "GSM6339930": "POP_Y2", "GSM6339931": "POP_Y3",
    "GSM6339932": "POP_Y4", "GSM6339933": "POP_Y5", "GSM6339934": "POP_Y6",
}
samples = list(sample_labels.keys())
counts_sub = counts_raw[samples]

keep = (counts_sub >= 10).sum(axis=1) >= 6
counts_filtered = counts_sub.loc[keep]
print(f"Genes retidos apos filtro: {counts_filtered.shape[0]} de {counts_raw.shape[0]}")

counts_t = counts_filtered.T
counts_t.index.name = "sample"
condition = ["POP" if sample_labels[s].startswith("POP") else "Control" for s in counts_t.index]
metadata = pd.DataFrame({"condition": condition}, index=counts_t.index)

dds = DeseqDataSet(counts=counts_t, metadata=metadata, design="~condition", refit_cooks=True)
dds.deseq2()

stat = DeseqStats(dds, contrast=["condition", "POP", "Control"], alpha=0.05)
stat.summary()
res = stat.results_df.copy()
res["GeneID"] = res.index
res = res.sort_values("padj")

res.to_csv("results/POP_Y_vs_Control_Y_DESeq2_completo.csv", index=False)
sig = res[(res["padj"] < 0.05) & (res["log2FoldChange"].abs() > 1)]
sig.to_csv("results/POP_Y_vs_Control_Y_DEG_sig.csv", index=False)

print(f"\nGenes testados: {res.shape[0]}")
print(f"DEGs significativos (padj<0.05, |log2FC|>1): {sig.shape[0]}")
print(f"  up em POP_Y: {(sig['log2FoldChange']>0).sum()}  down em POP_Y: {(sig['log2FoldChange']<0).sum()}")
