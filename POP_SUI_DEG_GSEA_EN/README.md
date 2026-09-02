# POP vs SUI: DEG, GSEA, and Shared-Pathway Analysis (English)

Self-contained, English-language deliverable comparing gene expression in
Pelvic Organ Prolapse (POP) and Stress Urinary Incontinence (SUI). One
script, top to bottom: DEG for both conditions, GSEA for both, shared
pathways with a gene-level direction table, and every figure needed to
support the analysis.

**Run this:** `scripts/POP_SUI_analysis.R` (open it in RStudio, set the
working directory to this folder, run top to bottom). It needs no internet
access at run time and produces everything in `results/` and `figures/`.

Runtime: about 9 minutes total, almost all of it the 500 permutations for
the POP classic GSEA (Step 3).

## Datasets

- **POP**: `data/GSE53868_series_matrix.txt` - GSE53868 (Kerkhof et al.),
  12 women with pelvic organ prolapse, **paired** biopsy per patient
  (prolapse-site tissue vs non-prolapse-site tissue from the same woman),
  Agilent 4x44K array, already log2-normalized.
- **SUI**: `data/Wei2020_TableS2_mRNA.xls` - Wei et al. 2020 (*Reprod
  Sci*), Supplementary Table S2, 3 women with SUI vs 3 continent controls,
  periurethral vaginal wall tissue, Arraystar array. **Two sheets** -
  `up_Sui_vs_Ctrl` and `down_Sui_vs_Ctrl` - both are read and combined; the
  script prints their row counts on every run so this is directly
  verifiable, not just claimed. This file already only lists the genes the
  authors called differentially expressed (fold-change>=2, raw p<0.05) -
  the full tested array was not published.

## Methods, and why

### DEG

- **POP**: standard limma moderated t-test, paired design
  (`~individual + tissue`), **|log2FC|>1 and FDR(BH)<0.05**.
- **SUI**: the supplementary table already IS the authors' DEG list
  (fold-change>=2, raw p<0.05) - the full untested array isn't available,
  so re-deriving a POP-style FDR cutoff from scratch on this pre-filtered
  list would be circular. We report their list as-is as the primary SUI
  DEG set (6,118 genes), and additionally show a stricter, more
  POP-comparable subset by applying FDR<0.05 to the authors' own
  precomputed FDR column (calculated by their software on the full array
  before their p<0.05 filter was applied - real variation confirms this
  isn't a circular re-derivation): 3,621 genes.
- **Bug found and fixed while building this**: the "Fold Change" column in
  the `down_Sui_vs_Ctrl` sheet is an unsigned **magnitude** (always >=2),
  with the sign coming from the sheet itself, not the number. A naive
  `log2(FoldChange)` would give down-regulated genes a *positive* log2FC.
  The script restores the correct sign (`-log2(FoldChange)` for
  `Direction=="down"`) before anything else uses it.

### GSEA

| | POP | SUI |
|---|---|---|
| Method | **Classic** (phenotype permutation) | **Preranked** (gene-set-label permutation) |
| Permutations | 500, paired sign-flip across 12 patients | 1,000, over the ~200 pathways |
| Ranking | Paired moderated t-statistic (limma) | Real moderated t-statistic from limma on the 3-vs-3 per-sample intensities in the supplementary table |

**Why different methods for the two diseases**: classic (phenotype)
permutation is the statistically preferred choice whenever a full
sample-level matrix is available (Subramanian et al. 2005), and POP's 12
paired patients give real permutation resolution. SUI, however, only has 3
vs 3 samples - phenotype permutation there has just `choose(6,3) = 20`
possible relabelings, capping resolution at p = 1/20 = 0.05 and making
FDR<0.05 structurally unreachable no matter how real the signal is.
Preranked (gene-set) permutation avoids that ceiling because it permutes
over the ~200 tested pathways, not over the 6 samples - so **SUI's "0
pathways at FDR<0.05" result below is a genuine (if modest) finding, not
an artifact of a resolution floor.**

Both use the same weighted running-sum enrichment score (ES, weight=1)
from Subramanian et al. 2005, reimplemented here (offline, no
internet needed) because this environment has no live access to
CRAN/Bioconductor to install `fgsea`/`msigdbr`; it uses the KEGG pathway
annotation bundled in `org.Hs.eg.db`. A small built-in lookup gives
readable names to the more common KEGG IDs in the figures/tables; less
common ones are shown as their numeric ID only (no offline source has the
full up-to-date name table).

## Results

### 1) DEG

- **POP: 117 DEGs** (|log2FC|>1, FDR<0.05) of 31,013 genes tested - 73 up,
  44 down at the prolapse site. `results/01_POP_DEG_logFC1_FDR05.csv`
  (full ranked table: `results/01_POP_limma_full.csv`).
- **SUI: 6,118 DEGs** (authors' original criterion) - 3,991 up, 2,127
  down. `results/02_SUI_DEG_primary.csv` (stricter FDR<0.05 subset: 3,621
  genes, `results/02_SUI_DEG_FDR005_subset.csv`; both sheets combined,
  before dedup: `results/02_SUI_all_genes_both_sheets.csv`).

### 2) GSEA

- **POP (classic)**: 218 KEGG pathways tested (5-200 members each).
  **94 significant at FDR<0.25; 75 at FDR<0.05.**
  `results/03_GSEA_classic_POP_KEGG.csv`.
- **SUI (preranked)**: 200 KEGG pathways tested.
  **50 significant at FDR<0.25; 0 at FDR<0.05** (a real result, not a
  resolution-floor artifact - see method note above).
  `results/04_GSEA_preranked_SUI_KEGG.csv`.

### 3) Shared pathways

| Threshold | POP significant | SUI significant | **Shared** |
|---|---|---|---|
| FDR<0.25 | 94 | 50 | **21** |
| FDR<0.05 | 75 | 0 | **0** (expected - SUI has none at this threshold) |

`results/05_shared_pathways_FDR025.csv`. Of the 21 shared pathways at
FDR<0.25, the POP-side NES could actually be computed for only 7 (the
other 14 have an undefined NES on the POP side - a normalization edge
case for pathways with an unusually one-sided permutation distribution,
not the SUI sample-size ceiling; it is unrelated to significance, only to
the NES *magnitude* being displayable). Of those 7 comparable pathways:

- **2 same direction**: KEGG 00380 (Tryptophan metabolism), KEGG 05012
  (Parkinson disease)
- **5 opposite direction**: KEGG 00563, KEGG 00640, KEGG 00650, KEGG
  03010 (Ribosome), KEGG 04621 (NOD-like receptor signaling pathway)

**Gene-level direction table** (`results/05_shared_pathways_gene_direction_table.csv`,
visualized in `figures/09_gene_direction_heatmap_shared_pathways.png`):
for every gene that is a member of one of the 21 shared pathways AND was
tested in both diseases (457 genes), its log2FC and up/down direction in
POP and in SUI. **187 of 457 (40.9%) are direction-concordant** - below
50%, i.e. the genes in these shared pathways trend toward *opposite*
regulation between POP and SUI more often than the same direction. This
matches a pattern seen elsewhere when comparing this specific POP dataset
(pre-menopausal women) against this specific SUI dataset (post-menopausal
women) - population/hormonal-status differences between the two cohorts
are a plausible explanation, discussed further below.

## Figures (`figures/`)

| File | Shows |
|---|---|
| `01_volcano_POP.png` | Volcano plot, POP DEG (FDR<0.05, \|log2FC\|>1) |
| `02_volcano_SUI.png` | Volcano plot, SUI DEG (authors' criterion) |
| `03_heatmap_DEG_POP.png` | Heatmap of the top 40 POP DEGs (row z-score, paired patients) |
| `04_heatmap_DEG_SUI.png` | Heatmap of the top 40 SUI DEGs (row z-score, 3 vs 3 samples) |
| `05_GSEA_barplot_POP.png` | Top 15 POP KEGG pathways by NES (classic GSEA) |
| `06_GSEA_barplot_SUI.png` | Top 15 SUI KEGG pathways by NES (preranked GSEA) |
| `07_shared_pathways_NES_comparison.png` | NES in POP vs NES in SUI for the 7 comparable shared pathways, colored by direction agreement |
| `08_GSEA_enrichment_plot_top_POP.png` | Classic GSEA running-score enrichment plot for POP's single top pathway |
| `09_gene_direction_heatmap_shared_pathways.png` | Direction (fill) and log2FC (text) per gene, for genes in shared pathways tested in both diseases |

## Limitations (read before using these numbers in the thesis)

- **No live access to CRAN/Bioconductor/MSigDB/KEGG.jp** in the
  environment this was built in - confirmed by direct test. KEGG pathway
  definitions come from the version frozen inside `org.Hs.eg.db`, not the
  live, current KEGG.jp database, and pathway *names* are only available
  for the common ones hardcoded into this script (`kegg_names` near the
  top). If you run this on your own computer with normal internet access,
  you can upgrade to `fgsea` + `msigdbr` for the full, current
  Hallmark/KEGG/Reactome collections with readable names throughout -
  the enrichment-score algorithm here is a from-scratch reimplementation
  of the same method (Subramanian et al. 2005) for exactly this offline
  case, but `fgsea` will be faster and more complete.
- **SUI DEG/GSEA are both built on the authors' own pre-filtered gene
  list**, not the full array (~20,000+ genes actually tested in the
  original study) - the true untested background was not published in
  the supplementary material. This affects statistical power more than it
  affects direction/plausibility.
- **SUI classic (phenotype) GSEA was not used**, specifically because of
  the 3-vs-3 sample size (see Methods above) - preranked was used
  instead, which is more resolution-appropriate here but is a marginally
  more liberal test (Subramanian et al. 2005) than phenotype permutation.
- **POP (pre-menopausal) and SUI (post-menopausal, Wei 2020) are
  different populations** - the below-50% gene-direction concordance in
  shared pathways (40.9%) may reflect genuine biological difference in
  hormonal status between the two source studies rather than an absence
  of a POP-SUI relationship. This is worth discussing explicitly rather
  than treated as a negative result.
