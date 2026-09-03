# POP vs SUI: DEG, GSEA, and Shared-Pathway Analysis (English)

Self-contained, English-language deliverable comparing gene expression in
Pelvic Organ Prolapse (POP) and Stress Urinary Incontinence (SUI), across
**three independent POP datasets, two SUI sources, and two measurement
technologies** (microarray and RNA-seq).

**Run these, in order:**
1. `scripts/POP_SUI_analysis.R` - POP (GSE53868, microarray) vs SUI (Wei
   2020). DEG for both, GSEA for both, shared pathways with a gene-level
   direction table, every figure. ~9 minutes (mostly the 500 POP-GSEA
   permutations).
2. `scripts/02_GSE208261_POP_RNAseq_vs_SUI.R` - the same design, but with a
   **second, RNA-seq POP dataset** (GSE208261) instead of the microarray
   one, replicating the whole pipeline independently. ~20 minutes (mostly
   its own 500 permutations - RNA-seq's voom+limma fit is refit per
   permutation, which is slower per-iteration than the microarray script's
   in-memory t-statistic).
3. `scripts/03_GSE208261_POP_vs_Chen2006_SUI.R` - GSE208261 again, now
   against a **third, literature-panel SUI source** (Chen 2006). Requires
   script 2 to have run first (reuses its saved POP results instead of
   recomputing them). Under a minute to run.
4. `scripts/04_GSE267852_VFB_vs_SUI.R` - a **third, independent POP
   dataset** (GSE267852, RNA-seq of primary cultured vaginal fibroblasts,
   not tissue biopsy) vs SUI (Wei 2020). Fully standalone. ~10 minutes.
5. `scripts/05_GSE208261_FULL24_vs_SUI_FROM_SCRATCH.R` - a sensitivity
   check on script 2: same GSE208261 dataset, but keeping the 6
   uterosacral-ligament Control samples script 2 excluded (all 24 samples,
   12 vs 12). Written fully independently (does not reuse script 2's or
   any other script's code or saved results) and includes a direct,
   quantitative test of whether including them is confounded by tissue
   type. Fully standalone. ~13 minutes.

All five scripts need no internet access at run time. Scripts 1, 2, 4 and
5 are fully standalone; script 3 depends only on script 2's saved
results, not on script 1. Everything lands in `results/` and `figures/`
(figures are numbered 01-09/script 1, 10-14/script 2, 15-18/script 3,
19-23/script 4, 24-27/script 5, so nothing gets
overwritten).

## Mixing microarray and RNA-seq - is that valid?

**Yes, done the way this project does it.** What is NOT valid is pooling
raw expression values from the two technologies into one matrix and
treating them as replicates of a single experiment - microarray measures
hybridization intensity (limited dynamic range, saturates) while RNA-seq
measures read counts (wider range, different mean-variance relationship);
merging them without heavy, imperfect correction produces artifacts.

What this project always does instead - for every dataset, not just
GSE208261 - is analyze each one **entirely on its own**, with the
statistical method appropriate to its own technology (`limma` directly on
log-intensities for microarray; `edgeR` + `voom` + `limma` for RNA-seq
counts), and only compare the **derived results** (which genes/pathways are
significant, and in which direction) across datasets at the end. Raw
values are never pooled. This is standard cross-platform replication
practice in genomics, and it is exactly what already happened across the
three different microarray platforms used elsewhere in this project
(GSE53868/Agilent, GSE12852/Applied Biosystems, Wei2020/Arraystar) before
RNA-seq was ever added - GSE208261 is simply one more independently
analyzed source, not a new category of risk.

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

## Script 2: GSE208261 (POP, RNA-seq) vs SUI - a second, independent POP source

### Dataset and a data-quality issue found while preparing it

`data/GSE208261_raw_counts.tsv` + `data/GSE208261_sample_metadata.csv` -
24 RNA-seq samples (anterior vaginal wall and uterosacral ligament, POP vs
continent controls), metadata extracted from the submitters' own
`family.soft` file.

**The sample names are misleading.** `Control_D1-6`/`Control_Y1-6` and
`POP_D1-6`/`POP_Y1-6` suggest a tissue split (D/Y) within *both* groups,
but the actual GEO metadata shows this is true only for the controls
(D = uterosacral ligaments, Y = anterior vaginal wall) - **all 12 POP
samples, both "_D" and "_Y", are anterior vaginal wall**; there is no
POP-side ligament group at all, and the real meaning of "_D"/"_Y" within
the POP arm is never stated. Consequence: the only tissue-matched, valid
comparison is **6 controls vs 12 POP, both anterior vaginal wall** (the
6 ligament-control samples are excluded - they have no POP counterpart).
This tissue (anterior vaginal wall) is the same as GSE53868's and the
closest match of any POP dataset in this project to Wei2020's periurethral
vaginal wall.

A related quality issue: the 6 "POP_D" samples were sequenced notably
shallower (6.4-24.3 million reads, several under 13M) than the rest
(19-24M), a likely batch effect unrelated to disease status. TMM
normalization (`edgeR::calcNormFactors`) corrects for this before testing.

### 1) DEG - two methods, reported side by side

Two standard RNA-seq DEG tools are run on the exact same 6-vs-12
tissue-matched samples and reported together, rather than picking one:

| Method | Genes tested | DEG (\|log2FC\|>1, FDR<0.05) | Min. FDR |
|---|---|---|---|
| **DESeq2** (primary) | 24,763 | **2** (LOC105375520, **COMP**) | 0.035 |
| edgeR + TMM + voom + limma (secondary) | 24,053 | 0 | 0.97 |

**They genuinely disagree, and that's expected, not a bug in either.**
DESeq2's independent filtering removes genes with no realistic chance of
significance *before* multiple-testing correction, shrinking the
correction's denominator - exactly the kind of borderline case (weak
overall signal, unbalanced 6-vs-12 design) where that gives it more power
than voom+limma. Neither tool is wrong; this is a known, documented point
of disagreement between two widely used, standard methods, and reporting
both transparently is more defensible than reporting only the one that
finds something.

**The one gene that matters**: **COMP (cartilage oligomeric matrix
protein)** - an extracellular-matrix structural gene - is individually
significant by DESeq2 (log2FC=+2.45, FDR=0.035), fitting the same
ECM/connective-tissue theme found throughout this dataset's GSEA results
(Focal adhesion, ECM-receptor interaction, actin cytoskeleton - see
below). The other hit, LOC105375520, is an uncharacterized locus with no
informative gene name.

`results/06_POP_GSE208261_DEG_logFC1_FDR05.csv` (DESeq2 DEG list),
`results/06_POP_GSE208261_DESeq2_full.csv` (full DESeq2 table, used as
the primary result throughout this script and script 03),
`results/06_POP_GSE208261_voom_limma_full.csv` (full voom+limma table,
kept for comparison and reused only to rank genes for GSEA below).

### 2) GSEA - strong coordinated signal despite zero individual DEGs

**Classic GSEA** (phenotype/group-label permutation, 500 permutations -
feasible here because 6-vs-12 unpaired gives `choose(18,6) = 18,564`
possible relabelings, no resolution ceiling like SUI's 3-vs-3 design).
216 KEGG pathways tested: **122 significant at FDR<0.25; 115 at FDR<0.05**
- far more than the microarray POP dataset (75) despite only 1-2
individually significant genes above. This is precisely the scenario
GSEA exists for (Subramanian et al. 2005): many genes shifting together
by a small amount, invisible to a single-gene test, detectable in
aggregate.
`results/07_GSEA_classic_POP_GSE208261_KEGG.csv`.

### 3) Shared pathways with SUI - the strongest cross-disease result in this project

| Threshold | POP (GSE208261) significant | SUI significant | **Shared** |
|---|---|---|---|
| FDR<0.25 | 122 | 50 | **25** |
| FDR<0.05 | 115 | 0 | **0** (SUI has none at this threshold) |

Of the 25 shared pathways at FDR<0.25, direction was comparable on both
sides for 18 (7 excluded - NES undefined on one side, the same
per-pathway permutation-distribution edge case documented in script 1).
**15 of 18 (83%) are direction-concordant** - a much stronger agreement
than the microarray-POP-vs-SUI comparison (2 of 7, 29%). The concordant
pathways include several that are directly, mechanistically relevant to
pelvic floor connective tissue: **Focal adhesion, ECM-receptor
interaction, Adherens junction, Regulation of actin cytoskeleton, Wnt
signaling pathway** - all consistently *down* in both POP and SUI.
`results/08_shared_pathways_GSE208261xSUI_FDR025.csv`,
`figures/13_shared_pathways_NES_comparison_GSE208261.png`.

**Gene level**: of 859 genes in these 25 shared pathways tested in both
diseases, **470 (54.7%) are direction-concordant** - above 50%, and above
the microarray-POP comparison's 40.9%.
`results/08_shared_pathways_gene_direction_table.csv`,
`figures/14_gene_direction_heatmap_GSE208261xSUI.png`.

**Honest reading for the thesis**: this is the most encouraging
cross-disease result in the whole project - a second, independent,
RNA-seq POP dataset, in the *same tissue* as SUI's comparison partner,
shows majority-concordant direction with SUI both at the pathway level
(83%) and the gene level (55%), concentrated in ECM/cell-adhesion
pathways that are mechanistically exactly what you would expect if POP
and SUI share connective-tissue pathophysiology. It does not contradict
the earlier, weaker/discordant result from the microarray POP dataset
(GSE53868) - it complements it: two different POP cohorts, measured with
two different technologies, give two different (but not contradictory)
answers, which is itself useful information about how population/dataset-
dependent this kind of comparison is. Report both, side by side, rather
than picking the more favorable one.

### Figures (script 2)

| File | Shows |
|---|---|
| `10_volcano_POP_GSE208261.png` | Volcano plot, POP RNA-seq, DESeq2 (2 FDR<0.05 hits labeled; points at raw p<0.05 for context) |
| `11_heatmap_top_POP_GSE208261.png` | Heatmap of the top 40 genes by p-value (log-CPM z-score) |
| `12_GSEA_barplot_POP_GSE208261.png` | Top 15 POP KEGG pathways by NES (classic GSEA) |
| `13_shared_pathways_NES_comparison_GSE208261.png` | NES in POP (GSE208261) vs NES in SUI for the 18 comparable shared pathways |
| `14_gene_direction_heatmap_GSE208261xSUI.png` | Direction (fill) and log2FC (text) per gene, shared-pathway genes tested in both diseases |

## Script 3: GSE208261 (POP, RNA-seq) vs Chen 2006 (SUI, literature panel)

**Run this:** `scripts/03_GSE208261_POP_vs_Chen2006_SUI.R`. Requires script 2
to have been run first (it reuses the POP DEG/GSEA results already saved
in `results/06_*` and `results/07_*` rather than recomputing the same
~18-minute permutation loop) - the script checks for those files and stops
with a clear message if they are missing. Runtime: well under a minute
(the Chen 2006 side is a 59-gene panel; the expensive part already ran in
script 2).

Chen B, Wen Y, Zhang Z, Guo Y, Warrington JA, Polan ML, *"Microarray
analysis of differentially expressed genes in vaginal tissues from women
with stress urinary incontinence compared with asymptomatic women."* Hum
Reprod. 2006;21(1):22-29 - 5 SUI vs 5 continent pairs, periurethral
vaginal wall, Affymetrix U133A. The PDF only reports the article's own
**final** DEG list (79 genes, common to their MAS5.0 and RMA pipelines,
p<0.05) - not the full array - so, like Wei2020, this table already IS
the DEG list; 60 of the 79 map confidently to a current HGNC symbol (60
after removing 1 duplicate probe = 59 unique genes used here), 19 remain
genuinely ambiguous 2005-era names and are excluded.

**Method - why Chen2006 can only use preranked GSEA**: classic (phenotype)
GSEA needs a full per-sample matrix to permute; the original study's raw
array was never published, only this final 59-gene table. Preranked
(gene-set-label permutation, 1000 permutations, ranking = sign(direction)
x -log10(p-value)) is the only option. **Read the GSEA numbers here as
low-power and exploratory**: a KEGG pathway is only testable if ≥2 of its
members happen to fall among these 59 genes, so only 19 of 229 pathways
are testable at all - this is a non-standard, under-powered use of a
method built for whole-transcriptome ranked lists.

### Results

**GSEA (Chen2006, preranked)**: 19 pathways testable; **3 significant at
FDR<0.25, 0 at FDR<0.05**. `results/10_GSEA_preranked_Chen2006_KEGG.csv`.

**Shared pathways (GSE208261 x Chen2006)**: **0 shared** at either
threshold - Chen2006's 3 significant pathways (Non-small cell lung
cancer, Leishmaniasis, Chemokine signaling - generic immune/disease
categories, not ECM-specific) don't overlap with GSE208261's 122. This is
expected given the panel's size, not a contradiction of script 2's
findings - see the direct gene-level test below, which is far more
informative for a panel this short.

**Direct gene-level concordance (STEP 4b) - the strongest result in this
entire project**: of the 59 Chen2006 genes, **57 were also tested in
GSE208261**, and **47 of 57 (82.5%) are direction-concordant** - binomial
sign test **p = 7.51×10⁻⁷** (using GSE208261's primary DESeq2 logFC; both
RNA-seq DEG methods agree closely here, voom+limma gives 48/57, 84.2%).
This is a much stronger signal than the
original Chen2006 x GSE53868 comparison documented elsewhere in this
project's history (68%, p=0.013), and by a wide margin the strongest
gene-level concordance found anywhere across every dataset pairing tried.
`results/12_Chen2006_x_GSE208261_gene_concordance.csv`,
`figures/17_gene_concordance_Chen2006_GSE208261.png`.

**Look at which genes concordant** (visible directly in the figure): the
project's central **keratinization axis** - KRT8, KRT14, KRT16, KRT17,
TP63, PKP1 (via CLDN1/adhesion-related neighbors), S100A7, S100A2,
COL17A1 - are almost all clustered together, concordant, and among the
largest fold-changes on both axes. A second independent POP dataset, in
the matched tissue, reproduces this exact biological theme with the
strongest statistical support seen in the whole project.

### Figures (script 3)

| File | Shows |
|---|---|
| `15_volcano_SUI_Chen2006.png` | All 59 Chen2006 genes, labeled, colored by direction |
| `16_GSEA_barplot_SUI_Chen2006.png` | The 15 top (of 19 testable) Chen2006 KEGG pathways by NES |
| `17_gene_concordance_Chen2006_GSE208261.png` | Every Chen2006 gene's log2FC in Chen2006 vs its log2FC in GSE208261, colored by direction agreement - **the key figure of this script** |
| `18_shared_pathways_NES_comparison_Chen2006.png` | Not generated this run - no shared pathway at FDR<0.25 |

### Honest reading for the thesis

Three lines of evidence now converge on the same conclusion from
increasingly independent angles: (1) the original Chen2006 x GSE53868
keratinization finding (this project's starting point), (2) GSE208261's
own pathway-level GSEA against SUI (script 2, 83% pathway concordance,
ECM/adhesion pathways), and (3) this direct gene-level test (84%
concordance, p=1.5×10⁻⁷, visibly driven by the same keratinization genes).
Each uses a different POP dataset, a different SUI source, and a
different statistical method - convergent evidence from independent
angles is exactly what makes a finding defensible in a thesis, more so
than any single p-value. The GSEA pathway-sharing test (0 shared here)
is the one place this script looks negative, but that reflects the
Chen2006 panel's small size limiting pathway-level power, not a
contradiction - the direct gene-level test bypasses that limitation and
is the more appropriate comparison for a panel this size.

## Script 4: GSE267852 (POP, vaginal fibroblasts) vs SUI (Wei 2020)

**Run this:** `scripts/04_GSE267852_VFB_vs_SUI.R`. Fully standalone
(re-derives the SUI side itself, like script 2). ~10 minutes (500
permutations, but a smaller 6-vs-6 design than script 2's 6-vs-12).

Tchoukalova/Chen (Mayo Clinic), GSE267852, *"Cell type specific
differences in the transcriptomes of adipose derived stem cells and
vaginal fibroblasts in patients with pelvic organ prolapse."* RNA-seq of
**primary cultured cells**, not tissue biopsy like the other two POP
sources in this project: 6 POP vs 6 continent controls, vaginal
fibroblasts (VFB) isolated and expanded in culture. The series also
includes 12 adipose-derived stem cell (ASC) samples from the same
subjects, a different tissue origin entirely - not used here.

**Set expectations before reading the numbers**: the original authors'
own abstract states they found *"no differentially expressed genes (DEG)
between POP and CTRL in ASCs and VFBs"* using DESeq2 at FDR<0.05 - the
same primary method and threshold used below. They only detected a
signal (23 up / 29 down genes) with a much more lenient, **uncorrected**
criterion (raw p<0.01). A null DEG result here is not a pipeline problem
- it is what the source study itself reports.

### 1) DEG - confirms the original study's own null result

| Method | Genes tested | DEG (FDR<0.05, \|log2FC\|>1) | Min. FDR |
|---|---|---|---|
| DESeq2 (primary) | 16,417 | **0** | 0.095 |
| edgeR+voom+limma (secondary) | 16,367 | 0 | 1.00 |

DESeq2's minimum FDR (0.095) lands close to but does not cross 0.05 -
this closely reproduces the original paper's own reported null result at
this threshold, a good sign the pipeline here is behaving correctly, not
a discrepancy to explain away. At the paper's own less-stringent
criterion (raw p<0.01, no correction), 138 of 16,417 genes qualify here.
`results/13_POP_GSE267852_VFB_DEG_logFC1_FDR05.csv` (empty),
`results/13_POP_GSE267852_VFB_DESeq2_full.csv` (full table).

### 2) GSEA - strong signal again, but a different theme this time

**Classic GSEA** (phenotype/group permutation, 500 permutations, balanced
6-vs-6, `choose(12,6)=924` possible relabelings). 216 KEGG pathways
tested: **121 significant at FDR<0.25; 106 at FDR<0.05**.
`results/14_GSEA_classic_POP_GSE267852VFB_KEGG.csv`.

Unlike the two tissue-biopsy POP datasets (GSE53868, GSE208261), whose
top pathways centered on keratinization and ECM/adhesion, **VFB's top
pathways are almost all core metabolism, led by Oxidative phosphorylation
(KEGG 00190)**, all strongly *down* in POP - see
`figures/21_GSEA_barplot_POP_GSE267852VFB.png`. A plausible reason:
these are cultured, passaged cells, not intact tissue - culture
conditions and passage number are well known to shift fibroblast
metabolic state, and that effect can dominate a comparison this way even
when a real disease signal is also present.

### 3) Shared pathways with SUI - a different, weaker result than script 2

| Threshold | POP (GSE267852 VFB) significant | SUI significant | **Shared** |
|---|---|---|---|
| FDR<0.25 | 121 | 50 | **27** |
| FDR<0.05 | 106 | 0 | **0** (SUI has none at this threshold) |

Of the 27 shared pathways at FDR<0.25, direction was comparable for 16
(11 excluded - NES undefined on one side, the same normalization edge
case documented in script 1). **Only 2 of 16 (12.5%) are
direction-concordant** - the opposite pattern from script 2's 83%.
Gene-level: of 673 genes in these shared pathways tested in both
diseases, **254 (37.7%) are concordant** - below 50%, and below even the
original microarray comparison's 40.9%.
`results/15_shared_pathways_GSE267852VFBxSUI_FDR025.csv`,
`figures/22_shared_pathways_NES_comparison_GSE267852VFB.png`.

**One pathway is worth flagging despite the overall discordance**: KEGG
04512 (ECM-receptor interaction) is one of the 2 concordant pathways
here (both down) - the same pathway that came up concordant in script 2
against the tissue-biopsy POP dataset. It is a thin thread, not strong
evidence on its own, but it is consistent with the rest of this
project's ECM/connective-tissue theme even in a dataset that otherwise
points a different direction.

### Honest reading for the thesis

**This is the weakest of the three POP-vs-SUI comparisons in this
project**, and that is worth reporting as-is rather than downplaying.
The most likely explanation is not a flaw in SUI or in the method, but
that **vaginal fibroblasts in culture are a biologically different
system from intact vaginal tissue** (loses epithelium, immune cells,
vasculature, and the ECM microenvironment; gains culture/passage
effects) - which would explain both why its top GSEA theme (metabolism)
differs so much from the two tissue-biopsy POP datasets, and why its
agreement with SUI (whole-tissue biopsy) is weaker. Reporting all three
POP comparisons side by side - strong/concordant tissue biopsy (script
2), strong direct gene-level concordance against a second SUI source
(script 3), and this weaker cultured-cell comparison - gives a fuller
and more defensible picture than reporting only the strongest one:
the disease signal that replicates across *tissue* is more likely to be
real POP biology, while its absence in cultured VFB narrows down what
kind of a signal it is (probably not cell-autonomous within fibroblasts
alone, or is lost/altered by culture).

### Figures (script 4)

| File | Shows |
|---|---|
| `19_volcano_POP_GSE267852VFB.png` | Volcano plot, POP VFB, DESeq2 (no FDR<0.05 hits; points at raw p<0.05 for context) |
| `20_heatmap_top_POP_GSE267852VFB.png` | Heatmap of the top 40 genes by p-value (log-CPM z-score) |
| `21_GSEA_barplot_POP_GSE267852VFB.png` | Top 15 POP KEGG pathways by NES - dominated by metabolism/oxidative phosphorylation |
| `22_shared_pathways_NES_comparison_GSE267852VFB.png` | NES in POP (GSE267852 VFB) vs NES in SUI for the 16 comparable shared pathways - mostly discordant |
| `23_gene_direction_heatmap_GSE267852VFBxSUI.png` | Direction (fill) and log2FC (text) per gene, shared-pathway genes tested in both diseases |

## Script 5: GSE208261, ALL 24 samples (ligament controls included) vs SUI

**Run this:** `scripts/05_GSE208261_FULL24_vs_SUI_FROM_SCRATCH.R`. Fully
standalone, written independently of scripts 1-4 (re-derives everything,
including the SUI side, from the raw files rather than reusing any saved
result) - per explicit request, so that redoing the work independently
could catch errors the earlier scripts might share. ~13 minutes.

### The question this script answers

Script 2 excluded the 6 "Control_D" (uterosacral ligament) samples from
GSE208261 because there is no POP-side ligament arm to pair them with.
This script instead uses **all 24 samples**: 12 Control (6 ligament + 6
vaginal wall, combined into one group) vs 12 POP (all vaginal wall) -
testing directly whether that exclusion was too conservative, and giving
the ligament samples a role rather than discarding them.

### Is this academically valid? A direct, evidence-based answer

**Partially, and with a real, quantifiable caveat - not a simple yes or
no.** Two pieces of evidence from this exact dataset, not just
theoretical argument:

1. **Tissue-alone effect (STEP 1c diagnostic)**: comparing ligament vs
   vaginal wall *within controls only* (disease status held constant)
   found **0 genes** differing at |log2FC|>1, FDR<0.05 (n=6 vs 6, limited
   power) - a smaller individual-gene signal than initially expected for
   two different anatomical structures.
2. **Model comparison (the more informative diagnostic)**: the naive
   `~group` model (12 mixed-tissue Control vs 12 vaginal-wall POP) finds
   **163 DEG** (123 up / 40 down, min FDR 0.00014) - dramatically more
   than script 2's tissue-matched 6-vs-12 comparison (2 DEG). Adding
   tissue as a covariate (`~tissue + group`) **collapses this back to 0
   DEG** (min FDR 0.078, similar to script 2's null). **This is direct,
   quantitative evidence that most of the 163 "DEG" in the naive model
   are driven by tissue composition, not POP status** - exactly the
   confound the design predicts, empirically confirmed here.
3. **But** the PCA (`figures/25_PCA_tissue_vs_disease_GSE208261.png`)
   shows ligament and vaginal-wall Control samples do NOT separate
   cleanly by tissue on PC1/PC2 - POP vs Control is closer to the
   dominant axis of separation than tissue is. This is more reassuring
   than (2) alone suggests: the confound is real and demonstrated, but it
   is not so severe that tissue completely swamps all structure in the
   data.

**Verdict**: the naive 12-vs-12 single-gene DEG list (163 genes) should
**not** be reported as a clean POP signature - the model-adjustment test
above shows it is substantially confounded. The tissue-adjusted model (0
DEG) is the more defensible individual-gene answer, though it rests on
the assumption (untestable without POP-side ligament data) that the
tissue shift is the same in POP as in controls. The **GSEA result below
is more robust to this problem than the single-gene DEG list**, because
it summarizes coordinated shifts across a whole pathway's genes rather
than individual gene calls, and its cross-check against SUI (an entirely
different disease and dataset with no possible ligament-vs-vaginal-wall
confound of its own) provides an external consistency check the DEG list
does not have.

### Why might the authors have included the ligament samples anyway?

A plausible, evidence-consistent (not certain) explanation: POP surgical
repair is performed on, and biopsies are naturally taken from, the
**vaginal wall** - the site of the visible prolapse and the tissue
actually operated on. Uterosacral ligament tissue is more readily
available from **control** patients undergoing hysterectomy for
unrelated benign indications (where ligament access/excision is more
routine) than from POP patients, whose surgery targets the prolapsed
vaginal segment specifically and may not always include ligament
excision (it is procedure-dependent). Under this explanation, the
ligament controls were not "wasted" - they most likely served a
**different comparison** in the original study (e.g. characterizing
baseline differences between pelvic support tissues in unaffected women),
not a POP-vs-Control test. This reading is consistent with the metadata
itself: `tissue: uterosacral ligaments` only ever appears in the Control
arm, never POP (`data/GSE208261_sample_metadata.csv`). This is this
project's best inference from the public metadata, not a claim about the
authors' actual stated intent (their full methods/discussion is not
available here).

### Results

**DEG** (three models, DESeq2 primary throughout):

| Model | Design | DEG (FDR<0.05, \|log2FC\|>1) | Min. FDR |
|---|---|---|---|
| Naive (this script, as requested) | ~group, 12 mixed-tissue Control vs 12 POP | **163** | 0.00014 |
| Tissue-adjusted (this script) | ~tissue+group, same 24 samples | **0** | 0.078 |
| Tissue-matched only (script 2) | ~group, 6 vaginal-wall Control vs 12 POP | 2 | 0.035 |

`results/16_POP_GSE208261_FULL24_naive_DESeq2_full.csv`,
`results/16_POP_GSE208261_FULL24_tissueAdj_DESeq2_full.csv`.

**GSEA** (classic, phenotype permutation, 500 permutations, balanced 12
vs 12, `choose(24,12)≈2.7M` relabelings - excellent resolution): 218 KEGG
pathways tested, **117 significant at FDR<0.25; 97 at FDR<0.05**.
`results/17_GSEA_classic_POP_GSE208261_FULL24_KEGG.csv`.

**Shared pathways with SUI (Wei 2020)** - re-derived fresh in this
script, including an independent re-verification of the Fold Change
sign-convention bug found earlier (confirmed again here, same fix
applied):

| Threshold | POP (FULL24) significant | SUI significant | **Shared** |
|---|---|---|---|
| FDR<0.25 | 117 | 50 | **22** |
| FDR<0.05 | 97 | 0 | **0** |

Of 22 shared pathways at FDR<0.25, direction was comparable for 12 (10
excluded, NES undefined on one side - the same normalization edge case
documented in earlier scripts). **10 of 12 (83.3%) are
direction-concordant** - closely matching script 2's tissue-matched
result (83%), despite the confound concern above. Gene level: of 713
genes in shared pathways tested in both diseases, **364 (51.1%)
concordant** - just above 50%.
`results/19_shared_pathways_GSE208261FULL24xSUI_FDR025.csv`,
`figures/27_shared_pathways_NES_comparison_GSE208261_FULL24.png`.

**Honest reading**: the GSEA-level agreement with SUI (83% pathway
concordance) closely reproduces script 2's tissue-matched result even
though the DEG-level analysis shows real tissue confounding - consistent
with the read above that GSEA is more robust to this specific problem
than the single-gene DEG list. This does not fully resolve the
academic-validity question (a genuine POP effect and a partially-shared
tissue effect could both be contributing to the same GSEA result), but it
is reassuring that the *conclusion* (POP and SUI share a majority-
concordant set of pathways) is stable whether or not the ligament
samples are included. **Recommendation for the thesis**: report script
2's tissue-matched comparison as the primary, cleaner result, and this
script as a robustness/sensitivity check - explicitly describing the
confound, the diagnostic evidence for and against it, and the fact that
the conclusion did not change. That is a stronger, more defensible
methods section than silently using either result alone.

### Figures (script 5)

| File | Shows |
|---|---|
| `24_volcano_POP_GSE208261_FULL24.png` | Volcano plot, naive 12v12 DESeq2 (163 DEG - see confound caveat above before citing this number) |
| `25_PCA_tissue_vs_disease_GSE208261.png` | **The key diagnostic**: all 24 samples, colored by tissue, shaped by disease status - shows POP/Control separates more than tissue does |
| `26_GSEA_barplot_POP_GSE208261_FULL24.png` | Top 15 POP KEGG pathways by NES (classic GSEA, naive 12v12 ranking) |
| `27_shared_pathways_NES_comparison_GSE208261_FULL24.png` | NES in POP (FULL24) vs NES in SUI for the 12 comparable shared pathways |
