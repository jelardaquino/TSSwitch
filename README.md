# TSSwitch

**Transcriptomic Splicing Switch**

Code for *Isoform Remodeling and Promoter Switching in Alzheimer's Disease Frontal Cortex: An Exploratory Long-Read RNA-seq Analysis* (Aquino, Lee, Zheng, Ku, Han, Chen).

---

## Overview

TSSwitch analyzes Oxford Nanopore long-read RNA-seq from postmortem frontal cortex (6 Alzheimer's disease [AD], 6 controls [CT]) to ask whether AD-associated transcriptome changes extend beyond gene-level abundance to isoform proportions, promoter selection, and transcript structure. Findings are then tested in an independent short-read cohort of 41 prefrontal cortex donors.

Analyses:

1. **Differential isoform usage (DIU)** with IsoformSwitchAnalyzeR (DEXSeq, three surrogate variables), with satuRn and CPM-based sensitivity analyses
2. **Promoter switching**: promoter usage from transcript TPM, tested by beta regression
3. **Transcription factor motif enrichment** at differentially used promoters (PWMEnrich)
4. **Targeted event-level splicing**: SUPPA2 PSI for the DIU genes, tested by beta regression (additive and condition-by-sex models)
5. **Gene- and transcript-level differential expression** of the DIU genes (DESeq2)
6. **Independent short-read replication** in E-MTAB-14167 (20 AD, 21 CT), quantified against the long-read Bambu annotation

### Main results

| Analysis | Result |
|---|---|
| DIU (FDR < 0.05, \|dIF\| > 0.1) | 109 isoform-switching events, 96 isoforms, 86 genes (e.g. *CLSTN3*, *SNCA*, *USP10*, *NSD2*) |
| NMD | PTC status assigned for 93 of 96 switching isoforms; 15 predicted NMD-sensitive |
| Promoter switching (BH-adjusted p ≤ 0.05) | 11 genes; 46 differentially used promoters (31 increased, 15 decreased in AD) |
| TF motifs (nominal p ≤ 0.05) | 179 motifs, 138 TFs; IRF9/IRF8/IRF7 and STAT2::STAT1 top-ranked at increased-usage promoters, NF-κB2/NF-κB1 at decreased-usage promoters |
| Event-level PSI (FDR < 0.05, \|dPSI\| ≥ 0.10) | 14 events in the additive model (4 in *PLSCR3*); 14 in the condition-by-sex model (e.g. *KDSR*, *CXXC1*, *C11orf80*) |
| Gene-level DE | 1 of 85 tested DIU genes is a DEG (*LSP1*); 70 of 86 show neither gene- nor transcript-level DE |
| satuRn check | 72 of the 86 DEXSeq DIU genes also significant with satuRn |
| CPM check | 74/86 DIU genes, 83/96 isoforms, 35/46 promoters, 9/11 promoter-switching genes, 11/14 PSI events per model remain significant |
| Short-read replication (41 donors) | Reciprocal *CLSTN3* switch and one *COL4A2* isoform switch replicated (same isoform, same direction); no promoter or PSI target reached FDR < 0.05 |

These results are exploratory. Only the *CLSTN3* and *COL4A2* switches replicated, so the remaining findings need larger long-read cohorts.

---

## Repository structure

```
TSSwitch/
├── syn52047893/                          # Long-read discovery analysis (6 AD vs 6 CT)
│   ├── tx_counts_to_tpm.r                    # Bambu counts -> TPM
│   ├── pca_tpm.r                             # PCA of transcript TPM (Fig. 2A)
│   ├── isoformswitchanalyzer_DEXSeq.r        # DIU: IsoformSwitchAnalyzeR + DEXSeq (Fig. 2B-C, Supp. Table 4)
│   ├── isoformswitchanalyzer_satuRn.r        # DIU with satuRn, same import/filter/SVs (Supp. Table 1)
│   ├── generate_all_significant_switch_plots_ORF_annotated.R  # Switch plots (Fig. 2D-E)
│   ├── custom_switch_plot_helpers.R          # ggplot2 helpers for switch plots
│   ├── example_custom_plots.R                # Example use of the helpers
│   ├── volcano_isoform_switches.r            # DIU volcano plot
│   ├── bar_plot_diu_switch_features.r        # Counts of switching genes/isoforms
│   ├── promoter_usage_analysis_86_genes_12samples.r  # Promoter usage beta regression (Fig. 3A-B, Supp. Table 6)
│   ├── circos_heatmap_atac_plots_86_genes_12samples.r # PWMEnrich motif enrichment + TF-promoter network (Fig. 3D-E, Supp. Table 7)
│   ├── run_regtools.sh                       # Splice junction extraction (RegTools)
│   ├── splicing_graph_tx_counts.r            # Splicing graph of promoter/transcript/junction counts (Fig. 3C)
│   ├── suppa_generate_events_bambu_gtf.sh    # SUPPA2 events and PSI from TPM
│   ├── diff_beta_reg_sig_diu_tx.r            # PSI beta regression, additive model (Supp. Table 3a)
│   ├── betareg_sex.r                         # PSI beta regression, condition-by-sex model (Fig. 4, Supp. Table 3b)
│   ├── gprofiler2_dotplot.r                  # g:Profiler GO enrichment (not reported in the manuscript)
│   ├── cpm_sensitivity/                      # All analyses repeated with CPM in place of TPM (Supp. Table 2)
│   │   ├── tx_counts_to_cpm.r                    # Bambu counts -> CPM
│   │   ├── suppa_psi_cpm.sh                      # SUPPA2 PSI from CPM (same IOE events)
│   │   ├── 01_diu_dexseq_cpm.r                   # DIU on CPM with the TPM-run surrogate variables; ORF/PTC calls
│   │   ├── 02a_psi_betareg_additive_cpm.r        # PSI, additive model
│   │   ├── 02b_psi_betareg_condition_x_sex_cpm.r # PSI, condition-by-sex model
│   │   ├── 03_promoter_usage_cpm.r               # Promoter usage
│   │   ├── 04_motif_tf_cpm.r                     # Motif enrichment
│   │   ├── 05_promoter_isoform_merge_cpm.r       # Link switching promoters to switching isoforms
│   │   ├── 06_compare_tpm_vs_cpm.r               # TPM vs CPM overlap report
│   │   ├── 07_supplementary_table_tpm_vs_cpm.r   # TPM vs CPM supplementary table
│   │   ├── diag_diu_abundance_vs_batch.r         # Diagnostic: abundance measure vs surrogate-variable effect
│   │   └── run_cpm_pipeline.sh                   # Runs steps 01-06
│   └── gene_tx_DE/                           # DIU genes vs gene- and transcript-level DE (Section 2.8)
│       ├── 08_reanalysis_deg_txdeg_diu_venn.R    # Gene DEG / transcript DEG / DIU Venn
│       ├── 09_supp_figure_deg_txdeg_diu.R        # Supplementary Figure 2
│       ├── 10_supp_table9_gene_tx_DE.R           # Per-gene table (Supplementary Table 5)
│       └── 10b_supp_table9_xlsx.py               # Writes that table to .xlsx
│
└── E-MTAB-14167/                         # Independent short-read replication (20 AD vs 21 CT)
    ├── salmon_bambu_index_quant.sh           # Salmon index from the Bambu extended annotation
    ├── 01_salmon_bambu_quant_emtab14167.sh   # Salmon quantification of all 41 PFC samples
    ├── 02_targeted_replication.R             # Promoter and PSI targets retested by beta regression (Supp. Table 9)
    ├── 03_diu_bambu_annotation_emtab14167.R  # Short-read DIU and long- vs short-read comparison (Supp. Table 8)
    ├── 04_replicated_feature_plots.R         # Per-sample plots of replicated features (Supp. Fig. 5)
    ├── 05_diu_replication_summary_figure.R   # DIU overlap and dIF concordance figure
    └── deconvolution/
        └── 01_deconvolution.R                # Neuronal fraction estimates read by 02 (sensitivity sets only)
```

Some output and script names keep earlier labels: `_9samples` files from the promoter script contain all 12 samples, and `10_supp_table9_*` writes what is now Supplementary Table 5.

---

## Datasets

| Folder | Accession | Description | Platform | Samples |
|---|---|---|---|---|
| `syn52047893/` | Synapse [syn52047893](https://www.synapse.org/#!Synapse:syn52047893); SRA SRP456327 | Postmortem frontal cortex, Heberle *et al.* 2025 (*Nat Biotechnol*) | Oxford Nanopore cDNA | 6 AD, 6 CT |
| `E-MTAB-14167/` | ArrayExpress [E-MTAB-14167](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-14167) (ENA ERP161086) | Postmortem prefrontal cortex | Illumina 150-bp paired-end | 20 AD (10 F, 10 M; age 82.3 ± 7.4), 21 CT (10 F, 11 M; age 88.6 ± 5.2) |

> ⚠️ **Data access:** The long-read data require an approved Synapse account and data use agreement through the [AD Knowledge Portal](https://adknowledgeportal.synapse.org/), or can be obtained from SRA (SRP456327). No raw or processed data are included in this repository.

---

## Methods summary

### Long-read processing
Reads were oriented and filtered with Pychopper v2.7.2, aligned to GRCh38 (no alternative contigs, with ERCC sequences) with minimap2, and filtered to MAPQ ≥ 10 with SAMtools v1.6, following Heberle *et al.* All 12 BAM files were merged for a single joint transcript discovery with Bambu, which gave an expanded annotation used throughout, including for the short-read index. Bambu counts were converted to TPM (`tx_counts_to_tpm.r`).

### Differential isoform usage
`isoformswitchanalyzer_DEXSeq.r` runs IsoformSwitchAnalyzeR v2.8.0 with `isoformSwitchTestDEXSeq`. Three surrogate variables estimated with sva during import are included as covariates. Significance: BH FDR < 0.05 and |dIF| > 0.1. Isoforms are called NMD-sensitive when the stop codon lies more than 50 nt upstream of the last exon-exon junction (PTC); isoforms without an assigned ORF are excluded from NMD summaries. The ORF/PTC step (`analyzeORF` on BSgenome.Hsapiens.UCSC.hg38) is the block at the end of `cpm_sensitivity/01_diu_dexseq_cpm.r`.

Sensitivity analyses:
- `isoformswitchanalyzer_satuRn.r`: satuRn v1.16.0 via `isoformSwitchTestSatuRn`, same filter and surrogate variables, all tested genes kept.
- `cpm_sensitivity/`: every analysis repeated with CPM in place of TPM, with the surrogate variables fixed to those of the TPM run.

### Promoter switching
Promoters span 1,000 bp upstream to 100 bp downstream of each transcript TSS. Promoter TPM is the sum of TPM over transcripts starting in the promoter, and promoter usage is promoter TPM divided by gene TPM. Usage is modeled per promoter by beta regression (`usage ~ condition`, logit link) with a Wald test on the condition term and BH adjustment (adjusted p ≤ 0.05). A gene is promoter-switching when it has at least two significant promoters with opposite directions.

### Motif enrichment
Significant promoters are split by direction (logit fold change > 0 is increased usage in AD). Strand-aware ±2 kb windows around each TSS are scanned with PWMEnrich v4.44.0 against the `PWMLogn.hg19.MotifDb.Hsap` log-normal background, and motifs are summarized per TF family with `groupReport`. Enriched TFs are intersected with the AD-associated TF list of Hill and Gammie (2022) to build the TF-promoter network.

> **Fix in this version:** the earlier `circos_heatmap_atac_plots_86_genes_12samples.r` swapped the direction labels, so promoters with increased usage in AD were labelled "Down in AD". The current script labels them "Up in AD", which matches the manuscript.

### Event-level splicing
SUPPA2 v2.3 defines SE, MX, RI, A5, A3, AF and AL events from the Bambu annotation, and PSI is computed from transcript TPM. Testing is restricted to events linked to transcripts of the 86 DIU genes. PSI is squeezed into (0, 1) as (ψ(n − 1) + 0.5)/n and modeled by beta regression with female as the sex reference. Events with fewer than three non-missing values in either group are skipped.
- Additive model (`diff_beta_reg_sig_diu_tx.r`): `PSI ~ condition + sex` vs `PSI ~ sex`, likelihood ratio test with 1 df.
- Condition-by-sex model (`betareg_sex.r`): `PSI ~ condition * sex` vs `PSI ~ sex`, likelihood ratio test with 2 df.

Significance: BH FDR < 0.05 and |dPSI| ≥ 0.10.

### Gene- and transcript-level DE
DESeq2 runs on Bambu gene and transcript counts with design `~ condition`. Features need ≥ 10 counts in ≥ 3 samples of either group. Log2 fold changes are apeglm-shrunken, and features are DE at FDR < 0.05 and |log2FC| > 1. The DESeq2 result tables are produced upstream of `gene_tx_DE/`. The scripts in that folder compare them with the 86 DIU genes by Ensembl ID, using a two-sided Fisher's exact test against the genes tested at each level.

### Short-read replication (`E-MTAB-14167/`)
`salmon_bambu_index_quant.sh` builds a Salmon index from the Bambu extended annotation (gffread transcript FASTA, k = 31). Its quantification step covers only the 14-donor subset that `02_targeted_replication.R` uses as a consistency check. `01_salmon_bambu_quant_emtab14167.sh` quantifies all 41 PFC samples with Salmon v1.10.3 (`-l A --gcBias --validateMappings`); median mapping rate is 89.4%.

- **DIU** (`03_diu_bambu_annotation_emtab14167.R`): uses the same import, filtering and DEXSeq settings as the long-read run, with age and sex added to the design alongside the surrogate variables. A long-read switch counts as replicated when the same isoform is significant (FDR < 0.05, |dIF| > 0.1) in the same direction. Gene overlap is tested with a one-sided Fisher's exact test, and dIF agreement with Spearman correlation.
- **Promoters and PSI events** (`02_targeted_replication.R`): the 46 promoters and 14 PSI events are pre-specified targets. They are recomputed from short-read TPM, tested by beta regression with `~ condition + age + sex`, and BH-corrected within each target class. Directional agreement is tested against 50% with a one-sided binomial test. The script also runs sensitivity sets that the manuscript does not report (APOE ε4, APOE3/3 only, donors outside the 7 vs 7 subset, neuronal fraction from `deconvolution/01_deconvolution.R`).

---

## Usage

> ⚠️ **Paths:** Every script uses hard-coded absolute paths from the analysis server (`/home/AD.UNLV.EDU/Shared_Data/...`). Edit the path variables at the top of each script before running.

### Long-read (`syn52047893/`)

```bash
Rscript tx_counts_to_tpm.r
Rscript pca_tpm.r

# DIU and sensitivity checks
Rscript isoformswitchanalyzer_DEXSeq.r
Rscript isoformswitchanalyzer_satuRn.r
Rscript generate_all_significant_switch_plots_ORF_annotated.R
Rscript volcano_isoform_switches.r
Rscript bar_plot_diu_switch_features.r

# Promoter switching and motif enrichment
Rscript promoter_usage_analysis_86_genes_12samples.r
Rscript circos_heatmap_atac_plots_86_genes_12samples.r

# Splicing graph (junctions from the MAPQ >= 10 BAMs)
bash run_regtools.sh
Rscript splicing_graph_tx_counts.r

# Event-level splicing
bash suppa_generate_events_bambu_gtf.sh
Rscript diff_beta_reg_sig_diu_tx.r
Rscript betareg_sex.r

# CPM sensitivity analysis
Rscript cpm_sensitivity/tx_counts_to_cpm.r   # also called by suppa_psi_cpm.sh
bash cpm_sensitivity/suppa_psi_cpm.sh
bash cpm_sensitivity/run_cpm_pipeline.sh --with-diu
Rscript cpm_sensitivity/07_supplementary_table_tpm_vs_cpm.r

# DIU vs gene/transcript DE (needs the DESeq2 result tables)
Rscript gene_tx_DE/09_supp_figure_deg_txdeg_diu.R      # sources 08_*
Rscript gene_tx_DE/10_supp_table9_gene_tx_DE.R
python3 gene_tx_DE/10b_supp_table9_xlsx.py             # needs openpyxl
```

### Short-read replication (`E-MTAB-14167/`)

```bash
bash salmon_bambu_index_quant.sh
bash 01_salmon_bambu_quant_emtab14167.sh
Rscript deconvolution/01_deconvolution.R    # required input for 02
Rscript 03_diu_bambu_annotation_emtab14167.R
Rscript 02_targeted_replication.R
Rscript 04_replicated_feature_plots.R
Rscript 05_diu_replication_summary_figure.R
```

---

## Software

| Tool | Version | Use |
|---|---|---|
| Pychopper | 2.7.2 | Read orientation and filtering |
| minimap2 | — | Long-read alignment |
| SAMtools | 1.6 | MAPQ filtering |
| Bambu | 1.0.0 (as reported in the manuscript) | Joint transcript discovery and quantification |
| IsoformSwitchAnalyzeR | 2.8.0 | DIU (DEXSeq, satuRn), ORF/PTC |
| satuRn | 1.16.0 | DIU sensitivity analysis |
| SUPPA2 | 2.3 | Splicing events and PSI |
| RegTools | 1.0.0 | Splice junction counts |
| GenomicRanges | 1.60.0 | Promoter intervals |
| BSgenome.Hsapiens.UCSC.hg38 | 1.4.5 | Promoter sequences, ORFs |
| PWMEnrich | 4.44.0 | Motif enrichment |
| Salmon | 1.10.3 | Short-read quantification |
| DESeq2 + apeglm | — | Gene/transcript DE |

R packages (CRAN): `betareg`, `broom`, `data.table`, `dplyr`, `tidyr`, `tidyverse`, `stringr`, `forcats`, `ggplot2`, `ggrepel`, `patchwork`, `gridExtra`, `circlize`, `matrixStats`, `nnls`, `Matrix`, `gprofiler2`.
Bioconductor: `IsoformSwitchAnalyzeR` (with its `DEXSeq`, `satuRn` and `sva` backends), `DESeq2`, `limma`, `tximport`, `GenomicRanges`, `GenomicFeatures`, `rtracklayer`, `BSgenome.Hsapiens.UCSC.hg38`, `PWMEnrich`, `PWMEnrich.Hsapiens.background`, `MotifDb`, `TFBSTools`, `SplicingGraphs`, `rhdf5`, `dtangle`.
Python: `openpyxl` (for `10b_supp_table9_xlsx.py`).

Reference: GRCh38 with Ensembl release 107 annotation plus ERCC spike-ins (`Homo_sapiens.GRCh38.107_ERCC.gtf`).

---

## Contact

For questions about this pipeline, please open a [GitHub issue](https://github.com/jelardaquino/TSSwitch/issues).
