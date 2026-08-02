# TSSwitch

**Transcriptomic Splicing Switch**

---

## Overview

TSSwitch is a bioinformatics analysis pipeline designed to investigate transcriptomic dysregulation and alternative splicing alterations in Alzheimer's Disease (AD) using long-read (Oxford Nanopore Technologies, ONT) RNA sequencing data from human (Homo sapiens, GRCh38) brain samples, with independent cross-platform validation using short-read (Illumina) RNA-seq.

The pipeline integrates multiple complementary analyses to characterize changes in RNA processing and transcript regulation between AD and healthy control (CT) samples, including:

- Differential transcript usage (DTU) via long-read sequencing
- Isoform switching analysis
- Targeted, gene-set-restricted alternative splicing (SUPPA2 + beta regression) on genes with significant DTU
- Differential promoter/transcription start site usage
- Transcription factor (TF) and motif enrichment
- Independent short-read RNA-seq concordance analysis of the significant DTU gene set

The name TSSwitch reflects the central theme of the project: identifying transcriptomic and splicing "switches" associated with Alzheimer's disease pathology. By combining long-read transcript reconstruction with a targeted, hypothesis-driven splicing analysis and an independent short-read concordance check, TSSwitch enables discovery of known and novel transcript isoforms, promoter shifts, and splicing events that may contribute to neurodegeneration and transcriptional reprogramming in AD — while directly assessing how reproducible these findings are across sequencing platforms.

Key analyses include:
- DTU testing with [IsoformSwitchAnalyzeR](https://bioconductor.org/packages/release/bioc/html/IsoformSwitchAnalyzeR.html) (DEXSeq backend, long-read)
- Alternative splicing event quantification (SUPPA2)
- **Targeted** differential splicing analysis (beta regression) restricted to transcripts belonging to genes with significant differential isoform usage (DIU)
- Differential promoter usage (beta regression)
- Motif enrichment at switching promoters (PWMEnrich)
- GO/pathway enrichment (g:Profiler2)
- Independent short-read concordance analysis of AD-vs-CT effect sizes for the significant DIU gene set
- Publication-quality visualization (volcano plots, switch plots, PCA, circos/heatmaps)

---

## Repository Structure

```
TSSwitch/
├── syn52047893/                 # Long-read Oxford Nanopore RNA-seq pipeline (Synapse dataset)
│   ├── tx_counts_to_tpm.r           # TPM conversion from Bambu count output
│   ├── pca_tpm.r                    # PCA of TPM expression matrix
│   ├── isoformswitchanalyzer_DEXSeq.r  # DTU testing with DEXSeq via IsoformSwitchAnalyzeR
│   ├── generate_all_significant_switch_plots_ORF_annotated.R  # Batch switch plot generation
│   ├── custom_switch_plot_helpers.R # Reusable ggplot2 helper functions for switch plots
│   ├── example_custom_plots.R       # Example usage of custom_switch_plot_helpers.R
│   ├── volcano_isoform_switches.r   # Volcano plots of DTU results and splicing summaries
│   ├── bar_plot_diu_switch_features.r # Bar plots summarizing switching feature counts
│   ├── promoter_usage_analysis_86_genes_12samples.r  # Promoter usage analysis restricted to the 86 significant DIU genes
│   ├── suppa_generate_events_bambu_gtf.sh  # SUPPA2 event generation and PSI calculation
│   ├── diff_beta_reg_sig_diu_tx.r   # Targeted beta regression: differential splicing (PSI ~ condition + sex) restricted to transcripts of the 86 significant DIU genes
│   ├── betareg_sex.r                # Same targeted approach, testing a condition × sex interaction
│   ├── gprofiler2_dotplot.r         # GO enrichment dot plots with g:Profiler2 on targeted splicing results
│   ├── splicing_graph_tx_counts.r   # Splicing graph construction and transcript count visualization
│   └── circos_heatmap_atac_plots_86_genes_12samples.r  # Circos/heatmaps for the 86 significant DIU genes
│
└── ERP161086/                   # Independent short-read RNA-seq concordance analysis
    ├── fastqc.sh                    # Raw read quality control
    ├── salmon_quant.sh               # Transcript quantification with Salmon (v0.14.2)
    ├── salmon_tx_tpm_counts.r        # TPM matrix preparation from Salmon output
    └── comparig_86_sig_genes_LR_SR_counts.R  # Long-read vs. short-read concordance analysis for the 86 significant DIU genes
```

---

## Datasets

| Folder | Accession | Description | Sequencing | Samples |
|---|---|---|---|---|
| `syn52047893/` | [syn52065646](https://www.synapse.org/#!Synapse:syn52065646) | Long-read RNA-seq (AD Knowledge Portal, controlled access) — primary DTU/DIU discovery dataset | Oxford Nanopore | 6 AD, 6 CT |
| `ERP161086/` | [ERP161086](https://www.ebi.ac.uk/ena/browser/view/ERP161086) | Short-read RNA-seq (European Nucleotide Archive) — independent postmortem prefrontal cortex cohort used for cross-platform concordance | Illumina paired-end, 150 bp | 7 AD, 7 CT |

> ⚠️ **Data Access:** The long-read dataset (`syn52047893`) requires an approved Synapse account and data use agreement through the [AD Knowledge Portal](https://adknowledgeportal.synapse.org/). It can also be accessed at the [Sequence Read Archive](https://www.ncbi.nlm.nih.gov/sra) under SRP456327. Raw data files are **not** included in this repository.

---

## Independent Short-Read RNA-seq Concordance Analysis (`ERP161086`)

### Methods

Transcript-level TPM estimates from the long-read dataset were compared with an independent short-read bulk RNA-seq dataset from postmortem prefrontal cortex (European Nucleotide Archive: ERP161086), consisting of 150-bp Illumina paired-end reads from AD (n = 7) and CT (n = 7) individuals. Reads were aligned to the GRCh38 reference genome, and TPM was quantified using Salmon v0.14.2.

Analyses were restricted to transcripts assigned to the significant DIU gene set (86 genes identified from the long-read DTU analysis) and present in both the long-read and short-read TPM matrices. TPM values were log-transformed as $\log_2(\text{TPM} + 1)$. For each shared transcript, mean expression was calculated separately in AD and CT samples, and an AD-vs-CT effect size was computed as:

$$\Delta_{AD-CT} = \text{mean}\ \log_2(\text{TPM}+1)_{AD} - \text{mean}\ \log_2(\text{TPM}+1)_{CT}$$

This produced a long-read effect size and a short-read effect size for each transcript. Transcript-level concordance between platforms was evaluated by comparing the direction and magnitude of these AD-vs-CT effects.

- **Directional concordance:** Transcripts were considered directionally concordant if the long-read and short-read effect sizes had the same sign, indicating that both platforms showed higher expression in AD or lower expression in AD relative to CT.
- **Large effect prioritization:** Transcripts were prioritized if they showed a large AD-vs-CT difference in both platforms, defined as an absolute log2 effect size above a selected threshold in both long- and short-read data.
- **Minimum expression filter:** To avoid prioritizing transcripts with low or absent expression, transcripts were also required to meet a minimum TPM expression threshold in at least one platform.
- **Effect-size ratio:** Similarity in effect magnitude was quantified using an effect-size ratio, calculated as the smaller absolute effect size divided by the larger absolute effect size:

$$\text{effect size ratio} = \frac{\min(|\Delta_{LR}|,\ |\Delta_{SR}|)}{\max(|\Delta_{LR}|,\ |\Delta_{SR}|)}$$

Values closer to 1 indicate more similar AD-vs-CT effect magnitudes between platforms.

- **Final classification:** Transcripts were classified as having large, similar AD-vs-CT differences if they were expressed above the minimum TPM threshold, showed large absolute AD-vs-CT effects in both platforms, changed in the same direction, and exceeded the minimum effect-size ratio threshold.

### Script

`ERP161086/comparig_86_sig_genes_LR_SR_counts.R` implements this workflow: it loads the long-read and short-read TPM matrices and metadata, restricts to the 86 significant DIU genes/transcripts shared between platforms, computes per-condition log2(TPM+1) means, AD-vs-CT effect sizes, directional concordance, and the effect-size ratio, and writes summary tables and scatter plots comparing long-read vs. short-read effect sizes.

---

## Targeted Differential Splicing Analysis (`syn52047893/diff_beta_reg_sig_diu_tx.r`, `betareg_sex.r`)

Rather than testing alternative splicing events genome-wide, the current pipeline takes a **targeted approach**: it restricts SUPPA2 event testing to transcripts belonging to the **86 genes with significant differential isoform usage (DIU)** identified by `isoformswitchanalyzer_DEXSeq.r`. This directly asks whether genes already implicated by DTU/DIU also show detectable differential alternative splicing (skipped exon, retained intron, mutually exclusive exon, alternative 5'/3' splice site, and alternative first/last exon events).

Workflow:
1. Read the list of significant DIU transcripts (`sig_transcripts.txt`) and match them to SUPPA2 `.ioe` event annotations for each event type (SE, RI, MX, AL, AF, A5, A3).
2. Run event-level beta regression of PSI ~ condition (+ sex, or condition × sex in `betareg_sex.r`) only for events linked to the significant DIU transcripts.
3. Filter results by BH-adjusted p-value (≤ 0.05) and effect size (|ΔPSI| ≥ 0.10).
4. Generate compact, human-readable plot labels from SUPPA2 event coordinates (e.g., gene | event type | exon/intron coordinates) and PSI-by-condition/sex figures for significant events.
5. Write per-event-type and combined summary tables (`all_beta_results.tsv`, `all_sig_beta_events.tsv`, `summary_by_event_type.tsv`) for downstream enrichment via `gprofiler2_dotplot.r`.

This replaces the previous whole-genome-wide beta regression on all detected SUPPA2 events, which has been removed from the repository.

---

## Dependencies

### System Requirements
- Linux or macOS
- R ≥ 4.2
- Bioconductor ≥ 3.16
- Python ≥ 3.8
- [GNU Parallel](https://www.gnu.org/software/parallel/) (for shell scripts)

### Short-Read Concordance Pipeline Tools (`ERP161086/`)
| Tool | Version | Purpose |
|---|---|---|
| [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) | ≥ 0.11 | Raw read QC |
| [Salmon](https://combine-lab.github.io/salmon/) | 0.14.2 | Transcript quantification |

### Long-Read Pipeline Tools (`syn52047893/`)
| Tool | Version | Purpose |
|---|---|---|
| [Minimap2](https://github.com/lh3/minimap2) | ≥ 2.24 | Long-read alignment |
| [Bambu](https://github.com/GoekeLab/bambu) (R) | ≥ 3.0 | Transcript discovery and quantification |
| [SUPPA2](https://github.com/comprna/SUPPA) | ≥ 2.3 | Alternative splicing event quantification |

### R Packages
```r
# Bioconductor
BiocManager::install(c(
  "IsoformSwitchAnalyzeR",  # Core DTU analysis framework
  "DEXSeq",                 # DTU statistical test (long-read)
  "bambu",                  # Long-read quantification
  "tximport",               # Salmon output import
  "GenomicRanges",          # Genomic interval operations
  "GenomicFeatures",        # TxDb construction
  "rtracklayer",            # GTF/GFF import/export
  "BSgenome.Hsapiens.UCSC.hg38",  # Genome sequence for ORF analysis
  "PWMEnrich",              # Transcription factor motif enrichment
  "PWMEnrich.Hsapiens.background",
  "SplicingGraphs",         # Splicing graph visualization
  "org.Hs.eg.db",           # Human gene ID mapping
  "biomaRt"                 # Ensembl annotation queries
))

# CRAN
install.packages(c(
  "ggplot2",      # Visualization
  "ggrepel",      # Non-overlapping text labels
  "ggpattern",    # Patterned fills for concordance figures
  "dplyr",        # Data manipulation
  "tidyr",        # Data reshaping
  "tidyverse",    # Data wrangling / readr / purrr
  "data.table",   # High-performance tabular data
  "betareg",      # Beta regression for PSI/IF values
  "broom",        # Tidy model outputs (betareg)
  "gprofiler2",   # g:Profiler GO enrichment
  "patchwork",    # Multi-panel figure layout
  "gridExtra",    # Grid-based figure layout
  "matrixStats",  # Row/column statistics
  "forcats",      # Factor manipulation
  "stringr"       # String manipulation
))
```

---

## Usage

> ⚠️ **Path Configuration:** All scripts currently use hard-coded absolute paths pointing to the original analysis server (e.g., `/home/AD.UNLV.EDU/Shared_Data/...`). Before running, update all path variables at the top of each script (`BASE_DIR`, `SUPPA_DIR`, `OUTPUT_DIR`, etc.) to match your local environment. A shared `config.R` approach is recommended for multi-script workflows.

### Long-Read Pipeline (`syn52047893/`)

```r
# 1. Convert Bambu counts to TPM
Rscript tx_counts_to_tpm.r

# 2. PCA of expression
Rscript pca_tpm.r

# 3. DTU testing with DEXSeq
Rscript isoformswitchanalyzer_DEXSeq.r

# 4. Generate switch plots
Rscript generate_all_significant_switch_plots_ORF_annotated.R

# 5. Volcano plots and splicing summaries
Rscript volcano_isoform_switches.r

# 6. Bar plots of switching feature counts
Rscript bar_plot_diu_switch_features.r

# 7. Differential promoter usage (86 significant DIU genes)
Rscript promoter_usage_analysis_86_genes_12samples.r
```

```bash
# 8. SUPPA2 event generation and PSI calculation
bash suppa_generate_events_bambu_gtf.sh
```

```r
# 9. Targeted beta regression on PSI values, restricted to the 86 significant DIU genes
Rscript diff_beta_reg_sig_diu_tx.r

# 9b. (Optional) Same targeted analysis with a condition x sex interaction term
Rscript betareg_sex.r

# 10. GO enrichment dot plots on the targeted splicing results
Rscript gprofiler2_dotplot.r

# 11. Circos plots, heatmaps, motif enrichment (86 significant DIU genes)
Rscript circos_heatmap_atac_plots_86_genes_12samples.r
```

### Independent Short-Read Concordance Analysis (`ERP161086/`)

```bash
# 1. Quality control
bash fastqc.sh

# 2. Transcript quantification (Salmon v0.14.2, requires Salmon index)
bash salmon_quant.sh
```

```r
# 3. Build TPM matrix
Rscript salmon_tx_tpm_counts.r

# 4. Compare long-read vs. short-read TPM for the 86 significant DIU genes:
#    effect sizes, directional concordance, effect-size ratio
Rscript comparig_86_sig_genes_LR_SR_counts.R
```

---

## Reference Genome

Both pipelines use **Ensembl GRCh38 release 107** with ERCC spike-in sequences:

```
Homo_sapiens.GRCh38.107_ERCC.gtf
```

The Salmon index for short-read quantification was built against Ensembl GRCh38 cDNA sequences.

---

## Key Outputs

| Output | Description |
|---|---|
| `switchList_with_gene_name.rds` | IsoformSwitchAnalyzeR object with DTU results (RDS checkpoint) |
| `IsoformSwitchAnalyzeR_top_switches.tsv` | Table of top significant isoform switches |
| `bambu_tx_tpm.txt` | Long-read transcript-level TPM matrix |
| `sig_transcripts.txt` | Transcript IDs from the 86 genes with significant DIU |
| `promoter_usage_betareg_results.tsv` | Beta regression results for differential promoter usage (86 significant DIU genes) |
| `all_beta_results.tsv` | All targeted PSI beta regression results (transcripts linked to the 86 significant DIU genes) |
| `all_sig_beta_events.tsv` | Significant targeted splicing events (adj. p ≤ 0.05, \|ΔPSI\| ≥ 0.10) |
| `summary_by_event_type.tsv` | Count of tested/significant events by SUPPA2 event type |
| `gprofiler_results.tsv` | GO enrichment results table for the targeted splicing analysis |
| `lr_sr_gene_level_log2fc_similarity.csv` | Per-gene long-read vs. short-read AD-vs-CT effect sizes and directional concordance |
| `lr_sr_log2fc_similarity_summary.csv` | Summary statistics (correlation, % same direction) for the concordance analysis |
| `figures/` | All publication-quality PNG figures (300 DPI) |

---

## Contact

For questions about this pipeline, please open a [GitHub Issue](https://github.com/jelardaquino/TSSwitch/issues).
