# TSSwitch

**Isoform Switching and Differential Transcript Usage Analysis in Alzheimer's Disease**

[![R](https://img.shields.io/badge/R-%3E%3D4.2-blue.svg)](https://www.r-project.org/)
[![Bioconductor](https://img.shields.io/badge/Bioconductor-%3E%3D3.16-green.svg)](https://bioconductor.org/)

---

## Overview

TSSwitch (Transcriptomic Splicing Switch) is a bioinformatics analysis pipeline designed to investigate transcriptomic dysregulation and alternative splicing alterations in Alzheimer’s Disease (AD) using both short-read (Illumina) and long-read (Oxford Nanopore Technologies, ONT) RNA sequencing data from human (Homo sapiens, GRCh38) brain samples.

The pipeline integrates multiple complementary analyses to characterize changes in RNA processing and transcript regulation between AD and healthy control (CT) samples, including:

Differential transcript usage (DTU)
Isoform switching analysis
Event-level alternative splicing detection
Differential promoter/transcription start site usage
Transcription factor (TF) and motif enrichment
Cross-platform transcriptomic comparisons between short- and long-read sequencing technologies

The name TSSwitch reflects the central theme of the project: identifying transcriptomic and splicing “switches” associated with Alzheimer’s disease pathology. By combining long-read transcript reconstruction with short-read statistical frameworks, TSSwitch enables the discovery of both known and novel transcript isoforms, promoter shifts, and splicing events that may contribute to neurodegeneration and transcriptional reprogramming in AD.

This pipeline provides a unified framework for studying isoform diversity, transcript regulation, and alternative splicing dynamics across sequencing platforms.

Key analyses include:
- DTU testing with [IsoformSwitchAnalyzeR](https://bioconductor.org/packages/release/bioc/html/IsoformSwitchAnalyzeR.html)
- Alternative splicing event quantification (SUPPA2)
- Differential splicing analysis (beta regression)
- Differential promoter usage (beta regression)
- Motif enrichment at switching promoters (PWMEnrich)
- GO/pathway enrichment (g:Profiler2, clusterProfiler)
- Publication-quality visualization (volcano plots, switch plots, PCA, circos/heatmaps)

---

## Repository Structure

```
TSSwitch/
├── PRJNA1206164/               # Short-read Illumina RNA-seq pipeline (public dataset)
│   ├── fastqc.sh                   # Raw read quality control
│   ├── trimmomatic.sh               # Adapter trimming and quality filtering
│   ├── salmon_quant.sh              # Transcript quantification with Salmon
│   ├── salmon_tx_tpm_counts.r       # TPM matrix preparation and PCA from Salmon output
│   ├── isoformswitchanalyzer_satuRn.r  # DTU testing with satuRn via IsoformSwitchAnalyzeR
│   ├── generate_all_significant_switch_plots_ORF_annotated.R  # Batch switch plot generation
│   ├── custom_switch_plot_helpers.R # Reusable ggplot2 helper functions for switch plots
│   ├── promoter_usage_analysis.r    # Differential promoter usage with beta regression
│   └── circos_heatmap_atac_plots.r  # Circos plots, heatmaps, and PWMEnrich motif analysis
│
└── syn52047893/                # Long-read Oxford Nanopore RNA-seq pipeline (Synapse dataset)
    ├── tx_counts_to_tpm.r           # TPM conversion from Bambu count output
    ├── pca_tpm.r                    # PCA of TPM expression matrix
    ├── isoformswitchanalyzer_DEXSeq.r  # DTU testing with DEXSeq via IsoformSwitchAnalyzeR
    ├── generate_all_significant_switch_plots_ORF_annotated.R  # Batch switch plot generation
    ├── custom_switch_plot_helpers.R # Reusable ggplot2 helper functions for switch plots
    ├── example_custom_plots.R       # Example usage of custom_switch_plot_helpers.R
    ├── volcano_isoform_switches.r   # Volcano plots of DTU results and splicing summaries
    ├── bar_plot_diu_switch_features.r # Bar plots summarizing switching feature counts
    ├── promoter_usage_analysis_86_genes_12samples.r  # Promoter usage analysis (long-read)
    ├── suppa_generate_events_bambu_gtf.sh  # SUPPA2 event generation and PSI calculation
    ├── beta_regression_psi.r        # Beta regression on PSI values for differential splicing
    ├── gprofiler2_dotplot.r         # GO enrichment dot plots with g:Profiler2
    ├── splicing_graph_tx_counts.r   # Splicing graph construction and transcript count visualization
    └── circos_heatmap_atac_plots_86_genes_12samples.r  # Circos/heatmaps for long-read dataset
```

---

## Datasets

| Folder | Accession | Description | Sequencing | Samples |
|---|---|---|---|---|
| `PRJNA1206164/` | [PRJNA1206164](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1206164) | Short-read RNA-seq (public, NCBI SRA) | Illumina paired-end | 7 AD vs. 8 CT |
| `syn52047893/` | [syn52065646](https://www.synapse.org/#!Synapse:syn52065646) | Long-read RNA-seq (AD Knowledge Portal, controlled access) | Oxford Nanopore | 6 AD, 6 CT |

> ⚠️ **Data Access:** The long-read dataset (`syn52047893`) requires an approved Synapse account and data use agreement through the [AD Knowledge Portal](https://adknowledgeportal.synapse.org/). But you can also access it at the [Sequence Read Archive](https://www.ncbi.nlm.nih.gov/sra) under SRP456327. Raw data files are **not** included in this repository.

---

## Dependencies

### System Requirements
- Linux or macOS
- R ≥ 4.2
- Bioconductor ≥ 3.16
- Python ≥ 3.8 (for SUPPA2)
- [GNU Parallel](https://www.gnu.org/software/parallel/) (for shell scripts)

### Short-Read Pipeline Tools
| Tool | Version | Purpose |
|---|---|---|
| [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) | ≥ 0.11 | Raw read QC |
| [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) | ≥ 0.39 | Adapter trimming |
| [Salmon](https://combine-lab.github.io/salmon/) | ≥ 1.9 | Transcript quantification |

### Long-Read Pipeline Tools
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
  "DEXSeq",                 # DTU statistical test (used in long-read)
  "satuRn",                 # DTU statistical test (used in short-read)
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
  "clusterProfiler",        # GO/KEGG enrichment
  "biomaRt"                 # Ensembl annotation queries
))

# CRAN
install.packages(c(
  "ggplot2",      # Visualization
  "ggrepel",      # Non-overlapping text labels
  "dplyr",        # Data manipulation
  "tidyr",        # Data reshaping
  "data.table",   # High-performance tabular data
  "betareg",      # Beta regression for PSI/IF values
  "gprofiler2",   # g:Profiler GO enrichment
  "patchwork",    # Multi-panel figure layout
  "gridExtra",    # Grid-based figure layout
  "matrixStats",  # Row/column statistics
  "glmmTMB",      # Mixed-effects beta regression
  "broom.mixed",  # Tidy model outputs
  "forcats",      # Factor manipulation
  "stringr"       # String manipulation
))
```

---

## Usage

> ⚠️ **Path Configuration:** All scripts currently use hard-coded absolute paths pointing to the original analysis server (e.g., `/home/AD.UNLV.EDU/Shared_Data/...`). Before running, update all path variables at the top of each script to match your local environment. A shared `config.R` approach is recommended for multi-script workflows.

### Short-Read Pipeline (`PRJNA1206164/`)

Run scripts in the following order from your working data directory:

```bash
# 1. Quality control
bash fastqc.sh

# 2. Trimming
bash trimmomatic.sh

# 3. Transcript quantification (requires Salmon index)
bash salmon_quant.sh
```

```r
# 4. Build TPM matrix and run PCA
Rscript salmon_tx_tpm_counts.r

# 5. DTU testing with satuRn
Rscript isoformswitchanalyzer_satuRn.r

# 6. Generate switch plots for significant isoforms
Rscript generate_all_significant_switch_plots_ORF_annotated.R

# 7. Differential promoter usage
Rscript promoter_usage_analysis.r

# 8. Circos plots, heatmaps, and motif enrichment
Rscript circos_heatmap_atac_plots.r
```

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

# 7. Differential promoter usage
Rscript promoter_usage_analysis_86_genes_12samples.r
```

```bash
# 8. SUPPA2 event generation and PSI calculation
bash suppa_generate_events_bambu_gtf.sh
```

```r
# 9. Beta regression on PSI values
Rscript beta_regression_psi.r

# 10. GO enrichment dot plots
Rscript gprofiler2_dotplot.r

# 11. Circos plots, heatmaps, motif enrichment
Rscript circos_heatmap_atac_plots_86_genes_12samples.r
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
| `transcript_tpm.txt` / `bambu_tx_tpm.txt` | Transcript-level TPM matrices |
| `promoter_usage_betareg_results.tsv` | Beta regression results for differential promoter usage |
| `suppa_beta_regression_all_results.tsv` | Beta regression results for SUPPA2 PSI events |
| `gprofiler_results.tsv` | GO enrichment results table |
| `figures/` | All publication-quality PNG figures (300 DPI) |

---

## Contact

For questions about this pipeline, please open a [GitHub Issue](https://github.com/jelardaquino/TSSwitch/issues).
