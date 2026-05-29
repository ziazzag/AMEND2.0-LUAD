# Active Module Identification in Lung Adenocarcinoma using AMEND 2.0

**INFO-F439 : Vrije Universiteit Brussel / Université Libre de Bruxelles**  
**Author:** Zachary Iazzag

Reimplementation of the AMEND 2.0 pipeline (Boyd et al., 2025) applied to TCGA lung adenocarcinoma (LUAD) data instead of the original kidney clear cell carcinoma (KIRC) dataset.

---

## Overview

AMEND (Active Module Identification using Experimental data and Network Diffusion) identifies a connected subgraph in a multiplex-heterogeneous molecular interaction network whose nodes have large experimental signal. It alternates between a diffusion step (RWR-MH) and a subgraph optimisation step (MWCS via heinz) until convergence.

This study applies AMEND 2.0 to 572 TCGA-LUAD samples (515 tumour, 57 matched normal) across three omic layers : mRNA expression, miRNA expression, and DNA methylation, integrated into a seven-layer network.

---

## Repository Structure

```
.
├── Code/
│   ├── TCGA_LUAD_Preprocessing.R   # Data loading, DE analysis, survival modelling
│   ├── TCGA_LUAD_networks.R        #Network construction (7 layers)
│   ├── TCGA_LUAD_amend_run.R       #AMEND runs (3 configurations)
│   └── TCGA_LUAD_analysis.R        #GO/DO/DGN enrichment, figures, Jaccard
├── Data/
│   └── LUAD/                       #Raw TCGA data 
    └── Interaction Networks/       #interaction networks
├── Results/                        #Generated outputs (CSV, figures)
```

---

## Data

 Downloaded from the **Broad GDAC Firehose** archive:

- URL: http://gdac.broadinstitute.org/runs/stddata__2016_01_28/data/LUAD/20160128/
- Run: `stddata__2016_01_28`

Files required (`Data/LUAD/`):

| File | Content |
|------|---------|
| `LUAD.uncv2.mRNAseq_RSEM_normalized_log2.txt` | mRNA expression (RSEM, log2) |
| `LUAD.miRseq_mature_RPM_log2.txt` | miRNA expression (RPM, log2) |
| `LUAD.meth.by_mean.data.txt` | DNA methylation (beta values) |
| `All_CDEs.txt` | Clinical data |

Interaction networks required ( `Data/Interaction Networks/`):

| File | Source |
|------|--------|
| `9606.protein.physical.links.v12.0.txt.gz` | [STRING v12](https://string-db.org) |
| `PPT-Ohmnet_tissues-combined.txt` | [OhmNet](http://snap.stanford.edu/ohmnet/) |
| `miRTarBase_SE_WR.xlsx` | [miRTarBase](https://mirtarbase.cuhk.edu.cn) |

The ARACNE gene regulatory network (`aracne_el_luad.txt`) is generated automatically from the `aracne.networks` R package (`regulonluad` object).

---

## Installation

**Critical version constraints:**
```r
# igraph 2.x breaks AMEND — must use 1.6.0
remotes::install_version("igraph", version = "1.6.0")

# AMEND from GitHub
remotes::install_github("samboyd0/AMEND")
```

**Full R package list:**

```r
# CRAN
install.packages(c("data.table", "dplyr", "ggplot2", "Matrix",
                   "survival", "missForest", "limma", "readxl", "remotes"))

# Bioconductor
BiocManager::install(c("fgsea", "DOSE", "org.Hs.eg.db", "AnnotationDbi",
                       "GO.db", "graphite", "aracne.networks"))
```



---

## Usage

Run the four scripts in order from the project root:

```r
source("Code/TCGA_LUAD_Preprocessing.R")
source("Code/TCGA_LUAD_networks.R")
source("Code/TCGA_LUAD_amend_run.R")
source("Code/TCGA_LUAD_analysis.R")
```

All outputs are written to `Results/`

---

## Run Configurations

| | Run 1 | Run 2 | Run 3 |
|---|---|---|---|
| Biased Random Walk | Yes | No | Yes |
| IN degree bias | No | No | Yes |
| Multiplex aggregation | No | No | Yes |
| Module size | 47 nodes | 56 nodes | 59 nodes |

Run 3 reproduces the original KIRC setup from Boyd et al. (2025).

---

## Key Results

The Run 3 module (59 nodes: 57 mRNA + 2 miRNA) is enriched in mitotic spindle checkpoint and chromosomal segregation processes. Key hub genes: **MYBL2** (degree 16), **TOP2A** (15), **BUB1B** (11), **NDC80** (11), **PLK1** (9).
---
