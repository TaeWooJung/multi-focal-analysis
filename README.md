# 🧬 Multifocal Prostate Cancer Cohort Analysis Using RNA-Seq Data

This repository contains code, data references, and workflows supporting the manuscript:

**"Multifocal cohort analysis unveils cell types associated with regional lymph node seeding in prostate cancer"**  
Louise de Schaetzen van Brienen, Taewoo Jung, et al.  
📅 Submitted to Genome Medicine, 2026  
🔗 [DOI]

---

## 📜 Abstract

> **(Background)** Understanding the molecular features that underlie metastatic prostate cancer (PCa) is essential to develop prognostic markers and improve treatment decisions. However, such studies are hampered by substantial intratumor and interpatient heterogeneity. 
**(Methods)** To cope with this heterogeneity, we propose a unique study design that leverages the statistical power of multifocal bulk transcriptome profiling with the resolution of a single-cell analysis to identify processes and cell types associated with regional metastatic lymph node seeding in PCa.
**(Results)** Elaborate analysis of these data allowed identifying a metric to distinguish, based on the multifocal expression data between lesions with high and low potential for regional metastatic lymph node seeding. Subsequently comparing the expression profiles of these lesions with respectively high and low metastatic potential identified an aggressiveness signature.  Overlaying this signature with single cell data identified proliferative luminal cells, an adipose derived cancer-associated fibroblast (CAF) state and a specific subtype of arterial endothelial cells. Assessing the prognostic value of these cell states in an independent dataset (TCGA-PRAD) confirmed their association with regional metastatic lymph node seeding and progression free survival and unveiled a complementary role for the proliferative luminal cells and the adipose derived CAF state in driving regional metastatic lymph node seeding. 
**(Conclusion)** Based on our analysis we hypothesize that lesions with high potential for regional metastatic lymph node seeding are mostly characterized by the presence of highly proliferative luminal cells and a transitioning towards an aggressive adipose derived CAF state.  Markers associated with these cell states largely explain the prognostic signal of currently used commercial signatures in PCa, further supporting the role of the identified cell states in driving regional metastatic lymph node seeding and providing an in depth understanding of the success of the currently used commercial signatures.


---

## 📁 Repository Contents

```bash
.
├── Data/                  # Processed data or links to raw datasets
├── Scripts/               # RNA-seq and single-cell analysis scripts (preprocessing, DEGs, clustering)
├── Results/               # Key figures and tables
├── Figures/               # Key files to generate figures in the paper
├── FigShare/              # Key figures and tables
├── SupplementaryTables/   # Key tables
└── README.md              # This file
```
--

## 📦 Setup Instructions
1. Clone the Repository
```bash
git clone https://github.com/TaeWooJung/multi-focal-analysis.git
cd multi-focal-analysis
```
2. Prerequisites
- conda
- R

--

## 🛠️ Run Analysis
To reproduce the full pipeline:

### Variant calling and seeding lesion identification
```bash
# Perform VarScan2 on tNGS data from locally advanced cohort
Scripts/tNGS_variant_calling/0_Varscan_all.sh                               # Perform VarScan2     
Scripts/tNGS_variant_calling/1_locally_advanced_varscan_variants.ipynb      # Pre-processing VarScan2 output
Scripts/tNGS_variant_calling/2_VarScanFiltering.ipynb                       # Filtering VarScan2 output
# Variant summary of locally advanced and de novo cohorts
Scripts/tNGS_variant_calling/3_tNGS_variant_call_analysis.ipynb             
# Main analysis
Scripts/1_bulk_seq_analysis.Rmd
```

### Singe cell preprocessing and annotations
```bash
# Preprocessing raw data from single cell dataset from Chen et al. (2021)
Scripts/single_cell_preprocessing/0_chenRawData.R
# Run inferCNV to classify malignant cells
Scripts/single_cell_preprocessing/1_inferCNV_analysis.R
```

### Aggressive cell identification and downstream analysis
```bash
Scripts/2_sc_analysis_luminal.Rmd            # Luminal cells
Scripts/3_sc_analysis_fibroblast.Rmd         # Fibroblast
Scripts/4_sc_analysis_endothelial.Rmd        # Endothelial cells
Scripts/5_sc_analysis_plots.Rmd              # Generating heatmaps
```

### Survival analysis on TCGA-PRAD
```bash
Scripts/6_survival_analysis.Rmd
```

--

## 🧪 Data Access
### Raw Data
RNA-seq and tNGS datasets from Multi-focal cohorts are available on European Genome-Phenome Archive (EGA):
- EGAS00001006466 [2](#ref2)
- EGAS00001006715 [3](#ref3)
TCGA-PRAD dataset from National Cancer Institute and GDC Data Portal:
- Bulk RNA-seq [4](#ref4)
- Genome annotation [5](#ref5)
- Clinical data [6](#ref6)

--

## 📘 References
1. <a id="ref1"></a> Chen S, Zhu G, Yang Y, Wang F, Xiao YT, Zhang N et al. Single-cell analysis reveals transcriptomic remodellings in distinct cell types that contribute to human prostate cancer progression. Nat Cell Biol 2021; 23:87-98. https://doi.org/10.1038/s41556-020-00613-6.
2. <a id="ref2"></a> Wyatt A. Sequencing data for the manuscript "Multi-focal sampling of de novo metastatic prostate cancer reveals complex polyclonality and enables accurate clinical genotyping". European Genome-Phenome Archive. https://ega-archive.org/studies/EGAS00001006466 (2022).
3. <a id="ref3"></a> Marchal K. Sequencing data for the manuscript "Multifocal cohort analysis unveils cell types associated with regional lymph node seeding in prostate cancer". European Genome-Phenome Archive. https://ega-archive.org/studies/EGAS00001006715 (2026).
4. <a id="ref4"></a> National Cancer Institute. TCGA-PRAD bulk RNA-seq dataset. Xenahub Portal. https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-PRAD.htseq_counts.tsv.gz (2023).
5. <a id="ref5"></a> Luo Y. xena-GDC-ETL. GitHub. https://github.com/ucscXena/xena-GDC-ETL/blob/master/xena_gdc_etl/resources/gencode.v22.annotation.gene.probeMap (2019).
6. <a id="ref6"></a> GDC Data Portal. TCGA Prostate Adenocarcinoma. https://portal.gdc.cancer.gov/projects/TCGA-PRAD. Accessed 17 Nov 2021.