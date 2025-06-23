# 🧬 Multifocal Prostate Cancer Cohort Analysis Using RNA-Seq Data

This repository contains code, data references, and workflows supporting the manuscript:

**"Multifocal cohort analysis unveils cell types associated with regional lymph node seeding in prostate cancer"**  
Louise de Schaetzen van Brienen, Taewoo Jung, et al.  
📅 Submitted to [Journal Name], 2025  
🔗 [DOI / preprint link]

---

## 📜 Abstract

> **(Background)** Understanding the molecular features that underlie metastatic prostate cancer is essential to develop prognostic markers and improve treatment decisions. However, such studies are hampered by substantial intratumor and interpatient heterogeneity. 
**(Methods)** To cope with this heterogeneity, we propose a unique study design that leverages the statistical power of multifocal bulk transcriptome profiling with the resolution of a single-cell analysis to identify processes and cell types associated with regional metastatic lymph node seeding in PCa.
**(Results)** Elaborate analysis of these data allowed identifying a metric to distinguish, based on the multifocal expression data between lesions with high and low potential for regional metastatic lymph node seeding. Subsequently comparing the expression profiles of these lesions with respectively high and low metastatic potential identified an aggressiveness signature.  Overlaying this signature with single cell data identified proliferative luminal cells, an adipose derived CAF state and a specific subtype of arterial endothelial cells. Assessing the prognostic value of these cell states in an independent dataset (TCGA-PRAD) confirmed their association with regional metastatic lymph node seeding and progression free survival and unveiled a complementary role for the proliferative luminal cells and the adipose derived CAF state in driving regional metastatic lymph node seeding. 
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
