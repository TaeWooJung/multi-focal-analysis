# List of packages you want to install
packages_to_install <- c("DESeq2", "tidyverse", "edgeR", "stringr", "Seurat", "ggplot2", "svglite", "readr", "limma", "gridExtra", "ComplexHeatmap", "survival", "survminer", "fgsea", "org.Hs.eg.db", "DT", "patchwork", "stats")

# Loop through each package
for (package_name in packages_to_install) {
  # Check if the package is already installed
  if (!requireNamespace(package_name, quietly = TRUE)) {
    # If not installed, install the package
    install.packages(package_name)
  }
}

# Load all the required packages
library(DESeq2)
library(tidyverse)
library(edgeR)
library(stringr)
library(ggplot2)
library(svglite)
library(readr)
library(limma)
library(gridExtra)
library(fgsea)
library(org.Hs.eg.db)
library(DT)
library(patchwork)
library(stats)
library(Seurat)
library(survival)
library(survminer)
library(ComplexHeatmap)

# Path information
PATH_SOURCE <- "/Users/taewoojung/Documents/PhD/Papers/GenomeMedicine/" # Change this to your local path
PATH_DATA <- paste0(PATH_SOURCE, "Data/")
PATH_DESEQ2 <- paste0(PATH_SOURCE, "Results/DESEq2/")
PATH_SIGNATURE <- paste0(PATH_SOURCE, "Results/Signatures/")
PATH_GSEA <- paste0(PATH_SOURCE, "Results/GSEA/")
PATH_FIG_SIGNATURE <- paste0(PATH_SOURCE, "Figures/Signatures/")
PATH_FIG_SEEDING_IDENTIFICATION <- paste0(PATH_SOURCE, "Figures/SeedingLesionIdentification/")
PATH_FIG_SCMARKERS <- paste0(PATH_SOURCE, "Figures/SingleCellMarkers/")
PATH_FIG_HEATMAPS <- paste0(PATH_SOURCE, "Figures/Heatmaps/")
PATH_FIG_SURVIVAL <- paste0(PATH_SOURCE, "Figures/SurvivalCurves/")
PATH_FIG_GSEA <- paste0(PATH_SOURCE, "Figures/GSEA/")
PATH_FIG_DOTPLOT <- paste0(PATH_SOURCE, "Figures/DotPlots/")
PATH_FIG_SC <- paste0(PATH_SOURCE, "Figures/SingleCell/")
PATH_FIGSHARE <- paste0(PATH_SOURCE, "FigShare/")

# Commercial signatures
prostadiag = c('REPS2', 'ACADL', 'SLC15A2', 'SLC22A3', 'ANO7', 'FMOD', 'HGD', 'CD38', 'AFF3', 'GREB1', 'AZGP1', 'ANPEP', 'NCAPD3', 'ASPN', 'CHRNA2', 'PAK1IP1', 'KHDRBS3', 'MGP', 'COL1A2', 'MOXD1', 'SFRP2', 'FRZB', 'SPARC', 'CDH11', 'SULF1', 'COL3A1', 'MS4A6A', 'COMP', 'COL10A1', 'NOX4', 'CXCL14',  'THBS2', 'VCAN', 'COL8A1', 'SFRP4', 'ASPN', 'COL1A1')
decipher = c('CAMK2N1', 'PBX1', 'LASP1', 'RABGAP1', 'IQGQP3', 'NFIB', 'THBS2', 'S1PR4', 'ANO7', 'PCDH7', 'MYBPC1', 
             'EPPK1', 'TSBP1', 'NUSAP1', 'ZWILCH', 'UBE2C', 
             'PCAT.32', 'GLYATL1P4', 'PCAT.80', 'TNFRSF19')
oncotypeDX = c('SFRP4', 'BGN', 'COL1A1', 'KLK2', 'SRD5A2', 'FAM13C', 'AZGP1', 'GSN', 'GSTM2', 'TPM2', 'FLNC', 'TPX2')
polaris = c('OXM1', 'CDC20', 'CDKN3', 'CDC2', 'KIF11', 'KIAA0101', 'NUSAP1', 'CENPF', 'ASPM', 'BUB1B', 'RRM2', 'DLGAP5', 'BIRC5', 'KIF20A',
            'PLK1', 'TOP2A', 'TK1', 'PBK', 'ASF1B', 'C18orf24', 'RAD54L', 'PTTG1', 'KIF4A', 'CDCA3', 'MCM10', 'PRC1', 'DTL', 'CEP55',
            'RAD51', 'CENPM', 'CDCA8', 'OIP5', 'SHCBP1', 'ORC6L', 'CCNB1', 'CHEK1', 'TACC3', 'MCM4', 'FANCI', 'KIF15', 'PLK4', 'APOBEC3B', 'NCAPG',
            'TRIP13', 'KIF23', 'NCAPH', 'TYMS', 'GINS1', 'STMN1', 'ZWINT', 'BLM', 'TTK', 'CDC6', 'KIF2C', 'RAD51AP1', 'NCAPG2')
# Luminal markers
proliferative <- c('MKI67', 'TOP2A',  'CDC20', 'CCNB1', 'CENPF', 'PTTG1')
luminal <- c('KRT8', 'KLK4', 'KLK3', 'MSMB', 'MT1E', 'ACPP', 'PA2G2A')
# Fibroblast markers
myofibroblasts <- c('MYH11', 'ACTA2', 'PLN', 'ACTG2', 'HRH2', 'SORBS1')
preadipocyte <- c('CFD', 'APOD', 'MGP', 'CXCL14')
ASCS <- c('PI16', 'COMP', 'ELN', 'EFEMP1', 'POSTN', 'MFAP5', 'IGFBP3', 'INHBA', 'MGP', 'FNDC5', 'MYH10', 'SFRP4', 'F2R', 'IGFBP5', 'VCAN', 'CCDC80', 'TIMP1', 'PTGIS', 'SERPINE2')
general_cafs <- c('FAP', 'S100A4', 'VIM', 'PDGFRB', 'PDPN', 'TNC')
matrix_cafs <- c('PDGFRA', 'DCN', 'LUM', 'VCAN', 'MFAP5', 'POSTN', 'FBLN1', 'FBLN2', 'LOX', 'LOXL1', 'CXCL14')
cafs_col11a1 <- c('FN1', 'CTHRC1', 'COMP', 'SFRP4', 'COL10A1', 'COL11A1', 'INHBA', 'THBS2')
# Endothelial markers
aec <- c('FBLN5', 'GJA5', 'FN1')
peri <- c('NDUFA4L2', 'RGS5', 'HIGD1B') # 'ACTA2'
hev <-  c('ACKR1', 'SELP')
lec <- c('PROX1')
tip <- c('ESM1', 'INSR')
