
#packages
library(infercnv)
library(Seurat)

#filepaths
rawCounts <- "~/Documents/prostateMaarten/single_cell_Chen/GSM4203181_data.raw.matrix.txt"
geneOrder <- "~/Documents/prostateMaarten/single_cell_Chen/geneInfo_dedup.tsv"
cellTypeAnnot <- "~/Documents/prostateMaarten/single_cell_Chen/cellType_seurat.tsv"
outPath <- "~/Documents/prostateMaarten/single_cell_Chen/inferCNV"

dir.create(outPath)

#read counts and celltypes
rawCounts <- read.table(rawCounts)
cellTypeAnnot <- read.table(cellTypeAnnot, row.names = 1)
rownames(cellTypeAnnot) <- gsub("-", "\\.", rownames(cellTypeAnnot))

#get patients
cells <- colnames(rawCounts)
cells <- as.data.frame(cbind(cells, gsub(".*\\.", "", cells)))
colnames(cells) <- c("cell", "patient")
cells$patient[cells$patient == '6'] <- '4' #sample 4&6 are from the same patient

patients <- unique(cells$patient)
for (patient in patients){
  print(patient)

  #subset data
  patientCells <- cells$cell[cells$patient == patient]
  patientCount <- rawCounts[,colnames(rawCounts) %in% patientCells]
  patientAnnot <- cellTypeAnnot[rownames(cellTypeAnnot) %in% patientCells, ]
  names(patientAnnot) <- rownames(cellTypeAnnot)[rownames(cellTypeAnnot) %in% patientCells]

  #remove celltypes with only 1 cell
  cellTypeCounts <- table(patientAnnot)
  cellTypesToRemove <- names(cellTypeCounts)[table(patientAnnot) < 2]
  cellsToRemove <- names(patientAnnot)[patientAnnot %in% cellTypesToRemove]
  patientCells <- patientCells[! patientCells %in% cellsToRemove]
  patientCount <- patientCount[! rownames(patientCount) %in% cellsToRemove,]
  patientAnnot <- patientAnnot[! names(patientAnnot) %in% cellsToRemove]


  #write patient Annotation to file
  dir.create(paste0(outPath, "/patient", patient))
  write.table(patientAnnot,
              file = paste0(outPath, "/patient", patient, "/celltypes.tsv"),
              sep = "\t", col.names = F)

  #get "normal" cells
  normalCells <- unique(patientAnnot)
  normalCells <- normalCells[! normalCells %in% c("Luminal", "Basal/Intermediate", "unknown")]

  #create the CNV object
  cnv <- CreateInfercnvObject(
    raw_counts_matrix = patientCount,
    gene_order_file = geneOrder,
    annotations_file = paste0(outPath, "/patient", patient, "/celltypes.tsv"),
    ref_group_names = normalCells,
    delim = "\t"
  )

  cnv <- infercnv::run(cnv,
                       cutoff=0.1,                  #0.1 for 10X data (minimum average read counts per gene for ref cells)
                       out_dir=paste0(outPath, "/patient", patient),
                       cluster_by_groups=TRUE,
                       denoise=TRUE,
                       HMM=TRUE,                    #predict CNV level
                       num_threads=10,
                       no_plot = TRUE)
  plot_cnv(cnv,
           out_dir = paste0(outPath, "/patient", patient),
           cluster_by_groups = TRUE,
           )

  #writes metadata to map_metadata_from_infercnv.txt
  seurat_obj = infercnv::add_to_seurat(infercnv_output_path=paste0(outPath, "/patient", patient),
                                       top_n=10)
  #read file in again
  cnvScores <- read.table(file = paste0(outPath, "/patient", patient, "/map_metadata_from_infercnv.txt"))
  cnvScores <- cnvScores[, grepl("proportion_scaled_cnv_.*", colnames(cnvScores))]
  cnvScores <- rowSums(cnvScores)

  #write scores to file
  write.table(cnvScores, file = paste0(outPath, "/patient", patient, "/cnvScores.txt"))

  #get top 10% of cells with highest score
  top10 <- floor(0.1 * length(cnvScores))
  top10 <- names(cnvScores)[order(cnvScores, decreasing = T)[1:top10]]

  resExpr <- cnv@expr.data
  top10Avg <- rowMeans(resExpr[, top10])
  corr <- cor(resExpr, top10Avg)

  CNAtable <- as.data.frame(cbind(cnvScores, corr))
  colnames(CNAtable) <- c("cnvScore", "correlation")
  CNAtable$celltype <- patientAnnot
  CNAtable$color <- "black"
  CNAtable$color[CNAtable$cnvScore < 1 & CNAtable$correlation < 0.4] <- "blue"
  CNAtable$color[CNAtable$cnvScore > 1 & CNAtable$correlation > 0.4] <- "red"

  pdf(paste0(outPath, "/patient", patient, "/malignantClassification.pdf"))
  plot(CNAtable$cnvScore, CNAtable$correlation,
       col = CNAtable$color,
       xlab = "CNA score",
       ylab = "CNA correlation to 10% most malignant",
       main = paste0("Patient ", patient))
  dev.off()

  write.table(CNAtable, file = paste0(outPath, "/patient", patient, "/CNAtable.tsv"),
              sep = "\t", row.names = T, col.names = T)

}

#combine with Seurat information
load("~/Documents/prostateMaarten/single_cell_Chen/clustering/chenSeuratObj.rds")
chen[["patientID"]] = chen$sampleID
chen$patientID[chen$patientID == 6] = 4 #sample 4 is the metastatic sample from the same patient as primary sample 6
allCNV <- NULL
for (patient in patients){
  patientChen <- subset(x = chen, subset = patientID == patient)
  
  patientCNAtable <- read.table(paste0(outPath, "/patient", patient, "/CNAtable.tsv"), sep = "\t", header = T)
  patientCNAtable$malignant <- "unresolved"
  patientCNAtable$malignant[patientCNAtable$color == "red"] <- "malignant"
  patientCNAtable$malignant[patientCNAtable$color == "blue"] <- "non-malignant"
  patientChen[["malignant"]] <- patientCNAtable[, "malignant", drop = FALSE]
  
  pdf(paste0(outPath, "/patient", patient, "/malignant_umap.pdf"))
  print(DimPlot(patientChen, reduction = "umap",  cols = c("red", "blue", "black"), group.by='malignant'))
  dev.off()
  
  pdf(paste0(outPath, "/patient", patient, "/malignant_tsne.pdf"))
  print(DimPlot(patientChen, reduction = "tsne", cols = c("red", "blue", "black"), group.by='malignant'))
  dev.off()
  
  allCNV <- rbind(allCNV, patientCNAtable)
}

#figure with all patients
chen[["malignant"]] <-  allCNV[, "malignant", drop = FALSE]

table(chen$malignant, Idents(chen))
write.table(cbind(chen$malignant, as.character(Idents(chen))), paste0(outPath, "/malignantAnnot.tsv"), sep  ="\t",
            row.names = T, col.names = F, quote = F)

pdf(paste0(outPath, "/malignant_umap.pdf"))
DimPlot(chen, reduction = "umap",  cols = c("red", "blue", "black"), group.by='malignant')
dev.off()

pdf(paste0(outPath, "/malignant_tsne.pdf"))
DimPlot(chen, reduction = "tsne", cols = c("red", "blue", "black"), group.by='malignant')
dev.off()

#additional figures for patient 4 (both metastatic & primary sample)
patient = "4"
patientChen <- subset(x = chen, subset = patientID == patient)

patientCNAtable <- read.table(paste0(outPath, "/patient", patient, "/CNAtable.tsv"), sep = "\t", header = T)
patientCNAtable$malignant <- "unresolved"
patientCNAtable$malignant[patientCNAtable$color == "red"] <- "malignant"
patientCNAtable$malignant[patientCNAtable$color == "blue"] <- "non-malignant"
patientChen[["malignant"]] <- patientCNAtable[, "malignant", drop = FALSE]
patientChen[["metastatic"]] <- "primary"
patientChen$metastatic[patientChen$sampleID == "4"] <- "metastatic"

pdf(paste0(outPath, "/patient4_metastatic_umap.pdf"))
DimPlot(patientChen, reduction = "umap", group.by = 'metastatic')
dev.off()

pdf(paste0(outPath, "/patient4_malignant_umap.pdf"))
DimPlot(patientChen, reduction = "umap", cols = c("red", "blue", "black"), group.by = "malignant")
dev.off()

#create the CNV object
# cnv <- CreateInfercnvObject(
#   raw_counts_matrix = rawCounts,
#   gene_order_file = geneOrder,
#   annotations_file = cellTypeAnnot,
#   ref_group_names = 
#   delim = "\t",
#   
# )
# 
# cnv <- infercnv::run(cnv,
#                      cutoff=0.1,                  #0.1 for 10X data (minimum average read counts per gene for ref cells)
#                      out_dir=paste0(outPath, "/all"),
#                      cluster_by_groups=TRUE,      #TODO where does it get patient annotation??
#                      denoise=TRUE,
#                      HMM=TRUE,                    #run HMM to predict CNV level
#                      HMM_report_by = 'cell',       #report HMM prediction per cell
#                      num_threads=10)



