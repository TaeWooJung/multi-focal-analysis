## PACKAGES --------------------------------------------------------------------
library("Seurat")
library("ggplot2")
library("dplyr")
library("tibble")
library("openxlsx")
library("harmony")

## PART 1: EXPLORATION OF PROVIDED SEURAT OBJECT -------------------------------
seuratObj <- readRDS("chenSeuratObj.rds")
cellannot <- read.csv("Documentation/cellannotation.csv")

seuratObj$cell.type <- seuratObj@active.ident
seuratObj$sample.id <- cellannot$sampleID
seuratObj$cell.type.paper <- cellannot$type
seuratObj$status <- cellannot$status
seuratObj$sample.site <- "Primary"
seuratObj@meta.data[seuratObj@meta.data$sample.id=="SC172", "sample.site"] <- "Metastatic lymph node"


## PART 2: BATCH EFFECT CORRECTION WITH HARMONY --------------------------------
seuratObj <- RunHarmony(seuratObj, group.by.vars = "sample.id", theta = 2)
seuratObj <- RunUMAP(seuratObj, dims = 1:30, reduction = "harmony", verbose = F)
seuratObj <- FindNeighbors(seuratObj, dims = 1:30, reduction = "harmony", verbose = F)
seuratObj <- FindClusters(seuratObj, resolution = 0.5, verbose = F)


## PART 3: ANNOTATION OF NEW CLUSTERING ----------------------------------------
seuratObj <- RenameIdents(
  seuratObj,
  `0` = "E_Luminal",
  `1` = "E_Luminal",
  `2` = "E_Luminal",
  `3` = "L_T/NK cells",
  `4` = "E_Luminal",
  `5` = "En_Vascular (and lymph) 1",
  `6` = "M_Macrophages/Monocytes",
  `7` = "F_Myofibroblasts/Pericytes",
  `8` = "M_Mast cells",
  `9` = "E_Luminal",
  `10` = "E_Luminal",
  `11` = "E_Luminal",
  `12` = "L_B/Plasma cells",
  `13` = "X_Unknown",
  `14` = "E_Luminal",
  `15` = "En_Vascular (and lymph) 2",
  `16` = "En_Vascular (and lymph) 3",
  `17` = "E_Luminal",
  `18` = "F_Aggressive",
  `19` = "E_Luminal",
  `20` = "L_T proliferating"
)

seuratObj@meta.data$cell.type.jef.L1 <- seuratObj@active.ident

saveRDS(seuratObj, "main.rds")

## PART 4: SUBCLUSTERING + ANNOTATION ------------------------------------------
### LUMINAL
#### Subset Epithelial cells
epithelial <- subset(
  seuratObj,
  idents = "E_Luminal"
)

#### Normalization, batch effect correction and clustering 
epithelial <- NormalizeData(epithelial, verbose = F)
epithelial <- FindVariableFeatures(epithelial, nfeatures = 1000, verbose = F)
epithelial <- ScaleData(epithelial, features = VariableFeatures(epithelial), verbose = F)
epithelial <- RunPCA(epithelial, reduction.name = "pca", verbose = F)

epithelial <- RunHarmony(epithelial, group.by.vars = "sample.id", theta = 2)
epithelial <- RunUMAP(epithelial, dims = 1:30, reduction = "harmony", verbose = F)
epithelial <- FindNeighbors(epithelial, dims = 1:30, reduction = "harmony", verbose = F)
epithelial <- FindClusters(epithelial, resolution = 0.5, verbose = F) 

#### Annotation
epithelial <- RenameIdents(
  epithelial,
  `0` = "Luminal",
  `1` = "Luminal",
  `2` = "Luminal",
  `3` = "Luminal",
  `4` = "Luminal",
  `5` = "Luminal",
  `6` = "Luminal",
  `7` = "Basal",
  `8` = "Luminal",
  `9` = "Luminal (Proliferating)",
  `10` = "Luminal",
  `11` = "Luminal (Ribo, translation)",
  `12` = "Hematopoietic"
)

epithelial@meta.data$cell.type.jef.L2 <- epithelial@active.ident

#### Only retain endothelial cells annotated in paper
epithelial <- subset(
  epithelial,
  cell.type.paper %in% c("Basal/Intermediate",
                         "Luminal",
                         "Luminal_agg")
)

saveRDS(epithelial, "epithelial.rds")


### FIBROBLASTS
#### Subset fibroblast cells
fibroblasts <- subset(
  seuratObj,
  idents = c("F_Myofibroblasts/Pericytes", "F_Aggressive")
)

#### Normalization, batch effect correction and clustering 
fibroblasts <- NormalizeData(fibroblasts, verbose = F)
fibroblasts <- FindVariableFeatures(fibroblasts, nfeatures = 1000, verbose = F)
fibroblasts <- ScaleData(fibroblasts, features = VariableFeatures(fibroblasts), verbose = F)
fibroblasts <- RunPCA(fibroblasts, reduction.name = "pca", verbose = F)

fibroblasts <- RunHarmony(fibroblasts, group.by.vars = "sample.id", theta = 2)
fibroblasts <- RunUMAP(fibroblasts, dims = 1:30, reduction = "harmony", verbose = F)
fibroblasts <- FindNeighbors(fibroblasts, dims = 1:30, reduction = "harmony", verbose = F)
fibroblasts <- FindClusters(fibroblasts, resolution = 0.5, verbose = F) 

#### Annotation
fibroblasts <- RenameIdents(
  fibroblasts,
  `0` = "Myofibroblasts",
  `1` = "Pericytes 1",
  `2` = "Myofibroblasts",
  `3` = "Myofibroblasts",
  `4` = "Pericytes 1",
  `5` = "Adipogenic fibroblasts",
  `6` = "Epithelial",
  `7` = "Unknown",
  `8` = "CAFs (ECM remodeling)",
  `9` = "Pericytes 2",
  `10` = "Immune"
)

fibroblasts@meta.data$cell.type.jef.L2 <- fibroblasts@active.ident

#### Only retain fibroblasts annotated in paper
fibroblasts <- subset(
  fibroblasts,
  cell.type.paper %in% c("Fibroblast",
                         "Fibroblast_agg")
)

saveRDS(fibroblasts, "fibroblasts.rds")


### FIBROBLASTS
#### Subset fibroblast cells
endothelial <- subset(
  seuratObj,
  idents = c("En_Vascular (and lymph) 1", "En_Vascular (and lymph) 2",
             "En_Vascular (and lymph) 3")
)

#### Normalization, batch effect correction and clustering
endothelial <- NormalizeData(endothelial, verbose = F)
endothelial <- FindVariableFeatures(endothelial, nfeatures = 1000, verbose = F)
endothelial <- ScaleData(endothelial, features = VariableFeatures(endothelial), verbose = F)
endothelial <- RunPCA(endothelial, reduction.name = "pca", verbose = F)

endothelial <- RunHarmony(endothelial, group.by.vars = "sample.id", theta = 2)
endothelial <- RunUMAP(endothelial, dims = 1:30, reduction = "harmony", verbose = F)
endothelial <- FindNeighbors(endothelial, dims = 1:30, reduction = "harmony", verbose = F)
endothelial <- FindClusters(endothelial, resolution = 0.5, verbose = F) 

#### Annotation
emb <- endothelial@reductions$umap@cell.embeddings
cls <- emb[which(emb[,"UMAP_1"] > 2 & 
                   emb[,"UMAP_2"] < 0.75 & 
                   emb[,"UMAP_1"] > 2 & 
                   emb[,"UMAP_2"] < 0),] %>% rownames()
cls <- intersect(cls, endothelial@meta.data[endothelial@meta.data$seurat_clusters==0,] %>% rownames())
Idents(endothelial, cells = cls) <- "0.1"

endothelial <- RenameIdents(
  endothelial,
  `0` = "HEVs and Venous ECs",
  `0.1` = "Tip cells",
  `1` = "Tip cells",
  `2` = "Arterial ECs",
  `3` = "HEVs and Venous ECs",
  `4` = "Endothelial-Pericyte",
  `5` = "HEVs and Venous ECs",
  `6` = "HEVs and Venous ECs",
  `7` = "HEVs and Venous ECs",
  `8` = "Lymphatic ECs",
  `9` = "ECs (Unfolded prot resp)"
)

endothelial@meta.data$cell.type.jef.L2 <- endothelial@active.ident

#### Only retain endothelial cells annotated in paper
endothelial <- subset(
  endothelial,
  cell.type.paper %in% c("Endothelial",
                         "Endothelial_agg")
)

saveRDS(endothelial, "endothelial.rds")