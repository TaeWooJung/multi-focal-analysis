# Load signatures
# Progression signature (225 genes)
progression_signature_df <- read.csv(paste0(PATH_SIGNATURE, "progression_signature_top225_PC1", ".csv"))
progression_signature <- progression_signature_df$Genes

# Aggressive signature (427 genes)
signature_info <- read.csv(paste0(PATH_SIGNATURE, "aggressive_signature_DU.csv"))
signature <- signature_info$Gene
signature_down <-  signature_info$Gene[signature_info$DU=="Down"]
signature_up <- signature_info$Gene[signature_info$DU=="Up"]

signature_df = data.frame(gene = signature)
signature_df$down = signature %in% signature_down
signature_df$up = signature %in% signature_up
signature_df$total = rowSums(signature_df[c('down', 'up')])
signature_df$color <- ifelse(signature_df$down, "blue", "red")
rownames(signature_df) <- signature

# Load Chen et al. (2021) dataset
chen <- readRDS(paste0(PATH_DATA, "chenSeuratObj.rds"))
cellannot <- read.csv(paste0(PATH_DATA, "cell_annotation.csv"))
cellannot["cellids"] = gsub("-", '.', cellannot$cellID)
cellannot["type"][cellannot["type"]=="Monolytic"] <- "Monocyte"
cellannot["type"][cellannot["type"]=="B cel"] <- "B"
cellannot["celltype"] = cellannot["type"]
cellannot["celltype"][cellannot["celltype"]=="Cell Cycle"] <- "Luminal"
chen$sampleID <- cellannot$sampleID
chen$celltype <- cellannot$celltype
print(unique(chen$celltype))
chen$status <- cellannot$status
Idents(chen) <- chen@meta.data$celltype
chen <- subset(chen, subset = celltype != "unknown")
chen_info <- chen@meta.data
sorted_idents <- c("Luminal", "Basal/Intermediate", 'Fibroblast', 'Endothelial', 'Monocyte', 'Mast', 'T', 'B')
Idents(chen) <- factor(x =Idents(chen), levels = sorted_idents)
chen_norm <- NormalizeData(chen)
chen_primaries <- subset(chen_norm, subset = sampleID %in% c("SC171", "SC172"), invert=T)

chen_primaries_info <- chen_primaries@meta.data
table(chen_primaries$sampleID)

# 380 genes out of 427 genes (aggressive signature) overlap 
signature_chen <- intersect(signature, rownames(chen_primaries))

set.seed(30)
random_signature <- sample(rownames(chen_primaries), length(signature_chen))

# plot <- DimPlot(chen_primaries)
# ggsave(paste0(PATH_FIG_SC, "SupplementaryFigure8.1_chen_UMAP_RP", ".svg"), plot=plot , width = 8, height = 5)
# write.csv(chen_primaries@reductions$umap@cell.embeddings, paste0(PATH_FIGSHARE, "SupplementaryFigure8.1_chen_UMAP_RP", ".csv"))

# Markers of specific cells
# for (cellname in c("Luminal", "Basal/Intermediate", "Fibroblast", "Endothelial", "Mast", "Monocyte", "T", "B")){
#   cells <- rownames(chen_primaries_info[which(chen_primaries_info$celltype == cellname), ])
#   markers  <- FindMarkers(chen_primaries, ident.1=cells,
#                           features=signature_chen,
#                           logfc.threshold=0,
#                           min.pct=0,
#                           min.cells.group=1,
#                           only.pos=TRUE,
#                           random.seed=1)
#   sign_markers <- markers[which(markers$p_val_adj<0.05),]
#   print(cellname)
#   DotPlot(chen_primaries, features=rownames(sign_markers)) + RotatedAxis() + coord_flip()
#   if (cellname == "Basal/Intermediate"){
#     cellname="Basal"
#   }
#   write.csv(sign_markers, paste0(PATH_SCMARKERS, "Markers_", cellname, ".csv"), row.names=TRUE)
# }

# epithelial_cells <- c("Basal/Intermediate", "Luminal")
# stromal_cells <- c("Fibroblast", "Endothelial")
# immune_cells <- c("B", "T", "Mast", "Monocyte")
# cellnames <- c("Epithelium", "Stromal", "Immune")
# 
# celltypes <- list(epithelial_cells, stromal_cells, immune_cells)
# n <- 1
# 
# for (celltype in celltypes){
#   cells <- rownames(chen_primaries_info[which(chen_primaries_info$celltype %in% celltype), ])
#   markers  <- FindMarkers(chen_primaries, ident.1=cells,
#                             features=signature_chen,
#                             only.pos=TRUE)
#   sign_markers <- markers[which(markers$p_val_adj<0.05),]
#   print(celltype)
#   print(sign_markers)
#   write.csv(sign_markers, paste0(PATH_SCMARKERS, "Markers_", cellnames[n], ".csv"), row.names=TRUE)
#   n <- n + 1
# }
