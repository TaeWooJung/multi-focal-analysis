#####################################
# analysis of Chen single cell data #
#####################################

#export LD_LIBRARY_PATH=/home/bioit/mhecke/anaconda3/envs/Seurat_agro/lib:$LD_LIBRARY_PATH
#export PKG_CONFIG_PATH=~/anaconda3/envs/Seurat_agro/lib/pkgconfig
#ulimit -s unlimited
library(Seurat)
library(fossil)
library(ggplot2)
library(pheatmap)
library(infercnv)
library(grid)
library(png)
library(ggplot2)
library(gridExtra)

#reproducability
set.seed(45)

outDir <- "~/Documents/prostateMaarten/single_cell_Chen/clustering"
dir.create(outDir)

#input data
#downloaded from GEO (GSE141445)
counts <- read.table("~/Documents/prostateMaarten/single_cell_Chen/GSM4203181_data.raw.matrix.txt")
chenClusters <- readRDS(file="/Users/taewoojung/Documents/PhD/Papers/GenomeMedicine/Data/epiTumor.basal.clusterID.rds") #only 16 epithelial clusters (see fig2a)

chen <- CreateSeuratObject(counts = counts)
chen

genes <- rownames(chen)

#cells already filtered
VlnPlot(chen, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
FeatureScatter(chen, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

#check clustering
names(chenClusters) <- gsub("-", ".", names(chenClusters))
chen[["chenCluster"]] <- chenClusters
#add additional cluster for NA's (= non epithelial cells)
levels(chen$chenCluster) <- c(levels(chen$chenCluster),"nonEpithelial")
chen$chenCluster[is.na(chen$chenCluster)] <- "nonEpithelial"

#normalize data
chen <- NormalizeData(chen, normalization.method = "LogNormalize", scale.factor = 10000)

#perform PCA
chen <- FindVariableFeatures(chen, selection.method = "vst", nfeatures = 2000)
VariableFeaturePlot(chen)
chen <- ScaleData(chen)
chen <- RunPCA(chen, features = VariableFeatures(chen))
DimPlot(chen, reduction = "pca")

chen <- JackStraw(chen, num.replicate = 100, dims = 40)
chen <- ScoreJackStraw(chen, dims = 1:40)
JackStrawPlot(chen, dims = 1:40)
ElbowPlot(chen, ndims = 30)

#use nr of PCs = 10 based on elbow plot
numPC = 10

#generating our own clusters
chen <- FindNeighbors(chen, dims = 1:numPC)
chen <- FindClusters(chen, resolution = 0.5)
#we found 20 clusters

#visualize using UMAP
chen <- RunUMAP(chen, dims = 1:numPC)
pdf(paste0(outDir, "/ownClusters_umap.pdf"))
DimPlot(chen, reduction = "umap") + ggtitle("own Seurat clusters")
dev.off()

#visualize using tsne
chen <- RunTSNE(chen, dims = 1:numPC)
pdf(paste0(outDir, "/ownClusters_tsne.pdf"))
DimPlot(chen, reduction = "tsne") + ggtitle("own Seurat clusters")
dev.off()

#check similarity with clusters from paper
pdf(paste0(outDir, "/chenClusters_umap.pdf"))
DimPlot(chen, reduction = "umap", group.by='chenCluster')
dev.off()

#compute rand index 
adj.rand.index(as.numeric(chen$seurat_clusters), as.numeric(as.factor(chen$chenCluster)))

#get marker genes for each cluster
chen.markers <- FindAllMarkers(chen, min.pct = 0.25, logfc.threshold = 0.25 )
clusterMarkers <- chen.markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC)

head(clusterMarkers)
write.table(clusterMarkers, file = "~/Documents/prostateMaarten/single_cell_Chen/markerGenes.tsv", row.names = F)

#get marker genes from paper chen
#chenMarkers = read.table("~/Documents/prostateMaarten/single_cell_Chen/markerGenesChen.tsv", sep = '\t', header = T)
orderedMarkers = read.table("~/Documents/prostateMaarten/single_cell_Chen/chenMarkersOrder.tsv")
orderedMarkers[,1] = rev(orderedMarkers[,1])

pdf(paste0(outDir, "/markerExpression.pdf"))
DotPlot(chen, features = orderedMarkers) + RotatedAxis() + coord_flip()
dev.off()

#based on fig 1B, we name our clusters
chen[["seuratNumbers"]] <- Idents(object = chen)

# Rename classes + cellMarker for cluster 16
chen <- RenameIdents(object = chen, 
                              `0` = "Luminal", `1` = "Luminal", `2` = "T",
                              `3` = "Luminal", `4` = "Luminal", `5` = "Luminal",
                              `6` = "Endothelial", `7` = "Luminal",`8` = "Luminal",
                              `9` = "Luminal", `10` = "Monolytic", `11` = "Fibroblast",
                              `12` = "Luminal", `13` = "Mast", `14` = "Endothelial",
                              `15` = "Basal/Intermediate", `16` = "B cel", `17` = "unknown",
                              `18` = "Endothelial", `19` = "Luminal", `20` = "Luminal")
pdf(paste0(outDir, "/celltypes_umap.pdf"))
DimPlot(chen, reduction = "umap")
dev.off()

pdf(paste0(outDir, "/celltypes_tsne.pdf"))
DimPlot(chen, reduction = "tsne")
dev.off()

#check with annotation paper -> very similar
par(mfrow=c(1,2))
DimPlot(chen, reduction = "umap")
DimPlot(chen, reduction = "umap",  group.by='chenCluster')

#convert "." to "-" again
seuratClusters = Idents(chen)
names(seuratClusters) = gsub( "\\.", "-", names(seuratClusters))

#write to file
write.table(seuratClusters, file = "~/Documents/prostateMaarten/single_cell_Chen/cellType_seurat.tsv", sep = "\t", col.names = F)

#check with sampleID's
chen[["sampleID"]] <- gsub( ".*\\.", "", colnames(chen))
DimPlot(chen, reduction = "umap",  group.by='sampleID')

DotPlot(chen, features = c("KLK4", "KLK2", "AR")) + RotatedAxis() + coord_flip()

#save seurat object
save(chen, file = paste0(outDir, "/chenSeuratObj.rds"))














#infer CNV for epithelial cells
chen.epithelial <- subset(x = chen, idents = c("Luminal", "Basal/Intermediate"))
chen.epithelial

#subcluster eptihelial cells
chen.epithelial <- FindVariableFeatures(chen.epithelial, selection.method = "vst", nfeatures = 2000)
VariableFeaturePlot(chen.epithelial)

chen.epithelial <- RunPCA(chen.epithelial, features = VariableFeatures(chen.epithelial))
DimPlot(chen.epithelial, reduction = "pca")

chen.epithelial <- JackStraw(chen.epithelial, num.replicate = 100, dims = 40)
chen.epithelial <- ScoreJackStraw(chen.epithelial, dims = 1:40)
JackStrawPlot(chen.epithelial, dims = 1:40)
ElbowPlot(chen.epithelial, ndims = 30)

numPC = 20

chen.epithelial <- FindNeighbors(chen.epithelial, dims = 1:numPC)
chen.epithelial <- FindClusters(chen.epithelial, resolution = 0.5)

DimPlot(chen.epithelial, reduction = "umap")
DimPlot(chen.epithelial, reduction = "umap",  group.by='sampleID')

#most clusters correspond to samples
table(Idents(chen.epithelial), chen.epithelial$sampleID)
pheatmap(table(Idents(chen.epithelial), chen.epithelial$sampleID), display_numbers = T)

#clusters don't really overlap with clusters found in Chen
table(chen.epithelial$seurat_clusters, chen.epithelial$chenCluster)
pheatmap(table(chen.epithelial$seurat_clusters, chen.epithelial$chenCluster), display_numbers = T)

#convert "." to "-" again
epithelialClusters = chen.epithelial$seurat_clusters
names(epithelialClusters) = gsub( "\\.", "-", names(epithelialClusters))

#write new clustering to file
write.table(chen.epithelial$seurat_clusters, file = "~/Documents/prostateMaarten/single_cell_Chen/epithelialClusters.tsv",
            sep = "\t", col.names = F)

#infer CNV
#reproducability
set.seed(45)

#filter duplicate genenames
gene_info = read.table("~/Documents/prostateMaarten/single_cell_Chen/geneInfo.tsv")
gene_info = gene_info[!duplicated(gene_info[,1]),]
write.table(gene_info, "~/Documents/prostateMaarten/single_cell_Chen/geneInfo_dedup.tsv", sep = "\t", col.names = F, row.names = F)

outdir = "~/Documents/prostateMaarten/single_cell_Chen/CNVAnalysis"
infercnv_obj = CreateInfercnvObject(
  raw_counts_matrix="~/Documents/prostateMaarten/single_cell_Chen/GSM4203181_data.raw.matrix.txt",
  annotations_file="~/Documents/prostateMaarten/single_cell_Chen/cellType_seurat.tsv",
  delim="\t",
  gene_order_file="~/Documents/prostateMaarten/single_cell_Chen/geneInfo_dedup.tsv",
  ref_group_names=c("T","Monolytic", "Fibroblast", "Endothelial", "B cel", "Mast"))

infercnv_obj_default = infercnv::run(
  infercnv_obj,
  cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
  out_dir=outdir,
  cluster_by_groups=TRUE, 
  plot_steps=FALSE,
  denoise=TRUE,
  HMM=FALSE,
  no_prelim_plot=TRUE,
  png_res=60
)

#chen did not use ref_groups, but set it to NULL
outdir_noRef = "~/Documents/prostateMaarten/single_cell_Chen/CNVAnalysis_noRef"
infercnv_obj = CreateInfercnvObject(
  raw_counts_matrix="~/Documents/prostateMaarten/single_cell_Chen/GSM4203181_data.raw.matrix.txt",
  annotations_file="~/Documents/prostateMaarten/single_cell_Chen/cellType_seurat.tsv",
  delim="\t",
  gene_order_file="~/Documents/prostateMaarten/single_cell_Chen/geneInfo_dedup.tsv",
  ref_group_names=c(NULL))

infercnv_obj_default = infercnv::run(
  infercnv_obj,
  cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
  out_dir=outdir_noRef,
  cluster_by_groups=FALSE, 
  plot_steps=FALSE
)
#output is a blank png
#following : https://github.com/broadinstitute/infercnv/issues/362
#something wrong with Cairo backend, disabling rasterisation should solve the problem (but is slower)
plot_cnv(infercnv_obj_default,
         out_dir = outdir_noRef,
         useRaster = F)


# per patient CNV figure
chen[["patientID"]] = chen$sampleID
chen$patientID[chen$patientID == 4] = 6 #sample 4 is the metastatic sample from the same patient as primary sample 6

annot <- read.table("~/Documents/prostateMaarten/single_cell_Chen/cellType_seurat_dot.tsv", sep = "\t")
for (patient in unique(chen$patientID)){
  write.table(annot[annot[,1] %in% Cells(chen[,chen$patientID == patient]),], 
              file = paste0("~/Documents/prostateMaarten/single_cell_Chen/cellType_seurat", patient, ".tsv"),
              sep = "\t", 
              row.names = F,
              col.names = F)
}

normalCells <- c("T","Monolytic", "Fibroblast", "Endothelial", "B cel", "Mast")
for (patient in unique(chen$patientID)){
  print(patient)
  
  #get "normal celltypes"
  celltypes <- read.table(paste0("~/Documents/prostateMaarten/single_cell_Chen/cellType_seurat", patient, ".tsv"))
  normalCellsPatient <- normalCells[normalCells %in% celltypes[,2]]
  
  outdir_patient = paste0("~/Documents/prostateMaarten/single_cell_Chen/CNVAnalysis_", patient)
  
  infercnv_obj = CreateInfercnvObject(
    raw_counts_matrix = counts[,colnames(counts) %in% Cells(chen[,chen$patientID == patient])],
    annotations_file=paste0("~/Documents/prostateMaarten/single_cell_Chen/cellType_seurat", patient, ".tsv"),
    delim="\t",
    gene_order_file="~/Documents/prostateMaarten/single_cell_Chen/geneInfo_dedup.tsv",
    ref_group_names=normalCellsPatient)
  
  infercnv_obj_default = infercnv::run(
    infercnv_obj,
    cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
    out_dir=outdir_patient,
    cluster_by_groups=TRUE, 
    plot_steps=FALSE,
    denoise=TRUE,
    HMM=FALSE,
    no_prelim_plot=TRUE,
    png_res=500
  )
}

patients <- unique(chen$patientID)
patients <- patients[!patients %in% c("2", "8", "11")] #error for patient 2, 8, 10,11??
for (patient in patients){
  print(patient)
  
  outdir_patient = paste0("~/Documents/prostateMaarten/single_cell_Chen/CNVAnalysisNoRef_", patient)
  
  infercnv_obj = CreateInfercnvObject(
    raw_counts_matrix = counts[,colnames(counts) %in% Cells(chen[,chen$patientID == patient])],
    annotations_file=paste0("~/Documents/prostateMaarten/single_cell_Chen/cellType_seurat", patient, ".tsv"),
    delim="\t",
    gene_order_file="~/Documents/prostateMaarten/single_cell_Chen/geneInfo_dedup.tsv",
    ref_group_names=NULL)
  
  infercnv_obj_default = infercnv::run(
    infercnv_obj,
    cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
    out_dir=outdir_patient,
    cluster_by_groups=TRUE, 
    plot_steps=FALSE,
    denoise=TRUE,
    HMM=FALSE,
    no_prelim_plot=TRUE,
    png_res=500
  )
}

patients <- as.character(1:13)

plots <- paste0("~/Documents/prostateMaarten/single_cell_Chen/CNVAnalysis_", patients, "/infercnv.png")
plots <- plots[file.exists(plots)]
plots <- lapply(ll <- plots,function(x){
  img <- as.raster(readPNG(x))
  rasterGrob(img, interpolate = FALSE)
})
#aggregate plots
ggsave("~/Documents/prostateMaarten/single_cell_Chen/CNVplot_ref.pdf", 
       marrangeGrob(grobs = plots, nrow=4, ncol=3,top=NULL))


# for non_ref files
plots <- paste0("~/Documents/prostateMaarten/single_cell_Chen/CNVAnalysisNoRef_", patients, "/infercnv.png")
plots <- plots[file.exists(plots)]
plots <- lapply(ll <- plots,function(x){
  img <- as.raster(readPNG(x))
  rasterGrob(img, interpolate = FALSE)
})
#aggregate plots
ggsave("~/Documents/prostateMaarten/single_cell_Chen/CNVplot_nonRef.pdf", 
       marrangeGrob(grobs = plots, nrow=3, ncol=3,top=NULL))

















genesLouise = c("ABRA", "AC011525.4", "AC053503.12", "ACOT11", "ACRV1", "ACTA2", "ACTC1", "ACTG2", "ACTN2", "ADAM33", "ADAMDEC1", "ADAMTS9-AS1", "ADAMTS9-AS2", "ADCY5", "ADCYAP1R1", "ADRA1A", "ADRA1D", "ADTRP", "AF001548.3", "AF001548.4", "AF131217.1", "AHNAK2", "AKR1B10", "ANO5", "ANO7", "AOC3", "AP000892.6", "APOB", "ARHGDIG", "ARHGEF25", "ARL4D", "ASB2", "ASB5", "ASPM", "ATP1A2", "ATP8A2", "BANK1", "BDNF", "BHMT2", "BNC2", "C10orf82", "C10orf95", "C1QL1", "C20orf166-AS1", "C21orf62", "C2orf40", "C3orf70", "C5orf38", "CA14", "CABP1", "CAMK2A", "CAP2", "CASQ2", "CDH20", "CES1", "CES1P1", "CHRDL1", "CHRDL2", "CHRM2", "CHST6", "CIDEC", "CLVS2", "CMTM5", "CNN1", "CNTFR", "CNTN1", "COL10A1", "COL19A1", "COL4A6", "COL6A1", "COL6A6", "COLGALT2", "CORO6", "COX7A1", "CPEB1", "CPNE7", "CRYAB", "CSPG4", "CSRP1", "CST1", "CST2", "CTB-1144G6.5", "CTF1", "CXCL9", "CYP11A1", "CYP2A6", "CYP4B1", "DES", "DNAJB5", "DNAJC5B", "DPYSL3", "DRD2", "EEF1A2", "ELF5", "EVX2", "F2RL2", "FAM155B", "FAM3D", "FAM46B", "FAM83B", "FAT3", "FBXL21", "FBXL22", "FENDRR", "FGF10", "FGFR2", "FHL1", "FLNA", "FLNC", "FOXF1", "FRMD7", "FXYD1", "GADL1", "GAL", "GATA5", "GDF15", "GLDN", "GLP2R", "GLYATL1", "GNAZ", "GPR1", "GPR133", "GPR156", "GPR158", "GSTM5", "HAPLN2", "HIST1H2AL", "HIST1H2BM", "HLF", "HPD", "HPSE2", "HSPB3", "HSPB6", "HSPB7", "HSPB8", "IGSF1", "ILDR2", "IP6K3", "IQCJ-SCHIP1", "ITGA5", "ITGB3", "ITIH6", "JAZF1-AS1", "JPH2", "JPH4", "KANK2", "KANK4", "KCNAB1", "KCNB1", "KCNH2", "KCNH5", "KCNIP3", "KCNMA1", "KCNMB1", "KCNQ4", "KRT16P6", "KRT2", "KY", "LDB3", "LEF1", "LGALS9C", "LGR6", "LINGO2", "LMOD1", "LRRTM1", "LURAP1", "LY6G6D", "LYNX1", "LYVE1", "MAL", "MAMDC2", "MAOB", "MAP1B", "MARCO", "MBNL1-AS1", "MCAM", "MCOLN2", "MEIS2", "MGARP", "MIR1-1HG", "MIR143HG", "MKRN2OS", "MLC1", "MPP2", "MRGPRF", "MROH7", "MRVI1", "MS4A8", "MSMB", "MSRB3", "MTTP", "MYBL2", "MYBPC1", "MYH11", "MYH2", "MYH6", "MYL9", "MYLK", "MYO3A", "MYOC", "MYOCD", "MYZAP", "NACAD", "NCAM1", "NCS1", "NECAB1", "NEXN", "NKX2-1", "NPAS4", "NPY6R", "NT5C1A", "ODF3L1", "OLR1", "OPTC", "OR2L1P", "OR2L2", "OSR1", "OVOL1", "PALLD", "PATE2", "PCA3", "PCED1B-AS1", "PCP4", "PCSK2", "PDLIM7", "PDZRN4", "PENK", "PGM5", "PGPEP1L", "PHGR1", "PHYHIP", "PI16", "PITX1", "PLA2G2C", "PLD5", "PNCK", "PPARGC1A", "PPP1R14A", "PPP1R1B", "PPP1R3C", "PRDM6", "PRIMA1", "PRKCB", "PRRG3", "PRSS50", "PTCHD1", "PTGIS", "RASL12", "RBFOX3", "RBM24", "RBPMS", "RBPMS2", "RCAN2", "REEP1", "RFPL2", "RIMS4", "RND2", "RNF112", "RP11-1334A24.5", "RP11-150O12.6", "RP11-309L24.6", "RP11-394O4.5", "RP11-395B7.2", "RP11-446H18.1", "RP11-46A10.6", "RP11-569G13.3", "RP11-627G23.1", "RP11-6O2.3", "RPE65", "RSPO2", "SBSPON", "SCUBE3", "SEMA3A", "SFTPA2", "SGCA", "SGCG", "SH3BGR", "SH3GL2", "SIX2", "SLC2A4", "SLC6A11", "SLC6A14", "SLC8A1", "SLC8A2", "SLCO1B3", "SLITRK6", "SMOC1", "SMTN", "SMTNL2", "SORBS1", "SOX14", "SPEG", "SPINK1", "STAC", "STMN4", "SULT4A1", "SYNC", "SYNM", "SYNPO2", "SYT13", "TACR2", "TAGLN", "TARID", "TBX4", "TBX5", "TCEAL2", "TGFB1I1", "TMC5", "TMEM158", "TMEM35", "TNS1", "TOX3", "TPM1", "TPM2", "UBXN10", "USP2", "VCL", "VIT", "WFDC1", "WISP2")

#filter out genes that are not expressed
PrctCellExpringGene <- function(object, genes, group.by = "all"){
  if(group.by == "all"){
    prct = unlist(lapply(genes,calc_helper, object=object))
    result = data.frame(Markers = genes, Cell_proportion = prct)
    return(result)
  }
  
  else{        
    list = SplitObject(object, group.by)
    factors = names(list)
    
    results = lapply(list, PrctCellExpringGene, genes=genes)
    for(i in 1:length(factors)){
      results[[i]]$Feature = factors[i]
    }
    combined = do.call("rbind", results)
    return(combined)
  }
}

calc_helper <- function(object,genes){
  counts = object[['RNA']]@counts
  ncells = ncol(counts)
  if(genes %in% row.names(counts)){
    sum(counts[genes,]>0)/ncells
  }else{return(NA)}
}

expressedGenes <- PrctCellExpringGene(chen, genesLouise, group.by = "all")

#filter non-expressed genes (24 genes)
expressedGenes <- expressedGenes[! is.na(expressedGenes$Cell_proportion),]

#filter lowly expressed genes (< 1% of samples express that gene -> 52)
expressedGenes <- expressedGenes[! expressedGenes$Cell_proportion < 0.01,]
nrow(expressedGenes)


p <- FeaturePlot(chen, features = expressedGenes$Markers, combine = FALSE)
for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + NoLegend() + NoAxes()
}
cowplot::plot_grid(plotlist = p)



#visual inspection of genes
for (i in 0:(length(genesLouise):10)){
  lowLim = i*10+1
  upLim = min((i+1)*10, length(genesLouise))
  sigPlot <- FeaturePlot(chen, features = genesLouise[lowLim:upLim])
  pdf(paste0("~/Documents/prostateMaarten/single_cell_Chen/signatureExpression", i, ".pdf"))
  print(sigPlot)
  dev.off()
}

#approach to test for signatures Maarten
# cell_typeA_marker_gene_list <- list(c("Gene1", "Gene2", "Gene3", "Gene4"))
# object <- AddModuleScore(object = object, features = cell_typeA_marker_gene_list, name = "cell_typeA_score")
# FeaturePlot(object = object, features = "cell_typeA_score1")



