# Load TCGA cohort
DN_counts <- read.table(paste0(PATH_DATA, "de_novo_counts_47.csv"), header=T, row.name=1, sep=",")
DN_counts <- subset(DN_counts, select = -M1RP_ID23_RP4)
colnames(DN_counts) <- ifelse(colnames(DN_counts) != "M1RP_ID23_RP44",
                              colnames(DN_counts), "M1RP_ID23_RP4") 

LA_counts <- read.table(paste0(PATH_DATA, "locally_advanced_counts.csv"), header=T, row.name=1, sep=",")
colnames(LA_counts) <- ifelse(colnames(LA_counts) != "HR_ID5_MLNL", colnames(LA_counts), "HR_ID5_NL") #one sample was wrongly labeled

# Expression 
tcga <- read.csv(paste0(PATH_DATA, "TCGA-PRAD.htseq_counts_nd.tsv.gz"), sep="\t", row.names = 1)
dim(tcga)

# RP vs. NL
sample = read.csv(paste0(PATH_DATA, "clinical/biospecimen.project-TCGA-PRAD.2021-11-17/sample.tsv"), sep="\t")
sample = sample[,c("sample_submitter_id", "sample_type")]
sample = unique(sample$sample_submitter_id[sample$sample_type=="Primary Tumor"])
sample = gsub("-", '.', sample)

tcga = tcga[,colnames(tcga)[colnames(tcga) %in% sample]]

# LN status
tcga_info <- read.csv(paste0(PATH_FIGSHARE, "tcga_RP_info", ".csv"), row.names = 1)

keep <- filterByExpr(tcga, min.prop = 0.1)
tcga <- tcga[keep,]
tcga <- tcga[apply(tcga, 1, var) >= 0.1,]
tcga <- log2(tcga + 1)
tcga <- as.data.frame(scale(tcga, center = colMeans(tcga), scale =  apply(tcga, 2, sd)))
