#Install packages
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")
library(DESeq2)

install.packages("gplots")
library(gplots)

install.packages("ggplot2")
library(ggplot2)

BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)

#SAMPLES
#GC number       Stage  
#071732          embryo stage XI - sample 1  
#071733          embryo stage XI - sample 2  
#071734          embryo stage XI - sample 3  
#071735          embryo stage XIV - sample 1  
#071736          embryo stage XIV - sample 2  
#071737          embryo stage XIV - sample 3  
#071738          hatchling - sample 1  
#071739          hatchling - sample 2  
#071740          hatchling - sample 3

FeatureCountsdata = read.table("counts.txt", header=TRUE, row.names=1)
FeatureCountsdata <- FeatureCountsdata[,-1]
FeatureCountsdata <- FeatureCountsdata[,-1]
FeatureCountsdata <- FeatureCountsdata[,-1]
FeatureCountsdata <- FeatureCountsdata[,-1]
FeatureCountsdata <- FeatureCountsdata[,-1]

condition <- data.frame(read.table("condition.txt"))

ddsQuantSeq <- DESeqDataSetFromMatrix(countData = FeatureCountsdata, colData = condition, design=~Stage)

smallestGroupSize <- 6
keep <- rowSums(counts(ddsQuantSeq) >= 10) >= smallestGroupSize # keep genes with at least 6 samples with a count of 10 or higher
ddsQuantSeq <- ddsQuantSeq[keep,]

colData(ddsQuantSeq)$Stage<-factor(colData(ddsQuantSeq)$Stage, levels=c('XI','XIV','XX.2'))
colData(ddsQuantSeq)

dds <- DESeq(ddsQuantSeq)
res <- results(dds, alpha = 0.05)

sizeFactors(dds) #which number is used for normalisation
dds@colData #each sample with its size factor per condition

normalized_counts_dds <- counts (dds, normalized=TRUE)
write.csv(normalized_counts_dds, file = "normalized_counts_dds.csv")

resOrdered <- res[order(res$padj),]
head(resOrdered)
write.table(resOrdered, file="resOrdered.txt")
summary(resOrdered)

##to check the fit of the dispersion estimates to the curve:
plotDispEsts(dds)

#count outliers
par(mar=c(8,5,2,2))
boxplot(log10(assays(dds)[["cooks"]]), range=0, las=2)
colnames()<- c("XI_R1","XI_R2","XI_R3","XIV_R1","XIV_R2","XIV_R3","XX_R1", "XX_R2", "XX_R3")


###check for variation sources
rld <- rlog (dds, blind= FALSE) #data transformation for visualization: log transformation, blind = FALSE, which means that differences between cell lines and treatment (the variables in the design) will not contribute to the expected variance-mean trend of the experiment. 
plotPCA (rld, intgroup="Stage") #results are moved to the vector rld


##making a heatmap to visualise within sample variation: lighter colours: HIGHER CORRELATION. Red : lower correlation 
rld_mat <- assay (rld)
rld_cor <- cor(rld_mat)
col = c("#FF0000", "#FF0000", "#FF0000")
par(mar=c(7,4,4,2)+0.1) 
png(filename='withinsamplevariation.png', width=800, height=750)
heatmap.2(rld_cor, col=redgreen(50), scale="row",
          key=TRUE, symkey=FALSE, density.info="none",cexRow=1,cexCol=1,margins=c(35,35),trace="none")
graphics.off()

#significants - getting rid of the outliers & low counts
sigs <- na.omit(resOrdered)
summary(sigs)
write.csv(sigs, file = "deseq_results.csv")

#Pairwise Comparison 

stageXIV_XI <- lfcShrink(dds , contrast=c("Stage", "XIV","XI"), type = "normal")
summary (stageXIV_XI)
sigs_stageXIV_XI <- na.omit(stageXIV_XI)
summary(sigs_stageXIV_XI)
sigs_stageXIV_XI <- sigs_stageXIV_XI[sigs_stageXIV_XI$padj < 0.05,]
write.csv(sigs_stageXIV_XI, file = "deseq_results_sigs_stageXIV_XI.csv")

stageXX_XI <- lfcShrink(dds , contrast=c("Stage", "XX.2", "XI"), type = "normal")
summary (stageXX_XI)
sigs_stageXX_XI <- na.omit(stageXX_XI)
summary(sigs_stageXX_XI)
sigs_stageXX_XI <- sigs_stageXX_XI[sigs_stageXX_XI$padj < 0.05,]
write.csv(sigs_stageXX_XI, file = "deseq_results_sigs_stageXX_XI.csv")

stageXX_XIV <- lfcShrink(dds , contrast=c("Stage", "XX.2", "XIV"), type = "normal")
summary (stageXX_XIV)
sigs_stageXIV_XX <- na.omit(stageXX_XIV)
summary(sigs_stageXX_XIV)
sigs_stageXX_XIV <- sigs_stageXIV_XX[sigs_stageXX_XIV$padj < 0.05,]
write.csv(sigs_stageXX_XIV, file = "deseq_results_sigs_stageXX_XIV.csv")

#To visualize as heatmap the results (adapted from: https://github.com/mousepixels/sanbomics_scripts/blob/main/tutorial_complex_Heatmap.Rmd)

df_sigs <- as.data.frame(sigs)

OctVul_Map <- read.csv('oct_mus_anno_1B.csv', header = FALSE)

keys <- OctVul_Map$V1
values <- OctVul_Map$V2

l <- list()
for (i in 1:length(keys)){
  l[keys[i]] <- values[i]
}

#for non-mapped labels
no_values <- setdiff(rownames(sigs), keys)
for (i in 1:length(no_values)){
  l[no_values[i]] <- 'NA'
}

df_sigs$symbol <- unlist(l[rownames(df_sigs)], use.names = FALSE)

mat<-assay(rld)[rownames(df.top), rownames(condition)] #sig genes x samples
colnames(mat) <- rownames(condition)
base_mean <- rowMeans(mat)
mat.scaled <- t(apply(mat, 1, scale)) #center and scale each column (Z-score) then transpose
colnames(mat.scaled)<- c("XI_R1","XI_R2","XI_R3","XIV_R1","XIV_R2","XIV_R3","XX_R1", "XX_R2", "XX_R3")

specific_gene_names <- c("OCTVUL_1B029026", "OCTVUL_1B026195", "OCTVUL_1B025607", "OCTVUL_1B010475", "OCTVUL_1B001263", "OCTVUL_1B031766", "OCTVUL_1B021910", "OCTVUL_1B010083", "OCTVUL_1B022853", "OCTVUL_1B008109", "OCTVUL_1B028282", "OCTVUL_1B027308","OCTVUL_1B027203", "OCTVUL_1B023609", "OCTVUL_1B031461", "OCTVUL_1B023804", "OCTVUL_1B015743", "OCTVUL_1B024495")
rows_keep <- which(rownames(mat.scaled) %in% specific_gene_names,nrow(mat.scaled))
rows_keep

h <- Heatmap(mat.scaled[rows_keep,], cluster_rows = F, 
              column_labels = colnames(mat.scaled),row_labels = df.top$symbol[rows_keep], name="Z-score",
              cluster_columns = F)
#row_labels = df.top$symbol[rows_keep] is the part which labels the row names, if you want to keep the OctVul number as row names, remove it from the script above.
h 

png("./heatmap.png", res = 300, width = 4000, height = 3500)
print(h)
dev.off()

