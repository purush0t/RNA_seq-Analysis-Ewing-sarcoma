if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install('org.Hs.eg.db')
BiocManager::install('AnnotationDbi')
BiocManager::install("dplyr")
BiocManager::install("SummarizedEzperiment")
BiocManager::install("DESeq2")
BiocManager::install("pheatmap")
BiocManager::install("EnhancedVolcano")
BiocManager::install("clusterProfiler")
BiocManager::install("enrichplot")

library(dplyr)
library(SummarizedExperiment)
library(DESeq2)
library(pheatmap)
library(EnhancedVolcano)
library(org.Hs.eg.db)
library(kableExtra)
library(ggplot2)
library(clusterProfiler)
library(enrichplot)
library(DOSE)
library(ggridges)
library(tidyverse)
library(pathview)
library(DESeq2)

# data is of object type:rangedsummarizedexperiment data 
rse <- readRDS("EwS.rds")
rse

#assay function extracts the primary assay data from rds
primary_countsdata <- assay(rse)
primary_countsdata
dim_primary_countsdata <- dim(primary_countsdata)
#print(paste("dimension of p_countsdata:",dim_primary_countsdata))
dim_primary_countsdata

#colData(rse)
# colData accessing all the meta data
meta_data <- data.frame(colData(rse))

#getiin all the column names from meta data

SRA_accesion <- rownames(meta_data)
class(SRA_accesion)


#condition_data <- data.frame("SRA"=SRA_accesion,"condition"=meta_data$condition)
condition_data <- data.frame(condition = meta_data$condition,
                             row.names = meta_data$run)
condition_data
##type(condition_data)

#removing entries where the read counts are 0
filtered_counts_data <- primary_countsdata[which(rowSums(primary_countsdata)>0),]
head(filtered_counts_data )

dim(filtered_counts_data)

#check if all the sam
all(rownames(condition_data) %in% colnames(filtered_counts_data))


dseq_object <- DESeqDataSetFromMatrix(countData = filtered_counts_data,colData = condition_data, design = ~condition)

dseq_object
#removing entries with low read counts
dseq_temp <- rowSums(counts(dseq_object)) >=100

filter_dseqobj <- dseq_object[dseq_temp]

# the control sample "shCTR is used as reference for Dseq"
filter_dseqobj$Treatment <- relevel(filter_dseqobj$condition,ref = "shCTR")
filter_dseqobj$Treatment

#running the Dseq2

filter_dseqobj <- DESeq(filter_dseqobj)

dseq_result <- results(filter_dseqobj)
dseq_result
summary(dseq_result)

# using pipes to , add new cloumns, remove version no.from ensembl id's , mapping to respective gene name's
df_dseq <- data.frame(dseq_result) %>%
  add_column(gene_id = 0, version = rownames(.))%>%
  mutate(ENSEMBL=gsub("\\..*","",rownames(.)))%>%
  mutate(gene_id = mapIds(org.Hs.eg.db,keys = .$ENSEMBL,keytype = "ENSEMBL",column = "SYMBOL"))%>%
  `rownames<-`(.$ENSEMBL)
#colnames(df_dseq)
#rownames(df_dseq)
head(df_dseq)

#creatin a copy of deseq , with new rownames (so the ensembl id dosent have version no.)
ens_dseq <- filter_dseqobj
rownames(ens_dseq) <- df_dseq$ENSEMBL


#variance stabilization(vst)
#For RNA-seq counts,the expected variance grows with the mean.so we do a variance stabilization transformation ,inside Dseq2 package
normal_dseq <- vst(filter_dseqobj,blind = FALSE)
#normal_dseq


#PCA(principle component analysis) PLOT!

pca_results <-
  plotPCA(
    normal_dseq,
    intgroup = c("condition"),
    returnData = TRUE # This argument tells R to return the PCA values
  )

#plot using ggplot() function
annotated_pca_plot<-ggplot(
  pca_results,
  aes(
    x = PC1,
    y = PC2,
    # plot points with different colors for each `condition` group
    shape=condition,
    color= condition,
    size = 1.5
    # plot points with different shapes for each `condition` group
  )
) +
  geom_point()
annotated_pca_plot


#MA plot
MAplot_results <- plotMA(dseq_result,ylim=c(-2,2),main ="MA Plot",returnData = TRUE)

#The plotMA() function returns a list containing two elements:
# M: A vector of log2 fold change values.
#A: A vector of average log2 expression levels.

MAplot_results$lfc_category <- ifelse(MAplot_results$lfc < 0, "red",
                              ifelse(MAplot_results$lfc == 0,"black", "purple"))

colnames(MAplot_results)

ggplot(MAplot_results,aes(x=mean,y=lfc))+
  geom_point(color = MAplot_results$lfc_category) +
  xlab("Mean Expression") +
  ylab("Log2 Fold Change") +
  ggtitle("MA Plot with Colored Points based on lfc")+
  scale_x_continuous(limits = c(0, 1000), breaks = seq(0, 40000, by = 40000)) +
  scale_y_continuous(limits = c(-3, 3), breaks = seq(-3, 3, by = 1))




