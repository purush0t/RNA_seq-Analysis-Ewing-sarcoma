#packages required!
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
  BiocManager::install("apeglm")
library(apeglm)
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

#LFC shrinkage MA PLOT
resLFC <- lfcShrink(filter_dseqobj, coef="condition_shEF1_vs_shCTR", type="apeglm")
plotMA(resLFC,ylim=c(-2,2), main ="LFC Shrinkage MA Plot" )



#volcano plot
#Volcano Plot for EWSR1-FLI1 knock-down vs control
EnhancedVolcano(df_dseq,
                lab = df_dseq$gene_id,
                x = 'log2FoldChange',
                y = 'pvalue',title = ' Volcano Plot for EWSR1-FLI1 knock-down vs control',
                labSize = 2, pointSize = 3.0)

#heatmap for upregulated genes
upregular1 <- df_dseq %>%
  mutate(gene_expression = df_dseq$log2FoldChange >0.1 & df_dseq$padj <0.05)%>%
  filter(gene_expression == TRUE) %>%
  arrange(desc(.$padj)) %>%
  distinct(gene_id, .keep_all = TRUE)%>%
  drop_na(gene_id) %>%
  .[1:10,]
rownames(upregular1)

temp1_dseqobj <- filter_dseqobj
rownames(temp1_dseqobj)<- df_dseq$ENSEMBL
temp1_dseqobj
nr_counts <- counts(temp1_dseqobj,normalized = T)[rownames(upregular1),]
nr_counts2 <- t(apply(nr_counts,1,scale))
colnames(nr_counts2) <- rownames(condition_data)

up_genes <- merge(nr_counts2,upregular1$gene_id, by = 0) %>%
  `rownames<-`(.$y)
#silicing the sample columns
up_genes <- up_genes[2:8]
colnames(up_genes)

#plotting upregulated heat map
pheatmap(up_genes,scale = 'row',main = "Upregulated Genes",cluster_rows = F,cluster_cols = F )


#heatmap for downregulated genes

downregular1 <- df_dseq %>%
  mutate(gene_expression = df_dseq$log2FoldChange < 0.1 & df_dseq$padj <0.05)%>%
  filter(gene_expression == TRUE) %>%
  arrange(.$padj) %>%
  distinct(gene_id, .keep_all = TRUE)%>%
  drop_na(gene_id) %>%
  .[1:10,]

downregular1

nr_counts <- counts(temp1_dseqobj,normalized = T)[rownames(downregular1),]
nr_counts2 <- t(apply(nr_counts,1,scale))
colnames(nr_counts2) <- rownames(condition_data)

down_genes <- merge(nr_counts2,downregular1$gene_id, by = 0) %>%
  `rownames<-`(.$y)
down_genes <- down_genes[2:8]
#plotting upregulated heat map
pheatmap(down_genes,scale = 'row',main = "Downregulated Genes",cluster_rows = F,cluster_cols = F )



#gene enrichment analysis using gseGO function from clusterprofiler

EA <- data.frame(log2FoldChange = df_dseq$log2FoldChange,gene_id =  df_dseq$gene_id,ENSEMBL= df_dseq$ENSEMBL) %>%
  `rownames<-`(.$ENSEMBL_ids) %>%
  .[!duplicated(.[c("gene_id")]),] %>%
  na.omit(.) %>%
  `rownames<-`(.$gene_id)

all(rownames(EA)) %in% all(EA$gene_id)

#creating a vector for gse go function
gid_lfc <- c("log2FoldChange" = EA[, 1])
names(gid_lfc) <- EA[,2]
gid_lfc <-sort(gid_lfc,decreasing = TRUE)

gse <- gseGO(geneList=gid_lfc, 
             ont = "ALL", 
             keyType =  "SYMBOL", 
             nPerm = 10000, 
             minGSSize = 3, 
             maxGSSize = 1000, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = org.Hs.eg.db, 
             pAdjustMethod = "none")

gse_sum <- as.data.frame(summary(gse))
head(gse_sum)

#dot plot using the gene enrichment analysis data
dotplot(gse,showCategory=6, title = "",
        font.size = 10,
        label_format = 30,
        split=".sign")+facet_grid(.~.sign)
