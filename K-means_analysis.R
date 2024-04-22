#K-means clustering analysis was performed on genes 
library(factoextra)
library(Mfuzz)

#root_fpkm_filter--root_fpkm_filter table with genes as rows and samples as columns;

# 1. K-means clustering analysis ----------------------------------------
#(1) Take the mean and standardize the dataset -----------
root_fpkm_filter_mean <- root_fpkm_filter %>% rowMeans %>% t %>% scale%>% t 

#(2) Finding the best cluster value K -----------------------
fviz_nbclust(root_fpkm_filter_mean,method = 'silhouette')#K=5

#(3) Building an expressionset object ----------------------
root_trans_express<-new ('ExpressionSet',exprs=root_fpkm_filter_mean %>%data.matrix)

#(4) Calculate the value of m----------------------------------
mestimate(root_trans_express)#1.72

#(5) cluster-----------------------------------------------------------
root_trans_cluster <- mfuzz(root_trans_express, c = 5, m = 1.72)
 
KEGG pathway enrichment analysis 
library(clusterProfiler)
library(ggplot2)


#KEGG_database—The first column is “gene_id”, and the second column is “KEGG pathway”; The third column is “Description”;

# KEGG enrichment analysis----------------------------------------------------
cluster_kegg_rich<-enricher(gene = root_trans_cluster$gene_id, 
                              TERM2GENE = KEGG_database [c('KEGG pathway', 'gene_id')] 
                              TERM2NAME = KEGG_database [c('KEGG pathway', 'Description')], 
                              pAdjustMethod = 'BH',
                              pvalueCutoff = 0.05,
                              qvalueCutoff = 0.2)


