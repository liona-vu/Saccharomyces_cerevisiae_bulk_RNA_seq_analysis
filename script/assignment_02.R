# binf110 assignment 02

# by: liona vu

# ====================================================
# Load in libraries
# ====================================================

library(tximport)
library(GenomicFeatures)
library(DESeq2)
library(ggplot2)
library(gridExtra)
library(apeglm)
library(pheatmap)
library(ggrepel)
library(BiocManager)
library(fontBitstreamVera) # Required for clusterProfiler
library(fontLiberation) # Required for clusterProfiler
library(tweenr) # Required for clusterProfiler
library(clusterProfiler)
library(org.Sc.sgd.db) # for S. cerevisiae
library(enrichplot)
library(tidyverse)
library(viridis)

# ====================================================
# Quality check for summary statistics
# ====================================================

#Read in summary statistics data from multiqc 
summary <- read.csv("multiqc_data/multiqc_general_stats.txt", sep = "\t")

#plot percent mapped
p1 <- ggplot(data = summary, aes(x = Sample, y = salmon.percent_mapped, fill = Sample)) +
  geom_bar(stat = "identity", colour = "black") +
  theme_minimal() +
  theme(legend.position="none") +
  labs(x = "Accession Numbers", y = "Percent (%) Mapped ") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  coord_cartesian(ylim = c(0,100))

#plot mapped reads in millions
p2 <- ggplot(data = summary, aes(x = Sample, y = salmon.num_mapped, fill = Sample)) +
  geom_bar(stat = "identity", colour = "black") +
  theme_minimal() +
  theme(legend.position="none") +
  labs(x = "Accession Numbers", y = "Mapped reads (Millions) ") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  coord_cartesian(ylim = c(0,10))

# plot duplicates
p3 <- ggplot(data = summary, aes(x = Sample, y = fastqc.percent_duplicates, fill = Sample)) +
  geom_bar(stat = "identity", colour = "black") +
  theme_minimal() +
  theme(legend.position="none") +
  labs(x = "Accession Numbers", y = "Duplicates reads (%) ") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  coord_cartesian(ylim = c(0,100))

#plot total reads in millions
p4 <- ggplot(data = summary, aes(x = Sample, y = fastqc.total_sequences, fill = Sample)) +
  geom_bar(stat = "identity", colour = "black") +
  theme_minimal() +
  theme(legend.position="none") +
  labs(x = "Accession Numbers", y = "Total sequences (Millions) ") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  coord_cartesian(ylim = c(0,10))


summary_plots <- grid.arrange(p1, p2, p3, p4)

# ====================================================
# Differential gene expression analysis - set up 
# ====================================================

#import data using tximport
txdb <- txdbmaker::makeTxDbFromGFF("ncbi_dataset/ncbi_dataset/data/GCF_000146045.2/genomic.gff")
k <-keys(txdb, keytype = "TXNAME")
tx2gene <- biomaRt::select(txdb, k, "GENEID", "TXNAME")

# create metadata of all SRR accession numbers
srr_samples <- c("SRR10551657", "SRR10551658", "SRR10551659", 
                 "SRR10551660", "SRR10551661", "SRR10551662",
                 "SRR10551663", "SRR10551664", "SRR10551665")

# Makes a metadata with conditions/stage
sample_table <- data.frame(run = srr_samples,
                           stage = factor(rep(c("Mature", "Thin", "Early"), each = 3))) 

#make rownames the same as the accession number
rownames(sample_table) <- sample_table$run

#Finds files paths that have the SRR accession numbers and the quant.sf file produces from Salmon
files <- file.path("quants", srr_samples,"quant.sf")

# import transcript data from Salmon
txi <- tximport(files, type = "salmon", tx2gene = tx2gene)
dds <- DESeqDataSetFromTximport(txi, colData = sample_table, design = ~stage)
dds_2 <- DESeq(dds, test = "Wald")

#Check comparison names which will be used in results
resultsNames(dds_2)
#Generate a results table for stages of yeast and comparing between thin and early biofilms
#res_thin_vs_early <- results(dds_2, name="stage_Thin_vs_Early")
res_thin_vs_early <- results(dds_2, contrast=c("stage", "Thin", "Early"))

#Generate a results table comparing between mature and early biofilms
#res_mature_vs_early <- results(dds_2, name="stage_Mature_vs_Early")
res_mature_vs_early <- results(dds_2, contrast = c("stage", "Mature", "Early"))

#Generate a results table comparing between mature and thin biofilms
res_mature_vs_thin <- results(dds_2, contrast= c("stage", "Mature", "Thin"))                   

# Look at tables for all comparisons
res_mature_vs_early
res_thin_vs_early
res_mature_vs_thin

#Shrink LFC for better and accurate gene counts
resLFC_thin_vs_early <- lfcShrink(dds_2, coef="stage_Thin_vs_Early", type="apeglm")
resLFC_mature_vs_early <- lfcShrink(dds_2, coef ="stage_Mature_vs_Early", type="apeglm")

# Relevel factor to get the other comparison, thin vs mature oh my Goooddddddd
dds$stage <- relevel(dds$stage, ref = "Thin") #Use Thin stage as the reference

#Run DESeq again
dds_3 <- DESeq(dds)

resultsNames(dds_3) #Check if Thin vs Mature is here

# Now can run LFC Shrinkage
resLFC_mature_vs_thin <- lfcShrink(dds_3, coef = "stage_Mature_vs_Thin", type="apeglm")

# ====================================================
# Plotting MA plots and checking structure
# ====================================================

# MA plot without LFC shrinkage for quality control to see if normalization occured
par(mfrow = c(2, 3))
plotMA(res_mature_vs_early, ylim=c(-10,10))
plotMA(res_thin_vs_early, ylim=c(-10, 10))
plotMA(res_mature_vs_thin, ylim =c(-10,10))
#THIS FUCKS UP RSTUDIO
#idx <- identify(res_mature_vs_early$baseMean, res_mature_vs_early$log2FoldChange)
#rownames(res_mature_vs_early)[idx]

# MA with apeglm applied
#par(mfrow = c(1, 3))
plotMA(resLFC_mature_vs_early, ylim=c(-10,10))
plotMA(resLFC_thin_vs_early, ylim=c(-10, 10))
plotMA(resLFC_mature_vs_thin, ylim=c(-10,10))

#remove na values
resLFC_mature_vs_early <- na.omit(resLFC_mature_vs_early)
resLFC_thin_vs_early <- na.omit(resLFC_thin_vs_early)
resLFC_mature_vs_thin <- na.omit(resLFC_mature_vs_thin)

# Volcano plot

# need to mark genes as upregulated, downregulated, or not significant to colour them in ggplot
# Exclude genes with under 2-fold change (log2FoldChange < 1)
# Create function to plot volcano plot and avoid repeating code
volcano_plot_fn <- function(lfc_deseq2, title_input) {
  lfc_df <- as.data.frame(lfc_deseq2)
  lfc_df$gene <- rownames(lfc_df)
  lfc_df$significant <- ifelse(lfc_df$padj < 0.05 & abs(lfc_df$log2FoldChange) > 1, 
                                                  ifelse(lfc_df$log2FoldChange > 0, "Up", "Down"), "Not Sig")
  lfc_df <- na.omit(lfc_df)
  
  # Here's some ggplot code for a volcano plot
  # We plot log2foldchange against -log10(adjusted p value)
  
  plot <- ggplot(lfc_df, aes(x = log2FoldChange, y = -log10(pvalue), color = significant)) +
    geom_point() +
    scale_color_manual(values = c("Down" = "blue", "Not Sig" = "gray", "Up" = "red")) +
    labs(x = "Log2 Fold Change", y = "-Log10 p-value", 
         title = paste0("Volcano Plot ", title_input)) +
    theme(legend.position = "right") #+ 
#    theme_minimal()
  
  return(plot)
}

# Generate volcano plot early to mature
vp_1 <- volcano_plot_fn(resLFC_mature_vs_early, "Early to Mature")
vp_2 <- volcano_plot_fn(resLFC_thin_vs_early, "Early to Thin")
vp_3 <- volcano_plot_fn(resLFC_mature_vs_thin, "Thin to Mature")

grid.arrange(vp_2, vp_3, vp_1, ncol = 3)
#top_genes_mature_vs_early <- head(order(abs(resLFC_mature_vs_early$padj), decreasing = FALSE), 20)
#gene_names_mature_vs_early <- rownames(resLFC_mature_vs_early)[top_genes_mature_vs_early]

# Make function to take LFC shrink results, filter by padj< 0.05, and select the top 20 genes
top_genes <- function(dataframe) {
  filtered_df <- as.data.frame(dataframe) %>%
    filter(padj < 0.05) %>%
    arrange(desc(abs(log2FoldChange)))
  
  filtered_LFC <- DESeqResults(filtered_df)
  genes_names <- head(rownames(filtered_LFC), n = 20)
  return(genes_names)
}

#filtered_df <- as.data.frame(resLFC_mature_vs_early) %>%
 # filter(padj < 0.5) %>%
  #arrange(desc(abs(log2FoldChange)))

#length(row.names(x))
#head(x, n = 20)

#resLFC_mature_vs_early_filtered <- DESeqResults(filtered_df)
#class(resLFC_mature_vs_early_filtered)
#gene_names_mature_vs_early <- head(rownames(resLFC_mature_vs_early_filtered), n =20)

#top_genes_thin_vs_early <- head(order(abs(resLFC_thin_vs_early$padj), decreasing = FALSE), n = 20)
#gene_names_thin_vs_early <- rownames(resLFC_thin_vs_early)[top_genes_thin_vs_early]

#top_genes_mature_vs_thin <- head(order(abs(resLFC_mature_vs_thin$padj), decreasing = FALSE), n = 20)
#gene_names_mature_vs_thin <- rownames(resLFC_mature_vs_thin)[top_genes_mature_vs_thin]

# Extract transformed & normalized counts with a variance stabilizing transformation
vsd <- vst(dds_2)

# Store counts in a matrix for the heatmap
gene_names_mature_vs_early <- top_genes(resLFC_mature_vs_early)
mat_mature_vs_early <- assay(vsd)[gene_names_mature_vs_early, ]

# Create heatmap visualizing the top 20 differential expressed genes for mature vs early
p5 <- pheatmap(mat_mature_vs_early, 
               color = inferno(100),
         #color = colorRampPalette(c("green", "white", "red"))(100),
         scale = "row",
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         sample_table = sample_table,
         show_rownames = TRUE,
         show_colnames = FALSE,
         main = "Early vs Mature")

gene_names_thin_vs_early <- top_genes(resLFC_thin_vs_early)
mat_thin_vs_early <- assay(vsd)[gene_names_thin_vs_early, ]

p6 <- pheatmap(mat_thin_vs_early, 
               color = inferno(100),
        # color = colorRampPalette(c("steelblue", "white", "salmon"))(100),
         scale = "row",
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         sample_table = sample_table,
         show_rownames = TRUE,
         show_colnames = TRUE,
         main = "Early vs Thin")

gene_names_mature_vs_thin <- top_genes(resLFC_mature_vs_thin)
mat_mature_vs_thin <- assay(vsd)[gene_names_mature_vs_thin,]

p7 <- pheatmap(mat_mature_vs_thin,
               color = inferno(100),
              # color = colorRampPalette(c("steelblue", "white", "salmon"))(100),
               scale = "row",
               cluster_rows = TRUE,
               cluster_cols = TRUE,
               sample_table = sample_table,
               show_rownames = TRUE,
               show_colnames = TRUE,
               main = "Thin vs Mature")

grid.arrange(p6[[4]], p7[[4]], p5[[4]], ncol = 3)
#PCA Plot

# Get the coordinates using plotPCA from DESeq2
pca_data <- plotPCA(vsd, intgroup = "stage", returnData = TRUE)

# Get percent variance explained by the top two principal components
percentVar <- round(100 * attr(pca_data, "percentVar"))

# GGplot code to display cell lines by colour, and treatment by shape
ggplot(pca_data, aes(x = PC1, y = PC2, color = stage, shape = stage, label = name)) +
  geom_point(size = 4) +
  theme_minimal() +
 # geom_text(hjust = 0.1, vjust = 0.1) +
  geom_label_repel(box.padding = 0.75) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("PCA Plot of Flor Yeast Samples") +
  theme(plot.title = element_text(size=15, hjust = 0.5, face = "bold"))

# ====================================================
# ORA - Gene analysis for early to/vs thin stage
# ====================================================

#keytypes(org.Sc.sgd.db)

#Create function for ORA 
gene_expression_compare <- function(lfc_results, ont) {
  #Create data frame and take names
  as_dataframe <- as.data.frame(lfc_results)
  ensembl_ids <- rownames(as_dataframe)
  
  gene_map <- bitr(ensembl_ids, 
                   fromType = "GENENAME", 
                   toType = c("ENTREZID", "ENSEMBL"),
                   OrgDb = org.Sc.sgd.db)
  
  as_dataframe$GENENAME <- rownames(as_dataframe)
  dataframe_2 <- merge(as_dataframe, gene_map, by = "GENENAME", all.x = TRUE)
  
  #Grabs upregulated genes
  sig_genes_up <- dataframe_2 %>%
    filter(pvalue<0.05) %>%
    filter(log2FoldChange > 1) %>%
    na.omit() %>%
    pull(ENTREZID) %>%
    unique()
  
  #Grabs downregulated genes
  sig_genes_down <- dataframe_2 %>%
    filter(pvalue<0.05) %>%
    filter(log2FoldChange < -1) %>%
    na.omit() %>%
    pull(ENTREZID) %>%
    unique()
  
  #Creates all genes for universe for compareCluster function
  all_genes <- dataframe_2 %>%
    pull(ENTREZID) %>%
    na.omit() %>%
    unique()
  
  #Compare upregulated and downregulated genes
  compare <- compareCluster(geneCluster = list(Upregulated = sig_genes_up,
                                               Downregulated = sig_genes_down),
                                          fun = "enrichGO",
                                          OrgDb = org.Sc.sgd.db,
                                          ont = ont,
                                          pvalueCutoff = 0.05,
                                          pAdjustMethod = "BH",
                                          universe = all_genes, # Benjamini Hochberg
                                          qvalueCutoff = 0.02,
                                          readable = FALSE)
  return(compare)
}

# apply function for ORA for thin vs early stages
compare_thin_vs_early <- gene_expression_compare(resLFC_thin_vs_early, "BP")
test <- clusterProfiler::simplify(x = compare_thin_vs_early, 
                                  cutoff = 0.7,
                                  by = "p.adjust",
                                  select_fun = min,
                                  measure = "Wang")
dp <- dotplot(test, showC)

dot_plot_1 <- dotplot(compare_thin_vs_early, showCategory = 10, title = "GO Biological Process \n Early vs Thin Stage")

# ORA - Gene analysis for thin to/vs mature stage

compare_mature_vs_thin <- gene_expression_compare(resLFC_mature_vs_thin, "BP")
dot_plot_2 <- dotplot(compare_mature_vs_thin, showCategory = 10, 
                      title = "GO Biological Processes \n Thin vs Mature Stage")

# ORA - Gene analysis for early to/vs mature stage
compare_early_vs_mature <- gene_expression_compare(resLFC_mature_vs_early, "BP")
dot_plot_3 <- dotplot(compare_early_vs_mature, showCategory = 10, 
                      title = "GO Biological Processes \n Thin vs Mature Stage")

grid.arrange(dot_plot_1, dot_plot_2, dotplot_3, ncol= 3, nrow = 1)

# GSEA analysis 
# Make a function
gsea_analysis <- function (LFC_results, ont, keytype) {
  gene_list_df <- as.data.frame(LFC_results)
  gene_list_df <- na.omit(gene_list_df)
  
  gene_list <- gene_list_df$log2FoldChange
  
  gene_list_df$NAMES <- rownames(gene_list_df)
  
  names(gene_list) <- gene_list_df$NAMES
  
  gene_list <- sort(gene_list, decreasing = TRUE)
  
  gsea_results <- gseGO(geneList = gene_list,
                        ont = "BP",
                        OrgDb = org.Sc.sgd.db, 
                        keyType = "COMMON",
                        pvalueCutoff = 0.05,
                        verbose =  TRUE, 
                        eps = 0 ) #P values less than 1.0 e-10
return(gsea_results)
                        
}

gsea_mature_vs_early <- gsea_analysis(resLFC_mature_vs_early, "BP", "COMMON")

gsea_plot_mature_vs_early<- gseaplot2(gsea_mature_vs_early, geneSetID = c(1,2,3), 
                                      color = c("#E495A5", "#86B875", "#7DB0DD"),
                                      pvalue_table = TRUE, 
                                      title = "GSEA of Early to Mature biofilm formation")


gsea_early_vs_thin <- gsea_analysis(resLFC_thin_vs_early, "BP", "COMMON")

gsea_plot_early_vs_thin <- gseaplot2(gsea_early_vs_thin, geneSetID = c(1:3),
                                         color = c("#E495A5", "#86B875", "#7DB0DD"),
                                         pvalue_table = TRUE, 
                                         title = "GSEA of Early to Thin biofilm formation")



gsea_thin_vs_mature <- gsea_analysis(resLFC_mature_vs_thin, "BP", "COMMON")

gsea_plot_mature_vs_thin <- gseaplot2(gsea_thin_vs_mature, geneSetID = c(1:3), 
                                         color = c("#E495A5", "#86B875", "#7DB0DD"),
                                         pvalue_table = TRUE, 
                                         title = "GSEA of Early to Thin biofilm formation")

plot_list(gsea_plot_mature_vs_early, gsea_plot_early_vs_thin, gsea_plot_mature_vs_thin)



gseaplot(gsea_results, by = "all", geneSetID = 1)
head(gsea_results)
