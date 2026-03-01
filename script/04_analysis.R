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
  scale_y_continuous(breaks = seq(0,100, 10)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  coord_cartesian(ylim = c(0,100))

#plot mapped reads in millions
p2 <- ggplot(data = summary, aes(x = Sample, y = salmon.num_mapped, fill = Sample)) +
  geom_bar(stat = "identity", colour = "black") +
  theme_minimal() +
  theme(legend.position="none") +
  labs(x = "Accession Numbers", y = "Mapped reads (Millions) ") +
  scale_y_continuous(breaks = seq(0,10, 1)) +
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


summary_plots <- grid.arrange(p1, p2, ncol = 2)

# ======================================================
# DESeq2 Differential gene expression analysis
# ======================================================

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
res_thin_vs_early <- results(dds_2, contrast=c("stage", "Thin", "Early"))

#Generate a results table comparing between mature and early biofilms
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

# Relevel factor to get the other comparison, thin vs mature
dds$stage <- relevel(dds$stage, ref = "Thin") #Use Thin stage as the reference

#Run DESeq again
dds_3 <- DESeq(dds)

resultsNames(dds_3) #Check if Thin vs Mature is here

# Now can run LFC Shrinkage
resLFC_mature_vs_thin <- lfcShrink(dds_3, coef = "stage_Mature_vs_Thin", type="apeglm")

# =======================================================
# Plotting MA plots, volcano plots, and checking structure
# =======================================================

# MA plot without LFC shrinkage for quality control to see if normalization occured
par(mfrow = c(2, 3))
plotMA(res_mature_vs_early, ylim=c(-10,10))
plotMA(res_thin_vs_early, ylim=c(-10, 10))
plotMA(res_mature_vs_thin, ylim =c(-10,10))

# MA with apeglm applied
plotMA(resLFC_mature_vs_early, ylim=c(-10,10))
plotMA(resLFC_thin_vs_early, ylim=c(-10, 10))
plotMA(resLFC_mature_vs_thin, ylim=c(-10,10))

#remove na values
resLFC_mature_vs_early <- na.omit(resLFC_mature_vs_early)
resLFC_thin_vs_early <- na.omit(resLFC_thin_vs_early)
resLFC_mature_vs_thin <- na.omit(resLFC_mature_vs_thin)

# Volcano plot

# Create function to plot volcano plot and avoid repeating code for categorizing genes as upregulated, downregulated, or not significant to colour them in ggplot
# Exclude genes under 2-fold change (log2FoldChange < 1)
volcano_plot_fn <- function(lfc_deseq2, title_input) {
  lfc_df <- as.data.frame(lfc_deseq2)
  lfc_df$gene <- rownames(lfc_df)
  lfc_df$significant <- ifelse(lfc_df$padj < 0.05 & abs(lfc_df$log2FoldChange) > 1, 
                                                  ifelse(lfc_df$log2FoldChange > 0, "Up", "Down"), "Not Sig")
  lfc_df <- na.omit(lfc_df)
  
  # plot log2foldchange against -log10(adjusted p value)
  plot <- ggplot(lfc_df, aes(x = log2FoldChange, y = -log10(pvalue), color = significant)) +
    geom_point() +
    scale_color_manual(values = c("Down" = "blue", "Not Sig" = "gray", "Up" = "red")) +
    labs(x = "Log2 Fold Change", y = "-Log10 p-value", 
         title = paste0("Volcano Plot ", title_input)) +
    theme(legend.position = "right")
  
  return(plot)
}

# Generate volcano plot for all stages
vp_1 <- volcano_plot_fn(resLFC_mature_vs_early, "Early to Mature")
vp_2 <- volcano_plot_fn(resLFC_thin_vs_early, "Early to Thin")
vp_3 <- volcano_plot_fn(resLFC_mature_vs_thin, "Thin to Mature")

#put all volcano plots next to each other
grid.arrange(vp_2, vp_3, vp_1, ncol = 2)

# Make function to take LFC shrink results, filter by padj< 0.05, and select the top 20 genes
top_genes <- function(dataframe) {
  filtered_df <- as.data.frame(dataframe) %>%
    filter(padj < 0.05) %>%
    arrange(desc(abs(log2FoldChange)))
  
  filtered_LFC <- DESeqResults(filtered_df)
  genes_names <- head(rownames(filtered_LFC), n = 20)
  return(genes_names)
}

# Extract transformed & normalized counts with a variance stabilizing transformation
vsd <- vst(dds_2)

# Store counts in a matrix for the heatmap
gene_names_mature_vs_early <- top_genes(resLFC_mature_vs_early)
mat_mature_vs_early <- assay(vsd)[gene_names_mature_vs_early, ]

# Create heatmap visualizing the top 20 differential expressed genes for mature vs early
annotate_colours <- data.frame(Stage = factor(rep(c("Mature", "Thin", "Early"), 
                                                  each = 3)))

#Creating annotation data to demarcate which samples belong to which biofilm stage,instead of constantly checking metadata
rownames(annotate_colours) <- colnames(mat_mature_vs_early)
ann_colours <- list(Stage = c("Mature" = "#E495A5",
                              "Thin" = "#86B875",
                              "Early" = "#7DB0DD"))

#plot heatmap
p5 <- pheatmap(mat_mature_vs_early, 
               color = inferno(100),
         scale = "row",
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         sample_table = sample_table,
         show_rownames = TRUE,
         show_colnames = TRUE,
         main = "Early to Mature",
         annotation_col = annotate_colours,
         annotation_colors = ann_colours)

#repeat with other stages
gene_names_thin_vs_early <- top_genes(resLFC_thin_vs_early)
mat_thin_vs_early <- assay(vsd)[gene_names_thin_vs_early, ]

p6 <- pheatmap(mat_thin_vs_early, 
               color = inferno(100),
         scale = "row",
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         sample_table = sample_table,
         show_rownames = TRUE,
         show_colnames = TRUE,
         main = "Early to Thin",
         annotation_col = annotate_colours,
         annotation_colors = ann_colours)

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
               main = "Thin to Mature",
              annotation_col = annotate_colours,
              annotation_colors = ann_colours)

#Plot all heatmap together
grid.arrange(p6[[4]], p7[[4]], p5[[4]], ncol = 3)

# Generate PCA Plot
# Get the coordinates using plotPCA from DESeq2
pca_data <- plotPCA(vsd, intgroup = "stage", returnData = TRUE)

# Get percent variance explained by the top two principal components
percentVar <- round(100 * attr(pca_data, "percentVar"))

# GGplot code to display cell lines by colour, and treatment by shape
ggplot(pca_data, aes(x = PC1, y = PC2, color = stage, shape = stage, label = name)) +
  geom_point(size = 4) +
  theme_minimal() +
  geom_label_repel(box.padding = 0.75) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("PCA Plot of Flor Yeast Samples") +
  theme(plot.title = element_text(size=15, hjust = 0.5, face = "bold"))

# ====================================================
# ORA - Gene analysis for all three stages comparison
# ====================================================

#check keytypes for database
keytypes(org.Sc.sgd.db)

#Create function for ORA analysis that takes in lcf shrinkage results and ontology type (i.e. BP, CC)
gene_expression_compare <- function(lfc_results, ont) {
  #Create data frame and take names
  as_dataframe <- as.data.frame(lfc_results)
  ensembl_ids <- rownames(as_dataframe)
  
  #Convert biological IDs from genename to entreid and ensembl
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
                                          pAdjustMethod = "BH",# Benjamini Hochberg
                                          universe = all_genes, 
                                          qvalueCutoff = 0.02,
                                          readable = FALSE)
  #Remove redundant GO terms
  compare_2 <- clusterProfiler::simplify(x = compare, 
                            cutoff = 0.7,
                            by = "p.adjust",
                            select_fun = min)
  return(compare_2)
}

# apply function for ORA for thin vs early stages
compare_thin_vs_early <- gene_expression_compare(resLFC_thin_vs_early, "BP")

#plot ORA results
dot_plot_1 <- dotplot(compare_thin_vs_early, showCategory = 10, title = "GO Biological Process \n Early vs Thin Stage")

# ORA - Gene analysis for thin to/vs mature stage
compare_mature_vs_thin <- gene_expression_compare(resLFC_mature_vs_thin, "BP")
dot_plot_2 <- dotplot(compare_mature_vs_thin, showCategory = 10, 
                      title = "GO Biological Processes \n Thin vs Mature Stage")

# ORA - Gene analysis for early to/vs mature stage
compare_early_vs_mature <- gene_expression_compare(resLFC_mature_vs_early, "BP")
dot_plot_3 <- dotplot(compare_early_vs_mature, showCategory = 10, 
                      title = "GO Biological Processes \n Thin vs Mature Stage")

#plot all dotplots together
grid.arrange(dot_plot_1, dot_plot_2, dotplot_3, ncol= 3, nrow = 1)

# ====================================================
# kegg analysis
# ====================================================

# make function for kegg analysis
kegg_analysis <- function(LFC_results) {
  
  as_dataframe <- as.data.frame(LFC_results)
  ensembl_ids <- rownames(as_dataframe)
  
  #Convert biological IDs from genename to entreid and ensembl
  gene_map <- bitr(ensembl_ids, 
                   fromType = "GENENAME", 
                   toType = c("ENTREZID", "ENSEMBL"),
                   OrgDb = org.Sc.sgd.db)
  
  as_dataframe$GENENAME <- rownames(as_dataframe)
  
  dataframe_2 <- merge(as_dataframe, gene_map, by = "GENENAME", all.x = TRUE)
  
  sig_up_genes_ensembl <-  dataframe_2%>%
    filter(padj < 0.05 & log2FoldChange > 1) %>%
    pull(ENSEMBL) %>%
    na.omit() %>%
    unique()
  
  sig_down_genes_ensembl <-  dataframe_2%>%
    filter(padj < 0.05 & log2FoldChange < - 1) %>%
    pull(ENSEMBL) %>%
    na.omit() %>%
    unique()
  
  #SetReadable is from DOSE which is for human data? Probably not needed here...
  #kegg_enrich <- setReadable(kegg_enrich, 
  #          OrgDb = org.Sc.sgd.db, 
  #      keyType = "ENTREZID")
  
  #kegg analysis
  kegg_enrich <- compareCluster(geneClusters = list(Upregulated = sig_up_genes_ensembl,
                                                    Downregulated = sig_down_genes_ensembl),
                                fun = "enrichKEGG",
                                organism = 'sce',
                                pvalueCutoff = 0.05,
                                qvalueCutoff = 0.2)
  return(kegg_enrich)
}

keg_early_to_thin <- kegg_analysis(resLFC_thin_vs_early)
kegg_thin_to_mature <- kegg_analysis(resLFC_mature_vs_thin)
kegg_early_to_mature <- kegg_analysis(resLFC_mature_vs_early)
# Dot plots

kegg_dotplot_1 <- dotplot(keg_early_to_thin, showCategory = 15, title = "KEGG Pathway Enrichment \n Early to Thin")
kegg_dotplot_2 <- dotplot(kegg_thin_to_mature, showCategory = 15, title = "KEGG Pathway Enrichment\nThin to Mature")
kegg_dotplot_3 <- dotplot(kegg_early_to_mature, showCategory = 15, title = "KEGG Pathway Enrichment\nEarly to Mature")

plot_list(kegg_dotplot_1, kegg_dotplot_2, kegg_dotplot_3)

# ====================================================
# Optional: GSEA analysis for all three stages (cuz I was curious about GSEA)
# ====================================================

# Analyis tutorial can be accessed from the link below:
# https://learn.gencore.bio.nyu.edu/rna-seq-analysis/gene-set-enrichment-analysis/

# Make a function that takes in LFC shrinkage results, ontology analysis, and keytype. also does data wrangling to get the proper gene names
gsea_analysis <- function (LFC_results, ont, keytype) {
  gene_list_df <- as.data.frame(LFC_results)
  gene_list_df <- na.omit(gene_list_df)
  
  gene_list <- gene_list_df$log2FoldChange
  
  gene_list_df$NAMES <- rownames(gene_list_df)
  
  names(gene_list) <- gene_list_df$NAMES
  
  gene_list <- sort(gene_list, decreasing = TRUE)
  
  #run GSEA 
  gsea_results <- gseGO(geneList = gene_list,
                        ont = ont,
                        OrgDb = org.Sc.sgd.db, 
                        keyType = keytype,
                        pvalueCutoff = 0.05,
                        pAdjustMethod = "BH",
                        verbose =  TRUE, 
                        eps = 0, #P values less than 1.0 e-10
                        by = "fgsea")
  
  #remove redundant go terms
  #https://github.com/YuLab-SMU/clusterProfiler/issues/28
  gsea_results_2 <- clusterProfiler::simplify(x = gsea_results, 
                                         cutoff = 0.7,
                                         by = "p.adjust",
                                         select_fun = min)
return(gsea_results_2)
}

#Run function and plot enrichment plot for all three stages
# For plotting, check the following tutorial:
# https://yulab-smu.top/biomedical-knowledge-mining-book/enrichplot.html
gsea_mature_vs_early <- gsea_analysis(resLFC_mature_vs_early, "BP", "COMMON")

gsea_plot_mature_vs_early<- gseaplot2(gsea_mature_vs_early, geneSetID = c(1,2,3), 
                                      color = c("#E495A5", "#86B875", "#7DB0DD"),
                                      pvalue_table = FALSE, 
                                      title = "GSEA of Early to Mature Biofilm Formation")


gsea_early_vs_thin <- gsea_analysis(resLFC_thin_vs_early, "BP", "COMMON")

gsea_plot_early_vs_thin <- gseaplot2(gsea_early_vs_thin, geneSetID = c(1:3),
                                         color = c("#E495A5", "#86B875", "#7DB0DD"),
                                         pvalue_table = FALSE, 
                                         title = "GSEA of Early to Thin Biofilm Formation")



gsea_thin_vs_mature <- gsea_analysis(resLFC_mature_vs_thin, "BP", "COMMON")

gsea_plot_mature_vs_thin <- gseaplot2(gsea_thin_vs_mature, geneSetID = c(1:3), 
                                         color = c("#E495A5", "#86B875", "#7DB0DD"),
                                         pvalue_table = FALSE, 
                                         title = "GSEA of Thin to Mature Biofilm Formation")

plot_list(gsea_plot_early_vs_thin, gsea_plot_mature_vs_thin, gsea_plot_mature_vs_early)
