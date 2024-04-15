#!/bin/env/R 
# Script for plotting Volcano Plots from DESEQ2 isoform transcript count data. 
# Author: Daniela Zarate, PhD.
# Adapted from DESEQ2 / TRINITY R scripts. 
# Built and run under R.4.3.3 - Angel Food Cake
#_______________________________________________________________________________#

# Load some packages using BiocManager
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.18")

BiocManager::install(c("DESeq2"))
BiocManager::install(c("edgeR"))

library(edgeR)
library(DESeq2)
library(gridExtra) 
library(cowplot)
library(ggplot2)
library(ggpubr)
#_______________________________________________________________________________#

# Set working directory
setwd("~/Desktop/AntTranscriptomics_project/new_isoform_volcanoes/")

#_______________________________________________________________________________#

# Read in the genes counts matrix 
data = read.table("Trinity.v2.isoform.counts.matrix", header=T, row.names=1, com='')

# Pull out all the columns with sm_sm and then all those with sm_sp: 
mmmp_col_ordering = c(1,2,3,7,8,9,17,18,23,28,29,37,4,6,5,10,13,14,
                      19,20,24,25,30,31,32,35,38,39)

# For sm_sm vs sp_sp 
mmpp_col_ordering = c(1,2,3,7,8,9,17,18,23,28,29,37,11,12,15,16,21,
                      22,26,27,33,34,36,40)

# For sm_sp vs sp_sp 
mppp_col_ordering = c(4,6,5,10,13,14,19,20,24,25,30,31,32,35,38,39,
                      11,12,15,16,21,22,26,27,33,34,36,40)

# Create a matrix, pulling out those columns that we indicated above:
mmmp_rnaseqMatrix = data[,mmmp_col_ordering]
mmpp_rnaseqMatrix = data[,mmpp_col_ordering]
mppp_rnaseqMatrix = data[,mppp_col_ordering]

# Round the matrix numbers to nearest whole number 
mmmp_rnaseqMatrix = round(mmmp_rnaseqMatrix)
mmpp_rnaseqMatrix = round(mmpp_rnaseqMatrix)
mppp_rnaseqMatrix = round(mppp_rnaseqMatrix)

# Generates count per million data per matrix row, if the value per cell is greater \
# than one, keep it and if the rowsum is greater than or equal to 2, keep the row. 
mmmp_rnaseqMatrix = mmmp_rnaseqMatrix[rowSums(cpm(mmmp_rnaseqMatrix) > 1) >= 2,]
mmpp_rnaseqMatrix = mmpp_rnaseqMatrix[rowSums(cpm(mmpp_rnaseqMatrix) > 1) >= 2,]
mppp_rnaseqMatrix = mppp_rnaseqMatrix[rowSums(cpm(mppp_rnaseqMatrix) > 1) >= 2,]

# Create a list of the replicates, 12 of sm_sm and 16 of sm_sp 
mmmp_conditions = data.frame(conditions=factor(c(rep("sm_sm", 12), rep("sm_sp", 16))))

# Create a list of the replicates, 12 of sm_sm and 12 of sp_sp 
mmpp_conditions = data.frame(conditions=factor(c(rep("sm_sm", 12), rep("sp_sp", 12))))

# Create a list of the replicates, 12 of sm_sp and 16 of sp_sp 
mppp_conditions = data.frame(conditions=factor(c(rep("sm_sp", 16), rep("sp_sp", 12))))

# replace conditions with colnames in the matrix (did it work?)
rownames(mmmp_conditions) = colnames(mmmp_rnaseqMatrix)
rownames(mmpp_conditions) = colnames(mmpp_rnaseqMatrix)
rownames(mppp_conditions) = colnames(mppp_rnaseqMatrix)

# Covert counts to integer mode 
mmmp_ddsFullCountTable <- DESeqDataSetFromMatrix(
  countData = mmmp_rnaseqMatrix,
  colData = mmmp_conditions,
  design = ~ conditions)


mmpp_ddsFullCountTable <- DESeqDataSetFromMatrix(
  countData = mmpp_rnaseqMatrix,
  colData = mmpp_conditions,
  design = ~ conditions)

mppp_ddsFullCountTable <- DESeqDataSetFromMatrix(
  countData = mppp_rnaseqMatrix,
  colData = mppp_conditions,
  design = ~ conditions)

# Run DESeq2: estimating size factors, dispersions, gene-wide dispersion estimates \
# mean-dispersion relationship, final dispersion estimates, fitting model and testing etc. 
mmmp_dds = DESeq(mmmp_ddsFullCountTable)
mmpp_dds = DESeq(mmpp_ddsFullCountTable)
mppp_dds = DESeq(mppp_ddsFullCountTable)

# Create a vector with names: 
mmmp_contrast=c("conditions","sm_sm","sm_sp")
mmpp_contrast=c("conditions","sm_sm","sp_sp")
mppp_contrast=c("conditions","sm_sp","sp_sp")

# create an object combining both contrast and dds:
mmmp_res = results(mmmp_dds, mmmp_contrast)
mmpp_res = results(mmpp_dds, mmpp_contrast)
mppp_res = results(mppp_dds, mppp_contrast)

# Calculate beaseMeans for both conditions 
mmmp_baseMeanA <- rowMeans(counts(mmmp_dds, normalized=TRUE)[,colData(mmmp_dds)$conditions == "sm_sm"])
mmmp_baseMeanB <- rowMeans(counts(mmmp_dds, normalized=TRUE)[,colData(mmmp_dds)$conditions == "sm_sp"])

mmpp_baseMeanA <- rowMeans(counts(mmpp_dds, normalized=TRUE)[,colData(mmpp_dds)$conditions == "sm_sm"])
mmpp_baseMeanB <- rowMeans(counts(mmpp_dds, normalized=TRUE)[,colData(mmpp_dds)$conditions == "sp_sp"])

mppp_baseMeanA <- rowMeans(counts(mppp_dds, normalized=TRUE)[,colData(mppp_dds)$conditions == "sm_sp"])
mppp_baseMeanB <- rowMeans(counts(mppp_dds, normalized=TRUE)[,colData(mppp_dds)$conditions == "sp_sp"])


# combine baseMean data from both genotypes
mmmp_res = cbind(mmmp_baseMeanA, mmmp_baseMeanB, as.data.frame(mmmp_res))
mmmp_res = cbind(sampleA="sm_sm", sampleB="sm_sp", as.data.frame(mmmp_res))

mmpp_res = cbind(mmpp_baseMeanA, mmpp_baseMeanB, as.data.frame(mmpp_res))
mmpp_res = cbind(sampleA="sm_sm", sampleB="sp_sp", as.data.frame(mmpp_res))

mppp_res = cbind(mppp_baseMeanA, mppp_baseMeanB, as.data.frame(mppp_res))
mppp_res = cbind(sampleA="sm_sp", sampleB="sp_sp", as.data.frame(mppp_res))


#  Order by p-value (smallest to lowest, if its < 1)
mmmp_res$padj[is.na(mmmp_res$padj)]  <- 1
mmmp_res = as.data.frame(mmmp_res[order(mmmp_res$pvalue),])

mmpp_res$padj[is.na(mmpp_res$padj)]  <- 1
mmpp_res = as.data.frame(mmpp_res[order(mmpp_res$pvalue),])

mppp_res$padj[is.na(mppp_res$padj)]  <- 1
mppp_res = as.data.frame(mppp_res[order(mppp_res$pvalue),])


# Add column to res df that prints sm_sm if sm_sm biased or sm_sp if sm_sp biased
mmmp_res$GENOBIAS <- ifelse(mmmp_res$mmmp_baseMeanA>mmmp_res$mmmp_baseMeanB, "sm_sm", "sm_sp")
mmpp_res$GENOBIAS <- ifelse(mmpp_res$mmpp_baseMeanA>mmpp_res$mmpp_baseMeanB, "sm_sm", "sp_sp")
mppp_res$GENOBIAS <- ifelse(mppp_res$mppp_baseMeanA>mppp_res$mppp_baseMeanB, "sm_sp", "sp_sp")

# PADJ top 1% values:
mmmp_padj=3.385693e-13
mmpp_padj=1.526518e-58
mppp_padj=8.753579e-45

#_______________________________________________________________________________#
                           # Volcano Plot Functions: 
#_______________________________________________________________________________#


## MM-MP 

mmmp_res$title <- "MM vs MP"
mm_mp_plot <- ggplot(mmmp_res, aes(x=log2FoldChange, y=-1*log10(padj))) +
  geom_point(size = 4.5, color="black", stroke = 0.1) +
  xlab("log FC") + theme_bw() + geom_point(aes(colour = 
    ifelse(3.385693e-13<padj & padj<=0.05 & GENOBIAS=="sm_sm","olivedrab3" ,
    ifelse(3.385693e-13<padj & padj<=0.05 & GENOBIAS=="sm_sp", "goldenrod",
    ifelse(padj<=3.385693e-13 & GENOBIAS=="sm_sm", "red",
    ifelse(padj<=3.385693e-13 & GENOBIAS=="sm_sp", "purple", "black")))))) +
  scale_color_identity(name = "DEGs",
      breaks = c("black", "goldenrod", "olivedrab3", "purple", "red"),
    labels = c("ns", "MP DEGs", "MM DEGs", "Top 1% MP DEGs", "Top 1% MM DEGs"),
       guide="legend") + 
  facet_grid(. ~ title) + theme(strip.text = element_text(size = 15)) +
   theme(legend.position = "bottom", 
         legend.title = element_text(size=20),
         legend.text = element_text(size=20)) +
   theme(axis.title = element_text(size=20)) +
  guides(colour = guide_legend(override.aes = list(size=5)))
mm_mp_plot
  
#_______________________________________________________________________________#

## MM-PP

mmpp_res$title <- "MM vs PP"
mm_pp_plot <- ggplot(mmpp_res, aes(x=log2FoldChange, y=-1*log10(padj))) +
  geom_point(size = 4.5, color="black", stroke = 0.1) +
  xlab("log FC") + theme_bw() + geom_point(aes(colour = 
    ifelse(1.526518e-58<padj & padj<=0.05 & GENOBIAS=="sm_sm","olivedrab3" ,
    ifelse(1.526518e-58<padj & padj<=0.05 & GENOBIAS=="sp_sp", "blue",
    ifelse(padj<=1.526518e-58 & GENOBIAS=="sm_sm", "red",
    ifelse(padj<=1.526518e-58 & GENOBIAS=="sp_sp", "cyan", "black")))))) +
   scale_color_identity(name = "DEGs", breaks = c("black", "blue", "olivedrab3", "red"),
    labels = c("ns", "PP DEGs", "MM DEGs", "Top 1% MM DEGs"), guide="legend") +
  facet_grid(. ~ title) + theme(strip.text = element_text(size = 15)) +
  theme(legend.position = "bottom", 
        legend.title = element_text(size=20),
        legend.text = element_text(size=20)) +
  theme(axis.title = element_text(size=20)) + 
   guides(colour = guide_legend(override.aes = list(size=5)))
mm_pp_plot

#_______________________________________________________________________________#


## MP-PP

mppp_res$title <- "MP vs PP"
mp_pp_plot <- ggplot(mppp_res, aes(x=log2FoldChange, y=-1*log10(padj))) +
  geom_point(size = 4.5, color="black", stroke = 0.1) +
  xlab("log FC") + theme_bw() + geom_point(aes(color = 
   ifelse(8.753579e-45<padj & padj<=0.05 & GENOBIAS=="sm_sp","goldenrod" ,
   ifelse(8.753579e-45<padj & padj<=0.05 & GENOBIAS=="sp_sp", "blue",
   ifelse(padj<=8.753579e-45 & GENOBIAS=="sm_sp", "purple",
   ifelse(padj<=8.753579e-45 & GENOBIAS=="sp_sp", "cyan", "black")))))) +
   scale_color_identity(name = "DEGs", breaks = c("black", "blue", "purple", "goldenrod"),
    labels = c("ns", "PP DEGs", "Top 1% MP DEGs", "MP DEGs"), guide="legend") +
   facet_grid(. ~ title) + theme(strip.text = element_text(size = 15)) +
  theme(legend.position = "bottom", 
        legend.title = element_text(size=20),
        legend.text = element_text(size=20)) +
  theme(axis.title = element_text(size=20)) + 
  guides(colour = guide_legend(override.aes = list(size=5)))
mp_pp_plot


#_______________________________________________________________________________#
                             # Volcano GridPlot # 
#_______________________________________________________________________________#

# Plot the three Volcanoes side by side without legends:
plot_grid(mm_mp_plot + theme(legend.position="none"),
          mm_pp_plot + theme(legend.position="none"),
          mp_pp_plot + theme(legend.position="none"),
          labels=c("A", "B", "C"), 
          ncol = 3, nrow = 1, align = "V") 

# Set the high-resolution png file 
png(("Isoforms.Volcanoes.Genotype.Pairwise.Comparisons.Big.png"),
    width = 15, 
    height = 15, 
    res = 300, 
    units = 'in', 
)
dev.off()

#_______________________________________________________________________________#



# Subset all padj < 0.05 for each PC 
mmmp_new_res <- subset(mmmp_res, mmmp_res$padj<=0.05)
nrow(mmmp_new_res)

mmpp_new_res <- subset(mmpp_res, mmpp_res$padj<=0.05)
nrow(mmpp_new_res)

mppp_new_res <- subset(mppp_res, mppp_res$padj<=0.05)
nrow(mppp_new_res)


# Subset all mm-biased from the padj < 0.05 subsection for each PC: 
mmmp_new_res2 <- subset(mmmp_new_res, mmmp_new_res$GENOBIAS=="sm_sm")
nrow(mmmp_new_res2)

mmpp_new_res2 <- subset(mmpp_new_res, mmpp_new_res$GENOBIAS=="sm_sm")
nrow(mmpp_new_res2)

mppp_new_res2 <- subset(mppp_new_res, mppp_new_res$GENOBIAS=="sm_sp")
nrow(mppp_new_res2)


# Subset all mm-biased from entire dataset 
mmmp_new_res3 <- subset(mmmp_res, mmmp_res$GENOBIAS=="sm_sm")
nrow(mmmp_new_res3)

mmpp_new_res3 <- subset(mmpp_res, mmpp_res$GENOBIAS=="sm_sm")
nrow(mmpp_new_res3)


mppp_new_res3 <- subset(mppp_res, mppp_res$GENOBIAS=="sm_sp")
nrow(mppp_new_res3)





#
write.table(res, file='Trinity.v2.isoform.counts.matrix.sm_sm_vs_sm_sp.DESeq2.DE_results', sep='        ', quote=FALSE)

write.table(rnaseqMatrix, file='Trinity.v2.isoform.counts.matrix.sm_sm_vs_sm_sp.DESeq2.count_matrix', sep='     ', quote=FALSE)

