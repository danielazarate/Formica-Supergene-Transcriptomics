#!/bin/env/R

# Script for plotting genome-wide FST calculated in sliding windows and with a window step. 
# FST calculated using VCFTOOLS and chromosomes stitched together using ScaffStitcher.py. 

# Script used in publication - \
# Zarate et al. In Prep. Differential Gene Expression of a supergene regulating \
# social behavior in an ant (Formica neoclara). 

# Authors: Daniela Zarate, PhD - UC Riverside 

install.packages("data.table")
install.packages("dplyr")
install.packages("reshape2")

library(data.table)
library(ggplot2)
library(dplyr)
library(reshape2)

# setwd
setwd("~/Desktop/AntTranscriptomics_project/")

# read in 3 files of data
fst<-fread("fst_pop1_pop2.txt.windowed.weir.fst.continuous",header=T)
fst<-fread("fst_pop1_pop3.txt.windowed.weir.fst.continuous",header=T)
fst<-fread("fst_pop2_pop3.txt.windowed.weir.fst.continuous",header=T)
fst<-fread("fst_pop4_pop5.txt.windowed.weir.fst.continuous",header=T)


# attach file to environment 
attach(fst)

# fix distance to be in Mb
fst$BIN_START=fst$BIN_START/1E6
fst$BIN_END=fst$BIN_END/1E6

# Weir's $FST$ can give values <0 so we reset to 0
fst$WEIGHTED_FST[which(fst$WEIGHTED_FST<0)]=0 

# use ggplot to plot weighted FST by Bin Start Position 
plot1 <- ggplot(fst, aes(BIN_START, WEIGHTED_FST)) + 
  geom_point() + 
  geom_smooth() + theme_bw()

plot1

# Calculate the median of the Weighted FST for each chromosome
axisdf <- fst%>% 
  group_by(CHROM) %>% 
  summarize(center=(max(BIN_START) + min(BIN_START) ) / 2 )

# add the chromosome median as chromosome breaks for the x axis 
plot2 <- plot1 + scale_x_continuous(name = "Genome Length", label = axisdf$CHROM, breaks = axisdf$center, limits=c(6000,184167025)
+ scale_y_continuous(name = "Weighted Fst", limits=c(0,1))) 

plot2

# Set the high-resolution png file 
png(("~/Desktop/AntTranscriptomics_project/Pop4.Pop5.Fst.jpg"),
    width = 26, 
    height = 7, 
    res = 300, 
    units = 'in',)
dev.off()


