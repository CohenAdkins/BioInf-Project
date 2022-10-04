# Script following video tutorial
library(dplyr)
library(tidyverse)
BiocManager::install("GEOquery")
library(GEOquery)
BiocManager::install("DESeq2")
library(DESeq2)

# Check working Directory
getwd()

# Imports Count Data
countData <- read.csv("combinedCountFile.csv", header=TRUE, sep=",", row.names = 1)

# Imports Meta Data
metaData <- read.delim("GSE212377_series_matrix.txt", sep="\t", header=TRUE)
# HERE need to extract and clean up a useful component from this metadata
# HOWEVER it appears that all the data here is homogenous and therefore of no use ??

# Install Homo Sapien genome annotations database for mapping gene names
BiocManager::install("org.Hs.eg.db") 
library("org.Hs.eg.db")
# EntrezID to Symbol for countData
transpCountData <- t(countData) #switches columns with rows (Samples are rows, genes are columns)
mapping <- select(org.Hs.eg.db, keys=colnames(transpCountData), column=c("SYMBOL"), keytype="ENTREZID")

#            Question 1
scaledCountData <- log2(countData) #log2 scale the data
range <- apply(X = scaledCountData, MARGIN = 1, FUN = range)
change <- range[2, ] - range[1, ]
plot(density(na.omit(change))) #Omits genes with no count values


#            Question 2
ddset <- DESeqDataSetFromMatrix( 
  countData = countData, # Here we supply non-normalized count data 
  colData = metaData, # Need to clean metadata first
  design = ~mutation_status ) # Need to select experimental variable to `design` 
deseq_object <- DESeq(ddset)

