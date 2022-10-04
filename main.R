# Script following video tutorial
library(dplyr)
library(tidyverse)
BiocManager::install("GEOquery")
library(GEOquery)
BiocManager::install("DESeq2")
library(DESeq2)

# Check working Directory
getwd()

# Imports Count Data and Re-formats
countData <- read.csv("combinedCountFile.csv", header=TRUE, sep=",")
colnames(countData) <- c(sapply(colnames(countData), function(x){substr(x,0, 10)}))

# Imports and Pre-processes Meta Data
metaData <- read.delim("RefinedMetaData.txt", sep="\t", header=FALSE)
table(sapply(metaData[1, ], function(x){substr(x, 13, nchar(x))}))
metaData[3, ] <- sapply(metaData[1, ], function(x){substr(x, 13, nchar(x))})
colnames(metaData) <- metaData[2,]
metaData <- t(metaData)
colnames(metaData) <- c('ExtHistology', 'SampleID', 'Histology')
metaData <- metaData[-1,]
metaData

# Install Homo Sapien genome annotations database for mapping gene names
BiocManager::install("org.Hs.eg.db") 
library("org.Hs.eg.db")
# EntrezID to Symbol for countData
transpCountData <- t(countData) #switches columns with rows (Samples are rows, genes are columns)
mapping <- select(org.Hs.eg.db, keys=transpCountData[1,]), column=c("SYMBOL"), keytype="ENTREZID")
transpCountData <- t(transpCountData)
table(transpCountData)
table(mapping)
mapping
transpCountData[1:5,1:5]
colnames(mapping) <- c('EntrezID', 'HugoID')

transpCountData <- transpCountData %>%
  left_join(., mapping, by = c("EntrezID" = "EntrezID"))

# x <- merge(transpCountData, mapping, by = c("EntrezID")
head(x)
countData[1:5,1:5]

#            Question 1
scaledCountData <- log2(countData) #log2 scale the data
range <- apply(X = scaledCountData, MARGIN = 1, FUN = range)
change <- range[2, ] - range[1, ]
plot(density(na.omit(change))) #Omits genes with no count values


#            Question 2
ddset <- DESeqDataSetFromMatrix( 
  countData = countData, # Here we supply non-normalized count data 
  colData = metaData, # Need to clean metadata first
  design = ~Histology ) # Need to select experimental variable to `design` 
deseq_object <- DESeq(ddset)

