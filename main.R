# Install Packages
install.packages("dplyr")
BiocManager::install("DESeq2")
install.packages("tidyverse")
BiocManager::install("GEOquery")
BiocManager::install("org.Hs.eg.db") 
BiocManager::install("M3C")
install.packages("ggplot2")
install.packages("devtools")
install.packages("ComplexHeatmap")
library("ggplot2")
library(M3C)
library("org.Hs.eg.db")
library(dplyr)
library(tidyverse)
library(GEOquery)
library("DESeq2")
library(devtools)
library(ComplexHeatmap)

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

# EntrezID to Symbol for countData
transpCountData <- t(countData) #switches columns with rows (Samples are rows, genes are columns)
mapping <- select(org.Hs.eg.db, keys=transpCountData[1,], column=c("SYMBOL"), keytype="ENTREZID")
transpCountData <- t(transpCountData)
colnames(mapping) <- c('EntrezID', 'HugoID')
countData <- merge(transpCountData, mapping, by = c("EntrezID")) # Adds HugoID to countData
countData <- countData %>%
  select(HugoID, everything()) # Moves HugoID column to front
countData$EntrezID <- NULL # Deletes EntrezID column
countData <- na.omit(countData) # Gets rid of NA values
rownames(countData) <- countData$HugoID # makes rownames HugoID
countData$HugoID <- NULL # Deletes HugoID because its now rowname

#            Question 1
scaledCountData <- log2(countData) #log2 scale the data
range <- apply(X = scaledCountData, MARGIN = 1, FUN = range)
change <- range[2, ] - range[1, ]
plot(density(na.omit(change))) #Omits genes with no count values


#            Question 2
metaData[order(row.names(metaData)), ] # Sorts metaData
countData <- countData[,order(colnames(countData))] # Sorts countData
ddset <- DESeqDataSetFromMatrix( 
  countData = countData, # Here we supply non-normalized count data 
  colData = metaData, # Need to clean metadata first
  design = ~ Histology) # Need to select experimental variable to `design` 
deseq_object <- DESeq(ddset)

vsd <- vst(ddset, blind=FALSE)
plotPCA(vsd, intgroup=c("Histology"))

#            Question 3
if (!("DESeq2" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  BiocManager::install("DESeq2", update = FALSE)
}
if (!("EnhancedVolcano" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  BiocManager::install("EnhancedVolcano", update = FALSE)
}
if (!("apeglm" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  BiocManager::install("apeglm", update = FALSE)
}

# Attach the DESeq2 library
library(DESeq2)
# Attach the ggplot2 library for plotting
library(ggplot2)
# We will need this so we can use the pipe: %>%
library(magrittr)

set.seed(12345)
expression_df <- vsd

deseq_results <- results(deseq_object)

# this is of class DESeqResults -- we want a data frame
deseq_df <- deseq_results %>%
  # make into data.frame
  as.data.frame() %>%
  # the gene names are row names -- let's make them a column for easy display
  tibble::rownames_to_column("Gene") %>%
  # add a column for significance threshold results
  dplyr::mutate(threshold = padj < 0.05) %>%
  # sort by statistic -- the highest values will be genes with
  # higher expression in RPL10 mutated samples
  dplyr::arrange(dplyr::desc(log2FoldChange))

# We'll assign this as `volcano_plot`
volcano_plot <- EnhancedVolcano::EnhancedVolcano(
  deseq_df,
  lab = deseq_df$Gene,
  x = "log2FoldChange",
  y = "padj",
  pCutoff = 0.01 # Loosen the cutoff since we supplied corrected p-values
)

# Print out plot here
volcano_plot

#            Question 4
Heatmap(metaData)










