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
install.packages("gprofiler2")
install.packages("ClusterR")
install.packages("cluster")
install.packages("cluster")
install.packages("factoextra")
library("ggplot2")
library(M3C)
library("org.Hs.eg.db")
library(dplyr)
library(tidyverse)
library(GEOquery)
library("DESeq2")
library(devtools)
library(ComplexHeatmap)
library(gprofiler2)
library(ClusterR)
library(cluster)
library(ggpubr)
library(factoextra)

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

#Second way of organizing metadata automatically
gse <- getGEO(GEO = 'GSE212377' , GSEMatrix = TRUE)
metaData2 <- pData(phenoData(gse[[1]]))
metaData2 <- select(metaData2, c(1,2))
#Picking out what columns we want from metadata automatically
metaData2<-metaData2 %>%
  select(1,2) %>%
  rename(Histology = title) %>%  #renaming columns
  rename(gene = geo_accession)   #renaming columns
#Pulling from countData, putting the counts in their own columns as well as the samples
count.long <- countData %>%
  gather(key = 'samples', value = 'Count', -EntrezID)
#combining two matrixes together and matching data by its GSM value
count.long <- count.long %>%
  left_join(., metaData2, by = c("samples"="gene"))

count.long$Count <- log10(count.long$Count) #logscaling counts
count.long$EntrezID <- as.character(count.long$EntrezID) #EntrezId to characters, heatmap needs this.
#picking grade 1-3 to use in plots
Interest <- c('MENI-TB100: WHO Grade-2 histology', 'MENI-TB015: WHO Grade-3 histology','MENI-TB097: WHO Grade-1 histology')

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

count.long %>% 
  filter(Histology  %in% Interest) %>% #filters for our hitologies
  ggplot(., aes(x = Count, fill = Histology))+geom_density() #plots our density

#Tsne plot
tsne(countData,colvec=c('gold'))


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
count.long %>% 
  filter(Histology %in% Interest) %>% #filters for our histologies
  ggplot(., aes(x = EntrezID , y =Histology, fill = Count))+ #plots a heatmap
  geom_tile()

#            Question 5 - gProfiler2 - Natalie
gostres <- gost(query = c("X:1000:1000000", "rs17396340", "GO:0005005", "ENSG00000156103", "NLRP1"), 
                organism = "hsapiens", ordered_query = FALSE, 
                multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                measure_underrepresentation = FALSE, evcodes = FALSE, 
                user_threshold = 0.05, correction_method = "g_SCS", 
                domain_scope = "annotated", custom_bg = NULL, 
                numeric_ns = "", sources = NULL, as_short_link = FALSE)
names(gostres)
head(gostres$result, 3)
names(gostres$metaData)
gostres_link <- gost(query = c("X:1000:1000000", "rs17396340", "GO:0005005", "ENSG00000156103", "NLRP1"), 
                     as_short_link = TRUE)
gostplot(gostres, capped = TRUE, interactive = TRUE)







#                     PROJECT 3

#            Question 2a
# Gets variance for each gene
allVariance <- as.data.frame(apply(scaledCountData, 1, FUN = var))
colnames(allVariance) <- c("Variance")

# Sorts gene variance in descending order
sortedVariance <- as.data.frame(allVariance[order(-allVariance$Variance), , drop = FALSE])

# 5000 Most Variable Genes
mostVar5000 <- merge(scaledCountData, head(sortedVariance, 5000), by=0)
row.names(mostVar5000) <- c(mostVar5000[,1]) # Sets row names
mostVar5000$Row.names <- NULL # Deletes Row.names column, no longer needed
mostVar5000$Variance <- NULL # Deletes Variance column, no longer needed

# 10 Most Variable Genes
mostVar10 <- merge(scaledCountData, head(sortedVariance, 10), by=0)
row.names(mostVar10) <- c(mostVar10[,1]) # Sets row names
mostVar10$Row.names <- NULL # Deletes Row.names column, no longer needed
mostVar10$Variance <- NULL # Deletes Variance column, no longer needed

# 100 Most Variable Genes
mostVar100 <- merge(scaledCountData, head(sortedVariance, 100), by=0)
row.names(mostVar100) <- c(mostVar100[,1]) # Sets row names
mostVar100$Row.names <- NULL # Deletes Row.names column, no longer needed
mostVar100$Variance <- NULL # Deletes Variance column, no longer needed

# 1000 Most Variable Genes
mostVar1000 <- merge(scaledCountData, head(sortedVariance, 1000), by=0)
row.names(mostVar1000) <- c(mostVar1000[,1]) # Sets row names
mostVar1000$Row.names <- NULL # Deletes Row.names column, no longer needed
mostVar1000$Variance <- NULL # Deletes Variance column, no longer needed

# 10000 Most Variable Genes
mostVar10000 <- merge(scaledCountData, head(sortedVariance, 10000), by=0)
row.names(mostVar10000) <- c(mostVar10000[,1]) # Sets row names
mostVar10000$Row.names <- NULL # Deletes Row.names column, no longer needed
mostVar10000$Variance <- NULL # Deletes Variance column, no longer needed


#            Question 2b-e Natalie

# 10
df10 <- as.data.frame(mostVar10)
hc10 <- hclust(dist(df10), method = "average")
plot(hc10)
rect.hclust(hc10, k = 3, border = 2:4)

# 100
df100 <- as.data.frame(mostVar100)
hc100 <- hclust(dist(df100), method = "average")
plot(hc100)
rect.hclust(hc100, k = 2, border = 2:3)

# 1000
df1000 <- as.data.frame(mostVar1000)
hc1000 <- hclust(dist(df1000), method = "average")
plot(hc1000, labels = FALSE)
rect.hclust(hc1000, k = 2, border = 2:3)

# 5000
df5000 <- as.data.frame(mostVar5000)
hc5000 <- hclust(dist(df5000), method = "average")
plot(hc5000, labels = FALSE)
rect.hclust(hc5000, k = 2, border = 2:3)

# 10000
df10000 <- as.data.frame(mostVar10000)
hc10000 <- hclust(dist(df10000), method = "average")
plot(hc10000, labels = FALSE)
rect.hclust(hc10000, k = 2, border = 2:3)


            #Question 2b- Cohen
set.seed(123)
#10
res.km10 <- kmeans(scale(t(mostVar10)[,-20]),5, nstart = 25)
fviz_cluster(res.km10, t(mostVar10))

#100
res.km100 <- kmeans(scale(t(mostVar100)[,-20]),5, nstart = 25)
fviz_cluster(res.km100, t(mostVar100))

#1000
res.km1000 <- kmeans(scale(t(mostVar1000)[,-20]),5, nstart = 25)
fviz_cluster(res.km1000, t(mostVar1000))

#5000
res.km5000 <- kmeans(scale(t(mostVar5000)[,-20]),3, nstart = 25)
fviz_cluster(res.km5000, t(mostVar5000))

#            Question 2b-e Avi
# PAM Plot 5000 Genes
res.pam <- pam(t(mostVar5000), 2)
fviz_cluster(res.pam, main = "PAM Cluster Plot", labelsize = 0, xlab = FALSE, ylab = FALSE)

# PAM Plot 10 Genes
res.pam <- pam(t(mostVar10), 3)
fviz_cluster(res.pam, main = "PAM Cluster Plot", labelsize = 0, xlab = FALSE, ylab = FALSE)

# PAM Plot 100 Genes
res.pam <- pam(t(mostVar100), 2)
fviz_cluster(res.pam, main = "PAM Cluster Plot", labelsize = 0, xlab = FALSE, ylab = FALSE)

# PAM Plot 1000 Genes
res.pam <- pam(t(mostVar1000), 3)
fviz_cluster(res.pam, main = "PAM Cluster Plot", labelsize = 0, xlab = FALSE, ylab = FALSE)

# PAM Plot 10000 Genes
res.pam <- pam(t(mostVar10000), 2)
fviz_cluster(res.pam, main = "PAM Cluster Plot", labelsize = 0, xlab = FALSE, ylab = FALSE)

# Alluvial Plot(IN PROGRESS)


#            Question 2b-e Gabriel

# Question  3a
#Heatmap of different places in the 5000 differently expressed genes
# 2 heatmaps of two different places

par(mfrow = c(1,2))
image(t(res.km)[,nrow(res.km):1], yaxt = "n",main = "orginal data")
set.seed(1234)
dataMatrix <- as.matrix(mostVariable)[sample(1:12),]
kmeansObj2 <- kmeans(dataMatrix, centers= 3)
par(mfrow = c(1,2), mar = c(2,4,.1,.1))
image(t(dataMatrix)[,nrow(dataMatrix):1], yaxt = "n")
image(t(dataMatrix)[, order(kmeansObj2$cluster)], yaxt = "n")











