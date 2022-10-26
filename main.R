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
install.packages("DataCombine")
install.packages("ggalluvial")
library("ggalluvial")
library("DataCombine")
install.packages("dendextend")
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
library(dendextend)

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

# Cleans metaData 10/26
histology <- as.data.frame(metaData)
histology$SampleID <- NULL
histology$ExtHistology <- NULL
Replaces <- data.frame(from = c("Not Classified", "WHO Grade-1 histology", "WHO Grade-2 histology", "WHO Grade-3 histology"), to = c("NA", "1", "2", "3"))
histology <- FindReplace(histology, Var = "Histology", Replaces, from = "from", to = "to", exact = TRUE, vector = FALSE)

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







#                     ASSIGNMENT 3

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
df10 <- as.data.frame(t(mostVar10))
hc10 <- hclust(dist(df10), method = "average")
plot(hc10)
rect.hclust(hc10, k = 4, border = 2:5)

# 100
df100 <- as.data.frame(t(mostVar100))
hc100 <- hclust(dist(df100), method = "average")
plot(hc100)
rect.hclust(hc100, k = 6, border = 2:7)

# 1000
df1000 <- as.data.frame(t(mostVar1000))
hc1000 <- hclust(dist(df1000), method = "average")
plot(hc1000)
rect.hclust(hc1000, k = 2, border = 2:3)

# 5000
df5000 <- as.data.frame(t(mostVar5000))
hc5000 <- hclust(dist(df5000), method = "average")
plot(hc5000)
rect.hclust(hc5000, k = 2, border = 2:3)

# 10000
df10000 <- as.data.frame(t(mostVar10000))
hc10000 <- hclust(dist(df10000), method = "average")
plot(hc10000)
rect.hclust(hc10000, k = 2, border = 2:3)

# Hclust Cluster Membership Vector for 5000 genes
HCclusMem5000 <- as.data.frame(cutree(hc5000, k = 2))
colnames(HCclusMem5000) <- c("HCCluster5000")

# Create a vector like 1 1 1 1 2 1 1 1 1 2 1 1 ... etc. based on the cluster group
hcGroups10 <- cutree(hc10, k = 4)
hcGroups100 <- cutree(hc100, k = 6)
hcGroups1000 <- cutree(hc1000, k = 2)
hcGroups5000 <- cutree(hc5000, k = 2)
hcGroups10000 <- cutree(hc10000, k = 2)

# Create a matrix of cluster data for alluvial
hcAlluvialMatrix <- matrix(hcGroups10, nrow = 102, ncol = 6)
hcAlluvialMatrix[,2] <- matrix(hcGroups100)
hcAlluvialMatrix[,3] <- matrix(hcGroups1000)
hcAlluvialMatrix[,4] <- matrix(hcGroups5000)
hcAlluvialMatrix[,5] <- matrix(hcGroups10000)

# Add WHO Histology to hcAlluvial
hcAlluvialMatrix[,6] <- metaData[,3]

# Alluvial
head(as.data.frame(hcAlluvialMatrix), n = 102)
is_alluvia_form(as.data.frame(hcAlluvialMatrix), axes = 1:5, silent = TRUE)
ggplot(as.data.frame(hcAlluvialMatrix),
       aes(axis1 = V1, axis2 = V2, axis3 = V3, axis4 = V4, axis5 = V5)) +
  geom_alluvium(aes(fill = V6), width = 1/12) +
  geom_stratum(width = 1/12, fill = "black", color = "grey") +
  geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("10 genes", "100 genes", "1000 genes", "5000 genes", "10000 genes"), expand = c(.05, .05)) +
  scale_fill_brewer(type = "qual", palette = "Set1") +
  ggtitle("Hierarchical Clustering by Different Gene Amounts")



            #Question 2b-e Cohen
set.seed(123)
#10
res.km10 <- kmeans(scale(t(mostVar10)[,-20]),3, nstart = 25)
fviz_cluster(res.km10, t(mostVar10), main = "K-Means plot", labelsize = 0)

#100
res.km100 <- kmeans(scale(t(mostVar100)[,-20]),4, nstart = 25)
fviz_cluster(res.km100, t(mostVar100), main = "K-Means plot", labelsize = 0)

#1000
res.km1000 <- kmeans(scale(t(mostVar1000)[,-20]),3, nstart = 25)
fviz_cluster(res.km1000, t(mostVar1000), main = "K-Means plot", labelsize = 0)

#5000
res.km5000 <- kmeans(scale(t(mostVar5000)[,-20]),3, nstart = 25)
fviz_cluster(res.km5000, t(mostVar5000), main = "K-Means plot", labelsize = 0)
KMclusMem5000 <- as.data.frame(res.km5000$cluster)
colnames(KMclusMem5000) <- c("KMCluster")

#10,000
res.km10000 <- kmeans(scale(t(mostVar10000)[,-20]),3, nstart = 25)
fviz_cluster(res.km10000, t(mostVar10000), main = "K-Means plot", labelsize = 0)


#            Question 2b-e Avi
# PAM Plot 5000 Genes
res.pam <- pam(t(mostVar5000), 2)
fviz_cluster(res.pam, main = "PAM Cluster Plot", labelsize = 0)
clusMem5000 <- as.data.frame(res.pam$clustering)
colnames(clusMem5000) <- c("PAMCluster5000")

# PAM Plot 10 Genes
res.pam <- pam(t(mostVar10), 3)
fviz_cluster(res.pam, main = "PAM Cluster Plot", labelsize = 0)
clusMem10 <- as.data.frame(res.pam$clustering)
colnames(clusMem10) <- c("PAMCluster10")

# PAM Plot 100 Genes
res.pam <- pam(t(mostVar100), 4)
fviz_cluster(res.pam, main = "PAM Cluster Plot", labelsize = 0)
clusMem100 <- as.data.frame(res.pam$clustering)
colnames(clusMem100) <- c("PAMCluster100")

# PAM Plot 1000 Genes
res.pam <- pam(t(mostVar1000), 3)
fviz_cluster(res.pam, main = "PAM Cluster Plot", labelsize = 0)
clusMem1000 <- as.data.frame(res.pam$clustering)
colnames(clusMem1000) <- c("PAMCluster1000")

# PAM Plot 10000 Genes
res.pam <- pam(t(mostVar10000), 2)
fviz_cluster(res.pam, main = "PAM Cluster Plot", labelsize = 0)
clusMem10000 <- as.data.frame(res.pam$clustering)
colnames(clusMem10000) <- c("PAMCluster10000")

# Alluvial Plot(IN PROGRESS)
PAMAlluvDF <- merge(histology, clusMem10, by=0)
row.names(PAMAlluvDF) <- c(PAMAlluvDF$Row.names)
PAMAlluvDF$Row.names <- NULL
PAMAlluvDF <- merge(PAMAlluvDF, clusMem100, by=0)
row.names(PAMAlluvDF) <- c(PAMAlluvDF$Row.names)
PAMAlluvDF$Row.names <- NULL
PAMAlluvDF <- merge(PAMAlluvDF, clusMem1000, by=0)
row.names(PAMAlluvDF) <- c(PAMAlluvDF$Row.names)
PAMAlluvDF$Row.names <- NULL
PAMAlluvDF <- merge(PAMAlluvDF, clusMem5000, by=0)
row.names(PAMAlluvDF) <- c(PAMAlluvDF$Row.names)
PAMAlluvDF$Row.names <- NULL
PAMAlluvDF <- merge(PAMAlluvDF, clusMem10000, by=0)
row.names(PAMAlluvDF) <- c(PAMAlluvDF$Row.names)


is_alluvia_form(as.data.frame(PAMAlluvDF), axes = 1:7, silent = TRUE)
c <- as.data.frame(Titanic)

ggplot(PAMAlluvDF,
       aes(y = Row.names,
           axis1 = Histology, axis2 = PAMCluster10, axis3 = PAMCluster100, axis4 = PAMCluster1000, axis5 = PAMCluster5000, axis6 = PAMCluster10000)) +
  geom_alluvium(aes(fill = Row.names),
                width = 0, knot.pos = 0, reverse = FALSE) +
  guides(fill = FALSE) +
  geom_stratum(width = 1/8, reverse = FALSE) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)),
            reverse = FALSE) +
  scale_x_continuous(breaks = 1:6, labels = c("Histology", "10 Genes", "100 Genes", "1000 Genes", "5000 Genes", "10000 Genes")) +
  ggtitle("PAM Alluvial Diagram")

#       Question 4
#       a)
# Chi squared test of PAM 100 Genes(4 clusters) vs Histology
x <- merge(clusMem100, histology, by=0)
row.names(x) <- c(x$Row.names)
x$Row.names <- NULL
clus100vsHist <- table(x$PAMCluster100, x$Histology)
chisq.test(clus100vsHist)

# Chi squared test of K-means 5000 Genes(3 clusters) vs Histology
x1 <- merge(KMclusMem5000, histology, by=0)
row.names(x1) <- c(x1$Row.names)
x1$Row.names <- NULL
KMclus5000vsHist <- table(x1$KMCluster, x1$Histology)
chisq.test(KMclus5000vsHist)

# Chi squared test of Hclust vs Histology
x2 <- merge(HCclusMem5000, histology, by=0)
row.names(x2) <- c(x2$Row.names)
x2$Row.names <- NULL
HCclus5000vsHist <- table(x2$HCCluster5000, x2$Histology)
chisq.test(HCclus5000vsHist)

#       b)
# Chi squared test of PAM 100 Genes(4 clusters) vs K-means
x3 <- merge(clusMem100, KMclusMem5000, by=0)
row.names(x3) <- c(x3$Row.names)
x3$Row.names <- NULL
clus100vsKMclus5000 <- table(x3$PAMCluster100, x3$KMCluster)
chisq.test(clus100vsKMclus5000)

# Chi squared test of PAM 100 Genes(4 clusters) vs Hclust
x4 <- merge(clusMem100, HCclusMem5000, by=0)
row.names(x4) <- c(x4$Row.names)
x4$Row.names <- NULL
clus100vsHCclus5000 <- table(x4$PAMCluster100, x4$HCCluster5000)
chisq.test(clus100vsHCclus5000)

# Chi squared test of K-means vs Hclust
x5 <- merge(KMclusMem5000, HCclusMem5000, by=0)
row.names(x5) <- c(x5$Row.names)
x5$Row.names <- NULL
KMclus5000vsHCclus5000 <- table(x5$KMCluster, x5$HCCluster5000)
chisq.test(KMclus5000vsHist)


# Adjust P-values

# Adjusted and Un-adjusted P-values Table
pTable <- data.frame (first_column  = c(chisq.test(clus100vsHist)$p.value, chisq.test(KMclus5000vsHist)$p.value, chisq.test(HCclus5000vsHist)$p.value, 
                                        chisq.test(clus100vsKMclus5000)$p.value, chisq.test(clus100vsHCclus5000)$p.value, chisq.test(KMclus5000vsHist)$p.value))


p.adjust(chisq.test(clus100vsKMclus5000)$p.value)


#         Question  3a
# Heatmap of different places in the 5000 differently expressed genes
# 2 heatmaps of two different places

# Hclust Heatmap

my_hclust_gene <- hclust(dist(mostVar5000)) #hclust clustering
my_gene_col <- cutree(tree = as.dendrogram(my_hclust_gene), k = 3) #Creating dendogram and cutting it for pheatmap
my_gene_col <- data.frame(cluster = ifelse(test = my_gene_col == 1, yes = "cluster 1", no = ifelse(test = my_gene_col ==2, yes = "cluster 2", no ="cluster 3")))#double Ifelse for 3 clusters

pheatmap(mostVar5000,  #heatmap for hclust with 3 rows/cols for dendogram
         annotation_row = my_gene_col,
         cutree_rows = 3,
         cutree_cols = 3,
         main="Hclust clustering for 5000 variable genes") #title

# K-Means Heatmap

pheatmap(t(mostVar5000),  #pheatfunction + transpose
         main="Kmeans clustering for 5000 variable genes", #title
         kmeans_k = 3,  #cluster count
         cutree_rows = 3, #dendogram rows
         cutree_cols = 3) #dendogram columns

# PAM heatmap
#so many packages to be able to create this one

install.packages("PAMhm")
install.packages("heatmapFlex")
install.packages("Heatplus")
install.packages("cluster")
install.packages("devtools")  
devtools::install_github("vfey/PAMhm")
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("Heatplus")
library(PAMhm)
library(heatmapFlex)
library(Biobase)
library(BiocGenerics)
library(parallel)
library(Heatplus)
library(cluster)
#Pam heatmap
mat <- matrix(t(mostVar5000),nrow=102) #102 sample columns
PAM.hm(mat, cluster.number = 3) # 3 clusters










