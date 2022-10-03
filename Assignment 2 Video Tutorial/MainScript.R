# Script following video tutorial

library(dplyr)
library(tidyverse)
library(GEOquery)

# Imports the data
dat <- read.csv(file = "data/GSE183947_fpkm.csv")

# Prints dimension, rows x columns
dim(dat)

# Get metadata
gse <- getGEO(GEO = "GSE183947", GSEMatrix = TRUE)
# (skipped something here, not sure if we need it?)

metadata <- pData(phenoData(gse[[1]]))

# Select specific columns from metadata
metadata %>%
  select(1,10,11,17) %>% # 1, 10, 11, 17 are hand picked
  head()
# (skipped renaming the columns)
# (skipped reformatting)

# Reshaping data
dat.long <- dat %>%
  rename(gene = X) %>%
  gather(key = 'samples', value = 'FPKM')

# Join dataframes (dat.long + metadata)
dat.long <- dat.long %>%
  left_join(., metadata, by = c("samples" = "description")) %>%

# Explore data


# ----- PLOTS -----

# Density plot
dat.long %<%
  filter(gene == 'BRCA1') %>%
  ggplot(., aes(x = FPKM, fill = tissue)) +
  geom_density()
  
  
  
  
