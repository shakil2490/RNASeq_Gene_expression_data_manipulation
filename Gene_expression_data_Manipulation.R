# load libraries
library(dplyr)
library(tidyverse)
library(GEOquery)

# read in the data ---------
data <- read.csv(file = "data/GSE183947_fpkm.csv")
dim(data)

# get metadata --------
gse <- getGEO(GEO = 'GSE183947', GSEMatrix = TRUE)

gse

metadata <- pData(phenoData(gse[[1]]))
head(metadata)

# select, mutate, rename ------------
metadata.modified <- metadata %>%
  select(1,10,11,17) %>%
  rename(tissue = characteristics_ch1) %>%
  rename(metastasis = characteristics_ch1.1) %>%
  mutate(tissue = gsub("tissue: ", "", tissue)) %>%
  mutate(metastasis = gsub("metastasis: ", "", metastasis))

head(data)

# reshaping data - from wide to long--------
data.long <- data %>%
  rename(gene = X) %>%
  gather(key = 'samples', value = 'FPKM', -gene)


# join dataframes = dat.long + metadata.modified
data.long <- data.long %>%
  left_join(., metadata.modified, by = c("samples" = "description")) 

# filter, group_by, summarize and arrange 
data.long %>%
  filter(gene == 'BRCA1' | gene == 'BRCA2') %>%
  group_by(gene, tissue.x) %>%
  summarize(mean_FPKM = mean(FPKM),
            median_FPKM = median(FPKM)) %>%
  arrange(mean_FPKM)
# barplot
data.long %>%
  filter(gene == 'BRCA1') %>%
  ggplot(., aes(x = samples, y = FPKM, fill = tissue.x)) +
  geom_col()

# density
data.long %>%
  filter(gene == 'BRCA1') %>%
  ggplot(., aes(x = FPKM, fill = tissue.x)) +
  geom_density(alpha = 0.3)

# boxplot 
data.long %>%
  filter(gene == 'BRCA1') %>%
  ggplot(., aes(x = metastasis.x, y = FPKM)) +
  geom_violin()


data.long %>%
  filter(gene == 'BRCA1') %>%
  ggplot(., aes(x = metastasis.x, y = FPKM)) +
  geom_boxplot()
#scatterplot
data.long %>%
  filter(gene == 'BRCA1' | gene == 'BRCA2') %>%
  spread(key = gene, value = FPKM) %>%
  ggplot(., aes(x = BRCA1, y = BRCA2, color = tissue.x)) +
  geom_point() +
  geom_smooth(method = 'lm', se = FALSE)

#heatmap
genes.of.interest <- c('BRCA1', 'BRCA2', 'TP53', 'ALK', 'MYCN')
data.long %>%
  filter(gene %in% genes.of.interest) %>%
  ggplot(., aes(x = samples, y = gene, fill = FPKM)) +
  geom_tile() +
  scale_fill_gradient(low = 'white', high = 'red')











