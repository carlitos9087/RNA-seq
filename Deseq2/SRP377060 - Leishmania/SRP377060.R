setwd("/Users/carlitos/Desktop/RNA-seq/")
# setwd("/Users/carlitos/Documents/")
library(dplyr)
library(ggrepel)
library(biomaRt)
library(dplyr)
library(tidyr)
library(dplyr)
library(biomaRt)
library(dplyr)
library(WGCNA)
library(DESeq2)
library(GEOquery)
library(tidyverse)
library(CorLevelPlot)
library(gridExtra)
library(stringr)
library(readxl)
library(openxlsx)
library(org.Hs.eg.db) 


tabular_dir  <- "Deseq2/SRP377060 - Leishmania/"
tabular_files <- list.files(path = tabular_dir, pattern = "\\.tabular$", full.names = TRUE)


read_tabular_file <- function(file) {
  df <- read.table(file, header = TRUE, stringsAsFactors = FALSE, sep="\t")
  colnames(df)[2] <- tools::file_path_sans_ext(basename(file))  # Renomeia a coluna "counts" para o nome do arquivo
  return(df)
}


tabular_dfs <- lapply(tabular_files, read_tabular_file)
combined_df <- Reduce(function(x, y) full_join(x, y, by = "Geneid"), tabular_dfs)
head(combined_df)
rownames(combined_df) = combined_df$Geneid
combined_df = combined_df[,-1]

data = combined_df
colnames(data)
# dados da amostra
phenoData  <-  read_excel("Deseq2/SRP377060 - Leishmania/SRP377060.xlsx", col_names = TRUE)
lista  <- phenoData$id
phenoData <- phenoData[,-1]
rownames(phenoData) <- lista


head(data)
phenoData

############################################################################################

# 2. QC - outlier detection ------------------------------------------------
# detect outlier genes

gsg <- goodSamplesGenes(t(data))
summary(gsg)
gsg$allOK

table(gsg$goodGenes)
table(gsg$goodSamples)

# remove genes that are detectd as outliers
data <- data[gsg$goodGenes == TRUE,]

# detect outlier samples - hierarchical clustering - method 1
htree <- hclust(dist(t(data)), method = "average")
plot(htree)


###########################################################
# pca - method 2
pca <- prcomp(t(data))
pca.dat <- pca$x

pca.var <- pca$sdev^2
pca.var.percent <- round(pca.var/sum(pca.var)*100, digits = 2)

pca.dat <- as.data.frame(pca.dat)

ggplot(pca.dat, aes(PC1, PC2)) +
  geom_point() +
  geom_text(label = rownames(pca.dat)) +
  labs(x = paste0('PC1: ', pca.var.percent[1], ' %'),
       y = paste0('PC2: ', pca.var.percent[2], ' %'))
###########################################################
colnames(data)
# exclude outlier samples
samples.to.be.excluded <- c("together_than_RNA_was_isolated_SRR1184508")
data.subset <- data[,!(colnames(data) %in% samples.to.be.excluded)]



pca <- prcomp(t(data.subset))
pca.dat <- pca$x
pca.var <- pca$sdev^2
pca.var.percent <- round(pca.var/sum(pca.var)*100, digits = 2)
pca.dat <- as.data.frame(pca.dat)
library(ggrepel)
ggplot(pca.dat, aes(PC1, PC2)) +
  geom_point() +
  geom_text_repel(aes(label = rownames(pca.dat)),
                  box.padding = 0.5,   # Ajuste o espaçamento ao redor do rótulo
                  point.padding = 0.5, # Ajuste o espaçamento entre o ponto e o rótulo
                  segment.color = 'grey50', # Cor das linhas conectando rótulos aos pontos
                  segment.size = 0.5) +    # Espessura das linhas
  labs(x = paste0('PC1: ', pca.var.percent[1], ' %'),
       y = paste0('PC2: ', pca.var.percent[2], ' %'))


#########################################################################################
log_counts <- log2(data.subset + 1)
# add a colorbar along the heatmap with sample condition
x = melt(as.matrix(log_counts))

colnames(x) = c('gene_id', 'sample', 'value')
ggplot(x, aes(x=value, color=sample)) + 
  geom_density() +
  guides(color = guide_legend(ncol = 2)) +
  theme( legend. = 2,
         legend.margin = margin(0, 0, 0, 0),  # Ajuste conforme necessário
         legend.position = "right",
         legend.text = element_text(size = 8),   # Ajuste o valor para o tamanho da fonte
         legend.key.size = unit(1.3, "lines")    # Mantém o tamanho dos ícones; ajuste conforme necessário
  )




#########################################################################################

# 3. Normalization ----------------------------------------------------------------------
# create a deseq2 dataset

row.names(phenoData)
colnames(data)
# exclude outlier samples
samples.to.be.excluded <- c("together_than_RNA_was_isolated_SRR1184508")
data.subset <- data[,!(colnames(data) %in% samples.to.be.excluded)]

# exclude outlier samples
colData <- phenoData %>% 
  filter(!row.names(.) %in% samples.to.be.excluded)


# fixing column names in colData
names(colData)
names(colData) <- gsub(':ch1', '', names(colData))
names(colData) <- gsub('\\s', '_', names(colData))
colData <- colData[1:2]
colnames(data.subset)
rownames(colData)<- colData$SampleName
rownames(colData)
# making the rownames and column names identical
all(rownames(colData) %in% colnames(data.subset))
all(rownames(colData) == colnames(data.subset))

# cont = 1
# for (i in rownames(colData)){
#   # print(rownames(colData)[cont])
#   a = colnames(data.subset)[cont]
#   b = rownames(colData)[cont]
#   print(a == b)
#   print(a)
#   print(b)
#   # print(cont)
#   cont = cont + 1
# }

(rownames(colData))[2] == (colnames(data.subset))[2]

# create dds
dds <- DESeqDataSetFromMatrix(countData = data.subset,
                              colData = colData,
                              design = ~ 1) # not spcifying model



## remove all genes with counts < 15 in more than 75% of samples (31*0.75=23.25)
## suggested by WGCNA on RNAseq FAQ

dds75 <- dds[rowSums(counts(dds) >= 15) >= 14,]
nrow(dds75) # 13284 genes

# perform variance stabilization
dds_norm <- vst(dds75)


# get normalized counts
norm.counts <- assay(dds_norm) %>% 
  t()



###########################################
# 4. Network Construction  ---------------------------------------------------
# Choose a set of soft-thresholding powers
power <- c(c(1:10), seq(from = 12, to = 50, by = 2))

# Call the network topology analysis function
sft <- pickSoftThreshold(norm.counts,
                         powerVector = power,
                         networkType = "signed",
                         verbose = 5)


sft.data <- sft$fitIndices

# visualization to pick power

a1 <- ggplot(sft.data, aes(Power, SFT.R.sq, label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.8, color = 'red') +
  labs(x = 'Power', y = 'Scale free topology model fit, signed R^2') +
  theme_classic()


a2 <- ggplot(sft.data, aes(Power, mean.k., label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  labs(x = 'Power', y = 'Mean Connectivity') +
  theme_classic()


grid.arrange(a1, a2, nrow = 2)









