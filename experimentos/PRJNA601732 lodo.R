#PRJNA290995
'''
if (!requireNamespace("dplyr", quietly = TRUE)) {
  install.packages("dplyr")
}
if (!requireNamespace("tidyr", quietly = TRUE)) {
  install.packages("tidyr")
}
if (!requireNamespace("biomaRt", quietly = TRUE)) {
  BiocManager::install("biomaRt")
}'''
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("clusterProfiler")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("pathview")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("apeglm")
# Instalar o pacote ggrepel se necessário
if (!requireNamespace("ggrepel", quietly = TRUE)) {
  install.packages("ggrepel")
}
# Instalar e carregar o pacote org.Mm.eg.db se necessário
if (!requireNamespace("org.Mm.eg.db", quietly = TRUE)) {
  install.packages("BiocManager")
  BiocManager::install("org.Mm.eg.db")
}

# Carregar o pacote ggrepel


library(org.Mm.eg.db)
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
library(dplyr) 


tabular_dir  <- "/Users/carlitos/Desktop/experimentos/exps fasta/fastas/Leishmania/PRJNA601732 - lodo//"
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


### NOTE: If there are batch effects observed, correct for them before moving ahead


#exclude outlier samples
samples.to.be.excluded <- c("Lama_Infected_72h_R3_SRR2163299")
data.subset <- data[,!(colnames(data) %in% samples.to.be.excluded)]
colnames(data.subset)
# for ( i in colnames(data.subset)){
#   print(i)
# }

# 3. Normalization ----------------------------------------------------------------------
# create a deseq2 dataset

phenoData  <-  read_excel("/Users/carlitos/Desktop/aaaaaaaaaaaaaaaaaaaaaaaaaaa - lodo.xlsx", col_names = TRUE)
lista  <- phenoData$id
phenoData <- phenoData[,-1]
rownames(phenoData) <- lista

# Visualizar o resultado
print(phenoData)


# exclude outlier samples
colData <- phenoData %>% 
  filter(!row.names(.) %in% samples.to.be.excluded)


# fixing column names in colData
names(colData)
names(colData) <- gsub(':ch1', '', names(colData))
names(colData) <- gsub('\\s', '_', names(colData))
# making the rownames and column names identical
all(rownames(colData) %in% colnames(data.subset))
all(rownames(colData) == colnames(data.subset))

cont = 1
for (i in rownames(colData)){
  # print(rownames(colData)[cont])
  a = colnames(data.subset)[cont]
  b = rownames(colData)[cont]
  print(a == b)
  print(a)
  print(b)
  # print(cont)
  cont = cont + 1
}

(rownames(colData))[2] == (colnames(data.subset))[2]

# create dds
dds <- DESeqDataSetFromMatrix(countData = data.subset,
                              colData = colData,
                              design = ~ 1) # not spcifying model



## remove all genes with counts < 15 in more than 75% of samples (31*0.75=23.25)
## suggested by WGCNA on RNAseq FAQ

dds75 <- dds[rowSums(counts(dds) >= 15) >= 15,]
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


# convert matrix to numeric
norm.counts[] <- sapply(norm.counts, as.numeric)

soft_power <- 12
temp_cor <- cor
cor <- WGCNA::cor


# memory estimate w.r.t blocksize
bwnet <- blockwiseModules(norm.counts,
                          maxBlockSize = 14000,
                          TOMType = "signed",
                          power = soft_power,
                          mergeCutHeight = 0.25,
                          numericLabels = FALSE,
                          randomSeed = 1234,
                          verbose = 3)

#######################################
# Obter a data e hora atual
current_time <- Sys.time()

# Formatar a data e hora
formatted_time <- format(current_time, "%Y-%m-%d %H:%M:%S")

# Imprimir a data e hora formatada
print(formatted_time)
###########################################
cor <- temp_cor

# save(bwnet, file = "/Users/carlitos/Desktop/bwnet_PRJNA601732.RData")
# load("/Users/carlitos/Desktop/bwnet.RData")

# 5. Module Eigengenes ---------------------------------------------------------
module_eigengenes <- bwnet$MEs


# Print out a preview
head(module_eigengenes)


# get number of genes for each module
color_table=table(bwnet$colors)
color_table
# write.csv(as.data.frame(color_table), file = "/Users/carlitos/Desktop/color_table.csv", row.names = TRUE)


# Plot the dendrogram and the module colors before and after merging underneath
plotDendroAndColors(bwnet$dendrograms[[1]], cbind(bwnet$unmergedColors, bwnet$colors),
                    c("unmerged", "merged"),
                    dendroLabels = FALSE,
                    addGuide = TRUE,
                    hang= 0.03,
                    guideHang = 0.05)




# grey module = all genes that doesn't fall into other modules were assigned to the grey module





# 6A. Relate modules to traits --------------------------------------------------
# module trait associations


# create traits file - binarize categorical variables
traits <- colData %>% 
  mutate(disease_state_bin = ifelse(grepl('Ldonovani', Treatment), 1, 0)) #%>% select(2)


# print(unique(traits))
# binarize categorical variables
unique(colData$Treatment)

colData$severity <- factor(colData$Treatment, levels = c(
  "Infected_Mouse_Ldonovani","Uninfected_Mouse"))

severity.out <- binarizeCategoricalColumns(colData$severity,
                                           includePairwise = FALSE,
                                           includeLevelVsAll = TRUE,
                                           minCount = 1)

traits <- cbind(traits, severity.out)
rownames(traits) <- traits$SampleName
traits <- traits %>% select(-1,-2,-3)

# Define numbers of genes and samples
nSamples <- nrow(norm.counts)
nGenes <- ncol(norm.counts)



module.trait.corr <- cor(module_eigengenes, traits, use = 'p')
module.trait.corr.pvals <- corPvalueStudent(module.trait.corr, nSamples)



# visualize module-trait association as a heatmap
nrow(module_eigengenes)
nrow(traits)

heatmap.data <- merge(module_eigengenes, traits, by = 'row.names')

head(heatmap.data)

heatmap.data <- heatmap.data %>% 
  column_to_rownames(var = 'Row.names')

colnames(heatmap.data)
# ordered_columns <- c(
#   "MEpurple", "MEgreen", "MEbrown", "MEyellow", "MEturquoise", "MEmagenta",
#   "MEpink", "MEred", "MEblack", "MEblue", "MEgreenyellow", "MEgrey",
#   "disease_state_bin", "data.Lama_Uninfected_24h.vs.all", "data.Lama_Uninfected_48h.vs.all",
#   "data.Lama_Uninfected_72h.vs.all", "data.Lama_Infected_4h.vs.all", "data.Lama_Infected_24h.vs.all",
#   "data.Lama_Infected_48h.vs.all", "data.Lama_Infected_72h.vs.all"
# )
# 
# # Reordenar as colunas
# heatmap.data <- heatmap.data %>% select(all_of(ordered_columns))


heatmap.data
colnames(heatmap.data)

colnames(traits)
colnames(heatmap.data)
CorLevelPlot(heatmap.data,
             x = names(heatmap.data)[10:11],
             y = names(heatmap.data)[1:9],rotLabX = 50,
             col = c("cyan", "white", "grey", "purple"))

# write.csv(as.data.frame(heatmap.data), file = "/Users/carlitos/Desktop/PRJNA290995_lama_heatmap.data.csv", row.names = TRUE)

module.gene.mapping <- as.data.frame(bwnet$colors)
module.gene.mapping %>% 
  filter(`bwnet$colors` == 'yellow') %>% 
  rownames()

class(bwnet$colors)

valores_interesse <- c('228994','213541','229096','231386','26383','56335','210529','268420','240255','69113','211064','66400','52463',
                       '66926','101867','52575','68789','328162','28114','13434','98956','21681','22608','107435','14534','18519','81601','244349',
                       '54169','217127','67773','269252','12914','73242','328572','433759','15182','15183','208727','15184',
                       '15185','56233','70315','79221','170787','93759','64383','64384','75387','68346','50721','209011')

bwnet$colors[valores_interesse]
table(bwnet$colors[valores_interesse])



# 6B. Intramodular analysis: Identifying driver genes ---------------

# Calculate the module membership and the associated p-values

# The module membership/intramodular connectivity is calculated as the correlation of the eigengene and the gene expression profile. 
# This quantifies the similarity of all genes on the array to every module.

module.membership.measure <- cor(module_eigengenes, norm.counts, use = 'p')
module.membership.measure.pvals <- corPvalueStudent(module.membership.measure, nSamples)


module.membership.measure.pvals[1:12,1:12]


# Calculate the gene significance and associated p-values

gene.signf.corr <- cor(norm.counts, severity.out$data.Lama_Infected_4h.vs.all, use = 'p')
gene.signf.corr.pvals <- corPvalueStudent(gene.signf.corr, nSamples)
gene.signf.corr=0
gene.signf.corr.pvals=0
gene.signf.corr.pvals %>% 
  as.data.frame() %>% 
  arrange(V1) %>% 
  head(25)


# Using the gene significance you can identify genes that have a high significance for trait of interest 
# Using the module membership measures you can identify genes with high module membership in interesting modules.

# 7. analisis enrich
green_genes <- module.gene.mapping %>% 
  filter(`bwnet$colors` == 'yellow') %>% 
  rownames()



# library(DESeq2)
# 
# colData <- colData[,-3]
# 
# dds_ana <- DESeqDataSetFromMatrix(countData = data.subset,
#                               colData = colData,
#                               design = ~Treatment)
# 
# 
# dds_ana <- DESeq(dds_ana)
# 
# #######
# 
# # Comparação personalizada entre Lama_Infected_4h e Lama_Uninfected_4h
# res <- results(dds_ana, contrast = c("Treatment", "Lama_Infected_72h", "Lama_Uninfected_72h"))
# 
# # Aplicar shrinkage do log2 fold change
# res <- lfcShrink(dds_ana, contrast = c("Treatment", "Lama_Infected_72h", "Lama_Uninfected_72h"),type = "normal")
# 
# # Visualizar resultados
# head(res)
# 
# # Filtrar genes com p-valor ajustado < 0.05 e |log2FoldChange| > 1
# diff_genes <- rownames(res[which(res$padj < 0.05 & abs(res$log2FoldChange) > 1), ])
# diff_genes
# library(clusterProfiler)
# library(org.Hs.eg.db)  # Use a base de dados apropriada para o seu organismo
# # Supondo que os IDs dos genes sejam Entrez IDs
# ego <- enrichGO(gene          = diff_genes,
#                 OrgDb         = org.Hs.eg.db,
#                 keyType       = 'ENTREZID',  # Ajuste o keyType conforme necessário
#                 ont           = "ALL",  # Pode ser "BP", "MF", "CC", ou "ALL"
#                 pAdjustMethod = "BH",
#                 qvalueCutoff  = 0.05,
#                 readable      = TRUE)
# 
# # Visualizar os resultados
# head(ego)
# library(enrichplot)
# 
# # Dotplot
# dotplot(ego)
# 
# # Barplot
# barplot(ego, showCategory = 10)
# 
# # EnrichMap
# emapplot(ego)

library(clusterProfiler)
library(org.Hs.eg.db)  # Use a base de dados apropriada para o seu organismo





# Supondo que os IDs dos genes sejam Entrez IDs
library(org.Mm.eg.db)

# Realizar a análise de enriquecimento GO
ego <- enrichGO(
  gene          = green_genes,
  OrgDb         = org.Mm.eg.db,
  keyType       = 'ENTREZID',  # Ajuste o keyType conforme necessário
  ont           = "ALL",  # Pode ser "BP", "MF", "CC", ou "ALL"
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

# Visualizar os resultados
head(ego)

# Visualize os resultados
head(ego)

library(enrichplot)

# Dotplot
dotplot(ego)

 # Barplot
barplot(ego, showCategory = 10)

# Supondo que o objeto ego já está criado
ego <- pairwise_termsim(ego)

library(ggplot2)

emap <- emapplot(ego)
# Ajustar o tamanho dos labels dos termos para um tamanho menor
emap + theme(
  legend.text = element_text(size = 8),
  legend.title = element_text(size = 8),
  plot.title = element_text(size = 10, face = "bold"),
  axis.text.x = element_text(size = 8),
  axis.text.y = element_text(size = 8)
)

ego_df <- as.data.frame(ego)
write.csv(ego_df, file = "/Users/carlitos/Desktop/enrichment_results_lodo_PRJNA601732.csv", row.names = FALSE)
