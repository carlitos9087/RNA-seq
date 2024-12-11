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
'''if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("Mus.musculus", "Homo.sapiens"))
'''


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



tabular_dir  <- "/Users/carlitos/Desktop/RNA-seq/experimentos/fastas/Aspergilus/GSE55648(homo sapiens)"
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
colnames(data) <- c("separably_were_united_for_RNA_isolation_SRR1184521","separably_were_united_for_RNA_isolation_SRR1184522","separably_were_united_for_RNA_isolation_SRR1184523","separably_were_united_for_RNA_isolation_SRR1184524",
                    "separably_were_united_for_RNA_isolation_SRR1184525","separably_were_united_for_RNA_isolation_SRR1184526","separably_were_united_for_RNA_isolation_SRR1184527","separably_were_united_for_RNA_isolation_SRR1184528",
                    "separably_were_united_for_RNA_isolation_SRR1184529","separably_were_united_for_RNA_isolation_SRR1184530","separably_were_united_for_RNA_isolation_SRR1184531","separably_were_united_for_RNA_isolation_SRR1184532",
                    "together_than_RNA_was_isolated_SRR1184505","together_than_RNA_was_isolated_SRR1184506","together_than_RNA_was_isolated_SRR1184507","together_than_RNA_was_isolated_SRR1184508",
                    "together_than_RNA_was_isolated_SRR1184509","together_than_RNA_was_isolated_SRR1184510","together_than_RNA_was_isolated_SRR1184511","together_than_RNA_was_isolated_SRR1184512",
                    "together_than_RNA_was_isolated_SRR1184513","together_than_RNA_was_isolated_SRR1184514","together_than_RNA_was_isolated_SRR1184515","together_than_RNA_was_isolated_SRR1184516",
                    "together_than_RNA_was_isolated_SRR1184517","together_than_RNA_was_isolated_SRR1184518","together_than_RNA_was_isolated_SRR1184519","together_than_RNA_was_isolated_SRR1184520")

colnames(data)
# dados da amostra
phenoData  <-  read_excel("/Users/carlitos/Desktop/RNA-seq/experimentos/GSE55648/GSE55648.xlsx", col_names = TRUE)
lista  <- phenoData$id
phenoData <- phenoData[,-1]
rownames(phenoData) <- lista





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

dds75 <- dds[rowSums(counts(dds) >= 15) >= 21,]
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

soft_power <- 14
temp_cor <- cor
cor <- WGCNA::cor
print("das")

# memory estimate w.r.t blocksize
bwnet <- blockwiseModules(norm.counts,
                          maxBlockSize = 14000,
                          TOMType = "signed",
                          power = soft_power,
                          mergeCutHeight = 0.25,
                          numericLabels = FALSE,
                          randomSeed = 1234,
                          verbose = 3)



cor <- temp_cor

save(bwnet, file = "/Users/carlitos/Desktop/bwnet_GSE55648.RData")
load("/Users/carlitos/Desktop/resultados/PRJNA290995/bwnet.RData")

# 5. Module Eigengenes ---------------------------------------------------------
module_eigengenes <- bwnet$MEs


# Print out a preview
head(module_eigengenes)


# get number of genes for each module
color_table=table(bwnet$colors)
color_table
color_table
barplot(color_table, 
        main = "Distribuição de Cores", 
        xlab = "Cores", 
        ylab = "Frequência", 
        col = names(color_table), 
        las = 2, 
        cex.names = 0.8)
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
  mutate(disease_state_bin = ifelse(grepl('_Infected_', Treatment), 1, 0)) #%>% select(2)


# print(unique(traits))
# binarize categorical variables
unique(colData$Treatment)

colData$severity <- factor(colData$Treatment, levels = c("separably" ,"together"))

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
module.trait.corr.pvals


# visualize module-trait association as a heatmap
nrow(module_eigengenes)
nrow(traits)

heatmap.data <- merge(module_eigengenes, traits, by = 'row.names')

head(heatmap.data)

heatmap.data <- heatmap.data %>% 
  column_to_rownames(var = 'Row.names')


# ordered_columns <- c(
#   "MEpurple", "MEgreen", "MEbrown", "MEyellow", "MEturquoise", "MEmagenta",
#   "MEpink", "MEred", "MEblack", "MEblue", "MEgreenyellow", "MEgrey",
#   "disease_state_bin", "data.Lama_Uninfected_24h.vs.all", "data.Lama_Uninfected_48h.vs.all",
#   "data.Lama_Uninfected_72h.vs.all", "data.Lama_Infected_4h.vs.all", "data.Lama_Infected_24h.vs.all",
#   "data.Lama_Infected_48h.vs.all", "data.Lama_Infected_72h.vs.all", "data.Lmj_Uninfected_4h.vs.all",
#   "data.Lmj_Uninfected_24h.vs.all", "data.Lmj_Uninfected_48h.vs.all", "data.Lmj_Uninfected_72h.vs.all",
#   "data.Lmj_Infected_4h.vs.all", "data.Lmj_Infected_24h.vs.all", "data.Lmj_Infected_48h.vs.all",
#   "data.Lmj_Infected_72h.vs.all"
# )

# Reordenar as colunas
#heatmap.data <- heatmap.data %>% select(all_of(ordered_columns))


heatmap.data
colnames(heatmap.data)

colnames(traits)
colnames(heatmap.data)
CorLevelPlot(heatmap.data,
             x = names(heatmap.data)[16:16],
             y = names(heatmap.data)[1:15],rotLabX = 50,
             col = c("blue1", "skyblue", "white", "pink", "red"))

# write.csv(as.data.frame(heatmap.data), file = "/Users/carlitos/Desktop/heatmap.data.csv", row.names = TRUE)

module.gene.mapping <- as.data.frame(bwnet$colors)
head(module.gene.mapping)
module.gene.mapping %>% 
  filter(`bwnet$colors` == 'yellow') %>% 
  rownames()



valores_interesse <- c("54915", "51441", "253943", "91746", "79068", "56339",
                       "57721", "54890", "64848", "221120", "8846", "84266",
                       "80312", "51605", "23378", "55006", "115708", "54888",
                       "1787", "55226", "10189", "4904", "8520", "2648",
                       "8850", "10524", "7994", "23522", "11143", "84148",
                       "9329", "1387", "79969", "2033", "3065", "3066",
                       "8841", "9759", "10014", "10013", "51564", "55869",
                       "9734", "83933", "23411", "22933", "23410", "23409",
                       "23408", "51548", "51547")

table(bwnet$colors[valores_interesse])

cores_interesse <- c("greenyellow", "green")
ids_filtrados <- valores_interesse[bwnet$colors[valores_interesse] %in% cores_interesse]
ids_filtrados
# 6B. Intramodular analysis: Identifying driver genes ---------------

# Calculate the module membership and the associated p-values

# The module membership/intramodular connectivity is calculated as the correlation of the eigengene and the gene expression profile. 
# This quantifies the similarity of all genes on the array to every module.

module.membership.measure <- cor(module_eigengenes, norm.counts, use = 'p')
module.membership.measure.pvals <- corPvalueStudent(module.membership.measure, nSamples)


module.membership.measure.pvals[1:12,1:12]


# Calculate the gene significance and associated p-values

gene.signf.corr <- cor(norm.counts, severity.out$data.together.vs.all, use = 'p')
gene.signf.corr.pvals <- corPvalueStudent(gene.signf.corr, nSamples)


gene.signf.corr.pvals %>% 
  as.data.frame() %>% 
  arrange(V1) %>% 
  head(25)
head(rownames(gene.signf.corr.pvals))
write.csv(as.data.frame(gene.signf.corr.pvals), file = "/Users/carlitos/Desktop/gene.signf.corr.pvals.csv", row.names = TRUE)
# Using the gene significance you can identify genes that have a high significance for trait of interest 
# Using the module membership measures you can identify genes with high module membership in interesting modules.

##################################################################################################################
library(Homo.sapiens)
keytypes(Homo.sapiens)
# Example query
gene_ids <- head(keys(Homo.sapiens, keytype='ENTREZID'), 2)

select(Homo.sapiens, keytype='ENTREZID', keys=gene_ids, 
       columns=c('ALIAS','SYMBOL', 'TXCHROM', 'TXSTART', 'TXEND'))
gene_ids
##################################################################################################################
head(module.gene.mapping)
module.gene.mapping %>% 
  filter(`bwnet$colors` == 'black') %>% 
  rownames()

# pick out a few modules of interest here




modules_of_interest = c("black")
modules_of_interest = c("black", "pink", "yellow")
modules_of_interest = c("black", "turquoise", "pink", "yellow")

module_df <- data.frame(
  gene_id = names(bwnet$colors),
  colors = labels2colors(bwnet$colors)
)

# Pull out list of genes in that module
submod = module_df %>%
  subset(colors %in% modules_of_interest)

row.names(module_df) = module_df$gene_id

# Get normalized expression for those genes

expr_normalized <- t(norm.counts)

colnames(expr_normalized)
lista  <- phenoData$Treatment
lista <- gsub("_4h", "_04h", lista)
lista

colnames(expr_normalized) <- lista
expr_normalized
subexpr = expr_normalized[submod$gene_id,]

submod_df = data.frame(subexpr) %>%
  mutate(
    gene_id = row.names(.)
  ) %>%
  pivot_longer(-gene_id) %>%
  mutate(
    module = module_df[gene_id,]$colors
  )

submod_df %>% ggplot(., aes(x=name, y=value, group=gene_id)) +
  geom_line(aes(color = module),
            alpha = 0.6) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90)
  ) +
  facet_grid(rows = vars(module)) +
  labs(x = "treatment",
       y = "normalized expression")
colnames(expr_normalized) <- rownames(norm.counts)

#################################################################################
#Enriquecimento
#################################################################################
# Carregar pacotes necessários
library(ggplot2)
library(enrichplot)
library(clusterProfiler)
library(org.Hs.eg.db)  # Use o banco de dados apropriado para o seu organismo
# Extrair genes do módulo preto
black_genes <- names(bwnet$colors[bwnet$colors == "black"])

ego <- enrichGO(gene          = black_genes,
                OrgDb         = org.Hs.eg.db,   # Substitua pelo banco de dados apropriado
                keyType       = "ENTREZID",       # Substitua por "ENTREZID" se usar IDs Entrez
                ont           = "ALL",
                library(beepr)          pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.05)

dotplot(ego, showCategory = 10)  # Dotplot dos termos enriquecidos

barplot(ego, showCategory = 10)

ego <- pairwise_termsim(ego)

emap <- emapplot(ego)
emap + theme(
  legend.text = element_text(size = 8),
  legend.title = element_text(size = 8),
  plot.title = element_text(size = 10, face = "bold"),
  axis.text.x = element_text(size = 8),
  axis.text.y = element_text(size = 8)
)

#####
ekegg <- enrichKEGG(gene         = black_genes,
                    organism     = 'hsa',       # Substitua pelo organismo apropriado
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.05,
                    qvalueCutoff = 0.05)

# Visualizar resultados
dotplot(ekegg)
#####################################################################################'
pink_genes <- names(bwnet$colors[bwnet$colors == "pink"])

ego <- enrichGO(gene          = pink_genes,
                OrgDb         = org.Hs.eg.db,   # Substitua pelo banco de dados apropriado
                keyType       = "ENTREZID",       # Substitua por "ENTREZID" se usar IDs Entrez
                ont           = "ALL",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.05)

dotplot(ego, showCategory = 10)  # Dotplot dos termos enriquecidos

barplot(ego, showCategory = 10)

ego <- pairwise_termsim(ego)

emap <- emapplot(ego)
emap + theme(
  legend.text = element_text(size = 8),
  legend.title = element_text(size = 8),
  plot.title = element_text(size = 7, face = "bold"),
  axis.text.x = element_text(size = 8),
  axis.text.y = element_text(size = 8)
)

#####
ekegg <- enrichKEGG(gene         = pink_genes,
                    organism     = 'hsa',       # Substitua pelo organismo apropriado
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.05,
                    qvalueCutoff = 0.05)

# Visualizar resultados
dotplot(ekegg, showCategory = 10)

barplot(ekegg, showCategory = 10)

ego <- pairwise_termsim(ekegg)

emap <- emapplot(ekegg)
emap + theme(
  legend.text = element_text(size = 8),
  legend.title = element_text(size = 8),
  plot.title = element_text(size = 7, face = "bold"),
  axis.text.x = element_text(size = 8),
  axis.text.y = element_text(size = 8)
)
#####################################################################################
modules_of_interest = c("black", "turquoise", "pink", "yellow")
turquoise_genes <- names(bwnet$colors[bwnet$colors == "turquoise"])

ego <- enrichGO(gene          = turquoise_genes,
                OrgDb         = org.Hs.eg.db,   # Substitua pelo banco de dados apropriado
                keyType       = "ENTREZID",       # Substitua por "ENTREZID" se usar IDs Entrez
                ont           = "ALL",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.05)

dotplot(ego, showCategory = 10)  # Dotplot dos termos enriquecidos

barplot(ego, showCategory = 10)

ego <- pairwise_termsim(ego)

emap <- emapplot(ego)
emap + theme(
  legend.text = element_text(size = 8),
  legend.title = element_text(size = 8),
  plot.title = element_text(size = 10, face = "bold"),
  axis.text.x = element_text(size = 8),
  axis.text.y = element_text(size = 8)
)

#####
ekegg <- enrichKEGG(gene         = turquoise_genes,
                    organism     = 'hsa',       # Substitua pelo organismo apropriado
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.05,
                    qvalueCutoff = 0.05)

# Visualizar resultados
dotplot(ekegg, showCategory = 10)

barplot(ekegg, showCategory = 10)

ego <- pairwise_termsim(ekegg)

emap <- emapplot(ekegg)
emap + theme(
  legend.text = element_text(size = 8),
  legend.title = element_text(size = 8),
  plot.title = element_text(size = 7, face = "bold"),
  axis.text.x = element_text(size = 8),
  axis.text.y = element_text(size = 8)
)

#####################################################################################
modules_of_interest = c("black", "turquoise", "pink", "yellow")
yellow_genes <- names(bwnet$colors[bwnet$colors == "yellow"])

ego <- enrichGO(gene          = yellow_genes,
                OrgDb         = org.Hs.eg.db,   # Substitua pelo banco de dados apropriado
                keyType       = "ENTREZID",       # Substitua por "ENTREZID" se usar IDs Entrez
                ont           = "ALL",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.05)

dotplot(ego, showCategory = 10)  # Dotplot dos termos enriquecidos

barplot(ego, showCategory = 10)

ego <- pairwise_termsim(ego)

emap <- emapplot(ego)
emap + theme(
  legend.text = element_text(size = 8),
  legend.title = element_text(size = 8),
  plot.title = element_text(size = 10, face = "bold"),
  axis.text.x = element_text(size = 8),
  axis.text.y = element_text(size = 8)
)

#####
ekegg <- enrichKEGG(gene         = yellow_genes,
                    organism     = 'hsa',       # Substitua pelo organismo apropriado
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.05,
                    qvalueCutoff = 0.05)

# Visualizar resultados
dotplot(ekegg, showCategory = 10)

barplot(ekegg, showCategory = 10)

ego <- pairwise_termsim(ekegg)

emap <- emapplot(ekegg)
emap + theme(
  legend.text = element_text(size = 8),
  legend.title = element_text(size = 8),
  plot.title = element_text(size = 7, face = "bold"),
  axis.text.x = element_text(size = 8),
  axis.text.y = element_text(size = 8)
)
#####################################################################################

modules_of_interest = c("black", "turquoise", "pink", "yellow")

module_df <- data.frame(
  gene_id = names(bwnet$colors),
  colors = labels2colors(bwnet$colors)
)

module_df


# Extrair a lista de genes pertencentes aos módulos de interesse
submod <- module_df %>%
  subset(colors %in% modules_of_interest)

# Definir os nomes das linhas como os IDs dos genes
row.names(module_df) <- module_df$gene_id

# Obter a expressão normalizada para esses genes
expr_normalized <- t(norm.counts)

# Ajustar os nomes das colunas para incluir zeros à esquerda em "4h"
lista <- gsub("_4h", "_04h", phenoData$Treatment)
colnames(expr_normalized) <- lista

# Selecionar as expressões dos genes dos módulos de interesse
subexpr <- expr_normalized[submod$gene_id,]

# Transformar a expressão em um formato longo para visualização
submod_df <- data.frame(subexpr) %>%
  mutate(gene_id = row.names(.)) %>%
  pivot_longer(-gene_id) %>%
  mutate(module = module_df[gene_id,]$colors)

# Visualizar a expressão normalizada por tratamento
submod_df %>% 
  ggplot(aes(x = name, y = value, group = gene_id, color = module)) +
  geom_line(alpha = 0.6) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90)) +
  facet_grid(rows = vars(module)) +
  labs(x = "Treatment", y = "Normalized Expression")

# Selecionar genes de interesse
genes_of_interest <- module_df %>%
  subset(colors %in% modules_of_interest)
#############################################
#selected_genes <- module.gene.mapping %>%
#  filter(`bwnet$colors` %in% c("black", "turquoise", "pink", "yellow"))
#selected_genes
#class(bwnet$unmergedColors)
###########################################

# Obter a expressão normalizada desses genes
expr_of_interest <- expr_normalized[genes_of_interest$gene_id,]



cordist <- function(dat) {
  cor_matrix  <- cor(t(dat))
  
  dist_matrix <- as.matrix(dist(dat, diag=TRUE, upper=TRUE))
  dist_matrix <- log1p(dist_matrix)
  dist_matrix <- 1 - (dist_matrix / max(dist_matrix))
  
  sign(cor_matrix) * ((abs(cor_matrix) + dist_matrix)/ 2)
}


sim_matrix <- cordist(expr_of_interest)

heatmap_indices <- sample(nrow(sim_matrix))

heatmap.2(t(sim_matrix[heatmap_indices, heatmap_indices]),
          col=redgreen(75),
          labRow=NA, labCol=NA, 
          trace='none', dendrogram='row',
          xlab='Gene', ylab='Gene',
          main='Similarity matrix',
          density.info='none', revC=TRUE)

adj_matrix <- adjacency.fromSimilarity(sim_matrix, power=20, type='signed')
rm(sim_matrix)
gc()

# Convert to matrix
gene_ids <- rownames(adj_matrix)

adj_matrix <- matrix(adj_matrix, nrow=nrow(adj_matrix))
rownames(adj_matrix) <- gene_ids
colnames(adj_matrix) <- gene_ids

heatmap.2(t(adj_matrix[heatmap_indices, heatmap_indices]),
          col=redgreen(75),
          labRow=NA, labCol=NA, 
          trace='none', dendrogram='row',
          xlab='Gene', ylab='Gene',
          main='Adjacency matrix',
          density.info='none', revC=TRUE)

# Cluster gene expression profiles; the flashClust function from
# the authors of WGCNA is another options for larger datasets.
# For input, we use the reciprocal of the adjacency matrix; hierarchical
# clustering works by comparing the _distance_ between objects instead of the
# _similarity_.
gene_tree <- hclust(as.dist(1 - adj_matrix), method="average")

# we will use the cuttreeDynamicTree method to break apart the hc dendrogram
# into separate modules
module_labels <- cutreeDynamicTree(dendro=gene_tree, minModuleSize=15,
                                   deepSplit=TRUE)

# assign a color to each module for easier visualization and referencing
module_colors <- labels2colors(module_labels)
module_colors


#' An igraph graph object representing the exported graph.
export_network_to_graphml <- function (adj_mat, filename=NULL, weighted=TRUE,
                                       threshold=0.5, max_edge_ratio=3,
                                       nodeAttr=NULL, nodeAttrDataFrame=NULL,
                                       edgeAttributes=NULL, verbose=FALSE) {
  library('igraph')
  
  # Determine filename to use
  if (is.null(filename)) {
    filename <- 'network.graphml'
  }
  
  max_edges <- max_edge_ratio * nrow(adj_mat)
  
  edge_to_total_ratio <- max_edges / length(adj_mat)
  edge_limit_cutoff <- as.numeric(quantile(abs(adj_mat), 1 - edge_to_total_ratio))
  
  # Ensure at least some edges are left
  min_threshold <- as.numeric(quantile(abs(adj_mat), 0.9999))
  threshold <- min(min_threshold, max(threshold, edge_limit_cutoff))
  
  # Remove edges below the threshold
  adj_mat[abs(adj_mat) < threshold] <- 0
  
  # Drop orphaned nodes (nodes with no connections)
  orphaned <- (colSums(adj_mat) == 0)
  adj_mat <- adj_mat[!orphaned, !orphaned]
  
  # Remove orphaned annotations
  if (!is.null(nodeAttr)) {
    nodeAttr <- nodeAttr[!orphaned]
  }
  if (!is.null(nodeAttrDataFrame)) {
    nodeAttrDataFrame <- nodeAttrDataFrame[!orphaned, ]
  }
  
  # Normalize edges to range [0, 1]
  is_zero <- adj_mat == 0
  is_negative <- adj_mat < 0
  adj_mat <- (abs(adj_mat) - threshold) / (max(adj_mat) - threshold)
  adj_mat[is_zero] <- 0
  adj_mat[is_negative] <- -adj_mat[is_negative]
  
  if (verbose) {
    message(sprintf("Outputting matrix with %d nodes and %d edges", 
                    nrow(adj_mat), sum(adj_mat > 0)))
  }
  
  # Create a weighted or unweighted graph
  if (weighted) {
    g <- graph.adjacency(adj_mat, mode='undirected', weighted=TRUE, diag=FALSE)
  } else {
    adj_mat[adj_mat != 0] <- 1
    g <- graph.adjacency(adj_mat, mode='undirected', diag=FALSE)
  }
  
  # Add node attributes from a vector
  if (!is.null(nodeAttr)) {
    g <- set_vertex_attr(g, "attr", value=nodeAttr)
  }
  
  # Add node attributes from a data frame
  if (!is.null(nodeAttrDataFrame)) {
    for (colname in colnames(nodeAttrDataFrame)) {
      g <- set_vertex_attr(g, colname, value=nodeAttrDataFrame[, colname])
    }
  }
  
  # Save graph to GraphML
  write.graph(g, filename, format='graphml')
  
  # Return the igraph object
  return(g)
}

# Obter anotações de genes
gene_info <- AnnotationDbi::select(Homo.sapiens, keytype='ENTREZID', 
                                   keys=rownames(expr_of_interest),
                                   columns=c('TXCHROM', 'TXSTRAND', 'GENENAME'))

# Renomear as colunas e adicionar informações do módulo
colnames(gene_info) <- c('gene_id', 'chr', 'strand', 'description')
gene_info <- gene_info[!duplicated(gene_info$gene_id), ]
gene_info <- cbind(gene_info, module=module_colors)

# Adicionar cores RGB para o Cytoscape
gene_info$color_rgb <- col2hex(gene_info$module)

# Exportar a rede para GraphML
g <- export_network_to_graphml(adj_matrix, 
                               filename='/Users/carlitos/Desktop/networkprj955.graphml',
                               threshold=0.04, nodeAttrDataFrame=gene_info)

#############################
# Calcular a TOM (Topological Overlap Matrix) para os módulos de interesse
TOM <- TOMsimilarityFromExpr(t(expr_of_interest), power = soft_power)

# Definir nomes das linhas e colunas como IDs dos genes
row.names(TOM) <- row.names(expr_of_interest)
colnames(TOM) <- row.names(expr_of_interest)

# Criar a lista de arestas para visualização em ferramentas como Cytoscape
# Remover prefixo "X" nos nomes dos genes da coluna gene2
edge_list <- data.frame(TOM) %>%
  mutate(gene1 = row.names(.)) %>%
  pivot_longer(-gene1, names_to = "gene2", values_to = "correlation") %>%
  mutate(gene2 = gsub("^X", "", gene2)) %>%  # Remove o "X" dos nomes dos genes
  unique() %>%
  subset(!(gene1 == gene2)) %>%
  mutate(
    module1 = module_df[gene1,]$colors,
    module2 = module_df[gene2,]$colors
  )

# Visualizar a lista de arestas corrigida
head(edge_list)


# Exportar a lista de arestas para um arquivo TSV
write_delim(edge_list, file = "edgelist.tsv", delim = "\t")


##################################################################################
gene_info <- AnnotationDbi::select(Homo.sapiens, keys = edge_list$gene1, 
                                   columns = c("SYMBOL", "TXCHROM", "TXSTART", "TXEND"), 
                                   keytype = "ENTREZID")
view(gene_info)
valores_interesse
edge_list$gene1
unique(edge_list$gene1)

valores_interesse[valores_interesse %in% edge_list$gene1]

# use OrganismDb to retrieve gene annotations
gene_info <- AnnotationDbi::select(Homo.sapiens, keytype='ENTREZID', keys= module_df$gene_id,
                                   columns=c('TXCHROM', 'TXSTRAND', 'GENENAME'))

export_network_to_graphml <- function (adj_mat, filename=NULL, weighted=TRUE,
                                       threshold=0.5, max_edge_ratio=3,
                                       nodeAttr=NULL, nodeAttrDataFrame=NULL,
                                       edgeAttributes=NULL, verbose=FALSE) {
  library('igraph')
  
  # Determine filename to use
  if (is.null(filename)) {
    filename='/Users/carlitos/Desktop/network.graphml'
  }
  
  # TODO 2015/04/09
  # Add option to rescale correlations for each module before applying
  # threshold (this is simpler than the previous approach of trying to
  # determine a different threshold for each module)
  #
  # Still, modules with very low correlations should be given somewhat
  # less priority than those with very high correlations.
  
  #module_colors <- unique(nodeAttrDataFrame$color)
  #module_genes <- which(nodeAttrDataFrame$color == color)
  #module_adjmat <- adj_mat[module_genes,]
  #num_genes <- length(module_genes)
  
  # Adjust threshold if needed to limit remaining edges
  max_edges <- max_edge_ratio * nrow(adj_mat)
  
  edge_to_total_ratio <- max_edges / length(adj_mat)
  edge_limit_cutoff <- as.numeric(quantile(abs(adj_mat), 1 - edge_to_total_ratio))
  
  # Also choose a minimum threshold to make sure that at least some edges
  # are left
  min_threshold <- as.numeric(quantile(abs(adj_mat), 0.9999))
  
  threshold <- min(min_threshold, max(threshold, edge_limit_cutoff))
  
  # Remove edges with weights lower than the cutoff
  adj_mat[abs(adj_mat) < threshold] <- 0
  
  # Drop any genes with no edges (TODO: Make optional)
  orphaned <- (colSums(adj_mat) == 0)
  adj_mat <- adj_mat[!orphaned, !orphaned]
  
  # Also remove annotation entries
  if (!is.null(nodeAttr)) {
    nodeAttr <- nodeAttr[!orphaned]
  }
  if (!is.null(nodeAttrDataFrame)) {
    nodeAttrDataFrame <- nodeAttrDataFrame[!orphaned,]
  }
  
  # Keep track of non-positive edges and rescale to range 0,1
  is_zero     <- adj_mat == 0
  is_negative <- adj_mat < 0
  
  adj_mat <- (abs(adj_mat) - threshold) / (max(adj_mat) - threshold)
  adj_mat[is_zero] <- 0
  adj_mat[is_negative] <- -adj_mat[is_negative]
  
  if (verbose) {
    message(sprintf("Outputting matrix with %d nodes and %d edges", 
                    nrow(adj_mat), sum(adj_mat > 0)))
  }
  
  # Create a new graph and add vertices
  # Weighted graph
  if (weighted) {
    g <- graph.adjacency(adj_mat, mode='undirected', weighted=TRUE, diag=FALSE)
  } else {
    adj_mat[adj_mat != 0] <- 1
    g <- graph.adjacency(adj_mat, mode='undirected', diag=FALSE)
  }
  
  # Add single node annotation from vector
  if (!is.null(nodeAttr)) {
    g <- set.vertex.attribute(g, "attr", value=nodeAttr)
  }
  
  # Add node one or more node annotations from a data frame
  if (!is.null(nodeAttrDataFrame)) {
    for (colname in colnames(nodeAttrDataFrame)) {
      g <- set.vertex.attribute(g, colname, value=nodeAttrDataFrame[,colname])
    }
  }
  
  edge_correlation_negative <- c()
  
  # neg_correlations[edge_list]
  edge_list <- get.edgelist(g)
  
  for (i in 1:nrow(edge_list)) {
    from <- edge_list[i, 1]    
    to   <- edge_list[i, 2]    
  }
  
  # Save graph to a file
  write.graph(g, filename, format='graphml')
  
  # return igraph
  return(g)
}


gene_tree <- hclust(as.dist(1 - TOM), method="average")

# we will use the cuttreeDynamicTree method to break apart the hc dendrogram
# into separate modules
module_labels <- cutreeDynamicTree(dendro=gene_tree, minModuleSize=15,
                                   deepSplit=TRUE)

# assign a color to each module for easier visualization and referencing

module_colors <- labels2colors(module_labels)
module_colors

colnames(gene_info) <- c('gene_id', 'description', 'chr', 'strand')

# for now, just grab the description for the first transcript
gene_info <- gene_info[!duplicated(gene_info$gene_id),]

gene_info <- cbind(gene_info, module=module_colors)

# Include RGB versions of module colors for better assignment in Cytoscape
gene_info$color_rgb <- col2hex(gene_info$module)

# first, it's a good idea to check the distribution of edges weights in our
# correlation matrix. This will help us choose a reasonable cutoff for
# exporting the network.
g <- export_network_to_graphml(adj_matrix, filename='~/network.graphml',
                               threshold=0.4, nodeAttrDataFrame=gene_info)



