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
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("clusterProfiler")
# 
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("pathview")
# 
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("apeglm")
# # Instalar o pacote ggrepel se necessário
# if (!requireNamespace("ggrepel", quietly = TRUE)) {
#   install.packages("ggrepel")
# }

if (!requireNamespace("openxlsx", quietly = TRUE)) {
  install.packages("openxlsx")
}

# Carregar o pacote ggrepel

setwd("/Users/carlitos/Desktop/RNA-seq/")
# setwd("/Users/carlitos/Documents/")

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

tabular_dir  <- "experimentos/fastas/Leishmania/PRJNA290995 --------lmj/"
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

# for ( i in colnames(data.subset)){
#   print(i)
# }

# 3. Normalization ----------------------------------------------------------------------
# create a deseq2 dataset

phenoData  <-  read_excel("experimentos/PRJNA290995_lmj/Phenodata lmj.xlsx", col_names = TRUE)
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

dds75 <- dds[rowSums(counts(dds) >= 15) >= 18,]
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

soft_power <- 18
temp_cor <- cor
cor <- WGCNA::cor


# memory estimate w.r.t blocksize
# bwnet <- blockwiseModules(norm.counts,
                          # maxBlockSize = 14000,
                          # TOMType = "signed",
                          # power = soft_power,
                          # mergeCutHeight = 0.25,
                          # numericLabels = FALSE,
                          # randomSeed = 1234,
                          # verbose = 3)


cor <- temp_cor

# save(bwnet, file = "PRJNA290995_lmj_bwnet.RData")
# load("/Users/carlitos/Desktop/resultados/PRJNA290995_lmj/bwnet.RData")

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
  mutate(disease_state_bin = ifelse(grepl('_Infected_', Treatment), 1, 0)) #%>% select(2)


# print(unique(traits))
# binarize categorical variables
unique(colData$Treatment)

colData$severity <- factor(colData$Treatment, levels = c(
  "Lmj_Uninfected_4h","Lmj_Uninfected_24h","Lmj_Uninfected_48h","Lmj_Uninfected_72h",
  "Lmj_Infected_4h","Lmj_Infected_24h","Lmj_Infected_48h","Lmj_Infected_72h"))

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
             x = names(heatmap.data)[21:27],
             y = names(heatmap.data)[1:20],rotLabX = 50,
             col = c("cyan", "white", "grey", "purple"))

# write.csv(as.data.frame(heatmap.data), file = "/Users/carlitos/Desktop/PRJNA290995_lama_heatmap.data.csv", row.names = TRUE)

module.gene.mapping <- as.data.frame(bwnet$colors)
module.gene.mapping %>% 
  filter(`bwnet$colors` == 'brown') %>% 
  rownames()

class(bwnet$colors)

genes344 = read_excel("/Users/carlitos/Desktop/acetylation_344.xlsx", sheet = 1)
valores_interesse344 = genes344$Entrez_ID
valores_interesse344

genes285 = read_excel("/Users/carlitos/Desktop/acetylation_344.xlsx", sheet = 2)
valores_interesse285 = genes285$Entrez_ID
valores_interesse285

valores_interesse <- c("54915", "51441", "253943", "91746", "79068", "56339",
                       "57721", "54890", "64848", "221120", "8846", "84266",
                       "80312", "51605", "23378", "55006", "115708", "54888",
                       "1787", "55226", "10189", "4904", "8520", "2648",
                       "8850", "10524", "7994", "23522", "11143", "84148",
                       "9329", "1387", "79969", "2033", "3065", "3066",
                       "8841", "9759", "10014", "10013", "51564", "55869",
                       "9734", "83933", "23411", "22933", "23410", "23409",
                       "23408", "51548", "51547")

colors_geral = bwnet$colors


df_table_geral <- as.data.frame(table(colors_geral), stringsAsFactors = FALSE)
colnames(df_table_geral) <- c("Colors_geral", "Frequency")
df_table_geral

colors_interesse = bwnet$colors[valores_interesse]
colors_interesse
bwnet$colors[valores_interesse]
table(bwnet$colors[valores_interesse])


df_table <- as.data.frame(table(colors_interesse), stringsAsFactors = FALSE)
colnames(df_table) <- c("Colors", "Frequency")

df_colors <- data.frame(genes = valores_interesse, Colors = colors_interesse)

################################################################################
colors_interesse344 = bwnet$colors[valores_interesse344]
colors_interesse344
bwnet$colors[valores_interesse344]
table(bwnet$colors[valores_interesse344])


df_table344 <- as.data.frame(table(colors_interesse344), stringsAsFactors = FALSE)
colnames(df_table344) <- c("Colors", "Frequency")

df_colors344 <- data.frame(genes = valores_interesse344, Colors = colors_interesse344)
################################################################################

colors_interesse285 = bwnet$colors[valores_interesse285]
colors_interesse285
bwnet$colors[valores_interesse285]
table(bwnet$colors[valores_interesse285])


df_table285 <- as.data.frame(table(colors_interesse285), stringsAsFactors = FALSE)
colnames(df_table285) <- c("Colors", "Frequency")

df_colors285 <- data.frame(genes = valores_interesse285, Colors = colors_interesse285)




library(writexl)
write_xlsx(list("colors geral conts" = df_table_geral, "Color Counts" = df_table, "Colors Info" = df_colors), path = "/Users/carlitos/Desktop/bwnet_colors.xlsx")
#############################################################
# Instale os pacotes, se necessário
if (!requireNamespace("openxlsx", quietly = TRUE)) {
  install.packages("openxlsx")
}

library(openxlsx)

# Criação do workbook
wb <- createWorkbook()

# Adiciona uma aba para "colors geral conts"
addWorksheet(wb, "colors geral conts")

# Adiciona a tabela "df_table_geral" na aba "colors geral conts"
writeData(wb, sheet = "colors geral conts", x = df_table_geral, startCol = 1, startRow = 1)

# Criação do gráfico em R
png("/Users/carlitos/Desktop/temp_plot.png", width = 800, height = 600) # Salva o gráfico como imagem temporária
cores_graf = df_table_geral$Colors_geral
barplot(
  df_table_geral$Frequency,
  names.arg = df_table_geral$Colors_geral,
  las = 2, # Rotação dos rótulos
  col = cores_graf,
  main = "Frequency of Colors",
  xlab = "Colors",
  ylab = "Frequency"
)
dev.off()

# Insere o gráfico como imagem na aba "colors geral conts"
insertImage(
  wb,
  sheet = "colors geral conts",
  file = "/Users/carlitos/Desktop/temp_plot.png",
  width = 10, height = 6,
  startCol = 5, startRow = 1
)

# Adiciona as outras tabelas em abas separadas
writeData(wb, sheet = addWorksheet(wb, "Color Counts"), x = df_table, startCol = 1, startRow = 1)
writeData(wb, sheet = addWorksheet(wb, "Colors Info"), x = df_colors, startCol = 1, startRow = 1)
writeData(wb, sheet = addWorksheet(wb, "Colors Interesse 344"), x = df_table344, startCol = 1, startRow = 1)
writeData(wb, sheet = addWorksheet(wb, "Colors Info 344"), x = df_colors344, startCol = 1, startRow = 1)
writeData(wb, sheet = addWorksheet(wb, "Colors Interesse 285"), x = df_table285, startCol = 1, startRow = 1)
writeData(wb, sheet = addWorksheet(wb, "Colors Info 285"), x = df_colors285, startCol = 1, startRow = 1)

# Salva o arquivo Excel no local especificado
saveWorkbook(wb, "/Users/carlitos/Desktop/bwnet_colors1.xlsx", overwrite = TRUE)

# Remove o arquivo temporário
file.remove("/Users/carlitos/Desktop/temp_plot.png")



#############################################################


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


gene.signf.corr.pvals %>% 
  as.data.frame() %>% 
  arrange(V1) %>% 
  head(25)


# Using the gene significance you can identify genes that have a high significance for trait of interest 
# Using the module membership measures you can identify genes with high module membership in interesting modules.

# 7. analisis enrich
green_genes <- module.gene.mapping %>% 
  filter(`bwnet$colors` == 'brown') %>% 
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
# 
# library(clusterProfiler)
# library(org.Hs.eg.db)  # Use a base de dados apropriada para o seu organismo
# 
# # Supondo que os IDs dos genes sejam Entrez IDs
# ego <- enrichGO(gene          = green_genes,
#                 OrgDb         = org.Hs.eg.db,
#                 keyType       = 'ENTREZID',  # Ajuste o keyType conforme necessário
#                 ont           = "ALL",  # Pode ser "BP", "MF", "CC", ou "ALL"
#                 pAdjustMethod = "BH",
#                 qvalueCutoff  = 0.05,
#                 readable      = TRUE)
# 
# # Visualize os resultados
# head(ego)
# 
# library(enrichplot)
# 
# # Dotplot
# dotplot(ego)
# 
# # Barplot
# barplot(ego, showCategory = 10)
# 
# 
# # Supondo que o objeto ego já está criado
# ego <- pairwise_termsim(ego)
# 
# library(ggplot2)
# 
# emap <- emapplot(ego)
# # Ajustar o tamanho dos labels dos termos para um tamanho menor
# emap + theme(
#   legend.text = element_text(size = 8),
#   legend.title = element_text(size = 8),
#   plot.title = element_text(size = 10, face = "bold"),
#   axis.text.x = element_text(size = 8),
#   axis.text.y = element_text(size = 8)
# )
# 
# 
# 
# ego_df <- as.data.frame(ego)
# write.csv(ego_df, file = "/Users/carlitos/Desktop/enrichment_results.csv", row.names = FALSE)



# Calcular o TOM (Topological Overlap Matrix)
# TOM <- TOMsimilarityFromExpr(norm.counts, power = soft_power, TOMType = "signed")

# save(TOM, file = "TOM_signed.RData")
load("TOM_signed.RData")

# Carregar as bibliotecas necessárias
library(Homo.sapiens)
library(AnnotationDbi)
library(org.Hs.eg.db)

# Suponha que 'geneIDs' seja uma lista de IDs de genes a partir do seu conjunto de dados
geneIDs <- colnames(norm.counts)

# Obter os símbolos dos genes e nomes completos
geneSymbols <- mapIds(org.Hs.eg.db, keys = geneIDs, column = "SYMBOL", keytype = "ENTREZID", multiVals = "first")
geneNames <- mapIds(org.Hs.eg.db, keys = geneIDs, column = "GENENAME", keytype = "ENTREZID", multiVals = "first")



# Criar uma tabela de nodos com cores dos módulos e anotações adicionais
nodeData <- data.frame(Node = geneIDs,
                       ModuleColor = bwnet$colors,
                       GeneSymbol = geneSymbols,
                       GeneName = geneNames)


# Remover quaisquer genes sem informações de símbolo ou nome, se necessário
nodeData <- na.omit(nodeData)

# Exportar para um arquivo que será usado no Cytoscape
write.table(nodeData, "CytoscapeInput-nodes1.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)


# Contar o número de conexões com TOM > 0.09
num_connections_above_threshold <- sum(TOM > 0.148)

cat("Número de conexões com TOM > 0.09:", num_connections_above_threshold, "\n")




# Carregar pacotes necessários
library(org.Hs.eg.db)
library(AnnotationDbi)

# Suponha que 'valores_interesse' contenha os IDs de genes relevantes
geneIDs <- valores_interesse

# Obter os símbolos dos genes com base nos IDs
geneSymbols <- mapIds(org.Hs.eg.db, keys = geneIDs, column = "SYMBOL", keytype = "ENTREZID", multiVals = "first")
geneNames <- mapIds(org.Hs.eg.db, keys = geneIDs, column = "GENENAME", keytype = "ENTREZID", multiVals = "first")

# Criar o data frame com cores e símbolos
geneData <- data.frame(
  GeneID = geneIDs,
  ModuleColor = bwnet$colors[geneIDs],  # Usando geneIDs em vez de valores_interesse
  GeneSymbol = geneSymbols,
  GeneName = geneNames
)

# Remover linhas com valores NA, se necessário
geneData <- na.omit(geneData)

# Visualizar o data frame resultante
print(geneData)

# Exportar para um arquivo
write.table(geneData, "genes_com_modulos.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)



# Dados fornecidos
color_table <- c(black = 95, blue = 1153, brown = 770, green = 120, greenyellow = 20, 
                 grey = 6345, magenta = 57, pink = 60, purple = 29, red = 100, 
                 turquoise = 3023, yellow = 502)

# Definir as cores correspondentes
colors <- c("black", "blue", "brown", "green", "greenyellow", "grey", "magenta", 
            "pink", "purple", "red", "turquoise", "yellow", "")

# Criar o gráfico de barras
barplot(color_table, 
        main = "Distribuição de Cores", 
        xlab = "Cores", 
        ylab = "Frequência", 
        col = colors, 
        las = 2, 
        cex.names = 0.8)
