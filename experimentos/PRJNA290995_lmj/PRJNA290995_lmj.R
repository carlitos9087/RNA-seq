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
if (!requireNamespace("rentrez", quietly = TRUE)) {
  install.packages("rentrez")
}




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


################################################################################
# Montando o conjunto de dados
################################################################################
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

################################################################################
# 2. QC - outlier detection 
# detect outlier genes
################################################################################
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




#exclude outlier samples
samples.to.be.excluded <- c("Lama_Infected_72h_R3_SRR2163299")
data.subset <- data[,!(colnames(data) %in% samples.to.be.excluded)]


################################################################################
# 3. Normalization ----------------------------------------------------------------------
################################################################################

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



################################################################################
# 4. Network Construction 
# Choose a set of soft-thresholding powers
################################################################################
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
load("/Users/carlitos/Desktop/resultados/PRJNA290995_lmj/bwnet.RData")

################################################################################
# 5. Module Eigengenes ---------------------------------------------------------
################################################################################


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




################################################################################
# 6A. Relate modules to traits --------------------------------------------------
# module trait associations
################################################################################

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
             x = names(heatmap.data)[24:31],
             y = names(heatmap.data)[1:20],rotLabX = 50,
             col = c("cyan", "white", "grey", "purple"))

# write.csv(as.data.frame(heatmap.data), file = "/Users/carlitos/Desktop/PRJNA290995_lama_heatmap.data.csv", row.names = TRUE)

# module.gene.mapping <- as.data.frame(bwnet$colors)
# module.gene.mapping %>% 
#   filter(`bwnet$colors` == 'brown') %>% 
#   rownames()
# 
# class(bwnet$colors)

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



df_colors
df_colors285
df_colors344

df_colors$genes <- as.character(df_colors$genes)

df_colors <- df_colors %>%
  mutate(GeneSymbol = mapIds( org.Hs.eg.db, keys = genes, column = "SYMBOL",
                            keytype = "ENTREZID", multiVals = "first" ))

df_colors <- dplyr::select(df_colors, GeneSymbol, genes, Colors)


df_colors285$genes <- as.character(df_colors285$genes)
df_colors285 <- df_colors285 %>%
  mutate(GeneSymbol = mapIds( org.Hs.eg.db, keys = genes, column = "SYMBOL",
                              keytype = "ENTREZID", multiVals = "first" ))
df_colors285 <- dplyr::select(df_colors285, GeneSymbol, genes, Colors)


df_colors344$genes <- as.character(df_colors344$genes)
df_colors344 <- df_colors344 %>%
  mutate(GeneSymbol = mapIds( org.Hs.eg.db, keys = genes, column = "SYMBOL",
                              keytype = "ENTREZID", multiVals = "first" ))
df_colors344 <- dplyr::select(df_colors344, GeneSymbol, genes, Colors)


#############################################################
library(openxlsx)

# Função para gerar gráficos e salvá-los como imagens
generate_plot <- function(data, colors_col, freq_col, output_file, main_title, x_label, y_label) {
  png(output_file, width = 800, height = 600) # Salva o gráfico como imagem
  par(mar = c(8, 4, 4, 2) + 0.5)
  par(mgp = c(5.5, 1, 0))
  barplot(
    data[[freq_col]],
    names.arg = data[[colors_col]],
    las = 2,
    col = data[[colors_col]],
    main = main_title,
    xlab = x_label,
    ylab = y_label
  )
  dev.off()
}

# Função para inserir gráfico no Excel
insert_plot_to_excel <- function(wb, sheet_name, image_file, start_col, start_row, width, height) {
  insertImage(
    wb,
    sheet = sheet_name,
    file = image_file,
    width = width, height = height,
    startCol = start_col, startRow = start_row
  )
}

# Dados e configurações para os gráficos e abas
plots_info <- list(
  list(
    data = df_table_geral,
    colors_col = "Colors_geral",
    freq_col = "Frequency",
    sheet_name = "colors geral conts",
    output_file = "/Users/carlitos/Desktop/temp_plot1.png",
    main_title = "Frequency of Colors (Geral)",
    x_label = "Colors",
    y_label = "Frequency"
  ),
  list(
    data = df_table,
    colors_col = "Colors",
    freq_col = "Frequency",
    sheet_name = "Color interesse Counts",
    output_file = "/Users/carlitos/Desktop/temp_plot2.png",
    main_title = "Frequency of Colors (Subset)",
    x_label = "Colors",
    y_label = "Frequency"
  ),
  list(
    data = df_table344,
    colors_col = "Colors",
    freq_col = "Frequency",
    sheet_name = "Colors Interesse 344",
    output_file = "/Users/carlitos/Desktop/temp_plot3.png",
    main_title = "Frequency of Colors (344)",
    x_label = "Colors",
    y_label = "Frequency"
  ),
  list(
    data = df_table285,
    colors_col = "Colors",
    freq_col = "Frequency",
    sheet_name = "Colors Interesse 285",
    output_file = "/Users/carlitos/Desktop/temp_plot4.png",
    main_title = "Frequency of Colors (285)",
    x_label = "Colors",
    y_label = "Frequency"
  )
)

# Criação do workbook e adição das abas
wb <- createWorkbook()
addWorksheet(wb, "colors geral conts")
addWorksheet(wb, "Color interesse Counts")
addWorksheet(wb, "Colors Info Interesse")
addWorksheet(wb, "Colors Interesse 344")
addWorksheet(wb, "Colors Info 344")
addWorksheet(wb, "Colors Interesse 285")
addWorksheet(wb, "Colors Info 285")

# Escreve os dados nas abas
writeData(wb, sheet = "colors geral conts", x = df_table_geral, startCol = 1, startRow = 1)
writeData(wb, sheet = "Color interesse Counts", x = df_table, startCol = 1, startRow = 1)
writeData(wb, sheet = "Colors Info Interesse", x = df_colors, startCol = 1, startRow = 1)
writeData(wb, sheet = "Colors Interesse 344", x = df_table344, startCol = 1, startRow = 1)
writeData(wb, sheet = "Colors Info 344", x = df_colors344, startCol = 1, startRow = 1)
writeData(wb, sheet = "Colors Interesse 285", x = df_table285, startCol = 1, startRow = 1)
writeData(wb, sheet = "Colors Info 285", x = df_colors285, startCol = 1, startRow = 1)

# Geração dos gráficos e inserção no Excel
for (plot_info in plots_info) {
  generate_plot(
    data = plot_info$data,
    colors_col = plot_info$colors_col,
    freq_col = plot_info$freq_col,
    output_file = plot_info$output_file,
    main_title = plot_info$main_title,
    x_label = plot_info$x_label,
    y_label = plot_info$y_label
  )
  insert_plot_to_excel(
    wb = wb,
    sheet_name = plot_info$sheet_name,
    image_file = plot_info$output_file,
    start_col = 5,
    start_row = 1,
    width = 8,
    height = 6
  )
}

# Salva o arquivo Excel no local especificado
saveWorkbook(wb, "/Users/carlitos/Desktop/bwnet_colors_PRJNA290995_lmj.xlsx", overwrite = TRUE)

# Remove os arquivos temporários
for (plot_info in plots_info) {
  file.remove(plot_info$output_file)
}


#############################################################


# 6B. Intramodular analysis: Identifying driver genes ---------------

# Calculate the module membership and the associated p-values

# The module membership/intramodular connectivity is calculated as the correlation of the eigengene and the gene expression profile. 
# This quantifies the similarity of all genes on the array to every module.

# module.membership.measure <- cor(module_eigengenes, norm.counts, use = 'p')
# module.membership.measure.pvals <- corPvalueStudent(module.membership.measure, nSamples)
# 
# 
# module.membership.measure.pvals[1:12,1:12]
# 
# 
# # Calculate the gene significance and associated p-values
# 
# gene.signf.corr <- cor(norm.counts, severity.out$data.Lama_Infected_4h.vs.all, use = 'p')
# gene.signf.corr.pvals <- corPvalueStudent(gene.signf.corr, nSamples)
# 
# 
# gene.signf.corr.pvals %>% 
#   as.data.frame() %>% 
#   arrange(V1) %>% 
#   head(25)
# 
# 
# # Using the gene significance you can identify genes that have a high significance for trait of interest 
# # Using the module membership measures you can identify genes with high module membership in interesting modules.
# 
# # 7. analisis enrich
# green_genes <- module.gene.mapping %>% 
#   filter(`bwnet$colors` == 'brown') %>% 
#   rownames()



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

# Carregar as bibliotecas necessárias
library(EnsDb.Hsapiens.v86)
library(Homo.sapiens)
library(AnnotationDbi)
library(org.Hs.eg.db)

# Suponha que 'geneIDs' seja uma lista de IDs de genes a partir do conjunto de dados
geneIDs <- colnames(norm.counts)

# Obter os símbolos dos genes e nomes completos
geneSymbols <- mapIds(org.Hs.eg.db, keys = geneIDs, column = "SYMBOL", keytype = "ENTREZID", multiVals = "first")
geneNames <- mapIds(org.Hs.eg.db, keys = geneIDs, column = "GENENAME", keytype = "ENTREZID", multiVals = "first")

# Criar uma tabela de nodos com cores dos módulos e anotações adicionais
nodeData <- data.frame(
  Node = geneIDs,
  ModuleColor = bwnet$colors,
  GeneSymbol = geneSymbols,
  GeneName = geneNames
)

# Remover genes sem informações de símbolo ou nome, se necessário
nodeData <- na.omit(nodeData)

# Exportar tabela de nodos para uso no Cytoscape
write.table(nodeData, "./experimentos/PRJNA290995_lmj/CytoscapeNodeFile-PRJNA290995_lmj.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)




# TOM = TOMsimilarityFromExpr(norm.counts, power = 16)
# save(TOM, file = "/Users/carlitos/Desktop/resultados/PRJNA290995_lmj/TOM_power16.RData")
load(file = "/Users/carlitos/Desktop/resultados/PRJNA290995_lmj/TOM_power16.RData")


# Definir limiar para TOM
threshold <- 0.17

sum(TOM > threshold)

# Obter os nomes dos genes
genes <- colnames(norm.counts)

# Configurar matriz topológica com nomes dos genes
dimnames(TOM) <- list(genes, genes)

# Identificar conexões que atendem ao limiar
TOM_indices <- which(TOM > threshold, arr.ind = TRUE)
TOM_indices <- TOM_indices[TOM_indices[, 1] < TOM_indices[, 2], ]  # Apenas conexões superiores à diagonal principal

# Criar tabela de arestas para o Cytoscape
edgeData <- data.frame(
  fromNode = genes[TOM_indices[, 1]],
  toNode = genes[TOM_indices[, 2]],
  weight = TOM[TOM_indices],
  direction = "undirected"
)

# Adicionar os símbolos dos genes correspondentes
edgeData$fromAltName <- nodeData$GeneSymbol[match(edgeData$fromNode, nodeData$Node)]
edgeData$toAltName <- nodeData$GeneSymbol[match(edgeData$toNode, nodeData$Node)]

# Exportar tabela de arestas para o Cytoscape
write.table(edgeData, "./experimentos/PRJNA290995_lmj/CytoscapeEdgeFile-PRJNA290995_lmj.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

# Mensagem de conclusão
cat("Arquivos 'CytoscapeEdgeFile-PRJNA290995.txt' e 'CytoscapeNodeFile-PRJNA290995.txt' foram gerados com sucesso.\n")







# 
# 
# # Carregar pacotes necessários
# library(org.Hs.eg.db)
# library(AnnotationDbi)
# 
# # Suponha que 'valores_interesse' contenha os IDs de genes relevantes
# geneIDs <- valores_interesse
# 
# # Obter os símbolos dos genes com base nos IDs
# geneSymbols <- mapIds(org.Hs.eg.db, keys = geneIDs, column = "SYMBOL", keytype = "ENTREZID", multiVals = "first")
# geneNames <- mapIds(org.Hs.eg.db, keys = geneIDs, column = "GENENAME", keytype = "ENTREZID", multiVals = "first")
# 
# # Criar o data frame com cores e símbolos
# geneData <- data.frame(
#   GeneID = geneIDs,
#   ModuleColor = bwnet$colors[geneIDs],  # Usando geneIDs em vez de valores_interesse
#   GeneSymbol = geneSymbols,
#   GeneName = geneNames
# )
# 
# # Remover linhas com valores NA, se necessário
# geneData <- na.omit(geneData)
# 
# # Visualizar o data frame resultante
# print(geneData)
# 
# # Exportar para um arquivo
# write.table(geneData, "genes_com_modulos_PRJNA290995.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)




