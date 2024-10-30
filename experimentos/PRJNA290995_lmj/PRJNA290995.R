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


library(knitr)
library(limma)
library(reshape2)
library(RColorBrewer)
library(ggplot2)
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


tabular_dir  <- "/Users/carlitos/Desktop/experimentos/exps fasta/fastas/SARS-CoV-2/GSE163959(homo sapiens)/"
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
phenoData  <-  read_excel("/Users/carlitos/Desktop/GSE163959 sarscov.xlsx", col_names = TRUE)
lista  <- phenoData$SampleName
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

# exclude outlier samples
samples.to.be.excluded <- c("Lama_Infected_72h_R3_SRR2163299")
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
log_counts <- log2(data + 1)
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

soft_power <- 24
temp_cor <- cor
cor <- WGCNA::cor
print("das")

'''
# memory estimate w.r.t blocksize
bwnet <- blockwiseModules(norm.counts,
                          maxBlockSize = 14100,
                          TOMType = "signed",
                          power = soft_power,
                          mergeCutHeight = 0.25,
                          numericLabels = FALSE,
                          randomSeed = 1234,
                          verbose = 3)
'''
bwnet = blockwiseModules(norm.counts, 
                         maxBlockSize = 14100,
                         power = soft_power, 
                         TOMType = "signed", 
                         minModuleSize = 100,
                         mergeCutHeight = 0.25,
                         numericLabels = FALSE,
                         saveTOMs = TRUE,
                         randomSeed = 1234,
                         saveTOMFileBase = "SpodopteraTOM-blockwise3",
                         verbose = 3)


save(bwnet, file = "/Users/carlitos/Desktop/resultados/GSE163959/GSE163959_3_bwnet.RData")
print('aaa')
#Loading the results pf single-block analysis
load(file = "SpodopteraTOM-blockwise2-block.1.RData")
# Re-labeling blockwise modules
#bwLabels = matchLabels(bwnet$colors, moduleLabels)

head(TOM)
#Converting the labels in colors to plot
#bwModuleColors = labels2colors(bwnet$colors)
# Substitui os números pelos nomes das cores diretamente no objeto bwnet
#bwnet$colors = labels2colors(bwnet$colors)
#save(bwnet, file = "/Users/carlitos/Desktop/resultados/GSE163959/GSE163959_2_bwnet.RData")
cor <- temp_cor

print("blabla")
#load("/Users/carlitos/Desktop/resultados/GSE163959/GSE163959_2_bwnet.RData")

# 5. Module Eigengenes ---------------------------------------------------------

load("/Users/carlitos/Desktop/resultados/GSE163959/GSE163959_3_bwnet.RData")
unique(bwnet$colors)


module_eigengenes <- bwnet$MEs
module_eigengenes

# Print out a preview
head(module_eigengenes)

bwnet$colors
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


# 6A. Relate modules to traits --------------------------------------------------
# module trait associations


# create traits file - binarize categorical variables
traits <- colData %>% 
  mutate(disease_state_bin = ifelse(grepl("SARS-", Treatment), 1, 0))

# print(unique(traits))
# binarize categorical variables
unique(colData$Treatment)

colData$severity <- factor(colData$Treatment, levels = c("SARS-CoV-2", "Uninfected" ))

severity.out <- binarizeCategoricalColumns(colData$severity,
                                           includePairwise = FALSE,
                                           includeLevelVsAll = TRUE,
                                           minCount = 1)

traits <- cbind(traits, severity.out)
rownames(traits) <- traits$SampleName
traits <- traits %>% select(-1,-2,-3)
#traits <- traits %>% select(-2,-3)

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


heatmap.data
colnames(heatmap.data)

colnames(traits)
colnames(heatmap.data)
CorLevelPlot(heatmap.data,
             x = names(heatmap.data)[14:15],
             y = names(heatmap.data)[1:13],rotLabX = 0,
             col = c("blue1", "skyblue", "white", "pink", "red"))

# write.csv(as.data.frame(heatmap.data), file = "/Users/carlitos/Desktop/heatmap.data.csv", row.names = TRUE)
###############################################################
# Calcular correlação e valores-p
module.trait.corr <- cor(module_eigengenes, traits, use = 'p')
module.trait.corr.pvals <- corPvalueStudent(module.trait.corr, nSamples)
# Preparar a matriz de texto para exibir no heatmap
textMatrix <- paste(signif(module.trait.corr, 2), "\n(",
                    signif(module.trait.corr.pvals, 1), ")", sep = "")

dim(textMatrix) <- dim(module.trait.corr)
# Definir as margens do gráfico
par(mar = c(6, 8.5, 3, 3))

# Gerar o heatmap
labeledHeatmap(Matrix = module.trait.corr,
               xLabels = names(traits),
               yLabels = names(module_eigengenes),
               ySymbols = names(module_eigengenes),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),  # Escala de cores
               textMatrix = textMatrix,    # Valores de correlação e p-valores
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),             # Limites da escala de cores
               main = "Module-Trait Relationships")


# Ajustar margens do gráfico
par(mar = c(8, 12, 2, 3))  # Aumenta a margem inferior para acomodar labels giradas

# Gerar o heatmap sem girar as labels do eixo X
labeledHeatmap(Matrix = module.trait.corr,
               xLabels = names(traits),
               yLabels = names(module_eigengenes),
               ySymbols = names(module_eigengenes),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.7,
               cex.lab = 1.2,
               zlim = c(-1, 1),
               main = "Module-Trait Relationships",
               xLabelsAngle = 0)  # Define o ângulo das labels como 0 (omitir labels no eixo X)



# Função para converter p-valores em asteriscos
pval_to_stars <- function(pval) {
  if (pval < 0.001) {
    return("***")
  } else if (pval < 0.01) {
    return("**")
  } else if (pval < 0.05) {
    return("*")
  } else {
    return("")
  }
}

# Aplicar a função aos p-valores para gerar a matriz de asteriscos
asteriskMatrix <- apply(module.trait.corr.pvals, c(1, 2), pval_to_stars)

# Combinar valores de correlação com os asteriscos no texto do heatmap
textMatrix <- paste(signif(module.trait.corr, 2), asteriskMatrix, sep = " ")
dim(textMatrix) <- dim(module.trait.corr)

# Ajustar margens do gráfico para acomodar labels giradas
par(mar = c(8, 12, 2, 3))  

# Gerar o heatmap com asteriscos indicando os p-valores
labeledHeatmap(Matrix = module.trait.corr,
               xLabels = names(traits),
               yLabels = names(module_eigengenes),
               ySymbols = names(module_eigengenes),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,  # Matriz de texto com correlações e asteriscos
               setStdMargins = FALSE,
               cex.text = 1,           # Tamanho do texto dentro do heatmap
               cex.lab = 1.2,            # Tamanho das labels dos eixos
               zlim = c(-1, 1),
               main = "Module-Trait Relationships",
               xLabelsAngle = 0)  # Gira as labels do eixo X em 0 graus (default)

# Adicionar labels do eixo X rotacionadas
xLabels <- names(traits)
axis(1, at = 1:length(xLabels), labels = FALSE)  # Cria a base das labels sem texto
text(x = 1:length(xLabels), y = par("usr")[3] - 0.5, srt = 45, adj = 1, 
     labels = xLabels, xpd = TRUE, cex = 1.2)  # Adiciona labels com rotação de 45 graus

###############################################################
uninfecteds = as.data.frame(traits$data.Uninfected.vs.all)
names(uninfecteds) = "Unfecteds"


modNames = substring(names(bwnet$MEs), 3)

geneModuleMembership = as.data.frame(cor(norm.counts, bwnet$MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))

names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")




geneTraitSignificance = as.data.frame(cor(norm.counts, uninfecteds, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))

names(geneTraitSignificance) = paste("GS.", names(uninfecteds), sep="")
names(GSPvalue) = paste("p.GS.", names(uninfecteds), sep="")

module = "greenyellow"
column = match(module, modNames)
moduleGenes = bwnet$colors == module
bwnet$colors == module


moduleGenes[1]

geneModuleMembership[moduleGenes, column]
(geneModuleMembership)

verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for Uninfecteds",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

length(rownames(geneModuleMembership))
# Ajustar margens
par(mfrow = c(1, 1), mar = c(5, 5, 4, 2))

# Criar o scatter plot básico com eixos e título
plot(abs(geneModuleMembership[moduleGenes, column]),
     abs(geneTraitSignificance[moduleGenes, 1]),
     xlab = paste("Module Membership in", module, "module"),
     ylab = "Gene significance for Peso 10 dias",
     main = "Module membership vs. gene significance",
     col = adjustcolor(module, alpha.f = 0.5),  # Ajustar a transparência
     pch = 21,         # Tipo de símbolo (círculo preenchido)
     bg = adjustcolor("greenyellow", alpha.f = 0.5),  # Cor de preenchimento
     cex = 1,          # Tamanho dos pontos
     lwd = 1)          # Largura das bordas (grossura)

# Personalização extra de círculos
points(abs(geneModuleMembership[moduleGenes, column]),
       abs(geneTraitSignificance[moduleGenes, 1]),
       pch = 21,                             # Tipo de símbolo: círculo preenchido
       col = "darkgreen",                     # Cor do contorno do círculo
       bg = adjustcolor("lightgreen", 0.6),   # Cor de preenchimento com transparência
       cex = 1,                            # Tamanho dos círculos
       lwd = 1)                              # Largura do contorno



geneInfo <- data.frame(
  Gene = rownames(geneModuleMembership)[moduleGenes],  # Nomes dos genes
  ModuleMembership = geneModuleMembership[moduleGenes, column],  # MM no módulo
  MMPvalue = MMPvalue[moduleGenes, column],  # P-valor da MM
  GeneSignificance = geneTraitSignificance[moduleGenes, 1],  # GS para 'Uninfecteds'
  GSPvalue = GSPvalue[moduleGenes, 1]  # P-valor da GS
)
rownames(geneModuleMembership)[moduleGenes]

geneModuleMembership[moduleGenes, column]
geneInfo
# Exibindo os primeiros genes
head(geneInfo)

# Definindo critérios de relevância: alta filiação ao módulo e alta significância gênica
importantGenes <- geneInfo %>%
  filter(abs(ModuleMembership) > 0.7 & abs(GeneSignificance) > 0.5) %>%  # Ajuste os thresholds conforme necessário
  arrange(desc(abs(GeneSignificance)))  # Ordena pelos genes mais significativos

# Exibindo os genes mais importantes
head(importantGenes)

# Selecionar e visualizar os genes mais relevantes (top 10)
topGenes <- importantGenes %>% head(10)
topGenes


mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")


entrez_ids <- topGenes$Gene

# Obtenha os nomes dos genes
infos <- getBM(attributes = c("entrezgene_id", "hgnc_symbol"),
                   filters = "entrezgene_id",
                   values = entrez_ids,
                   mart = mart)

# Combine os dados
topGenes <- merge(topGenes, infos, by.x = "Gene", by.y = "entrezgene_id", all.x = TRUE)

# Renomeie a coluna para algo mais intuitivo
names(topGenes)[names(topGenes) == "hgnc_symbol"] <- "GeneName"



# Exibir a lista dos genes mais importantes
print(topGenes)

# Scatter plot dos genes mais importantes
library(ggplot2)

ggplot(topGenes, aes(x = abs(ModuleMembership), y = abs(GeneSignificance))) +
  geom_point(color = "blue", size = 4, alpha = 0.7) +
  geom_text(aes(label = GeneName), vjust = -1, size = 3.5) +  # Rótulos com os nomes dos genes
  labs(x = "Module Membership",
       y = "Gene Significance",
       title = "Top Genes: Module Membership vs. Gene Significance") +
  theme_minimal()

##############################################################
genes = colnames(norm.counts)
moduleColors = bwnet$colors
#if you want export specific colors, substitute the second modulecolors by above modules

inModule = is.finite(match(moduleColors, moduleColors))
modGenes = genes[inModule]

dist(TOM)
modTOM = TOM[inModule, inModule]
dimnames(modTOM) = list(modGenes, modGenes)
modTOMSignificantes = which(modTOM>0.4)

###############################################################


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

###########################################
# Supondo que você tenha os IDs dos genes em um vetor chamado geneIDs
geneIDs <- names(bwnet$colors)  # ou qualquer vetor de IDs que você tenha

# Crie o dataframe
df <- data.frame(ID = geneIDs, Color = bwnet$colors)

# Visualize o dataframe
head(df)

mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Liste todos os atributos disponíveis
#attributes <- listAttributes(mart)
#head(attributes)


entrez_ids <- df$ID

# Obtenha os nomes dos genes
gene_info <- getBM(attributes = c("entrezgene_id", "hgnc_symbol"),
                   filters = "entrezgene_id",
                   values = entrez_ids,
                   mart = mart)

# Combine os dados
df <- merge(df, gene_info, by.x = "ID", by.y = "entrezgene_id", all.x = TRUE)

# Renomeie a coluna para algo mais intuitivo
names(df)[names(df) == "hgnc_symbol"] <- "GeneName"

# Visualize o dataframe atualizado
head(df)

df_filtrado <- df[df$ID %in% valores_interesse, ]

# Visualize o dataframe filtrado
head(df_filtrado)
# Supondo que seu dataframe se chame df_filtrado
write.csv(df_filtrado, "/Users/carlitos/Desktop/resultados/GSE163959/genes_de_filtrado.csv", row.names = FALSE)

################################################



#Exporting the network to a cytoscape format
#Recalculating topological overlap, if necessary
TOM = 0
TOM = TOMsimilarityFromExpr(norm.counts, power = 16)
save(TOM, file = "/Users/carlitos/Desktop/resultados/GSE163959/TOM.RData")
load(file = "/Users/carlitos/Desktop/resultados/GSE163959/TOM.RData")
#Select the modules
#modules = c("greenyellow"); #chose modules that u want to export
#Select the gene modules
modules = c("blue", "greenyellow")
genes = colnames(norm.counts)




# Carregar as cores dos módulos para cada gene
moduleColors <- bwnet$colors  # Supondo que bwnet já contenha os dados de módulos

# Selecionar genes que pertencem ao módulo "greenyellow"
genes = colnames(norm.counts)
inModule = moduleColors == "blue"
modGenes = genes[inModule]

# Selecionar o TOM correspondente aos genes do módulo
modTOM = TOM[inModule, inModule]
dimnames(modTOM) = list(modGenes, modGenes)

# Verificar a dimensão do TOM extraído
dim(modTOM)
modTOM

# Exportar a rede para o formato do Cytoscape
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = "CytoscapeEdgeFile.txt",
                               nodeFile = "CytoscapeNodeFile.txt",
                               weighted = TRUE,
                               threshold = 0.04,  # Ajustar o threshold conforme necessário
                               nodeNames = modGenes,
                               nodeAttr = moduleColors[inModule])




# Definir o limiar de similaridade topológica (por exemplo, 0.4)
threshold = 0.04

# Filtrar o TOM aplicando o threshold (deixando apenas valores acima do limiar)
filteredTOM = TOM
filteredTOM[filteredTOM < threshold] = 0  # Zerar valores abaixo do limiar

# Obter os nomes dos genes
genes = colnames(norm.counts)

# Exportar a rede para o formato do Cytoscape, usando o TOM filtrado
cyt = exportNetworkToCytoscape(filteredTOM,
                               edgeFile = "FilteredCytoscapeEdgeFile.txt",
                               nodeFile = "FilteredCytoscapeNodeFile.txt",
                               weighted = TRUE,
                               threshold = threshold,  # Limiar já aplicado, mas pode ser ajustado aqui também
                               nodeNames = genes,
                               nodeAttr = moduleColors)


print("cabouse")
################################################

cores_interesse <- c( "turquoise")
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

gene.signf.corr <- cor(norm.counts, severity.out$data.Mock.vs.all, use = 'p')
gene.signf.corr.pvals <- corPvalueStudent(gene.signf.corr, nSamples)


gene.signf.corr.pvals %>% 
  as.data.frame() %>% 
  arrange(V1) %>% 
  head(25)


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
modules_of_interest = c( "turquoise", "brown", "blue")

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
black_genes <- names(bwnet$colors[bwnet$colors == "turquoise"])

ego <- enrichGO(gene          = black_genes,
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
# Converte o resultado do enriquecimento para um data frame
enrichment_results <- as.data.frame(ego)

# Salva como CSV
write.csv(enrichment_results, file = "/Users/carlitos/Desktop/turquoise kegg enrichment_results.csv", row.names = FALSE)

#####
ekegg <- enrichKEGG(gene         = black_genes,
                    organism     = 'hsa',       # Substitua pelo organismo apropriado
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.05,
                    qvalueCutoff = 0.05)

# Visualizar resultados
dotplot(ekegg)
barplot(ekegg, showCategory = 10)

ego <- pairwise_termsim(ekegg)

emap <- emapplot(ekegg)
emap + theme(
  legend.text = element_text(size = 8),
  legend.title = element_text(size = 8),
  plot.title = element_text(size = 10, face = "bold"),
  axis.text.x = element_text(size = 8),
  axis.text.y = element_text(size = 8)
)
kegg_results <- as.data.frame(ekegg)

# Salva os resultados em um arquivo CSV
write.csv(kegg_results, file = "/Users/carlitos/Desktop/turquoise kegg_enrichment_results.csv", row.names = FALSE)
#####################################################################################'
black_genes <- names(bwnet$colors[bwnet$colors == "brown"])

ego <- enrichGO(gene          = black_genes,
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
# Converte o resultado do enriquecimento para um data frame
enrichment_results <- as.data.frame(ego)

# Salva como CSV
write.csv(enrichment_results, file = "/Users/carlitos/Desktop/brown enrichment_results.csv", row.names = FALSE)

#####
ekegg <- enrichKEGG(gene         = black_genes,
                    organism     = 'hsa',       # Substitua pelo organismo apropriado
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.05,
                    qvalueCutoff = 0.05)

# Visualizar resultados
dotplot(ekegg)
barplot(ekegg, showCategory = 10)

ekegg <- pairwise_termsim(ekegg)

emap <- emapplot(ekegg)
emap + theme(
  legend.text = element_text(size = 8),
  legend.title = element_text(size = 8),
  plot.title = element_text(size = 10, face = "bold"),
  axis.text.x = element_text(size = 8),
  axis.text.y = element_text(size = 8)
)
kegg_results <- as.data.frame(ekegg)

# Salva os resultados em um arquivo CSV
write.csv(kegg_results, file = "/Users/carlitos/Desktop/brown kegg_enrichment_results.csv", row.names = FALSE)
#####################################################################################
blue_genes <- names(bwnet$colors[bwnet$colors == "blue"])

ego <- enrichGO(gene          = blue_genes,
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
# Converte o resultado do enriquecimento para um data frame
enrichment_results <- as.data.frame(ego)

# Salva como CSV
write.csv(enrichment_results, file = "/Users/carlitos/Desktop/blue enrichment_results.csv", row.names = FALSE)

#####
ekegg <- enrichKEGG(gene         = blue_genes,
                    organism     = 'hsa',       # Substitua pelo organismo apropriado
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.05,
                    qvalueCutoff = 0.05)

# Visualizar resultados
dotplot(ekegg)
barplot(ekegg, showCategory = 10)

ekegg <- pairwise_termsim(ekegg)

emap <- emapplot(ekegg)
emap + theme(
  legend.text = element_text(size = 8),
  legend.title = element_text(size = 8),
  plot.title = element_text(size = 10, face = "bold"),
  axis.text.x = element_text(size = 8),
  axis.text.y = element_text(size = 8)
)
kegg_results <- as.data.frame(ekegg)

# Salva os resultados em um arquivo CSV
write.csv(kegg_results, file = "/Users/carlitos/Desktop/blue kegg_enrichment_results.csv", row.names = FALSE)

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



