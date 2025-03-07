# install.packages("pheatmap")
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
library(ggrepel)


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



################################################################

library(DESeq2)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
run_deseq_analysis <- function(data, phenoData, group1, group2, output_file) {
  # Criar diretório se não existir
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  
  # Filtrando apenas as amostras relevantes
  samples <- phenoData$SampleName[phenoData$Treatment %in% c(group1, group2)]
  
  # Criando subconjunto de dados
  counts <- data[, samples]
  pheno <- phenoData[phenoData$SampleName %in% samples, ]
  
  # Criando objeto DESeqDataSet
  dds <- DESeqDataSetFromMatrix(countData = counts,
                                colData = pheno,
                                design = ~ Treatment)
  
  # Rodando a análise DESeq2
  dds <- DESeq(dds)
  
  # Extraindo resultados
  res <- results(dds, contrast = c("Treatment", group2, group1))
  
  # Salvando resultados em formato tabular
  write.table(as.data.frame(res), file = output_file, sep = "\t", quote = FALSE, row.names = TRUE)
  
  return(res)
}



run_pca_analysis <- function(dds, output_file) {
  rld <- rlog(dds, blind = TRUE)
  pca_data <- plotPCA(rld, intgroup = "Treatment", returnData = TRUE)
  
  # Calculando a variância explicada
  pca_variance <- attr(pca_data, "percentVar")
  
  # Criando o gráfico com rótulos ajustados
  p <- ggplot(pca_data, aes(PC1, PC2, color = Treatment, label = rownames(pca_data))) +
    geom_point(size = 3) +
    geom_text_repel(size = 3, max.overlaps = 10) +  # Evita sobreposição e corte de nomes
    xlab(paste0("PC1: ", round(pca_variance[1] * 100, 2), "% variance")) +
    ylab(paste0("PC2: ", round(pca_variance[2] * 100, 2), "% variance")) +
    theme_minimal()
  
  ggsave(output_file, p, width = 6, height = 5)  # Ajuste no tamanho do gráfico
  
  return(p)
}


run_sample_distances <- function(dds, output_file) {
  rld <- rlog(dds, blind = TRUE)
  sampleDists <- dist(t(assay(rld)))
  sampleDistMatrix <- as.matrix(sampleDists)
  rownames(sampleDistMatrix) <- colnames(dds)
  colnames(sampleDistMatrix) <- colnames(dds)
  
  colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
  pheatmap(sampleDistMatrix, col = colors, filename = output_file)
}


# Executando a função para as comparações desejadas
# res_nf_ilp <- run_deseq_analysis(data, phenoData, "Non-Infected", "Infected_Live_Parasites", "./Deseq2/SRP377060 - Leishmania/results/DESeq2_NonInfected_vs_InfectedLive.tabular")
# res_nf_ifp <- run_deseq_analysis(data, phenoData, "Non-Infected", "Infected_Fixed_Parasites", "./Deseq2/SRP377060 - Leishmania/results/DESeq2_NonInfected_vs_InfectedFixed.tabular")

dds_nf_ilp <- run_deseq_analysis(data, phenoData, "Non-Infected", "Infected_Live_Parasites", "./Deseq2/SRP377060 - Leishmania/results/DESeq2_NonInfected_vs_InfectedLive.tabular")
dds_nf_ifp <- run_deseq_analysis(data, phenoData, "Non-Infected", "Infected_Fixed_Parasites", "./Deseq2/SRP377060 - Leishmania/results/DESeq2_NonInfected_vs_InfectedFixed.tabular")
# Visualizando os primeiros resultados
head(dds_nf_ilp)
head(dds_nf_ifp)



# Executando PCA
pca_nf_ilp <- run_pca_analysis(dds_nf_ilp, "./Deseq2/SRP377060 - Leishmania/results/PCA_NonInfected_vs_InfectedLive.jpeg")
pca_nf_ifp <- run_pca_analysis(dds_nf_ifp, "./Deseq2/SRP377060 - Leishmania/results/PCA_NonInfected_vs_InfectedFixed.jpeg")

# Executando Sample-to-Sample Distances
run_sample_distances(dds_nf_ilp, "./Deseq2/SRP377060 - Leishmania/results/SampleDistances_NonInfected_vs_InfectedLive.jpeg")
run_sample_distances(dds_nf_ifp, "./Deseq2/SRP377060 - Leishmania/results/SampleDistances_NonInfected_vs_InfectedFixed.jpeg")

