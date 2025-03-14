# Define o diretório de trabalho
setwd("/Users/carlitos/Desktop/RNA-seq/")

# Carregando os pacotes necessários
library(dplyr)
library(ggrepel)
library(biomaRt)
library(tidyr)
library(WGCNA)
library(DESeq2)
library(GEOquery)
library(tidyverse)
library(CorLevelPlot)
library(gridExtra)
library(stringr)
library(readxl)
library(openxlsx)
library(pheatmap)
library(RColorBrewer)

###########################################
# 0. Remoção de amostras indesejadas
###########################################
# Lista de amostras que serão removidas
samples_to_remove <- c("DC_Infected_Live_Parasites_3_SRR19400242", 
                       "DC_Infected_Fixed_Parasites_3_SRR19400239")

###########################################
# 1. Leitura e combinação dos arquivos tabulares
###########################################
tabular_dir  <- "Deseq2/SRP377060 - Leishmania/"
tabular_files <- list.files(path = tabular_dir, pattern = "\\.tabular$", full.names = TRUE)

read_tabular_file <- function(file) {
  df <- read.table(file, header = TRUE, stringsAsFactors = FALSE, sep="\t")
  # Renomeia a coluna de contagem para o nome do arquivo (sem extensão)
  colnames(df)[2] <- tools::file_path_sans_ext(basename(file))
  return(df)
}

tabular_dfs <- lapply(tabular_files, read_tabular_file)
combined_df <- Reduce(function(x, y) full_join(x, y, by = "Geneid"), tabular_dfs)
head(combined_df)
rownames(combined_df) <- combined_df$Geneid
combined_df <- combined_df[,-1]

# Converte os dados para uma matriz numérica
data <- as.matrix(combined_df)
storage.mode(data) <- "numeric"

# Remove as amostras indesejadas da matriz de contagem
data <- data[, !(colnames(data) %in% samples_to_remove)]
print(colnames(data))

###########################################
# 2. Leitura dos dados fenotípicos (phenoData)
###########################################
# Supondo que o arquivo Excel contenha as colunas "id", "SampleName" e "Treatment"
phenoData <- read_excel("Deseq2/SRP377060 - Leishmania/SRP377060.xlsx", col_names = TRUE)
lista <- phenoData$id
phenoData <- phenoData[,-1]    # Remove a coluna "id"
rownames(phenoData) <- lista    # Define os nomes das linhas como os ids

# Remove as amostras indesejadas do phenoData
phenoData <- phenoData[!(phenoData$SampleName %in% samples_to_remove), ]
head(phenoData)

###########################################
# 3. QC - Detecção de outliers com WGCNA
###########################################
gsg <- goodSamplesGenes(t(data))
summary(gsg)
print(gsg$allOK)
table(gsg$goodGenes)
table(gsg$goodSamples)

# Remove os genes identificados como outliers
data <- data[gsg$goodGenes, ]

###########################################
# 4. Funções para análise DESeq2, PCA, Sample-to-Sample Distances,
#    Dispersion estimates, Histogram of p-values e MA-plot
###########################################
# Função DESeq2: análise e salvamento dos resultados
run_deseq_analysis <- function(data, phenoData, group1, group2, output_file) {
  # Cria o diretório de saída, se não existir
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  
  # Seleciona as amostras que pertencem aos grupos de interesse e que também estejam nos dados
  valid_samples <- intersect(phenoData$SampleName[phenoData$Treatment %in% c(group1, group2)], colnames(data))
  if(length(valid_samples) == 0){
    stop("Nenhuma amostra encontrada para os grupos especificados.")
  }
  
  # Subconjunto dos dados de contagem e dos metadados
  counts <- data[, valid_samples, drop = FALSE]
  pheno <- phenoData[phenoData$SampleName %in% valid_samples, ]
  
  # Cria o objeto DESeqDataSet
  dds <- DESeqDataSetFromMatrix(countData = counts,
                                colData = pheno,
                                design = ~ Treatment)
  
  # Executa a análise DESeq2
  dds <- DESeq(dds)
  
  # Extrai os resultados do contraste especificado
  res <- results(dds, contrast = c("Treatment", group2, group1))
  
  # Salva os resultados em formato tabular
  write.table(as.data.frame(res), file = output_file, sep = "\t", quote = FALSE, row.names = TRUE)
  
  # Retorna o objeto DESeqDataSet para análises subsequentes
  return(dds)
}

# Função para PCA com alta resolução (usando ggsave com dpi maior)
run_pca_analysis <- function(dds, output_file) {
  rld <- rlog(dds, blind = TRUE)
  pca_data <- plotPCA(rld, intgroup = "Treatment", returnData = TRUE)
  pca_variance <- attr(pca_data, "percentVar")
  
  p <- ggplot(pca_data, aes(PC1, PC2, color = Treatment, label = rownames(pca_data))) +
    geom_point(size = 3) +
    geom_text_repel(size = 3, max.overlaps = 10) +
    xlab(paste0("PC1: ", round(pca_variance[1] * 100, 2), "% variance")) +
    ylab(paste0("PC2: ", round(pca_variance[2] * 100, 2), "% variance")) +
    theme_minimal()
  
  ggsave(output_file, p, width = 6, height = 5, dpi = 300)
  return(p)
}

# Função para Sample-to-Sample Distances com alta resolução
run_sample_distances <- function(dds, output_file) {
  rld <- rlog(dds, blind = TRUE)
  sampleDists <- dist(t(assay(rld)))
  sampleDistMatrix <- as.matrix(sampleDists)
  rownames(sampleDistMatrix) <- colnames(dds)
  colnames(sampleDistMatrix) <- colnames(dds)
  
  colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
  pheatmap(sampleDistMatrix, col = colors, filename = output_file, 
           width = 8, height = 6, units = "in", res = 300)
}

# Função para o gráfico de Dispersion Estimates com alta resolução
run_dispersion_plot <- function(dds, output_file) {
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  png(filename = output_file, width = 8, height = 6, units = "in", res = 300)
  plotDispEsts(dds)
  dev.off()
}

# Função para gerar o Histograma de p-values com alta resolução
run_pval_histogram <- function(dds, output_file) {
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  res <- results(dds)
  png(filename = output_file, width = 8, height = 6, units = "in", res = 300)
  hist(res$pvalue,
       breaks = 50,
       col = "skyblue",
       main = "Histogram of p-values",
       xlab = "p-value")
  dev.off()
}

# Função para gerar o MA-plot com alta resolução e escala de -8 a 8
run_ma_plot <- function(dds, output_file) {
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  res <- results(dds)
  png(filename = output_file, width = 8, height = 6, units = "in", res = 300)
  plotMA(res, main = "MA-plot", ylim = c(-8, 8))
  dev.off()
}


###########################################
# 5. Execução das análises
###########################################

# Executa DESeq2 para as comparações desejadas
dds_nf_ilp <- run_deseq_analysis(data, phenoData, "Non-Infected", "Infected_Live_Parasites", 
                                 "./Deseq2/SRP377060 - Leishmania/results/DESeq2_NonInfected_vs_InfectedLive.tabular")
dds_nf_ifp <- run_deseq_analysis(data, phenoData, "Non-Infected", "Infected_Fixed_Parasites", 
                                 "./Deseq2/SRP377060 - Leishmania/results/DESeq2_NonInfected_vs_InfectedFixed.tabular")

# Visualiza as contagens (primeiras linhas) dos objetos DESeqDataSet
head(assay(dds_nf_ilp))
head(assay(dds_nf_ifp))

# Executa e salva os gráficos de PCA
pca_nf_ilp <- run_pca_analysis(dds_nf_ilp, "./Deseq2/SRP377060 - Leishmania/results/PCA_NonInfected_vs_InfectedLive.jpeg")
pca_nf_ifp <- run_pca_analysis(dds_nf_ifp, "./Deseq2/SRP377060 - Leishmania/results/PCA_NonInfected_vs_InfectedFixed.jpeg")

# Executa e salva os mapas de distância entre amostras
run_sample_distances(dds_nf_ilp, "./Deseq2/SRP377060 - Leishmania/results/SampleDistances_NonInfected_vs_InfectedLive.jpeg")
run_sample_distances(dds_nf_ifp, "./Deseq2/SRP377060 - Leishmania/results/SampleDistances_NonInfected_vs_InfectedFixed.jpeg")

# Executa e salva os gráficos de Dispersion, Histograma de p-values e MA-plot
run_dispersion_plot(dds_nf_ilp, "./Deseq2/SRP377060 - Leishmania/results/Dispersion_NonInfected_vs_InfectedLive.jpeg")
run_dispersion_plot(dds_nf_ifp, "./Deseq2/SRP377060 - Leishmania/results/Dispersion_NonInfected_vs_InfectedFixed.jpeg")

run_pval_histogram(dds_nf_ilp, "./Deseq2/SRP377060 - Leishmania/results/Histogram_pvalues_NonInfected_vs_InfectedLive.jpeg")
run_pval_histogram(dds_nf_ifp, "./Deseq2/SRP377060 - Leishmania/results/Histogram_pvalues_NonInfected_vs_InfectedFixed.jpeg")

run_ma_plot(dds_nf_ilp, "./Deseq2/SRP377060 - Leishmania/results/MAplot_NonInfected_vs_InfectedLive.jpeg")
run_ma_plot(dds_nf_ifp, "./Deseq2/SRP377060 - Leishmania/results/MAplot_NonInfected_vs_InfectedFixed.jpeg")

