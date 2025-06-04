# install.packages("VennDiagram")
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
library(gridExtra)
library(stringr)
library(readxl)
library(openxlsx)
library(pheatmap)
library(RColorBrewer)
library(ggplot2)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(VennDiagram)


###########################################
# 0. Remoção de amostras indesejadas
###########################################
# Lista de amostras que serão removidas
samples_to_remove <- c("DC_Infected_Live_Parasites_3_SRR19400242", 
                       "DC_Infected_Fixed_Parasites_3_SRR19400239")

###########################################
# 1. Leitura e combinação dos arquivos tabulares
###########################################
tabular_dir  <- "Deseq2/SRP185421 - M tuberculosis/"
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
phenoData <- read_excel("Deseq2/SRP185421 - M tuberculosis/SRP185421.xlsx", col_names = TRUE)
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

run_pca_analysis <- function(dds, output_file) {
  rld <- rlog(dds, blind = TRUE)
  pca_data <- plotPCA(rld, intgroup = "Treatment", returnData = TRUE)
  pca_variance <- attr(pca_data, "percentVar")
  
  write.xlsx(pca_data, file = sub(".jpeg", ".xlsx", output_file))
  
  p <- ggplot(pca_data, aes(PC1, PC2, color = Treatment, label = rownames(pca_data))) +
    geom_point(size = 3) +
    geom_text_repel(size = 3, max.overlaps = 10) +
    xlab(paste0("PC1: ", round(pca_variance[1] * 100, 2), "% variance")) +
    ylab(paste0("PC2: ", round(pca_variance[2] * 100, 2), "% variance")) +
    theme_minimal()
  
  ggsave(output_file, p, width = 6, height = 5, dpi = 300)
  return(p)
}

# Função para detectar o organismo com base nos IDs
detect_organism <- function(entrez_ids) {
  # Verifica se os IDs pertencem ao organismo humano ou camundongo
  hs_ids <- keys(org.Hs.eg.db, keytype = "ENTREZID")
  mm_ids <- keys(org.Mm.eg.db, keytype = "ENTREZID")
  
  if (any(entrez_ids %in% hs_ids)) {
    return("Homo sapiens")
  } else if (any(entrez_ids %in% mm_ids)) {
    return("Mus musculus")
  } else {
    stop("Organismo não reconhecido para os IDs fornecidos.")
  }
}

# Função para mapear IDs Entrez para símbolos de genes, ignorando erros
convert_entrez_to_symbol <- function(entrez_ids) {
  organism <- tryCatch(detect_organism(entrez_ids), error = function(e) NA)
  
  if (is.na(organism)) {
    warning("Nenhum organismo reconhecido. IDs serão ignorados.")
    return(rep(NA, length(entrez_ids)))
  }
  
  gene_symbols <- tryCatch({
    if (organism == "Homo sapiens") {
      mapIds(org.Hs.eg.db, keys = as.character(entrez_ids), 
             column = "SYMBOL", keytype = "ENTREZID", multiVals = "first")
    } else {
      mapIds(org.Mm.eg.db, keys = as.character(entrez_ids), 
             column = "SYMBOL", keytype = "ENTREZID", multiVals = "first")
    }
  }, error = function(e) {
    warning("Erro ao mapear IDs. Alguns genes podem estar ausentes.")
    return(rep(NA, length(entrez_ids)))
  })
  
  return(gene_symbols)
}


# run_volcano_plot <- function(dds, output_file, pval_threshold = 0.05, log2fc_threshold = 1) {
#   dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
#   
#   res <- results(dds)
#   res_df <- na.omit(as.data.frame(res))
#   res_df$GeneID <- rownames(res_df)  # Mantém os IDs Entrez
#   
#   # Converte IDs Entrez para símbolos de genes
#   res_df$Gene <- convert_entrez_to_symbol(res_df$GeneID)
#   
#   res_df$Significance <- "Not Significant"
#   res_df$Significance[res_df$pvalue < pval_threshold & abs(res_df$log2FoldChange) >= log2fc_threshold] <- "Significant"
#   res_df$Color <- ifelse(res_df$Significance == "Significant", ifelse(res_df$log2FoldChange > 0, "red", "blue"), "grey")
#   
#   write.xlsx(res_df, file = sub(".jpeg", ".xlsx", output_file), rowNames = FALSE)
#   
#   top_genes <- res_df[res_df$Significance == "Significant" & res_df$pvalue < 0.01, ]
#   
#   p <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(pvalue), color = Color, label = Gene)) +
#     geom_point(alpha = 0.7, size = 2) +
#     scale_color_manual(values = c("red" = "red", "blue" = "blue", "grey" = "grey")) +
#     geom_hline(yintercept = -log10(pval_threshold), linetype = "dashed", color = "black") +
#     geom_vline(xintercept = c(-log2fc_threshold, log2fc_threshold), linetype = "dashed", color = "black") +
#     xlab("Log2 Fold Change") +
#     ylab("-Log10 P-value") +
#     ggtitle("Volcano Plot") +
#     theme_minimal() +
#     theme(legend.position = "none") +
#     geom_text_repel(data = top_genes, aes(label = Gene), size = 3, max.overlaps = 10)
#   
#   ggsave(output_file, p, width = 8, height = 6, dpi = 300)
#   return(p)
# }
# 
# run_ma_plot <- function(dds, output_file) {
#   dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
#   
#   res <- results(dds)
#   res_df <- as.data.frame(res)
#   res_df$GeneID <- rownames(res_df)  # Mantém os IDs Entrez
#   
#   # Converte IDs Entrez para símbolos de genes
#   res_df$Gene <- convert_entrez_to_symbol(res_df$GeneID)
#   
#   res_df$Significance <- "Not Significant"
#   res_df$Significance[res_df$padj < 0.05] <- "Significant"
#   res_df$Color <- "grey"
#   res_df$Color[res_df$padj < 0.05 & res_df$log2FoldChange > 0] <- "red"
#   res_df$Color[res_df$padj < 0.05 & res_df$log2FoldChange < 0] <- "blue"
#   
#   write.xlsx(res_df, file = sub(".jpeg", ".xlsx", output_file), rowNames = FALSE)
#   
#   top_genes <- res_df[abs(res_df$log2FoldChange) > 2 & res_df$padj < 0.05, ]
#   
#   p <- ggplot(res_df, aes(x = log10(baseMean), y = log2FoldChange, color = Color, label = Gene)) +
#     geom_point(alpha = 0.7, size = 2) +
#     scale_color_manual(values = c("red" = "red", "blue" = "blue", "grey" = "grey")) +
#     geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
#     xlab("Log10 Mean Expression") +
#     ylab("Log2 Fold Change") +
#     ggtitle("MA-plot") +
#     theme_minimal() +
#     theme(legend.position = "none") +
#     geom_text_repel(data = top_genes, aes(label = Gene), size = 3, max.overlaps = 10)
#   
#   ggsave(output_file, p, width = 8, height = 6, dpi = 300)
#   return(p)
# }

# 1. Volcano Plot usando padj (FDR) em vez de pvalue bruto
run_volcano_plot <- function(dds, output_file, padj_threshold = 0.05, log2fc_threshold = 1) {
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  
  res <- results(dds)
  res_df <- na.omit(as.data.frame(res))
  res_df$GeneID <- rownames(res_df)  # Mantém os IDs Entrez
  
  # Converte IDs Entrez para símbolos de genes
  res_df$Gene <- convert_entrez_to_symbol(res_df$GeneID)
  
  # Define significância baseada em FDR e magnitude de FC
  res_df$Significance <- "Not Significant"
  sig_idx <- res_df$padj < padj_threshold & abs(res_df$log2FoldChange) >= log2fc_threshold
  res_df$Significance[sig_idx] <- "Significant"
  
  # Define cores
  res_df$Color <- "grey"
  res_df$Color[res_df$Significance == "Significant" & res_df$log2FoldChange > 0] <- "red"
  res_df$Color[res_df$Significance == "Significant" & res_df$log2FoldChange < 0] <- "blue"
  
  # Salva tabela com resultados e anotações
  write.xlsx(res_df, file = sub("\\.jpeg$", ".xlsx", output_file), rowNames = FALSE)
  
  # Seleciona top genes para rotular (mais rigoroso no FDR)
  top_genes <- subset(res_df, padj < padj_threshold/10 & abs(log2FoldChange) >= log2fc_threshold)
  
  # Gera o gráfico
  p <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = Color, label = Gene)) +
    geom_point(alpha = 0.7, size = 2) +
    scale_color_manual(values = c("red" = "red", "blue" = "blue", "grey" = "grey")) +
    geom_hline(yintercept = -log10(padj_threshold), linetype = "dashed", color = "black") +
    geom_vline(xintercept = c(-log2fc_threshold, log2fc_threshold), linetype = "dashed", color = "black") +
    xlab("Log2 Fold Change") +
    ylab("-Log10 Adjusted P-value (FDR)") +
    ggtitle("Volcano Plot (FDR + log2FC)") +
    theme_minimal() +
    theme(legend.position = "none") +
    geom_text_repel(data = top_genes, aes(label = Gene), size = 3, max.overlaps = 10)
  
  ggsave(output_file, p, width = 8, height = 6, dpi = 300)
  return(p)
}


# 2. MA-plot já utiliza padj, apenas ajustamos thresholds parametrizáveis
run_ma_plot <- function(dds, output_file, padj_threshold = 0.05, log2fc_label = 2) {
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  
  res <- results(dds)
  res_df <- as.data.frame(res)
  res_df$GeneID <- rownames(res_df)
  
  # Converte IDs Entrez para símbolos de genes
  res_df$Gene <- convert_entrez_to_symbol(res_df$GeneID)
  
  # Define significância via FDR
  res_df$Significance <- ifelse(res_df$padj < padj_threshold, "Significant", "Not Significant")
  
  # Define cores baseadas em direção da FC
  res_df$Color <- "grey"
  res_df$Color[res_df$padj < padj_threshold & res_df$log2FoldChange > 0] <- "red"
  res_df$Color[res_df$padj < padj_threshold & res_df$log2FoldChange < 0] <- "blue"
  
  # Salva tabela
  write.xlsx(res_df, file = sub("\\.jpeg$", ".xlsx", output_file), rowNames = FALSE)
  
  # Genes para rotular: magnitude alta e FDR significativo
  top_genes <- subset(res_df, padj < padj_threshold & abs(log2FoldChange) >= log2fc_label)
  
  # Gera o MA-plot
  p <- ggplot(res_df, aes(x = log10(baseMean), y = log2FoldChange, color = Color, label = Gene)) +
    geom_point(alpha = 0.7, size = 2) +
    scale_color_manual(values = c("red" = "red", "blue" = "blue", "grey" = "grey")) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    xlab("Log10 Mean Expression") +
    ylab("Log2 Fold Change") +
    ggtitle("MA-plot (FDR + log2FC)") +
    theme_minimal() +
    theme(legend.position = "none") +
    geom_text_repel(data = top_genes, aes(label = Gene), size = 3, max.overlaps = 10)
  
  ggsave(output_file, p, width = 8, height = 6, dpi = 300)
  return(p)
}

phenoData
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


# Função para criar os diagramas de Venn e salvar os resultados com a tabela contendo IDs e nomes de genes
run_venn_analysis <- function(dds1, dds2, comparison1, comparison2, output_prefix) {
  # Extrai os resultados e converte para data.frame
  res1 <- as.data.frame(results(dds1))
  res2 <- as.data.frame(results(dds2))
  
  # Filtra genes upregulados e downregulados, removendo os casos com NA
  up1 <- rownames(res1)[!is.na(res1$log2FoldChange) & !is.na(res1$padj) & (res1$log2FoldChange > 1) & (res1$padj < 0.05)]
  up2 <- rownames(res2)[!is.na(res2$log2FoldChange) & !is.na(res2$padj) & (res2$log2FoldChange > 1) & (res2$padj < 0.05)]
  
  down1 <- rownames(res1)[!is.na(res1$log2FoldChange) & !is.na(res1$padj) & (res1$log2FoldChange < -1) & (res1$padj < 0.05)]
  down2 <- rownames(res2)[!is.na(res2$log2FoldChange) & !is.na(res2$padj) & (res2$log2FoldChange < -1) & (res2$padj < 0.05)]
  
  # Interseção e exclusividade
  up_intersect <- intersect(up1, up2)
  up_only1 <- setdiff(up1, up2)
  up_only2 <- setdiff(up2, up1)
  
  down_intersect <- intersect(down1, down2)
  down_only1 <- setdiff(down1, down2)
  down_only2 <- setdiff(down2, down1)
  
  # Ajuste dos nomes nos diagramas de Venn
  venn_up <- venn.diagram(
    x = list(Comp1 = up1, Comp2 = up2),
    category.names = c(comparison1, comparison2),
    filename = paste0(output_prefix, "_Venn_Up.jpeg"),
    output = TRUE,
    fill = c("red", "blue"),
    alpha = 0.5,
    cex = 1.5,
    cat.cex = 0.8,     # Diminui o tamanho do texto dos círculos
    cat.pos = c(-30, 30),
    cat.dist = c(0.05, 0.05)
  )
  
  venn_down <- venn.diagram(
    x = list(Comp1 = down1, Comp2 = down2),
    category.names = c(comparison1, comparison2),
    filename = paste0(output_prefix, "_Venn_Down.jpeg"),
    output = TRUE,
    fill = c("red", "blue"),
    alpha = 0.5,
    cex = 1.5,
    cat.cex = 0.8,
    cat.pos = c(-30, 30),
    cat.dist = c(0.05, 0.05)
  )
  
  # Função auxiliar para gerar um data.frame com Entrez ID e Gene_Name
  get_gene_info <- function(gene_list) {
    data.frame(
      Entrez_ID = gene_list,
      Gene_Name = convert_entrez_to_symbol(gene_list)
    )
  }
  
  # Criar planilha Excel com as listas
  wb <- createWorkbook()
  
  addWorksheet(wb, "Upregulated_Intersect")
  addWorksheet(wb, "Upregulated_Only_Comp1")
  addWorksheet(wb, "Upregulated_Only_Comp2")
  addWorksheet(wb, "Downregulated_Intersect")
  addWorksheet(wb, "Downregulated_Only_Comp1")
  addWorksheet(wb, "Downregulated_Only_Comp2")
  
  writeData(wb, "Upregulated_Intersect", get_gene_info(up_intersect))
  writeData(wb, "Upregulated_Only_Comp1", get_gene_info(up_only1))
  writeData(wb, "Upregulated_Only_Comp2", get_gene_info(up_only2))
  writeData(wb, "Downregulated_Intersect", get_gene_info(down_intersect))
  writeData(wb, "Downregulated_Only_Comp1", get_gene_info(down_only1))
  writeData(wb, "Downregulated_Only_Comp2", get_gene_info(down_only2))
  
  saveWorkbook(wb, paste0(output_prefix, "_VennData.xlsx"), overwrite = TRUE)
  
  message("Diagramas de Venn e arquivos Excel gerados com sucesso!")
}


run_deseq_up_down_analysis <- function(dds, up_threshold, down_threshold, padj_cutoff, output_file) {
  # Executa a análise diferencial
  res <- results(dds)
  
  # Filtra genes upregulados
  up_genes <- subset(res, log2FoldChange > up_threshold & padj < padj_cutoff)
  up_genes <- up_genes[!is.na(up_genes$padj), ]
  
  # Filtra genes downregulados
  down_genes <- subset(res, log2FoldChange < down_threshold & padj < padj_cutoff)
  down_genes <- down_genes[!is.na(down_genes$padj), ]
  
  # Converte rownames (genes) para uma coluna separada
  up_genes$Entrez_ID <- rownames(up_genes)
  down_genes$Entrez_ID <- rownames(down_genes)
  
  # Salva os resultados em Excel
  wb <- createWorkbook()
  addWorksheet(wb, "Upregulated")
  addWorksheet(wb, "Downregulated")
  writeData(wb, "Upregulated", up_genes)
  writeData(wb, "Downregulated", down_genes)
  saveWorkbook(wb, output_file, overwrite = TRUE)
  
  return(list(up = up_genes, down = down_genes))
}

###########################################
# 5. Execução das análises
###########################################


# Antes de criar o DESeqDataSet, force Treatment como fator:
phenoData$Treatment <- factor(phenoData$Treatment, levels = c("Uninfected", "Infected"))

# Agora execute a análise DESeq2:
dds_nf_ilp <- run_deseq_analysis(
  data, 
  phenoData, 
  "Uninfected", 
  "Infected",
  "./Deseq2/SRP185421 - M tuberculosis/results/DESeq2_Uninfected_vs_Infected.tabular"
)





intersect(phenoData$SampleName[phenoData$Treatment %in% c("Uninfected", "Infected")], colnames(data))
colnames(data)
phenoData


# # Executa DESeq2 para as comparações desejadas
# dds_nf_ilp <- run_deseq_analysis(data, phenoData, "Uninfected", "Infected",
#                                  "./Deseq2/SRP185421 - M tuberculosis/results/DESeq2_Uninfected_vs_Infected.tabular")
# # 
# dds_nf_ifp <- run_deseq_analysis(data, phenoData, "Non-Infected", "Infected_Fixed_Parasites",
#                                  "./Deseq2/SRP377060 - Leishmania/results/DESeq2_NonInfected_vs_InfectedFixed.tabular")

# Visualiza as contagens (primeiras linhas) dos objetos DESeqDataSet
head(assay(dds_nf_ilp))
# head(assay(dds_nf_ifp))

# Executa e salva os gráficos de PCA
pca_nf_ilp <- run_pca_analysis(dds_nf_ilp, "./Deseq2/SRP185421 - M tuberculosis/results/PCA_Uninfected_vs_Infected.jpeg")
# pca_nf_ifp <- run_pca_analysis(dds_nf_ifp, "./Deseq2/SRP377060 - Leishmania/results/PCA_NonInfected_vs_InfectedFixed.jpeg")

# Executa e salva os mapas de distância entre amostras
run_sample_distances(dds_nf_ilp, "./Deseq2/SRP185421 - M tuberculosis/results/SampleDistances_Uninfected_vs_Infected.jpeg")
# run_sample_distances(dds_nf_ifp, "./Deseq2/SRP377060 - Leishmania/results/SampleDistances_NonInfected_vs_InfectedFixed.jpeg")

# Executa e salva os gráficos de Dispersion, Histograma de p-values e MA-plot
run_dispersion_plot(dds_nf_ilp, "./Deseq2/SRP185421 - M tuberculosis/results/Dispersion_Uninfected_vs_Infected.jpeg")
# run_dispersion_plot(dds_nf_ifp, "./Deseq2/SRP377060 - Leishmania/results/Dispersion_NonInfected_vs_InfectedFixed.jpeg")

run_pval_histogram(dds_nf_ilp, "./Deseq2/SRP185421 - M tuberculosis/results/Histogram_pvalues_Uninfected_vs_Infected.jpeg")
# run_pval_histogram(dds_nf_ifp, "./Deseq2/SRP377060 - Leishmania/results/Histogram_pvalues_NonInfected_vs_InfectedFixed.jpeg")

ma1 = run_ma_plot(dds_nf_ilp, "./Deseq2/SRP185421 - M tuberculosis/results/MAplot_Uninfected_vs_Infected.jpeg")
# ma2 = run_ma_plot(dds_nf_ifp, "./Deseq2/SRP377060 - Leishmania/results/MAplot_NonInfected_vs_InfectedFixed.jpeg")
head(ma1$data)

# Executa e salva os Volcano Plots
run_volcano_plot(dds_nf_ilp, "./Deseq2/SRP185421 - M tuberculosis/results/Volcano_Uninfected_vs_Infected.jpeg")
# run_volcano_plot(dds_nf_ifp, "./Deseq2/SRP377060 - Leishmania/results/Volcano_NonInfected_vs_InfectedFixed.jpeg")

# a= run_volcano_plot(dds_nf_ifp, "./Deseq2/SRP377060 - Leishmania/results/Volcano_NonInfected_vs_InfectedFixed.jpeg")
# head(a$data)



# run_venn_analysis(dds_nf_ilp, dds_nf_ifp, 
#                   "NonInfected_vs_InfectedLive", "NonInfected_vs_InfectedFixed",
#                   "./Deseq2/SRP377060 - Leishmania/results/Venn_Analysis")


############################

# library(openxlsx)
# library(clusterProfiler)
# library(org.Hs.eg.db)
# library(org.Mm.eg.db)
# library(ggplot2)
# library(enrichplot)
# 
# # Função para detectar o organismo com base nos IDs
# detect_organism <- function(entrez_ids) {
#   hs_ids <- keys(org.Hs.eg.db, keytype = "ENTREZID")
#   mm_ids <- keys(org.Mm.eg.db, keytype = "ENTREZID")
#   
#   if (any(entrez_ids %in% hs_ids)) {
#     return("Homo sapiens")
#   } else if (any(entrez_ids %in% mm_ids)) {
#     return("Mus musculus")
#   } else {
#     stop("Organismo não reconhecido para os IDs fornecidos.")
#   }
# }
# 
# # Define o arquivo de entrada (o Excel com as abas de Venn)
# input_file <- "./Deseq2/SRP377060 - Leishmania/results/Venn_Analysis_VennData.xlsx"
# 
# # Obter os nomes das abas no arquivo Excel
# sheet_names <- getSheetNames(input_file)

# Loop por cada aba para realizar enriquecimento
# for(sheet in sheet_names) {
#   # Ler a aba
#   gene_df <- read.xlsx(input_file, sheet = sheet)
#   
#   # Verifica se a coluna "Entrez_ID" existe
#   if(!("Entrez_ID" %in% colnames(gene_df))) {
#     stop("A coluna 'Entrez_ID' não foi encontrada na aba ", sheet)
#   }
#   
#   # Pega a lista de IDs, removendo NAs ou strings vazias
#   gene_list <- as.character(gene_df$Entrez_ID)
#   gene_list <- gene_list[!is.na(gene_list) & gene_list != ""]
#   
#   # Detecta o organismo usando a função detect_organism
#   organism <- tryCatch(detect_organism(gene_list), error = function(e) NA)
#   if(is.na(organism)) {
#     message("Nenhum organismo reconhecido para a aba ", sheet, ". Pulando enriquecimento.")
#     next
#   }
#   
#   # Seleciona o OrgDb com base no organismo detectado
#   orgDb <- if(organism == "Homo sapiens") org.Hs.eg.db else org.Mm.eg.db
#   
#   # Executa enriquecimento GO para Biological Process (BP)
#   ego <- enrichGO(gene         = gene_list,
#                   OrgDb        = orgDb,
#                   keyType      = "ENTREZID",
#                   ont          = "BP",
#                   pAdjustMethod = "BH",
#                   pvalueCutoff  = 0.05,
#                   qvalueCutoff  = 0.2,
#                   readable     = TRUE)
#   
#   # Se não houver genes mapeados, pula a aba
#   if(is.null(ego) || nrow(as.data.frame(ego)) == 0) {
#     message("Nenhum gene mapeado para enriquecimento na aba ", sheet, ". Pulando.")
#     next
#   }
#   
#   # Salva a tabela de enriquecimento em um arquivo Excel
#   output_excel <- paste0("./Deseq2/SRP377060 - Leishmania/results/Enrichment_", sheet, ".xlsx")
#   write.xlsx(as.data.frame(ego), file = output_excel, overwrite = TRUE)
#   
#   # Gera e salva um dotplot (como PDF) para visualização
#   output_pdf <- paste0("./Deseq2/SRP377060 - Leishmania/results/Enrichment_", sheet, ".pdf")
#   pdf(output_pdf, width = 10, height = 8)
#   print(dotplot(ego, showCategory = 20) + ggtitle(paste("GO Enrichment for", sheet)))
#   dev.off()
#   
#   message("Enrichment analysis for sheet '", sheet, "' completed!")
# }

######################################



library(DESeq2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(openxlsx)
library(ggplot2)
library(enrichplot)

# Função para detectar o organismo com base nos IDs
detect_organism <- function(entrez_ids) {
  hs_ids <- keys(org.Hs.eg.db, keytype = "ENTREZID")
  mm_ids <- keys(org.Mm.eg.db, keytype = "ENTREZID")
  
  if (any(entrez_ids %in% hs_ids)) {
    return("Homo sapiens")
  } else if (any(entrez_ids %in% mm_ids)) {
    return("Mus musculus")
  } else {
    stop("Organismo não reconhecido para os IDs fornecidos.")
  }
}

# Função para selecionar o banco de dados correto
select_orgDb <- function(gene_ids) {
  organism <- detect_organism(gene_ids)
  if (organism == "Homo sapiens") {
    return(org.Hs.eg.db)
  } else {
    return(org.Mm.eg.db)
  }
}

# Função principal para análise de enriquecimento
run_deseq_up_down_enrichment <- function(dds, 
                                         up_threshold = 1, 
                                         down_threshold = -1, 
                                         padj_cutoff = 0.05,
                                         ont = "BP",
                                         pAdjustMethod = "BH",
                                         pvalueCutoff = 0.05,
                                         qvalueCutoff = 0.2,
                                         output_file_prefix) {
  
  # 1. Extrai DEGs e salva os resultados em Excel
  deg_results <- run_deseq_up_down_analysis(dds, 
                                            up_threshold = up_threshold, 
                                            down_threshold = down_threshold, 
                                            padj_cutoff = padj_cutoff,
                                            output_file = paste0(output_file_prefix, "_DEGs.xlsx"))
  
  # Obtém os vetores de IDs (como caracteres) para up e down
  up_ids <- as.character(deg_results$up$Entrez_ID)
  down_ids <- as.character(deg_results$down$Entrez_ID)
  
  # Remove possíveis strings vazias
  up_ids <- up_ids[up_ids != ""]
  down_ids <- down_ids[down_ids != ""]
  
  # Determina o OrgDb a partir dos IDs
  all_ids <- unique(c(up_ids, down_ids))
  orgDb <- select_orgDb(all_ids)
  
  # 2. Realiza enriquecimento GO para os genes upregulados e downregulados
  ego_up <- enrichGO(gene = up_ids,
                     OrgDb = orgDb,
                     keyType = "ENTREZID",
                     ont = ont,
                     pAdjustMethod = pAdjustMethod,
                     pvalueCutoff = pvalueCutoff,
                     qvalueCutoff = qvalueCutoff,
                     readable = TRUE)
  
  ego_down <- enrichGO(gene = down_ids,
                       OrgDb = orgDb,
                       keyType = "ENTREZID",
                       ont = ont,
                       pAdjustMethod = pAdjustMethod,
                       pvalueCutoff = pvalueCutoff,
                       qvalueCutoff = qvalueCutoff,
                       readable = TRUE)
  
  # 3. Salva os resultados de enriquecimento em uma planilha Excel com duas abas
  wb <- createWorkbook()
  addWorksheet(wb, "Up_Enrichment")
  addWorksheet(wb, "Down_Enrichment")
  writeData(wb, "Up_Enrichment", as.data.frame(ego_up))
  writeData(wb, "Down_Enrichment", as.data.frame(ego_down))
  saveWorkbook(wb, paste0(output_file_prefix, "_Enrichment.xlsx"), overwrite = TRUE)
  
  # 4. Gera e salva os dotplots de enriquecimento (PDF)
  output_pdf_up <- paste0(output_file_prefix, "_Enrichment_Up.pdf")
  output_pdf_down <- paste0(output_file_prefix, "_Enrichment_Down.pdf")
  
  pdf(output_pdf_up, width = 10, height = 8)
  print(dotplot(ego_up, showCategory = 20) + ggtitle("GO Enrichment for Upregulated Genes"))
  dev.off()
  
  pdf(output_pdf_down, width = 10, height = 8)
  print(dotplot(ego_down, showCategory = 20) + ggtitle("GO Enrichment for Downregulated Genes"))
  dev.off()
  
  message("Enrichment analysis for DEGs completed! Files saved with prefix: ", output_file_prefix)
  
  return(list(ego_up = ego_up, ego_down = ego_down))
}

# Exemplo de uso:


run_deseq_up_down_enrichment <- function(dds,
                                         up_threshold   = 1,
                                         down_threshold = -1,
                                         padj_cutoff    = 0.05,
                                         ontologies     = c("BP","MF"),
                                         pAdjustMethod  = "BH",
                                         pvalueCutoff   = 0.05,
                                         qvalueCutoff   = 0.2,
                                         output_prefix) {
  
  # (1) Extrai DEGs e mantém em memória
  degs <- run_deseq_up_down_analysis(dds,
                                     up_threshold   = up_threshold,
                                     down_threshold = down_threshold,
                                     padj_cutoff    = padj_cutoff,
                                     output_file    = paste0(output_prefix, "_DEGs.xlsx"))
  up_ids   <- na.omit(as.character(degs$up$Entrez_ID))
  down_ids <- na.omit(as.character(degs$down$Entrez_ID))
  
  # Determina o OrgDb
  orgDb <- select_orgDb(unique(c(up_ids, down_ids)))
  
  # Resultado em lista
  enrich_results <- list()
  
  # Inicia workbook para salvar todos os resultados
  wb <- createWorkbook()
  
  for(ont in ontologies) {
    # up-regulated
    ego_up <- enrichGO(gene          = up_ids,
                       OrgDb         = orgDb,
                       keyType       = "ENTREZID",
                       ont           = ont,
                       pAdjustMethod = pAdjustMethod,
                       pvalueCutoff  = pvalueCutoff,
                       qvalueCutoff  = qvalueCutoff,
                       readable      = TRUE)
    enrich_results[[paste0("up_", ont)]] <- ego_up
    
    sheet_up <- paste0("Up_", ont)
    addWorksheet(wb, sheet_up)
    writeData(wb, sheet_up, as.data.frame(ego_up))
    pdf(paste0(output_prefix, "_", sheet_up, ".pdf"), width = 10, height = 8)
    print(dotplot(ego_up, showCategory = 20) + ggtitle(paste("GO", ont, "Enrichment for Upregulated Genes")))
    dev.off()
    
    # down-regulated
    ego_down <- enrichGO(gene          = down_ids,
                         OrgDb         = orgDb,
                         keyType       = "ENTREZID",
                         ont           = ont,
                         pAdjustMethod = pAdjustMethod,
                         pvalueCutoff  = pvalueCutoff,
                         qvalueCutoff  = qvalueCutoff,
                         readable      = TRUE)
    enrich_results[[paste0("down_", ont)]] <- ego_down
    
    sheet_down <- paste0("Down_", ont)
    addWorksheet(wb, sheet_down)
    writeData(wb, sheet_down, as.data.frame(ego_down))
    pdf(paste0(output_prefix, "_", sheet_down, ".pdf"), width = 10, height = 8)
    print(dotplot(ego_down, showCategory = 20) + ggtitle(paste("GO", ont, "Enrichment for Downregulated Genes")))
    dev.off()
  }
  
  # Salva workbook com todas as abas
  saveWorkbook(wb, paste0(output_prefix, "_Enrichment_GO.xlsx"), overwrite = TRUE)
  message("GO enrichment (BP & MF) completed! Files saved with prefix: ", output_prefix)
  
  # Retorna lista com objetos enrichResult
  return(enrich_results)
}



# # Supondo que 'dds_nf_ilp' seja o objeto DESeq2 resultante da comparação desejada.
# up_down_enrich_results <- run_deseq_up_down_enrichment(dds = dds_nf_ilp,
#                                                        up_threshold = 1,
#                                                        down_threshold = -1,
#                                                        padj_cutoff = 0.05,
#                                                        ont = "BP",
#                                                        pAdjustMethod = "BH",
#                                                        pvalueCutoff = 0.05,
#                                                        qvalueCutoff = 0.2,
#                                                        output_file_prefix = "./Deseq2/SRP185421 - M tuberculosis/results/DESeq2_UpDown")

up_down_enrich_results <- run_deseq_up_down_enrichment(
  dds             = dds_nf_ilp,
  up_threshold    = 1,
  down_threshold  = -1,
  padj_cutoff     = 0.05,
  ontologies      = c("BP","MF"),
  pAdjustMethod   = "BH",
  pvalueCutoff    = 0.05,
  qvalueCutoff    = 0.2,
  output_prefix   = "./Deseq2/SRP185421 - M tuberculosis/results/DESeq2_UpDown"
)

up_down_enrich_results


#################
# 
# 
# run_deseq_up_down_analysis <- function(dds, up_threshold, down_threshold, padj_cutoff, output_file) {
#   library(AnnotationDbi)
#   
#   # Executa a análise diferencial
#   res <- results(dds)
#   
#   # Filtra genes upregulados
#   up_genes <- subset(res, log2FoldChange > up_threshold & padj < padj_cutoff)
#   up_genes <- up_genes[!is.na(up_genes$padj), ]
#   
#   # Filtra genes downregulados
#   down_genes <- subset(res, log2FoldChange < down_threshold & padj < padj_cutoff)
#   down_genes <- down_genes[!is.na(down_genes$padj), ]
#   
#   # Converte rownames (genes) para uma coluna separada
#   up_genes$Entrez_ID <- rownames(up_genes)
#   down_genes$Entrez_ID <- rownames(down_genes)
#   
#   # Seleciona o banco de dados correto com base nos IDs Entrez
#   orgDb <- select_orgDb(up_genes$Entrez_ID)
#   
#   # Converte Entrez ID para Symbol
#   up_gene_symbols <- mapIds(orgDb, keys = up_genes$Entrez_ID, 
#                             column = "SYMBOL", keytype = "ENTREZID", multiVals = "first")
#   down_gene_symbols <- mapIds(orgDb, keys = down_genes$Entrez_ID, 
#                               column = "SYMBOL", keytype = "ENTREZID", multiVals = "first")
#   
#   # Adiciona os símbolos aos data frames
#   up_genes$Gene_Symbol <- up_gene_symbols
#   down_genes$Gene_Symbol <- down_gene_symbols
#   
#   # Salva os resultados em Excel
#   wb <- createWorkbook()
#   addWorksheet(wb, "Upregulated")
#   addWorksheet(wb, "Downregulated")
#   writeData(wb, "Upregulated", up_genes)
#   writeData(wb, "Downregulated", down_genes)
#   saveWorkbook(wb, output_file, overwrite = TRUE)
#   
#   return(list(up = up_genes, down = down_genes))
# }
# 
# up_down_enrich_results <- run_deseq_up_down_enrichment(dds = dds_nf_ilp,
#                                                        up_threshold = 1,
#                                                        down_threshold = -1,
#                                                        padj_cutoff = 0.05,
#                                                        ont = "BP",
#                                                        pAdjustMethod = "BH",
#                                                        pvalueCutoff = 0.05,
#                                                        qvalueCutoff = 0.2,
#                                                        output_file_prefix = "./Deseq2/SRP185421 - M tuberculosis/results/DESeq2_UpDown")
# 


# head(up_down_enrich_results)
# 
# 
# keywords <- c("acetyltransferase", "acetylation", "desacetylation",
#               "histone lysine","histone", 
#               "methyltransferase",
#               "deacetylase",
#               "histone deacetylase",
#               "lysine acetyltransferase",
#               "lysine")
# 
# # monta uma regex que busca qualquer um deles (case insensitive)
# pattern  <- paste(keywords, collapse = "|")
# # 
# res_list <- up_down_enrich_results
# # 
# library(dplyr)
# filtered <- lapply(names(res_list), function(nm) {
#   df <- as.data.frame(res_list[[nm]])
#   df %>%
#     filter(grepl(pattern, Description, ignore.case = TRUE)) %>%
#     mutate(category = nm)
# })
# filtered_df <- bind_rows(filtered)
# filtered_df
# 
# 
# 
# # library(openxlsx)
# #
# wb <- createWorkbook()
# for(nm in unique(filtered_df$category)) {
#   sheet_df <- filtered_df %>% filter(category == nm)
#   addWorksheet(wb, nm)
#   writeData(wb, nm, sheet_df)
# }
# saveWorkbook(wb, "./Deseq2/SRP185421 - M tuberculosis/results/filtered_enrichment_keywords.xlsx", overwrite = TRUE)
# 
# 
# write.csv(filtered_df,
#           "./Deseq2/SRP185421 - M tuberculosis/results/filtered_enrichment_keywords_all.csv",
#           row.names = FALSE)
# 
# 
# 
# 
# library(dplyr)
# library(openxlsx)
# 
# # 1) lista de genes de interesse
# targets <- c("NAT10","HAT1",
#              "KAT2A","KAT2B","KAT5","KAT6A","KAT6B","KAT7","KAT8","KAT12",
#              "GTF3C4","CREBBP","aTAT1","p300","HDAC1","HDAC2","HDAC3",
#              "HDAC4","HDAC5","HDAC6","HDAC7","HDAC8","HDAC9","HDAC10",
#              "SIRT1","SIRT2","SIRT3","SIRT4","SIRT5","SIRT6","SIRT7")
# 
# # monta regex: palavra exata, ignorando case
# pattern_genes <- paste0("\\b(", paste(targets, collapse="|"), ")\\b")
# 
# # 2) percorre cada categoria de enrichment, filtra pelo geneID
# filtered_by_gene <- lapply(names(up_down_enrich_results), function(cat) {
#   df <- as.data.frame(up_down_enrich_results[[cat]])
#   df %>%
#     # geneID é uma string do tipo "GeneA/GeneB/…"
#     filter(grepl(pattern_genes, geneID, ignore.case = TRUE)) %>%
#     mutate(category = cat)
# })
# 
# # une em um só
# filtered_genes_df <- bind_rows(filtered_by_gene)
# 
# # 3) salva em Excel, uma aba por categoria
# wb2 <- createWorkbook()
# for(cat in unique(filtered_genes_df$category)) {
#   sheet_df <- filtered_genes_df %>% filter(category == cat)
#   addWorksheet(wb2, cat)
#   writeData(wb2, cat, sheet_df)
# }
# saveWorkbook(wb2, "./Deseq2/SRP185421 - M tuberculosis/results/filtered_by_genes_enrichment_targets.xlsx", overwrite = TRUE)
# 
# # 4) opcional: salva um CSV único
# write.csv(filtered_genes_df,
#           "./Deseq2/SRP185421 - M tuberculosis/results/filtered_by_genes_enrichment_all.csv",
#           row.names = FALSE)
# 
# # Em memória você fica com:
# #   filtered_genes_df        # todas as linhas que contêm algum gene alvo
# #   filtered_by_gene$up_BP   # etc., por categoria
# 
# 
# ######################################################################################################################
# 
# 
# 
# library(dplyr)
# library(stringr)
# 
# # 1) sua lista de genes-alvo em uppercase
# targets_up <- toupper(targets)
# 
# # 2) função que extrai apenas os targets, retornando NA se não houver
# extract_targets <- function(geneID_string) {
#   genes <- str_split(geneID_string, "/", simplify = TRUE)
#   found <- genes[ toupper(genes) %in% targets_up ]
#   if (length(found) == 0) return(NA_character_)
#   paste(found, collapse = "/")
# }
# 
# # 3) função que reordena colocando targets na frente
# reorder_with_targets_first <- function(geneID_string) {
#   genes <- str_split(geneID_string, "/", simplify = TRUE)
#   is_tgt <- toupper(genes) %in% targets_up
#   # targets primeiro, depois os demais na ordem original
#   new_order <- c(genes[is_tgt], genes[!is_tgt])
#   paste(new_order, collapse = "/")
# }
# 
# # 4) aplica ao seu data.frame
# df2 <- filtered_genes_df %>%
#   mutate(
#     only_targets   = sapply(geneID, extract_targets),
#     geneID_reorder = sapply(geneID, reorder_with_targets_first)
#   )
# 
# # Exemplo de uso:
# df2 %>% select(category, geneID, only_targets, geneID_reorder) %>% head()
# 
# write.xlsx(df2, "./Deseq2/SRP185421 - M tuberculosis/results/with_targets_and_reordered.xlsx", overwrite = TRUE)
# 

###################################################
###################################################


# # --- 1. Pacotes e Paths ------------------------------------------------------
# library(dplyr)
# library(stringr)
# library(purrr)
# library(openxlsx)
# 
# out_dir       <- "./Deseq2/SRP185421 - M tuberculosis/results/Enrichment_Analysis/"
# terms_file    <- file.path(out_dir, "enrichment_terms.xlsx")
# targets_file  <- file.path(out_dir, "enrichment_targets.xlsx")
# combined_file <- file.path(out_dir, "enrichment_combined.xlsx")
# dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
# 
# # --- 2. Definições -----------------------------------------------------------
# targets   <- c(
#   "NAT10","HAT1","KAT2A","KAT2B","KAT5","KAT6A","KAT6B",
#   "KAT7","KAT8","KAT12","GTF3C4","CREBBP","aTAT1","p300",
#   "HDAC1","HDAC2","HDAC3","HDAC4","HDAC5","HDAC6","HDAC7",
#   "HDAC8","HDAC9","HDAC10","SIRT1","SIRT2","SIRT3","SIRT4",
#   "SIRT5","SIRT6","SIRT7"
# )
# keywords  <- c(
#   "acetyltransferase","acetylation","desacetylation",
#   "histone lysine","histone","methyltransferase",
#   "deacetylase","histone deacetylase",
#   "lysine acetyltransferase","lysine"
# )
# 
# pat_genes    <- paste0("\\b(", paste(targets, collapse="|"), ")\\b")
# pat_keywords <- paste(keywords, collapse="|")
# targets_up   <- toupper(targets)
# 
# # --- 3. Funções Auxiliares ---------------------------------------------------
# 
# # Extrai apenas os genes-alvo de um geneID “A/B/C…”
# extract_targets <- function(geneID) {
#   parts <- str_split(geneID, "/", simplify = TRUE)
#   found <- parts[toupper(parts) %in% targets_up]
#   if (length(found) == 0) NA_character_ else paste(found, collapse = "/")
# }
# 
# # Reordena geneID colocando os alvos na frente
# reorder_genes <- function(geneID) {
#   parts <- str_split(geneID, "/", simplify = TRUE)
#   tgt    <- parts[toupper(parts) %in% targets_up]
#   other  <- parts[!toupper(parts) %in% targets_up]
#   paste(c(tgt, other), collapse = "/")
# }
# 
# # Filtra e anexa colunas comuns
# filter_and_augment <- function(results_list, pattern, on_desc = TRUE) {
#   map_dfr(names(results_list), function(cat) {
#     df <- as.data.frame(results_list[[cat]])
#     df$original_row <- seq_len(nrow(df))
#     df %>%
#       mutate(category = cat) %>%
#       filter(
#         if (on_desc) str_detect(Description, regex(pattern, ignore_case = TRUE))
#         else         str_detect(geneID,      regex(pattern, ignore_case = TRUE))
#       ) %>%
#       mutate(
#         only_targets   = vapply(geneID, extract_targets, FUN.VALUE = character(1)),
#         geneID_reorder = vapply(geneID, reorder_genes,   FUN.VALUE = character(1))
#       )
#   })
# }
# 
# # Escreve abas de um data.frame dividido por “category”
# write_by_category <- function(df, path) {
#   wb <- createWorkbook()
#   df %>% split(.$category) %>%
#     iwalk(~{
#       addWorksheet(wb, .y)
#       writeData(wb, .y, .x)
#     })
#   saveWorkbook(wb, path, overwrite = TRUE)
# }
# 
# # --- 4. Geração das três planilhas ------------------------------------------
# 
# # (A) Enriquecimento apenas por termos
# df_terms   <- filter_and_augment(up_down_enrich_results, pat_keywords, TRUE)
# write_by_category(df_terms,       terms_file)
# 
# # (B) Enriquecimento apenas por genes-alvo
# df_targets <- filter_and_augment(up_down_enrich_results, pat_genes,    FALSE)
# write_by_category(df_targets,     targets_file)
# 
# # (C) União (termos ∪ targets), sem duplicar entradas por category+ID
# df_combined <- bind_rows(df_terms, df_targets) %>%
#   distinct(category, ID, .keep_all = TRUE)
# write_by_category(df_combined,    combined_file)

# --- 1. Pacotes e Paths ------------------------------------------------------
# library(dplyr)
# library(stringr)
# library(purrr)
# library(openxlsx)
# 
# out_dir       <- "./Deseq2/SRP377060 - Leishmania/results/Enrichment_Analysis/"
# terms_file    <- file.path(out_dir, "enrichment_terms.xlsx")
# targets_file  <- file.path(out_dir, "enrichment_targets.xlsx")
# combined_file <- file.path(out_dir, "enrichment_combined.xlsx")
# dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
# 
# # --- 2. Definições -----------------------------------------------------------
# targets   <- c(
#   "NAT10","HAT1","KAT2A","KAT2B","KAT5","KAT6A","KAT6B",
#   "KAT7","KAT8","KAT12","GTF3C4","CREBBP","aTAT1","p300",
#   "HDAC1","HDAC2","HDAC3","HDAC4","HDAC5","HDAC6","HDAC7",
#   "HDAC8","HDAC9","HDAC10","SIRT1","SIRT2","SIRT3","SIRT4",
#   "SIRT5","SIRT6","SIRT7"
# )
# keywords  <- c(
#   "acetyltransferase","acetylation","desacetylation",
#   "histone lysine","histone","methyltransferase",
#   "deacetylase","histone deacetylase",
#   "lysine acetyltransferase","lysine"
# )
# 
# pat_genes    <- paste0("\\b(", paste(targets, collapse="|"), ")\\b")
# pat_keywords <- paste(keywords, collapse="|")
# targets_up   <- toupper(targets)
# 
# # --- 3. Funções Auxiliares ---------------------------------------------------
# 
# # Extrai apenas os genes-alvo de um geneID “A/B/C…”
# extract_targets <- function(geneID) {
#   parts <- str_split(geneID, "/", simplify = TRUE)
#   found <- parts[toupper(parts) %in% targets_up]
#   if (length(found) == 0) NA_character_ else paste(found, collapse = "/")
# }
# 
# # Reordena geneID colocando os alvos na frente
# reorder_genes <- function(geneID) {
#   parts <- str_split(geneID, "/", simplify = TRUE)
#   tgt    <- parts[toupper(parts) %in% targets_up]
#   other  <- parts[!toupper(parts) %in% targets_up]
#   paste(c(tgt, other), collapse = "/")
# }
# 
# # Filtra e anexa colunas comuns
# filter_and_augment <- function(results_list, pattern, on_desc = TRUE) {
#   map_dfr(names(results_list), function(cat) {
#     df <- as.data.frame(results_list[[cat]])
#     df$original_row <- seq_len(nrow(df))
#     df %>%
#       mutate(category = cat) %>%
#       filter(
#         if (on_desc) str_detect(Description, regex(pattern, ignore_case = TRUE))
#         else         str_detect(geneID,      regex(pattern, ignore_case = TRUE))
#       ) %>%
#       mutate(
#         only_targets   = vapply(geneID, extract_targets, FUN.VALUE = character(1)),
#         geneID_reorder = vapply(geneID, reorder_genes,   FUN.VALUE = character(1))
#       )
#   })
# }
# 
# # Escreve abas de um data.frame dividido por “category”
# write_by_category <- function(df, path) {
#   wb <- createWorkbook()
#   df %>% split(.$category) %>%
#     iwalk(~{
#       addWorksheet(wb, .y)
#       writeData(wb, .y, .x)
#     })
#   saveWorkbook(wb, path, overwrite = TRUE)
# }
# 
# # --- 4. Geração das três planilhas ------------------------------------------
# 
# # (A) Enriquecimento apenas por termos
# df_terms   <- filter_and_augment(up_down_enrich_results, pat_keywords, TRUE)
# write_by_category(df_terms,       terms_file)
# 
# # (B) Enriquecimento apenas por genes-alvo
# df_targets <- filter_and_augment(up_down_enrich_results, pat_genes,    FALSE)
# write_by_category(df_targets,     targets_file)
# 
# # (C) União (termos ∪ targets), sem duplicar entradas por category+ID
# df_combined <- bind_rows(df_terms, df_targets) %>%
#   distinct(category, ID, .keep_all = TRUE)
# write_by_category(df_combined,    combined_file)

# --- 1. Pacotes e Paths ------------------------------------------------------
# library(dplyr)
# library(stringr)
# library(purrr)
# library(openxlsx)
# 
# out_dir       <- "./Deseq2/SRP185421 - M tuberculosis/results/Enrichment_Analysis/"
# terms_file    <- file.path(out_dir, "enrichment_terms.xlsx")
# targets_file  <- file.path(out_dir, "enrichment_targets.xlsx")
# combined_file <- file.path(out_dir, "enrichment_combined.xlsx")
# dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
# 
# # --- 2. Definições -----------------------------------------------------------
# targets   <- c(
#   "NAT10","HAT1","KAT2A","KAT2B","KAT5","KAT6A","KAT6B",
#   "KAT7","KAT8","KAT12","GTF3C4","CREBBP","aTAT1","p300",
#   "HDAC1","HDAC2","HDAC3","HDAC4","HDAC5","HDAC6","HDAC7",
#   "HDAC8","HDAC9","HDAC10","SIRT1","SIRT2","SIRT3","SIRT4",
#   "SIRT5","SIRT6","SIRT7"
# )
# keywords  <- c(
#   "acetyltransferase","acetylation","desacetylation",
#   "histone lysine","histone","methyltransferase",
#   "deacetylase","histone deacetylase",
#   "lysine acetyltransferase","lysine"
# )
# 
# pat_genes    <- paste0("\\b(", paste(targets, collapse="|"), ")\\b")
# pat_keywords <- paste(keywords, collapse="|")
# targets_up   <- toupper(targets)
# 
# # --- 3. Funções Auxiliares ---------------------------------------------------
# 
# # Extrai apenas os genes-alvo de um geneID “A/B/C…”
# extract_targets <- function(geneID) {
#   parts <- str_split(geneID, "/", simplify = TRUE)
#   found <- parts[toupper(parts) %in% targets_up]
#   if (length(found) == 0) NA_character_ else paste(found, collapse = "/")
# }
# 
# # Reordena geneID colocando os alvos na frente
# reorder_genes <- function(geneID) {
#   parts <- str_split(geneID, "/", simplify = TRUE)
#   tgt    <- parts[toupper(parts) %in% targets_up]
#   other  <- parts[!toupper(parts) %in% targets_up]
#   paste(c(tgt, other), collapse = "/")
# }
# 
# # Filtra e anexa colunas comuns
# filter_and_augment <- function(results_list, pattern, on_desc = TRUE) {
#   map_dfr(names(results_list), function(cat) {
#     df <- as.data.frame(results_list[[cat]])
#     df$original_row <- seq_len(nrow(df))
#     df %>%
#       mutate(category = cat) %>%
#       filter(
#         if (on_desc) str_detect(Description, regex(pattern, ignore_case = TRUE))
#         else         str_detect(geneID,      regex(pattern, ignore_case = TRUE))
#       ) %>%
#       mutate(
#         only_targets   = vapply(geneID, extract_targets, FUN.VALUE = character(1)),
#         geneID_reorder = vapply(geneID, reorder_genes,   FUN.VALUE = character(1))
#       )
#   })
# }
# 
# # Escreve abas de um data.frame dividido por “category”
# write_by_category <- function(df, path) {
#   wb <- createWorkbook()
#   df %>% split(.$category) %>%
#     iwalk(~{
#       addWorksheet(wb, .y)
#       writeData(wb, .y, .x)
#     })
#   saveWorkbook(wb, path, overwrite = TRUE)
# }
# 
# # --- 4. Geração das três planilhas ------------------------------------------
# 
# # (A) Enriquecimento apenas por termos
# df_terms   <- filter_and_augment(up_down_enrich_results, pat_keywords, TRUE)
# write_by_category(df_terms,       terms_file)
# 
# # (B) Enriquecimento apenas por genes-alvo
# df_targets <- filter_and_augment(up_down_enrich_results, pat_genes,    FALSE)
# write_by_category(df_targets,     targets_file)
# 
# # (C) União (termos ∪ targets), sem duplicar entradas por category+ID
# df_combined <- bind_rows(df_terms, df_targets) %>%
#   distinct(category, ID, .keep_all = TRUE)
# write_by_category(df_combined,    combined_file)



