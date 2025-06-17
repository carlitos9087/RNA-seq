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
base_dir <- "Deseq2/GSE148796(homo sapiens)- HIV/"
results_dir <- file.path(base_dir, "results")

tabular_dir <- base_dir
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
phenoData <- read_excel("./Deseq2/GSE148796(homo sapiens)- HIV/GSE148796.xlsx", col_names = TRUE)
lista <- phenoData$id
# phenoData <- phenoData[,-1]
phenoData <- phenoData[,-1]# Remove a coluna "id"
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


run_volcano_plot <- function(dds, output_file, pval_threshold = 0.05, log2fc_threshold = 1) {
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  
  res <- results(dds)
  res_df <- na.omit(as.data.frame(res))
  res_df$GeneID <- rownames(res_df)  # Mantém os IDs Entrez
  
  # Converte IDs Entrez para símbolos de genes
  res_df$Gene <- convert_entrez_to_symbol(res_df$GeneID)
  
  res_df$Significance <- "Not Significant"
  res_df$Significance[res_df$pvalue < pval_threshold & abs(res_df$log2FoldChange) >= log2fc_threshold] <- "Significant"
  res_df$Color <- ifelse(res_df$Significance == "Significant", ifelse(res_df$log2FoldChange > 0, "red", "blue"), "grey")
  
  write.xlsx(res_df, file = sub(".jpeg", ".xlsx", output_file), rowNames = FALSE)
  
  top_genes <- res_df[res_df$Significance == "Significant" & res_df$pvalue < 0.01, ]
  
  p <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(pvalue), color = Color, label = Gene)) +
    geom_point(alpha = 0.7, size = 2) +
    scale_color_manual(values = c("red" = "red", "blue" = "blue", "grey" = "grey")) +
    geom_hline(yintercept = -log10(pval_threshold), linetype = "dashed", color = "black") +
    geom_vline(xintercept = c(-log2fc_threshold, log2fc_threshold), linetype = "dashed", color = "black") +
    xlab("Log2 Fold Change") +
    ylab("-Log10 P-value") +
    ggtitle("Volcano Plot") +
    theme_minimal() +
    theme(legend.position = "none") +
    geom_text_repel(data = top_genes, aes(label = Gene), size = 3, max.overlaps = 10)
  
  ggsave(output_file, p, width = 8, height = 6, dpi = 300)
  return(p)
}

run_ma_plot <- function(dds, output_file) {
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  
  res <- results(dds)
  res_df <- as.data.frame(res)
  res_df$GeneID <- rownames(res_df)  # Mantém os IDs Entrez
  
  # Converte IDs Entrez para símbolos de genes
  res_df$Gene <- convert_entrez_to_symbol(res_df$GeneID)
  
  res_df$Significance <- "Not Significant"
  res_df$Significance[res_df$padj < 0.05] <- "Significant"
  res_df$Color <- "grey"
  res_df$Color[res_df$padj < 0.05 & res_df$log2FoldChange > 0] <- "red"
  res_df$Color[res_df$padj < 0.05 & res_df$log2FoldChange < 0] <- "blue"
  
  write.xlsx(res_df, file = sub(".jpeg", ".xlsx", output_file), rowNames = FALSE)
  
  top_genes <- res_df[abs(res_df$log2FoldChange) > 2 & res_df$padj < 0.05, ]
  
  p <- ggplot(res_df, aes(x = log10(baseMean), y = log2FoldChange, color = Color, label = Gene)) +
    geom_point(alpha = 0.7, size = 2) +
    scale_color_manual(values = c("red" = "red", "blue" = "blue", "grey" = "grey")) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    xlab("Log10 Mean Expression") +
    ylab("Log2 Fold Change") +
    ggtitle("MA-plot") +
    theme_minimal() +
    theme(legend.position = "none") +
    geom_text_repel(data = top_genes, aes(label = Gene), size = 3, max.overlaps = 10)
  
  ggsave(output_file, p, width = 8, height = 6, dpi = 300)
  return(p)
}

# phenoData
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

intersect(phenoData$SampleName[phenoData$Treatment %in% c("Mock", "HIV-1_infected")], colnames(data))
colnames(data)
phenoData

# Executa DESeq2 para a comparação Mock vs HIV-1_infected
dds_Mock_vs_HIV1_infected <- run_deseq_analysis(
  data, phenoData, 
  "Mock",            # Grupo controle
  "HIV-1_infected",  # Grupo experimental
  file.path(results_dir, "DESeq2_Mock_vs_HIV-1_infected.tabular")
)

# Análises de PCA
pca_Mock_vs_HIV1 <- run_pca_analysis(
  dds_Mock_vs_HIV1_infected, 
  file.path(results_dir, "PCA_Mock_vs_HIV-1_infected.jpeg")
)

# Distâncias entre amostras
run_sample_distances(
  dds_Mock_vs_HIV1_infected, 
  file.path(results_dir, "SampleDistances_Mock_vs_HIV-1_infected.jpeg")
)

# Gráficos de dispersão
run_dispersion_plot(
  dds_Mock_vs_HIV1_infected, 
  file.path(results_dir, "Dispersion_Mock_vs_HIV-1_infected.jpeg")
)

# Histogramas de p-values
run_pval_histogram(
  dds_Mock_vs_HIV1_infected, 
  file.path(results_dir, "Histogram_pvalues_Mock_vs_HIV-1_infected.jpeg")
)

# MA Plot
ma_Mock_vs_HIV1 <- run_ma_plot(
  dds_Mock_vs_HIV1_infected, 
  file.path(results_dir, "MAplot_Mock_vs_HIV-1_infected.jpeg")
)

# Volcano Plot
run_volcano_plot(
  dds_Mock_vs_HIV1_infected, 
  file.path(results_dir, "Volcano_Mock_vs_HIV-1_infected.jpeg")
)


############################

library(openxlsx)
library(clusterProfiler)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
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
                                         up_threshold   = 1,
                                         down_threshold = -1,
                                         padj_cutoff    = 0.05,
                                         ontologies     = c("BP", "MF"),
                                         pAdjustMethod  = "BH",
                                         pvalueCutoff   = 0.05,
                                         qvalueCutoff   = 0.2,
                                         output_prefix) {
  
  # (1) Extrai DEGs e salva em arquivo
  degs <- run_deseq_up_down_analysis(dds,
                                     up_threshold   = up_threshold,
                                     down_threshold = down_threshold,
                                     padj_cutoff    = padj_cutoff,
                                     output_file    = paste0(output_prefix, "_DEGs.xlsx"))
  
  up_ids   <- na.omit(as.character(degs$up$Entrez_ID))
  down_ids <- na.omit(as.character(degs$down$Entrez_ID))
  
  if (length(up_ids) == 0 && length(down_ids) == 0) {
    stop("Nenhum gene upregulated ou downregulated encontrado para enriquecimento.")
  }
  
  # Determina o OrgDb
  orgDb <- select_orgDb(unique(c(up_ids, down_ids)))
  
  # Inicializa lista de resultados e workbook
  enrich_results <- list()
  wb <- createWorkbook()
  
  for (ont in ontologies) {
    
    # UPREGULATED
    if (length(up_ids) > 0) {
      ego_up <- tryCatch({
        enrichGO(gene          = up_ids,
                 OrgDb         = orgDb,
                 keyType       = "ENTREZID",
                 ont           = ont,
                 pAdjustMethod = pAdjustMethod,
                 pvalueCutoff  = pvalueCutoff,
                 qvalueCutoff  = qvalueCutoff,
                 readable      = TRUE)
      }, error = function(e) NULL)
      
      if (!is.null(ego_up) && nrow(as.data.frame(ego_up)) > 0) {
        enrich_results[[paste0("up_", ont)]] <- ego_up
        
        sheet_up <- paste0("Up_", ont)
        addWorksheet(wb, sheet_up)
        writeData(wb, sheet_up, as.data.frame(ego_up))
        pdf(paste0(output_prefix, "_", sheet_up, ".pdf"), width = 10, height = 8)
        print(dotplot(ego_up, showCategory = 20) + ggtitle(paste("GO", ont, "Enrichment for Upregulated Genes")))
        dev.off()
      } else {
        message("Nenhum termo enriquecido encontrado para genes upregulated (", ont, ").")
      }
    }
    
    # DOWNREGULATED
    if (length(down_ids) > 0) {
      ego_down <- tryCatch({
        enrichGO(gene          = down_ids,
                 OrgDb         = orgDb,
                 keyType       = "ENTREZID",
                 ont           = ont,
                 pAdjustMethod = pAdjustMethod,
                 pvalueCutoff  = pvalueCutoff,
                 qvalueCutoff  = qvalueCutoff,
                 readable      = TRUE)
      }, error = function(e) NULL)
      
      if (!is.null(ego_down) && nrow(as.data.frame(ego_down)) > 0) {
        enrich_results[[paste0("down_", ont)]] <- ego_down
        
        sheet_down <- paste0("Down_", ont)
        addWorksheet(wb, sheet_down)
        writeData(wb, sheet_down, as.data.frame(ego_down))
        pdf(paste0(output_prefix, "_", sheet_down, ".pdf"), width = 10, height = 8)
        print(dotplot(ego_down, showCategory = 20) + ggtitle(paste("GO", ont, "Enrichment for Downregulated Genes")))
        dev.off()
      } else {
        message("Nenhum termo enriquecido encontrado para genes downregulated (", ont, ").")
      }
    }
  }
  
  # Salva workbook apenas se tiver conteúdo
  if (length(wb$worksheets) > 0) {
    saveWorkbook(wb, paste0(output_prefix, "_Enrichment_GO.xlsx"), overwrite = TRUE)
  } else {
    message("Nenhum resultado de enriquecimento foi salvo no Excel.")
  }
  
  message("GO enrichment (BP & MF) finalizado. Prefixo: ", output_prefix)
  return(enrich_results)
}



up_down_enrich_dds_Mock_vs_HIV1_infected <- run_deseq_up_down_enrichment(
  dds             = dds_Mock_vs_HIV1_infected,
  up_threshold    = 1,
  down_threshold  = -1,
  padj_cutoff     = 0.05,
  ontologies      = c("BP","MF"),
  pAdjustMethod   = "BH",
  pvalueCutoff    = 0.05,
  qvalueCutoff    = 0.2,
  output_prefix = file.path(results_dir, "DESeq2_UpDown_dds_Mock_vs_HIV1_infected")
  
)

results_dir

resultados <- list(
  dds_Mock_vs_HIV1_infected = up_down_enrich_dds_Mock_vs_HIV1_infected
)

# Diretórios de saída correspondentes
diretorios_saida <- list(
  dds_Mock_vs_HIV1_infected = file.path(results_dir, "Enrichment_Analysis_UpDown_dds_Mock_vs_HIV1_infected")
  
)

# Alvos de interesse
targets <- c(
  "NAT10","HAT1","KAT2A","KAT2B","KAT5","KAT6A","KAT6B",
  "KAT7","KAT8","KAT12","GTF3C4","CREBBP","aTAT1","p300",
  "HDAC1","HDAC2","HDAC3","HDAC4","HDAC5","HDAC6","HDAC7",
  "HDAC8","HDAC9","HDAC10","SIRT1","SIRT2","SIRT3","SIRT4",
  "SIRT5","SIRT6","SIRT7"
)

# Palavras-chave para filtragem
keywords <- c(
  "acetyltransferase","acetylation","desacetylation",
  "histone lysine","histone","methyltransferase",
  "deacetylase","histone deacetylase",
  "lysine acetyltransferase","lysine"
)

# Padrões regex para filtragem
pat_genes <- paste0("\\b(", paste(targets, collapse = "|"), ")\\b")
pat_keywords <- paste(keywords, collapse = "|")
targets_up <- toupper(targets)

library(dplyr)
library(stringr)
library(purrr)
library(openxlsx)

# Extrai apenas os genes-alvo de um geneID “A/B/C…”
extract_targets <- function(geneID) {
  parts <- str_split(geneID, "/", simplify = TRUE)
  found <- parts[toupper(parts) %in% targets_up]
  if (length(found) == 0) NA_character_ else paste(found, collapse = "/")
}

# Reordena geneID colocando os alvos na frente
reorder_genes <- function(geneID) {
  parts <- str_split(geneID, "/", simplify = TRUE)
  tgt <- parts[toupper(parts) %in% targets_up]
  other <- parts[!toupper(parts) %in% targets_up]
  paste(c(tgt, other), collapse = "/")
}

# Filtra e anexa colunas comuns
filter_and_augment <- function(results_list, pattern, on_desc = TRUE) {
  map_dfr(names(results_list), function(cat) {
    df <- as.data.frame(results_list[[cat]])
    df$original_row <- seq_len(nrow(df))
    df %>%
      mutate(category = cat) %>%
      filter(
        if (on_desc) str_detect(Description, regex(pattern, ignore_case = TRUE))
        else         str_detect(geneID, regex(pattern, ignore_case = TRUE))
      ) %>%
      mutate(
        only_targets = vapply(geneID, extract_targets, FUN.VALUE = character(1)),
        geneID_reorder = vapply(geneID, reorder_genes, FUN.VALUE = character(1))
      )
  })
}

# Escreve abas de um data.frame dividido por “category”
write_by_category <- function(df, path) {
  wb <- createWorkbook()
  df %>% split(.$category) %>%
    iwalk(~{
      addWorksheet(wb, .y)
      writeData(wb, .y, .x)
    })
  saveWorkbook(wb, path, overwrite = TRUE)
}

# Itera sobre cada conjunto de resultados
walk(names(resultados), function(nome) {
  resultado <- resultados[[nome]]
  out_dir <- diretorios_saida[[nome]]
  
  # Cria o diretório de saída, se não existir
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Define os caminhos dos arquivos de saída
  terms_file <- file.path(out_dir, "enrichment_terms.xlsx")
  targets_file <- file.path(out_dir, "enrichment_targets.xlsx")
  combined_file <- file.path(out_dir, "enrichment_combined.xlsx")
  
  # (A) Enriquecimento apenas por termos
  df_terms <- filter_and_augment(resultado, pat_keywords, TRUE)
  write_by_category(df_terms, terms_file)
  
  # (B) Enriquecimento apenas por genes-alvo
  df_targets <- filter_and_augment(resultado, pat_genes, FALSE)
  write_by_category(df_targets, targets_file)
  
  # (C) União (termos ∪ targets), sem duplicar entradas por category+ID
  df_combined <- bind_rows(df_terms, df_targets) %>%
    distinct(category, ID, .keep_all = TRUE)
  write_by_category(df_combined, combined_file)
})







