library(shiny)
library(DESeq2)
library(readxl)
library(dplyr)
library(pheatmap)
library(data.table)
library(ggplot2)
library(plotly)
library(htmlwidgets)
library(RColorBrewer) 
library(ggrepel)
library(tidyverse)
library(reshape2)
library(DOSE)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(enrichR)
library(biomaRt)
library(RWeka)
library(mice)
library(Hmisc)
library(VIM)
library(UpSetR)
library(sva)
library(writexl)
library(preprocessCore)



setwd("/Users/mortezaabyadeh/Documents/npc-fty/published papers")
npc <- read_excel("8w.xlsx")


dim(npc)
colnames(npc)[1] <- "Gene.names"

ctrl_names <- paste0("CTRL-", seq_len(3))

npc_names <- paste0("NPC-", seq_len(3))

colnames(npc)[2:4] <- npc_names
colnames(npc)[5:7] <- ctrl_names

colnames(npc)
head(npc)

boxplot(npc[, -1], main = "Boxplot of Proteome Intensities", 
        las = 2, outline = FALSE, col = "lightblue",
        ylab = "Protein Intensity")

pca <- prcomp(t(npc[, -1]), scale. = TRUE)
plot(pca$x[,1:2], col = rainbow(65), pch = 16, main = "PCA")





#npc_log2 <- npc[, -1]  
#npc_log2 <- log2(npc_log2)  

#npc_log2 <- cbind(Gene.names = npc$Gene.names, npc_log2)


#head(npc_log2)


#pca <- prcomp(t(npc[, -1]), scale. = TRUE)
#plot(pca$x[,1:2], col = rainbow(65), pch = 16, main = "PCA")






ctrl_columns <- grep("^CTRL", colnames(npc), value = TRUE)
npc_columns <- grep("^NPC", colnames(npc), value = TRUE)

geomean_ctrl <- apply(npc[, ctrl_columns], 1, function(x) exp(mean(log(x[x > 0]), na.rm = TRUE)))
geomean_npc <- apply(npc[, npc_columns], 1, function(x) exp(mean(log(x[x > 0]), na.rm = TRUE)))

npc$Geomean_CTRL <- geomean_ctrl
npc$Geomean_NPC <- geomean_npc

npc$Fold_Change <- npc$Geomean_NPC / npc$Geomean_CTRL

npc$Log2_FC <- log2(npc$Fold_Change)

p_values <- apply(npc[, c(ctrl_columns, npc_columns)], 1, function(x) {
  t_test_result <- t.test(x[1:length(ctrl_columns)], x[(length(ctrl_columns) + 1):length(x)])
  return(t_test_result$p.value)
})

npc$p_value <- p_values

npc$p_adj <- p.adjust(npc$p_value, method = "BH")

head(npc[, 5:ncol(npc)])

filtered_data <- npc[npc$p_adj < 0.05, ]
npc_data <- npc[npc$p_value < 0.05, ]
dim(npc_data)
# Get the dimensions of the filtered data
dim(filtered_data)
colnames(npc)[1] <- "Gene"
colnames(npc)[11] <- "log2FC"

write.xlsx(npc, "npc_datapw8.xlsx")


pca_maker <- function(df){
  df_filtered <- df[df$`p_adj` < 0.05, ]
  numberic_column_dataset1 <- df_filtered[,2:7]
  dataset1_numeric <- data.frame(lapply(numberic_column_dataset1, function(x) if(is.numeric(x)) as.integer(x) else x))
  
  dataset1_numeric <- na.omit(dataset1_numeric)
  colnames(dataset1_numeric)
  #gr <- sub("\\..*", "", colnames(dataset1_numeric))
  gr <- sub("\\-.*", "", colnames(dataset1_numeric))
  gr <- as.factor(gr)
  gr
  levels(gr)
  dataset1_numeric <- na.omit(dataset1_numeric)
  colnames(dataset1_numeric) <- gr
  #rownames(colData) <- colnames(dataset1_numeric)
  colnames(dataset1_numeric) <- gr
  dim(dataset1_numeric)
  cds <- DESeqDataSetFromMatrix(countData = dataset1_numeric, colData <- data.frame(Group = gr), design = ~ Group)
  #cds <- DESeq(cds)
  data.norm <- log2(1+counts(cds, normalized=F))
  data.mean.center <- t(scale(t(data.norm), scale = F))
  pc <- prcomp(data.mean.center)
  pcr <- data.frame(pc$r)
  pcr$group <- gr
  pc_var <- pc$sdev^2 / sum(pc$sdev^2)
  
  ggplot(pcr, aes(PC1, PC2, color=group)) + geom_point(size=10, alpha=0.6, shape = 19) + theme_bw() + 
    theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
          axis.text = element_text(size = 13),
          axis.title = element_text(size = 15, face = "bold"),
          legend.title = element_text(size = 15, face = "bold"),
          legend.text = element_text(size = 15),
          legend.position = "right") + 
    labs(
      x = paste0("PC1 (", round(pc_var[1] * 100, 2), "% variance)"),
      y = paste0("PC2 (", round(pc_var[2] * 100, 2), "% variance)"),
      title = "PCA Plot"
    )
}

heatmap_maker <- function(df){
  df <- df
  df_filtered <- df[df$`p_adj` < 0.05, ]
  df_filtered <- df_filtered[complete.cases(df_filtered), ]
  i_subset <- df_filtered[,2:7]
  i_matrix <- as.matrix(i_subset)
  i_log <- log2(i_matrix + 1)
  i <- t(scale(t(i_log), scale = FALSE))
  i_annotation_col <- data.frame(Group= c(sub("-.*", "", colnames(i))))
  rownames(i_annotation_col) <- colnames(i)
  i_head <- pheatmap(i,
                     cluster_cols = TRUE,  
                     scale = "row",  
                     color = colorRampPalette(c("green", "black", "red"))(100), 
                     main = "Heatmap of Normalized Counts",  
                     fontsize = 13,  
                     fontsize_row = 12,  
                     fontsize_col = 13,  
                     border_color = NA,  
                     #angle_col = 45,
                     annotation_col = i_annotation_col,
                     width = 7,
                     height = 7,
                     show_rownames = FALSE
  )
}


KEGG_maker_down <- function(df){
  df$diffexpressed <- ifelse(df$`p_adj` < 0.05 & df$`log2FC` > 0, "Up",
                             ifelse(df$`p_adj` < 0.05 & df$`log2FC` < 0, "Down", "No"))
  df$delabel <- df$Gene
  df <- df[complete.cases(df), ]
  df <- df[df$diffexpressed=="Down", ]
  gene <- df$Gene
  gene <- as.data.frame(gene)
  gene <- as.list(gene)
  class(gene)
  listEnrichrDbs()
  dbs <- listEnrichrDbs()
  db_names <- dbs$libraryName
  selected_dbs <- c("KEGG_2021_Human")
  enriched <- enrichr(df$Gene, selected_dbs)
  kegg_result <- as.data.frame(enriched[["KEGG_2021_Human"]])
  kegg_filtered <- kegg_result[kegg_result$Adjusted.P.value < 0.05, ]
  kegg_sorted <- kegg_filtered[order(kegg_filtered$Adjusted.P.value), ]
  top_10_kegg <- head(kegg_sorted, 20)
  ggplot(top_10_kegg, aes(x = reorder(Term, -Adjusted.P.value), y = -log10(Adjusted.P.value))) +
    geom_bar(stat = "identity", fill = "steelblue") +
    coord_flip() +
    xlab("KEGG Pathways") +
    ylab("-log10(Adjusted P-value)") +
    ggtitle("KEGG Pathway Enrichment Analysis") +
    theme_minimal() +
    theme(axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14),
          axis.text.y = element_text(size = 12))
}

KEGG_maker_up <- function(df){
  df$diffexpressed <- ifelse(df$`p_adj` < 0.05 & df$`log2FC` > 0, "Up",
                             ifelse(df$`p_adj` < 0.05 & df$`log2FC` < 0, "Down", "No"))
  df$delabel <- df$Gene
  df <- df[complete.cases(df), ]
  df <- df[df$diffexpressed=="Up", ]
  gene <- df$Gene
  gene <- as.data.frame(gene)
  gene <- as.list(gene)
  class(gene)
  listEnrichrDbs()
  dbs <- listEnrichrDbs()
  db_names <- dbs$libraryName
  selected_dbs <- c("KEGG_2021_Human")
  enriched <- enrichr(df$Gene, selected_dbs)
  kegg_result <- as.data.frame(enriched[["KEGG_2021_Human"]])
  kegg_filtered <- kegg_result[kegg_result$Adjusted.P.value < 0.05, ]
  kegg_sorted <- kegg_filtered[order(kegg_filtered$Adjusted.P.value), ]
  top_10_kegg <- head(kegg_sorted, 20)
  ggplot(top_10_kegg, aes(x = reorder(Term, -Adjusted.P.value), y = -log10(Adjusted.P.value))) +
    geom_bar(stat = "identity", fill = "red") +
    coord_flip() +
    xlab("KEGG Pathways") +
    ylab("-log10(Adjusted P-value)") +
    ggtitle("KEGG Pathway Enrichment Analysis") +
    theme_minimal() +
    theme(axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14),
          axis.text.y = element_text(size = 12))
}



volcano_maker <- function(df){
  df <- df[!df$p_adj == 0, ]
  df$diffexpressed <- ifelse(df$`p_adj` < 0.05 & df$`log2FC` > 0, "Up",
                             ifelse(df$`p_adj` < 0.05 & df$`log2FC` < 0, "Down", "No"))
  df$delabel <- df$Gene
  df <- df[complete.cases(df), ]
  counts <- table(df$diffexpressed)
  max_y <- round(max(-log10(df$p_adj)))
  max_x <- round(max(df$log2FC))
  downregulated_count <- counts["Down"]
  not_significant_count <- counts["No"]
  upregulated_count <- counts["Up"]
  
  ggplot(data = df, aes(x = log2FC, y = -log10(p_adj), col = diffexpressed)) +
    geom_vline(xintercept = c(-0.6, 0.6), col = "gray", linetype = 'dashed') +
    geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') + 
    geom_point(size = 3) + 
    scale_color_manual(values = c("#00AFBB", "grey", "#bb0c00"), 
                       labels = c(paste("Downregulated (", downregulated_count, ")", sep = ""),
                                  paste("Not significant (", not_significant_count, ")", sep = ""),
                                  paste("Upregulated (", upregulated_count, ")", sep = ""))) +
    coord_cartesian(ylim = c(0, max_y), xlim = c(-max_x, max_x)) +
    scale_x_continuous(breaks = seq(-max_x, max_x, 1)) + 
    ggtitle('Volcano plot') +
    theme_bw() + 
    theme(legend.title = element_blank(),
          legend.text = element_text(size = 15)) 
}







head(npc)
heatmap_maker(npc)
pca_maker(npc)
volcano_maker(npc)

KEGG_maker_down(npc)
KEGG_maker_up(npc)




top_kegg_up <- KEGG_maker_up(npc)

write_xlsx(top_kegg_up, path = "KEGG_Upregulated_Top20.xlsx")



top_kegg_down <- KEGG_maker_down(npc)

write_xlsx(top_kegg_down, path = "KEGG_Downregulated_Top20.xlsx")


























### normalize ###### dont need to normalize 

log_expr <- log2(as.matrix(npc[, -1]) + 1)

norm_expr <- normalize.quantiles(log_expr)

batch <- rep(1:(ncol(norm_expr)/4), each = 4)

#### batch effect correction ###### same batch and pca showed that dont need batch effcet correction

combat_expr <- ComBat(dat = norm_expr, batch = batch, par.prior = TRUE)

#Combine with gene names
npc_corrected <- cbind(Gene.names = npc$Gene.names, as.data.frame(combat_expr))



head(npc_corrected)

expr_matrix1 <- as.matrix(npc_corrected[, -1])

annotation_col <- data.frame(
  Batch = batch
)
rownames(annotation_col) <- colnames(expr_matrix1)

pheatmap(
  cor(expr_matrix1),                  
  annotation_col = annotation_col,  
  main = "Batch Effect Check",
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  show_rownames = TRUE,
  show_colnames = TRUE
)



pca <- prcomp(t(npc_corrected[, -1]), scale. = TRUE)
plot(pca$x[,1:2], col = rainbow(65), pch = 16, main = "PCA")


###### batch effect correction without normalization


expr_matrix <- as.matrix(npc[, -1])
ncol(expr_matrix)
head(expr_matrix)

batch <- rep(1:(ncol(expr_matrix)/4), each = 4)
log_expr <- log2(expr_matrix + 1)  

combat_expr <- ComBat(dat = log_expr, batch = batch, par.prior = TRUE, prior.plots = FALSE)
npc_corrected <- cbind(Gene.names = npc$Gene.names, as.data.frame(combat_expr))

pca <- prcomp(t(npc_corrected[, -1]), scale. = TRUE)
plot(pca$x[,1:2], col = rainbow(65), pch = 16, main = "PCA")

###heatmap to check batch effcet

annotation_col <- data.frame(
  Batch = batch
)
rownames(annotation_col) <- colnames(expr_matrix)

pheatmap(
  cor(expr_matrix),                  # correlation heatmap of samples
  annotation_col = annotation_col,  # annotate columns with batch info
  main = "Batch Effect Check",
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  show_rownames = TRUE,
  show_colnames = TRUE
)









colnames(npc)


pca_maker <- function(df){
  df_filtered <- df[df$`p_adj` < 0.05, ]
  numberic_column_dataset1 <- df_filtered[,2:65]
  dataset1_numeric <- data.frame(lapply(numberic_column_dataset1, function(x) if(is.numeric(x)) as.integer(x) else x))
  
  dataset1_numeric <- na.omit(dataset1_numeric)
  colnames(dataset1_numeric)
  #gr <- sub("\\..*", "", colnames(dataset1_numeric))
  gr <- sub("\\-.*", "", colnames(dataset1_numeric))
  gr <- as.factor(gr)
  gr
  levels(gr)
  dataset1_numeric <- na.omit(dataset1_numeric)
  colnames(dataset1_numeric) <- gr
  #rownames(colData) <- colnames(dataset1_numeric)
  colnames(dataset1_numeric) <- gr
  dim(dataset1_numeric)
  cds <- DESeqDataSetFromMatrix(countData = dataset1_numeric, colData <- data.frame(Group = gr), design = ~ Group)
  #cds <- DESeq(cds)
  data.norm <- log2(1+counts(cds, normalized=F))
  data.mean.center <- t(scale(t(data.norm), scale = F))
  pc <- prcomp(data.mean.center)
  pcr <- data.frame(pc$r)
  pcr$group <- gr
  pc_var <- pc$sdev^2 / sum(pc$sdev^2)
  
  ggplot(pcr, aes(PC1, PC2, color=group)) + geom_point(size=10, alpha=0.6, shape = 19) + theme_bw() + 
    theme(plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
          axis.text = element_text(size = 13),
          axis.title = element_text(size = 15, face = "bold"),
          legend.title = element_text(size = 15, face = "bold"),
          legend.text = element_text(size = 15),
          legend.position = "right") + 
    labs(
      x = paste0("PC1 (", round(pc_var[1] * 100, 2), "% variance)"),
      y = paste0("PC2 (", round(pc_var[2] * 100, 2), "% variance)"),
      title = "PCA Plot"
    )
}



#### differentially expressed anlayis after normalization



npc_log2 <- npc[, -1]  # Exclude the Gene.names column (assuming it's the first column)
npc_log2 <- log2(npc_log2)  # Apply log2 transformation

# Add the Gene.names column back
npc_log2 <- cbind(Gene.names = npc$Gene.names, npc_log2)

# View the normalized data
head(npc_log2)


#expr_matrix <- as.matrix(npc_log2[, -1])

#batch <- rep(1:(ncol(expr_matrix)/4), each = 4)
#log_expr <- log2(expr_matrix + 1)  

#combat_expr <- ComBat(dat = log_expr, batch = batch, par.prior = TRUE, prior.plots = FALSE)
#npc_corrected <- cbind(Gene.names = npc$Gene.names, as.data.frame(combat_expr))

pca <- prcomp(t(npc[, -1]), scale. = TRUE)
plot(pca$x[,1:2], col = rainbow(65), pch = 16, main = "PCA")






ctrl_columns <- grep("^CTRL", colnames(npc), value = TRUE)
npc_columns <- grep("^NPC", colnames(npc), value = TRUE)

geomean_ctrl <- apply(npc[, ctrl_columns], 1, function(x) exp(mean(log(x[x > 0]), na.rm = TRUE)))
geomean_npc <- apply(npc[, npc_columns], 1, function(x) exp(mean(log(x[x > 0]), na.rm = TRUE)))

npc$Geomean_CTRL <- geomean_ctrl
npc$Geomean_NPC <- geomean_npc

npc$Fold_Change <- npc$Geomean_NPC / npc$Geomean_CTRL

npc$Log2_FC <- log2(npc$Fold_Change)

p_values <- apply(npc[, c(ctrl_columns, npc_columns)], 1, function(x) {
  t_test_result <- t.test(x[1:length(ctrl_columns)], x[(length(ctrl_columns) + 1):length(x)])
  return(t_test_result$p.value)
})

npc$p_value <- p_values

npc$p_adj <- p.adjust(npc$p_value, method = "BH")

head(npc)

filtered_data <- npc[npc$p_adj < 0.05, ]

# Get the dimensions of the filtered data
dim(filtered_data)
pca_maker(npc)



heatmap_maker <- function(df){
  df <- df
  df_filtered <- df[df$`p_adj` < 0.05, ]
  df_filtered <- df_filtered[complete.cases(df_filtered), ]
  i_subset <- df_filtered[,2:65]
  i_matrix <- as.matrix(i_subset)
  i_log <- log2(i_matrix + 1)
  i <- t(scale(t(i_log), scale = FALSE))
  i_annotation_col <- data.frame(Group= c(sub("-.*", "", colnames(i))))
  rownames(i_annotation_col) <- colnames(i)
  i_head <- pheatmap(i,
                     cluster_cols = TRUE,  
                     scale = "row",  
                     color = colorRampPalette(c("green", "black", "red"))(100), 
                     main = "Heatmap of Normalized Counts",  
                     fontsize = 13,  
                     fontsize_row = 12,  
                     fontsize_col = 13,  
                     border_color = NA,  
                     #angle_col = 45,
                     annotation_col = i_annotation_col,
                     width = 7,
                     height = 7,
                     show_rownames = FALSE
  )
}
heatmap_maker(npc)





### get the DAP after batch effects correction ####

data <- read_excel("Proteome data-sci-npc1.xlsx")

selected_columns <- c(
  "Gene.names",
  "WC0004", "WC0001", "WC0002", "WC0003",
  "WC0019", "WC0017", "WC0018", "WC0020",
  "WC0035", "WC0033", "WC0034", "WC0036",
  "WC0050", "WC0049", "WC0051", "WC0052",
  "WC0068", "WC0065", "WC0066", "WC0067",
  "WC0081", "WC0084", "WC0083", "WC0082",
  "WC0098", "WC0097", "WC0099", "WC0100",
  "WC0113", "WC0115", "WC0116", "WC0114",
  "WC0130", "WC0129", "WC0132", "WC0131",
  "WC0148", "WC0146", "WC0147", "WC0145",
  "WC0164", "WC0162", "WC0163", "WC0161",
  "WC0183", "WC0181", "WC0182", "WC0184",
  "WC0208", "WC0207", "WC0206", "WC0205",
  "WC0231", "WC0230", "WC0229", "WC0232",
  "WC0255", "WC0253", "WC0256", "WC0254",
  "WC0093", "WC0094", "WC0096", "WC0095"
)


npc <- data[, selected_columns]
dim(npc)


colnames(npc)[1] <- "Gene.names"

ctrl_names <- paste0("CTRL-", seq_len(60))

npc_names <- paste0("NPC-", seq_len(4))

colnames(npc)[2:61] <- ctrl_names
colnames(npc)[62:65] <- npc_names


library(openxlsx)
write.xlsx(npc, "npc_data.xlsx") 




### batch effcet
expr_matrix <- as.matrix(npc[, -1])
ncol(expr_matrix)
head(expr_matrix)

batch <- rep(1:(ncol(expr_matrix)/4), each = 4)
log_expr <- log2(expr_matrix + 1)  

combat_expr <- ComBat(dat = log_expr, batch = batch, par.prior = TRUE, prior.plots = FALSE)
npc_corrected <- cbind(Gene.names = npc$Gene.names, as.data.frame(combat_expr))

pca <- prcomp(t(npc_corrected[, -1]), scale. = TRUE)
plot(pca$x[,1:2], col = rainbow(65), pch = 16, main = "PCA")

head(npc_corrected)


















expr_matrix <- as.matrix(npc[, -1])
ncol(expr_matrix)
head(expr_matrix)

batch <- rep(1:(ncol(expr_matrix)/4), each = 4)
log_expr <- log2(expr_matrix + 1)  

combat_expr <- ComBat(dat = log_expr, batch = batch, par.prior = TRUE, prior.plots = FALSE)
npc_corrected <- cbind(Gene.names = npc$Gene.names, as.data.frame(combat_expr))

pca <- prcomp(t(npc_corrected[, -1]), scale. = TRUE)
plot(pca$x[,1:2], col = rainbow(65), pch = 16, main = "PCA")







boxplot(npc_corrected[, -1], main = "Boxplot of Proteome Intensities", 
        las = 2, outline = FALSE, col = "lightblue",
        ylab = "Protein Intensity")



ctrl_columns <- grep("^CTRL", colnames(npc_corrected), value = TRUE)
npc_columns <- grep("^NPC", colnames(npc_corrected), value = TRUE)

geomean_ctrl <- apply(npc_corrected[, ctrl_columns], 1, function(x) exp(mean(log(x[x > 0]), na.rm = TRUE)))
geomean_npc <- apply(npc_corrected[, npc_columns], 1, function(x) exp(mean(log(x[x > 0]), na.rm = TRUE)))

npc_corrected$Geomean_CTRL <- geomean_ctrl
npc_corrected$Geomean_NPC <- geomean_npc

npc_corrected$Fold_Change <- npc_corrected$Geomean_NPC / npc_corrected$Geomean_CTRL

npc_corrected$Log2_FC <- log2(npc_corrected$Fold_Change)

p_values <- apply(npc_corrected[, c(ctrl_columns, npc_columns)], 1, function(x) {
  t_test_result <- t.test(x[1:length(ctrl_columns)], x[(length(ctrl_columns) + 1):length(x)])
  return(t_test_result$p.value)
})

npc_corrected$p_value <- p_values

npc_corrected$p_adj <- p.adjust(npc_corrected$p_value, method = "BH")

head(npc_corrected)

filtered_data <- npc_corrected[npc_corrected$p_adj < 0.05, ]

# Get the dimensions of the filtered data
dim(filtered_data)
pca_maker(npc)




