# Author: Komal S. Rathi
# Function: PCA plot of voom normalized data
# Date: 04/03/2020

library(Rtsne)
library(ggpubr)

tsne.plot <- function(voomData, gene_list = NULL, topVar = 500, meta, fname, plx, color_var, shape_var){
  
  # subset voom normalized data using meta file
  rownames(meta) <- meta$sample
  voomData <- voomData[,rownames(meta)]
  
  # check if meta and expression are compatible
  if(identical(rownames(meta), colnames(voomData))) {
    print("Proceed")
  } else {
    break
  }
  
  # t-SNE before combat adjustment
  set.seed(42)
  tsneOut <- Rtsne(t(voomData), initial_dims = 50, perplexity = plx, max_iter = 1000)
  tsneOut <- data.frame(tsneOut$Y, meta)
  p <- ggplot(tsneOut, aes(X1, X2)) +
    geom_point(size = 5, alpha = 0.5, aes_string(color = color_var, shape = shape_var)) +
    geom_text(aes(label = sample), size = 2) + 
    theme_bw() +
    ggtitle("All expressed genes") +
    theme_Publication2() + xlab("PC1") + ylab("PC2") 
  res <- p
  
  pdf(file = fname, width = 12, height = 5, onefile = F)
  # if gene_list is present use that
  if(!is.null(gene_list)){
    voomData <- voomData[which(rownames(voomData) %in% gene_list),]
    tsneOut <- Rtsne(t(voomData), initial_dims = 50, perplexity = plx, max_iter = 1000)
    tsneOut <- data.frame(tsneOut$Y, meta)
    q <- ggplot(tsneOut, aes(X1, X2)) +
      geom_point(size = 5, alpha = 0.5, aes_string(color = color_var, shape = shape_var)) +
      geom_text(aes(label = sample), size = 2) + 
      theme_bw() +
      ggtitle("Housekeeping genes") +
      theme_Publication2() + xlab("PC1") + ylab("PC2")
    res <- ggarrange(p, q, common.legend = T)
  } 
  
  # if topVar is present use that
  if(!is.null(gene_list)){
    # use the top variable genes
    rv <- rowVars(voomData)
    select <- order(rv, decreasing=TRUE)[seq_len(topVar)]
    voomData <- voomData[select,]
    tsneOut <- Rtsne(t(voomData), initial_dims = 50, perplexity = plx, max_iter = 1000)
    tsneOut <- data.frame(tsneOut$Y, meta)
    q <- ggplot(tsneOut, aes(X1, X2)) +
      geom_point(size = 5, alpha = 0.5, aes_string(color = color_var, shape = shape_var)) +
      geom_text(aes(label = sample), size = 2) + 
      theme_bw() +
      ggtitle("Most variable genes") +
      theme_Publication2() + xlab("PC1") + ylab("PC2")
    res <- ggarrange(p, q, common.legend = T)
  } 
  
  # print plot
  print(res)
  dev.off()
}

pca.plot <- function(voomData, gene_list = NULL, topVar = 500, meta, fname, color_var, shape_var){
  
  # subset voom normalized data using meta file
  rownames(meta) <- meta$sample
  voomData <- voomData[,rownames(meta)]
  
  # check if meta and expression are compatible
  if(identical(rownames(meta), colnames(voomData))) {
    print("Proceed")
  } else {
    break
  }
  
  # pca
  prData <- prcomp(voomData)
  pca.data <- prData$rotation
  pca.data <- data.frame(pca.data)[1:4]
  pca.data <- data.frame(pca.data, meta)
  pdf(file = fname, width = 12, height = 8, onefile = F)
  p <- ggplot(pca.data, aes(PC1, PC2)) +
    geom_point(size = 5, alpha = 0.5, aes_string(color = color_var, shape = shape_var)) +
    geom_text(aes(label = sample), size = 2) + 
    theme_bw() +
    ggtitle("All expressed genes") +
    theme_Publication2()  
  q <- ggplot(pca.data, aes(PC3, PC4)) +
    geom_point(size = 5, alpha = 0.5, aes_string(color = color_var, shape = shape_var)) +
    geom_text(aes(label = sample), size = 2) + 
    theme_bw() +
    ggtitle("All expressed genes") +
    theme_Publication2() 
  res <- ggarrange(p, q, ncol = 2, nrow = 2, common.legend = T)
  
  # if gene_list is present use that
  if(!is.null(gene_list)){
    voomData <- voomData[which(rownames(voomData) %in% gene_list),]
    prData <- prcomp(voomData)
    pca.data <- prData$rotation
    pca.data <- data.frame(pca.data)[1:4]
    pca.data <- data.frame(pca.data, meta)
    r <- ggplot(pca.data, aes(PC1, PC2)) +
      geom_point(size = 5, alpha = 0.5, aes_string(color = color_var, shape = shape_var)) +
      geom_text(aes(label = sample), size = 2) + 
      theme_bw() +
      ggtitle("Housekeeping genes") +
      theme_Publication2()  
    s <- ggplot(pca.data, aes(PC3, PC4)) +
      geom_point(size = 5, alpha = 0.5, aes_string(color = color_var, shape = shape_var)) +
      geom_text(aes(label = sample), size = 2) + 
      theme_bw() +
      ggtitle("Housekeeping genes") +
      theme_Publication2() 
    res <- ggarrange(p, q, r, s, ncol = 2, nrow = 2, common.legend = T)
  }
  
  # if topVar is present use that
  if(!is.null(topVar)){
    
    # use the top variable genes
    rv <- rowVars(voomData)
    select <- order(rv, decreasing=TRUE)[seq_len(topVar)]
    voomData <- voomData[select,]
    prData <- prcomp(voomData)
    pca.data <- prData$rotation
    pca.data <- data.frame(pca.data)[1:4]
    pca.data <- data.frame(pca.data, meta)
    r <- ggplot(pca.data, aes(PC1, PC2)) +
      geom_point(size = 5, alpha = 0.5, aes_string(color = color_var, shape = shape_var)) +
      geom_text(aes(label = sample), size = 2) + 
      theme_bw() +
      ggtitle("Most variable genes") +
      theme_Publication2()  
    s <- ggplot(pca.data, aes(PC3, PC4)) +
      geom_point(size = 5, alpha = 0.5, aes_string(color = color_var, shape = shape_var)) +
      geom_text(aes(label = sample), size = 2) + 
      theme_bw() +
      ggtitle("Most variable genes") +
      theme_Publication2() 
    res <- ggarrange(p, q, r, s, ncol = 2, nrow = 2, common.legend = T)
  }
  
  # print plot
  print(res)
  dev.off()
}
