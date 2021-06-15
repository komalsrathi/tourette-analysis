# Author: Komal S. Rathi
# Differential expression (Human)
# Date: 09/02/2020

options(java.parameters = "- Xmx1024m")
library(optparse)
library(tidyverse)
library(limma)
library(edgeR)
library(sva)
library(readxl)
library(xlsx)
library(matrixStats)

option_list <- list(
  make_option(c("--counts_matrix"), type = "character",
              help = "RData object of counts"),
  make_option(c("--fpkm_matrix"), type = "character",
              help = "RData object of fpkm"),
  make_option(c("--meta_file"), type = "character",
              help = "Metadata file for samples (.tsv)"),
  make_option(c("--gene_list"), type = "character",
              help = "House keeping gene list for tSNE/PCA plots"),
  make_option(c("--var_filter"), type = "character", default = TRUE,
              help = "Variance filter: T or F"),
  make_option(c("--type"), type = "character",
              help = "Type of comparison"),
  make_option(c("--col"), type = "character",
              help = "column for comparison"),
  make_option(c("--fc"), type = "character", default = 0,
              help = "fold-change for output files"),
  make_option(c("--plx"), type = "character", default = 0,
              help = "perplexity for t-SNE"),
  make_option(c("--prefix"), type = "character",
              help = "Prefix for output files"),
  make_option(c("--excel"), type = "character",
              default = TRUE,
              help = "Write to excel file?"),
  make_option(c("--text"), type = "character",
              default = TRUE,
              help = "Write to txt file?"),
  make_option(c("--outdir"), type = "character",
              help = "Output directory path")
)

# parse parameters
opt <- parse_args(OptionParser(option_list = option_list))
counts_matrix <- opt$counts_matrix
fpkm_matrix <- opt$fpkm_matrix
meta_file <- opt$meta_file
gene_list <- opt$gene_list
var_filter <- opt$var_filter
type <- opt$type
col <- opt$col
fc <- opt$fc
plx <- opt$plx
prefix <- opt$prefix
excel <- opt$excel
text <- opt$text
outdir <- opt$outdir

# source functions
source('R/design_pairs.R')
source('R/annotate_limma.R')
source('R/pubTheme.R')
source('R/PCA-plot.R')
source('R/Volcano-plot.R')
source('R/filterExpr.R')
source('R/add_fpkm_expression.R')

# read data
load(counts_matrix)
load(fpkm_matrix)
expr.counts.mat <- expr.counts[[1]]
expr.counts.annot <- expr.counts[[2]]
expr.fpkm <- expr.fpkm[[1]]
meta <- read.delim(meta_file, stringsAsFactors = F)
gene_list <- read.delim(gene_list,  stringsAsFactors = F)
type <- trimws(strsplit(type,",")[[1]]) 
plx <- as.numeric(plx)
fc <- as.numeric(fc)

# make output directories
pca.out <- file.path(outdir, 'pca')
dir.create(pca.out, showWarnings = F, recursive = TRUE)
tsne.out <- file.path(outdir, 'tsne')
dir.create(tsne.out, showWarnings = F, recursive = TRUE)
volcano.out <- file.path(outdir, 'volcano')
dir.create(volcano.out, showWarnings = F, recursive = TRUE)
summary.out <- file.path(outdir, 'summary')
dir.create(summary.out, showWarnings = F, recursive = TRUE)

# Analysis:
# 1. All comparisons between cell lines

# generalized function
diff.expr <- function(expr, expr.fpkm, meta, annot, type = NULL, 
                      gene_list,
                      var = var, fc = 0, plx = 7, 
                      fname = "plot_name", var_filter = TRUE,
                      write_to_excel = TRUE,  write_to_text = TRUE){
  
  # filter meta file by type
  meta <- meta %>%
    mutate(tmp = !!as.name(var)) %>%
    filter(tmp %in% type) %>%
    dplyr::select(-c(tmp)) %>%
    arrange(sample)
  
  # filter gene expression matrix by samples
  rownames(meta) <- meta$sample
  expr <- expr[,rownames(meta)]
  
  # filter gene expression for low expression
  print(var_filter)
  expr <- filterExpr(expr.counts.mat = expr, design = NULL, group = meta[,var], var.filter = var_filter)
  
  if(identical(rownames(meta), colnames(expr))) {
    print("Proceed")
  } else {
    break
  }
  
  # subset fpkm to count matrix samples
  expr.fpkm <- expr.fpkm[,colnames(expr)]
  
  # create design
  var <- factor(meta[,var])
  design <- model.matrix(~0+var)
  colnames(design) <- levels(var)
  rownames(design) <- meta$sample
  print(dim(design))
  
  # voom normalize
  y <- DGEList(counts = as.matrix(expr))
  y <- calcNormFactors(y)
  v <- voom(counts = y, design = design, plot = FALSE)
  voomData <- v$E
  
  # tsne and pca plot of voom normalized data
  pca.fname <- paste0(fname, '.pdf')
  pca.plot(voomData = voomData, topVar = 500, meta = meta, fname = file.path(pca.out, pca.fname), color_var = 'label', shape_var = 'treat')
  tsne.plot(voomData = voomData, topVar = 500, meta = meta, fname = file.path(tsne.out, pca.fname), plx = plx, color_var = 'label', shape_var = 'treat')
  
  # all levels of design
  all.pairs <- design.pairs(levels = colnames(design))
  
  # differential expression
  fit <- lmFit(voomData, design)
  fit2 <- contrasts.fit(fit, all.pairs)
  fit2 <- eBayes(fit2)
  
  # create list of all combinations
  t <- ncol(fit2$coefficients)
  newList <- list()
  for(i in 1:t){
    comp <- colnames(fit2$coefficients)[i] # comparison
    comp <- gsub(' ','_',comp)
    outputLimma <- topTable(fit2, coef = i, number = Inf)
    
    # plot volcano
    volc.fname <- paste0(fname, '-volcano.pdf')
    plotVolcano(result = outputLimma, fname = file.path(volcano.out, volc.fname), yaxis = "P.Value", lfcutoff = 0, pvalcutoff = 0.05)
    
    if(nrow(outputLimma) > 0){
      outputLimma[,paste0(comp,'_logFC')] <- outputLimma$logFC
      rev.comp <- paste(rev(strsplit(comp, "-")[[1]]), collapse = "-")
      outputLimma[,paste0(rev.comp,'_logFC')] <- -1*(outputLimma$logFC)
      
      # annotate and filter output
      outputLimma <- annotate.limma(x = outputLimma, foldchange = fc, annot)
      
      # add fpkm output
      outputLimma <- add.fpkm(x = outputLimma, fpkm_mat = expr.fpkm)

      print(paste0("Output dimensions:", dim(outputLimma)))
    } else {
      print("No significant difference")
    }
    newList[[i]] <- outputLimma
  }
  
  # write to excel file
  if(write_to_excel == TRUE){
    print("Writing output to excel...")
    xls.fname <- file.path(summary.out, paste0(fname,'.xlsx'))
    for (i in 1:length(newList)){
      if(nrow(newList[[i]]) > 1){
        write.xlsx(x = newList[[i]], file = xls.fname, row.names = FALSE, sheetName = "all_comparisons", append=T)
        gc()
      }
    }
  }
  
  # write to text file
  if(write_to_text == TRUE){
    print("Writing output to text...")
    txt.fname <- file.path(summary.out, paste0(fname, '.txt'))
    for(i in 1:length(newList)){
      if(nrow(newList[[i]]) > 1){
        tmp <- newList[[i]] %>%
          filter(adj.P.Val < 0.05)
        if(nrow(tmp) > 0) {
          write.table(tmp, file = txt.fname, quote = F, row.names = F, sep = "\t")
        }
      }
    }
  }
}

# call function
diff.expr(expr = expr.counts.mat, expr.fpkm = expr.fpkm, 
          meta = meta, gene_list = gene_list$genes, 
          annot = expr.counts.annot, 
          type = type, var = col, fc = fc, plx = plx, var_filter = var_filter,
          fname = prefix, write_to_excel = excel, write_to_text = text)
