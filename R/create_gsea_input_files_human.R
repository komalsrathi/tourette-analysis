# Author: Komal S. Rathi
# Function: Create GSEA input files

# libraries
library(tidyverse)
library(reshape2)
library(biomaRt)
library(optparse)

option_list <- list(
  make_option(c("--count_matrix"), type = "character",
              help = "Count matrix (.RData)"),
  make_option(c("--meta_file"), type = "character",
              help = "Metadata file"),
  make_option(c("--output_dir"), type = "character",
              help = "Output directory"),
  make_option(c("--prefix"), type = "character",
              help = "Prefix for output files")
)
# parse parameters
opt <- parse_args(OptionParser(option_list = option_list))
count_matrix <- opt$count_matrix
meta_file <- opt$meta_file
output_dir <- opt$output_dir
prefix <- opt$prefix
dir.create(output_dir, showWarnings = F, recursive = T)

# load count matrix
load(count_matrix)
expr.counts.mat <- expr.counts[[1]]
expr.counts.mat <- expr.counts.mat[apply(expr.counts.mat!=0, 1, all),]

# metadata
meta_file <- read.delim(meta_file, stringsAsFactors = F)

# function
gsea.input <- function(counts_collapsed, meta, groups, gct_file, cls_file) {
  
  # add group to meta file
  meta <- meta %>%
    filter(label %in% groups) %>%
    mutate(tmp = sample) %>%
    arrange(label) %>%
    column_to_rownames('tmp') 
  
  # order matrix
  counts_collapsed <- counts_collapsed[,rownames(meta)]
  
  # print dimensions
  print(groups)
  print(dim(meta))
  print(dim(counts_collapsed))
  print(colnames(counts_collapsed))
  
  # gct file
  gct <- counts_collapsed
  add <- data.frame(NAME = c("#1.2", nrow(gct), "NAME"), 
                    Description = c('', ncol(gct), "Description"))
  total.cols <- ncol(gct) + 2
  add[,3:total.cols] <- ''
  colnames(add)[3:total.cols] <- colnames(gct)
  add[3, 3:total.cols] <- colnames(gct)
  annot <- data.frame(NAME = rownames(gct), Description = 'na')
  annot <- merge(annot, gct, by.x = 'NAME', by.y = 'row.names')
  add <- rbind(add, annot)
  write.table(add, file = gct_file, quote = F, sep = "\t", col.names = F, row.names = F)
  
  # phenotype file
  groups <- unique(meta$label)
  ngroups <- length(groups)
  ph <- matrix(nrow = 3, ncol = ncol(gct))
  # first row
  ph[1,1] <- ncol(gct)
  ph[1,2] <- ngroups
  ph[1,3] <- 1
  # second row
  ph[2,1:ngroups] <- groups
  ph[2,1] <- paste0('# ', ph[2,1])
  # third row
  ph[3,] <- meta$label
  ph <- as.data.frame(ph)
  write.table(ph, file = cls_file, quote = F, sep = " ", na = "", col.names = F, row.names = F)
}

# 399-2 vs 370-1
gsea.input(counts_collapsed = expr.counts.mat, 
           meta = meta_file, 
           groups = c("399_2", "370_1"), 
           gct_file = file.path(output_dir, paste0(prefix, '_gsea_370_1.gct')),
           cls_file = file.path(output_dir, paste0(prefix, '_gsea_370_1.cls')))

# 399-2 vs 370-5
gsea.input(counts_collapsed = expr.counts.mat, 
           meta = meta_file, 
           groups = c("399_2", "370_5"), 
           gct_file = file.path(output_dir, paste0(prefix, '_gsea_370_5.gct')),
           cls_file = file.path(output_dir, paste0(prefix, '_gsea_370_5.cls')))


# by group
gsea.input_bygroup <- function(counts_collapsed, meta, groups, gct_file, cls_file) {
  
  # add group to meta file
  meta <- meta %>%
    filter(treat %in% groups) %>%
    mutate(tmp = sample) %>%
    arrange(treat) %>%
    column_to_rownames('tmp') 
  
  # order matrix
  counts_collapsed <- counts_collapsed[,rownames(meta)]
  
  # print dimensions
  print(groups)
  print(dim(meta))
  print(dim(counts_collapsed))
  print(colnames(counts_collapsed))
  
  # gct file
  gct <- counts_collapsed
  add <- data.frame(NAME = c("#1.2", nrow(gct), "NAME"), 
                    Description = c('', ncol(gct), "Description"))
  total.cols <- ncol(gct) + 2
  add[,3:total.cols] <- ''
  colnames(add)[3:total.cols] <- colnames(gct)
  add[3, 3:total.cols] <- colnames(gct)
  annot <- data.frame(NAME = rownames(gct), Description = 'na')
  annot <- merge(annot, gct, by.x = 'NAME', by.y = 'row.names')
  add <- rbind(add, annot)
  write.table(add, file = gct_file, quote = F, sep = "\t", col.names = F, row.names = F)
  
  # phenotype file
  groups <- unique(meta$treat)
  ngroups <- length(groups)
  ph <- matrix(nrow = 3, ncol = ncol(gct))
  # first row
  ph[1,1] <- ncol(gct)
  ph[1,2] <- ngroups
  ph[1,3] <- 1
  # second row
  ph[2,1:ngroups] <- groups
  ph[2,1] <- paste0('# ', ph[2,1])
  # third row
  ph[3,] <- meta$treat
  ph <- as.data.frame(ph)
  write.table(ph, file = cls_file, quote = F, sep = " ", na = "", col.names = F, row.names = F)
}

# control vs treat
gsea.input_bygroup(counts_collapsed = expr.counts.mat, 
           meta = meta_file, 
           groups = c("Control", "Treat"), 
           gct_file = file.path(output_dir, paste0(prefix, '_gsea_control_vs_treat.gct')),
           cls_file = file.path(output_dir, paste0(prefix, '_gsea_control_vs_treat.cls')))
