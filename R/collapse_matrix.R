# Author: Komal S. Rathi
# Function: To create collapsed matrices of gene symbols and samples
# Date: 04/02/2020

library(optparse)
library(tidyverse)
library(data.table)

option_list <- list(
  make_option(c("--input"), type = "character",
              help = "RData object with merged RSEM"),
  make_option(c("--annot"), type = "character",
              help = "Annotation in tab delim format"),
  make_option(c("--prefix"), type = "character",
              help = "Prefix for output files"),
  make_option(c("--outdir"), type = "character",
              help = "Output directory path")
)

# parse parameters
opt <- parse_args(OptionParser(option_list = option_list))
input <- opt$input
annot <- opt$annot
prefix <- opt$prefix
outdir <- opt$outdir

# load RSEM data (untrimmed)
load(input)
# load('data/skeletal_muscle_mm10.RData')

# annotation
annot <- data.table::fread(annot)
# annot <- data.table::fread('data/gencode.vM17.annotation.txt')
annot <- annot %>%
  select(gene_id,  gene_symbol, biotype) %>%
  unique()

# merge annotation with expression
collapse.mat <- function(expr, geneAnnot){
  
  # add annotation
  expr <- geneAnnot %>%
    select(gene_id, gene_symbol)  %>%
    inner_join(expr, by = c('gene_id'))
  
  # collapse to gene symbols
  expr <- expr %>%
    dplyr::mutate(means = rowMeans(select(.,-gene_id, -gene_symbol))) %>%
    arrange(desc(means)) %>%
    distinct(gene_symbol, .keep_all = TRUE) %>% 
    dplyr::select(-c(means)) %>%
    column_to_rownames(var = "gene_symbol")
  
  # subset annotation
  geneAnnot <- geneAnnot %>%
    filter(gene_id %in% expr$gene_id)
  
  # remove gene_id from matrix
  expr <- expr %>%
    select(-c(gene_id))
  
  newList <- list(expr, geneAnnot)
  return(newList)
}

# create collapsed count matrix
expr.counts <- collapse.mat(expr = expr.counts, geneAnnot = annot)
expr.counts.annot <- expr.counts[[2]]
expr.counts.mat <- expr.counts[[1]]

# create collapsed fpkm matrix
expr.fpkm <- collapse.mat(expr = expr.fpkm, geneAnnot = annot)
expr.fpkm.annot <- expr.fpkm[[2]]
expr.fpkm.mat <- expr.fpkm[[1]]

file1 <- file.path(outdir, paste0(prefix,'_collapsed_counts_matrix.RData'))
file2 <- file.path(outdir, paste0(prefix,'_collapsed_fpkm_matrix.RData'))

save(expr.counts, file = file1)
save(expr.fpkm, file = file2)
