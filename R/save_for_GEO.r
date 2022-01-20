library(tidyverse)

# counts
load('human_collapsed_counts_matrix.RData')
expr_counts <- expr.counts[[1]]
colnames(expr_counts)
dat <- read.delim('human-meta-data-filtered.txt')
dat$sample <- as.character(dat$sample)
expr_counts <- expr_counts %>%
  dplyr::select(dat$sample)
expr_counts <- expr_counts %>%
  rownames_to_column("Gene_Symbol")
readr::write_tsv(expr_counts, file = 'tourette_expected_counts_matrix.tsv')

# fpkm
load('human_collapsed_fpkm_matrix.RData')
expr_fpkm <- expr.fpkm[[1]]
colnames(expr_fpkm)
dat <- read.delim('human-meta-data-filtered.txt')
dat$sample <- as.character(dat$sample)
expr_fpkm <- expr_fpkm %>%
  dplyr::select(dat$sample)
expr_fpkm <- expr_fpkm %>%
  rownames_to_column("Gene_Symbol")
readr::write_tsv(expr_fpkm, file = 'tourette_fpkm_matrix.tsv')
