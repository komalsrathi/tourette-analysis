# merge fpkm output to excel file

add.fpkm <- function(x, fpkm_mat){
  x <- x %>%
    inner_join(fpkm_mat %>%
                 as.data.frame() %>%
                 rownames_to_column('gene_symbol'), by = "gene_symbol")
  return(x)
}
