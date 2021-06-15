# Author: Komal S. Rathi
# Function: Annotate and filter limma output
# Date: 04/02/2020

annotate.limma <- function(x, foldchange, annot) {
  
  # select specific columns
  x <- x %>%
    rownames_to_column("gene_symbol") %>%
    dplyr::select(gene_symbol, matches("_logFC"), P.Value, adj.P.Val, AveExpr) %>%
    inner_join(annot %>% dplyr::select(-c(gene_id)), by = 'gene_symbol')

  # filter by logFC
  varname1 <- paste0(gsub('_logFC','',colnames(x)[2]),'_DEGAnnot')
  varname2 <- paste0(gsub('_logFC','',colnames(x)[3]),'_DEGAnnot')
  x <- x %>%
    mutate(!!varname1 := ifelse(!!as.symbol(colnames(x)[2]) > foldchange, "Up", "Down"),
           !!varname2 := ifelse(!!as.symbol(colnames(x)[3]) > foldchange, "Up", "Down"))
  
  # return
  return(x)
}