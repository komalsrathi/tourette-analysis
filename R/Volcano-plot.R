# Author: Komal S. Rathi
# Function: Volcano plot of limma output
# Date: 04/28/2020
library(scales)

# accessory for volcano plot
reverselog_trans <- function(base = exp(1)) {
  trans <- function(x) -log(x, base)
  inv <- function(x) base^(-x)
  trans_new(paste0("reverselog-", format(base)), trans, inv,
            log_breaks(base = base),
            domain = c(1e-100, Inf))
}

# volcano plot, takes in limma analysis
plotVolcano <- function(result, fname, yaxis = c('P.Value','adj.P.Val'), title = "Volcano Plot (DEGs)", lfcutoff = 0, pvalcutoff = 0.05) {

  result$yaxis <- result[,yaxis]
  if(!is.null(lfcutoff)){
    result$DEGAnnot <- 'Other'
    result$DEGAnnot[result$yaxis < pvalcutoff & result$logFC > lfcutoff] <- "Up"
    result$DEGAnnot[result$yaxis < pvalcutoff & result$logFC < -(lfcutoff)] <- "Down"
  }

  result <- result %>%
    group_by(DEGAnnot) %>%
    mutate(freq = n()) %>%
    ungroup() %>%
    mutate(DEGAnnot = paste0(DEGAnnot, ' (n = ', freq, ')'))

  p <- ggplot(result, aes(x = logFC, y = yaxis, color = DEGAnnot)) +
    geom_point() + scale_y_continuous(trans = reverselog_trans(10)) +
    ggtitle(title) + scale_colour_manual(values = c("darkblue","gray","darkred")) +
    theme_Publication() + labs(color = "DEG") + ylab(yaxis)
  ggsave(filename = fname, plot = p, device = "pdf", width = 7, height = 5)
}
