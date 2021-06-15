library(tidyverse)
library(xlsx)

# with variance filter
diffexpr.399.370.1 <- read.delim('results/human/diffexpr-399_2-vs-370_1/summary/diffexpr-399_2-vs-370_1.txt', check.names = F)
diffexpr.399.370.5 <- read.delim('results/human/diffexpr-399_2-vs-370_5/summary/diffexpr-399_2-vs-370_5.txt', check.names = F)

merged.output  <- diffexpr.399.370.1 %>%
  mutate(`370_1-399_2_P.Value` = P.Value,
         `370_1-399_2_adj.P.Val` = adj.P.Val) %>%
  select(gene_symbol, `370_1-399_2_logFC`, `370_1-399_2_P.Value`, `370_1-399_2_adj.P.Val`, `370_1-399_2_DEGAnnot`) %>%
  full_join(diffexpr.399.370.5 %>%
              mutate(`370_5-399_2_P.Value` = P.Value,
                     `370_5-399_2_adj.P.Val` = adj.P.Val) %>%
               select(gene_symbol, `370_5-399_2_logFC`, `370_5-399_2_P.Value`, `370_5-399_2_adj.P.Val`, `370_5-399_2_DEGAnnot`),
             by = 'gene_symbol')
write.xlsx(merged.output, file = 'results/human/diffexpr_combined.xlsx', row.names = F)

# no variance filter
diffexpr.399.370.1 <- read.delim('results/human/diffexpr-399_2-vs-370_1_novarfilter/summary/diffexpr-399_2-vs-370_1.txt', check.names = F)
diffexpr.399.370.5 <- read.delim('results/human/diffexpr-399_2-vs-370_5_novarfilter/summary/diffexpr-399_2-vs-370_5.txt', check.names = F)

merged.output  <- diffexpr.399.370.1 %>%
  mutate(`370_1-399_2_P.Value` = P.Value,
         `370_1-399_2_adj.P.Val` = adj.P.Val) %>%
  select(gene_symbol, `370_1-399_2_logFC`, `370_1-399_2_P.Value`, `370_1-399_2_adj.P.Val`, `370_1-399_2_DEGAnnot`) %>%
  full_join(diffexpr.399.370.5 %>%
              mutate(`370_5-399_2_P.Value` = P.Value,
                     `370_5-399_2_adj.P.Val` = adj.P.Val) %>%
              select(gene_symbol, `370_5-399_2_logFC`, `370_5-399_2_P.Value`, `370_5-399_2_adj.P.Val`, `370_5-399_2_DEGAnnot`),
            by = 'gene_symbol')
write.xlsx(merged.output, file = 'results/human/diffexpr_combined_novarfilter.xlsx', row.names = F)
