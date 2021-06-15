# Collapse RNA-seq
Rscript R/collapse_matrix.R \
--input data/human/human_hg38.RData \
--annot data/gencode.v23.annotation.txt \
--prefix human \
--outdir data/human

# Remove filter by variance
# 399_2 vs 370-1
Rscript R/diffexpr_human.R \
--counts_matrix data/human/human_collapsed_counts_matrix.RData \
--fpkm_matrix data/human/human_collapsed_fpkm_matrix.RData \
--meta_file data/human/human-meta-data-filtered.txt \
--gene_list data/gene-lists/housekeeping-genes-human.txt \
--var_filter FALSE \
--type '399_2, 370_1' \
--col 'label' \
--fc 0 \
--plx 3 \
--prefix 'diffexpr-399_2-vs-370_1' \
--excel TRUE \
--text TRUE \
--outdir results/human_filtered/diffexpr-399_2-vs-370_1_novarfilter

# 399_2 vs 370_5
Rscript R/diffexpr_human.R \
--counts_matrix data/human/human_collapsed_counts_matrix.RData \
--fpkm_matrix data/human/human_collapsed_fpkm_matrix.RData \
--meta_file data/human/human-meta-data-filtered.txt \
--gene_list data/gene-lists/housekeeping-genes-human.txt \
--var_filter FALSE \
--type '399_2, 370_5' \
--col 'label' \
--fc 0 \
--plx 3 \
--prefix 'diffexpr-399_2-vs-370_5' \
--excel TRUE \
--text TRUE \
--outdir results/human_filtered/diffexpr-399_2-vs-370_5_novarfilter

# 370_1 vs 370_5
Rscript R/diffexpr_human.R \
--counts_matrix data/human/human_collapsed_counts_matrix.RData \
--fpkm_matrix data/human/human_collapsed_fpkm_matrix.RData \
--meta_file data/human/human-meta-data-filtered.txt \
--gene_list data/gene-lists/housekeeping-genes-human.txt \
--var_filter FALSE \
--type '370_1, 370_5' \
--col 'label' \
--fc 0 \
--plx 3 \
--prefix 'diffexpr-370_1-vs-370_5' \
--excel TRUE \
--text TRUE \
--outdir results/human_filtered/diffexpr-370_1-vs-370_5_novarfilter

# remove variance filter and do one comparison tumor vs normals
Rscript R/diffexpr_human.R \
--counts_matrix data/human/human_collapsed_counts_matrix.RData \
--fpkm_matrix data/human/human_collapsed_fpkm_matrix.RData \
--meta_file data/human/human-meta-data-filtered.txt \
--gene_list data/gene-lists/housekeeping-genes-human.txt \
--var_filter FALSE \
--type 'Treat, Control' \
--col 'treat' \
--fc 0 \
--plx 3 \
--prefix 'diffexpr_treat_vs_control' \
--excel TRUE \
--text TRUE \
--outdir results/human_filtered/diffexpr_treat_vs_control_novarfilter
