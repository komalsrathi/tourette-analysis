# human
Rscript --vanilla R/create_gsea_input_files_human.R \
--count_matrix 'data/human/human_collapsed_counts_matrix.RData' \
--meta_file 'data/human/human-meta-data-filtered.txt' \
--output_dir 'results/human_filtered/gsea' \
--prefix 'human'

