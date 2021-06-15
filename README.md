# Tourette Analysis

Authors: Komal S. Rathi (@komalsrathi)

## Introduction

Various steps in the analysis:

1. QC
2. Process using STAR-RSEM pipeline
2. Differential gene expression
3. Pathway analysis using GSEA
4. Figures

## Methods

18 Human Tourett RNA-sequencing Fastq files, consisting of 6 controls and 12 treatment samples, were processed using the STAR alignment tool and subsequently normalized using the RSEM package based upon the hg38 reference genome and the Gencode version 23 gene annotation.

Differential gene expression analysis was performed by comparing each gene in the control group vs the treatment group. The voom procedure was used to normalize the RSEM generated expected counts followed by differential expression testing using R package limma to obtain P-values and LogFC. Specifically, a total of 58581 genes were tested for differential expression between the control and treatment samples. Pathway enrichment was performed on control vs treatment samples using Gene Set Enrichment Analysis (GSEA) version 4.1.0 using a weighted scoring scheme and Hallmark and C2 CP genesets.


