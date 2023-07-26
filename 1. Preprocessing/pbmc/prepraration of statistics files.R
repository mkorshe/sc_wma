bulk=read.table('/groups/umcg-franke-scrna/tmp01/projects/sc-eqtlgen-consortium-pipeline/ongoing/wg3-QTL-mapping/WMA_round_II/pseudobulk/input/Input_data_5ds.tsv',sep='\t',header=T)

library(dplyr)

bulk$abs_zscore <- abs(bulk$OverallZScore.eqtlgenng)
bulk_gene <- bulk %>% group_by(feature_idng) %>%  arrange(abs_zscore) %>% distinct(feature_idng, .keep_all= TRUE)

write.table(bulk_gene, '/groups/umcg-franke-scrna/tmp01/projects/sc-eqtlgen-consortium-pipeline/ongoing/wg3-QTL-mapping/WMA_round_II/pseudobulk/input/gene_Input_data_5ds.tsv',sep='\t',col.names = T, row.names = F, quote = F)

mono=read.table('/groups/umcg-franke-scrna/tmp01/projects/sc-eqtlgen-consortium-pipeline/ongoing/wg3-QTL-mapping/WMA_round_II/monocyte/input/Input_data_5ds.tsv',sep='\t',header=T)


mono$abs_zscore <- abs(mono$OverallZScore.eqtlgenng)
mono_gene <- mono %>% group_by(feature_idng) %>%  arrange(abs_zscore) %>% distinct(feature_idng, .keep_all= TRUE)

write.table(mono_gene, '/groups/umcg-franke-scrna/tmp01/projects/sc-eqtlgen-consortium-pipeline/ongoing/wg3-QTL-mapping/WMA_round_II/monocyte/input/gene_Input_data_5ds.tsv',sep='\t',col.names = T, row.names = F, quote = F)
mono_gene$id <- 'mono'
bulk_gene$id <- 'bulk'

bulk_mono <- rbind(mono_gene,bulk_gene )
write.table(bulk_mono, '/groups/umcg-franke-scrna/tmp01/projects/sc-eqtlgen-consortium-pipeline/ongoing/wg3-QTL-mapping/WMA_round_II/gene_stats_bulk_mono.tsv',sep='\t',col.names = T, row.names = F, quote = F)


# expracting cell proportions 


library(Seurat)
library(Matrix)

library(tidyverse)
library(liana)
library(nichenetr)
library(Seurat)
library(ggrepel)
library(dplyr)
library(cowplot)
library(readr)
library(stringr)
library(cowplot)
library(matrixStats)


STEMI_sequrat <- '/groups/umcg-franke-scrna/tmp01/releases/blokland-2020/v1/seurat/cardio.integrated.20210301.rds'
NG_Seurat <- '/groups/umcg-franke-scrna/tmp01/releases/wijst-2018-hg19/v1/clustering/pilot3_seurat3_200420_sct_azimuth.rds'
v2_1m <- '/groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/seurat_preprocess_samples/objects/1M_v2_mediumQC_ctd_rnanormed_demuxids_20210905.rds'
v3_1m <- '/groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/seurat_preprocess_samples/objects/1M_v3_mediumQC_ctd_rnanormed_demuxids_20210905.rds'


setwd('/groups/umcg-franke-scrna/tmp01/projects/sc-eqtlgen-consortium-pipeline/ongoing/wg3-QTL-mapping/WMA_round_II/')

v2 <- readRDS(v2_1m)

v2 <- v2[, !is.na(v2@meta.data$timepoint)]
v2 <- v2[, !is.na(v2@meta.data$assignment)]
v2 <- v2[, !is.na(v2@meta.data$cell_type_lowerres)]
v2 <- v2[, v2@meta.data$timepoint=='UT']

# v3 <- readRDS(v3_1m)
# stemi <- readRDS(STEMI_sequrat)
# v2 <- readRDS(NG_Seurat)
# 

freq_cell_type <- as.data.frame(table(v2$cell_type_lowerres))
identity_ds <- 'v2_1m'
write.table(freq_cell_type, paste0(identity_ds, "_cell_proportions_lowres.tsv"),sep='\t', quote = F, row.names = F, col.names = T)



stemi <- readRDS(STEMI_sequrat)
stemi <- stemi[,stemi@meta.data$timepoint_ll =='t8w']

v2 <-stemi[,stemi@meta.data$chem=='V2']              
# v2 <- v2[, !is.na(v2@meta.data$timepoint)]
v2 <- v2[, !is.na(v2@meta.data$assignment)]
v2 <- v2[, !is.na(v2@meta.data$cell_type_lowerres)]
# v2 <- v2[, v2@meta.data$timepoint=='UT']

# v3 <- readRDS(v3_1m)
# stemi <- readRDS(STEMI_sequrat)
# v2 <- readRDS(NG_Seurat)
# 

freq_cell_type <- as.data.frame(table(v2$cell_type_lowerres))
identity_ds <- 'stemi_v2'
write.table(freq_cell_type, paste0(identity_ds, "_cell_proportions_lowres.tsv"),sep='\t', quote = F, row.names = F, col.names = T)


v2 <-stemi[,stemi@meta.data$chem=='V3']              
# v2 <- v2[, !is.na(v2@meta.data$timepoint)]
v2 <- v2[, !is.na(v2@meta.data$assignment)]
v2 <- v2[, !is.na(v2@meta.data$cell_type_lowerres)]
# v2 <- v2[, v2@meta.data$timepoint=='UT']

# v3 <- readRDS(v3_1m)
# stemi <- readRDS(STEMI_sequrat)
# v2 <- readRDS(NG_Seurat)
# 

freq_cell_type <- as.data.frame(table(v2$cell_type_lowerres))
identity_ds <- 'stemi_v3'
write.table(freq_cell_type, paste0(identity_ds, "_cell_proportions_lowres.tsv"),sep='\t', quote = F, row.names = F, col.names = T)


v2 <- readRDS(NG_Seurat)

# v2 <- v2[, !is.na(v2@meta.data$timepoint)]
v2 <- v2[, !is.na(v2@meta.data$snumber)]
v2 <- v2[, !is.na(v2@meta.data$cell_type)]
# v2 <- v2[, v2@meta.data$timepoint=='UT']

# v3 <- readRDS(v3_1m)
# stemi <- readRDS(STEMI_sequrat)
# v2 <- readRDS(NG_Seurat)
# 

freq_cell_type <- as.data.frame(table(v2$cell_type))
identity_ds <- 'ng'
write.table(freq_cell_type, paste0(identity_ds, "_cell_proportions_lowres.tsv"),sep='\t', quote = F, row.names = F, col.names = T)




# ipsc
/groups/umcg-franke-scrna/tmp01/projects/sc-eqtlgen-consortium-pipeline/ongoing/wg3-QTL-mapping/LIMIX/run_limix/ipsc/readlines_ipcs/updated_top_results_filtered

file:
/groups/umcg-franke-scrna/tmp01/projects/sc-eqtlgen-consortium-pipeline/ongoing/wg3-QTL-mapping/LIMIX/run_limix/ipsc/readlines_ipcs/updated_top_results_filtered/overlaped_eqtls/top_gene_features.tsv


setwd('/groups/umcg-franke-scrna/tmp01/users/umcg-mkorshevniuk/LIMIX/cfdr_limix/limix_qtl_and_expr_stat/all_genes_eqtlgen/paralel_alternative_220907/pairwise/pseudobulk_readlines/updated_strategy_20221021/input_for_visualization/')

