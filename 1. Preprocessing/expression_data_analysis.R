library(dplyr)

eqtltop <- read.table('/groups/umcg-franke-scrna/tmp01/projects/sc-eqtlgen-consortium-pipeline/ongoing/wg3-QTL-mapping/LIMIX/run_limix/limix_like_eQTLGen_top_key_coulmns.txt', sep='\t', header=T)
colnames(eqtltop) <- paste(colnames(eqtltop) , '.eqtlgen',sep='')

ng_path <- '/groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/LIMIX/wma_limix/old_conda_output/ng/'
stemiv2_path <- '/groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/LIMIX/wma_limix/old_conda_output/stemiv2/'
stemiv3_path <- '/groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/LIMIX/wma_limix/old_conda_output/stemiv3/'
onemil_v2_path <- '/groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/LIMIX/cfdr_limix/output/v2/bulk/'
onemil_v3_path <- '/groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/LIMIX/cfdr_limix/output/v3/bulk/'

stemi_genomes <-'/groups/umcg-franke-scrna/tmp01/releases/blokland-2020/v1/genotype/stemi_all_merged_plink1.bim' 
ng_genomes <-'/groups/umcg-franke-scrna/tmp01/releases/wijst-2018-hg19/v1/genotype/plink.bim' 
m1_genomes <-'/groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/genotype/1M_LL/1M_LL_genotypes.bim' 

bim_files_path <- '/groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/LIMIX/cfdr_limix/bim_files_annotated/'
path <- ng_genomes
m1 <- read.table(paste0(bim_files_path, '1m_snp_ids_annotated.txt'), sep='\t', header=T)
ng <- read.table(paste0(bim_files_path, 'ng_snp_ids_annotated.txt'), sep='\t', header=T)
stemi <- read.table(paste0(bim_files_path, 'stemi_snp_ids_annotated.txt'), sep='\t', header=T)


# ng_path - file with limix_oupput  
change_tab <- function(ng_path, mark, bim_file){
  
  tab_full <- read.table(paste0(ng_path, 'qtl_results_all.txt'), sep='\t', header=T)
  tab <- merge(tab_full, bim_file[,c("V2",'snp','V5','V6' )], by.x ='snp_id', by.y='V2')
  
  tab$snp_gene <- paste(tab$snp, tab$feature_id, sep='_')
  
  tab <- tab[tab$snp_gene %in% eqtltop$snp_gene.eqtlgen, ]
  #colnames(tab)[colnames(tab) == 'V5'] <- 'assessed_allele'
  colnames(tab)[colnames(tab) == 'V5'] <- 'alt_allele'
  dim(tab)
  
  # tab <- subset(tab, select=c('feature_id', 'snp_id','effect_allele',"alt_allele", 'p_value','beta', 'beta_se','empirical_feature_p_value','n_samples'))
  tab$OverallZScore <- tab$beta/tab$beta_se
  tab <- tab[!duplicated(tab$snp_gene),]
  
  tab <- merge(tab, eqtltop[,c('snp_gene.eqtlgen', 'p_value.eqtlgen','assessed_allele.eqtlgen','alt_allele.eqtlgen','OverallZScore.eqtlgen')], by.x='snp_gene', by.y='snp_gene.eqtlgen')
  dim(tab)
  
  #tab_full2 <- tab_full %>% group_by(feature_id) %>%  arrange(p_value)  %>% distinct(feature_id, .keep_all= TRUE)
  
  tab$corrected_zscore <- 'a'
  tab <- tab[tab$assessed_allele == tab$assessed_allele.eqtlgen | tab$assessed_allele == tab$alt_allele.eqtlgen,  ]
  tab <- tab[tab$alt_allele == tab$assessed_allele.eqtlgen | tab$alt_allele == tab$alt_allele.eqtlgen,  ]
  
  tab <- transform(tab,  corrected_zscore = ifelse(tab$assessed_allele == tab$assessed_allele.eqtlgen, OverallZScore, OverallZScore*(-1)))
  dim(tab)
  head(tab[,c('snp_gene','assessed_allele.eqtlgen','alt_allele.eqtlgen', 'assessed_allele', 'alt_allele','corrected_zscore')])
  
  
  
  write.table(tab, paste0(ng_path, 'limix_like_eQTLGen_top_Genes_loc.txt'), sep='\t', quote=F, row.names = F)
  colnames(tab) <- paste(mark, colnames(tab), sep="_")
  colnames(tab)[colnames(tab) == paste(mark, 'snp_gene', sep="_")] <- 'snp_gene'
  
  return(tab)
  
}

ng_processed <- change_tab(ng_path,'ng', ng)

bim_file = m1
ng_path = onemil_v2_path
change_tab_2 <- function(ng_path, mark, bim_file){
  tab_full <- read.table(paste0(ng_path, 'qtl_results_all.txt'), sep='\t', header=T)
  tab <- merge(tab_full, bim_file, by.x='snp_id', by.y='V2')
  tab$snp_id <- tab_full$snp_id.y
  
  tab$snp_gene <- paste(tab$snp, tab$feature_id, sep='_')
  
  tab <- tab[tab$snp_gene %in% eqtltop$snp_gene.eqtlgen, ]
  dim(tab)
  #colnames(tab)[colnames(tab) == 'V5'] <- 'assessed_allele'
  colnames(tab)[colnames(tab) == 'V5'] <- 'alt_allele'
  dim(tab)
  
  # tab <- subset(tab, select=c('feature_id', 'snp_id','effect_allele',"alt_allele", 'p_value','beta', 'beta_se','empirical_feature_p_value','n_samples'))
  tab$OverallZScore <- tab$beta/tab$beta_se
  tab <- tab[!duplicated(tab$snp_gene),]
  
  tab <- merge(tab, eqtltop[,c('snp_gene.eqtlgen', 'p_value.eqtlgen','assessed_allele.eqtlgen','alt_allele.eqtlgen','OverallZScore.eqtlgen')], by.x='snp_gene', by.y='snp_gene.eqtlgen')
  dim(tab)
  #tab_full2 <- tab_full %>% group_by(feature_id) %>%  arrange(p_value)  %>% distinct(feature_id, .keep_all= TRUE)
  
  tab$corrected_zscore <- 'a'
  tab <- tab[tab$assessed_allele == tab$assessed_allele.eqtlgen | tab$assessed_allele == tab$alt_allele.eqtlgen,  ]
  tab <- tab[tab$alt_allele == tab$assessed_allele.eqtlgen | tab$alt_allele == tab$alt_allele.eqtlgen,  ]
  
  tab <- transform(tab,  corrected_zscore = ifelse(tab$assessed_allele == tab$assessed_allele.eqtlgen, OverallZScore, OverallZScore*(-1)))
  
  
  write.table(tab, paste0(ng_path, 'limix_like_eQTLGen_top_Genes_loc.txt'), sep='\t', quote=F, row.names = F)
  # colnames(tab) <- paste(mark, colnames(tab), sep="_")
  # colnames(tab)[colnames(tab) == paste(mark, 'snp_gene', sep="_")] <- 'snp_gene'
  # 
  return(tab)
  
}

stemi_v2_processed <- change_tab_2(stemiv2_path,'stemiv2', stemi)
stemi_v3_processed <- change_tab_2(stemiv3_path,'stemiv3',stemi)

m1_v2_processed <- change_tab_2(onemil_v2_path,'v2_1M',m1)
m1_v3_processed <- change_tab_2(onemil_v3_path,'v3_1M', m1)
dim(m1_v3_processed)


#### adding expression info ####


path_for_expression_stat_output <- '/groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/LIMIX/cfdr_limix/limix_qtl_and_expr_stat/update_220407/'

file_SD_donor<- read.table('/groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/GRN_reconstruction/Seurat_features/donor_SD/summary_mean_SD_for_all_ds.tsv',sep='\t',header=T)
colnames(file_SD_donor)[colnames(file_SD_donor) == 'v2'] <- 'v2_1m'
colnames(file_SD_donor)[colnames(file_SD_donor) == 'v3'] <- 'v3_1m'

path_donor <- '/groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/GRN_reconstruction/Seurat_features/aggragated_by_donor_pseudobulk/'

path_cell <- '/groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/GRN_reconstruction/Seurat_features/cell_pseudobulk/'
dataset_id <- 'v2_1m'
path_with_qtls <- onemil_v2_path


expr_stats <- function(dataset_id, path_with_qtls, path_donor, path_cell, file_SD_donor, path_for_expression_stat_output){
  
  donor_tab <- read.table(paste(path_donor, dataset_id, 'stat_donor_aggregated.tsv', sep=''), sep='\t', header=T)
  colnames(donor_tab)<- paste0(colnames(donor_tab), '.donor')
  colnames(donor_tab)[colnames(donor_tab) == 'gene.donor'] <- 'feature_id'
  
  cell_tab <- read.table(paste(path_cell, dataset_id, '.tsv', sep=''), sep='\t', header=T)
  colnames(cell_tab)<- paste0(colnames(cell_tab), '.cell')
  colnames(cell_tab)[colnames(cell_tab) == 'feature_id.cell'] <- 'feature_id'
  
  qtl_stat <- read.table(paste(path_with_qtls,  'limix_like_eQTLGen_top_Genes_loc.txt', sep=''), sep='\t', header=T)
  
  qtl_stat_merged <- merge(qtl_stat,donor_tab, by ='feature_id')
  qtl_stat_merged <- merge(qtl_stat_merged, cell_tab, by ='feature_id')
  qtl_stat_merged <- merge(qtl_stat_merged,file_SD_donor[,c('gene', dataset_id)], by.x ='feature_id',by.y='gene')
  colnames(qtl_stat_merged)[colnames(qtl_stat_merged) == dataset_id] <- 'avr_SD_donor'
  
  columns_to_remove <- c('feature_chromosome','feature_start', 'feature_end', 'n_e_samples','snp_chromosome', 'snp_position', 'call_rate', "V1", "V3", "V4", 'assessed_allele.1','feature_id.eqtlgen', 'snp_loc.eqtlgen', 'zeros_log.donor', 'hgnc.cell', 'median_log.cell', 'length.cell','gc.cell',  'length_log10.cell', 'length_log2.cell', 'zeros.donor')
  dif <- setdiff(colnames(qtl_stat_merged),columns_to_remove)
  qtl_stat_merged_clean <- subset(qtl_stat_merged, select =dif)
  
  
  write.table(qtl_stat_merged, paste(path_for_expression_stat_output,dataset_id, '_all.tsv', sep=''), sep='\t')
  write.table(qtl_stat_merged_clean, paste(path_for_expression_stat_output,dataset_id, '_clean.tsv', sep=''), sep='\t')
  
}


expr_stats(dataset_id='v2_1m', path_with_qtls=onemil_v2_path, path_donor, path_cell, file_SD_donor, path_for_expression_stat_output)

expr_stats(dataset_id='v3_1m', path_with_qtls=onemil_v3_path, path_donor, path_cell, file_SD_donor, path_for_expression_stat_output)

expr_stats(dataset_id='stemi_v2', path_with_qtls=stemiv2_path, path_donor, path_cell, file_SD_donor, path_for_expression_stat_output)

expr_stats(dataset_id='stemi_v3', path_with_qtls=stemiv3_path, path_donor, path_cell, file_SD_donor, path_for_expression_stat_output)

expr_stats(dataset_id='ng', path_with_qtls=ng_path, path_donor, path_cell, file_SD_donor, path_for_expression_stat_output)


#### leaving overlap 

path_for_expression_stat_output_overlap <- '/groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/LIMIX/cfdr_limix/limix_qtl_and_expr_stat/update_220407/overlap/'

file_1m_v2 <- read.table(paste(path_for_expression_stat_output, dataset_id, '_clean.tsv', sep=''), sep='\t', header=T)

eqtls_overlap <- unique(file_1m_v2$snp_gene)
ds_ids <- c('v2_1m','v3_1m','stemi_v2','stemi_v3','ng')

for(ds_id in ds_ids){
  file_stat_qtl <- read.table(paste(path_for_expression_stat_output, ds_id, '_clean.tsv', sep=''), sep='\t', header=T)
  eqtls_tab <- unique(file_stat_qtl$snp_gene)
  length(eqtls_tab)
  eqtls_overlap <- intersect(eqtls_overlap, eqtls_tab)
}

eqtls_overlap <- unique(eqtls_overlap)
length(eqtls_overlap)

for(ds_id in ds_ids){
  file_stat_qtl <- read.table(paste(path_for_expression_stat_output, ds_id, '_clean.tsv', sep=''), sep='\t', header=T)
  file_stat_qtl <- file_stat_qtl[file_stat_qtl$snp_gene %in% eqtls_overlap,]
  write.table(file_stat_qtl, paste(path_for_expression_stat_output_overlap,ds_id, '.tsv', sep=''), sep='\t')
}
