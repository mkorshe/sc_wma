# iPSC processing 

library(dplyr)
di_with_aggregated_significant_results <- ('/groups/umcg-franke-scrna/tmp01/projects/sc-eqtlgen-consortium-pipeline/ongoing/wg3-QTL-mapping/LIMIX/run_limix/ipsc/readlines_ipcs/updated_top_results_filtered/wma_output_recalculated_SS/aggregated/')


# cd_ss <- read.table(paste0(di_with_aggregated_significant_results, '/intermediate/CDNumber_of_donors_by_bulk.tsv'),sep='\t',header=T)


all_files_in_dir <- list.files('/groups/umcg-franke-scrna/tmp01/projects/sc-eqtlgen-consortium-pipeline/ongoing/wg3-QTL-mapping/LIMIX/run_limix/ipsc/readlines_ipcs/updated_top_results_filtered/wma_output_recalculated_SS/aggregated/intermediate')

all_files_in_dir_by_bulk <- all_files_in_dir[grep( 'bulk.tsv',all_files_in_dir)]


dataset_combos <- c('EF', "DF","CD","CE")

dataset_combo <- 'EF'
file_wma <- 'EFTotal_reads_per_cohortvariance_log.cell_by_bulk.tsv'

for(dataset_combo in dataset_combos){
  all_files_in_dir_by_bulk_dataset <- all_files_in_dir_by_bulk[grep( dataset_combo,all_files_in_dir_by_bulk)]
  
  significant_effects <- c()
  significant_effects_eqtl <- c()
  
  for(file_wma in all_files_in_dir_by_bulk_dataset){
    wma_output_per_weight <- read.table(paste0(di_with_aggregated_significant_results, '/intermediate/',file_wma),sep='\t',header=T)
    # wma_output_per_weight_bulk_and_wma <- wma_output_per_weight[wma_output_per_weight$significant_in_A ==1 & wma_output_per_weight$significant_in_wma ==1, ]
    
    print(file_wma)

    wma_output_per_weight_bulk_and_wma <- wma_output_per_weight[ wma_output_per_weight$significant_in_wma ==1, ]
    
    
    min(top_QTL$global_FDR)
    print(length(which(top_QTL$global_FDR<0.1)))
    
    
    significant_effects <- c(significant_effects, wma_output_per_weight_bulk_and_wma$feature_id_B)
    significant_effects <- unique(significant_effects)
    
    significant_effects_eqtl <- c(significant_effects_eqtl, wma_output_per_weight_bulk_and_wma$snp_gene)
    significant_effects_eqtl <- unique(significant_effects_eqtl)
    
  }
  write.table(significant_effects, paste0(di_with_aggregated_significant_results,'UPDATE_MTC/',dataset_combo, '_TP_gene_effects_redefined.tsv'),sep='\t', row.names = F, col.names = F, quote = F)
  
  write.table(significant_effects_eqtl, paste0(di_with_aggregated_significant_results,'UPDATE_MTC/',dataset_combo, '_TP_eqtl_effects_redefined.tsv'),sep='\t', row.names = F, col.names = F, quote = F)
  
  
}

# Extra checks: 
# length(significant_effects)
# significant_effects_u <- unique(significant_effects)
# length(significant_effects_u)
# 
# # 
# A_bulk_sample_beta_param <- read.table('/groups/umcg-franke-scrna/tmp01/projects/sc-eqtlgen-consortium-pipeline/ongoing/wg3-QTL-mapping/LIMIX/run_limix/ipsc/readlines_ipcs/updated_top_results_filtered/A/Stats_qtl_results_polished.txt',sep='\t',header = T)
# 
#   length(intersect(A_bulk_sample_beta_param[A_bulk_sample_beta_param$significant_A ==1,]$snp_gene,significant_effects_u ))