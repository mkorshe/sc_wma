

library("data.table")   
library("dplyr")  
library("stringr")  
library("readr")  
dataset_combo <- 'stemi_v2_vs_stemi_v3'
dataset_combos <- c('5ds','ng_vs_stemi_v2','ng_vs_stemi_v3','ng_vs_v2_1m','ng_vs_v3_1m','stemi_v3_vs_v3_1m','v2_1m_vs_stemi_v2','v2_1m_vs_stemi_v3','v2_1m_vs_v3_1m','v3_1m_vs_stemi_v2','stemi_v2_vs_stemi_v3')

for(dataset_combo in dataset_combos[11]){
# dataset_combo <- 'EF'
dir_aggregated <- '/groups/umcg-franke-scrna/tmp01/projects/sc-eqtlgen-consortium-pipeline/ongoing/wg3-QTL-mapping/WMA_round_II/backup/WMA_round_II/monocyte/output/MTC_corrected/'

dir_with_WMA_results <- '/groups/umcg-franke-scrna/tmp01/projects/sc-eqtlgen-consortium-pipeline/ongoing/wg3-QTL-mapping/WMA_round_II/backup/WMA_round_II/monocyte/output/WMA/aggregated/'
setwd(dir_aggregated)

# extracting weights 

colnames_of_complete_output <- read.table('/groups/umcg-franke-scrna/tmp01/projects/sc-eqtlgen-consortium-pipeline/ongoing/wg3-QTL-mapping/WMA_round_II/backup/WMA_round_II/pseudobulk/output/WMA/colnames.tsv' ,sep='\t', header = T)

# # 
# singificant_coeqtls_readlines <- read_lines(paste0(dir_with_WMA_results, dataset_combo, ".tsv"), skip=0, n_max=10)
# 
# eqtl <-  as.data.frame(str_split_fixed(singificant_coeqtls_readlines, "\t", 463)) 
# colnames(eqtl) <- colnames_of_complete_output$x

'%!in%' <- function(x,y)!('%in%'(x,y))

zeros_cols <-    colnames_of_complete_output$x[grep('zeros', colnames_of_complete_output$x)]
colnames_of_complete_output <- colnames_of_complete_output[colnames_of_complete_output$x %!in%zeros_cols, ]

# A_bulk_sample_beta_param <- read.table('/groups/umcg-franke-scrna/tmp01/projects/sc-eqtlgen-consortium-pipeline/ongoing/wg3-QTL-mapping/LIMIX/run_limix/ipsc/readlines_ipcs/updated_top_results_filtered/A/Stats_qtl_results_polished.txt',sep='\t',header = T)

v2_1m_stats <- read.table('/groups/umcg-franke-scrna/tmp01/projects/sc-eqtlgen-consortium-pipeline/ongoing/wg3-QTL-mapping/WMA_round_II/backup/WMA_round_II/pseudobulk/input/v2_beta_param_eqtlgen_zscore_fdr.tsv',sep='\t',header = T)


# 
 colnames_of_complete_output[colnames_of_complete_output=='ZW_weight_SE'] <- 'ZW_SE'
colnames_of_complete_output[colnames_of_complete_output=='pvalue_ZW_SE'] <- 'pvalue_SE'


weight_columns <-colnames_of_complete_output[grep('pvalue_',colnames_of_complete_output)]
weight_columns <- gsub("pvalue_", "", weight_columns)


####### recalculating TP eGenes and eQTLs measurements using the new 

TP_table_genes <- read.table(paste0(dir_aggregated, dataset_combo, '_TP_gene_effects_redefined.tsv'), sep = '\t', header=F)
TP_table_eqtl <- read.table(paste0(dir_aggregated, dataset_combo, '_TP_eqtl_effects_redefined.tsv'), sep = '\t', header=F)


sum_tab <- data.frame(c('eGenes','eQTls','cor_all','cor_sign','TP_eqtls','TN_eqtls','FP_eqtls','FN_eqtls','Number_of_eqtls','OE_sign_both','OE_sign_ZW', 'total_number_of_tests','TP_gene','TN_gene','FP_gene','FN_gene'))


colnames(sum_tab) <- 'param'
sum_tab_A <- sum_tab
'%!in%' <- function(x,y)!('%in%'(x,y))


summarizing_wma_pbmc <- function(id, qtl_tab_wma,qtlsA_sig,top_QTL,sum_tab,TP_table_genes,TP_table_eqtl){
  # top_QTL_intermediate <- table_qtl %>% group_by(feature_idv2_1m) %>%  arrange(pval)  %>% distinct(feature_idv2_1m, .keep_all= TRUE)
  
  vector_of_results <- c(length(which(top_QTL$global_FDR<0.1)),#   eGenes
                         length(unique(qtlsA_sig$snp_gene)), # eQTLs 
                         
                         cor.test(qtl_tab_wma$OverallZScore.eqtlgenv2_1m, qtl_tab_wma$zscore)[[4]][[1]], #cor.all
                         cor.test(qtlsA_sig$OverallZScore.eqtlgenv2_1m, qtlsA_sig$zscore)[[4]][[1]], #cor.sign
                         
                         length(unique(qtlsA_sig[qtlsA_sig$snp_gene %in% TP_table_eqtl$V1 ,]$snp_gene)),  # TP
                         length(unique(qtl_tab_wma[qtl_tab_wma$snp_gene %!in% TP_table_eqtl$V1  & qtl_tab_wma$snp_gene %!in% qtlsA_sig$snp_gene  ,]$snp_gene)),  # TN
                         length(unique(qtlsA_sig[qtlsA_sig$snp_gene %!in% TP_table_eqtl$V1 ,]$snp_gene)),  # FP
                         length(unique(qtl_tab_wma[qtl_tab_wma$snp_gene %in% TP_table_eqtl$V1  & qtl_tab_wma$snp_gene %!in% qtlsA_sig$snp_gene  ,]$snp_gene)),  # FN
                         
                         dim(qtl_tab_wma)[1],
                         length(unique(qtl_tab_wma[qtl_tab_wma$significant ==1& qtl_tab_wma$significant_A  ==1 & qtl_tab_wma$oposite_effect  ==1,]$snp_gene)),
                         length(unique(qtl_tab_wma[qtl_tab_wma$significant ==1& qtl_tab_wma$significant_A  =='0' & qtl_tab_wma$oposite_effect  ==1,]$snp_gene)),
                         dim(qtl_tab_wma)[1],
                         
                         # Gene 
                         length(unique(qtlsA_sig[qtlsA_sig$feature_idv2_1m %in% TP_table_genes$V1 ,]$feature_idv2_1m)),  # TP
                         length(unique(qtl_tab_wma[qtl_tab_wma$feature_idv2_1m %!in% TP_table_genes$V1  & qtl_tab_wma$feature_idv2_1m %!in% qtlsA_sig$feature_idv2_1m  ,]$feature_idv2_1m)),  # TN
                         length(unique(qtlsA_sig[qtlsA_sig$feature_idv2_1m %!in% TP_table_genes$V1 ,]$feature_idv2_1m)),  # FP
                         length(unique(qtl_tab_wma[qtl_tab_wma$feature_idv2_1m %in% TP_table_genes$V1  & qtl_tab_wma$feature_idv2_1m %!in% qtlsA_sig$feature_idv2_1m  ,]$feature_idv2_1m))  # FN
  )
  
  sum_tab[,id] <- vector_of_results
  
  return(sum_tab)
}

# 

for(i in weight_columns){
  
  zcol <-  print(which(colnames_of_complete_output == paste0('ZW_',i)))
  pval_col <-   print(which(colnames_of_complete_output == paste0('pvalue_',i)))
  
  qtl_tab_wma <- fread(paste0(dir_with_WMA_results, dataset_combo, ".tsv"),select = c(1, zcol,pval_col))
  colnames(qtl_tab_wma) <- c('snp_gene','zscore','pval')
  qtl_tab_wma <- qtl_tab_wma[!(duplicated(qtl_tab_wma)), ]
  qtl_tab_wma <- qtl_tab_wma[complete.cases(qtl_tab_wma), ]
  
  # qtl_tab_wma <- merge(qtl_tab_wma, B_stats, by='snp_gene')
  qtl_tab_wma <- merge(qtl_tab_wma, v2_1m_stats,  by='snp_gene')
  #qtl_tab_wma <- merge(qtl_tab_wma, B_stats,  by='snp_gene')
  
  qtl_tab_wma <- qtl_tab_wma[!(duplicated(qtl_tab_wma)), ]
  
  top_QTL <- qtl_tab_wma %>% group_by(feature_idv2_1m) %>%  arrange(pval)  %>% distinct(feature_idv2_1m, .keep_all= TRUE)
  
  # by sample A 
  
  top_QTL$empirical_feature_pvalue <- top_QTL$beta_paramv2_1m * top_QTL$pval
  top_QTL$empirical_feature_pvalue[top_QTL$empirical_feature_pvalue>1] <-  1
  #top_QTL["global_FDR"] =  qvalue::qvalue(top_QTL$empirical_feature_pvalue)$qvalues
  top_QTL["global_FDR"] =  qvalue::qvalue(top_QTL$empirical_feature_pvalue, pi0 = 1)$qvalues
  
  min(top_QTL$global_FDR)
  print(length(which(top_QTL$global_FDR<0.1)))
  
  
  qtl_tab_wma$empirical_feature_pvalue <- qtl_tab_wma$beta_paramv2_1m * qtl_tab_wma$pval
  qtl_tab_wma$empirical_feature_pvalue[qtl_tab_wma$empirical_feature_pvalue>1] <-  1
  qtl_tab_wma["global_FDR"] =  qvalue::qvalue(qtl_tab_wma$empirical_feature_pvalue)$qvalues
  
  qtlsA_sig = qtl_tab_wma[which(qtl_tab_wma$empirical_feature_pvalue<=max(top_QTL$empirical_feature_pvalue[which(top_QTL$global_FDR<0.1)])),]
  qtl_tab_wma$significant <-  ifelse(qtl_tab_wma$snp_gene%in% qtlsA_sig$snp_gene, 1,0)
  qtl_tab_wma$oposite_effect<- ifelse(qtl_tab_wma$OverallZScore.eqtlgenv2_1m * qtl_tab_wma$zscore > 0, 0,1)
  # 
  sum_tab_A <-   summarizing_wma_pbmc(id=i, qtl_tab_wma,qtlsA_sig, top_QTL,sum_tab_A,TP_table_genes,TP_table_eqtl)
  
  write.table(sum_tab_A, paste(dir_aggregated, 'output_sumamry/','sign_by_B_BH_bulk_',dataset_combo, ".tsv", sep=''), sep = '\t', col.names = T, row.names = F)
  
  # write.table(qtlsA_sig, paste(dir_aggregated, '/intermediate/','sign_by_BH_bulk_',dataset_combo, '_',i,".tsv", sep=''), sep = '\t', col.names = T, row.names = F)
  print(paste(dataset_combo, i, 'DONE'))
}

write.table(sum_tab_A, paste(dir_aggregated, '/output_sumamry/','1_sign_by_B_BH_bulk_',dataset_combo, ".tsv", sep=''), sep = '\t', col.names = T, row.names = F)

}