# MTC for iPSC analysis 

library("data.table")   
library("dplyr")  
library("stringr")  
library("readr")  
dataset_combo <- 'EF'
# dataset_combo <- 'EF'
dir_aggregated <- '/groups/umcg-franke-scrna/tmp01/projects/sc-eqtlgen-consortium-pipeline/ongoing/wg3-QTL-mapping/LIMIX/run_limix/ipsc/readlines_ipcs/updated_top_results_filtered/wma_output_recalculated_SS/aggregated/'
setwd(dir_aggregated)

# extracting weights 

colnames_of_complete_output <- read.table('colnames_ZW.tsv' ,sep='\t')

# 

A_bulk_sample_beta_param <- read.table('/groups/umcg-franke-scrna/tmp01/projects/sc-eqtlgen-consortium-pipeline/ongoing/wg3-QTL-mapping/LIMIX/run_limix/ipsc/readlines_ipcs/updated_top_results_filtered/A/Stats_qtl_results_polished.txt',sep='\t',header = T)

'%!in%' <- function(x,y)!('%in%'(x,y))
# 

weight_columns <-colnames_of_complete_output$V1[grep('pvalue_',colnames_of_complete_output$V1)]
weight_columns <- gsub("pvalue_", "", weight_columns)

# weight_columns <-weight_columns[-grep('zeros_log',weight_columns)]

# ref data  
# B_stats <- fread('/groups/umcg-franke-scrna/tmp01/projects/sc-eqtlgen-consortium-pipeline/ongoing/wg3-QTL-mapping/LIMIX/run_limix/ipsc/readlines_ipcs/updated_top_results_filtered/overlaped_eqtls/all_effects_with_stats.tsv',select = c('snp_gene', 'OveralZscore_B','feature_id_B','beta_param_B','significant_B'))


sum_tab <- data.frame(c('eGenes','eQTls','cor_all','cor_sign','TP_SS','TN_SS','FP_SS','FN_SS','OE_sign_both','OE_sign_ZW', 'total_number_of_tests','TP_SS_gene','TN_SS_gene','FP_SS_gene','FN_SS_gene'))


colnames(sum_tab) <- 'param'
sum_tab_B <- sum_tab
sum_tab_A <- sum_tab
'%!in%' <- function(x,y)!('%in%'(x,y))


summarizing_wma_ipsc_A <- function(id, qtl_tab_wma,qtlsA_sig,top_QTL,sum_tab,TP_table_genes,TP_table_eqtl){
  # top_QTL_intermediate <- table_qtl %>% group_by(feature_id_A) %>%  arrange(pval)  %>% distinct(feature_id_A, .keep_all= TRUE)
  
  vector_of_results <- c(length(which(top_QTL$global_FDR<0.1)),#   eGenes
                         length(unique(qtlsA_sig$snp_gene)), # eQTLs 
                         
                         cor.test(qtl_tab_wma$OveralZscore_A, qtl_tab_wma$zscore)[[4]][[1]], #cor.all
                         cor.test(qtlsA_sig$OveralZscore_A, qtlsA_sig$zscore)[[4]][[1]], #cor.sign
                         
                         length(unique(qtlsA_sig[qtlsA_sig$snp_gene %in% TP_table_eqtl$V1 ,]$snp_gene)),  # TP
                         length(unique(qtl_tab_wma[qtl_tab_wma$snp_gene %!in% TP_table_eqtl$V1  & qtl_tab_wma$snp_gene %!in% qtlsA_sig$snp_gene  ,]$snp_gene)),  # TN
                         length(unique(qtlsA_sig[qtlsA_sig$snp_gene %!in% TP_table_eqtl$V1 ,]$snp_gene)),  # FP
                         length(unique(qtl_tab_wma[qtl_tab_wma$snp_gene %in% TP_table_eqtl$V1  & qtl_tab_wma$snp_gene %!in% qtlsA_sig$snp_gene  ,]$snp_gene)),  # FN
                         
                         length(unique(qtl_tab_wma[qtl_tab_wma$significant ==1& qtl_tab_wma$significant_A  ==1 & qtl_tab_wma$oposite_effect  ==1,]$snp_gene)),
                         length(unique(qtl_tab_wma[qtl_tab_wma$significant ==1& qtl_tab_wma$significant_A  =='0' & qtl_tab_wma$oposite_effect  ==1,]$snp_gene)),
                         dim(qtl_tab_wma)[1],
                         
                         # Gene 
                         length(unique(qtlsA_sig[qtlsA_sig$feature_id_A %in% TP_table_genes$V1 ,]$feature_id_A)),  # TP
                         length(unique(qtl_tab_wma[qtl_tab_wma$feature_id_A %!in% TP_table_genes$V1  & qtl_tab_wma$feature_id_A %!in% qtlsA_sig$feature_id_A  ,]$feature_id_A)),  # TN
                         length(unique(qtlsA_sig[qtlsA_sig$feature_id_A %!in% TP_table_genes$V1 ,]$feature_id_A)),  # FP
                         length(unique(qtl_tab_wma[qtl_tab_wma$feature_id_A %in% TP_table_genes$V1  & qtl_tab_wma$feature_id_A %!in% qtlsA_sig$feature_id_A  ,]$feature_id_A))  # FN
  )
  
  
  sum_tab[,id] <- vector_of_results
  
  return(sum_tab)
}


for(i in weight_columns[-1]){
  
  zcol <-  print(which(colnames_of_complete_output == paste0('ZW_',i)))
  pval_col <-   print(which(colnames_of_complete_output == paste0('pvalue_',i)))
  
  TP_table_genes <- read.table(paste0(dir_aggregated,'/UPDATE_MTC/', dataset_combo, '_TP_gene_effects_redefined.tsv'), sep = '\t', header=F)
  TP_table_eqtl <- read.table(paste0(dir_aggregated,'/UPDATE_MTC/', dataset_combo, '_TP_eqtl_effects_redefined.tsv'), sep = '\t', header=F)
  
  # zcol <- i
  # pval_col <- zcol+286
  # qtl_tab_wma2 <- read_lines(paste0(dataset_combo, ".tsv"), n_max = 2)
  # qtl_tab_wma2_df <- as.data.frame(str_split_fixed(qtl_tab_wma2, "\t", 577)) 
  # 
  # 
  
  qtl_tab_wma <- fread(paste0(dataset_combo, ".tsv"),select = c(1, zcol,pval_col))
  colnames(qtl_tab_wma) <- c('snp_gene','zscore','pval')
  qtl_tab_wma <- qtl_tab_wma[!(duplicated(qtl_tab_wma)), ]
  qtl_tab_wma <- qtl_tab_wma[complete.cases(qtl_tab_wma), ]
  
  # qtl_tab_wma <- merge(qtl_tab_wma, B_stats, by='snp_gene')
  qtl_tab_wma <- merge(qtl_tab_wma, A_bulk_sample_beta_param,  by='snp_gene')
  #qtl_tab_wma <- merge(qtl_tab_wma, B_stats,  by='snp_gene')
  
  qtl_tab_wma <- qtl_tab_wma[!(duplicated(qtl_tab_wma)), ]
  
  top_QTL <- qtl_tab_wma %>% group_by(feature_id_A) %>%  arrange(pval)  %>% distinct(feature_id_A, .keep_all= TRUE)
  
  # by sample A 
  
  top_QTL$empirical_feature_pvalue <- top_QTL$beta_param_A * top_QTL$pval
  top_QTL$empirical_feature_pvalue[top_QTL$empirical_feature_pvalue>1] <-  1
  #top_QTL["global_FDR"] =  qvalue::qvalue(top_QTL$empirical_feature_pvalue)$qvalues
  top_QTL["global_FDR"] =  qvalue::qvalue(top_QTL$empirical_feature_pvalue, pi0 = 1)$qvalues
  
  min(top_QTL$global_FDR)
  print(length(which(top_QTL$global_FDR<0.1)))
  
  
  qtl_tab_wma$empirical_feature_pvalue <- qtl_tab_wma$beta_param_A * qtl_tab_wma$pval
  qtl_tab_wma$empirical_feature_pvalue[qtl_tab_wma$empirical_feature_pvalue>1] <-  1
  qtl_tab_wma["global_FDR"] =  qvalue::qvalue(qtl_tab_wma$empirical_feature_pvalue)$qvalues
  
  qtlsA_sig = qtl_tab_wma[which(qtl_tab_wma$empirical_feature_pvalue<=max(top_QTL$empirical_feature_pvalue[which(top_QTL$global_FDR<0.1)])),]
  qtl_tab_wma$significant <-  ifelse(qtl_tab_wma$snp_gene%in% qtlsA_sig$snp_gene, 1,0)
  qtl_tab_wma$oposite_effect<- ifelse(qtl_tab_wma$OveralZscore_A * qtl_tab_wma$zscore > 0, 0,1)
  
  sum_tab_A <- summarizing_wma_ipsc_A(id=i, qtl_tab_wma,qtlsA_sig, top_QTL,sum_tab_A,TP_table_genes,TP_table_eqtl)
  
  write.table(sum_tab_A, paste(dir_aggregated, '/UPDATE_MTC/intermediate/','sign_by_B_BH_bulk_',dataset_combo, ".tsv", sep=''), sep = '\t', col.names = T, row.names = F)
  
  write.table(qtlsA_sig, paste(dir_aggregated, '/UPDATE_MTC/intermediate/','sign_by_B_BH_bulk_',dataset_combo, i,".tsv", sep=''), sep = '\t', col.names = T, row.names = F)
  print(paste(dataset_combo, i, 'DONE'))
}

i ='ZW_weight_SE'

zcol <-  print(which(colnames_of_complete_output == paste0(i)))
pval_col <-   print(which(colnames_of_complete_output == paste0('pvalue_','ZW_SE')))


qtl_tab_wma <- fread(paste0(dataset_combo, ".tsv"),select = c(1, zcol,pval_col))
colnames(qtl_tab_wma) <- c('snp_gene','zscore','pval')
qtl_tab_wma <- qtl_tab_wma[!(duplicated(qtl_tab_wma)), ]
qtl_tab_wma <- qtl_tab_wma[complete.cases(qtl_tab_wma), ]

# qtl_tab_wma <- merge(qtl_tab_wma, B_stats, by='snp_gene')
qtl_tab_wma <- merge(qtl_tab_wma, A_bulk_sample_beta_param,  by='snp_gene')
#qtl_tab_wma <- merge(qtl_tab_wma, B_stats,  by='snp_gene')

qtl_tab_wma <- qtl_tab_wma[!(duplicated(qtl_tab_wma)), ]

top_QTL <- qtl_tab_wma %>% group_by(feature_id_A) %>%  arrange(pval)  %>% distinct(feature_id_A, .keep_all= TRUE)

# by sample A 

top_QTL$empirical_feature_pvalue <- top_QTL$beta_param_A * top_QTL$pval
top_QTL$empirical_feature_pvalue[top_QTL$empirical_feature_pvalue>1] <-  1
#top_QTL["global_FDR"] =  qvalue::qvalue(top_QTL$empirical_feature_pvalue)$qvalues
top_QTL[,"global_FDR"] =  qvalue::qvalue(top_QTL$empirical_feature_pvalue, pi0 = 1)$qvalues

min(top_QTL$global_FDR)
print(length(which(top_QTL$global_FDR<0.1)))


qtl_tab_wma$empirical_feature_pvalue <- qtl_tab_wma$beta_param_A * qtl_tab_wma$pval
qtl_tab_wma$empirical_feature_pvalue[qtl_tab_wma$empirical_feature_pvalue>1] <-  1
qtl_tab_wma[,"global_FDR"] =  qvalue::qvalue(qtl_tab_wma$empirical_feature_pvalue)$qvalues

qtlsA_sig = qtl_tab_wma[which(qtl_tab_wma$empirical_feature_pvalue<=max(top_QTL$empirical_feature_pvalue[which(top_QTL$global_FDR<0.1)])),]
qtl_tab_wma$significant <-  ifelse(qtl_tab_wma$snp_gene%in% qtlsA_sig$snp_gene, 1,0)
qtl_tab_wma$oposite_effect<- ifelse(qtl_tab_wma$OveralZscore_A * qtl_tab_wma$zscore > 0, 0,1)

sum_tab_A <- summarizing_wma_ipsc_A(id=i, qtl_tab_wma,qtlsA_sig, top_QTL,sum_tab_A,TP_table_genes,TP_table_eqtl)

# writing the results 

write.table(sum_tab_A, paste(dir_aggregated, '/UPDATE_MTC/intermediate/','sign_by_B_BH_bulk_',dataset_combo, ".tsv", sep=''), sep = '\t', col.names = T, row.names = F)

write.table(qtlsA_sig, paste(dir_aggregated, '/UPDATE_MTC/intermediate/','sign_by_B_BH_bulk_',dataset_combo, i,".tsv", sep=''), sep = '\t', col.names = T, row.names = F)

write.table(sum_tab_A, paste(dataset_combo, "_summary_A_bulk_BH.tsv", sep=''), sep = '\t', col.names = T, row.names = F)



# extracting se and ss info
dataset_combos <-c('EF','CE','CD','DF')
ss_weights <- c('SE','Number_of_donors', 'Total_reads_per_cohort','Total_number_of_cells','Counts_per_person','Counts_Reads_per_cell','Avr_number_of_cell_per_donor')





### reanalysing SE data 

dataset_combos <-c('EF','CE','CD','DF')
dataset_combo <- dataset_combos[1]


A_bulk_sample_beta_param <- read.table('/groups/umcg-franke-scrna/tmp01/projects/sc-eqtlgen-consortium-pipeline/ongoing/wg3-QTL-mapping/LIMIX/run_limix/ipsc/readlines_ipcs/updated_top_results_filtered/A/Stats_qtl_results_polished.txt',sep='\t',header = T)


for(dataset_combo in dataset_combos){
  dataset_df <- read.table( paste(dir_aggregated, '/UPDATE_MTC/intermediate/','sign_by_B_BH_bulk_',dataset_combo, ".tsv", sep=''), sep='\t',header=T)
  
  
  TP_table_genes <- read.table(paste0(dir_aggregated,'/UPDATE_MTC/', dataset_combo, '_TP_gene_effects_redefined.tsv'), sep = '\t', header=F)
  TP_table_eqtl <- read.table(paste0(dir_aggregated,'/UPDATE_MTC/', dataset_combo, '_TP_eqtl_effects_redefined.tsv'), sep = '\t', header=F)
  
  
  i ='ZW_weight_SE'
  
  zcol <-  print(which(colnames_of_complete_output == paste0(i)))
  pval_col <-   print(which(colnames_of_complete_output == paste0('pvalue_','ZW_SE')))
  
  
  qtl_tab_wma <- fread(paste0(dataset_combo, ".tsv"),select = c(1, zcol,pval_col))
  colnames(qtl_tab_wma) <- c('snp_gene','zscore','pval')
  qtl_tab_wma <- qtl_tab_wma[!(duplicated(qtl_tab_wma)), ]
  qtl_tab_wma <- qtl_tab_wma[complete.cases(qtl_tab_wma), ]
  
  # qtl_tab_wma <- merge(qtl_tab_wma, B_stats, by='snp_gene')
  qtl_tab_wma <- merge(qtl_tab_wma, A_bulk_sample_beta_param,  by='snp_gene')
  #qtl_tab_wma <- merge(qtl_tab_wma, B_stats,  by='snp_gene')
  
  qtl_tab_wma <- qtl_tab_wma[!(duplicated(qtl_tab_wma)), ]
  
  top_QTL <- qtl_tab_wma %>% group_by(feature_id_A) %>%  arrange(pval)  %>% distinct(feature_id_A, .keep_all= TRUE)
  
  # by sample A 
  
  top_QTL$empirical_feature_pvalue <- top_QTL$beta_param_A * top_QTL$pval
  top_QTL$empirical_feature_pvalue[top_QTL$empirical_feature_pvalue>1] <-  1
  #top_QTL["global_FDR"] =  qvalue::qvalue(top_QTL$empirical_feature_pvalue)$qvalues
  top_QTL[,"global_FDR"] =  qvalue::qvalue(top_QTL$empirical_feature_pvalue, pi0 = 1)$qvalues
  
  min(top_QTL$global_FDR)
  print(length(which(top_QTL$global_FDR<0.1)))
  
  
  qtl_tab_wma$empirical_feature_pvalue <- qtl_tab_wma$beta_param_A * qtl_tab_wma$pval
  qtl_tab_wma$empirical_feature_pvalue[qtl_tab_wma$empirical_feature_pvalue>1] <-  1
  qtl_tab_wma[,"global_FDR"] =  qvalue::qvalue(qtl_tab_wma$empirical_feature_pvalue)$qvalues
  
  qtlsA_sig = qtl_tab_wma[which(qtl_tab_wma$empirical_feature_pvalue<=max(top_QTL$empirical_feature_pvalue[which(top_QTL$global_FDR<0.1)])),]
  qtl_tab_wma$significant <-  ifelse(qtl_tab_wma$snp_gene%in% qtlsA_sig$snp_gene, 1,0)
  qtl_tab_wma$oposite_effect<- ifelse(qtl_tab_wma$OveralZscore_A * qtl_tab_wma$zscore > 0, 0,1)
  
  sum_tab <- data.frame(c('eGenes','eQTls','cor_all','cor_sign','TP_SS','TN_SS','FP_SS','FN_SS','OE_sign_both','OE_sign_ZW', 'total_number_of_tests','TP_SS_gene','TN_SS_gene','FP_SS_gene','FN_SS_gene'))
  
  
  colnames(sum_tab) <- 'param'
  sum_tab_B <- sum_tab
  sum_tab_A <- sum_tab
  
  sum_tab_A <- summarizing_wma_ipsc_A(id=i, qtl_tab_wma,qtlsA_sig, top_QTL,sum_tab_A,TP_table_genes,TP_table_eqtl)
  
  
  dataset_df$SE <- sum_tab_A$ZW_weight_SE
  
  write.table(dataset_df, paste(dir_aggregated, '/UPDATE_MTC/intermediate/','sign_by_B_BH_bulk_',dataset_combo, ".tsv", sep=''), sep = '\t', col.names = T, row.names = F)
  
  
}

