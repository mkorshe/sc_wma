
  library("data.table")   
  library("dplyr")  
  library("stringr")  
  library("readr")  
  dataset_combo <- 'DF'
  # dataset_combo <- 'EF'
  dir_aggregated <- '/groups/umcg-franke-scrna/tmp01/projects/sc-eqtlgen-consortium-pipeline/ongoing/wg3-QTL-mapping/LIMIX/run_limix/ipsc/readlines_ipcs/updated_top_results_filtered/wma_output_recalculated_SS/aggregated_2/'
  setwd(dir_aggregated)
  
  # extracting weights 
  
  colnames_of_complete_output <- read.table('colnames.tsv' ,sep='\t')
  
  # 
  
  A_bulk_sample_beta_param <- read.table('/groups/umcg-franke-scrna/tmp01/projects/sc-eqtlgen-consortium-pipeline/ongoing/wg3-QTL-mapping/LIMIX/run_limix/ipsc/readlines_ipcs/updated_top_results_filtered/A/Stats_qtl_results_polished.txt',sep='\t',header = T)
  
  '%!in%' <- function(x,y)!('%in%'(x,y))
  # 
  # colnames_of_complete_output[colnames_of_complete_output=='ZW_SE'] <- 'SE'
  colnames_of_complete_output[colnames_of_complete_output=='pvalue_ZW_SE'] <- 'pvalue_SE'
  
  
  weight_columns <-colnames_of_complete_output$x[grep('pvalue_',colnames_of_complete_output$x)]
  weight_columns <- gsub("pvalue_", "", weight_columns)


  for(i in weight_columns){
    
    zcol <-  print(which(colnames_of_complete_output == paste0('ZW_',i)))
    pval_col <-   print(which(colnames_of_complete_output == paste0('pvalue_',i)))
    
    # TP_table_genes <- read.table(paste0(dir_aggregated,'/UPDATE_MTC/', dataset_combo, '_TP_gene_effects_redefined.tsv'), sep = '\t', header=F)
    # TP_table_eqtl <- read.table(paste0(dir_aggregated,'/UPDATE_MTC/', dataset_combo, '_TP_eqtl_effects_redefined.tsv'), sep = '\t', header=F)
    # 
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
    # qtl_tab_wma$significant <-  ifelse(qtl_tab_wma$snp_gene%in% qtlsA_sig$snp_gene, 1,0)
    # qtl_tab_wma$oposite_effect<- ifelse(qtl_tab_wma$OveralZscore_A * qtl_tab_wma$zscore > 0, 0,1)
    # 
    #sum_tab_A <- summarizing_wma_ipsc_A(id=i, qtl_tab_wma,qtlsA_sig, top_QTL,sum_tab_A,TP_table_genes,TP_table_eqtl)
    
   # write.table(sum_tab_A, paste(dir_aggregated, '/UPDATE_MTC/intermediate/','sign_by_B_BH_bulk_',dataset_combo, ".tsv", sep=''), sep = '\t', col.names = T, row.names = F)
    
    write.table(qtlsA_sig, paste(dir_aggregated, '/intermediate/','sign_by_BH_bulk_',dataset_combo, '_',i,".tsv", sep=''), sep = '\t', col.names = T, row.names = F)
    print(paste(dataset_combo, i, 'DONE'))
  }
  
####### getting TP eGenes and eQTLs
  
  all_files_in_dir <- list.files('/groups/umcg-franke-scrna/tmp01/projects/sc-eqtlgen-consortium-pipeline/ongoing/wg3-QTL-mapping/LIMIX/run_limix/ipsc/readlines_ipcs/updated_top_results_filtered/wma_output_recalculated_SS/aggregated_2/intermediate')
  
  all_files_in_dir_by_bulk <- all_files_in_dir[grep( 'bulk',all_files_in_dir)]

  all_files_in_dir_by_bulk_dataset <- all_files_in_dir_by_bulk[grep( dataset_combo,all_files_in_dir_by_bulk)]
  
  significant_effects <- c()
  significant_effects_eqtl <- c()
  
  for(file_wma in all_files_in_dir_by_bulk_dataset){
    wma_output_per_weight <- read.table(paste0(dir_aggregated, '/intermediate/',file_wma),sep='\t',header=T)
    # wma_output_per_weight_bulk_and_wma <- wma_output_per_weight[wma_output_per_weight$significant_in_A ==1 & wma_output_per_weight$significant_in_wma ==1, ]
    print(file_wma)
  
    # wma_output_per_weight_bulk_and_wma <- wma_output_per_weight[ wma_output_per_weight$significant_in_wma ==1, ]
    min(top_QTL$global_FDR)
    print(length(which(top_QTL$global_FDR<0.1)))
    
    significant_effects <- c(significant_effects, wma_output_per_weight$feature_id_A)
    significant_effects <- unique(significant_effects)

    significant_effects_eqtl <- c(significant_effects_eqtl, wma_output_per_weight$snp_gene)
    significant_effects_eqtl <- unique(significant_effects_eqtl)
    
  }
  
  significant_effects <- significant_effects[significant_effects %in%A_bulk_sample_beta_param[A_bulk_sample_beta_param$significant_A ==1,]$feature_id_A ]
  
  significant_effects_eqtl <- significant_effects_eqtl[significant_effects_eqtl %in%A_bulk_sample_beta_param[A_bulk_sample_beta_param$significant_A ==1,]$snp_gene ]
  
  write.table(significant_effects, paste0(dir_aggregated,'MTC/',dataset_combo, '_TP_gene_effects_redefined.tsv'),sep='\t', row.names = F, col.names = F, quote = F)
  
  write.table(significant_effects_eqtl, paste0(dir_aggregated,'MTC/',dataset_combo, '_TP_eqtl_effects_redefined.tsv'),sep='\t', row.names = F, col.names = F, quote = F)
  
  
####### recalculating TP eGenes and eQTLs measurements using the new 
    
    TP_table_genes <- read.table(paste0(dir_aggregated,'/MTC/', dataset_combo, '_TP_gene_effects_redefined.tsv'), sep = '\t', header=F)
    TP_table_eqtl <- read.table(paste0(dir_aggregated,'/MTC/', dataset_combo, '_TP_eqtl_effects_redefined.tsv'), sep = '\t', header=F)
    
    
    sum_tab <- data.frame(c('eGenes','eQTls','cor_all','cor_sign','TP_eqtls','TN_eqtls','FP_eqtls','FN_eqtls','Number_of_eqtls','OE_sign_both','OE_sign_ZW', 'total_number_of_tests','TP_gene','TN_gene','FP_gene','FN_gene'))
    
    
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
                             
                             dim(qtl_tab_wma)[1],
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
    
    
    
    for(i in weight_columns){
      
      zcol <-  print(which(colnames_of_complete_output == paste0('ZW_',i)))
      pval_col <-   print(which(colnames_of_complete_output == paste0('pvalue_',i)))
      
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
      # 
      sum_tab_A <- summarizing_wma_ipsc_A(id=i, qtl_tab_wma,qtlsA_sig, top_QTL,sum_tab_A,TP_table_genes,TP_table_eqtl)
      
      write.table(sum_tab_A, paste(dir_aggregated, '/MTC/intermediate/','sign_by_B_BH_bulk_',dataset_combo, ".tsv", sep=''), sep = '\t', col.names = T, row.names = F)
      
      # write.table(qtlsA_sig, paste(dir_aggregated, '/intermediate/','sign_by_BH_bulk_',dataset_combo, '_',i,".tsv", sep=''), sep = '\t', col.names = T, row.names = F)
      print(paste(dataset_combo, i, 'DONE'))
    }
    
    write.table(sum_tab_A, paste(dir_aggregated, '/MTC/intermediate/','1_sign_by_B_BH_bulk_',dataset_combo, ".tsv", sep=''), sep = '\t', col.names = T, row.names = F)
    