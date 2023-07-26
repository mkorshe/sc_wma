# MTC and significant QTL extraction


library("data.table")   
library("dplyr")  
library("stringr")  
library("readr")  
library("qvalue")  
# cd /groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/LIMIX/WMA_round_II/scripts/MTC

dataset_combos <- c('5ds','ng_vs_stemi_v2','ng_vs_stemi_v3','ng_vs_v2_1m','ng_vs_v3_1m', 'stemi_v3_vs_v3_1m','v2_1m_vs_stemi_v2','v2_1m_vs_stemi_v3','v2_1m_vs_v3_1m','v3_1m_vs_stemi_v2')
dataset_combo_1 <- c('5ds','ng_vs_stemi_v2','ng_vs_stemi_v3','ng_vs_v2_1m','ng_vs_v3_1m')
dataset_combo_2 <- c('stemi_v3_vs_v3_1m','v2_1m_vs_stemi_v2','v2_1m_vs_stemi_v3','v2_1m_vs_v3_1m','v3_1m_vs_stemi_v2')


out_path <- ('/groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/LIMIX/WMA_round_II/monocyte/output/MTC_corrected/')
path_with_files <- '/groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/LIMIX/WMA_round_II/monocyte/output/WMA/aggregated/'
# setwd(out_path)
# 
# test  <- read.table('/groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/LIMIX/WMA_round_II/pseudobulk/output/WMA/5ds/group_2/6805.tsv', sep='\t', header=T) # file which contains information baout weighting identifiers
# 
# colnames(test) <- colnames_for_df


colnames_for_df  <- read.table('/groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/LIMIX/WMA_round_II/pseudobulk/output/WMA/colnames.tsv', sep='\t', header=T) # file which contains information baout weighting identifiers
colnames_for_df <- c(colnames_for_df$x)

colnames_of_complete_output <- colnames_for_df
# colnames_of_complete_output <- colnames_of_complete_output$x
# # 
# colnames_of_complete_output <- colnames_of_complete_output[-1]

weight_columns <-colnames_for_df[grep('pvalue_',colnames_for_df)]
weight_columns <- gsub("pvalue_", "", weight_columns)


# ref data  


v2_1m_stats <- read.table('/groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/LIMIX/WMA_round_II/pseudobulk/input/v2_beta_param_eqtlgen_zscore_fdr.tsv',sep='\t',header = T)
# v2_1m_stats <- v2_1m_stats[,c('snp_gene', 'feature_idv2_1m','OverallZScore.eqtlgenv2_1m','BonferroniP.eqtlgenv2_1m','beta_paramv2_1m')] # 268641
# write.table(v2_1m_stats,'/groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/LIMIX/WMA_round_II/pseudobulk/input/v2_beta_param_eqtlgen_zscore_fdr.tsv', sep = '\t', col.names = T, row.names = F, quote=F)

v2_1m_stats$significant_A <-   ifelse(v2_1m_stats$BonferroniP.eqtlgenv2_1m  < 0.05, 1,0)

sum_tab <- data.frame(c('eGenes','eQTLs','correlation of all effects','correlation of significant effects','OE in significant eQTLs significant in bulk', 'Number of eQTLs'))
colnames(sum_tab) <- 'param'
# sum_tab <- data.frame(c('eGenes','eQTLs','correlation of all effects','correlation of significant effects','TP','TN','FP','FN','OE in significant eQTLs','OE in significant eQTLs significant in bulk', 'Number of eQTLs','TP gene','TN gene','FP gene','FN gene', 'Precision','Sensitivity'))


'%!in%' <- function(x,y)!('%in%'(x,y))

sum_tab_v2_1m <- sum_tab

summarizing_wma_pbmc <- function(id, qtl_tab_wma,top_QTL,qtlsA_sig,sum_tab){
  vector_of_results <- c(length(which(top_QTL$global_FDR <0.1 )),#   eGenes
                         length(unique(qtlsA_sig$snp_gene)), # eQTLs 
                         
                         cor.test(qtl_tab_wma$OverallZScore.eqtlgenv2_1m, qtl_tab_wma$zscore)[[4]][[1]], #cor.all
                         cor.test(qtlsA_sig$OverallZScore.eqtlgenv2_1m, qtlsA_sig$zscore)[[4]][[1]], #cor.sign
                         
                         # 
                         
                         length(unique(qtl_tab_wma[qtl_tab_wma$global_FDR <0.1 & qtl_tab_wma$significant_A  =='1' & qtl_tab_wma$oposite_effect  ==1,]$snp_gene)),
                         dim(qtl_tab_wma)[1]
                         
                         
  )
  
  
  sum_tab[,id] <- vector_of_results
  
  return(sum_tab)
  
}

# bulk2 <-read.table(      paste0("/groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/LIMIX/WMA_round_II/pseudobulk/output/WMA/aggregated/",dataset_combo,".tsv"),sep='\t',header=T)


i <- weight_columns[4]
dataset_combo <- '5ds'
# dataset_combo <-'ng_vs_v2_1m'
snps_TP <- c()


#mtc_pbmc <- function(dataset_combo){
  
  path_to_file <- paste0(path_with_files,dataset_combo,".tsv")
  # 
  # singificant_coeqtls_readlines <- read_lines(path_to_file, skip=0, n_max=100)
  # # # 
  # eqtl <-  as.data.frame(str_split_fixed(singificant_coeqtls_readlines, "\t", 463)) 
  # # 
  # colnames(eqtl) <- colnames_of_complete_output
  # 
   for(i in  weight_columns[2:length(weight_columns)]){
  #   
     zcol <-  print(which(colnames_of_complete_output == paste0('ZW_',i)))
     pval_col <-   print(which(colnames_of_complete_output == paste0('pvalue_',i)))
  # 
  # qtl_tab_wma <- fread(path_to_file,select = c(1,zcol,pval_col))

    qtl_tab_wma <-as.data.frame(read_tsv(file = path_to_file,col_select =c(1,zcol,pval_col),col_names=F))
  
    colnames(qtl_tab_wma) <- c('snp_gene','zscore','pval')
    # qtl_tab_wma$snp_gene <- gsub("\\\\","", qtl_tab_wma$snp_gene ,ignore.case = TRUE)
    
    qtl_tab_wma <- qtl_tab_wma[!(duplicated(qtl_tab_wma)), ]
    qtl_tab_wma <- qtl_tab_wma[complete.cases(qtl_tab_wma), ]
    #eqtltop$snp_gene.eqtlgen
    qtl_tab_wma <- merge(qtl_tab_wma, v2_1m_stats, by='snp_gene')
    qtl_tab_wma <- qtl_tab_wma[!(duplicated(qtl_tab_wma)), ]
    
    top_QTL <- qtl_tab_wma %>% group_by(feature_idv2_1m) %>%  arrange(pval)  %>% distinct(feature_idv2_1m, .keep_all= TRUE)#10865
    
    top_QTL$empirical_feature_pvalue <- top_QTL$beta_paramv2_1m * top_QTL$pval
    top_QTL$empirical_feature_pvalue[top_QTL$empirical_feature_pvalue>1] <-  1
    top_QTL[,"global_FDR"] =  qvalue::qvalue(top_QTL$empirical_feature_pvalue)$qvalues
    
    top_QTL$significant_in_wma <-  ifelse(top_QTL$global_FDR <0.1, 1,0)
    
    min(top_QTL$global_FDR)
    print(length(which(top_QTL$global_FDR<0.1)))
    
    qtl_tab_wma$empirical_feature_pvalue <- qtl_tab_wma$beta_paramv2_1m * qtl_tab_wma$pval
    qtl_tab_wma$empirical_feature_pvalue[qtl_tab_wma$empirical_feature_pvalue>1] <-  1
    qtl_tab_wma[,"global_FDR"] =  qvalue::qvalue(qtl_tab_wma$empirical_feature_pvalue)$qvalues
    
    qtlsA_sig = qtl_tab_wma[which(qtl_tab_wma$empirical_feature_pvalue<=max(top_QTL$empirical_feature_pvalue[which(top_QTL$global_FDR < 0.1)])),]
    qtlsA_sig$significant_in_wma <-  ifelse(qtlsA_sig$global_FDR <0.1, 1,0)
    
    qtlsA_sig$significant_in_A <-  ifelse(qtlsA_sig$significant_in_wma ==1 & (qtlsA_sig$zscore*qtlsA_sig$OverallZScore.eqtlgenv2_1m) > 0 , 1,0)
    
    print(length(which(qtlsA_sig$global_FDR<0.1)))
    print(length(which(qtlsA_sig$significant_in_A==1)))
    
    qtl_tab_wma$oposite_effect<- ifelse(qtl_tab_wma$OverallZScore.eqtlgenv2_1m * qtl_tab_wma$zscore > 0, 0,1)
    qtl_tab_wma$significant_in_wma <-  ifelse(qtl_tab_wma$snp_gene %in%  qtlsA_sig$snp_gene, 1,0)
    sum_tab_v2_1m <- summarizing_wma_pbmc(id=i, qtl_tab_wma,top_QTL,qtlsA_sig,sum_tab_v2_1m)
    
    print(i)
    snps_TP <- unique(snps_TP,qtlsA_sig[qtlsA_sig$significant_in_A ==1,]$snp_gene )
    write.table(qtlsA_sig, paste(out_path,'/intermediate/',dataset_combo,i,"_by_bulk.tsv", sep=''), sep = '\t', col.names = T, row.names = F)
    write.table(sum_tab_v2_1m, paste(out_path,'/intermediate/','1_',dataset_combo,"_by_bulk.tsv", sep=''), sep = '\t', col.names = T, row.names = F)
    
  }
  
  i <- 'SE'
  
  qtl_tab_wma <-as.data.frame(read_tsv(file = path_to_file,col_select =c(1,114,115),col_names=F))
  
  colnames(qtl_tab_wma) <- c('snp_gene','zscore','pval')

  qtl_tab_wma <- qtl_tab_wma[!(duplicated(qtl_tab_wma)), ]
  qtl_tab_wma <- qtl_tab_wma[complete.cases(qtl_tab_wma), ]
  #eqtltop$snp_gene.eqtlgen
  qtl_tab_wma <- merge(qtl_tab_wma, v2_1m_stats, by='snp_gene')
  qtl_tab_wma <- qtl_tab_wma[!(duplicated(qtl_tab_wma)), ]
  top_QTL <- qtl_tab_wma %>% group_by(feature_idv2_1m) %>%  arrange(pval)  %>% distinct(feature_idv2_1m, .keep_all= TRUE)#10865
  
  top_QTL$empirical_feature_pvalue <- top_QTL$beta_paramv2_1m * top_QTL$pval
  top_QTL$empirical_feature_pvalue[top_QTL$empirical_feature_pvalue>1] <-  1
  top_QTL[,"global_FDR"] =  qvalue::qvalue(top_QTL$empirical_feature_pvalue)$qvalues
  top_QTL$significant_in_wma <-  ifelse(top_QTL$global_FDR <0.1, 1,0)
  
  min(top_QTL$global_FDR)
  print(length(which(top_QTL$global_FDR<0.1)))
  qtl_tab_wma$empirical_feature_pvalue <- qtl_tab_wma$beta_paramv2_1m * qtl_tab_wma$pval
  qtl_tab_wma$empirical_feature_pvalue[qtl_tab_wma$empirical_feature_pvalue>1] <-  1
  qtl_tab_wma[,"global_FDR"] =  qvalue::qvalue(qtl_tab_wma$empirical_feature_pvalue)$qvalues
  
  qtlsA_sig = qtl_tab_wma[which(qtl_tab_wma$empirical_feature_pvalue<=max(top_QTL$empirical_feature_pvalue[which(top_QTL$global_FDR < 0.1)])),]
  qtlsA_sig$significant_in_wma <-  ifelse(qtlsA_sig$global_FDR <0.1, 1,0)
  qtlsA_sig$significant_in_A <-  ifelse(qtlsA_sig$significant_in_wma ==1 & (qtlsA_sig$zscore*qtlsA_sig$OverallZScore.eqtlgenv2_1m) > 0 , 1,0)
  print(length(which(qtlsA_sig$global_FDR<0.1)))
  print(length(which(qtlsA_sig$significant_in_A==1)))
  qtl_tab_wma$oposite_effect<- ifelse(qtl_tab_wma$OverallZScore.eqtlgenv2_1m * qtl_tab_wma$zscore > 0, 0,1)
  qtl_tab_wma$significant_in_wma <-  ifelse(qtl_tab_wma$snp_gene %in%  qtlsA_sig$snp_gene, 1,0)
  sum_tab_v2_1m <- summarizing_wma_pbmc(id=i, qtl_tab_wma,top_QTL,qtlsA_sig,sum_tab_v2_1m)
  print(i)
  snps_TP <- unique(snps_TP,qtlsA_sig[qtlsA_sig$significant_in_A ==1,]$snp_gene )
  
  write.table(qtlsA_sig, paste(out_path,'/intermediate/',dataset_combo,i,"_by_bulk.tsv", sep=''), sep = '\t', col.names = T, row.names = F)
  write.table(sum_tab_v2_1m, paste(out_path,'/intermediate/','2_',dataset_combo,"_by_bulk.tsv", sep=''), sep = '\t', col.names = T, row.names = F)
  write.table(snps_TP, paste(out_path,'/intermediate/','TP_',dataset_combo,".tsv", sep=''), sep = '\t', col.names = T, row.names = F)
  print(paste(dataset_combo, 'DONE'))
#}

#mtc_pbmc(dataset_combo='5ds')


# 
# 
# for(dataset_combo in dataset_combo_2){
#   mtc_pbmc(dataset_combo)
# }
# 





