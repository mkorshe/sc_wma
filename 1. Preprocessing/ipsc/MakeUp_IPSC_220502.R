# making up iPSC files 

# Setting information about expression and making up expression tab 
 
path_with_ipsc_stat <- '/groups/umcg-franke-scrna/tmp01/projects/sc-eqtlgen-consortium-pipeline/ongoing/wg3-QTL-mapping/LIMIX/run_limix/ipsc/expression_stat/ipsc_stat_qtl_filtered/'

samples_ids <- c('B', "C","D", "E", "F")

path_for_output <- '/groups/umcg-franke-scrna/tmp01/projects/sc-eqtlgen-consortium-pipeline/ongoing/wg3-QTL-mapping/LIMIX/run_limix/ipsc/expression_stat/ipsc_stat_qtl_merged_220502/'

making_expr_stat <- function(samples_id,path_for_output ){
  table_qtl<- read.
}



# # Expression tables #
# # expr1 <- read.table('/groups/umcg-franke-scrna/tmp01/projects/sc-eqtlgen-consortium-pipeline/ongoing/wg3-QTL-mapping/LIMIX/run_limix/ipsc/exp/ss2.tsv', header = T, se
# p ='\t', fill =T)
# # #ss2_metadata <-  read.table('/Users/m/Downloads/wma/iPSC/expression/SingleCell_ENSEMBL75_featureCounts_counts.tsv/CellSampleMetadata.txt', header = F, sep ='\t', fil
# l =T)
# # expr2 <- read.table('/groups/umcg-franke-scrna/tmp01/projects/sc-eqtlgen-consortium-pipeline/ongoing/wg3-QTL-mapping/LIMIX/run_limix/ipsc/exp/10x.txt', header = T, se
# p ='\t', fill =T)
# # 
# # 


#/groups/umcg-franke-scrna/tmp01/projects/sc-eqtlgen-consortium-pipeline/ongoing/wg3-QTL-mapping/LIMIX/run_limix/ipsc/ipsc_29_11_2021

# Loading QTL Output
limix_qtl_ss2_87_raw <- read.table('/groups/umcg-franke-scrna/tmp01/projects/sc-eqtlgen-consortium-pipeline/ongoing/wg3-QTL-mapping/LIMIX/run_limix/ipsc/ipsc_29_11_2021/WeigthedMetaAnalysis/b_SS2_all_Samples_dS/qtl_results_all_filteredDs.txt.gz', header = T)

limix_qtl_ss2_67_raw <- read.table('/groups/umcg-franke-scrna/tmp01/projects/sc-eqtlgen-consortium-pipeline/ongoing/wg3-QTL-mapping/LIMIX/run_limix/ipsc/ipsc_29_11_2021/WeigthedMetaAnalysis/c_SS2_Other_Samples_dS/qtl_results_all_filteredDs.txt.gz', header = T)

limix_qtl_ss2_25_raw <- read.table('/groups/umcg-franke-scrna/tmp01/projects/sc-eqtlgen-consortium-pipeline/ongoing/wg3-QTL-mapping/LIMIX/run_limix/ipsc/ipsc_29_11_2021/WeigthedMetaAnalysis/d_SS2_26_Samples_dS/qtl_results_all_filteredDs.txt.gz', header = T)

limix_qtl_10x_25_raw <- read.table('/groups/umcg-franke-scrna/tmp01/projects/sc-eqtlgen-consortium-pipeline/ongoing/wg3-QTL-mapping/LIMIX/run_limix/ipsc/ipsc_29_11_2021/WeigthedMetaAnalysis/e_10X_26_Samples_dS/qtl_results_all_filtered.txt.gz', header = T)

limix_qtl_ss2_25_2_raw <- read.table('/groups/umcg-franke-scrna/tmp01/projects/sc-eqtlgen-consortium-pipeline/ongoing/wg3-QTL-mapping/LIMIX/run_limix/ipsc/ipsc_29_11_2021/WeigthedMetaAnalysis/f_SS2_26_SamplesR_dS/qtl_results_all_filteredDs.txt.gz', header = T)


# 

out_dir <- '/groups/umcg-franke-scrna/tmp01/projects/sc-eqtlgen-consortium-pipeline/ongoing/wg3-QTL-mapping/LIMIX/run_limix/ipsc/batch2/20220315/'


prepare_data_limix <- function(object_path, dataset_id){
    object <- read.table(object_path, sep = '\t', header=T )
    object$snp_gene <- paste(object$snp_id, object$feature_id, sep='_')
    object$OverallZScore <- object$beta/object$beta_se
    object <- object[object$snp_gene %in% limix_qtl_ss2_87_top_sub$snp_gene, ]
    object <- subset(object, select = c( 'snp_gene', "feature_id","p_value","snp_id", "beta_se","assessed_allele", "OverallZScore",'empirical_feature_p_value'))
    
     object <- object %>% group_by(feature_id) %>%  arrange(empirical_feature_p_value)
    multTestOut = multtest::mt.rawp2adjp(c("BH","Bonferroni"),rawp = object$empirical_feature_p_value)$adjp
    # # add BH and Bonferroni corrections
    object["BH"] = multTestOut[,2] 
    # object <- object[object$BH < 0.1, ]
    # object <- object %>% distinct(feature_id, .keep_all= TRUE)
    colnames(object) <- paste(colnames(object), dataset_id, sep="_")
    colnames(object)[colnames(object) ==paste('snp_gene',dataset_id, sep="_")] <- 'snp_gene'
  
  return(object)
}

limix_qtl_ss2_67 <- prepare_data_limix(limix_qtl_ss2_67_path, "SS2_67_C")
limix_qtl_ss2_25 <- prepare_data_limix(limix_qtl_ss2_25_path,'SS2_25_D')
limix_qtl_10x_25 <- prepare_data_limix(limix_qtl_10x_25_path,'10X_25_E')
limix_qtl_ss2_25_2 <- prepare_data_limix(limix_qtl_ss2_25_2_path, 'SS2_25_F')

limix_qtl_ss2_87_top_sub <- merge(limix_qtl_ss2_87_top_sub, limix_qtl_ss2_67, by='snp_gene')
limix_qtl_ss2_87_top_sub <- merge(limix_qtl_ss2_87_top_sub, limix_qtl_ss2_25, by='snp_gene')
limix_qtl_ss2_87_top_sub <- merge(limix_qtl_ss2_87_top_sub, limix_qtl_10x_25, by='snp_gene')
limix_qtl_ss2_87_top_sub <- merge(limix_qtl_ss2_87_top_sub, limix_qtl_ss2_25_2, by='snp_gene')


limix_qtl_ss2_87_top_sub_backup_merged <- limix_qtl_ss2_87_top_sub
table0 = limix_qtl_ss2_87_top_sub_backup

limix_qtl_ss2_87_top_sub <- transform(limix_qtl_ss2_87_top_sub,  OverallZScore_SS2_67_C = ifelse(limix_qtl_ss2_87_top_sub$assessed_allele_SS2_87_B == limix_qtl_ss2_87_top_sub$assessed_allele_SS2_67_C, OverallZScore_SS2_67_C, OverallZScore_SS2_67_C*(-1)))

limix_qtl_ss2_87_top_sub <- transform(limix_qtl_ss2_87_top_sub,  OverallZScore_SS2_25_D = ifelse(limix_qtl_ss2_87_top_sub$assessed_allele_SS2_87_B == limix_qtl_ss2_87_top_sub$assessed_allele_SS2_25_D, OverallZScore_SS2_25_D, OverallZScore_SS2_25_D*(-1)))

limix_qtl_ss2_87_top_sub <- transform(limix_qtl_ss2_87_top_sub,  OverallZScore_10X_25_E = ifelse(limix_qtl_ss2_87_top_sub$assessed_allele_SS2_87_B == limix_qtl_ss2_87_top_sub$assessed_allele_10X_25_E, OverallZScore_10X_25_E, OverallZScore_10X_25_E*(-1)))

limix_qtl_ss2_87_top_sub <- transform(limix_qtl_ss2_87_top_sub,  OverallZScore_SS2_25_F = ifelse(limix_qtl_ss2_87_top_sub$assessed_allele_SS2_87_B == limix_qtl_ss2_87_top_sub$assessed_allele_SS2_25_F, OverallZScore_SS2_25_F, OverallZScore_SS2_25_F*(-1)))

write.table(limix_qtl_ss2_87_top_sub, paste0(path_out, 'merged_input.tsv', sep=''), sep='\t', col.names = T)

# get BH number #
merged_dt <- limix_qtl_ss2_87_top_sub
BH_table <- data.frame(c('id', "BH"))
colnames(BH_table) <- 'value'
dataset_ids <- c('SS2_67_C','SS2_25_D','10X_25_E','SS2_25_F')
  for (dataset_id in dataset_ids){
    BH_numbers <- merged_dt[,paste0('BH_',dataset_id)]
    BH_n <- length(BH_numbers[BH_numbers <0.05])
    BH_vector <- c(dataset_id, BH_n)
    BH_table[,dataset_id] <- BH_vector
  }
BH_table

#limix_qtl_ss2_87_top_sub$r2_SS2_87_B <- zToCor(z=limix_qtl_ss2_87_top_sub$OverallZScore_SS2_87_B, df = 87)
limix_qtl_ss2_87_top_sub$r2_SS2_67_C <- zToCor(z=limix_qtl_ss2_87_top_sub$OverallZScore_SS2_67_C, df = 67)
limix_qtl_ss2_87_top_sub$r2_SS2_25_D <- zToCor(z=limix_qtl_ss2_87_top_sub$OverallZScore_SS2_25_D, df = 25)
limix_qtl_ss2_87_top_sub$r2_10X_25_E <- zToCor(z=limix_qtl_ss2_87_top_sub$OverallZScore_10X_25_E, df = 25)
limix_qtl_ss2_87_top_sub$r2_SS2_25_F <- zToCor(z=limix_qtl_ss2_87_top_sub$OverallZScore_SS2_25_F, df = 25)

sum_tab <- data.frame(c('BH_SS', "BH_SE",'BH_SS_rev', "BH_SE_rev", 'cor_SS', 'cor_SE', 'cor_SS_rev', 'cor_SE_rev'))

SS1=67
SS2=25

weighting_ids <- c('mean.donor','mean.cell','mean_log.donor','mean_log.cell','sd.donor','sd.cell','sd_log.donor','sd_log.cell','variance.donor','variance.cell', 'cv.donor','cv.cell', 'cv_log.donor','cv_log.cell',  'median.donor','median.cell','sd_prop')


WMA_processing_limix <- function(table0, table1,table2, sample_size1, sample_size2,mark, out_dir){
  
  
  merged_dt <- merge(table1,table2, by="snp_gene")
  #merged_dt <- as.matrix(merged_dt)
  merged_dt <- merge(merged_dt, table0, by="snp_gene" )
  
  head(merged_dt)
  dim(merged_dt)
  
  merged_dt <- merged_dt[order(merged_dt$p_value.x),]
  bh_x <- merged_dt %>% distinct(feature_id, .keep_all= TRUE)
  multTestOut_table1 = multtest::mt.rawp2adjp(c("BH","Bonferroni"),rawp = bh_x$p_value.x)$adjp
  bh_x$BH.x <- multTestOut_table1[,2]
  bh_x$Bonferroni.x <- multTestOut_table1[,3]
  
  bh_x_top <- bh_x[bh_x$BH.x <0.05,]
  
  merged_dt <- merged_dt[order(merged_dt$p_value.y),]
  bh_y <- merged_dt %>% distinct(feature_id, .keep_all= TRUE)
  multTestOut_table2 = multtest::mt.rawp2adjp(c("BH","Bonferroni"),rawp = bh_y$p_value.y)$adjp
  bh_y$BH.y <- multTestOut_table2[,2]
  bh_y$Bonferroni.y <- multTestOut_table2[,3]
  bh_y_top <- bh_y[bh_y$BH.y <0.05,]
  
  merged_dt <-merged_dt[order(merged_dt$p_value),]
  bh <- merged_dt %>% distinct(feature_id, .keep_all= TRUE)
  multTestOut_table3 = multtest::mt.rawp2adjp(c("BH","Bonferroni"),rawp = bh$p_value)$adjp
  bh$BH <- multTestOut_table3[,2]
  bh$Bonferroni <- multTestOut_table3[,3]
  
  bh_top <- bh[bh$BH <0.05,]
  
  
  ### WMA on sample size
  
  for(i in 1:nrow(merged_dt)) {
    eqtl <- merged_dt[i,]
    weighted_z_score_ss <- get_weighted_z_score(
      zScores=c(eqtl$OverallZScore.x, eqtl$OverallZScore.y),
      weights=c(1,1),
      sampleSizes = c(sample_size1, sample_size2))
    merged_dt[i,"ZW_SS"] <- weighted_z_score_ss                              
  }
  merged_dt$pvalue_ss <- convert.z.score(merged_dt$ZW_SS)
  merged_dt<- merged_dt[order(merged_dt$pvalue_ss),]
  
  ZW_SS <- merged_dt %>% distinct(feature_id, .keep_all= TRUE)
  multTestOut_ss <- multtest::mt.rawp2adjp(ZW_SS$pvalue_ss, proc = c("BH","Bonferroni"))$adjp
  ZW_SS$cfdr_ss_BH<- multTestOut_ss[,2]
  ZW_SS$cfdr_Bonferroni.ss <- multTestOut_ss[,3]
  
  ZW_SS_top <- ZW_SS[ZW_SS$cfdr_ss_BH <0.05,] 
  
  
  for(i in 1:nrow(merged_dt)) {
    eqtl <- merged_dt[i,]
    weighted_z_score <- get_weighted_z_score_SE(
      zScores=c(eqtl$OverallZScore.x, eqtl$OverallZScore.y),
      weights=c(1,1),
      SEs = c(eqtl$beta_se.x, eqtl$beta_se.y)
    )
    merged_dt$ZW_SE[i] <- weighted_z_score          
  }
  merged_dt$pvalue_se <- convert.z.score(merged_dt$ZW_SE)
  merged_dt <- merged_dt[order(merged_dt$pvalue_se),]
  ZW_SE <- merged_dt %>% distinct(feature_id, .keep_all= TRUE)
  
  multTestOut_se <- multtest::mt.rawp2adjp(ZW_SE$pvalue_se, proc = c("BH","Bonferroni"))$adjp
  ZW_SE$cfdr_se_BH<- multTestOut_se[,2]
  ZW_SE$cfdr_Bonferroni.se <- multTestOut_se[,3]
  
  ZW_SE_top <- ZW_SE[ZW_SE$cfdr_se_BH <0.05,] 
  
  #snps_in_the_table <- intersect(bh$snp_gene, intersect(bh.x$snp_gene,bh.y$snp_gene),intersect(ZW_SS$snp_gene,ZW_SE$snp_gene))
  snps_in_the_table <-c(bh_top$snp_gene, bh_x_top$snp_gene,bh_y_top$snp_gene,ZW_SS$snp_gene, ZW_SE$snp_gene)
  snps_in_the_table <- snps_in_the_table[!duplicated(snps_in_the_table)]
  length(snps_in_the_table)
  merged_dt <- merged_dt[merged_dt$snp_gene %in% snps_in_the_table,]
  
  summary_tab <- data.frame(c("BH_tab1","BH_tab2", "BH_tab","BH_SS","BH_SE", 
                              "corr_BH_tab1","corr_BH_tab2","corr_BH_SS","corr_BH_SE",
                              "OE_tab1", "OE_tab2", "OE_SS", "OE_SE" ),
                            c(dim(bh_x_top)[1], dim(bh_y_top)[1], dim(bh_top)[1],dim(ZW_SS_top)[1],dim(ZW_SE_top)[1],
                              
                              cor.test(bh_x_top$OverallZScore.x, bh_x_top$OverallZScore, method =  c("pearson"))[[4]],
                              cor.test(bh_y_top$OverallZScore.y, bh_y_top$OverallZScore, method =  c("pearson"))[[4]],
                              cor.test(ZW_SS_top$ZW_SS, ZW_SS_top$OverallZScore, method =  c("pearson"))[[4]],
                              cor.test(ZW_SE_top$ZW_SE, ZW_SE_top$OverallZScore, method =  c("pearson"))[[4]], 
                              
                              dim(bh_x_top[bh_x_top$OverallZScore.x/bh_x_top$OverallZScore <0,])[1],
                              dim(bh_y_top[bh_y_top$OverallZScore.x/bh_y_top$OverallZScore <0,])[1],
                              dim(ZW_SS_top[ZW_SS_top$OverallZScore.x/ZW_SS_top$OverallZScore <0,])[1],
                              dim(ZW_SE_top[ZW_SE_top$OverallZScore.x/ZW_SE_top$OverallZScore <0,])[1]
                            )
                            
                            
  )
  summary_tab$id <- paste(mark)
  colnames(summary_tab) <- NULL 
  
  name_for_output2 <- paste(out_dir, mark, "_summary.tsv", sep='')
  write.table(summary_tab, name_for_output2, col.names = T, quote=F, sep ='\t')
  
  
  colnames(bh_x_top) <- paste(colnames(bh_x_top), "_bh_x", sep="")
  colnames(bh_y_top) <- paste(colnames(bh_y_top), "_bh_y", sep="")
  colnames(bh_top) <- paste(colnames(bh_top), "_bh", sep="")
  colnames(ZW_SS) <- paste(colnames(ZW_SS), "_ZW_SS", sep="")
  colnames(ZW_SE) <- paste(colnames(ZW_SE), "_ZW_SE", sep="")
  
  merged_dt_final <- merge(merged_dt, bh_x_top, by.x="snp_gene",by.y="snp_gene_bh_x", all.x=T)
  merged_dt_final <- merge(merged_dt_final, bh_y_top, by.x="snp_gene",by.y="snp_gene_bh_y", all.x=T)
  merged_dt_final <- merge(merged_dt_final, bh_top,by.x="snp_gene",by.y="snp_gene_bh", all.x=T)
  merged_dt_final <- merge(merged_dt_final, ZW_SS,by.x="snp_gene",by.y="snp_gene_ZW_SS", all.x=T)
  merged_dt_final <- merge(merged_dt_final, ZW_SE,by.x="snp_gene",by.y="snp_gene_ZW_SE", all.x=T)
  
  
  
  name_for_output <- paste(out_dir, mark, ".tsv", sep='')
  write.table(merged_dt_final, name_for_output, col.names = T, quote=F, sep ='\t')
  
}

# 
# intersected_snps <- intersect(intersect(limix_qtl_v2$snp_gene,limix_qtl_v3$snp_gene),intersect(limix_qtl_v2s$snp_gene,eqtlgen$snp_gene))
# length(intersected_snps)
# 
# limix_qtl_v2_intersect <- limix_qtl_v2[limix_qtl_v2$snp_gene %in%intersected_snps, ]
# limix_qtl_v2s_intersect <- limix_qtl_v2s[limix_qtl_v2s$snp_gene %in%intersected_snps, ]
# limix_qtl_v3_intersect <- limix_qtl_v3[limix_qtl_v3$snp_gene %in%intersected_snps, ]
# eqtlgen_intersect <- eqtlgen[eqtlgen$snp_gene %in%intersected_snps, ]
# 
# table0 = eqtlgen
# table1 = limix_qtl_v2
# table2 = limix_qtl_v3
# table3 = limix_qtl_v2s

# 
# WMA_processing_limix(table0=eqtlgen, table1=limix_qtl_v2,table2=limix_qtl_v3, sample_size1=72, sample_size2=32,mark="v2_v3",out_dir)
# WMA_processing_limix(table0=eqtlgen, table1=limix_qtl_v2s,table2=limix_qtl_v3, sample_size1=72, sample_size2=32,mark="v2s_v3",out_dir)
# 
# #
# WMA_processing_limix(table0=eqtlgen_intersect, table1=limix_qtl_v2_intersect,table2=limix_qtl_v3_intersect, sample_size1=72, sample_size2=32,mark="v2_v3_int",out_dir)
# WMA_processing_limix(table0=eqtlgen_intersect, table1=limix_qtl_v2s_intersect,table2=limix_qtl_v3_intersect, sample_size1=72, sample_size2=32,mark="v2s_v3_int",out_dir)
# 

#### Grid search ####
#weights <- seq(0,2,by=0.75)

weights <- seq(0,2,by=0.025)
weights

#out_dir_grid <- out_dir

wma_grid_search <- function(table0,table1,table2, weights, out_dir_grid, mark){
  
  #### prep ####
  names_for_columns <-  c("weight", "n_eqtl_ss", "n_eqtl_se", "zscore_cor_ss","zscore_cor_se", "intersected_eqtls", 'corr_ss_se', "lm_ss_se_intersect", "lm_ss_se_coef", "zscore_cor_ss_overlap","zscore_cor_se_overlap")
  
  output_table <- as.data.frame(t(names_for_columns))
  colnames(output_table) <- names_for_columns
  
  # output_table_se <- data.frame("weight", "n_eqtl_se","zscore_cor_se")
  # colnames(output_table_se) <- c("weight", "n_eqtl_se","zscore_cor_se")
  # 
  merged_dt <- merge(table1, table2, by="snp_gene")
  merged_dt <- merge(merged_dt, table0, by="snp_gene" )
  
  #### weights ####
  for(weight1 in weights){
    
    for(i in 1:nrow(merged_dt)) {
      eqtl <- merged_dt[i,]
      
      weighted_z_score <- get_weighted_z_score(
        zScores=c(eqtl$OverallZScore.x, eqtl$OverallZScore.y),
        weights=c(weight1,1),
        sampleSizes = c(1, 1))
      
      merged_dt[i,"ZW"] <- weighted_z_score   
      
      weighted_z_score_SE <- get_weighted_z_score_SE(
        zScores=c(eqtl$OverallZScore.x, eqtl$OverallZScore.y),
        weights=c(weight1,1),
        SEs = c(eqtl$beta_se.x, eqtl$beta_se.y)
      )
      merged_dt$ZW_SE[i] <- weighted_z_score_SE 
      
    }
    merged_dt$pvalue_ZW <- convert.z.score(merged_dt$ZW)
    merged_dt$pvalue_ZW_SE <- convert.z.score(merged_dt$ZW_SE)
    
    
    tab1 <- merged_dt[order(merged_dt$pvalue_ZW),]
    tab1 <- tab1 %>% distinct(feature_id.x, .keep_all= TRUE)
    multTestOut_pvalue_ZW <- multtest::mt.rawp2adjp(tab1$pvalue_ZW, proc = c("BH"))$adjp
    tab1$cfdr_ZW <- multTestOut_pvalue_ZW[,2]
    wma_eqtls <- length(tab1[tab1$cfdr_ZW <0.05, ]$snp_gene)
    
    tab2 <- merged_dt[order(merged_dt$pvalue_ZW_SE),]
    tab2 <- tab2 %>% distinct(feature_id.x, .keep_all= TRUE)
    multTestOut_pvalue_ZW_SE <- multtest::mt.rawp2adjp(tab2$pvalue_ZW_SE, proc = c("BH"))$adjp
    tab2$cfdr_ZW_SE <- multTestOut_pvalue_ZW_SE[,2]
    wma_eqtls_se <- length(tab2[tab2$cfdr_ZW_SE <0.05, ]$snp_gene)
    
    #tab_ss <- merge(tab1, tab2, by="snp_gene", suffixes = c('tab1', "tab2"))
    tab <- merge(tab1, tab2, by="snp_gene", suffixes = c('.tab1', ".tab2"))
    
    
    wma_output_tab <- data.frame(paste(weight1),wma_eqtls,wma_eqtls_se, cor.test(tab1$ZW, tab1$OverallZScore, method =  c("pearson"))[[4]], cor.test(tab2$ZW, tab2$OverallZScore, method =  c("pearson"))[[4]], dim(tab)[1], cor.test(tab$ZW.tab1, tab$ZW_SE.tab1, method =  c("pearson"))[[4]] , 
                                 cor.test(tab$ZW.tab1, tab$OverallZScore.tab1, method =  c("pearson"))[[4]],
                                 cor.test(tab$ZW_SE.tab1, tab$OverallZScore.tab1, method =  c("pearson"))[[4]],
                                 lm(tab$ZW.tab1 ~ tab$ZW_SE.tab1)[[1]][[1]],lm(tab$ZW.tab1 ~ tab$ZW_SE.tab1)[[1]][[2]]
    )
    
    colnames(wma_output_tab) <- c("weight", "n_eqtl_ss", "n_eqtl_se", "zscore_cor_ss","zscore_cor_se", "intersected_eqtls", 'corr_ss_se',  "zscore_cor_ss_overlap","zscore_cor_se_overlap","lm_ss_se_intersect", "lm_ss_se_coef" )
    
    output_table<- rbind(output_table,wma_output_tab)
    
    merged_dt_out <- merge( tab1[,c("snp_gene", "cfdr_ZW")],tab2[,c("snp_gene", "cfdr_ZW_SE")], by="snp_gene", all = T)
    merged_dt_out <- merge(merged_dt_out, merged_dt, by="snp_gene", all = T)
    merged_dt_out[mapply(is.na, merged_dt_out)] <- 0
    
    name_for_output_weight <- paste(out_dir_grid, weight1, mark, ".tsv", sep='')
    write.table(merged_dt_out, name_for_output_weight, col.names = T, quote=F, sep ='\t')
    
  }
  
  
  #}
  #output_table <- merge(output_table_ss,output_table_se, by="weight")
  
  # output_table <- output_table[2:nrow(output_table),]
  # output_table$n_eqtl <- as.numeric(output_table$n_eqtl)
  # output_table$proportion_from_ds <- output_table$n_eqtl/dim(merged_dt)[1]
  # output_table$max <- output_table$n_eqtl/max(output_table$n_eqtl)
  # 
  # output_table$proportion_from_sign_x <- output_table$n_eqtl/(dim(merged_dt[merged_dt$BH.x <0.05,])[1])
  # output_table$proportion_from_sign_y <- output_table$n_eqtl/(dim(merged_dt[merged_dt$BH.y <0.05,])[1])
  output_table[mapply(is.infinite, output_table)] <- 0
  
  name_for_output2 <- paste(out_dir_grid, mark, "_summary.tsv", sep='')
  write.table(output_table, name_for_output2, col.names = T, quote=F, sep ='\t')
  
  #return(output_table)
  
}


out_dir_grid = "/groups/umcg-franke-scrna/tmp01/projects/sc-eqtlgen-consortium-pipeline/ongoing/wg3-QTL-mapping/LIMIX/run_limix/ipsc/output_29_11/grid/"


# correlation with bulk 
WMA_processing_limix(table0=limix_qtl_bulk, table1=limix_qtl_ss2_67,table2=limix_qtl_10x_25, sample_size1=67, sample_size2=26,mark="bulk_cor_ss2_67_10x_26_ACE", out_dir)
WMA_processing_limix(table0=limix_qtl_ss2_87, table1=limix_qtl_ss2_67,table2=limix_qtl_10x_25, sample_size1=67, sample_size2=26,mark="ss2_87_cor_ss2_67_10x_26_BCE", out_dir)

WMA_processing_limix(table0=limix_qtl_bulk, table1=limix_qtl_ss2_67,table2=limix_qtl_ss2_25, sample_size1=67, sample_size2=26,mark="bulk_cor_ss2_67_ss2_26_ACD", out_dir)
WMA_processing_limix(table0=limix_qtl_ss2_87, table1=limix_qtl_ss2_67,table2=limix_qtl_ss2_25, sample_size1=67, sample_size2=26,mark="ss2_87_cor_ss2_67_ss2_26_BCD", out_dir)

WMA_processing_limix(table0=limix_qtl_bulk, table1=limix_qtl_ss2_25_2,table2=limix_qtl_ss2_25, sample_size1=26, sample_size2=26,mark="bulk_cor_ss2_26_ss2_26_AEF", out_dir)
WMA_processing_limix(table0=limix_qtl_ss2_87, table1=limix_qtl_ss2_25_2,table2=limix_qtl_ss2_25, sample_size1=26, sample_size2=26,mark="ss2_87_cor_ss2_26_ss2_26_BEF", out_dir)


