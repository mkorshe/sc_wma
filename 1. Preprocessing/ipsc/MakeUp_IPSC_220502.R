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
