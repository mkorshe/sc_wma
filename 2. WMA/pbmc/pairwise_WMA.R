# #!/usr/bin/env bash

DS1=$1
DS2=$2
NGROUP=$3
EXPR_STAT=$4
OUTDIR=$5
NEQTLS=$6

cat  << EOF

# paralel with readlines check

library(dplyr)
library(multtest)
library(stringr)
library(readr)
#ds_ids <- 'ng_vs_stemi_v2'
# # 
# #
 # ds_id1='ng'
 # ds_id2='stemi_v2'
 # N=3
 # expression_stats='/groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/LIMIX/WMA_round_II/pseudobulk/input/Input_data_5ds.tsv'
 # 
 # output_path_for_summaries='/groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/LIMIX/WMA_round_II/pseudobulk/output/WMA/'
 # nrow_datatable=680377
# #

# expression_stats='/groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/LIMIX/WMA_round_II/monocyte/input/Input_data_5ds.tsv'
# expression_stats_bulk=read.table('/groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/LIMIX/WMA_round_II/pseudobulk/input/Input_data_5ds.tsv',sep='\t',header=T)
# colnames(expression_stats_bulk) <-  singificant_coeqtls_readlines_conames[1,]

# write.table(expression_stats_bulk, '/groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/LIMIX/WMA_round_II/pseudobulk/input/Input_data_5ds.tsv',sep='\t', col.names = T, quote = F, row.names = F)

ds_id1 <- '$DS1' #v2_1m
ds_id2 <- '$DS2' #v3_1m
N <- as.numeric('$NGROUP')
nrow_datatable <- as.numeric('$NEQTLS')#680377
expression_stats <- '$EXPR_STAT'
output_path_for_summaries <- '$OUTDIR'



dataset_idents <- paste0(ds_id1, '_vs_', ds_id2)
n_of_groups <-  100

#chunk = $CHUNK

#dir.create('/Users/korshe/createtestdir', showWarnings = TRUE, recursive = FALSE, mode = "0777")
#dataset_combinations <- c('v2_1m_vs1m_v3', 'v2_1m_vsstemi_v2','v2_1m_vsstemi_v3' ,'ng_vsstemi_v2','ng_vsstemi_v3','ng_vs1m_v2','ng_vs1m_v3', 'stemi_v3_vs1m_v3', 'stemi_v2_vs1m_v3', 'stemi_v2_vsstemi_v3')

weights_ss <- read.table('/groups/umcg-franke-scrna/tmp01/projects/sc-eqtlgen-consortium-pipeline/ongoing/wg3-QTL-mapping/LIMIX/run_limix/wma/pseudobulk_1000_eTQLGen/PBMC_dataset_specifications.csv', sep=',', header=T)

samplesize_ids <- colnames(weights_ss)[2:7]

weights_ss[,'ids_wma'] <- c('v2_1m','v3_1m', 'ng', 'stemi_v2', 'stemi_v3')

weighting_ids <- read.table('/groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/LIMIX/WMA_round_II/weights_list.tsv', sep='\t', header=T)

weighting_ids <- weighting_ids[,'x']
weighting_ids <- weighting_ids[-33]
# mean_donor_L_variance_cell
#### Functions #####

zToCor <- function(z, df){
  t <- qt(pnorm(-abs(z),log.p = T), df, log.p=T)
  r <- (t/sqrt(df+t^2))^2
  # r <- -log10(r)
  return(r)
}

get_weighted_z_score <- function(zScores, sampleSizes, weights) {
  # set starting values
  weightSum <- 0
  zSum <- 0
  
  for(i in 1:length(zScores)) {
    # add to the Z
    w <- sqrt(sampleSizes[i]) * weights[i]
    zSum <- zSum + (w * zScores[i])
    # add to the total sum
    weightSum <- weightSum + (w*w)
  }
  weightedZ <- zSum / sqrt(weightSum)
  return(weightedZ)
}

get_weighted_z_score_no_ss <- function(zScores, weights) {
  # set starting values
  weightSum <- 0
  zSum <- 0
  
  for(i in 1:length(zScores)) {
    # add to the Z
    w <- 1 * weights[i]
    zSum <- zSum + (w * zScores[i])
    # add to the total sum
    weightSum <- weightSum + (w*w)
  }
  weightedZ <- zSum / sqrt(weightSum)
  return(weightedZ)
}

get_weighted_z_score_SE <- function(zScores,weights, SEs) {
  # set starting values
  weightSum <- 0
  zSum <- 0
  
  for(i in 1:length(zScores)) {
    # add to the Z
    w <- 1/((SEs[i])^2) *weights[i]
    zSum <- zSum + (w * zScores[i])
    # add to the total sum
    weightSum <- weightSum + (w*w)
  }
  weightedZ <- zSum / sqrt(weightSum)
  return(weightedZ)
}

convert.z.score<-function(z, one.sided=NULL) {
  if(is.null(one.sided)) {
    pval = pnorm(-abs(z));
    pval = 2 * pval
  } else if(one.sided=="-") {
    pval = pnorm(z);
  } else {
    pval = pnorm(-z);
  }
  return(pval);
}

#### Calculation ZWeighted 


#dataset_5ds <- read.table(paste(expression_stats), sep='\t', header= T)

#head(dataset_5ds)

#dir.create(paste(output_path_for_summaries,  '5ds_wma', sep=''), showWarnings = TRUE, recursive = FALSE, mode = "0777")

#nrow_datatable <- nrow(merged_dt)
#write.table(colnames(dataset_5ds), '/groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/LIMIX/cfdr_limix/limix_qtl_and_expr_stat/all_genes_eqtlgen/monocytes_limix/expression_files/expression_stat/qtls_and_expression/colnames_5ds.tsv', sep='\t', row.names = F, quote = F, col.names=F)
colnames_for_output <- read.table('/groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/LIMIX/WMA_round_II/colnames_for_input.tsv', sep='\t', header = T)
colnames_for_output <- colnames_for_output[,'x']
#rm(merged_dt)

#### stage 2 ####

f = 1

nrow_datatable_s <-  round(nrow_datatable/n_of_groups,0)
samplesize_id <- samplesize_ids[1]
#N <- 2

#n_of_groups = 4
#for(N in 1:n_of_groups){

dir.create(paste(output_path_for_summaries, dataset_idents ,'/group_', N, sep=''), showWarnings = TRUE, recursive = FALSE, mode = "0777")

n_min <-  ((N-1)*nrow_datatable_s) +1
n_max <-  N*nrow_datatable_s 
sequence_of_values <- seq(n_min, n_max,1000) 

for(f in sequence_of_values){
  singificant_coeqtls_readlines_conames <- read_lines(paste(expression_stats), skip=0, n_max=1)
  singificant_coeqtls_readlines_conames <-  as.data.frame(str_split_fixed(singificant_coeqtls_readlines_conames, "\t", 278)) 
  
  singificant_coeqtls_readlines <- read_lines(paste(expression_stats), skip=f, n_max=1000)
  eqtl <-  as.data.frame(str_split_fixed(singificant_coeqtls_readlines, "\t", 278)) 
  colnames(eqtl) <- singificant_coeqtls_readlines_conames[1,]
  
  
  eqtl_ds1 <- colnames(eqtl)[grep(paste0(ds_id1), colnames(eqtl))]
  eqtl_ds2 <- colnames(eqtl)[grep(paste0(ds_id2), colnames(eqtl))]
  eqtl <- eqtl %>%  select(one_of(c('snp_gene',eqtl_ds1,eqtl_ds2)))
  colnames(eqtl) <- gsub(ds_id1, "ds1", colnames(eqtl) )
  colnames(eqtl) <- gsub(ds_id2, "ds2", colnames(eqtl) )
  
  
  # wma analysis 
  # SE-based
  for(i in 1:nrow(eqtl)){
    
    weighted_z_score_se <- get_weighted_z_score_SE(zScores= as.numeric(c(eqtl[i,'corrected_zscoreds1'], eqtl[i,'corrected_zscoreds2'])), weights=c(1,1),SEs=as.numeric(c(eqtl[i,'beta_seds1'], eqtl[i,'beta_seds2'])))
    
    eqtl[i,"ZW_weight_SE"] <- weighted_z_score_se
  }
  eqtl[,'pvalue_ZW_SE'] <- convert.z.score(eqtl[,'ZW_weight_SE'])
  
  for (id in weighting_ids){
  #  id = substring(id,1, nchar(id)-1)
    
    weight1 <- paste(id,"ds1",sep='')
    weight2 <- paste(id,"ds2",sep='')
 
    for(i in 1:nrow(eqtl)){
      
      weighted_z_score_se <- get_weighted_z_score_SE(zScores= as.numeric(c(eqtl[i,'corrected_zscoreds1'], eqtl[i,'corrected_zscoreds2'])),
                                                     weights=as.numeric(c(eqtl[i,weight1][[1]],eqtl[i,weight2][[1]])),
                                                     SEs=as.numeric(c(eqtl[i,'beta_seds1'], eqtl[i,'beta_seds2'])))
      
      eqtl[i,paste0("ZW_SE", id)] <- weighted_z_score_se
    }
    eqtl[,paste0("pvalue_SE", id)] <- convert.z.score(eqtl[,paste0("ZW_SE", id)])
    
  }
  
  
  for(samplesize_id in samplesize_ids){
    
    weights_SS_counts <-  weights_ss[weights_ss[,'ids_wma'] %in% c(ds_id1, ds_id2) ,  samplesize_id]
    
    SS1 <-   weights_SS_counts[1]
    SS2 <-   weights_SS_counts[2]
    
    sampleSizes = c(SS1,SS2)
   
     for(i in 1:nrow(eqtl)){
      
      weighted_z_score_ss <- get_weighted_z_score(zScores= as.numeric(c(eqtl[i,'corrected_zscoreds1'], eqtl[i,'corrected_zscoreds2'])),sampleSizes = sampleSizes,  weights=c(1,1))
      
      eqtl[i,paste0("ZW_",samplesize_id)] <- weighted_z_score_ss
    }
    eqtl[,paste0("pvalue_",samplesize_id)] <- convert.z.score(eqtl[,paste0("ZW_",samplesize_id)])
    
    # doing WMA with weights 
    
    for (id in weighting_ids){
      print(id)
      #id = substring(id,1, nchar(id)-1)
      
      weight1 <- paste(id,"ds1",sep='')
      weight2 <- paste(id,"ds2",sep='')
      
      for(i in 1:nrow(eqtl)){
        
        weighted_z_score_ss <- get_weighted_z_score(zScores= as.numeric(c(eqtl[i,'corrected_zscoreds1'], eqtl[i,'corrected_zscoreds2'])), sampleSizes = sampleSizes, weights=c(as.numeric(eqtl[i,weight1][[1]]),as.numeric(eqtl[i,weight2][[1]])))
        
        eqtl[i,paste0("ZW_",samplesize_id, id)] <- weighted_z_score_ss
      }
      eqtl[,paste0("pvalue_",samplesize_id, id)] <- convert.z.score(eqtl[,paste0("ZW_",samplesize_id, id)])
      
      
      print(paste(samplesize_id, id, 'weight is in a done'))
      
    }
  }
  
  param_cols1 <- colnames(eqtl)[grep('ds1', colnames(eqtl))]
  param_cols2 <- colnames(eqtl)[grep('ds2', colnames(eqtl))]

  eqtl_polished <- eqtl %>%  select(-one_of(c(param_cols1,param_cols2)), snp_gene)
  
  # writing output
  # write.table(colnames(eqtl_polished), paste(output_path_for_summaries, "colnames_pairwise.tsv", sep=''), sep = '\t', col.names = T, row.names = F, quote=F)
  
  write.table(eqtl_polished, paste(output_path_for_summaries,dataset_idents ,'/group_', N ,'/',f, ".tsv", sep=''), sep = '\t', col.names = F, row.names = F)
  
  f <- f+1000
  #print(F)
  #i = i+nrow_datatable_s
  
}


EOF
