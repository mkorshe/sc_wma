setwd('/Users/korshe/Documents/Data_Groningen/WMA_sep2022/T_preprocessed/')


ef <- read.table('T_by_B_BH_bulk_EF.tsv',sep='\t', header=T)
df <- read.table('T_by_B_BH_bulk_DF.tsv',sep='\t', header=T)
cd <- read.table('T_by_B_BH_bulk_CD.tsv',sep='\t', header=T)
ce <- read.table('T_by_B_BH_bulk_CE.tsv',sep='\t', header=T)
#df <- ce

f1_score <- function(df){
  df$increase_eGenes_SS <- df$eGenes/df[df$weight == 'Number_of_donors','eGenes']
  #  df$increase_eGenes_SE<- df$eGenes/df[df$weight == 'ZW_weight_SE','eGenes']
  
  df$increase_eQTL_SS <- df$eQTls/df[df$weight == 'Number_of_donors','eQTls']
  #  df$increase_eQTL_SE<- df$eQTls/df[df$weight == 'ZW_weight_SE','eQTls']
  
  df$increase_Precision_SS <- df$Precision/df[df$weight == 'Number_of_donors','Precision']
  df$increase_Sensitivity_SS <- df$Sensitivity/df[df$weight == 'Number_of_donors','eGenes']
  
  df$F1_score <-  2*(df$Precision*df$Sensitivity)/(df$Sensitivity+df$Precision)
  df$increase_eGenes_SS <- df$eGenes/df[df$weight == 'Number_of_donors','eGenes']
#  df$increase_eGenes_SE<- df$eGenes/df[df$weight == 'ZW_weight_SE','eGenes']
  
  df$increase_eQTL_SS <- df$eQTls/df[df$weight == 'Number_of_donors','eQTls']
#  df$increase_eQTL_SE<- df$eQTls/df[df$weight == 'ZW_weight_SE','eQTls']
  
  df$increase_Precision_SS <- df$Precision/df[df$weight == 'Number_of_donors','Precision']
  df$increase_Sensitivity_SS <- df$Sensitivity/df[df$weight == 'Number_of_donors','eGenes']
  
  df <- df[order(df$F1_score, decreasing = T),]
  return(df)
}
ef<- f1_score(ef)
df<- f1_score(df)
cd<- f1_score(cd)
ce<- f1_score(ce)


ss2 <- merge(cd, df, by='weight')


# calculated the weighted mean 

E <- 25
C <- 67

for(i in 1:nrow(ss2)){
ss2[i,'weighted_mean_F1'] <- 
  weighted.mean(c(ss2[i,'F1_score.x'],
                  ss2[i,'F1_score.y']), 
                
                c((C+E),
                  (E+E)
                ))
}
# ss2$final_f1 <- (ss2$F1_score.x+ss2$F1_score.y)/2
# ss2_10x$final_f1 <- (ss2_10x$F1_score.x+ss2_10x$F1_score.y)/2


ss2 <- ss2[order(ss2$weighted_mean_F1, decreasing = T),]
head(ss2)



# ss2$final_f1 <- (ss2$F1_score.x+ss2$F1_score.y)/2
# ss2_10x$final_f1 <- (ss2_10x$F1_score.x+ss2_10x$F1_score.y)/2

ss2[1,c('weight','eGenes.x','eQTls.x','Precision.x','Sensitivity.x','eGenes.y','eQTls.y','Precision.y','Sensitivity.y','F1_score.x','F1_score.y','weighted_mean_F1')]
ss2[ss2$weight == 'Number_of_donors',c('weight','eGenes.x','eQTls.x','Precision.x','Sensitivity.x','eGenes.y','eQTls.y','Precision.y','Sensitivity.y','weighted_mean_F1')]

best_weight_genes <-ss2[ss2$weighted_mean_F1_genes==max(ss2$weighted_mean_F1_genes),]$weight[1]
best_weight_eqtls <-ss2[ss2$weighted_mean_F1_eqtls==max(ss2$weighted_mean_F1_eqtls),]$weight[1]
best_weight_genes_eqtl <-ss2[ss2$weighted_mean_F1_genes_eqtls==max(ss2$weighted_mean_F1_genes_eqtls),]$weight[1]


test_df <- as.data.frame(t(rbind(ss2[ss2$weight==best_weight_genes,],
                                 ss2[ss2$weight==best_weight_eqtls,],
                                 ss2[ss2$weight==best_weight_genes_eqtl,],
                                 ss2[ss2$weight == 'Number_of_donors',] )))
colnames(test_df) <- c('genes','eqtls','gene_eqtl','SS')
##


ss2_10x <-  merge(ce, ef, by='weight')

for(i in 1:nrow(ss2_10x)){
  ss2_10x[i,'weighted_mean_F1'] <- 
    weighted.mean(c(ss2_10x[i,'F1_score.x'],
                    ss2_10x[i,'F1_score.y']), 
                  
                  c((C+E),
                    (E+E)
                  ))
}

ss2_10x <- ss2_10x[order(ss2_10x$weighted_mean_F1, decreasing = T),]

ss2_10x[1,c('weight','eGenes.x','eQTls.x','Precision.x','Sensitivity.x','eGenes.y','eQTls.y','Precision.y','Sensitivity.y','F1_score.x','F1_score.y','weighted_mean_F1')]

ss2_10x[ss2_10x$weight == 'Number_of_donors',c('weight','eGenes.x','eQTls.x','Precision.x','Sensitivity.x','eGenes.y','eQTls.y','Precision.y','Sensitivity.y','weighted_mean_F1')]

##### Trying to reimplement metanalaysis 
setwd('/Users/korshe/Documents/Data_Groningen/WMA_sep2022/T_preprocessed/')


ds <- read.table('T_by_B_BH_bulk_CE.tsv',sep='\t', header=T)



ef <- read.table('T_by_B_BH_bulk_EF.tsv',sep='\t', header=T)
df <- read.table('T_by_B_BH_bulk_DF.tsv',sep='\t', header=T)
cd <- read.table('T_by_B_BH_bulk_CD.tsv',sep='\t', header=T)
ce <- read.table('T_by_B_BH_bulk_CE.tsv',sep='\t', header=T)
#df <- ce

get_new_set_pf_F1 <- function(ds){
  ds$f1_genes <- ds$TP_SS_gene/(ds$TP_SS_gene+0.5*(ds$FP_SS_gene+ds$FN_SS_gene))
  ds$increase_f1_genes <- ds$f1_genes/ds[ds$weight == 'Number_of_donors',]$f1_genes
  
  ds$f1_eQTLs  <- ds$TP_SS/(ds$TP_SS+0.5*(ds$FP_SS+ds$FN_SS))
  ds$increase_f1_eQTLs <- ds$f1_eQTLs/ds[ds$weight == 'Number_of_donors',]$f1_eQTLs
  
  ds$mean_f1_genes_eqtls <- (ds$f1_eQTLs+ds$f1_genes)/2
  ds$increase_f1_genes_eqtls <- ds$mean_f1_genes_eqtls/ds[ds$weight == 'Number_of_donors',]$mean_f1_genes_eqtls
  
  return(ds)
}

ef<- get_new_set_pf_F1(ef)
df<- get_new_set_pf_F1(df)
cd<- get_new_set_pf_F1(cd)
ce<- get_new_set_pf_F1(ce)

# meta-analysis on f1 from egenes only
E <- 25
C <- 67

ss2 <- merge(cd, df, by='weight')
# calculated the weighted mean 

for(i in 1:nrow(ss2)){
  ss2[i,'weighted_mean_F1_genes'] <- 
    weighted.mean(c(ss2[i,'f1_genes.x'],
                    ss2[i,'f1_genes.y']), 
                  
                  c((C+E),
                    (E+E)
                  ))
  ss2[i,'weighted_mean_F1_eqtls'] <- 
    weighted.mean(c(ss2[i,'f1_eQTLs.x'],
                    ss2[i,'f1_eQTLs.y']), 
                  
                  c((C+E),
                    (E+E)
                  ))
  ss2[i,'weighted_mean_F1_genes_eqtls'] <- 
    weighted.mean(c(ss2[i,'increase_f1_genes_eqtls.x'],
                    ss2[i,'increase_f1_genes_eqtls.y']), 
                  
                  c((C+E),
                    (E+E)
                  ))
}
# ss2$final_f1 <- (ss2$F1_score.x+ss2$F1_score.y)/2
# ss2_10x$final_f1 <- (ss2_10x$F1_score.x+ss2_10x$F1_score.y)/2

zeros_cols <-    ss2$weight[grep('zeros', ss2$weight)]
ss2<- ss2[ss2$weight %!in%zeros_cols, ]

ss2 <- ss2[order(ss2$weighted_mean_F1_genes, decreasing = T),]
head(ss2)


ss2[1,c('weight','eGenes.x','eQTls.x','Precision.x','Sensitivity.x','eGenes.y','eQTls.y','Precision.y','Sensitivity.y','f1_genes.x','f1_genes.y','weighted_mean_F1_genes')]
ss2[ss2$weight == 'Number_of_donors',c('weight','eGenes.x','eQTls.x','Precision.x','Sensitivity.x','eGenes.y','eQTls.y','f1_genes.x','f1_genes.y','weighted_mean_F1_genes')]


ss2_10x <-  merge(ce, ef, by='weight')

for(i in 1:nrow(ss2_10x)){
  ss2_10x[i,'weighted_mean_F1_genes'] <- 
    weighted.mean(c(ss2_10x[i,'f1_genes.x'],
                    ss2_10x[i,'f1_genes.y']), 
                  
                  c((C+E),
                    (E+E)
                  ))
  ss2_10x[i,'weighted_mean_F1_eqtls'] <- 
    weighted.mean(c(ss2_10x[i,'f1_eQTLs.x'],
                    ss2_10x[i,'f1_eQTLs.y']), 
                  
                  c((C+E),
                    (E+E)
                  ))
  ss2_10x[i,'weighted_mean_F1_genes_eqtls'] <- 
    weighted.mean(c(ss2_10x[i,'increase_f1_genes_eqtls.x'],
                    ss2_10x[i,'increase_f1_genes_eqtls.y']), 
                  
                  c((C+E),
                    (E+E)
                  ))
}


ss2_10x[1,c('weight','eGenes.x','eQTls.x','Precision.x','Sensitivity.x','eGenes.y','eQTls.y','Precision.y','Sensitivity.y','f1_genes.x','f1_genes.y','weighted_mean_F1_genes')]
ss2_10x[ss2_10x$weight == 'Number_of_donors',c('weight','eGenes.x','eQTls.x','Precision.x','Sensitivity.x','eGenes.y','eQTls.y','f1_genes.x','f1_genes.y','weighted_mean_F1_genes')]



# ss2$final_f1 <- (ss2$F1_score.x+ss2$F1_score.y)/2
# ss2_10x$final_f1 <- (ss2_10x$F1_score.x+ss2_10x$F1_score.y)/2




ss2_10x <- ss2_10x[order(ss2_10x$weighted_mean_F1_genes, decreasing = T),]

ss2_10x[1,c('weight','eGenes.x','eQTls.x','Precision.x','Sensitivity.x','eGenes.y','eQTls.y','Precision.y','Sensitivity.y','F1_score.x','F1_score.y','increase_f1_genes.x','weighted_mean_F1_genes')]

ss2_10x[ss2_10x$weight == 'Number_of_donors',c('weight','eGenes.x','eQTls.x','Precision.x','Sensitivity.x','eGenes.y','eQTls.y','Precision.y','Sensitivity.y','increase_f1_genes.x','weighted_mean_F1_genes')]
####

ss2_10x <- ss2_10x[order(ss2_10x$weighted_mean_F1_genes_eqtls, decreasing = T),]

ss2_10x[1,c('weight','eGenes.x','eQTls.x','Precision.x','Sensitivity.x','eGenes.y','eQTls.y','Precision.y','Sensitivity.y','F1_score.x','F1_score.y','increase_f1_genes.x','weighted_mean_F1_genes_eqtls')]

ss2_10x[ss2_10x$weight == 'Number_of_donors',c('weight','eGenes.x','eQTls.x','Precision.x','Sensitivity.x','eGenes.y','eQTls.y','Precision.y','Sensitivity.y','increase_f1_genes.x','weighted_mean_F1_genes_eqtls')]



#######



# dataset_combos <- c('CE', 'CD')
'%!in%' <- function(x,y)!('%in%'(x,y))

dataframe_ipsc$weight <- rownames(dataframe_ipsc)
zeros_columns <-    dataframe_ipsc$weight[grep('zeros', dataframe_ipsc$weight)]

# dataframe_ipsc <- dataframe_ipsc[dataframe_ipsc$weight %!in% c('A', )]

 prepare_agr_rank_ipsc <- function(dataframe_ipsc){
  dataframe_ipsc$weight <- rownames(dataframe_ipsc)
  dataframe_ipsc <- dataframe_ipsc[dataframe_ipsc$weight %!in% c('A',zeros_columns ),]
  
  dataframe_ipsc= dataframe_ipsc[order(dataframe_ipsc$Precision,dataframe_ipsc$Sensitivity,decreasing = c(T,T)),]
  dataframe_ipsc["rank_genes"] = rev(c(1:nrow(dataframe_ipsc)))
  
  # dataframe_ipsc$pair <- paste(dataframe_ipsc$weight, dataframe_ipsc$rank_genes, sep='_')
  return(dataframe_ipsc)
 }
 
 ce_rank <- prepare_agr_rank_ipsc(ce)
 cd_rank <- prepare_agr_rank_ipsc(cd)
 ef_rank <- prepare_agr_rank_ipsc(ef)
 df_rank <- prepare_agr_rank_ipsc(df)
 
 # SS2 and SS2
 
 SS2 <- merge(cd_rank, df_rank, by='weight')
 
 for(n in 1:nrow(SS2)){
   SS2$weighted_mean_rank[n] <-weighted.mean(c(SS2$rank_genes.x[n], SS2$rank_genes.y[n]), c(87,50))
   
 }
 SS2$mean_rank <- (SS2$rank_genes.x + SS2$rank_genes.y)/2
 SS2 <- SS2[order(SS2$weighted_mean_rank, decreasing = T),]
 
 # Summarising top weight:
 # paste('the most optimal weight (according to the weighted rank) is',SS2$weighted_mean_rank[1],', which leads to increase to number of eGenes in CD dataset to', round(SS2$eGenes.x[1]/SS2[SS2$weight == 'Number_of_donors',]$eGenes.x, digits = 2),'comparing to SS-based weighing and ' , round(SS2$eGenes.x[1]/SS2[SS2$weight == 'SE',]$eGenes.x, 'comparing to SE-based weighing.')
 
 
 summary_weighting <- as.data.frame(c(paste(SS2$weight[1], 'eGenes'), 
                                      paste(SS2$weight[1], 'eQTL'),
                                      paste(SS2$weight[1], 'Sensitivity'),
                                      paste(SS2$weight[1], 'Precision'),
                                      
                                      paste('SS', 'eGenes'), 
                                      paste('SS', 'eQTL'),
                                      paste('SS', 'Sensitivity'),
                                      paste('SS', 'Precision'),
                                      
                                      paste('SE', 'eGenes'), 
                                      paste('SE', 'eQTL'),
                                      paste('SE', 'Sensitivity'),
                                      paste('SE', 'Precision')
                                      ))
 
 colnames(summary_weighting) <- 'parameter'
 
 summary_weighting$CD_absolute_numbers <- c(SS2$eGenes.x[1],SS2$eQTls.x[1],SS2$Sensitivity.x[1],SS2$Precision.x[1],
                                            SS2[SS2$weight == 'Number_of_donors',]$eGenes.x,SS2[SS2$weight == 'Number_of_donors',]$eQTls.x,SS2[SS2$weight == 'Number_of_donors',]$Sensitivity.x,SS2[SS2$weight == 'Number_of_donors',]$Precision.x,
                                            SS2[SS2$weight == 'SE',]$eGenes.x,SS2[SS2$weight == 'SE',]$eQTls.x,SS2[SS2$weight == 'SE',]$Sensitivity.x,SS2[SS2$weight == 'SE',]$Precision.x)
 SS2_and_SS2 <- SS2
 summary_weighting$CD_increase_comparing_to_SS <- c(SS2$eGenes.x[1]/SS2[SS2$weight == 'Number_of_donors',]$eGenes.x,
                                                 SS2$eQTls.x[1]/SS2[SS2$weight == 'Number_of_donors',]$eQTls.x,
                                                 SS2$Sensitivity.x[1]/SS2[SS2$weight == 'Number_of_donors',]$Sensitivity.x,
                                                 SS2$Precision.x[1]/SS2[SS2$weight == 'Number_of_donors',]$Precision.x,
                                                 1,1,1,1,
                                                 
                                            SS2[SS2$weight == 'SE',]$eGenes.x/SS2[SS2$weight == 'Number_of_donors',]$eGenes.x,
                                            SS2[SS2$weight == 'SE',]$eQTls.x/SS2[SS2$weight == 'Number_of_donors',]$eQTls.x,
                                            SS2[SS2$weight == 'SE',]$Sensitivity.x/SS2[SS2$weight == 'Number_of_donors',]$Sensitivity.x,
                                            SS2[SS2$weight == 'SE',]$Precision.x/SS2[SS2$weight == 'Number_of_donors',]$Precision.x)
 
 
 ## for dataset 2
 
 summary_weighting$DF_absolute_numbers <- c(SS2$eGenes.y[1],SS2$eQTls.y[1],SS2$Sensitivity.y[1],SS2$Precision.y[1],
                                            SS2[SS2$weight == 'Number_of_donors',]$eGenes.y,SS2[SS2$weight == 'Number_of_donors',]$eQTls.y,SS2[SS2$weight == 'Number_of_donors',]$Sensitivity.y,SS2[SS2$weight == 'Number_of_donors',]$Precision.y,
                                            SS2[SS2$weight == 'SE',]$eGenes.y,SS2[SS2$weight == 'SE',]$eQTls.y,SS2[SS2$weight == 'SE',]$Sensitivity.y,SS2[SS2$weight == 'SE',]$Precision.y)
 
 summary_weighting$DF_increase_comparing_to_SS <- c(SS2$eGenes.y[1]/SS2[SS2$weight == 'Number_of_donors',]$eGenes.y,
                                                    SS2$eQTls.y[1]/SS2[SS2$weight == 'Number_of_donors',]$eQTls.y,
                                                    SS2$Sensitivity.y[1]/SS2[SS2$weight == 'Number_of_donors',]$Sensitivity.y,
                                                    SS2$Precision.y[1]/SS2[SS2$weight == 'Number_of_donors',]$Precision.y,
                                                    1,1,1,1,
                                                    
                                                    SS2[SS2$weight == 'SE',]$eGenes.y/SS2[SS2$weight == 'Number_of_donors',]$eGenes.y,
                                                    SS2[SS2$weight == 'SE',]$eQTls.y/SS2[SS2$weight == 'Number_of_donors',]$eQTls.y,
                                                    SS2[SS2$weight == 'SE',]$Sensitivity.y/SS2[SS2$weight == 'Number_of_donors',]$Sensitivity.y,
                                                    SS2[SS2$weight == 'SE',]$Precision.y/SS2[SS2$weight == 'Number_of_donors',]$Precision.y)
 
 summary_weighting$CD_increase_comparing_to_SS <- round(summary_weighting$CD_increase_comparing_to_SS,digits = 3)
 summary_weighting$DF_increase_comparing_to_SS <- round(summary_weighting$DF_increase_comparing_to_SS,digits = 3)
 
 summary_weighting_SS2 <- summary_weighting
 
 
 # SS2 + 10x 
 
 
 SS2 <- merge(ef_rank, ce_rank, by='weight')
 
 for(n in 1:nrow(SS2)){
   SS2$weighted_mean_rank[n] <-weighted.mean(c(SS2$rank_genes.x[n], SS2$rank_genes.y[n]), c(50,87))
   
 }
 SS2$mean_rank <- (SS2$rank_genes.x + SS2$rank_genes.y)/2
 SS2 <- SS2[order(SS2$weighted_mean_rank, decreasing = T),]
 
 # Summarising top weight:
 # paste('the most optimal weight (according to the weighted rank) is',SS2$weighted_mean_rank[1],', which leads to increase to number of eGenes in CD dataset to', round(SS2$eGenes.x[1]/SS2[SS2$weight == 'Number_of_donors',]$eGenes.x, digits = 2),'comparing to SS-based weighing and ' , round(SS2$eGenes.x[1]/SS2[SS2$weight == 'SE',]$eGenes.x, 'comparing to SE-based weighing.')
 
 
 summary_weighting <- as.data.frame(c(paste(SS2$weight[1], 'eGenes'), 
                                      paste(SS2$weight[1], 'eQTL'),
                                      paste(SS2$weight[1], 'Sensitivity'),
                                      paste(SS2$weight[1], 'Precision'),
                                      
                                      paste('SS', 'eGenes'), 
                                      paste('SS', 'eQTL'),
                                      paste('SS', 'Sensitivity'),
                                      paste('SS', 'Precision'),
                                      
                                      paste('SE', 'eGenes'), 
                                      paste('SE', 'eQTL'),
                                      paste('SE', 'Sensitivity'),
                                      paste('SE', 'Precision')
 ))
 
 colnames(summary_weighting) <- 'parameter'
 
 summary_weighting$EF_absolute_numbers <- c(SS2$eGenes.x[1],SS2$eQTls.x[1],SS2$Sensitivity.x[1],SS2$Precision.x[1],
                                            SS2[SS2$weight == 'Number_of_donors',]$eGenes.x,SS2[SS2$weight == 'Number_of_donors',]$eQTls.x,SS2[SS2$weight == 'Number_of_donors',]$Sensitivity.x,SS2[SS2$weight == 'Number_of_donors',]$Precision.x,
                                            SS2[SS2$weight == 'SE',]$eGenes.x,SS2[SS2$weight == 'SE',]$eQTls.x,SS2[SS2$weight == 'SE',]$Sensitivity.x,SS2[SS2$weight == 'SE',]$Precision.x)
 
 summary_weighting$EF_increase_comparing_to_SS <- c(SS2$eGenes.x[1]/SS2[SS2$weight == 'Number_of_donors',]$eGenes.x,
                                                    SS2$eQTls.x[1]/SS2[SS2$weight == 'Number_of_donors',]$eQTls.x,
                                                    SS2$Sensitivity.x[1]/SS2[SS2$weight == 'Number_of_donors',]$Sensitivity.x,
                                                    SS2$Precision.x[1]/SS2[SS2$weight == 'Number_of_donors',]$Precision.x,
                                                    1,1,1,1,
                                                    
                                                    SS2[SS2$weight == 'SE',]$eGenes.x/SS2[SS2$weight == 'Number_of_donors',]$eGenes.x,
                                                    SS2[SS2$weight == 'SE',]$eQTls.x/SS2[SS2$weight == 'Number_of_donors',]$eQTls.x,
                                                    SS2[SS2$weight == 'SE',]$Sensitivity.x/SS2[SS2$weight == 'Number_of_donors',]$Sensitivity.x,
                                                    SS2[SS2$weight == 'SE',]$Precision.x/SS2[SS2$weight == 'Number_of_donors',]$Precision.x)
 
 
 ## for dataset 2
 
 summary_weighting$CE_absolute_numbers <- c(SS2$eGenes.y[1],SS2$eQTls.y[1],SS2$Sensitivity.y[1],SS2$Precision.y[1],
                                            SS2[SS2$weight == 'Number_of_donors',]$eGenes.y,SS2[SS2$weight == 'Number_of_donors',]$eQTls.y,SS2[SS2$weight == 'Number_of_donors',]$Sensitivity.y,SS2[SS2$weight == 'Number_of_donors',]$Precision.y,
                                            SS2[SS2$weight == 'SE',]$eGenes.y,SS2[SS2$weight == 'SE',]$eQTls.y,SS2[SS2$weight == 'SE',]$Sensitivity.y,SS2[SS2$weight == 'SE',]$Precision.y)
 
 summary_weighting$CE_increase_comparing_to_SS <- c(SS2$eGenes.y[1]/SS2[SS2$weight == 'Number_of_donors',]$eGenes.y,
                                                    SS2$eQTls.y[1]/SS2[SS2$weight == 'Number_of_donors',]$eQTls.y,
                                                    SS2$Sensitivity.y[1]/SS2[SS2$weight == 'Number_of_donors',]$Sensitivity.y,
                                                    SS2$Precision.y[1]/SS2[SS2$weight == 'Number_of_donors',]$Precision.y,
                                                    1,1,1,1,
                                                    
                                                    SS2[SS2$weight == 'SE',]$eGenes.y/SS2[SS2$weight == 'Number_of_donors',]$eGenes.y,
                                                    SS2[SS2$weight == 'SE',]$eQTls.y/SS2[SS2$weight == 'Number_of_donors',]$eQTls.y,
                                                    SS2[SS2$weight == 'SE',]$Sensitivity.y/SS2[SS2$weight == 'Number_of_donors',]$Sensitivity.y,
                                                    SS2[SS2$weight == 'SE',]$Precision.y/SS2[SS2$weight == 'Number_of_donors',]$Precision.y)
 
 summary_weighting$EF_increase_comparing_to_SS <- round(summary_weighting$EF_increase_comparing_to_SS,digits = 3)
 summary_weighting$CE_increase_comparing_to_SS <- round(summary_weighting$CE_increase_comparing_to_SS,digits = 3)
 
 summary_weighting_SS2_10x <- summary_weighting
 SS2_and_10x <- SS2
 
 
 write.table(summary_weighting_SS2_10x, '10x_SS2_weighted_rank.tsv',sep='\t', quote = F, row.names = F)
 write.table(summary_weighting_SS2, 'SS2_SS2_weighted_rank.tsv',sep='\t', quote = F, row.names = F)
 
 
 write.table(SS2_and_10x, 'SS2_10x_weighted_rank_all.tsv',sep='\t', quote = F, row.names = F)
 write.table(SS2_and_SS2, 'SS2_SS2_weighted_rank_all.tsv',sep='\t', quote = F, row.names = F)
 
 