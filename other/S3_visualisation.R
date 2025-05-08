

library("data.table")   
library("dplyr")  
library("stringr")  
library("readr")  
library("qvalue")  
library(gridExtra)
library(tidyverse)
library(scales)


# output_path_for_summaries <- '/Users/mkorshe/Documents/WMA_input_files/monocyte/output/aggregated/'

output_path_for_summaries_TP <- '/Users/mkorshe/Documents/WMA_input_files/replication_revision/by_ref_TP/'

summary_weight_ref_all <- fread(paste0(output_path_for_summaries_TP,'all_weights_MTC_processing.tsv' ), data.table = F)

summary_weight_ref_all_eqtlgen <- summary_weight_ref_all[summary_weight_ref_all$ref =='eqtlgen',]
summary_weight_ref_all_onek1k <- summary_weight_ref_all[summary_weight_ref_all$ref =='onek1k',]

####
summary_weight_ref_all_M_cells_donor <- summary_weight_ref_all[summary_weight_ref_all$weight %in% summary_weight_ref_all$weight[grep('M_cells_donor',summary_weight_ref_all$weight )],]


summary_weight_ref_all_SS_SE_only <- summary_weight_ref_all[summary_weight_ref_all$weight %in% c('SS','SE','M_cells_donor','M_UMI_cell','M_UMI_donor','T_cells','T_UMI'),]

summary_weight_ref_all_SS_SE_only_eqtlgen <- summary_weight_ref_all_SS_SE_only[summary_weight_ref_all_SS_SE_only$ref == 'eqtlgen',]


summary_weight_ref_all_SS_SE_only_onek1k<- summary_weight_ref_all_SS_SE_only[summary_weight_ref_all_SS_SE_only$ref == 'onek1k',]


out_tab <- as.data.frame(c('SE','M_cells_donor','M_UMI_cell','M_UMI_donor','T_cells','T_UMI'))
colnames(out_tab) <- 'SS_weight'


out_tab_average_increase_egenes <- out_tab
# 
# for(ss_weight in out_tab$SS_weight){
#   
#   for(ref in c('eqtlgen', 'onek1k')){
#     
#     summary_weight_ref_all_SS_SE_only_ref <- summary_weight_ref_all_SS_SE_only[summary_weight_ref_all_SS_SE_only$ref == ref,]
#     
#     summary_weight_ref_all_SS_SE_only_eqtlgen_ds1 <- summary_weight_ref_all_SS_SE_only_ref[summary_weight_ref_all_SS_SE_only_ref$weight %in% c('SS'),]
#     
#     
#     summary_weight_ref_all_SS_SE_only_eqtlgen_ds2 <- summary_weight_ref_all_SS_SE_only_ref[summary_weight_ref_all_SS_SE_only_ref$weight %in% c(ss_weight),]
#     
#     out_tab[out_tab$SS_weight == ss_weight, ref] <-  sum(summary_weight_ref_all_SS_SE_only_eqtlgen_ds1$eGenes < summary_weight_ref_all_SS_SE_only_eqtlgen_ds2$eGenes)    
#     
#     out_tab_average_increase_egenes[out_tab_average_increase_egenes$SS_weight == ss_weight, ref] <-  mean(summary_weight_ref_all_SS_SE_only_eqtlgen_ds2$eGenes - summary_weight_ref_all_SS_SE_only_eqtlgen_ds1$eGenes)    
#     
#     
#   }
#   
# }
# 
# ###
for(ss_weight in out_tab$SS_weight){

  for(ref in c('eqtlgen', 'onek1k')){

    summary_weight_ref_all_SS_SE_only_ref <- summary_weight_ref_all_SS_SE_only[summary_weight_ref_all_SS_SE_only$ref == ref,]

    summary_weight_ref_all_SS_SE_only_eqtlgen_ds1 <- summary_weight_ref_all_SS_SE_only_ref[summary_weight_ref_all_SS_SE_only_ref$weight %in% c('SS'),]


    summary_weight_ref_all_SS_SE_only_eqtlgen_ds2 <- summary_weight_ref_all_SS_SE_only_ref[summary_weight_ref_all_SS_SE_only_ref$weight %in% c(ss_weight),]

    out_tab[out_tab$SS_weight == ss_weight, ref] <-  sum(summary_weight_ref_all_SS_SE_only_eqtlgen_ds1$eGenes <= summary_weight_ref_all_SS_SE_only_eqtlgen_ds2$eGenes)

    out_tab_average_increase_egenes[out_tab_average_increase_egenes$SS_weight == ss_weight, ref] <-  mean(summary_weight_ref_all_SS_SE_only_eqtlgen_ds2$eGenes - summary_weight_ref_all_SS_SE_only_eqtlgen_ds1$eGenes)


  }

}
out_tab
#### Plot 1#### 

# out_tab[out_tab$onek1k ==9,]$onek1k <- 7
df_long <- out_tab %>%
  pivot_longer(cols = c(eqtlgen, onek1k), names_to = "dataset", values_to = "value")

# Reorder SS_weight **within each dataset**

# df_long$new_id

df_long <- df_long %>%
  mutate(SS_weight = recode(SS_weight,
                            "SE" = "Standard error",
                            "M_cells_donor" = "Average number of cells",
                            "M_UMI_cell" = "Counts per cell",
                            "M_UMI_donor" = "Counts per donor",
                            "T_cells" = "Total number of cells",
                            "T_UMI" = "Total counts per cohort"
  ))

ss_colors <- c(
  "Standard error" = "#62C2A6",
  "Average number of cells" = "#D44757",
  "Counts per cell" = "#F3764F",
  "Counts per donor" = "#9A164B",
  "Total number of cells" = "#3589BC",
  "Total counts per cohort" = "#5C59A5"
)
# 
# df_long <- df_long %>%
#   group_by(dataset) %>%
#   mutate(SS_weight = fct_reorder(SS_weight, value))
df_long$value <- as.numeric(df_long$value)
# df_long <- df_long[df_long$dataset =='onek1k',]
df_long <- df_long %>%
  group_by(dataset) %>%
  mutate(SS_weight = fct_reorder(SS_weight, value,.desc = F)) %>%
  ungroup()

# Plot
pl_a <- ggplot(df_long[df_long$dataset == 'eqtlgen',], aes(y = SS_weight, x = value, fill = SS_weight)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = value), hjust = -0.2, size = 3.5) +
  # facet_wrap(~ dataset, ncol = 2, scales = "free_y") +
  scale_fill_manual(values = ss_colors) +
  theme_light() +
  labs(y = NULL, x = "Number of dataset combinations", title = "a") +
  theme(
    axis.text.y = element_text(size = 10),
    legend.position = "none"
  ) +
  coord_cartesian(xlim = c(0, max(df_long$value) + 1)) +
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank())

pl_a

pl_b <- ggplot(df_long[df_long$dataset == 'onek1k',], aes(y = SS_weight, x = value, fill = SS_weight)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = value), hjust = -0.2, size = 3.5) +
  # facet_wrap(~ dataset, ncol = 2, scales = "free_y") +
  scale_fill_manual(values = ss_colors) +
  theme_light() +
  labs(y = NULL, x = "Number of dataset combinations", title = "b") +
  theme(
    axis.text.y = element_text(size = 10),
    legend.position = "none"
  ) +
  coord_cartesian(xlim = c(0, max(df_long$value) + 1)) +
  theme(axis.line = element_line(color='black'),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank())

pl_b

####

##### panel 2 #####


############

summary_weight_ref_all_M_cells_donor <- summary_weight_ref_all[summary_weight_ref_all$weight %in% summary_weight_ref_all$weight[grep('M_cells_donor',summary_weight_ref_all$weight )],]


summary_weight_ref_all_M_cells_donor_out_for_ds <- as.data.frame(unique(summary_weight_ref_all_M_cells_donor$weight))
colnames(summary_weight_ref_all_M_cells_donor_out_for_ds) <- 'param'


summary_weight_ref_all_M_cells_donor_out_for_ds_onek1k <- summary_weight_ref_all_M_cells_donor_out_for_ds

summary_weight_ref_all_M_cells_donor_out_for_ds_onek1k_TP <- summary_weight_ref_all_M_cells_donor_out_for_ds


summary_weight_ref_all_M_cells_donor_out_for_ds_eqtlgen <- summary_weight_ref_all_M_cells_donor_out_for_ds

summary_weight_ref_all_M_cells_donor_out_for_ds_onek1k <- summary_weight_ref_all_M_cells_donor_out_for_ds
summary_weight_ref_all_M_cells_donor_out_for_ds_onek1k_TP <- summary_weight_ref_all_M_cells_donor_out_for_ds

summary_weight_ref_all_M_cells_donor_out_for_ds_eqtlgen <- summary_weight_ref_all_M_cells_donor_out_for_ds
summary_weight_ref_all_M_cells_donor_out_for_ds_eqtlgen_TP <- summary_weight_ref_all_M_cells_donor_out_for_ds


for(extra_weight in out_tab_m_cells_average_increase_egenes$extra_weight){
  
  # for(ref in c('eqtlgen', 'onek1k')){
  for(ds in unique(summary_weight_ref_all_SS_SE_only_ref$ds)){
    
    summary_weight_ref_ds_onek1k <- summary_weight_ref_all_M_cells_donor[summary_weight_ref_all_M_cells_donor$ds == ds & summary_weight_ref_all_M_cells_donor$ref =='onek1k',]
    
    summary_weight_ref_ds_onek1k$dif_TP <- summary_weight_ref_ds_onek1k$TP.eGene - summary_weight_ref_ds_onek1k[summary_weight_ref_ds_onek1k$weight == 'M_cells_donor',]$TP.eGene
    
    summary_weight_ref_ds_onek1k$dif_egene <- summary_weight_ref_ds_onek1k$eGenes - summary_weight_ref_ds_onek1k[summary_weight_ref_ds_onek1k$weight == 'M_cells_donor',]$eGenes
    
    summary_weight_ref_all_M_cells_donor_out_for_ds_onek1k[,ds] <- summary_weight_ref_ds_onek1k$dif_egene
    
    summary_weight_ref_all_M_cells_donor_out_for_ds_onek1k_TP[,ds] <- summary_weight_ref_ds_onek1k$dif_TP
    
    ##
    
    
    summary_weight_ref_ds_eqtlgen <- summary_weight_ref_all_M_cells_donor[summary_weight_ref_all_M_cells_donor$ds == ds & summary_weight_ref_all_M_cells_donor$ref =='eqtlgen',]
    
    summary_weight_ref_ds_eqtlgen$dif_TP <- summary_weight_ref_ds_eqtlgen$TP.eGene - summary_weight_ref_ds_eqtlgen[summary_weight_ref_ds_eqtlgen$weight == 'M_cells_donor',]$TP.eGene
    
    summary_weight_ref_ds_eqtlgen$dif_egene <- summary_weight_ref_ds_eqtlgen$eGenes - summary_weight_ref_ds_eqtlgen[summary_weight_ref_ds_eqtlgen$weight == 'M_cells_donor',]$eGenes
    
    summary_weight_ref_all_M_cells_donor_out_for_ds_eqtlgen[,ds] <- summary_weight_ref_ds_eqtlgen$dif_egene
    
    summary_weight_ref_all_M_cells_donor_out_for_ds_eqtlgen_TP[,ds] <- summary_weight_ref_ds_eqtlgen$dif_TP
    
    
    
  }
  
}


make_a_heatmap <- function(df, title_paste){
  
  # Long format
  df_long <- df %>%
    pivot_longer(-param, names_to = "comparison", values_to = "value")
  
  # Assign color category (green/red/white)
  df_long <- df_long %>%
    mutate(fill_color = case_when(
      value < 0 ~ "negative",
      value > 0 ~ "positive",
      TRUE ~ "zero"
    ))
  
  # Scale abs(value) to [0, 1] per column for opacity
  df_long <- as.data.frame(df_long)
  df_long <- df_long %>%
    group_by(comparison) %>%
    mutate(
      max_val = max(abs(value), na.rm = TRUE),
      opacity = ifelse(max_val == 0, 0, abs(value) / max_val)
    ) %>%
    select(-max_val) %>%
    ungroup()
  
  comparison_names <- c(
    "ng_vs_stemi_v2"    = "Wijst and Blokland (V2)",
    "ng_vs_stemi_v3"    = "Wijst and Blokland (V3)",
    "ng_vs_v2_1m"       = "Wijst and Oelen (V2)",
    "ng_vs_v3_1m"       = "Wijst and Oelen (V3)",
    "stemi_v3_vs_v3_1m" = "Blokland (V3) and Oelen (V3)",
    "v2_1m_vs_stemi_v2" = "Blokland (V2) and Oelen (V2)",
    "v2_1m_vs_stemi_v3" = "Blokland (V3) and Oelen (V2)",
    "v2_1m_vs_v3_1m"    = "Oelen (V2) and Oelen (V3)",
    "v3_1m_vs_stemi_v2" = "Oelen (V3) and Blokland (V2)",
    "5ds"               = "All datasets"
  )
  
  # Apply the renaming
  df_long <- df_long %>%
    mutate(comparison = recode(comparison, !!!comparison_names))
  
  # Map colors
  fill_colors <- c("negative" = "#D93C33", "positive" = "#1F9E5A", "zero" = "white")
  
  df_long$param <- gsub('M_cells_donor','',df_long$param)
  df_long$param[df_long$param=='maf_'] <- 'maf'
  # Plot
  ggplot(df_long[df_long$param != '',], aes(y = comparison, x = param)) +
    geom_tile(aes(fill = fill_color,alpha = opacity)) +
    geom_text(aes(label = value), size = 3.5) +
    scale_fill_manual(values = fill_colors, guide = "none") +
    scale_alpha(range = c(0.2, 1), guide = "none") +  # controls opacity range
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.title = element_blank(),
      panel.grid = element_blank() 
    )+ ggtitle(title_paste)
  
}


make_a_heatmap(summary_weight_ref_all_M_cells_donor_out_for_ds_eqtlgen_TP, 'TP_by_eqtlgen')
make_a_heatmap(summary_weight_ref_all_M_cells_donor_out_for_ds_eqtlgen, 'All eGenes by eqtlgen')

summary_weight_ref_all_M_cells_donor_out_for_ds_eqtlgen$param <- gsub('M_cells_donor', '', summary_weight_ref_all_M_cells_donor_out_for_ds_eqtlgen$param )

make_a_heatmap(summary_weight_ref_all_M_cells_donor_out_for_ds_onek1k_TP, 'TP eGenes by OneK1K')
make_a_heatmap(summary_weight_ref_all_M_cells_donor_out_for_ds_onek1k, 'c \n \n All eGenes by OneK1K')

pl_c <- make_a_heatmap(summary_weight_ref_all_M_cells_donor_out_for_ds_onek1k, 'c \n \n Difference of eGenes; by OneK1K')

##### Panel 3 #####



out_tab_extra_weight <- as.data.frame(unique(summary_weight_ref_all_M_cells_donor$weight))
colnames(out_tab_extra_weight) <- 'extra_weight'

out_tab_extra_weight$eqtlgen <- 0
out_tab_extra_weight$onek1k <- 0

# 
# out_tab_average_N_ds_egenes <- out_tab_extra_weight
out_tab_m_cells <- out_tab_extra_weight
out_tab_m_cells_average_increase_egenes <- out_tab_extra_weight


for(extra_weight in out_tab_extra_weight$extra_weight){
  
  for(ref in c('eqtlgen', 'onek1k')){
    
    summary_weight_ref_all_SS_SE_only_ref <- summary_weight_ref_all_M_cells_donor[summary_weight_ref_all_M_cells_donor$ref == ref,]
    
    summary_weight_ref_all_SS_SE_only_eqtlgen_ds1 <- summary_weight_ref_all_SS_SE_only_ref[summary_weight_ref_all_SS_SE_only_ref$weight %in% c('M_cells_donor'),]
    
    
    summary_weight_ref_all_SS_SE_only_eqtlgen_ds2 <- summary_weight_ref_all_SS_SE_only_ref[summary_weight_ref_all_SS_SE_only_ref$weight %in% c(extra_weight),]
    
    out_tab_m_cells[out_tab_m_cells$extra_weight == extra_weight, ref] <-  sum(summary_weight_ref_all_SS_SE_only_eqtlgen_ds1$eGenes < summary_weight_ref_all_SS_SE_only_eqtlgen_ds2$eGenes)    
    
    out_tab_m_cells_average_increase_egenes[out_tab_m_cells_average_increase_egenes$extra_weight == extra_weight, ref] <-  mean(summary_weight_ref_all_SS_SE_only_eqtlgen_ds2$eGenes - summary_weight_ref_all_SS_SE_only_eqtlgen_ds1$eGenes)    
    
    
  }
  
}
out_tab_m_cells




###



out_tab_m_cells <- as.data.frame(c(unique(summary_weight_ref_all_M_cells_donor$weight)))
colnames(out_tab_m_cells) <- 'extra_weight'



out_tab_m_cells_average_increase_egenes <- out_tab_m_cells
out_tab_m_cells_average_increase_egenes$eqtlgen <- 0
out_tab_m_cells_average_increase_egenes$onek1k <- 0

summary_weight_ref_all_M_cells_donor <- summary_weight_ref_all[summary_weight_ref_all$weight %in% summary_weight_ref_all$weight[grep('M_cells_donor',summary_weight_ref_all$weight )],]

summary_weight_ref_all_M_cells_donor <- summary_weight_ref_all_M_cells_donor[order(summary_weight_ref_all_M_cells_donor$weight),]


for(extra_weight in out_tab_m_cells_average_increase_egenes$extra_weight){
  
  for(ref in c('eqtlgen', 'onek1k')){
    
    vec1 <- summary_weight_ref_all_M_cells_donor[summary_weight_ref_all_M_cells_donor$weight == extra_weight & summary_weight_ref_all_M_cells_donor$ref == ref,]$eGenes
    vec2 <- summary_weight_ref_all[ summary_weight_ref_all$ref == ref& summary_weight_ref_all$weight =='M_cells_donor',]$eGenes
    
    
    
    out_tab_m_cells_average_increase_egenes[out_tab_m_cells_average_increase_egenes$extra_weight == extra_weight,ref] <- length(vec1[vec1 > vec2])
    
  }
  
}

out_tab_m_cells <- out_tab_m_cells[out_tab_m_cells$eqtlgen >0 & out_tab_m_cells$onek1k > 0,]

out_tab_m_cells$extra_weight <- gsub('M_cells_donor','',out_tab_m_cells$extra_weight)

df_sorted <- out_tab_m_cells %>%
  mutate(total = eqtlgen + onek1k) %>%
  arrange(desc(-total)) %>%
  mutate(extra_weight = factor(extra_weight, levels = extra_weight))  # Preserve order

# Convert to long format
df_long <- df_sorted %>%
  pivot_longer(cols = c(eqtlgen, onek1k), names_to = "reference", values_to = "count")

# Plot


pl_d <- ggplot(df_long, aes(x = extra_weight, y = count)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  coord_flip() +  # Flip for readability
  facet_grid(cols = vars(reference)) + 
  labs(x = "Extra weight", y = "Count (More eGenes than 'Average cells per donor')") +
  theme_minimal() +
  ggtitle(" d \n \n Number of combinations where extra weight is outperforming 'Average cells per donor'") +
  theme(axis.text.y = element_text(size = 10)) +  # controls opacity range
  theme_minimal() +
  theme(
    axis.text.x = element_text( hjust = 1),
    axis.title = element_blank(),
    panel.grid = element_blank() 
  )

pl_d


grid.arrange(pl_a, pl_b, ncol =2)
grid.arrange( pl_c,pl_d, ncol =1)

grid.arrange(glob_1, glob_2)

grid.arrange(
  grobs = gl,
  widths = c(2, 1, 1),
  layout_matrix = rbind(c(1, 2, NA),
                        c(3, 3, 4))
)
