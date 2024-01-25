## this script takes limma output tables from 
## 2_CTLH_proteomes_2023_limma_stats_final_230501

## plots heatmap with dendrogram using just sinificant proteins for that comparison

library(tidyverse)
library(magrittr)
library(pheatmap)

#################
#### OPTIONS ####
#################

options(max.print=100)


################
#### inputs ####
################

# specify root folder that contains the data
base = "C:/Users/imnot/OneDrive/Postdoc/Projects/Maitland_proteomes/Analysis_v3_separated_Apr23/"

# data with imputed MaxLFQ tables
dat.dir = paste0(base,"data/")

# specify the folder for limma output tables
stat_dir = paste0(base,"stats/")


# filenames of limma output tables
limma_files = list.files(stat_dir,pattern = "limma_significant.tsv")

# filenames of imputed MaxLFQ table
imputed_files = list.files(dat.dir,pattern = "_wide_imputed.tsv")

# df to loop over
to_loop =
  tibble(limma_files,
       imputed_files)

## add proteins to this df
to_loop %<>%
  separate(limma_files,remove=F,into=c("prot1","prot2")) %>%
  gather(key,protein,prot1:prot2) %>%
  dplyr::select(-key)


#################
#### outputs ####
#################

plot_dir = paste0(base,"plots/")



#############
#### PCA ####
#############


i=1
for(i in 1:nrow(to_loop)){
  
  # limma file
  lim_file = to_loop[i,1] %>% as.character()
  # imputed data file
  imp_file = to_loop[i,2] %>% as.character()
  #this protein
  this_prot = to_loop[i,3] %>% as.character()
  # name for this prot pair
  this_prot_pair = str_split(lim_file,pattern="_")[[1]][1]
  
  # load limma table
  ## filter this table on just the sifnificant proteins FOR THIS ONE PROTEIN
  sig_prots =
    read_tsv(paste0(stat_dir,lim_file)) %>%
    separate(contrast,into="prot") %>%
    filter(prot==!!this_prot)

  
  ## load imputed MaxLFQ table
  wide = read_tsv(paste0(dat.dir,imp_file))
  
  ## filter the wide table to just significant proteins
  wide %<>%
    filter(ID %in% sig_prots$ID) %>%
    dplyr::select(Gene,ID,contains(!!this_prot))
  
  ## get into numeric matrix
  mat =
    wide %>%
    dplyr::select(contains(!!this_prot)) %>%
    as.matrix()

  #create data frame for column annotations on heatmap
  anno =
    data.frame(sample=as.character(colnames(mat))) %>%
    mutate(sample2 = sample) %>%
    column_to_rownames("sample2") %>%
    separate(sample, into=c("Protein","condition","Rep")) %>%
    dplyr::select(-Rep,-Protein)
  
  ## change colors of annotation
  anno$condition=as.factor(anno$condition)
  
  ann_colors = list(
    condition = c(control="#377eb8ff", knockout="#e4211cff"))
  

  pdf(paste0(plot_dir,"Heatmap_",this_prot,".pdf"),2,1.5)
  pheatmap(mat,scale="row",
           color=colorRampPalette(c("navy", "white", "red"))(50),
           cutree_cols=2,
           cutree_rows=2,
           annotation_col = anno,
           show_rownames=F,
           show_colnames = F,
           clustering_method="complete",
           annotation_colors = ann_colors,
           border_color = NA,
           treeheight_row = 0.3,
           treeheight_col = 0.3
           )
  Sys.sleep(3)
  dev.off()
  
}

