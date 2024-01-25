## this script takes limma output tables from 
## 1_1_CTLH_proteomes_2023_data_cleaning_filtering_230428
## 2_CTLH_proteomes_2023_limma_stats_final_230501

## calculates median centered values from imputed data
## combines with non-imputed MaxLFQ
## saves a single long-format data table for use later (in plotting)

library(tidyverse)
library(magrittr)


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

# filenames of non-imputed MaxLFQ table (long format)
long_files = list.files(dat.dir,pattern = "msbooster_long.tsv")

# df to loop over
to_loop =
  tibble(limma_files,
         imputed_files,
         long_files)


#################
#### outputs ####
#################

# long out
long_comb_out_file = paste0(dat.dir,"Combined_protein_measurements_long.tsv")


#######################
#### MEDIAN CENTER ####
#######################


i=1
long_all = tibble()
for(i in 1:nrow(to_loop)){
  
  
  #### which samples are we operating on
  # limma file
  lim_file = to_loop[i,1] %>% as.character()
  # imputed data file
  imp_file = to_loop[i,2] %>% as.character()
  # non imputed data file
  long_file = to_loop[i,3] %>% as.character()
  #this protein
  #this_prot = to_loop[i,3] %>% as.character()
  # name for this prot pair
  this_prot_pair = str_split(lim_file,pattern="_")[[1]][1]
  # get a vector of each proteins name for filtering columns on
  these_prots = str_split(this_prot_pair,pattern="-")[[1]]
  
  
  ###################
  #### LOAD DATA ####
  ###################
  
  ## load imputed MaxLFQ table
  imp_dat = read_tsv(paste0(dat.dir,imp_file))
  
  ## load non-imputed MaxLFQ table
  ## get into nicer format for use later
  long = 
    read_tsv(paste0(dat.dir,long_file)) %>%
    dplyr::select(ID,Gene,sample,LFQ)
  
  colnames(long)[3] = "condition"
  
  ## get unique gene IDs mapping table
  gene_ids =
    long %>%
    dplyr::select(ID,Gene) %>%
    distinct()
  
  
  ##########################
  #### MEDIAN CENTERING ####
  ##########################
  
  ## do median centering for each protein ON THE IMPUTED VALUES
  ## center each protein expression by subtracting median (median-centered Log2 intensities)
  ## have to transpose the wide matrix first as scale/sweep operates on columns by default
  new_mat = 
    imp_dat %>%
    dplyr::select(ID,starts_with(these_prots)) %>%
    column_to_rownames("ID") %>%
    t()
  
  # get protein medians
  med.att <- apply(new_mat, 2, median)
  
  # perform centering, join to protein IDs
  cent_dat = 
    sweep(data.matrix(new_mat), 2, med.att) %>%
    t() %>%
    data.frame() %>%
    rownames_to_column("ID") %>%
    tibble()
  
  
  ## get into long format
  keycol <- "condition"
  valuecol <- "centd"
  gathercols <- colnames(dplyr::select(cent_dat,-ID))
  
  cent_dat %<>%
    gather_(keycol, valuecol, gathercols) %>%
    separate(condition, into=c("prot","type","rep"),sep="_",remove=F) %>%
    mutate(group = paste(prot,type,sep="_")) %>%
    distinct()
  
  # check all IDs to check how close to zero the protein medians are
  chk = 
    cent_dat %>%
    group_by(ID) %>%
    summarise(median_id = median(centd)) %>%
    arrange(median_id)
  
  ## tiny errors, proceed
  min(chk$median_id)
  max(chk$median_id)
  
  
  ##########################
  #### COMBINING TABLES ####
  ##########################
  
  ## combine the median centered data with the imputed data
  ## first get into long format
  imp_dat %<>%
    dplyr::select(ID,name,imputed,num_NAs,starts_with(these_prots))
  
  ## get into long format
  keycol <- "condition"
  valuecol <- "imputed_LFQ"
  gathercols <- colnames(dplyr::select(imp_dat,starts_with(these_prots)))

  imp_dat %<>%
    gather_(keycol, valuecol, gathercols)
  
  comb = 
    left_join(cent_dat,imp_dat,by=c("ID","condition")) %>%
    dplyr::select(ID,condition,prot,centd,imputed_LFQ,imputed,num_NAs) %>%
    left_join(gene_ids)
  
  comb %<>%
    left_join(long,by=c("ID","condition")) %>%
    dplyr::select(ID,Gene.x,prot,condition,imputed,num_NAs,LFQ,imputed_LFQ,centd) %>%
    mutate(LFQ=log2(LFQ)) %>%
    rename(Gene=Gene.x)
  
  ## add to main output table
  long_all %<>% rbind.data.frame(comb)
  
}


write.table(long_all,
            long_comb_out_file,
            col.names = T,
            row.names = F,
            quote=F,
            sep="\t")