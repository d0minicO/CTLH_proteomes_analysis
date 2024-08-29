## this script takes limma output tables from 
## 2_CTLH_proteomes_2023_limma_stats_final_230501

## performs overlaps of significant proteins between the proteins

library(tidyverse)
library(magrittr)
library(UpSetR)

#################
#### OPTIONS ####
#################

options(max.print=100)


################
#### inputs ####
################

# specify root folder that contains the data
base = "C:/Users/dowens/OneDrive/Postdoc/Projects/Maitland_proteomes/Analysis_v3_separated_Apr23/"

# specify the data folder for protein tables
dat.dir = paste0(base,"data/")


# specify the folder for limma output tables
stat_dir = paste0(base,"stats/")


# filenames of limma output tables
in_files = list.files(stat_dir,pattern = "limma_significant.tsv")



## CTLH complex members
CTLH_symbols =
  c("MKLN1",
    "RANBP9",
    "RMND5A",
    "RMND5B",
    "MAEA",
    "GID8",
    "WDR26",
    "ARMC8",
    "GID4")

CTLH_ids =
  c(
    "Q9UL63",
    "Q96S59",
    "Q9H871",
    "Q96G75",
    "Q7L5Y9",
    "Q9NWU2",
    "Q9H7D7",
    "Q8IUR7",
    "Q8IVV7"
  )

CTLH = 
  tibble(
    CTLH_symbols,
    CTLH_ids
  )





#################
#### outputs ####
#################

plot_dir = paste0(base,"plots/")



##################
#### OVERLAPS ####
##################


#######
## 1 ##
#######

## basic plots for figures

### ALL 4 DATA SETS TOGETHER


in_file = in_files[1]
sig_all = tibble()
for(in_file in in_files){
  
  temp=
    read_tsv(paste0(stat_dir,in_file)) %>%
    separate(contrast,into="prot")
  
  sig_all %<>% rbind.data.frame(temp)
  
}

#saveRDS(sig_all,paste0(base,"sig_all_data.Rds"))


#sig_all %>% filter(Gene=="MKLN1")


####### UPSET ########

## get into wide format for upsetr
## must be df not tbl_df, must be 0 not NA
upset_data = 
  sig_all %>%
  mutate(sig = 1) %>%
  dplyr::select(ID,prot,sig) %>%
  spread(key=prot,value=sig) %>%
  data.frame() %>%
  replace(is.na(.), 0)

# make the plot
upset_plot <- upset(upset_data, order.by = "freq", decreasing = T,text.scale = .8)

# save it
pdf(paste0(plot_dir,"Upset_all_DEPs.pdf"),3,2.5)
print(upset_plot)
dev.off()



####### 4-WAY VENN ########


## construct list to perform overlaps with ggvenndiagram
a = filter(sig_all,prot=="RMND5A") %>% pull(ID)
b = filter(sig_all,prot=="MAEA") %>% pull(ID)
c = filter(sig_all,prot=="MKLN1") %>% pull(ID)
d = filter(sig_all,prot=="WDR26") %>% pull(ID)




## overlap with all

sig_list = 
  list(
    a,b,c,d
  )


names(sig_list) = 
  c(
    "RMND5A",
    "MAEA",
    "MKLN1",
    "WDR26"
  )


## plot the venn for all
ggVennDiagram(sig_list)

ggsave(filename = paste0(plot_dir,"venn_sig_prots_all_overlaps.pdf"),
       width=4.5,height = 3.5)



### DUO DATA SETS SEPARATELY

in_file = in_files[1]
for(in_file in in_files){
  
  # name for this prot pair
  this_prot = str_split(in_file,pattern="_")[[1]][1]
  
  # load limma table
  ## filter this table on just the sifnificant proteins
  sig_prots =
    read_tsv(paste0(stat_dir,in_file)) %>%
    separate(contrast,into="prot")
    
  ### using UPSET
  
  ## get into wide format for upsetr
  ## must be df not tbl_df, must be 0 not NA
  test = 
    sig_prots %>%
    mutate(sig = 1) %>%
    dplyr::select(ID,prot,sig) %>%
    spread(key=prot,value=sig) %>%
    data.frame() %>%
    replace(is.na(.), 0)
  
  # make the plot
  upset_plot <- upset(test, order.by = "freq", decreasing = T)
  
  # save it
  pdf(paste0(plot_dir,"upset_all_DEPs_",this_prot,".pdf"),3.5,2.5)
  print(upset_plot)
  dev.off()
  
  
  ### now venn diagrams
  
  prot1 = str_split(this_prot,pattern="-")[[1]][1]
  prot2 = str_split(this_prot,pattern="-")[[1]][2]
  
  ## construct list to perform overlaps with ggvenndiagram
  a = filter(sig_prots,prot==prot1) %>% pull(ID)
  b = filter(sig_prots,prot==prot2) %>% pull(ID)
  
  ## overlap with all
  
  sig_list = 
    list(
      a,b
    )
  
  
  names(sig_list) = 
    c(
      prot1,
      prot2
    )
  
  
  ## plot the venn for all
  ggVennDiagram(sig_list)
  
  ggsave(filename = paste0(plot_dir,"venn_sig_prots_all_overlaps_",this_prot,".pdf"),
         width=4.5,height = 3.5)
  
}





#######
## 2 ##
#######

#### INTERSECT PREVIOUS CTLH DATA

## get count for each protein for how many datasets it was significant in
sig_counts =
  sig_all %>%
  group_by(Gene,ID) %>%
  dplyr::count() %>%
  arrange(desc(n))

# Correlate list of changed proteins in different KO with previous CTLH data (Gid4 BioID, AP-MS in Onea et al., RanBPM proteome and ubiquitome in Maitland et al., 2021)

## save output table of CTLH interactors that are significantly changed in CTLH proteomes 
top_genes =
  sig_counts %>%
  left_join(sig_all) %>%
  filter(prot %in% c("RMND5A","MAEA")) %>%
  filter(n>1)

# Gid4 BioID
# AP-MS in Onea et al.
# RanBPM proteome and ubiquitome in Maitland et al., 2021

## get Uniprot IDs in the biotype of choice to do overlaps with

## GID4 BioID Owens et al
bioid = 
  read_tsv(paste0(dat.dir,"GID4_bioid_combined.tsv")) %>%
  filter(GID4_SP_max>.59)

bioid_genes =
  bioid %>%
  filter(hgnc_symbol %in% top_genes$Gene) %>%
  rename(Gene=hgnc_symbol)

top_genes %>%
  filter(Gene %in% bioid_genes$Gene)

#tibble(Gene=prot_list) %>%
#  filter(Gene %in% bioid$hgnc_symbol)

#bioid %<>% pull(ID)

## construct list to perform overlaps with ggvenndiagram
#mk = sig_prots %>% filter(prot=="MKLN1") %>% pull(ID) %>% unique()
#wd = sig_prots %>% filter(prot=="WDR26") %>% pull(ID) %>% unique()

#sig_list = 
#  list(
#    bioid,mk,wd
#  )

#names(sig_list) = 
#  c(
#    "GID4 BioID HiConf",
#    "MKLN1",
#    "WDR26"
#  )

## plot the venn for all
#ggVennDiagram(sig_list)

#ggsave(filename = paste0(plot_dir2,"Combined_GID4_BioID_venn_sig_prots_overlaps.pdf"),
#       width=5,height = 4)




## GID4 + PFI-7 proteomes Owens et al
gid4_proteome = 
  read_tsv("C:/Users/imnot/OneDrive/Postdoc/Projects/GID4/Paper/Bioinformatics/Proteomes/Supp_table_5.txt") %>%
  filter(grepl("GID4",comparison))

gid4_proteome_genes =
  gid4_proteome %>%
  filter(hgnc_symbol %in% top_genes$Gene)

top_genes %>%
  filter(Gene %in% gid4_proteome_genes$hgnc_symbol)


#tibble(Gene=prot_list) %>%
#  filter(Gene %in% gid4_proteome$hgnc_symbol)


## construct list to overlap
#sig_list = 
#  list(
#    gid4_proteome,mk,wd
#  )

#names(sig_list) = 
#  c(
#    "GID4 proteome",
#    "MKLN1",
#    "WDR26"
#  )

## plot the venn for all
#ggVennDiagram(sig_list)

#ggsave(filename = paste0(plot_dir2,"Combined_GID4_proteome_venn_sig_prots_overlaps.pdf"),
#       width=5,height = 4)


pfi7_proteome =
  read_tsv("C:/Users/imnot/OneDrive/Postdoc/Projects/GID4/Paper/Bioinformatics/Proteomes/Supp_table_5.txt") %>%
  filter(grepl("PFI7",comparison))

pfi7_proteome_genes =
  pfi7_proteome %>%
  filter(hgnc_symbol %in% top_genes$Gene)

top_genes %>%
  filter(Gene %in% pfi7_proteome_genes$hgnc_symbol)


#tibble(Gene=prot_list) %>%
#  filter(Gene %in% pfi7_proteome$hgnc_symbol)

## construct list to overlap
#sig_list = 
#  list(
#    pfi7_proteome,mk,wd
#  )

#names(sig_list) = 
#  c(
#    "PFI7 proteome",
#    "MKLN1",
#    "WDR26"
#  )

## plot the venn for all
#ggVennDiagram(sig_list)

#ggsave(filename = paste0(plot_dir2,"Combined_PFI7_proteome_venn_sig_prots_overlaps.pdf"),
#       width=5,height = 4)



## RanBPM AP-MS Onea et al 2022
apms = read_tsv(paste0(dat.dir,"Onea_HCIP_trimmed2.txt"),col_names=c("Gene","Loc"))


apms_genes =
  apms %>%
  filter(Gene %in% top_genes$Gene)

top_genes %>%
  filter(Gene %in% apms_genes$Gene)


#tibble(Gene=prot_list) %>%
#  filter(Gene %in% apms$Gene)



#############################################
#### INTERSECT PREVIOUS CTLH INTERACTORS ####
#############################################

# load the full interactome database from interactR 2.0
ints = readRDS("C:/Users/imnot/OneDrive/Postdoc/Projects/interactR/data/InteractR_interactome_full.Rds")

## find CTLH interactors in interactR
CTLH_ints_main =
  ints %>%
  filter(A %in% CTLH_symbols | B %in% CTLH_symbols)

CTLH_ints = unique(c(CTLH_ints_main$A, CTLH_ints_main$B))
CTLH_ints = CTLH_ints[!CTLH_ints%in%CTLH_symbols]


CTLH_ints_genes =
  CTLH_ints_main %>%
  #filter(A %in% top_genes$Gene | B %in% top_genes$Gene) %>%
  filter(!(A %in% CTLH_symbols & B %in% CTLH_symbols)) %>%
  dplyr::select(A,B) %>%
  gather(key,Gene,A:B) %>%
  dplyr::select(-key) %>%
  filter(!Gene%in%CTLH_symbols) %>%
  filter(Gene %in% top_genes$Gene) %>%
  distinct()


top_genes %>%
  filter(Gene %in% CTLH_ints_genes$Gene)

#tibble(Gene=prot_list) %>%
#  filter(Gene %in% CTLH_ints_genes$Gene)


## search for a specific gene interacting with CTLH
to_search = "CKAP4"
CTLH_ints_main %>%
  filter(A==to_search | B==to_search)




################################
#### LOOK FOR PRO N DEGRONS ####
################################


proteome = readRDS(paste0(dat.dir,"2021_05_14_Homo_sapiens_UP000005640_degron_marked.Rds"))

perfect_degron_IDs =
  proteome %>%
  dplyr::select(ID,Degron) %>%
  filter(grepl("Perfect",Degron))

unique(proteome$Degron)


# find the increased proteins that have a perfect or perfect if trimmed degron
sig_all %>%
  filter(logFC<0) %>%
  left_join(proteome,by="ID") %>%
  filter(Degron=="Perfect") %>%
  write_tsv(paste0(stat_dir,'ProN-degron-increased-proteins.tsv'))
