## this script takes limma output tables from 
## 8_CTLH_proteomes_2023_median_centering_230503

## plots MaxLFQ intensities, imputed intensites, and median centred values
## for a specific set of proteins

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

# data with full long table
dat.dir = paste0(base,"data/")

# long combined table
long_comb_out_file = paste0(dat.dir,"Combined_protein_measurements_long.tsv")


# specify the folder for limma output tables
stat_dir = paste0(base,"stats/")
# filenames of limma output tables
limma_files = list.files(stat_dir,pattern = "limma_significant.tsv")


#################
#### outputs ####
#################

plot_dir = paste0(base,"plots/")

gene_dir = paste0(plot_dir,"genes/")
dir.create(gene_dir,showWarnings = F)


###################
#### LOAD DATA ####
###################

long = read_tsv(long_comb_out_file)

## get unique gene IDs mapping table
gene_ids =
  long %>%
  dplyr::select(ID,Gene) %>%
  distinct()



######################
#### GENE PLOTTER ####
######################


## find the full set of significant proteins to plot
## nice to have saved
sig_prots_all = tibble()
i=1
for(lim_file in limma_files){
  
  # load limma table
  ## filter this table on just the sifnificant proteins FOR THIS ONE PROTEIN
  sig_prots =
    read_tsv(paste0(stat_dir,lim_file)) %>%
    separate(contrast,into="prot")
  
  sig_prots_all %<>% rbind.data.frame(sig_prots)
  
}




#######
## 1 ##
#######


## plot each DEP in a single plot loop with the four different values possible
## with significance shown

prot_list = c(
  unique(sig_prots_all$Gene)
)


## get list of common chaned proteins in WDR26 and MKLN1
mkln =
  sig_prots_all %>%
  filter(prot=="MKLN1") %>%
  pull(Gene)

wdr26 =
  sig_prots_all %>%
  filter(prot=="WDR26") %>%
  pull(Gene)

## in common
prot_list = intersect(mkln,wdr26)

# Find items NOT shared between the vectors
unique_to_mkln <- setdiff(mkln, wdr26)
unique_to_wdr26 <- setdiff(wdr26, mkln)

# Combine the unique items into one vector
not_shared_items <- c(unique_to_mkln, unique_to_wdr26)
prot_list = not_shared_items

## find the opposite direction proteins significant in both
mkln_up =
  sig_prots_all %>%
  filter(prot=="MKLN1" & logFC<0) %>%
  pull(Gene)

mkln_down =
  sig_prots_all %>%
  filter(prot=="MKLN1" & logFC>0) %>%
  pull(Gene)

wdr26_up =
  sig_prots_all %>%
  filter(prot=="WDR26" & logFC<0) %>%
  pull(Gene)

wdr26_down =
  sig_prots_all %>%
  filter(prot=="WDR26" & logFC>0) %>%
  pull(Gene)

## in common
prot_list = intersect(mkln_up,wdr26_down)
prot_list = intersect(mkln_down,wdr26_up)


#prot_list = c("ANPEP","LRRC15", "FN1","SNCG", "LPL", "CPS1", "AKRC12") # c("PKM","ENO1")  

#this_prot = prot_list[1]

#prot_list = c("MKLN1","HMGCS1")
#this_prot = "MKLN1"

for(this_prot in prot_list){
  
  cat(this_prot,"\n")
  
  temp =
    long %>%
    filter(Gene==!!this_prot)
  
  
  ## find which prots were significant to add this info visually to the plot
  sig_temp =
    sig_prots_all %>%
    filter(Gene==!!this_prot) %>%
    dplyr::select(ID,Gene,prot) %>%
    mutate(sig="Sig.")
  
  ## join the sig temp info to the temp plotting df
  temp %<>%
    left_join(sig_temp,by=c("ID", "Gene", "prot"))
  
  
  ## get into long format over the different quantification types
  temp %<>%
    gather(quant_type,value,LFQ:centd)
  
  
  ## then fix the NAs in sig
  temp %<>%
    ungroup %>%
    mutate(sig=if_else(
      is.na(sig),
      "N.S.",
      sig))
  
  ## set up group for the plot
  temp %<>%
    separate(condition,into=c("prot","condition","rep")) %>%
    mutate(group=paste0(prot,"_",condition))
  
  ## set up a variable for the plot height based on how many exp_groups it was detected in
  ## not all proteins detected in both!!
  if(length(unique(temp$prot))==4){
    plot_height = 2
  } else {
    plot_height = 1.1
  }
  
  ## plot log2(LFQ), imputed LFQ, and median centered log2 imputed LFQ in one plot
  p =
    ggplot(temp,aes(group,value))+
    geom_boxplot(size=.1,outlier.shape = NA,aes(colour=condition))+
    geom_point(shape=21,size=.5,stroke=0,position=position_jitterdodge(.6), aes(fill=sig))+
    scale_colour_manual(values=c("#377eb8ff","#e4211cff"))+
    scale_fill_manual(values=c("gray60","red"))+
    ggtitle(this_prot)+
    facet_wrap(prot~quant_type,scales="free",ncol=6)+
    xlab(NULL)+
    ylab(NULL)+
    theme_bw()+
    theme(
      strip.background = element_rect(fill=NA,size=.1),
      strip.text = element_text(margin = margin(.1,.1,.1,.1, "mm"),size=4),
      axis.ticks=element_line(size=.1),
      axis.ticks.length = unit(.3,"mm"),
      panel.background = element_blank(),
      legend.key.size = unit(2,"mm"),
      legend.title=element_blank(),
      axis.text.x = element_blank(),#(size=5,angle=45,hjust=1,vjust=1),   
      plot.title = element_text(hjust = 0.5,vjust=-1.5),
      plot.subtitle = element_text(hjust = 0.5,vjust=-.5),
      panel.border=element_rect(size=.1),
      panel.grid=element_blank(),
      text=element_text(size=5))
  
  ggsave(p,
         filename = paste0(gene_dir,this_prot,".pdf"),
         width=4,
         height=plot_height)
  
  
}





#######
## 2 ##
#######

## nicer smaller plots of just median centered values, selected proteins


prot_list =  c("CKAP4","HMGCS1") #c("PKM","ENO1") #"LIMA1" #"CSNK2A2" #c("ANPEP","BAIAP2","MYO6","TLE") #c("MKLN1","HMGCS1")
#this_prot = "MKLN1"

for(this_prot in prot_list){
  
  cat(this_prot,"\n")
  
  temp =
    long %>%
    filter(Gene==!!this_prot)
  
  ## set up group for the plot
  temp %<>%
    separate(condition,into=c("prot","condition","rep")) %>%
    mutate(group=paste0(prot,"_",condition))
  
  ## set up an exp_group for the plot
  temp %<>%
    mutate(exp_group=if_else(prot %in% c("RMND5A","MAEA"),"RMND5A-MAEA","MKLN1-WRD26"))
  
  ## set up a variable for the plot height based on how many exp_groups it was detected in
  ## not all proteins detected in both!!
  if(length(unique(temp$prot))==4){
    plot_width = 2.3
  } else {
    plot_width = 1.5
  }
  
  temp$prot=factor(temp$prot,levels=c("RMND5A","MAEA","MKLN1","WDR26"))
  
  ## plot log2(LFQ), imputed LFQ, and median centered log2 imputed LFQ in one plot
  p =
    ggplot(temp,aes(prot,centd,fill=condition))+
    geom_boxplot(size=.1,outlier.shape = NA,colour="black",position=position_dodge2(width=.6))+
    geom_point(shape=21,size=.7,stroke=.2,position=position_dodge2(.6),fill=NA,aes(colour=condition))+
    scale_colour_manual(values=c("#377eb8ff","#e4211cff"))+
    scale_fill_manual(values=c("#377eb8ff","#e4211cff"))+
    facet_wrap(~exp_group,scales="free",ncol=6)+
    ggtitle(this_prot)+
    xlab(NULL)+
    ylab(NULL)+
    theme_bw()+
    theme(
      strip.background = element_rect(fill=NA,size=.1),
      strip.text = element_text(margin = margin(.1,.1,.1,.1, "mm"),size=4),
      axis.ticks=element_line(size=.1),
      axis.ticks.length = unit(.3,"mm"),
      panel.background = element_blank(),
      legend.key.size = unit(2,"mm"),
      legend.title=element_blank(),
      axis.text.x = element_text(size=4),   
      plot.title = element_text(hjust = 0.5,vjust=-1.5),
      plot.subtitle = element_text(hjust = 0.5,vjust=-.5),
      panel.border=element_rect(size=.1),
      panel.grid=element_blank(),
      text=element_text(size=5))
  
  ggsave(p,
         filename = paste0(gene_dir,this_prot,"_cleaner.pdf"),
         width=plot_width,
         height=1.1)
  
  
}

