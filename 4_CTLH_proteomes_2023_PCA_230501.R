## this script takes limma output tables from 
## 2_CTLH_proteomes_2023_limma_stats_final_230501

## performs overlaps of significant proteins between the proteins

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
  
  ## do the PCA on combined data
  pc = prcomp(t(mat), scale. = F, center=T)
  
  sig_prot_num_temp = length(unique(sig_prots$ID))
  
  ## set up values and sample labels to plot using raw ggplot
  pc_vals =
    pc$x %>%
    data.frame() %>%
    rownames_to_column("sample") %>%
    separate(sample, into=c("Protein","condition","Rep"))
  
  
  ## get the variance explained for each component
  var_explained <- pc$sdev^2/sum(pc$sdev^2)
  
  # plot PC1 and PC2
  pc_vals %>%
    ggplot(aes(x=PC1,y=PC2))+
    geom_point(size=1,aes(shape=Protein,colour=condition))+
    scale_colour_manual(values=c("#377eb8ff","#e4211cff"))+
    ggtitle("PCA plot of proteomics combined",paste0(sig_prot_num_temp," DEPs"))+
    labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"),
         y=paste0("PC2: ",round(var_explained[2]*100,1),"%")) +
    theme_bw()+
    theme(
      strip.background = element_rect(fill=NA,size=.1),
      axis.ticks=element_line(size=.1),
      #axis.ticks.x=element_blank(),
      axis.ticks.length=unit(.3,"mm"),
      panel.background = element_blank(),
      legend.key.size = unit(.2,"cm"),
      #axis.text.x = element_text(size=6,angle=45,hjust=1,vjust=1),   
      plot.title = element_text(hjust = 0.5,vjust=-1.5),
      plot.subtitle = element_text(hjust = 0.5,vjust=-.5),
      panel.border=element_rect(size=.1),
      panel.grid=element_blank(),
      text=element_text(size=6),
      strip.text.x = element_text(margin = margin(0,0,0,0, "cm")))
  
  
  ggsave(paste0(plot_dir,"PCA_combined_",this_prot,".pdf"),
         width=1.7,
         height=1.2)
  
  
  
  #pdf(paste0(plot_dir,"upset_all_DEPs_",this_prot,".pdf"),3.5,2.5)
  #print(upset_plot)
  #dev.off()
  

  #ggsave(filename = paste0(plot_dir,"venn_sig_prots_all_overlaps_",this_prot,".pdf"),
  #       width=4.5,height = 3.5)
  
}

