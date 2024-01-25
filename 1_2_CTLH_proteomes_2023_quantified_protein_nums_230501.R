## this script takes filtered LONG format tables from
## 1_CTLH_proteomes_2023

## performs QC plot of number of quantified proteins in each sample


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

# specify the data folder for protein tables
dat.dir = paste0(base,"data/")


# output of fragpipe combined_protein.tsv file
## give it a nice name as outputs will be based on it
in_files =
  c(
    "combined_protein_mkln_wdr_msbooster_long.tsv",
    "combined_protein_rmd_maea_msbooster_long.tsv"
  )


in_files = paste0(dat.dir,in_files)


#################
#### outputs ####
#################

plot_dir = paste0(base,"plots/")
dir.create(plot_dir,showWarnings = F)


#################################
#### LOAD PREVIOUS STEP DATA ####
#################################


## combining into one table
long = tibble()
for(in_file in in_files){
  
  long %<>% rbind.data.frame(read_tsv(in_file))
  
}


##################
#### PLOTTING ####
##################


## plot number of proteins in each sample
counts =
  long %>%
  group_by(sample,group) %>%
  dplyr::count() %>%
  separate(group,into=c("Protein","condition"))


counts$Protein=factor(counts$Protein,levels=c("RMND5A","MAEA","MKLN1","WDR26"))

ggplot(counts,aes(Protein,n,colour=condition))+
  stat_summary(fun = mean, geom = "bar",position=position_dodge(.7),width=.7,size=.1,fill=NA) +
  geom_point(position=position_dodge2(width=.5),shape=21,size=.5,stroke=.2,fill=NA)+
  scale_fill_manual(values=c("#377eb8ff","#e4211cff"))+
  scale_colour_manual(values=c("#377eb8ff","#e4211cff"))+
  ylab("Quantified proteins count")+
  xlab(NULL)+
  theme_bw()+
  theme(
    strip.background = element_rect(fill=NA,size=.1),
    axis.ticks=element_line(size=.1),
    axis.ticks.x=element_blank(),
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

ggsave(filename=paste0(plot_dir,"QC_quantified_proteins.pdf"),
       width=2.2,
       height=1.1)
