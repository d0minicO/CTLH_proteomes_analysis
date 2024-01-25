## this script takes limma output tables from 
## 2_CTLH_proteomes_2023_limma_stats_final_230501

## plots volcano plots off differential proteins ## highlighting MKLN1

library(tidyverse)
library(magrittr)
library(ggrastr)
library(ggrepel)
library(patchwork)

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


## the full significant proteins tables
in_files =
  c(
    "MKLN1-WDR26_msboost_limma_all.tsv",
    "RMND5A-MAEA_msboost_limma_all.tsv"
  )

## significance cutoffs
p_thresh = .05
lfc_thresh = log2(1.5)

#################
#### outputs ####
#################

plot_dir = paste0(base,"plots/")



######################
#### VOLCANO PLOT ####
######################

## get a combined stats table to plot volcano from
sig_prots_all = tibble()
for(in_file in in_files){

  temp =
    read_tsv(paste0(stat_dir,in_file)) %>%
    separate(contrast,into="prot") %>%
    mutate(abs_fc = abs(logFC))
  
  sig_prots_all %<>% rbind.data.frame(temp)
  
}


## prepare the df for volcano analysis
volc =
  sig_prots_all %>%
  mutate(minuslogp = -log10(adj.P.Val)) %>%
  arrange(desc(minuslogp))


## set up a column for significance
## can chose whether or not to include the log2 FC
volc %<>%
  mutate(sig = if_else(
    (adj.P.Val<=p_thresh & abs_fc > lfc_thresh),
    "sig",
    "notsig"
  ))


## swap the direction of the logFC to make more intuitive sense
## points will be to the right if they were
## higher in the knockout
volc %<>%
  mutate(logFC = -logFC)



## list of proteins to plot volcano for
prots =
  c(
    "RMND5A",
    "MAEA",
    "MKLN1",
    "WDR26"
  )

## for each plot, chose the proteins to label
## separate with : if multiple
to_label =
  c(
    "MKLN1:HMGCS1", # RMND5A labels
    "MKLN1:HMGCS1", # MAEA labels
    "HMGCS1", # MKLN1 labels
    "HMGCS1" # WDR26 labels
    )

to_loop =
  tibble(prots,
         to_label)



# set up universal features of the plots
boxpad = .5
pointpad = .5
minlength = .01
legpos = "bottom"


## loop throuh each protein at a time
plot_list = list()
i=1
for(i in 1:nrow(to_loop)){
  
  ## this protein to plot
  this_prot = to_loop[i,1] %>% as.character()
  ## proteins to label
  ## deal with multiple entries colon separated
  this_labels = 
    to_loop[i,2] %>%
    mutate(label_split = str_split(to_label,pattern=":")) %>%
    unnest(cols=c(label_split)) %>%
    pull(label_split)
  
  ## filter the main volc table on just this prot
  temp =
    volc %>%
    filter(prot==!!this_prot)
  
  ## set up label column to use in geom_text_repel
  temp %<>%
    mutate(label = if_else(
      Gene %in% this_labels,
      Gene,
      "nolabel"
    ))
  
  
  ## only plot a y axis for the first plot
  if(i==1){
    ylab_val = "-log10(p.adjusted)"
  } else {
    ylab_val = NULL
  }
  
  
  ## create dummy variable to get legend to look nice
  #temp %<>%
  #  mutate(allhits = "All proteins")
  
  ## get number of sig proteins (up or down) to annotate to upper RH corner of plot
  prot_nums =
    temp %>%
    filter(sig=="sig") %>%
    dplyr::select(logFC) %>%
    mutate(pos_or_neg = if_else(
      logFC>0,
      "up",
      "down"
    )) %>%
    group_by(pos_or_neg) %>%
    dplyr::count()
  
  
  ## find x and y values to position proteins labelling
  label_y = max(temp$minuslogp)*.9
  label_x_left = min(temp$logFC)*.9
  label_x_right = max(temp$logFC)*.9
  
  # custom volcano plot
  plot_list[[i]]=
    ggplot(temp,aes(logFC,minuslogp, label = label))+
    geom_hline(yintercept = -log10(p_thresh),linetype="dashed",size=.1)+
    geom_vline(xintercept = log2(lfc_thresh),linetype="dashed",size=.1)+
    geom_vline(xintercept = -log2(lfc_thresh),linetype="dashed",size=.1)+
    geom_point_rast(data=subset(temp, sig!="sig"),alpha=.4,raster.dpi=900,shape=".",colour="gray60")+ # rasterise as many points are slow
    geom_point_rast(data=subset(temp, sig=="sig"&label=="nolabel"),alpha=.8,raster.dpi=900,shape=".")+
    geom_point_rast(data=subset(temp, label!="nolabel"),alpha=.95,raster.dpi=900,shape=20,size=.1,colour="red")+
    # decreased proteins number label
    annotate(geom="text", x=label_x_left, y=label_y, label=prot_nums[1,2],color="black",size=1)+
    # increased proteins number label
    annotate(geom="text", x=label_x_right, y=label_y, label=prot_nums[2,2],color="black",size=1)+
    ggtitle(this_prot)+
    ylab(ylab_val)+
    xlab("Log2FC")+
    theme_bw()+
    theme(
      panel.grid = element_blank(),
      panel.border = element_rect(size=.1),
      axis.ticks = element_line(size=.1),
      text=element_text(size=5),
      legend.key.size = unit(5,"mm"),
      title=element_text(size=2.5),
      legend.position = "bottom"
    ) +
    geom_text_repel(data          = subset(temp, label!="nolabel"),
                    colour="black",
                    size          = 1,
                    box.padding   = boxpad,
                    point.padding = pointpad,
                    max.overlaps = 40,
                    force         = 100,
                    segment.size  = 0.1,
                    min.segment.length = minlength,
                    segment.color = "grey50",
                    #direction     = "x")
    )
  
}



wrap_plots(plot_list) +
  plot_layout(guides="collect",ncol=4) +
  plot_annotation(title = "CTLH_proteomes_volcano",
                  theme = 
                    theme(
                      plot.title = element_text(size = 5),
                      legend.position = legpos
                    ))

ggsave(paste0(plot_dir,"Volcano_all.pdf"),
       width=3.8,
       height=1.3)

