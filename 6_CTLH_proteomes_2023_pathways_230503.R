## this script takes limma output tables from 
## 2_CTLH_proteomes_2023_limma_stats_230501

## looks for pathway changes

library(tidyverse)
library(magrittr)
library(pheatmap)
library(biomaRt)
library(clusterProfiler)
library(org.Hs.eg.db)

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



##################
#### PATHWAYS ####
##################

## find the full set of significant proteins to plot
## nice to have saved
sig_prots_all = tibble()
i=1
for(in_file in in_files){
  
  # load limma table
  ## filter this table on just the sifnificant proteins FOR THIS ONE PROTEIN
  sig_prots =
    read_tsv(paste0(stat_dir,in_file)) %>%
    separate(contrast,into="prot")
  
  sig_prots_all %<>% rbind.data.frame(sig_prots)
  
}





prots =
  c(
    "WDR26",
    "MKLN1"
  )


path_out = tibble()
p=prots[2]
for(p in prots){
  
  # filter on just this prot
  temp =
    sig_prots_all %>%
    filter(prot==!!p)
  
  # get just this prots significant ones
  temp_sig =
    temp %>%
    mutate(abs_fc = abs(logFC)) %>%
    filter(adj.P.Val<p_thresh & abs_fc>lfc_thresh)
  
  # set your list of differentially expressed genes
  ## find the terms for specific direction (up or down reg)
  up =
    temp_sig %>%
    filter(logFC<0) %>%
    pull(ID)
  
  down =
    temp_sig %>%
    filter(logFC>0) %>%
    pull(ID)

  
  # the whole list of possible genes
  universe = unique(temp$ID)
  
  ## list keytypes for this db
  keytypes(org.Hs.eg.db)
  
  # run GO terms enrichment analysis using the 'enrichGO' function
  go_results_up = enrichGO(gene = up,
                        universe = universe,
                        keyType = "UNIPROT",
                        ont = "ALL",
                        pool=T,
                        pAdjustMethod = "BH",
                        OrgDb = "org.Hs.eg.db",
                        pvalueCutoff= 0.05)
  
  go_results_down = enrichGO(gene = down,
                           universe = universe,
                           keyType = "UNIPROT",
                           ont = "ALL",
                           pool=T,
                           pAdjustMethod = "BH",
                           OrgDb = "org.Hs.eg.db",
                           pvalueCutoff= 0.05)
  
  # extract the top enriched GO terms
  top_GO_terms_up <- head(go_results_up, n = 100)
  top_GO_terms_down <- head(go_results_down, n = 100)
  
  # plot the GO terms enrichment results, only if entries
  if(nrow(top_GO_terms_up)>0){
    dotplot(go_results_up)
    ggsave(paste0(plot_dir,"clusterProfiler",p,"_UP.pdf"),height=4,width=5)
    ## add to table
    temp_df =
      as_tibble(go_results_up) %>%
      mutate(group=paste0(p,"_","Increased"))
    path_out %<>% rbind.data.frame(temp_df)
  }
  
  if(nrow(top_GO_terms_down)>0){
    dotplot(go_results_down)
    ggsave(paste0(plot_dir,"clusterProfiler",p,"_DOWN.pdf"),height=4,width=5)
    ## add to table
    temp_df =
      as_tibble(go_results_down) %>%
      mutate(group=paste0(p,"_","Decreased"))
    path_out %<>% rbind.data.frame(temp_df)
  }
}
