## this script takes limma output tables from 
## 2_CTLH_proteomes_2023_limma_stats_final_230501

## looks for pathway changes ####### not finished

library(tidyverse)
library(magrittr)
library(pheatmap)
library(biomaRt)

# install and load the 'clusterProfiler' package
BiocManager::install("clusterProfiler")
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
    #"RMND5A",
    #"MAEA"
  )


path_out = tibble()
p=prots[1]
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


## plotting function
path_out %>%
  group_by(group) %>%              # Group data by 'group' column
  arrange(p.adjust, .by_group = TRUE) %>%   # Sort data by 'p_value' column within each group
  slice_head(n = 10) %>%
  ungroup() %>%                    # Ungroup the data
  arrange(GeneRatio) %>%  
  ggplot(aes(GeneRatio,Description,fill=p.adjust,size=Count))+
  geom_point()+
  facet_wrap(~group,scales="free")


path_out %>% filter(grepl("embryonic",Description))

## explore and combination plot




















## check whether proteins are cytoplasmic or nuclear

## use biomart for this

## get all the sig_prots IDs
#IDs = unique(sig_prots$ID)


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


}
  
  
  
  
  

ensembl <- useMart("ensembl")
ensembl_hsapiens <- useDataset("hsapiens_gene_ensembl", mart = ensembl)

#listAttributes(ensembl_hsapiens)

## get all GO terms associated with significant uniprot IDs
mapping =
  getBM(attributes = c("uniprotswissprot", "hgnc_symbol", "go_id"), 
        filters = "uniprotswissprot", 
        values = IDs, 
        mart = ensembl_hsapiens)

## alternative to biomart for mapping gene to GO ids
#gene_to_go <- select(org.Hs.eg.db, keys="ENSG00000198901", columns="GO", keytype="ENSEMBL")

## now map the GO:id to GO term name
## get the relevant GO term IDs for the main nuclear or cytoplasmic/non-nuclear GO terms

nuc_terms = c("nucleus","nucleolus","nucleoplasm","chromatin")
cyt_terms = c("cytosol","cytoplasm","plasma membrane","mitochondrion","endoplasmic reticulum","Golgi apparatus","cytoskeleton")

go_to_term = 
  select(GO.db, keys=mapping$go_id, columns="TERM") %>%
  as_tibble() %>%
  filter(TERM %in% c(nuc_terms,cyt_terms)) %>%
  distinct() %>%
  rename(go_id=GOID) %>%
  mutate(loc_group = case_when(
    TERM %in% nuc_terms ~ "nucleus",
    TERM %in% cyt_terms ~ "cytoplasm"
  ))


## now loop through each protein
ids = unique(mapping$uniprotswissprot)
ID=ids[1]
out.df=tibble()
for(ID in ids){
  
  
  cat(ID,"\n")
  
  # get just this ID
  temp=
    mapping %>%
    filter(uniprotswissprot==!!ID)
  
  ## find out if any of these go terms are the predefined location ones
  loc_group = 
    left_join(temp,go_to_term,by="go_id") %>%
    drop_na() %>%
    dplyr::select(loc_group) %>%
    distinct() %>%
    pull(loc_group) %>%
    paste0(collapse="_")
  
  out.df %<>% rbind.data.frame(tibble(ID,loc_group))
}


### now create a mapping table so that each ID gets one localisation value
id_loc =
  out.df %>%
  mutate(localization=case_when(
    loc_group==""~"neither",
    loc_group=="cytoplasm"~"cytoplasm",
    (loc_group=="cytoplasm_nucleus"|loc_group=="nucleus_cytoplasm")~"both",
    loc_group=="nucleus"~"nucleus"
  )) %>%
  dplyr::select(1,3)


## now plot this info with the directional and KO info

ud =
  sig_prots %>%
  mutate(up_down = if_else(logFC<0,"Increased","Decreased")) %>%
  left_join(id_loc) %>%
  mutate(localization=if_else(is.na(localization),"neither",localization)) %>%
  group_by(contrast,up_down,localization) %>%
  dplyr::count() %>%
  separate(contrast,into="Knockout",sep="_")

ud$localization=factor(ud$localization,levels=c("cytoplasm","nucleus","both","neither"))
ud$Knockout=factor(ud$Knockout,levels=c("RMND5A","MAEA","WDR26","MKLN1"))

ggplot(ud,aes(localization,n,fill=up_down))+
  geom_bar(stat="identity",position=position_dodge2())+
  facet_grid(~Knockout)+
  labs(fill="Protein",y="Number of significant proteins")+
  scale_fill_manual(values=c("#878787ff","#da0017ff"))+
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.border = element_rect(size=.1),
        strip.background = element_blank(),
        axis.ticks = element_blank(),
        plot.margin = unit(c(0,0,0,0), 'cm'),
        legend.key.size = unit(.2, 'cm'),
        text=element_text(size=6,colour="black"),
        axis.text.x = element_text(angle=45,hjust=1,vjust=1))






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
    filter(name %in% sig_prots$ID) %>%
    dplyr::select(name,ID,contains(!!this_prot))
  
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

