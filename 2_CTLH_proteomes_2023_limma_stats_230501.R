## this script takes filtered wide format tables from
## 1_1_CTLH_proteomes_2023

## saved environment image to load back later
img="C:/Users/imnot/OneDrive/Postdoc/Projects/Maitland_proteomes/Analysis_v3_separated_Apr23/CTLH_proteomes_workspace_230503.Rdata"
#save.image(img)
load(img)

## quickly get p values for a gene
sig_prots_main %>%
  filter(Gene=="MKLN1") # COL3A1

## performs additional filtering, imputation, and stats with limma


library(tidyverse)
library(magrittr)
library(DEP)
library(SummarizedExperiment)
library(limma)


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


## DATA FRAME WITH EXPERIMENT PARAMETERS TO LOOP OVER

# filenames of wide table (located in ./data)
in_files =
  c(
    #"combined_protein_mkln_wdr_wide.tsv",
    "combined_protein_mkln_wdr_msbooster_wide.tsv",
    #"combined_protein_rmd_maea_wide.tsv",
    "combined_protein_rmd_maea_msbooster_wide.tsv"
  )


# nice name for the experiment group / protein
proteins =
  c(
    #"MKLN1-WDR26_normal",
    "MKLN1-WDR26_msboost",
    #"RMND5A-MAEA_normal",
    "RMND5A-MAEA_msboost"
  )


to_analyze =
  tibble(in_file=in_files,
         protein=proteins)


## p threshold and logFC threshold for significant proteins
p_thresh = .05
lfc_thresh = log2(1.5)


#################
#### outputs ####
#################

plot_dir = paste0(base,"plots/")

stat_dir = paste0(base,"stats/")
dir.create(stat_dir,showWarnings = F)


# imputed wide out file
imp_out_files = gsub(".tsv$", "_imputed.tsv", in_files)



########################################
#### LOOP FOR EACH EXPERIMENT GROUP ####
########################################

## main output table from limma for all proteins / groups
sig_prots_main = tibble()
i = 2
for(i in 1:nrow(to_analyze)){
  
  # this experiment to analyze
  in_file = paste0(dat.dir,to_analyze[i,1])
  protein = as.character(to_analyze[i,2])
  
  cat("Now analyzing",protein,"\n")
  
  #################################
  #### LOAD PREVIOUS STEP DATA ####
  #################################
  
  wide = read_tsv(in_file)
  
  ## get a conversion table of gene name to ID to use later
  ## deduplicate gene names
  gene_id =
    wide %>%
    dplyr::select(ID,Gene) %>%
    distinct(Gene, .keep_all = T)
  
  #########################
  #### PREPARE FOR DEP ####
  #########################
  
  # Generate a SummarizedExperiment object using an experimental design
  LFQ_columns <- 3:ncol(wide)# get LFQ column numbers (for area / intensity output (MaxLFQ))
  
  # get experimental design using wide table column names
  exp_des = 
    tibble(label = colnames(wide)[LFQ_columns]) %>%
    separate(label,into=c("protein", "condition", "replicate"),sep="_",remove=F) %>%
    mutate(condition=paste0(protein,"_",condition)) %>%
    dplyr::select(label,condition,replicate)
  
  ## make sure all LFQ cols are numeric
  wide %<>% mutate(across(all_of(LFQ_columns),as.numeric))
  
  ## run make_unique to give necessary columns and check each protein is on its own row
  wide %<>% make_unique(names = "Gene",ids="ID")
  
  ## make the summarized experiment object
  data_se <- make_se(wide, LFQ_columns, exp_des)
  
  
  #################################
  #### FILTERING / NORMALIZING ####
  #################################
  
  ## according to Goeminne et al 2020, default DEP filtering is too stringent and can be dropped
  # https://doi.org/10.1021/acs.analchem.9b04375
  
  # Filter for proteins that are identified in all replicates of at least one condition
  #data_filt <- data_se # filter_missval(data_se, thr = 0)
  
  # Do a Less stringent filtering, though optional
  # Filter for proteins that are identified in 2 out of 4 replicates of at least one condition
  data_filt <- filter_missval(data_se, thr = 2)
  
  # no more filtering
  #data_filt <- data_se
  
  # Normalize the data
  data_filt <- normalize_vsn(data_filt)
  
  #meanSdPlot(data_filt)

  ########################
  #### DEP IMPUTATION ####
  ########################
  
  # set seed to have imputation be reproducible
  set.seed(123)
  
  ## performing mixed imputation as per DEP vignette and Goeminne et al 2020
  proteins_MNAR <- get_df_long(data_filt) %>%
    group_by(name, condition) %>%
    summarize(NAs = all(is.na(intensity))) %>% 
    filter(NAs) %>% 
    pull(name) %>% 
    unique()
  
  # Get a logical vector
  MNAR <- names(data_filt) %in% proteins_MNAR
  
  # Perform a mixed imputation
  set.seed(123)
  mixed_imputation <- impute(
    data_filt, 
    fun = "mixed",
    randna = !MNAR, # we have to define MAR which is the opposite of MNAR
    mar = "knn", # imputation function for MAR
    mnar = "QRILC") # imputation function for MNAR
  
  ## compare imputation options with a plot
  plot_imputation(data_filt,mixed_imputation)
  
  ## paste in a variable name here
  ggsave(paste0(plot_dir,"QC_imputations_",protein,".pdf"),
         width=6,
         height=6)
  
  
  ### export imputed wide tables for use later
  get_df_wide(mixed_imputation) %>%
    write.table(file = paste0(dat.dir,imp_out_files[i]),
                row.names = F,
                col.names = T,
                quote=F,
                sep="\t")
  
    
  
  #################################################
  #### TEST DIFFERENTIAL ABUNDANCE USING LIMMA ####
  #################################################
  
  ## make the design matrix for limma, derived from sample info in col names
  design =
    tibble(coldat = colnames(mixed_imputation)) %>%
    separate(coldat, into=c("prot","condition","rep"),remove=F) %>%
    mutate(group=paste0(prot,"_",condition)) %>%
    mutate(val=1) %>%
    spread(key=group,value=val) %>%
    replace(is.na(.),0) %>%
    dplyr::select(-coldat,-condition,-rep,-prot)

  
  # contrasts to use for limma (based on column names)
  ## contrasts got using the protein type
  if(grepl("MKLN1",protein)){
    contrast = c(
      "MKLN1_control-MKLN1_knockout",
      "WDR26_control-WDR26_knockout"
    )
  } else if(grepl("RMND5A",protein)){
    contrast = c(
      "RMND5A_control-RMND5A_knockout",
      "MAEA_control-MAEA_knockout"
    )
  }
  
  
  ## limma linear regression on mixed imputed data
  ## get imputed data as numeric matrix
  mat = assays(mixed_imputation)[[1]]
  fit = lmFit(mat, design)
  
  # get p values for each contrast
  sig_prots = tibble()
  #cont=contrasts[1]
  for(cont in contrast){
    
    cat("Testing contrast ",cont,"\n")
    
    # set up the contrasts
    cont.matrix = makeContrasts(contrasts = cont, levels=design)
    
    fit2 = contrasts.fit(fit, cont.matrix)
    
    # correct using empirical bayes fit
    fit2 = eBayes(fit2)
    
    # extract all proteins with BH adjusted p values and log2FC
    sig_prots_temp =
      topTable(fit2, adjust="none",n=Inf) %>%
      as.data.frame() %>%
      rownames_to_column("ID") %>%
      as_tibble() %>%
      mutate(contrast=cont)
    
    sig_prots %<>% rbind.data.frame(sig_prots_temp)
    
  }
  
  # adjust p values for multiple comparisons (Benjamini & Hochberg)
  sig_prots %<>%
    mutate(adj.P.Val = p.adjust(P.Value,method="BH"))
  
  
  ## add back the gene names to the table
  ## need to deduplicatea some gene names to allow matching to IDs properly
  sig_prots %<>%
    mutate(Gene=ID) %>%
    separate(Gene,into=c("Gene2","extra"),sep="\\.") %>%
    mutate(Gene=Gene2) %>%
    dplyr::select(-Gene2,-extra,-ID) %>%
    left_join(gene_id)
  
  ## format nicer for export
  sig_prots %<>%
    dplyr::select(Gene,ID,1:7)
  
  
  # export the limma results for this experiment group
  ## keep all the proteins, useful for checking p values of specific proteins later
  write.table(sig_prots,
              file = paste0(stat_dir,protein,"_limma_all.tsv"),
              row.names = F,
              col.names = T,
              quote=F,
              sep="\t")
  
  ## add to the big main for proteins table
  sig_prots_main %<>% rbind.data.frame(sig_prots)
  
  ## tidy up the table by filtering just on significant proteins
  sig_prots %<>%
    mutate(abs_fc = abs(logFC)) %>%
    filter(adj.P.Val<p_thresh & abs_fc>lfc_thresh)
  
  # export the significant proteins
  write.table(sig_prots,
              file = paste0(stat_dir,protein,"_limma_significant.tsv"),
              row.names = F,
              col.names = T,
              quote=F,
              sep="\t")
  
  
}


## how many significant proteins in each comparison?
sig_prots_main %>%
  mutate(abs_fc = abs(logFC)) %>%
  filter(adj.P.Val<p_thresh & abs_fc>lfc_thresh) %>%
  group_by(contrast) %>%
  dplyr::count()


## what about up and down?
sig_prots_main %>%
  mutate(abs_fc = abs(logFC)) %>%
  mutate(updown=if_else(logFC>0,"down","up")) %>%
  filter(adj.P.Val<p_thresh & abs_fc>lfc_thresh) %>%
  group_by(contrast,updown) %>%
  dplyr::count()

