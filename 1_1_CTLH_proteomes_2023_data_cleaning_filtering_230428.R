## this script takes a fragpipe combined_protein.tsv file
## and saves cleanded and filtered long format and wide format tables


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
dir.create(dat.dir,showWarnings = F)

# output of fragpipe combined_protein.tsv file
## give it a nice name as outputs will be based on it
in_files =
  c(
    "combined_protein_mkln_wdr.tsv",
    "combined_protein_mkln_wdr_msbooster.tsv",
    "combined_protein_rmd_maea.tsv",
    "combined_protein_rmd_maea_msbooster.tsv"
  )

f="combined_protein_mkln_wdr.tsv"

data_prepper = function(dat.dir,f){
  ## function that takes the data directory where input / output files go
  ## all files needed by the script must be in here
  
  # inputs (in dat.dir)
  
  # 1 # combined_protein.tsv files from Fragpipe
  # 2 # fragpipe-files.fp-manifest file from Fragpipe
  # 3 # maxquant contaminants file contaminants.fasta
  # 4 # Hao contaminants file Hao_contaminants.tsv
  
  # outputs  (saved into dat.dir)
  
  # 1 # a combined_protein_long.tsv file for each in_file
  # 2 # a combined_protein_wide.tsv file for each in_file
  
  
  ### SCRIPT START ###
  
  cat("exporting tables for", f, "\n")
  
  ################
  #### inputs ####
  ################
  # these must be in the data directory (dat.dir)
  
  # the input file
  in_file = paste0(dat.dir,f)
  
  # sample mappings from the fragpipe manifest
  sample_info_file = paste0(dat.dir,"fragpipe-files.fp-manifest")
  
  ## maxquant contaminants fasta db file
  con_file_mq = paste0(dat.dir,"contaminants.fasta")
  
  # Hao contaminants file ## https://pubs.acs.org/doi/10.1021/acs.jproteome.2c00145
  con_file = paste0(dat.dir,"Hao_contaminants.tsv")
  
  
  
  #################
  #### outputs ####
  #################
  
  # long out file
  long_out_file = gsub(".tsv$", "_long.tsv", in_file)
  
  # wide out file
  wide_out_file = gsub(".tsv$", "_wide.tsv", in_file)
  
  
  
  ###################
  #### DATA LOAD ####
  ###################
  
  ### FRAGPIPE DATA ###
  dat = read_tsv(in_file)
  
  
  # read the annotated likely contaminants Hao group
  con_dat = 
    read_tsv(con_file) %>%
    dplyr::rename(ID=`Uniprot ID`,Gene=`Gene names`,Source=`Source of Contamination`) %>%
    mutate(con="CON",db="Hao") %>%
    dplyr::select(ID,Gene,Source,con,db)
  
  
  # Read in fasta file MaxQuant contaminants DB and clean up
  fasta = readLines(con_file_mq)
  
  # Extract first part of each entry
  ids = grep("^>", fasta, value = TRUE)
  ids = sub("^>(\\S+).*", "\\1", ids)
  
  con_dat_mq = 
    tibble(ID=ids) %>%
    mutate(Gene=NA,Source=NA) %>%
    mutate(con="CON",db="MQ")
  
  
  ## combine Hao and MQ contaminants into a single DB
  con_dat %<>% rbind.data.frame(con_dat_mq)
  
  
  
  ############################
  #### FILTERING PROTEINS ####
  ############################
  
  ## only keep human proteins
  dat %<>%
    filter(Organism=="Homo sapiens OX=9606")
  
  ## see distribution of proteins based on total number of peptides and unique spectral counts
  dat %>%
    mutate(protein_group = case_when(
      `Combined Total Peptides`>1&`Combined Unique Spectral Count`<2 ~ ">1_pep::<2_Unique.SPC",
      `Combined Total Peptides`<2&`Combined Unique Spectral Count`>1 ~ "<2_pep::>1_Unique.SPC",
      `Combined Total Peptides`>1&`Combined Unique Spectral Count`>1 ~ ">1_pep::>1_Unique.SPC",
      `Combined Total Peptides`<2&`Combined Unique Spectral Count`<2 ~ "<2_pep::<2_Unique.SPC"
    )) %>%
    group_by(protein_group) %>%
    dplyr::count()
  
  
  ## only keep proteins with >1 peptide and >1 unique spectral count
  dat %<>%
    filter(`Combined Total Peptides`>1&`Combined Unique Spectral Count`>1)
  
  dat %<>%
    dplyr::select(2,4,contains("LFQ")) %>%
    dplyr::rename(ID=`Protein ID`)
  
  
  
  #############################
  #### REMOVE CONTAMINANTS ####
  #############################
  
  ## Explore number of contaminants in the data that were annotated by either Hao & MQ
  dat %>%
    inner_join(con_dat,by="ID") %>%
    group_by(db) %>%
    dplyr::count()
  
  
  ## remove the contaminants
  dat %<>%
    left_join(con_dat,by="ID") %>%
    filter(is.na(con))
  
  
  ##########################
  #### LONG FORMAT DATA ####
  ##########################
  
  
  ## part 1 ##
  ## continued in part 2 after saving wide again ##
  
  ## gather wide format programmatically to get long data
  keycol <- "sample"
  valuecol <- "LFQ"
  gathercols <- colnames(dat)[3:ncol(dat)]
  
  long = gather_(dat, keycol, valuecol, gathercols)
  
  
  ## clean up cols and match to samples
  long$sample = gsub(" MaxLFQ Intensity","", long$sample)
  
  # get a new col name to join to the sample IDs later
  long %<>%
    dplyr::rename(sample_ID = sample)
  
  # parse the sample info from FragPipe manifest
  sample_info = 
    read_tsv(sample_info_file,col_names = F) %>%
    mutate(sample = paste0(X2,"_",X3)) %>%
    dplyr::select(X1,sample) %>%
    mutate(X1=gsub("/scratch/dowens/CTLH_proteomes/mzml/","",X1)) %>%
    mutate(X1=gsub(".mzML","",X1))
  
  ## use separate to get more useful metadata cols
  ## cleaning up needed for BGC samples as naming is different
  sample_info %<>%
    mutate(X1 = if_else(
      grepl("BGC",X1),
      gsub("KO","KO_",X1),
      X1
    )) %>%
    mutate(X1 = if_else(
      grepl("BGC",X1),
      gsub("WT","WT_",X1),
      X1
    ))
  
  ## now separate sample name to get useful info
  sample_info %<>%
    separate(X1, into=c(NA,"date","protein","type","rep"),sep="_") %>%
    mutate(rep=parse_number(rep),
           date=gsub("Jan","01",date),
           date=gsub("Mar","03",date),
           protein=gsub("BGC","MKLN1",protein),
           protein=gsub("untreated","",protein),
           protein=gsub("WDR","WDR26",protein)
    ) %>%## now add whether it is a KO or control column
    mutate(sampleType = case_when(
      type %in% c("EV","WT") ~ "control",
      type=="KO" ~ "knockout"
    )) %>%
    dplyr::rename(sample_ID = sample) %>%
    mutate(sample=paste(protein,sampleType,rep,sep="_"))
  
  
  ## now join sample info to the long data table
  long %<>% 
    left_join(sample_info) %>%
    dplyr::rename(Gene=Gene.x)
  
  ## when technical replicates included that were run on different dates
  ## this is causing duplicate protein measurements to be included
  ## quick solution is to drop date and make distinct
  long %<>%
    dplyr::select(-date) %>%
    distinct()
  
  
  ##########################
  #### WIDE FORMAT DATA ####
  ##########################
  
  ## keep long table but also spread back into wide table for supplementary tables / export
  ## also make LFQ numeric
  wide = 
    long %>%
    dplyr::select(ID,Gene,sample,LFQ) %>%
    distinct() %>%
    mutate(LFQ=as.numeric(LFQ))
  
  # spread back to long format
  wide = spread(data=wide,key=sample,value=LFQ)
  
  # only keep rows with atleast one non-zero measurement (in a numeric col)
  wide %<>%
    filter_if(is.numeric, any_vars(. != 0)) %>%
    # mystery extra column to be deleted
    dplyr::select(-ncol(wide))
  
  
  ##########################
  #### LONG FORMAT DATA ####
  ##########################
  
  ## part 2 ##
  
  ## create a grouping variable
  long %<>%
    mutate(group=paste0(protein,"_",sampleType))
  
  ## filter out the zero quants to make table smaller
  long %<>%
    filter(LFQ>0)
  
  
  
  
  ######################
  #### SAVE OUTPUTS ####
  ######################
  
  
  ## outputs are wide and long format cleaned data
  
  wide %>%
    write.table(wide_out_file,
                quote=F,
                col.names = T,
                row.names = F,
                sep="\t")
  
  
  
  long %>%
    write.table(long_out_file,
                quote=F,
                col.names = T,
                row.names = F,
                sep="\t")
  
}


# run function on input files I have
for(f in in_files){
  
  data_prepper(dat.dir,f)
  
}