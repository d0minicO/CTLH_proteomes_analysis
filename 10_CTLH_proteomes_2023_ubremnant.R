## this script takes combined_modified_peptide.tsv from FragPipe
## performs cleaning, filtering and imputation
## exports QC plots and logFC plot


library(tidyverse)
library(magrittr)
library(DEP)
library(SummarizedExperiment)
library(limma)
library(ggrepel)
library(ggrastr)

#################
#### OPTIONS ####
#################

options(max.print=100)


################
#### inputs ####
################

# specify root folder that contains the data
base = "C:/Users/imnot/OneDrive/Postdoc/Projects/Maitland_proteomes/Analysis_v3_separated_Apr23/"


dat.dir = paste0(base,"data/")

## significance cutoffs
p_thresh = .05
lfc_thresh = log2(1.5)

# files must be in the data directory (dat.dir)

# the input file
in_file = paste0(dat.dir,'combined_modified_peptide_digly.tsv')

# sample mappings from the fragpipe manifest
sample_info_file = paste0(dat.dir,"experiment_annotation.tsv")

## maxquant contaminants fasta db file
con_file_mq = paste0(dat.dir,"contaminants.fasta")

# Hao contaminants file ## https://pubs.acs.org/doi/10.1021/acs.jproteome.2c00145
con_file = paste0(dat.dir,"Hao_contaminants.tsv")


# location of whole proteome fold changes from limma for RMND5A KO (located in ./data)
prot_files = "RMND5A-MAEA_msboost_limma_all.tsv"





#################
#### outputs ####
#################

plot_dir = paste0(base,"plots/")

stat_dir = paste0(base,"stats/")

###################
#### DATA LOAD ####
###################

dat =
  read_tsv(in_file)


#############################
#### REMOVE CONTAMINANTS ####
#############################


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

## only keep human peptides
dat %<>%
  filter(grepl('_HUMAN',`Entry Name`)) %>%
  dplyr::rename(ID=`Protein ID`)


## Explore number of contaminants in the data that were annotated by either Hao & MQ
dat %>%
  inner_join(con_dat,by="ID") %>%
  group_by(db) %>%
  dplyr::count()


## remove the contaminants
dat %<>%
  left_join(con_dat,by="ID") %>%
  filter(is.na(con))



#######################
#### ENRICHMENT QC ####
#######################


## count and plot how many peptides contain the KeGG mod
## K[114.0429]
qc =
  dat %>%
  mutate(is_remnant = if_else(grepl('K\\[114.0429\\]',`Modified Sequence`),T,F))

qc %>%
  group_by(is_remnant) %>%
  dplyr::count()


qc %<>%
  select(is_remnant, contains('MaxLFQ')) %>%
  pivot_longer(
    cols = contains('MaxLFQ') & !is_remnant,
    values_to = 'value',
    names_to = 'key'
  ) %>%
  filter(value>0) %>%
  group_by(is_remnant,key) %>%
  dplyr::count()


# get the sample names properly matched
samples =
  read_tsv(sample_info_file) %>%
  separate(file,into=c(NA,NA,NA,NA,NA,NA,NA,'sample')) %>%
  select(sample_name,sample) %>%
  separate(sample,into=c('condition','rep'),sep='n') %>%
  mutate(condition = gsub('EV','Control',condition)) %>%
  mutate(condition = gsub('R5KO','RMND5A KO',condition))

qc %<>%
  mutate(sample_name = gsub(' MaxLFQ Intensity','',key)) %>%
  select(-key) %>%
  left_join(samples)


ggplot(qc,aes(condition,n,colour=is_remnant))+
  geom_boxplot(position=position_dodge2(.5),width=.5,size=.1)+
  geom_point(position=position_dodge2(.5),shape=21,fill=NA,size=.8)+
  labs(x='Condition',y='Number of peptides',colour='Contains\nK-GG')+
  scale_colour_manual(values=rev(c("#377eb8ff","#e4211cff")))+
  theme_bw()+
  theme(
    strip.background = element_rect(fill=NA,size=.1),
    axis.ticks=element_line(size=.1),
    axis.ticks.length=unit(.3,"mm"),
    panel.background = element_blank(),
    legend.key.size = unit(.2,"cm"),
    plot.title = element_text(hjust = 0.5,vjust=-1.5),
    plot.subtitle = element_text(hjust = 0.5,vjust=-.5),
    panel.border=element_rect(size=.1),
    panel.grid=element_blank(),
    text=element_text(size=6),
    strip.text.x = element_text(margin = margin(0,0,0,0, "cm")))

ggsave(paste0(plot_dir,'digly_counts.pdf'),
       width=1.6,
       height=1.2)
  

##############################
#### FILTER ON UBREM PEPS ####
##############################


dat %<>%
  filter(grepl('K\\[114.0429\\]',`Modified Sequence`))



#################################
#### Get unique Ubiquitin ID ####
#################################

## get a mapping table of UBID and protein etc
## clean up the other mods as not interested in those
## M[15.9949] C[57.0215] n[42.0106]
## mask 'K\\[114.0429\\]' as B (which does not appear in peptide sequence) to allow position of K-Ub to be found more easily
new_dat =
  dat %>%
  mutate(Peptide = gsub('M\\[15.9949\\]','M',`Modified Sequence`)) %>%
  mutate(Peptide = gsub('C\\[57.0215\\]','C',Peptide)) %>%
  mutate(Peptide = gsub('n\\[42.0106\\]','',Peptide)) %>%
  mutate(Peptide = gsub('K\\[114.0429\\]','B',Peptide)) %>%
  select(ID,
         Gene.x,
         Peptide,
         Start,
         End,
         contains('MaxLFQ')) %>%
  # str_locate_all to locate positions of the modified K (aka B) residues
  # collapse to single comma separated character vector
  mutate(Positions = sapply(str_locate_all(Peptide, "B"), function(x) paste(x[,1], collapse = ","))) %>%
  rowwise() %>%
  mutate(ubid = {
    pos = str_split(Positions,pattern=',')[[1]] %>% as.numeric()
    pos = pos + Start - 1
    # format a ubid name
    paste0(Gene.x,'_',paste(paste0('K.ub.',pos),collapse='_'))
  }) %>%
  ungroup() %>%
  # reverse teh Kub masking
  mutate(Peptide = gsub('B','K.ub.',Peptide))

## now get just ubid and metadata
mapping =
  new_dat %>%
  select(Peptide,ID,Gene.x,ubid) %>%
  rename(Gene = Gene.x) %>%
  distinct()


## now get just data and ubid
## get into long format and group by ubid to collapse any duplicated peptides
long =
  new_dat %>%
  select(ubid,contains('MaxLFQ')) %>%
  pivot_longer(
    cols = contains('MaxLFQ') & !ubid,
    names_to = "key",
    values_to = "value"
  )

## collapse any duplicated entries using median while removing zeros
long %<>%
  filter(value>0) %>%
  group_by(ubid,key) %>%
  mutate(new_val = median(value)) %>%
  select(ubid,key,new_val) %>%
  distinct()

# get the sample names properly matched
samples =
  read_tsv(sample_info_file) %>%
  separate(file,into=c(NA,NA,NA,NA,NA,NA,NA,'sample')) %>%
  select(sample_name,sample) %>%
  separate(sample,into=c('condition','rep'),sep='n')

# join to the table
long %<>%
  mutate(sample_name = gsub(' MaxLFQ Intensity','',key)) %>%
  select(-key) %>%
  left_join(samples)


# now get into wide format again
wide =
  long %>%
  ungroup() %>%
  mutate(sampleID = paste0(condition,'_',rep)) %>%
  select(-sample_name,-condition,-rep,-key) %>%
  pivot_wider(
  names_from = sampleID,
  values_from = new_val) %>%
  select(ubid,starts_with('EV'),starts_with('R5'))



## data
wide

## mappings
mapping


#########################
#### PREPARE FOR DEP ####
#########################


# Generate a SummarizedExperiment object using an experimental design
LFQ_columns <- 2:ncol(wide)# get LFQ column numbers (for area / intensity output (MaxLFQ))

# get experimental design using wide table column names
exp_des = 
  tibble(label = colnames(wide)[LFQ_columns]) %>%
  separate(label,into=c("condition", "replicate"),sep="_",remove=F) %>%
  dplyr::select(label,condition,replicate)

## make sure all LFQ cols are numeric
wide %<>% mutate(across(all_of(LFQ_columns),as.numeric))

## run make_unique to give necessary columns and check each protein is on its own row
wide %<>% make_unique(names = "ubid",ids="ubid")

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

meanSdPlot(data_filt)

########################
#### DEP IMPUTATION ####
########################

# set seed to have imputation be reproducible
set.seed(123)

## performing mixed imputation as per DEP vignette and Goeminne et al 2020
## MNAR defined as missing in 3/4 replicates or 4/4
proteins_MNAR <- get_df_long(data_filt) %>%
  group_by(name, condition) %>%
  summarise(NA_count = sum(is.na(intensity), na.rm = FALSE)) %>%
  filter(NA_count>2) %>% 
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
full =
  get_df_wide(mixed_imputation) %>%
  select(ID,imputed,num_NAs,starts_with('EV'),starts_with('R5')) %>%
  rename(ubid = ID) %>%
  left_join(mapping)

write.table(full,
            file = paste0(dat.dir,'digly_imputed'),
            row.names = F,
            col.names = T,
            quote=F,
            sep="\t")



###########################################
#### ESTIMATE FOLD CHANGES USING LIMMA ####
###########################################

## make the design matrix for limma, derived from sample info in col names
design =
  tibble(coldat = colnames(mixed_imputation)) %>%
  separate(coldat, into=c("condition","rep"),remove=F) %>%
  mutate(group=paste0(condition)) %>%
  mutate(val=1) %>%
  spread(key=group,value=val) %>%
  replace(is.na(.),0) %>%
  dplyr::select(-coldat,-condition,-rep)


# contrasts to use for limma (based on column names)
cont = 'EV-R5KO'

## limma linear regression on mixed imputed data
## get imputed data as numeric matrix
mat = assays(mixed_imputation)[[1]]
fit = lmFit(mat, design)

# set up the contrasts
cont.matrix = makeContrasts(contrasts = cont, levels=design)

fit2 = contrasts.fit(fit, cont.matrix)

# correct using empirical bayes fit
fit2 = eBayes(fit2)

# extract all proteins and log2FC
sig_prots_all =
  topTable(fit2, adjust="none",n=Inf) %>%
  as.data.frame() %>%
  rownames_to_column("ID") %>%
  as_tibble() %>%
  mutate(contrast=cont)


## add mapping info before exports
sig_prots_all %<>%
  left_join(mapping,by=c('ID'='ubid')) %>%
  dplyr::rename(ubid=ID,ID=ID.y)

# export the limma results for this experiment group
## keep all the proteins, useful for checking p values of specific proteins later
write.table(sig_prots_all,
            file = paste0(stat_dir,"digly_limma_all.tsv"),
            row.names = F,
            col.names = T,
            quote=F,
            sep="\t")



######################################################
#### PLOT PROTEIN FOLD CHANGE VS K-UB FOLD CHANGE ####
######################################################


# load the proteome fold changes
prot_fc = 
  read_tsv(paste0(stat_dir,prot_files)) %>%
  filter(grepl('RMND5A',contrast)) %>%
  select(Gene,ID,logFC)


# decide if to just keep p<0.05
ub_fc =
  sig_prots_all %>%
  #filter(P.Value<0.05) %>%
  select(ubid,ID,Gene,logFC)

# combine
comb =
  full_join(prot_fc,ub_fc,by=c("Gene", "ID")) %>%
  dplyr::rename(ub_fc = logFC.y,
         prot_fc = logFC.x) %>%
  drop_na()


# colour MKLN1 in red
comb %<>%
  mutate(type = if_else(Gene=='MKLN1',Gene,'Other'))


## swap the direction of the logFC to make more intuitive sense
## have done this for all volcano plots
## points will be to the right if they were
## higher in the knockout
comb %<>%
  mutate(prot_fc = -prot_fc,
         ub_fc = -ub_fc)


## set up sites to label
comb %<>%
  mutate(label = if_else(prot_fc>1 & ub_fc<(-1), ubid, 'none'))


ggplot(comb,aes(prot_fc,ub_fc,colour=type,label=label))+
  geom_point_rast(data=subset(comb, type=="Other"),alpha=.5,raster.dpi=300,colour="gray60")+ # rasterise as many points are slow
  geom_point_rast(data=subset(comb, type!="Other"),alpha=1,raster.dpi=300,colour="red",size=2)+ # rasterise as many points are slow
  geom_vline(xintercept = 1,linetype='dashed',size=.1)+
  geom_vline(xintercept = -1,linetype='dashed',size=.1)+
  geom_hline(yintercept = 1,linetype='dashed',size=.1)+
  geom_hline(yintercept = -1,linetype='dashed',size=.1)+
  labs(x='Total protein log2FC RMND5A KO:Control',
       y='K-GG log2FC RMND5A KO:Control')+
  theme_bw()+
  theme(panel.grid=element_blank())+
  geom_text_repel(data          = subset(comb, label!="none"),
                  colour="black",
                  size          = 3,
                  #box.padding   = .1,
                  #point.padding = .1,
                  max.overlaps = 400,
                  max.iter = 200,
                  force         = 10,
                  segment.size  = 0.1,
                  segment.color = "black")

ggsave(paste0(plot_dir,"Volcano_digly.pdf"),
       width=5,
       height=5)

#############
#### PCA ####
#############




## get imputed MaxLFQ data
dat =
  full %>%
  select(ubid,contains('EV'),contains('R5KO'))

## filter the wide table to just significant proteins
dat %<>%
  filter(ubid %in% sig_prots_all$ubid)

## get into numeric matrix
mat =
  dat %>%
  dplyr::select(-ubid) %>%
  as.matrix()

## do the PCA on combined data
pc = prcomp(t(mat), scale. = F, center=T)

sig_prot_num_temp = length(unique(sig_prots_all$ubid))

## set up values and sample labels to plot using raw ggplot
pc_vals =
  pc$x %>%
  data.frame() %>%
  rownames_to_column("sample") %>%
  separate(sample, into=c("condition","Rep")) %>%
  mutate(condition = gsub('EV','Control',condition)) %>%
  mutate(condition = gsub('R5KO','RMND5A KO',condition))


## get the variance explained for each component
var_explained <- pc$sdev^2/sum(pc$sdev^2)

# plot PC1 and PC2
pc_vals %>%
  ggplot(aes(x=PC1,y=PC2))+
  geom_point(size=1,aes(colour=condition),shape=21,fill=NA)+
  scale_colour_manual(values=c("#377eb8ff","#e4211cff"))+
  ggtitle("PCA plot of K-GG proteomics",paste0(sig_prot_num_temp," K-GG peptides"))+
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


ggsave(paste0(plot_dir,"PCA_digly.pdf"),
       width=2,
       height=1.5)


