suppressPackageStartupMessages({
    library(optparse)
    library(tidyverse)
    library(data.table)
    library(reshape2)
    library(glue)
})
option_list = list(
  make_option("--in_dir", action = "store", default = NA, type = "character",
              help = "Full path to working directory containing majiq, rmats, leafcutter results"),
  make_option("--in_comp", action = "store", default = NA, type = "character",
              help = "The comparison name along with group header"),
  make_option("--out_dir", action = "store", default = NA, type = "character",
              help = "Full path to output directory."),
  make_option("--groups_file", action = "store", default = NA, type = "character",
              help = "Groups file"),  
  make_option("--out_prefix", action = "store", default = NA, type = "character",
              help = "The name of the output table.  File output will look like {out_prefix}_master_table.csv"),
  make_option("--majiq_cutoffs", action="store", default=NA, type = "character",
              help = "The cutoffs used for majiq_v1 and/or majiq_v2 results.  Separate by space, comma, or semicolon."),
  make_option("--include_tools", action="store", default=NA, type = "character",
              help = "Which tools to include in the master table.  Separate by space, comma, or semicolon.
              Values are: rmats rmats_turbo leafcutter majiq_v1 majiq_v2"),
  make_option("--gtf_path", action="store", default=NA, type = "character",
              help = "full path and file name of gtf file used for analysis.  Used in majiq_v2 hgnc gene symbol mapping.")
)

opt = parse_args(OptionParser(option_list=option_list))
# Rscript src/make_master_table_v3.0.1.r --in_dir /mnt/depts/dept04/compbio/users/zouyang/tools/SplicingPipe/test --in_comp Description1_g_Het_amiR_NT_noASO_vs_Het_NOamiR --out_dir /mnt/depts/dept04/compbio/users/zouyang/tools/SplicingPipe/test/master_table --groups_file /mnt/depts/dept04/compbio/users/zouyang/tools/SplicingPipe/test/leafcutter_output/tmp/Description1_g_Het_amiR_NT_noASO_vs_Het_NOamiR/grp_files.txt --out_prefix Description1_g_Het_amiR_NT_noASO_vs_Het_NOamiR --majiq_cutoffs 0.1,0.2,0.3,0.4 --include_tools rmats,leafcutter,majiq_v2 --gtf_path /edgehpc/dept/compbio/users/mzavodsz/genomes/chimeric/Mouse.GRCm38.vM25.hUNC13A_KI.GRCh38/Mouse.GRCm38.vM25.hUNC13A_KI.GRCh38.v3.gtf
#Debug
if(FALSE){
  opt = list(in_dir = "/edgehpc/dept/compbio/projects/splice_pipeline_BSSI/results/dev/TST11435/DH_ALS_DH_ctrl", 
             out_dir = "/edgehpc/dept/compbio/projects/splice_pipeline_BSSI/results/dev/TST11435/DH_ALS_DH_ctrl", 
             out_prefix = "DH_ALS_DH_ctrl",
             majiq_cutoffs = "0.1,0.2,0.3,0.4",
             include_tools = "rmats,rmats_turbo,leafcutter,majiq_v2",#,majiq_v1
             gtf_path = "/edgehpc/dept/compbio/projects/splice_pipeline_BSSI/annotation/HUman.GRCh38.v28/gencode.v28.annotation.gtf")
}

#Convert tools to vector
include_tools = unlist(strsplit(opt$include_tools, split = "[ ,;]+"))
print(opt$include_tools)
print(include_tools)

#Convert majiq cutoffs to a vector
if(!is.na(opt$majiq_cutoffs)){
majiq_cutoffs = unlist(strsplit(opt$majiq_cutoffs, split = "[ ,;]+"))
}

#Check that the folders are in the in_dir
#need.dirs <- c("majiq_v1","leafcutter","rmats")
if(is.na(opt$in_comp)){
    test.dirs <- file.path(opt$in_dir, include_tools)
    #Check that the files are in the in_dir
    need.files <- list(rmats=c(
        if("rmats" %in% include_tools){
            file.path(opt$in_dir, "rmats", 
                      c('A3SS.MATS.JCEC.txt','A5SS.MATS.JCEC.txt','MXE.MATS.JCEC.txt','RI.MATS.JCEC.txt','SE.MATS.JCEC.txt'))
        }),
        rmats_turbo=c(if("rmats_turbo" %in% include_tools){
            file.path(opt$in_dir, "rmats_turbo", 
                      c('A3SS.MATS.JCEC.txt','A5SS.MATS.JCEC.txt','MXE.MATS.JCEC.txt','RI.MATS.JCEC.txt','SE.MATS.JCEC.txt'))
        }),
        leafcutter=c(if("leafcutter" %in% include_tools){
            file.path(opt$in_dir, "leafcutter", 
                      c("leafcutter_ds_res_cluster_significance.txt","leafcutter_ds_res_effect_sizes.txt","leafviz.Rdata"))
        }),
        majiq_v1=c(if("majiq_v1" %in% include_tools){
            file.path(opt$in_dir, "majiq_v1", paste0("cutoff_",majiq_cutoffs,"_reference_alternative.deltapsi.tsv") )
        }),
        majiq_v2=c(if("majiq_v2" %in% include_tools){
            file.path(opt$in_dir, "majiq_v2", paste0("cutoff_",majiq_cutoffs,"_reference-alternative.deltapsi.voila.tsv") ) #v2 unexpectedly adds 'voila' to filename
        })
    )
}else{
    test.dirs <- file.path(opt$in_dir, paste0(gsub("_v2","",include_tools),"_output"),"tmp",opt$in_comp)
    #Check that the files are in the in_dir
    opt$out_prefix <- opt$in_comp
    strCom <- strsplit(opt$in_comp,"_g_")[[1]][2]
    need.files <- list(
        rmats=c(if("rmats" %in% include_tools){
            file.path(opt$in_dir, "rmats_output","tmp",opt$in_comp, 
                      c('A3SS.MATS.JCEC.txt','A5SS.MATS.JCEC.txt','MXE.MATS.JCEC.txt','RI.MATS.JCEC.txt','SE.MATS.JCEC.txt'))
        }),
        leafcutter=c(if("leafcutter" %in% include_tools){
            file.path(opt$in_dir, "leafcutter_output","tmp",opt$in_comp, 
                      paste0(strCom,"_",c("ds_res_cluster_significance.txt","ds_res_effect_sizes.txt","leafviz.Rdata")))
        }),
        majiq_v2=c(if("majiq_v2" %in% include_tools){
            file.path(opt$in_dir, "majiq_output","tmp",opt$in_comp,
                      paste0(strCom,"_","cutoff_",majiq_cutoffs,"_reference-alternative.deltapsi.voila.tsv") )
        })
    )
}
#print(test.dirs)
if(all(dir.exists(test.dirs))){
  message(paste0("all analysis directories are found in ", opt$in_dir))
} else{ stop(paste0("analysis directories are missing from ", opt$in_dir)) }

if(all(file.exists(unlist(need.files,use.names=F)))){
  message(paste0("all analysis files are found in ", opt$in_dir))
} else{ stop(paste0("analysis files are missing:", paste(unlist(need.files,use.names=F)[!file.exists(unlist(need.files,use.names=F))],collapse="\n"))) }

################## Prep rMATS portion of master table ###############

#rMATS helper functions

#make M
mcol <- function(df){
  #ID column is duplicated, force unique
  #names(df) <- make.unique(names(df))
  
  #Map function
  map.m <- function(x) {sum(as.numeric(unlist(strsplit(x, split = ",")))) / length(as.numeric(unlist(strsplit(x, split = ","))))}
  
  #map.m map for the collapsed columns
  df.i <- df %>% 
    select(ID, IJC_SAMPLE_1, IJC_SAMPLE_2, SJC_SAMPLE_1, SJC_SAMPLE_2) %>%
    mutate(i1 = map(IJC_SAMPLE_1, function(x) map.m(x)),
           i2 = map(IJC_SAMPLE_2, function(x) map.m(x)),
           s1 = map(SJC_SAMPLE_1, function(x) map.m(x)),
           s2 = map(SJC_SAMPLE_2, function(x) map.m(x))
    ) %>%
    select(ID, i1, i2, s1, s2) %>%
    unnest(c(i1,i2,s1,s2)) %>%
    mutate(M = pmax(i1,i2,s1,s2)) %>% 
    select(ID, M)
}

#make psi_1 / psi_2 averages
psicol <- function(df){
  #Map function
  map.m <- function(x) {mean(as.numeric(unlist(strsplit(x, split = ","))),na.rm = TRUE)}
  
  df.i <- df %>% 
    select(ID, PSI_1_persub, PSI_2_persub) %>%
    mutate(psi1 = map(PSI_1_persub, function(x) map.m(x)),
           psi2 = map(PSI_2_persub, function(x) map.m(x))
           ) %>%
    select(ID, psi1, psi2) %>%
    unnest(c(psi1, psi2)) %>%
    select(ID, PSI_1 = psi1, PSI_2 = psi2)
}

if("rmats" %in% include_tools){

#Names of the relevant rmats files within the working directory
rmats_files = need.files[['rmats']]
#file.path(opt$in_dir, "rmats", 
#                        c('A3SS.MATS.JCEC.txt','A5SS.MATS.JCEC.txt','MXE.MATS.JCEC.txt','RI.MATS.JCEC.txt','SE.MATS.JCEC.txt'))
rmats_names = c('A3SS', 'A5SS', 'MXE', 'RI', 'SE')

#Read in the files
rmats_list <- lapply(rmats_files, fread)
rmats_format <- list()
#Format the files
for(i in 1:length(rmats_list)){
  
  file0 <- as.data.frame(rmats_list[[i]])
  
  #ID column is duplicated, force unique
  names(file0) <- make.unique(names(file0))
  
  #Basic reformatting
  file1 <- file0 %>% 
    mutate(Event = rmats_names[i],
           Algorithm = "rmats") %>%
    rename("dPSI" = "IncLevelDifference",
           "PSI_1_persub" = "IncLevel1",
           "PSI_2_persub" = "IncLevel2") 
  #create m
  file.m <- mcol(file0)
  #Add M
  file2 <- file1 %>% 
    left_join(file.m)
  #create PSI_1 and PSI_2 (averages)
  file.psi <- psicol(file2)
  #Add PSI_1 and PSI_2
  file3 <- file2 %>% 
    left_join(file.psi)
  #Event info (collapse all columns not in the keep.col)
  keep.col <- c('chr','strand','geneSymbol', 'GeneID', 'Event', 'Algorithm', 'FDR',
                'PValue', 'dPSI', 'M', 'PSI_1', 'PSI_2')
  file.event <- file3 %>% select(-one_of(keep.col)) %>% select(-ID.1) %>% select(ID, everything())
  file.event[] <- Map(paste, names(file.event), file.event, sep=":")
  file.event <- file.event %>% 
    mutate(ID = as.integer( gsub("ID:","",ID)) ) %>%
    unite(col = Event_info, colnames(.)[-1],sep="|") %>%
    select(ID, Event_info)
  #Join back to file2
  rmats_format[[i]] <- file3 %>% left_join(file.event) %>% select(keep.col, Event_info)
}
#Make the rMATS portion of the table
rmats_format_df <- rbindlist(rmats_format)

}

################## Prep rMATS turbo portion of master table ###############

if("rmats_turbo" %in% include_tools){
  
  #Names of the relevant rmats files within the working directory
  rmats_turbo_files = need.files[['rmats_turbo']]
  #file.path(opt$in_dir, "rmats_turbo", 
  #                        c('A3SS.MATS.JCEC.txt','A5SS.MATS.JCEC.txt','MXE.MATS.JCEC.txt','RI.MATS.JCEC.txt','SE.MATS.JCEC.txt'))
  rmats_turbo_names = c('A3SS', 'A5SS', 'MXE', 'RI', 'SE')
  
  #Read in the files
  rmats_turbo_list <- lapply(rmats_turbo_files, fread)
  rmats_turbo_format <- list()
  #Format the files
  for(i in 1:length(rmats_turbo_list)){
    
    file0 <- as.data.frame(rmats_turbo_list[[i]])
    
    #ID column is duplicated, force unique
    names(file0) <- make.unique(names(file0))
    
    #Basic reformatting
    file1 <- file0 %>% 
      mutate(Event = rmats_turbo_names[i],
             Algorithm = "rmats_turbo") %>%
      rename("dPSI" = "IncLevelDifference",
             "PSI_1_persub" = "IncLevel1",
             "PSI_2_persub" = "IncLevel2") 
    #create m
    file.m <- mcol(file0)
    #Add M
    file2 <- file1 %>% 
      left_join(file.m)
    #create PSI_1 and PSI_2 (averages)
    file.psi <- psicol(file2)
    #Add PSI_1 and PSI_2
    file3 <- file2 %>% 
      left_join(file.psi)
    #Event info (collapse all columns not in the keep.col)
    keep.col <- c('chr','strand','geneSymbol', 'GeneID', 'Event', 'Algorithm', 'FDR',
                  'PValue', 'dPSI', 'M', 'PSI_1', 'PSI_2')
    file.event <- file3 %>% select(-one_of(keep.col)) %>% select(-ID.1) %>% select(ID, everything())
    file.event[] <- Map(paste, names(file.event), file.event, sep=":")
    file.event <- file.event %>% 
      mutate(ID = as.integer( gsub("ID:","",ID)) ) %>%
      unite(col = Event_info, colnames(.)[-1],sep="|") %>%
      select(ID, Event_info)
    #Join back to file2
    file4<-file3 %>% left_join(file.event) %>% select(keep.col, Event_info)
    #as.dataframe was reading these as logical characters in the case of sparse data, which was breaking rbind's ability to join an empty df
    if(nrow(file4)==0){
      file4=data.frame(chr=character(0), strand=character(0), geneSymbol=character(0), GeneID=character(0), Event=character(0), Algorithm=character(0), FDR=numeric(0), PValue=numeric(0), dPSI=numeric(0), M=numeric(0), PSI_1=numeric(0), PSI_2=numeric(0), Event_info=numeric(0))
    }
    rmats_turbo_format[[i]] <- file4
  }
  #Make the rMATS portion of the table
  rmats_turbo_format_df <- rbindlist(rmats_turbo_format)
  
  
}

############################# Prep leafcutter portion of master table ##################################

if("leafcutter" %in% include_tools){
#Names of the relevant leafcutter files within the working directory
leafcutter_files = need.files[['leafcutter']]
    #file.path(opt$in_dir, "leafcutter", 
    #                         c("leafcutter_ds_res_cluster_significance.txt","leafcutter_ds_res_effect_sizes.txt","leafviz.Rdata"))

# read in files:
#pvalues
table_cs <- fread(leafcutter_files[1])
#drop NA from table_cs (can keep track of it with the status)
table_cs.badstatus <- table_cs %>% filter(status != "Success")
table_cs <- table_cs %>% drop_na
#psis
table_es <- fread(leafcutter_files[2])
#anno table
load(leafcutter_files[3])
table_anno <- clusters 

#add a cluster column that matches table_cs to table_es (format: chr:clu_number)
table_es$cluster <- paste0(str_split_fixed(table_es$intron, ":", 5)[,1],":",str_split_fixed(table_es$intron, ":", 5)[,4])

#make the cluster column in table_anno match as well
table_anno$cluster <- paste0(str_split_fixed(table_anno$coord,":",2)[,1],":",table_anno$clusterID)

#drop unneeded columns from table_anno, format gene
table_anno <- table_anno %>% 
  mutate(gene = gsub("(<i>|</i>)","",gene) ) %>%
  select(-clusterID, -FDR)

#merge table_es and table_cs
table_es_cs_anno <- table_es %>% full_join(table_cs, by="cluster") %>% full_join(table_anno, by="cluster")

#Read groups from file
if(file.exists(opt$groups_file)){
    groups_file <- fread(opt$groups_file,header = F)
    groups <- unique(groups_file$V2)
}else{
    groups <- unlist(strsplit(opt$groups_file,","))
}


#format the columns
leafcutter_format_df.all <- table_es_cs_anno %>%
  select(geneSymbol = gene, 
         Event = cluster,
         FDR = p.adjust,
         PValue = p,
         dPSI = deltapsi,
         M = status,
         PSI_1 = paste(groups[1]),
         PSI_2 = paste(groups[2]),
         Annotation_leafcutter = annotation,
         everything()) %>%
  mutate(Algorithm = "leafcutter")

leafcutter_format_df.all$chr <- str_split_fixed(leafcutter_format_df.all$coord, ":", 2)[,1]
leafcutter_format_df.all$strand <- str_split_fixed(leafcutter_format_df.all$intron, "_", 3)[,3]

leafcutter_format_df.all <- leafcutter_format_df.all %>% filter(FDR < 1)

#make the event info (use intron as an index for joining)
keep.col <- c("chr","strand","geneSymbol", "Event", "FDR", "PValue", "dPSI", "M", "PSI_1", "PSI_2", "Annotation_leafcutter", "Algorithm")
leafcutter.format.event <- leafcutter_format_df.all %>% select(-one_of(keep.col)) %>% select(intron, everything())
leafcutter.format.event[] <- Map(paste, names(leafcutter.format.event), leafcutter.format.event, sep=":")
leafcutter.format.event <- leafcutter.format.event %>% 
  mutate(intron = gsub("intron:","",intron),
         intron.info = intron) %>%
  unite(col = Event_info, colnames(.)[-1],sep="|") %>%
  select(intron, Event_info)

#make the final leafcutter file
leafcutter_format_df.final <- leafcutter_format_df.all %>% select(intron, one_of(keep.col)) %>% full_join(leafcutter.format.event, by="intron") %>%
  select(-intron)
}

############################# Prep majiq (v1) portion of master table #######################

if("majiq_v1" %in% include_tools){

#Get the files
#Names of the relevant majiq files within the working directory
majiq_files = need.files[['majiq_v1']]
#file.path(opt$in_dir, "majiq_v1", paste0("cutoff_",majiq_cutoffs,"_reference_alternative.deltapsi.tsv"))
#Names of each cutoff (flexible)
majiq_cutoff_names <- str_split_fixed(basename(majiq_files),"_",4)[,2]

#Read in the files
majiq_list <- lapply(majiq_files, fread)
majiq_format <- list()

for(i in 1:length(majiq_list)){
  
  majiq0 <- majiq_list[[i]]
  
  #Start to format, several columns need to be un-collapsed
  majiq.format.all <- majiq0 %>% 
    select(geneSymbol = `#Gene Name`,
           GeneID = `Gene ID`,
           Event = `LSV ID`, #Use as index
           everything()
    ) %>%
    mutate(Algorithm = "majiq_v1",
           Event_type = case_when(A5SS & !A3SS & !ES ~ "A5SS",
                                  A3SS & !A5SS & !ES ~ "A3SS",
                                  ES & !A3SS & !A5SS ~ "ES",
                                  A5SS & A3SS & !ES ~ "A5SS_A3SS",
                                  A5SS & ES & !A3SS ~ "A5SS_ES",
                                  A3SS & ES & !A5SS ~ "A3SS_ES",
                                  A3SS & A5SS & ES ~ "A5SS_A3SS_ES",
                                  !A5SS & !A3SS & !ES ~ "MXE/RI") ) #look at any of these cases, anything with >4 exons too
  
  #These columns need to be un-collapsed
  majiq.expand <- majiq.format.all %>% select(Event, 
                                              dPSI = `E(dPSI) per LSV junction`,
                                              PSI_1 = `alternative E(PSI)`,
                                              PSI_2 = `reference E(PSI)`,
                                              Junc_coords = `Junctions coords`, 
                                              PdPSI = matches("^P\\(.*>") #regex flexible for different cutoffs
                                              ) 
  
  expand.cols <- c("dPSI", "PSI_1", "PSI_2", 
                   "Junc_coords", 
                   "PdPSI")
  expand.list <- list()
  for(j in 1:length(expand.cols)){
    expand.list[[j]] <- majiq.expand %>% select(Event, !!expand.cols[[j]]) %>% separate_rows(!!expand.cols[[j]], sep=";")
  }
  
  #Make the uncollapsed dataframe and remove Event column duplicates
  majiq.expand <- do.call(cbind, expand.list)
  majiq.expand <- majiq.expand[!duplicated(as.list(majiq.expand))]
  
  #Join expand to the table
  majiq.format.df1 <- majiq.format.all %>% dplyr::select(chr, strand, geneSymbol, GeneID, Event, Event_type, Algorithm) %>%
    full_join(majiq.expand)
  
  #Some columns to keep
  keep.col = c("chr","strand",'geneSymbol', "GeneID", 
               'Event_type', 'Algorithm', 
               'PdPSI', 'dPSI', 'PSI_1', 'PSI_2', 'Junc_coords', 'Event_info')
  
  #Make the event info (Event column as index)
  majiq.format.event <- majiq.format.all %>% dplyr::select(Event, `LSV Type`, `Exons coords`, `IR coords`)
  
  majiq.format.event[] <- Map(paste, names(majiq.format.event), majiq.format.event, sep=":")
  majiq.format.event <- majiq.format.event %>% 
    mutate(Event = gsub("Event:","",Event)) %>%
    unite(col = Event_info, colnames(.)[-1],sep="|") %>%
    select(Event, Event_info)
  
  majiq_format[[i]] <- majiq.format.df1 %>% full_join(majiq.format.event) 
}
names(majiq_format) <- majiq_cutoff_names
majiq_format_df.final <- rbindlist(majiq_format, idcol = "Cutoff_dPSI")

}

############################# Prep majiq (v2) portion of master table #######################

if("majiq_v2" %in% include_tools){
#Get the files
#Names of the relevant majiq files within the working directory
majiq_v2_files = need.files[['majiq_v2']]
#file.path(opt$in_dir, "majiq_v2", paste0("cutoff_",majiq_cutoffs,"_reference-alternative.deltapsi.voila.tsv")) #Unexpected voila added to filename in majiq v2
#Names of each cutoff (flexible)
majiq_v2_cutoff_names <- str_split_fixed(basename(majiq_v2_files),"_",4)[,2]

#Read in the files
majiq_v2_list <- lapply(majiq_v2_files, fread)
majiq_v2_format <- list()

#Read in gtf file for hgnc gene symbol mapping
gtf.tx <- fread(file.path(opt$gtf_path)) #load the gtf file w/tx info
colnames(gtf.tx) <- c("CHR", "source", "type", "start", "end", "score", "strand", "frame", "info")
gtf.tx$gene_name <- str_trim(gsub("(gene_name |\")", "",str_split_fixed(gtf.tx$info, pattern = ";", n = 5)[,4]),side="both")
gtf.tx$gene_id <- str_trim(gsub("(gene_id |\")","",str_split_fixed(gtf.tx$info, pattern = ";", n = 2)[,1]),side="both")
#Make a gene map
gene_map <- gtf.tx %>% dplyr::select(gene_name, gene_id, chr=CHR, strand) %>% distinct()

rm(gtf.tx)

for(i in 1:length(majiq_v2_list)){
  
  majiq_v2_0 <- majiq_v2_list[[i]]
  
  #Start to format, several columns need to be un-collapsed
  majiq_v2_format_all <- majiq_v2_0 %>% 
    left_join(gene_map) %>%
    select(geneSymbol = gene_name,
           GeneID = gene_id,
           Event = lsv_id, #Use as index
           everything()
    ) %>% mutate(Algorithm = "majiq_v2",
                 Event_type=NA)
  
  #These columns need to be un-collapsed
  majiq_v2_expand <- majiq_v2_format_all %>% select(Event, 
                                              dPSI = mean_dpsi_per_lsv_junction,
                                              PSI_1 = alternative_mean_psi,
                                              PSI_2 = reference_mean_psi,
                                              Junc_coords = junctions_coords, 
                                              PdPSI = probability_changing
  ) 
  
  expand.cols <- c("dPSI", "PSI_1", "PSI_2", 
                   "Junc_coords", 
                   "PdPSI")
  expand.list <- list()
  for(j in 1:length(expand.cols)){
    expand.list[[j]] <- majiq_v2_expand %>% select(Event, !!expand.cols[[j]]) %>% separate_rows(!!expand.cols[[j]], sep=";")
  }
  
  #Make the uncollapsed dataframe and remove Event column duplicates
  majiq.expand <- do.call(cbind, expand.list)
  majiq.expand <- majiq.expand[!duplicated(as.list(majiq.expand))]
  
  #Join expand to the table
  majiq.format.df1 <- majiq_v2_format_all %>% dplyr::select(chr, strand, geneSymbol, GeneID, Event, Event_type, Algorithm) %>%
    full_join(majiq.expand)
  
  #Some columns to keep
  keep.col = c("chr","strand",'geneSymbol', "GeneID", 
               'Event_type', 'Algorithm', 
               'PdPSI', 'dPSI', 'PSI_1', 'PSI_2', 'Junc_coords', 'Event_info')
  
  #Make the event info (Event column as index)
  majiq.format.event <- majiq_v2_format_all %>% dplyr::select(Event, lsv_type, exons_coords, ir_coords) 
  
  majiq.format.event[] <- Map(paste, names(majiq.format.event), majiq.format.event, sep=":")
  majiq.format.event <- majiq.format.event %>% 
    mutate(Event = gsub("Event:","",Event)) %>%
    unite(col = Event_info, colnames(.)[-1],sep="|") %>%
    select(Event, Event_info)
  
  majiq_v2_format[[i]] <- majiq.format.df1 %>% full_join(majiq.format.event) 
}
names(majiq_v2_format) <- majiq_v2_cutoff_names
majiq_format_df_v2.final <- rbindlist(majiq_v2_format, idcol = "Cutoff_dPSI")
}

############ Merge all tables and write out ##################

bind_list <- list()

if("rmats" %in% include_tools){
  bind_list <- append(bind_list, list(rmats_format_df))
}
if("rmats_turbo" %in% include_tools){
  bind_list <- append(bind_list, list(rmats_turbo_format_df))
}
if("leafcutter" %in% include_tools){
  bind_list <- append(bind_list, list(leafcutter_format_df.final))
}
if("majiq_v1" %in% include_tools){
  bind_list <- append(bind_list, list(majiq_format_df.final))
}
if("majiq_v2" %in% include_tools){
  bind_list <- append(bind_list, list(majiq_format_df_v2.final))
}

master_table <- rbindlist(bind_list, fill = TRUE )
fwrite(master_table, file.path(opt$out_dir, glue("{opt$out_prefix}_master_table.csv")), row.names = FALSE)
