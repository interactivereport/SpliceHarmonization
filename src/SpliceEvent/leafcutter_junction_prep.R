suppressPackageStartupMessages({
    library(tidyr)
    library(stringr)
    library(readr)
    library(dplyr)
    library(GenomicRanges)
    library(optparse)
    library(parallel)
})
num_cores <- detectCores() - 3

option_list <- list(
  make_option(
    c("-i", "--input"),
    type = "character",
    default = "leafcutter/*_cluster_significance.txt",
    help = "Path with glob character to Leafcutter result files.
    [default %default]",
    metavar = "character"
  ),
  make_option(
    c("-o", "--output"),
    type = "character",
    default = "leafcutter/leafcutter_junctions_prep.csv",
    help = "Path to output file [default %default]",
    metavar = "character"
  )
)

opt <- parse_args(OptionParser(option_list = option_list))
files <- Sys.glob(opt$input)

file_names <- str_split(files, "/", simplify = T)
file_names <- file_names[, ncol(file_names) - 1]
message("Comparison names: ", paste0(file_names, collapse = "\t"))

cluster_sig_file <- files

cluster_sig <- lapply(
  cluster_sig_file,
  read.table,
  sep = "\t",
  header = TRUE
)

names(cluster_sig) <- file_names
cluster_sig <- bind_rows(cluster_sig, .id = "comparison")

effec_size_files <- gsub(
  x = files,
  pattern = "cluster_significance",
  replacement = "effect_sizes"
)
es <- lapply(
  effec_size_files,
  read.table,
  sep = "\t",
  header = T,
  colClasses = c("character", "double", "double", "double", "double"),
  col.names=c( "intron", "logef", "reference", "alternative", "deltapsi")
)
names(es) <- file_names
es <- bind_rows(es, .id = "comparison")


# parse the intron column for merging
es$intron <- str_replace_all(es$intron, "_", ":")

intron <- str_split(es$intron, pattern = ":", simplify = T)
colnames(intron) <- c(
  "chr", "start", "end", "clu", "clu_number", "strand"
)
intron <- as_tibble(intron)

es <- bind_cols(es, intron)
# cluster will be the pivot for merging
es$cluster <- as.character(
  str_glue_data(es, "{chr}:{clu}_{clu_number}_{strand}")
)
es$cluster2 <- as.character(
  str_glue_data(es, "{chr}:{clu}_{clu_number}_{strand}")
)
message("Merging tables")


res <- inner_join(es, cluster_sig, by = c("comparison", "cluster"))
# res$chr <- gsub("chr", "", res$chr)

res <- select(res, -c("clu", "clu_number"))
# create a unique junction column for each row
message("Number of junctions output by Leafcutter ", nrow(res))
res <- res %>%
  filter(p.adjust < 1.1) %>%
  mutate(method = "Leafcutter") %>%
  arrange(p.adjust) %>%
  distinct(comparison, chr, start, end, strand, .keep_all = TRUE)

res$ID <- str_split(res$cluster, ":", simplify = TRUE)[, 2]
res$cluster <- res$cluster2
res$cluster2 <- NULL
message("Number of junctions after filtering ", nrow(res))


leafViz_files <- gsub(
    x=files, 
    pattern= 'ds_res_cluster_significance.txt', 
    replacement='leafviz.Rdata')

leafViz <- lapply(leafViz_files,function(file){
    load(file) 
    df <- clusters  %>%
    mutate(clusters = paste0(str_split_fixed(coord, ":", 2)[, 1], ":", clusterID)) %>%
    mutate(gene = str_remove_all(gene, "<i>|</i>")) %>%
    select(-clusterID, -FDR)
    return(df)
    })

names(leafViz) <- file_names
leafViz <- bind_rows(leafViz, .id = "comparison")

res <- res %>% 
    left_join(leafViz, by = c("cluster" = "clusters")) %>% 
    select(c('comparison.x', 'intron', 'reference',
                        'alternative', 'chr', 
                        'start', 'end', 
                        'strand', 'p.adjust', 
                        'cluster', 'method', 'gene', 'annotation', 'ID')) %>%
    dplyr::rename('comparison' = 'comparison.x', 'PSI_1' = 'reference', 'PSI_2' = 'alternative', 'geneSymbol' ='gene', 'type' = 'annotation') 

message("Number of junctions after filtering ", nrow(res))
write_csv(res, opt$output)