suppressPackageStartupMessages({
    library(tidyr)
    library(stringr)
    library(readr)
    library(dplyr)
    library(GenomicRanges)
    library(optparse)
})

option_list <- list(
  make_option(
    c("-i", "--input"),
    type = "character",
    default = "majiq/*.tsv",
    help = "Path with glob character to Leafcutter result files.
    [default %default]",
    metavar = "character"
  ),
  make_option(
    c("-o", "--output"),
    type = "character",
    default = "majiq/majiq_junctions_prep.csv",
    help = "Path to output file [default %default]",
    metavar = "character"
  )
)

opt <- parse_args(OptionParser(option_list = option_list))
files <- Sys.glob(opt$input)
file_names <- gsub(
    x = files,
    replacement = "majiq_\\1",
    pattern = ".*/.*cutoff_([0-9]\\.[0-9]).*")

read_majiq_out <- function(x) {
  read.table(
    x,
    sep = "\t",
    header = TRUE,
    col.names = c(
      "gene_name",
      "gene_id",
      "lsv_id",
      "mean_dpsi_per_lsv_junction",
      "probability_changing",
      "probability_non_changing",
      "ref_mean_psi",
      "alt_mean_psi",
      "lsv_type",
      "num_junctions",
      "num_exons",
      "de_novo_junctions",
      "chr",
      "strand",
      "junctions_coords",
      "exons_coords",
      "ir_coords",
      "ucsc_lsv_link"
    ),
    na.strings=c("\"\"","","na"),
    colClasses = c(
      "character",
      "character",
      "character",
      "character",
      "character",
      "character",
      "character",
      "character",
      "character",
      "integer",
      "integer",
      "character",
      "character",
      "character",
      "character",
      "character",
      "character",
      "character"
    )
  )
}

res <- lapply(files, read_majiq_out)
names(res) <- file_names
# res <- bind_rows(res, .id = "comparison") %>%
#   mutate(
#     lsv_type = substr(lsv_type, 3, nchar(lsv_type)),
#     lsv_type = str_replace_all(lsv_type, "\\|", ";")
#   ) %>% 
#   separate_rows(
#     "junctions_coords",
#     "mean_dpsi_per_lsv_junction",
#     "probability_changing",
#     "probability_non_changing",
#     "lsv_type",
#     "ref_mean_psi",
#     "alt_mean_psi",
#     "gene_name",
#     sep = ";",
#     convert = T
#   )

res <- bind_rows(res, .id = "comparison") %>%
  mutate(
    lsv_type = substr(lsv_type, 3, nchar(lsv_type)),
    lsv_type = str_replace_all(lsv_type, "\\|", ";")
  ) %>% 
  mutate(across(c(junctions_coords, mean_dpsi_per_lsv_junction, probability_changing,
                  probability_non_changing, lsv_type, ref_mean_psi, alt_mean_psi, gene_name),
                ~str_replace_all(., c("^na;|;na$" = "", ";na;|na" = "")))) %>%
  separate_rows(
    junctions_coords,
    mean_dpsi_per_lsv_junction,
    probability_changing,
    probability_non_changing,
    lsv_type,
    ref_mean_psi,
    alt_mean_psi,
    gene_name,
    sep = ";",
    convert = TRUE
  )

junction_pattern <- "(\\d+)-(\\d+)"
junctions_coords <- str_match(
  res$junctions_coords, junction_pattern
)[, c(2, 3)]

res["start"] <- junctions_coords[, 1]
res["end"] <- junctions_coords[, 2]

res <- res %>% 
  mutate(
    start = as.numeric(start),
    end = as.numeric(end),
    start = ifelse(lsv_type == "i", start - 1, start),
    end = ifelse(lsv_type == "i", end + 1, end)
  )

message("Number of junctions output by Majiq ", nrow(res))
res <- dplyr::select(res, c(
  chr,
  start,
  end,
  strand,
  comparison,
  probability_changing,
  probability_non_changing,
  lsv_id,
  lsv_type,
  ref_mean_psi,
  alt_mean_psi,
  mean_dpsi_per_lsv_junction,
  gene_name
)) %>%
  filter(probability_non_changing < 1.1) %>%
  mutate(method = "majiq") %>%
  arrange(probability_non_changing) #%>%
 # distinct(comparison, chr, start, end, strand, .keep_all = TRUE)

res <- res %>% 
    dplyr::rename('dpsi' = 'mean_dpsi_per_lsv_junction', 'PSI_1' = 'ref_mean_psi',  'PSI_2' = 'alt_mean_psi','PdPSI' = 'probability_changing' , 'geneSymbol' = 'gene_name')

write_csv(res, opt$output)