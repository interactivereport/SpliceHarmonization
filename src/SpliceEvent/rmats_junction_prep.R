suppressPackageStartupMessages({
    library(tidyr)
    library(stringr)
    library(readr)
    library(dplyr)
    library(GenomicRanges)
    library(optparse)
    library(purrr)
})

option_list <- list(
  make_option(
    c("-i", "--input"),
    type = "character",
    # default = "rmats/*.MATS.JC.txt",
    default = "rmats/*.MATS.JCEC.txt", # change it from JC to JCEC
    help = "Path with glob character to Leafcutter result files.
    [default %default]",
    metavar = "character"
  ),
  make_option(
    c("-o", "--output"),
    type = "character",
    default = "rmats/rmats_junctions_prep.csv",
    help = "Path to output file [default %default]",
    metavar = "character"
  ),

  make_option(
    c("-n", "--ref_name"),
    type = "character",
    default = "DMSO",
    help = "Specify the DMSO string [default]",
    metavar = "character"
  ),

  make_option(
    c("-M", "--M_cutoff"),
    type = 'double',
    default = 0,
    help = 'M_cutoff values',
    metavar = 'character'
  )
)

opt <- parse_args(OptionParser(option_list = option_list))
files <- Sys.glob(opt$input)
print(files)
map.m <- function(x) {
  sum(as.numeric(unlist(strsplit(x, split = ",")))) / length(as.numeric(unlist(strsplit(x, split = ","))))
}

process_RMATS_ass <- function(df, start, end , type, FDR = 1.1, m_cutoff=opt$M_cutoff) {
  df <- df %>%
    mutate(SJC1_mean = map_dbl(SJC_SAMPLE_1, function(x) map.m(x)),
         SJC2_mean = map_dbl(SJC_SAMPLE_2, function(x) map.m(x)),  # Corrected this line
         IJC1_mean = map_dbl(IJC_SAMPLE_1, function(x) map.m(x)),
         IJC2_mean = map_dbl(IJC_SAMPLE_2, function(x) map.m(x))) %>% 
    mutate(max_v = pmax(SJC1_mean, SJC2_mean, IJC1_mean, IJC2_mean, na.rm = TRUE)) %>%
    dplyr::filter(max_v >= m_cutoff) %>%
    dplyr::filter(FDR < !!FDR) 

    if (nrow(df) == 0) {
      return(df)
      }

    df<- df %>%
    dplyr::select(
      ID, chr, !!start, !!end, strand, comparison, FDR, IncLevel1, IncLevel2, max_v, geneSymbol
    ) %>%
    dplyr::rename(c(start = !!start, end = !!end, PSI_1 = IncLevel1, PSI_2 = IncLevel2, M = max_v)) %>%
    mutate(start = pmin(start, end), end = pmax(start, end)) 
    # %>%
    # mutate(PSI_1 = mean(as.numeric(unlist(strsplit(PSI_1, split = ","))))) %>%
    # mutate(PSI_2 = mean(as.numeric(unlist(strsplit(PSI_2, split = ",")))))
  df$type <- type
    
    df

}

#' @param df dataframe from rMATs
#' @param start name of the start column
#' @param end name of the end column
#' @param type flag for splice junction type
#' @param FDR is the FDR cutoff
#' @return GenomicRange of the selected SJ
#' @export
#'
process_RMATS <- function(df, start, end, type, FDR = 1.1, m_cutoff=opt$M_cutoff) {
  df <- df %>%
    mutate(SJC1_mean = map_dbl(SJC_SAMPLE_1, function(x) map.m(x)),
         SJC2_mean = map_dbl(SJC_SAMPLE_2, function(x) map.m(x)),  # Corrected this line
         IJC1_mean = map_dbl(IJC_SAMPLE_1, function(x) map.m(x)),
         IJC2_mean = map_dbl(IJC_SAMPLE_2, function(x) map.m(x))) %>% 
    mutate(max_v = pmax(SJC1_mean, SJC2_mean, IJC1_mean, IJC2_mean, na.rm = TRUE)) %>%
    dplyr::filter(max_v >= m_cutoff) %>%
    dplyr::filter(FDR < !!FDR) 

    if (nrow(df) == 0) {
      return(df)
      }

    df <- df %>% 
    dplyr::select(
      ID, chr, !!start, !!end, strand, comparison,
      FDR, IncLevel1, IncLevel2, max_v, geneSymbol) %>%
    dplyr::rename(c(start = !!start, end = !!end, PSI_1 = IncLevel1, PSI_2 = IncLevel2, M = max_v)) 


  df$type <- type
    df

}


get_rmats_coord <- function(.files, .comparison) {
  message("Processing files for ", .comparison)
  x <- lapply(.files, read.delim, "\t", header = TRUE)

  lapply(seq_along(.files), function(i) {
    message(.files[[i]], ": has ", nrow(x[[i]]), " entries")
  })
  .split_files <- strsplit(.files, "/", fixed = T)
  as_type <- str_replace(res[['rmats']], paste0(".*(?:", opt$ref_name, "_|", opt$ref_name, "-bridge_)(.*?)\\.txt"), "\\1")

  # names(x) <- as_type
  # this for loop need to be commented when no rmats in simulation
  names(x) <- as_type
  for (i in names(x)) {
    x[[i]]$comparison <- .comparison
      x[[]]
  } 
  # this for loop need to be commented when no rmats in simulation
    
  es_ssj <- process_RMATS(
      x$SE, "upstreamEE", "downstreamES", "ES_SJ",
      FDR = 1.1
    )
    es_isj1 <- process_RMATS(
      x$SE, "upstreamEE", "exonStart_0base", "EIL_SJ",
      FDR = 1.1
    )
      
    es_isj2 <- process_RMATS(
      x$SE, "exonEnd", "downstreamES", "EIR_SJ",
      FDR = 1.1
    )
 
    
    mxe_sj11 <- process_RMATS(
      x$MXE, "upstreamEE", "X1stExonStart_0base", "MXE1L_SJ",
      FDR = 1.1
    )
      
      
    mxe_sj12 <- process_RMATS(
      x$MXE, "X1stExonEnd", "downstreamES", "MXE1R_SJ",
      FDR = 1.1
    ) 
      
      
    mxe_sj21 <- process_RMATS(
      x$MXE, "upstreamEE", "X2ndExonStart_0base", "MXE2L_SJ",
      FDR = 1.1
    )
    mxe_sj22 <- process_RMATS(
      x$MXE, "X2ndExonEnd", "downstreamES", "MXE2R_SJ",
      FDR = 1.1
    ) 

    ir <- process_RMATS(
      x$RI, "upstreamEE", "downstreamES", "IR",
      FDR = 1.1
    )
    
    irL  <- process_RMATS_ass(
      x$RI, "upstreamEE", "riExonStart_0base", "IRL",
        FDR =1.1
    )
    
      
    irR  <- process_RMATS_ass(
      x$RI, "riExonEnd", "downstreamES", "IRR",
        FDR =1.1
    )  
  
  
  
    a5ss_SSJ_pos <- process_RMATS_ass(
              x$A5SS[x$A5SS['strand'] == "+", ], "longExonEnd", "flankingES", "A5SS_SSJ",
              FDR = 1.1
            )
      
    a5ss_SSJ_neg <- process_RMATS_ass(
              x$A5SS[x$A5SS['strand'] == "-", ], "flankingEE", "longExonStart_0base", "A5SS_SSJ",
              FDR = 1.1
            )
      
    a5ss_ISJ_pos <- process_RMATS_ass(
              x$A5SS[x$A5SS['strand'] == "+", ], "shortEE", "flankingES", "A5SS_ISJ",
              FDR = 1.1
            )
      
    a5ss_ISJ_neg <- process_RMATS_ass(
              x$A5SS[x$A5SS['strand'] == "-", ], "flankingEE", "shortES", "A5SS_ISJ",
              FDR = 1.1
            )

    
  
  a3ss_SSJ_neg <- process_RMATS_ass(
              x$A3SS[x$A3SS['strand'] == "-", ], "longExonEnd", "flankingES", "A3SS_SSJ",
              FDR = 1.1
            )
      
    a3ss_SSJ_pos <- process_RMATS_ass(
              x$A3SS[x$A3SS['strand'] == "+", ], "flankingEE", "longExonStart_0base", "A3SS_SSJ",
              FDR = 1.1
            )
      
    a3ss_ISJ_neg <- process_RMATS_ass(
              x$A3SS[x$A3SS['strand'] == "-", ], "shortEE", "flankingES", "A3SS_ISJ",
              FDR = 1.1
            )
      
    a3ss_ISJ_pos <- process_RMATS_ass(
              x$A3SS[x$A3SS['strand'] == "+", ], "flankingEE", "shortES", "A3SS_ISJ",
              FDR = 1.1
            )


  res <- bind_rows(
    lst(es_ssj, es_isj1, es_isj2, mxe_sj11, mxe_sj12, mxe_sj21, mxe_sj22, ir, irL, irR, 
        a5ss_SSJ_pos, a5ss_SSJ_neg, 
        a5ss_ISJ_pos, a5ss_ISJ_neg, a3ss_SSJ_pos, a3ss_SSJ_neg, 
        a3ss_ISJ_pos, a3ss_ISJ_neg
       )
      )
  res$method <- "rmats"
  
      res
}


split_path <- strsplit(files, "/", fixed = T)
comparison <- 'rmats'

res <- split(files, comparison)

message("Loading and procesing rMATs files")
res <- lapply(setNames(names(res), names(res)), function(x) {
  get_rmats_coord(res[[x]], x)
})
res <- bind_rows(res)

res$ID <- paste0(res$type, res$ID) 
res <- res[res$end - res$start > 1, ]
res<- res %>%
    arrange(FDR) #%>%
    #distinct(comparison, chr, start, end, strand,  .keep_all = TRUE)

message("Number of junctions after filtering ", nrow(res))
write_csv(res, opt$output)