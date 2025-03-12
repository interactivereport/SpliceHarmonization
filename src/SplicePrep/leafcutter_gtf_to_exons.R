#/edgehpc/dept/compbio/users/emarshal/splice_pipeline_testing/RMATS_test_data/leafcutter/leafcutter/scripts/gtf_to_exons.R
suppressPackageStartupMessages({
library(data.table)
library(dplyr)
})
main <- function(){
    args = commandArgs(trailingOnly=TRUE)
    if (length(args)<2) {
        stop("Usage is: Rscript gtf_to_exons_for_leafcutter.R input.gtf.gz output.txt.gz")
    }
    cat("Reading in ",args[1],"\n")
    gtf=fread(args[1], data.table = F, col.names=c("chr","source","feature","start","end","a","strand","b","dat"))
    
    cat("Processing...\n")
    gtf = gtf %>% filter( feature=="exon" )
    
    gn_where=regexpr("gene_name \"[^ ]+\"" , gtf$dat) # find gene_names in dat
    gn_where=gn_where + 11 # ignore "gene_name" label
    attr(gn_where,"match.length")=attr(gn_where,"match.length") - 11- 1 # cutoff trailing quote mark
    gtf$gene_name=regmatches(gtf$dat, gn_where )
    gtf = gtf %>% select( chr, start, end, strand, gene_name ) %>% distinct()
    cat("Saving exons to ",args[2],"\n")
    gz=gzfile(args[2],"w")
    write.table(gtf, gz, row.names = F, quote=F, sep="\t")
    close(gz)
}


main()