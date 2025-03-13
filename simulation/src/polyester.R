suppressPackageStartupMessages({
    library(polyester)
    library(Biostrings)
    library(optparse)
})

option_list <- list(
  make_option(
    "--ref_input", 
    type = "character",
    default = NULL,
    help = "Path with ref_fasta file",
    metavar = "REF_INPUT"
  ),
  make_option(
    "--alt_input", 
    type = "character",
    default = NULL,
    help = "Path with alt_fasta file",
    metavar = "ALT_INPUT"
  ),
  make_option(
    "--output_dir", 
    type = "character",
    default = NULL,
    help = "Path to output dir",
    metavar = "OUTPUT_DIR"
  ),
  make_option(
    "--num_reads", 
    type = "integer",
    default = NULL,
    help = "Total count seq read for the fasta",
    metavar = "NUM_READS"
  ),
  make_option(
    "--cntl_propo", 
    type = "double",
    default = NULL,
    help = "Proportional of control",
    metavar = "CNTL_PROPO"
  ),
  make_option(
    "--alt_propo", 
    type = "double",
    default = NULL,
    help = "Proportional of alternative",
    metavar = "ALT_PROPO"
  )
)

args <- commandArgs(trailingOnly = TRUE)
parser <- OptionParser(option_list = option_list)

tryCatch({
  opt <- parse_args(parser, args)
}, error = function(e) {
  cat("Error parsing options: ", e$message, "\n")
  quit(status = 1)
})

# if (!is.null(opt$ref_input)) {
#   print(opt$ref_input)
# } else {
#   print("No reference input provided.")
# }


# print(opt$cntl_propo)
# print(opt$alt_propo)
set.seed(123)

fasta_file_A <- opt$ref_input
fasta_file_B <- opt$alt_input

# print(fasta_file_A)

if (!is.null(opt$cntl_propo) && !is.null(opt$alt_propo) && opt$cntl_propo > 0 && opt$alt_propo > 0) {
  outdir <- paste0(opt$output_dir, '_num_' , opt$num_reads , '_cntlP_', opt$cntl_propo ,'_altP_' ,opt$alt_propo)
  
} else {
  outdir <- paste0(opt$output_dir, '_num_' , opt$num_reads )
}


cntl_outdir <- paste0(outdir, '/condition_cntl')
alt_outdir <- paste0(outdir , '/condition_alt')

if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE)
}

if (!dir.exists(cntl_outdir)) {
  dir.create(cntl_outdir, recursive = TRUE)
}

if (!dir.exists(alt_outdir)) {
  dir.create(alt_outdir, recursive = TRUE)
}

total_reads_per_file <- opt$num_reads

# control
simulate_experiment(
  fasta = fasta_file_A,
  num_reps = 3,
  readlen =100, 
  reads_per_transcript=total_reads_per_file, 
  fold_changes=3,
  outdir=cntl_outdir,
  paired = TRUE
)

# alternative
# print(opt$cntl_propo)

if (!is.null(opt$cntl_propo) && !is.null(opt$alt_propo) && opt$cntl_propo > 0 && opt$alt_propo > 0) {
  simulate_experiment(
    fasta = c(fasta_file_A, fasta_file_B),
    num_reps =3, 
    readlen =100, 
    reads_per_transcript = c(total_reads_per_file*opt$cntl_propo, total_reads_per_file*opt$alt_propo),
    fold_changes = c(3,3),
    outdir=alt_outdir,
    paired= TRUE
  )
} else {
  simulate_experiment(
    fasta = fasta_file_B,
    num_reps = 3,
    readlen =100, 
    reads_per_transcript=total_reads_per_file, 
    fold_changes=3,
    outdir=alt_outdir,
    paired = TRUE
)
}


# num_genes_A <- length(readDNAStringSet(fasta_file_A))
# num_genes_B <- length(readDNAStringSet(fasta_file_B))
# seqs_A <- readDNAStringSet(fasta_file_A)
# seqs_B <- readDNAStringSet(fasta_file_B)
# expr_levels <- round(runif(num_genes_A, opt$minC, opt$maxC))


# # fold_changes_changeA <- round(runif(num_genes_A, 3, 4))
# # fold_changes_changeB <- round(runif(num_genes_B, 3, 4))
# # Simulate reads for Condition A
# simulate_experiment(fasta=fasta_file_A, reads_per_transcript=expr_levels, num_reps=3, readlen =100,
#                     outdir=opt$output_dir+"/condition_cntl", paired=TRUE)

# # Simulate reads for Condition B using the same fold changes (no change)
# simulate_experiment(fasta=fasta_file_B, reads_per_transcript=expr_levels, num_reps=3, readlen =100,
#                     outdir=opt$output_dir+"/condition_alt", paired=TRUE)