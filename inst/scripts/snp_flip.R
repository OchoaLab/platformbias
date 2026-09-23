library(optparse) 
library(platformbias)

# terminal inputs
option_list = list(
  make_option(c( "-f", "--file"), type = "character", default = 'biased', 
              help = "input file with controls from two platforms", metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
# get values

input_file <- opt$f

###########################################
# identify significant SNPs in current iteration

# identify SNPs that can be flipped
remove_snps = readLines("remove_phase1.txt")  

bim <- read.table( paste0( input_file, '.bim' ), header = FALSE )
colnames( bim ) <- c('chr', 'id', 'posg', 'pos', 'alt', 'ref')
ids_revcomp = bim$id[ bim$ref == revcomp( bim$alt ) ]

sig_snps_flip = intersect( remove_snps, ids_revcomp )
sig_snps_remove = setdiff( remove_snps, ids_revcomp )

# write output and rerun saige
writeLines(sig_snps_flip, "phase2_init_flip.txt")
writeLines(sig_snps_remove, "phase2_init_remove.txt")
