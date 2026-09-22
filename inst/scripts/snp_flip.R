library(readr)
library(genio)
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
remove_snps = read_lines("remove_phase1.txt")  

bim <- read_bim(input_file)
ids_revcomp = bim$id[ bim$ref == revcomp( bim$alt ) ]

sig_snps_flip = intersect( remove_snps, ids_revcomp )
sig_snps_remove = setdiff( remove_snps, ids_revcomp )

# write output and rerun saige
write_lines(sig_snps_flip, "phase2_init_flip.txt")
write_lines(sig_snps_remove, "phase2_init_remove.txt")

print("phase 1 remove SNPs:")
print(length(remove_snps))
print("phase 2 flip SNPs:")
print(length(sig_snps_flip))
print("phase 2 remove SNPs:")
print(length(sig_snps_remove))


