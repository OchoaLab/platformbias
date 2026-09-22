library(tidyverse)
library(genio)
library(platformbias)
library(optparse) 

# terminal inputs
option_list = list(
  make_option(c( "-f", "--file"), type = "character", default = 'biased', 
              help = "input file with controls from two platforms", metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
# get values

input_file <- opt$f

# takes predictions (disorganized remove and flip lists) for a certain run, and produces a nice classification vector and other useful info, similar to preds.txt.gz for the simpler methods

# assumes we're inside the phase2 dir

# read BIM file to add data to
bim <- read_bim( input_file ) %>% select( chr, id, ref, alt )

# identify reverse complement cases, which are "flippable" (not strictly necessary, though these files I've made usually have them and allow for internal checks)
bim$revcomp <- bim$ref == revcomp( bim$alt )

### LMM DATA ###

# first read final "phase2" edits, including permanently removed and temporarily flipped locus
remove <- read_lines( "phase2_init_remove.txt" )
flip <- read_lines( "phase2_init_flip.txt" )

# check that remove and flip are unique and disjoint
stopifnot( length( remove ) == length( unique( remove ) ) )
stopifnot( length( flip ) == length( unique( flip ) ) )
stopifnot( length( intersect( remove, flip ) ) == 0 )

# phase 2 removals 
phase2_snps <- read_lines( "remove_phase2.txt" )
total_remove <- union( remove, phase2_snps )

# final text classification
# this is most cases
bim$category <- 'keep'
# "flip" will have cases overwritten into "remove" if they were subsequently removed (in phase 3), happens in next step
bim$category[ bim$id %in% flip ] <- 'flip'
bim$category[ bim$id %in% total_remove ] <- 'remove'

# data checks, should pass if all is good
# all indexes are present and ordered continuously, match IDs (big assumption above!)
#stopifnot( all( bim$id == 1 : nrow( bim ) ) )
# confirm that all loci that were flipped are actually flippable
stopifnot( all( bim$revcomp[ bim$category == 'flip' ] ) )

# save updated/extended `bim`, the key calculation!
write_tsv( bim, 'preds.txt.gz' )
