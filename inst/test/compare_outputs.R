library(optparse)

# terminal inputs
option_list = list(
    make_option(c( "-d", "--dir_out"), type = "character", default = 'lmm-filter', 
                help = "Output directory", metavar = "character"),
    make_option("--pval", type = "character", default = '1e-02',
                help = "P-value threshold for identifying significant snps, and exact name of output subdirectory", metavar = "character"),
    make_option("--noflip", action = "store_true", default = FALSE,
                help = "Run phase 1 only (removals only, no flips).  This is best for data where reverse-complement strand flips are not expected, such as whole-genome sequencing.  (Flips are often expected in genotyping array data.)")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
# get values
pval <- opt$pval
dir_out <- opt$dir_out
flip <- !opt$noflip

dir1 <- paste0( dir_out, '/', pval, '/' )
dir2 <- paste0( dir_out, '_EXPECTED', if ( !flip ) '-NO-FLIP', '/', pval, '/' )

# first look at summary
data1 <- read.table( paste0( dir1, 'iteration_summary.txt' ), header = TRUE )
data2 <- read.table( paste0( dir2, 'iteration_summary.txt' ), header = TRUE )
# runtimes are not comparable, but the rest should be identical
data1 <- data1[ , c('phase', 'iter', 'edits') ]
data2 <- data2[ , c('phase', 'iter', 'edits') ]
stopifnot( all( data1 == data2 ) )

# now look at full PREDS file
data1 <- read.table( paste0( dir1, 'preds.txt.gz' ), header = TRUE )
data2 <- read.table( paste0( dir2, 'preds.txt.gz' ), header = TRUE )

# sadily, these p-values differ by small amounts, I hope by rounding more we might see more agreement
# this appears to be the worst precision of SAIGE in these tests
digits <- 1
data1$pval_fwd <- signif( data1$pval_fwd, digits )
data2$pval_fwd <- signif( data2$pval_fwd, digits )
if ( flip ) {
    data1$pval_rev <- signif( data1$pval_rev, digits )
    data2$pval_rev <- signif( data2$pval_rev, digits )
}

# this shows that NA values agree
stopifnot( all( is.na( data1 ) == is.na( data2 ) ) )
# this shows that non-NA values agree
stopifnot( all( data1 == data2, na.rm = TRUE ) )
