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

# now look at full PREDS file
data1 <- read.table( paste0( dir_out, '/', pval, '/preds.txt.gz' ), header = TRUE )

# pval must be numeric now for comparisons to work
pval <- as.numeric( pval )

# check the p-values
# make a copy that excludes untestable loci (presumably fixed?)
bim2 <- data1[ !is.na( data1$pval_fwd ), ]

# keep should have insignificant FWD p-values
stopifnot( all( bim2$pval_fwd[ bim2$category == 'keep' ] >= pval ) )
if ( !flip ) {
    
    # otherwise they are remove, and must be significant
    stopifnot( all( bim2$pval_fwd[ bim2$category == 'remove' ] < pval ) )
    # confirm that there are no flips
    stopifnot( all( bim2$category != 'flip' ) )
    
} else {
    
    # otherwise they are remove or flip, and non-revcomp cases must be significant
    stopifnot( all( bim2$pval_fwd[ bim2$category != 'keep' & !bim2$revcomp ] < pval ) )

    # confirm that non-revcomp are never flip
    stopifnot( all( bim2$category[ ! bim2$revcomp ] != 'flip' ) )
    # and all their REV p-values are NA
    stopifnot( all( is.na( bim2$pval_rev[ ! bim2$revcomp ] ) ) )

    # now look at revcomps and their REV p-values
    bim2 <- bim2[ bim2$revcomp, ]
    # for the rest, only makes sense to make comparisons when both p-values are present (if we had thresholds, this wouldn't happen)
    bim2 <- bim2[ !is.na( bim2$pval_fwd ) & !is.na( bim2$pval_rev ), ]
    # removes all have significant REV p-values
    stopifnot( all( bim2$pval_rev[ bim2$category == 'remove' ] < pval ) )
    # keeps have bigger FWD p-values than REV, when both are present
    # keep also must have insignificant FWD p-values (already tested above for all cases)
    stopifnot( all( bim2$pval_fwd[ bim2$category == 'keep' ] >= bim2$pval_rev[ bim2$category == 'keep' ] ) )
    # the reverse is true for flip
    stopifnot( all( bim2$pval_fwd[ bim2$category == 'flip' ] < bim2$pval_rev[ bim2$category == 'flip' ] ) )
    stopifnot( all( bim2$pval_rev[ bim2$category == 'flip' ] >= pval ) )

}
