# now look at full PREDS file
data1 <- read.table( 'lmm-filter/1e-02/preds.txt.gz', header = TRUE )

# test p-values with default p-value threshold
pval <- 1e-2

# check the p-values
# make a copy that excludes untestable loci (presumably fixed?)
bim2 <- data1[ !is.na( data1$pval_fwd ), ]

# keep should have insignificant FWD p-values
stopifnot( all( bim2$pval_fwd[ bim2$category == 'keep' ] >= pval ) )
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
