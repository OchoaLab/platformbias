# first look at summary
data1 <- read.table( 'lmm-filter/1e-02/iteration_summary.txt', header = TRUE )
data2 <- read.table( 'lmm-filter_EXPECTED/1e-02/iteration_summary.txt', header = TRUE )
# runtimes are not comparable, but the rest should be identical
data1 <- data1[ , c('phase', 'iter', 'remove') ]
data2 <- data2[ , c('phase', 'iter', 'remove') ]
stopifnot( all( data1 == data2 ) )

# now look at full PREDS file
data1 <- read.table( 'lmm-filter/1e-02/preds.txt.gz', header = TRUE )
data2 <- read.table( 'lmm-filter_EXPECTED/1e-02/preds.txt.gz', header = TRUE )


# test p-values with default p-value threshold
pval <- 1e-2

# check the p-values since they're a new addition
# make a copy that excludes untestable loci (presumably fixed?)
bim2 <- data1[ !is.na( data1$pval_fwd ), ]

# keep should have insignificant FWD p-values
stopifnot( all( bim2$pval_fwd[ bim2$category == 'keep' ] >= pval ) )
# otherwise they are remove, and they must be significant
stopifnot( all( bim2$pval_fwd[ bim2$category == 'remove' ] < pval ) )

# confirm that non-revcomp are never flip
stopifnot( all( bim2$category[ ! bim2$revcomp ] != 'flip' ) )
# and all their REV p-values are NA
stopifnot( all( is.na( bim2$pval_rev[ ! bim2$revcomp ] ) ) )

# now look at revcomps and their REV p-values
bim2 <- bim2[ bim2$revcomp, ]
# now, all keeps have NA REV p-values
stopifnot( all( is.na( bim2$pval_rev[ bim2$category == 'keep' ] ) ) )
# removes all have significant REV p-values
stopifnot( all( bim2$pval_rev[ bim2$category == 'remove' ] < pval ) )
# flips have significant FWD but insignificant REV p-values
stopifnot( all( bim2$pval_fwd[ bim2$category == 'flip' ] < pval ) )
stopifnot( all( bim2$pval_rev[ bim2$category == 'flip' ] >= pval ) )


# now compare to previous outputs
# sadily, these p-values differ by small amounts, I hope by rounding more we might see more agreement
# this appears to be the precision of SAIGE in these tests
digits <- 2
data1$pval_fwd <- signif( data1$pval_fwd, digits )
data1$pval_rev <- signif( data1$pval_rev, digits )
data2$pval_fwd <- signif( data2$pval_fwd, digits )
data2$pval_rev <- signif( data2$pval_rev, digits )

# this shows that NA values agree
stopifnot( all( is.na( data1 ) == is.na( data2 ) ) )
# this shows that non-NA values agree
stopifnot( all( data1 == data2, na.rm = TRUE ) )

