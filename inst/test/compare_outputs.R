# first look at summary
data1 <- read.table( 'lmm-filter/1e-02/iteration_summary.txt', header = TRUE )
data2 <- read.table( 'lmm-filter_EXPECTED/1e-02/iteration_summary.txt', header = TRUE )
# runtimes are not comparable, but the rest should be identical
data1 <- data1[ , c('phase', 'iter', 'edits') ]
data2 <- data2[ , c('phase', 'iter', 'edits') ]
stopifnot( all( data1 == data2 ) )

# now look at full PREDS file
data1 <- read.table( 'lmm-filter/1e-02/preds.txt.gz', header = TRUE )
data2 <- read.table( 'lmm-filter_EXPECTED/1e-02/preds.txt.gz', header = TRUE )

# sadily, these p-values differ by small amounts, I hope by rounding more we might see more agreement
# this appears to be the worst precision of SAIGE in these tests
digits <- 1
data1$pval_fwd <- signif( data1$pval_fwd, digits )
data1$pval_rev <- signif( data1$pval_rev, digits )
data2$pval_fwd <- signif( data2$pval_fwd, digits )
data2$pval_rev <- signif( data2$pval_rev, digits )

# this shows that NA values agree
stopifnot( all( is.na( data1 ) == is.na( data2 ) ) )
# this shows that non-NA values agree
stopifnot( all( data1 == data2, na.rm = TRUE ) )
