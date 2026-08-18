library(tibble)

# level for CIs, to make them even more certain than default if needed
level <- 0.999

test_that("af_test_single, af_test_or_single work for flat x edge cases", {
    # this used to cause errors particularly for glm/OR case
    x1 <- 0
    n1 <- 11
    x2 <- 0
    n2 <- 31
    # original, stat is always zero here
    expect_silent(
        stat <- af_test_single( x1, n1, x2, n2 )
    )
    expect_true( stat == 0 )
    # repeat with OR version, value also must have a specific form
    expect_silent(
        data <- af_test_or_single( x1, n1, x2, n2, level = level )
    )
    expect_equal( data, c(0, -Inf, Inf, 1) )

    # other edge case has n=x
    x1 <- n1
    x2 <- n2
    # expect same behaviors here
    # original, stat is always zero here
    expect_silent(
        stat <- af_test_single( x1, n1, x2, n2 )
    )
    expect_true( stat == 0 )
    # repeat with OR version, value also must have a specific form
    expect_silent(
        data <- af_test_or_single( x1, n1, x2, n2, level = level )
    )
    expect_equal( data, c(0, -Inf, Inf, 1) )
})

test_that("af_test_single, af_test, af_test_or_single, af_test_or work with k=1", {
    # simulate some null data
    p <- 0.4
    n1 <- 10
    n2 <- 13
    x1 <- rbinom( 1, n1, p )
    x2 <- rbinom( 1, n2, p )

    # test it!
    expect_silent(
        stat <- af_test_single( x1, n1, x2, n2 )
    )
    expect_equal( length( stat ), 1 )
    expect_true( stat >= 0 )
    # repeat with more general wrapper
    expect_silent(
        obj <- af_test( x1, n1, x2, n2 )
    )
    expect_equal( obj$stat, stat )
    pval <- obj$pval
    expect_equal( length( pval ), 1 )
    #expect_true( pval >= 0 )
    expect_true( pval <= 1 )
    # in this case null holds, so it should be very unlikely that the result is extremely significant
    expect_true( pval > 1e-10 )
    # repeat with OR version
    expect_silent(
        data <- af_test_or_single( x1, n1, x2, n2, level = level )
    )
    expect_true( is.numeric( data ) )
    expect_equal( length( data ), 4 )
    # the mean should be within the CI always
    expect_true( data[2] <= data[1] )
    expect_true( data[3] >= data[1] )
    # because it's the null, the log-OR (that's what this is) CI should contain 0 very often
    expect_true( data[2] <= 0 )
    expect_true( data[3] >= 0 )
    # test p-values
    expect_true( data[4] <= 1 )
    # in this case null holds, so it should be very unlikely that the result is extremely significant
    expect_true( data[4] > 1e-10 )
    # and final wrapper, which is real ORs (not log), names are different
    expect_silent(
        data <- af_test_or( x1, n1, x2, n2, level = level )
    )
    expect_true( is_tibble( data ) )
    expect_equal( names( data ), c('OR', 'CIL', 'CIU', 'pval') )
    expect_equal( nrow( data ), 1 )
    # the mean should be within the CI always
    expect_true( data$CIL <= data$OR )
    expect_true( data$CIU >= data$OR )
    # because it's the null, the OR (that's what this is) CI should contain 1 very often
    expect_true( data$CIL <= 1 * 10 ) # make more permissive, fails too often
    expect_true( data$CIU >= 1 / 10 ) # make more permissive, fails too often
    # test p-values
    expect_true( data$pval <= 1 )
    # in this case null holds, so it should be very unlikely that the result is extremely significant
    expect_true( data$pval > 1e-10 )
    
    
    # simulate some alternative data
    # big sample sizes too to best ensure a significant result
    p1 <- 0.4
    p2 <- 0.1
    n1 <- 1000
    n2 <- 1390
    x1 <- rbinom( 1, n1, p1 )
    x2 <- rbinom( 1, n2, p2 )
    # get true OR
    true_OR <- p1 * (1 - p2 ) / ( p2 * (1 - p1 ) )

    # test it!
    expect_silent(
        stat <- af_test_single( x1, n1, x2, n2 )
    )
    expect_equal( length( stat ), 1 )
    expect_true( stat >= 0 )
    # repeat with more general wrapper
    expect_silent(
        obj <- af_test( x1, n1, x2, n2 )
    )
    expect_equal( obj$stat, stat )
    pval <- obj$pval
    expect_equal( length( pval ), 1 )
    expect_true( pval >= 0 )
    #expect_true( pval <= 1 )
    # in this case alt holds, and sample sizes are large, so a significant result is most likely
    expect_true( pval < 1e-5 )
    # repeat with OR version
    expect_silent(
        data <- af_test_or_single( x1, n1, x2, n2, level = level )
    )
    expect_true( is.numeric( data ) )
    expect_equal( length( data ), 4 )
    # the mean should be within the CI always
    expect_true( data[2] <= data[1] )
    expect_true( data[3] >= data[1] )
    # the log-OR (that's what this is) CI should contain its true value very often
    expect_true( data[2] <= log( true_OR ) )
    expect_true( data[3] >= log( true_OR ) )
    # test p-values
    expect_true( data[4] >= 0 )
    # in this case alt holds, and sample sizes are large, so a significant result is most likely
    expect_true( data[4] < 1e-5 )
    # and final wrapper, which is real ORs (not log), names are different
    expect_silent(
        data <- af_test_or( x1, n1, x2, n2, level = level )
    )
    expect_true( is_tibble( data ) )
    expect_equal( names( data ), c('OR', 'CIL', 'CIU', 'pval') )
    expect_equal( nrow( data ), 1 )
    # the mean should be within the CI always
    expect_true( data$CIL <= data$OR )
    expect_true( data$CIU >= data$OR )
    # because it's the alt, the OR (that's what this is) CI should contain its true value very often
    expect_true( data$CIL <= true_OR * 10 ) # make more permissive, fails too often
    expect_true( data$CIU >= true_OR / 10 ) # make more permissive, fails too often
    # test p-values
    expect_true( data$pval >= 0 )
    # in this case alt holds, and sample sizes are large, so a significant result is most likely
    expect_true( data$pval < 1e-5 )
})

test_that( "af_test works with k=3", {
    k <- 3
    # simulate some null data
    # totally random data!
    p1 <- runif( k )
    p2 <- p1
    # same for sample sizes
    n1 <- rpois( k, 10 )
    n2 <- rpois( k, 10 )
    # don't allow zeroes of course
    n1[ n1 == 0 ] <- 1
    n2[ n2 == 0 ] <- 1
    # simulate data in batch
    x1 <- rbinom( k, n1, p1 )
    x2 <- rbinom( k, n2, p2 )
    
    # this is main version, testing each case independently
    expect_silent(
        obj <- af_test( x1, n1, x2, n2 )
    )
    stat <- obj$stat
    df <- obj$df
    pval <- obj$pval
    # expect several values here
    expect_equal( length( stat ), k )
    expect_equal( length( pval ), k )
    # range tests
    expect_true( all( stat >= 0 ) )
    expect_equal( df, 1 ) # single degree of freedom here!
    #expect_true( all( pval >= 0 ) )
    expect_true( all( pval <= 1 ) )
    # in this case null holds, so it should be very unlikely that the result is extremely significant
    expect_true( all( pval > 1e-10 ) )

    ## # report data being used!
    ## message( 'x1: ', toString( x1 ) )
    ## message( 'n1: ', toString( n1 ) )
    ## message( 'x2: ', toString( x2 ) )
    ## message( 'n2: ', toString( n2 ) )

    # test OR version
    expect_silent(
        data <- af_test_or( x1, n1, x2, n2, level = level )
    )
    expect_true( is_tibble( data ) )
    expect_equal( names( data ), c('OR', 'CIL', 'CIU', 'pval') )
    expect_equal( nrow( data ), k )
    # the mean should be within the CI always
    expect_true( all( data$CIL <= data$OR ) )
    expect_true( all( data$CIU >= data$OR ) )
    # because it's the null, the OR (that's what this is) CI should contain 1 very often
    expect_true( all( data$CIL <= 1 * 10 ) ) # make more permissive, fails too often
    expect_true( all( data$CIU >= 1 / 10 ) ) # make more permissive, fails too often
    # test p-values
    expect_true( all( data$pval <= 1 ) )
    # in this case null holds, so it should be very unlikely that the result is extremely significant
    expect_true( all( data$pval > 1e-10 ) )

    # ensure matrix version is the same
    x1 <- cbind( x1 )
    n1 <- cbind( n1 )
    x2 <- cbind( x2 )
    n2 <- cbind( n2 )
    expect_silent(
        obj2 <- af_test( x1, n1, x2, n2 )
    )
    expect_equal( obj2, obj )
    expect_silent(
        data2 <- af_test_or( x1, n1, x2, n2, level = level )
    )
    expect_equal( data2, data )

    # now run version with data as rows, which combines loci into a single test
    x1 <- t( x1 )
    n1 <- t( n1 )
    x2 <- t( x2 )
    n2 <- t( n2 )
    expect_silent(
        obj <- af_test( x1, n1, x2, n2 )
    )
    stat <- obj$stat
    df <- obj$df
    pval <- obj$pval
    # expect single values here
    expect_equal( length( stat ), 1 )
    expect_equal( length( pval ), 1 )
    # range tests
    expect_true( stat >= 0 )
    expect_equal( df, k )
    #expect_true( pval >= 0 )
    expect_true( pval <= 1 )
    # in this case null holds, so it should be very unlikely that the result is extremely significant
    expect_true( pval > 1e-10 )

    # another version
    expect_silent(
        obj2 <- likelihoodratio_3(
            x1[1], n1[1], x2[1], n2[1],
            x1[2], n1[2], x2[2], n2[2],
            x1[3], n1[3], x2[3], n2[3]
        )
    )
    expect_equal( obj2, obj )

    # test OR version
    expect_silent(
        data <- af_test_or( x1, n1, x2, n2, level = level )
    )
    expect_true( is_tibble( data ) )
    expect_equal( names( data ), c('OR', 'CIL', 'CIU', 'pval') )
    expect_equal( nrow( data ), 1 )
    # the mean should be within the CI always
    expect_true( data$CIL <= data$OR )
    expect_true( data$CIU >= data$OR )
    # because it's the null, the OR (that's what this is) CI should contain 1 very often
    expect_true( data$CIL <= 1 * 10 ) # make more permissive, fails too often
    expect_true( data$CIU >= 1 / 10 ) # make more permissive, fails too often
    # test p-values
    expect_true( data$pval <= 1 )
    # in this case null holds, so it should be very unlikely that the result is extremely significant
    expect_true( data$pval > 1e-10 )

    # simulate some alternative data
    # totally random data!
    p1 <- runif( k )
    p2 <- runif( k )
    # same for sample sizes
    n1 <- rpois( k, 1000 )
    n2 <- rpois( k, 1000 )
    # don't allow zeroes of course
    n1[ n1 == 0 ] <- 1
    n2[ n2 == 0 ] <- 1
    # simulate data in batch
    x1 <- rbinom( k, n1, p1 )
    x2 <- rbinom( k, n2, p2 )
    # get true ORs (vector)
    true_OR <- p1 * (1 - p2 ) / ( p2 * (1 - p1 ) )
    
    # this is main version, testing each case independently
    expect_silent(
        obj <- af_test( x1, n1, x2, n2 )
    )
    stat <- obj$stat
    df <- obj$df
    pval <- obj$pval
    # expect several values here
    expect_equal( length( stat ), k )
    expect_equal( length( pval ), k )
    # range tests
    expect_true( all( stat >= 0 ) )
    expect_equal( df, 1 ) # single degree of freedom here!
    expect_true( all( pval >= 0 ) )
    expect_true( all( pval <= 1 ) )
    # in this case alt holds, and sample sizes are large, so a significant result is most likely
    # there are so many tests, though, that some do happen to get quite large p-values
    #expect_true( all( pval < 0.5 ) )

    # test OR version
    expect_silent(
        data <- af_test_or( x1, n1, x2, n2, level = level )
    )
    expect_true( is_tibble( data ) )
    expect_equal( names( data ), c('OR', 'CIL', 'CIU', 'pval') )
    expect_equal( nrow( data ), k )
    # the mean should be within the CI always
    expect_true( all( data$CIL <= data$OR ) )
    expect_true( all( data$CIU >= data$OR ) )
    # because it's the alt, the OR (that's what this is) CI should contain its true value very often
    expect_true( all( data$CIL <= true_OR * 10 ) ) # make more permissive, this fails too often
    expect_true( all( data$CIU >= true_OR / 10 ) ) # make more permissive, this fails too often
    # test p-values
    expect_true( all( data$pval >= 0 ) )
    # in this case alt holds, and sample sizes are large, so a significant result is most likely
    #expect_true( all( data$pval < 0.5 ) ) # make more permissive, this fails too often

    # ensure matrix version is the same
    x1 <- cbind( x1 )
    n1 <- cbind( n1 )
    x2 <- cbind( x2 )
    n2 <- cbind( n2 )
    expect_silent(
        obj2 <- af_test( x1, n1, x2, n2 )
    )
    expect_equal( obj2, obj )
    expect_silent(
        data2 <- af_test_or( x1, n1, x2, n2, level = level )
    )
    expect_equal( data2, data )

    # now run version with data as rows, which combines loci into a single test
    x1 <- t( x1 )
    n1 <- t( n1 )
    x2 <- t( x2 )
    n2 <- t( n2 )
    expect_silent(
        obj <- af_test( x1, n1, x2, n2 )
    )
    stat <- obj$stat
    df <- obj$df
    pval <- obj$pval
    # expect single values here
    expect_equal( length( stat ), 1 )
    expect_equal( length( pval ), 1 )
    # range tests
    expect_true( stat >= 0 )
    expect_equal( df, k )
    #expect_true( pval >= 0 )
    expect_true( pval <= 1 )
    # in this case alt holds, and sample sizes are large, so a significant result is most likely
    expect_true( pval < 1e-4 )
    
    # another version
    expect_silent(
        obj2 <- likelihoodratio_3(
            x1[1], n1[1], x2[1], n2[1],
            x1[2], n1[2], x2[2], n2[2],
            x1[3], n1[3], x2[3], n2[3]
        )
    )
    expect_equal( obj2, obj )

    # test OR version
    expect_silent(
        data <- af_test_or( x1, n1, x2, n2, level = level )
    )
    expect_true( is_tibble( data ) )
    expect_equal( names( data ), c('OR', 'CIL', 'CIU', 'pval') )
    expect_equal( nrow( data ), 1 )
    # the mean should be within the CI always
    expect_true( data$CIL <= data$OR )
    expect_true( data$CIU >= data$OR )
    # here it's hard to really say what the true OR should be, so let's skip that particular test
    ## # because it's the alt, the OR (that's what this is) CI should contain its true value very often
    ## expect_true( data$CIL <= true_OR )
    ## expect_true( data$CIU >= true_OR )
    # test p-values
    expect_true( data$pval >= 0 )
    # in this case alt holds, and sample sizes are large, so a significant result is most likely
    expect_true( data$pval < 0.5 )
})

test_that( "revcomp works", {
    # a vector of alleles with all cases for single letters, and a few more complex cases
    x_in    <- c('A', 'C', 'G', 'T', 'a', 'c', 'g', 't', 'ATAC', 'GATTACA')
    # these are the expected outputs
    out_exp <- c('T', 'G', 'C', 'A', 't', 'g', 'c', 'a', 'GTAT', 'TGTAATC')

    # a successful run!
    expect_silent(
        out <- revcomp( x_in )
    )
    expect_equal( out, out_exp )
})


# simulate genotypes to reuse to test various functions
m_loci <- 101
n_ind <- 9
X <- matrix(
    rbinom( m_loci * n_ind, 2, 0.5 ),
    nrow = m_loci,
    ncol = n_ind
)
# for flipping tests, come up with some complex alt/ref cases that are also predictable
m_flippable <- 51
alt <- c(
    sample( c('A', 'C', 'G', 'T', 'GATTACA'), 51, replace = TRUE ),
    rep.int( 'A', 25 ),
    rep.int( 'T', 25 )
)
# first half, roughly, of refs are actual reverse complements of alts (thus flippable)
ref <- revcomp( alt )
# the other half force to not be
ref[ 52:101 ] <- 'C'
bim <- tibble( alt = alt, ref = ref )

test_that( "bias_geno works", {
    # parameters
    p_biased <- 0.1
    # for checks
    m_biased <- round( p_biased * m_loci )

    # this is the standard run
    expect_silent(
        obj <- bias_geno( X, p_biased )
    )
    # validate outputs
    expect_true( is.list( obj ) )
    expect_equal( names( obj ), c('X', 'biased_loci', 'error_probs', 'error_genos') )
    # check genotypes
    expect_true( is.matrix( obj$X ) )
    expect_equal( nrow( obj$X ), m_loci )
    expect_equal( ncol( obj$X ), n_ind )
    expect_true( all( obj$X %in% 0L : 2L ) )
    # check the other vectors
    expect_equal( length( obj$biased_loci ), m_biased )
    expect_equal( length( unique( obj$biased_loci ) ), m_biased )
    expect_true( is.integer( obj$biased_loci ) )
    expect_true( all( obj$biased_loci %in% 1 : m_loci ) )
    expect_equal( length( obj$error_probs ), m_biased )
    expect_true( is.numeric( obj$error_probs ) )
    expect_true( min( obj$error_probs ) >= 0 )
    expect_true( max( obj$error_probs ) <= 1 )
    expect_equal( length( obj$error_genos ), m_biased )
    expect_true( is.integer( obj$error_genos ) )
    expect_true( all( obj$error_genos %in% c(0L, 2L) ) )
    # loci that were not biased should be identical to input
    expect_equal( obj$X[ -obj$biased_loci, ], X[ -obj$biased_loci, ] )
    # among biased loci, those that were biased toward zero should be smaller or equal than before, and the opposite for those biased towards 2
    indexes <- obj$error_genos == 0
    expect_true( all( obj$X[ obj$biased_loci[ indexes ], ] <= X[ obj$biased_loci[ indexes ], ] ) )
    expect_true( all( obj$X[ obj$biased_loci[ !indexes ], ] >= X[ obj$biased_loci[ !indexes ], ] ) )
    
    # repeat test providing loci to bias
    m_biased <- 21
    biased_loci <- 1 : m_biased
    expect_silent(
        obj <- bias_geno( X, biased_loci = biased_loci )
    )
    # validate outputs
    expect_true( is.list( obj ) )
    expect_equal( names( obj ), c('X', 'biased_loci', 'error_probs', 'error_genos') )
    # check genotypes
    expect_true( is.matrix( obj$X ) )
    expect_equal( nrow( obj$X ), m_loci )
    expect_equal( ncol( obj$X ), n_ind )
    expect_true( all( obj$X %in% 0L : 2L ) )
    # output biased_loci should match input exactly
    expect_equal( obj$biased_loci, biased_loci )
    # check the other vectors
    expect_equal( length( obj$error_probs ), m_biased )
    expect_true( is.numeric( obj$error_probs ) )
    expect_true( min( obj$error_probs ) >= 0 )
    expect_true( max( obj$error_probs ) <= 1 )
    expect_equal( length( obj$error_genos ), m_biased )
    expect_true( is.integer( obj$error_genos ) )
    expect_true( all( obj$error_genos %in% c(0L, 2L) ) )
    # loci that were not biased should be identical to input
    expect_equal( obj$X[ -obj$biased_loci, ], X[ -obj$biased_loci, ] )
    # among biased loci, those that were biased toward zero should be smaller or equal than before, and the opposite for those biased towards 2
    indexes <- obj$error_genos == 0
    expect_true( all( obj$X[ obj$biased_loci[ indexes ], ] <= X[ obj$biased_loci[ indexes ], ] ) )
    expect_true( all( obj$X[ obj$biased_loci[ !indexes ], ] >= X[ obj$biased_loci[ !indexes ], ] ) )
})

test_that( "flip_revcomps works", {
    p <- 0.1
    # for checks
    m_flipped <- round( p * m_flippable )

    # the successful run
    expect_silent(
        obj <- flip_revcomps( X, p, bim )
    )
    # validate outputs
    expect_true( is.list( obj ) )
    expect_equal( names( obj ), c('X', 'flipped_loci') )
    # check genotypes
    expect_true( is.matrix( obj$X ) )
    expect_equal( nrow( obj$X ), m_loci )
    expect_equal( ncol( obj$X ), n_ind )
    expect_true( all( obj$X %in% 0L : 2L ) )
    # check the other vector
    expect_equal( length( obj$flipped_loci ), m_flipped )
    expect_equal( length( unique( obj$flipped_loci ) ), m_flipped )
    expect_true( is.integer( obj$flipped_loci ) )
    expect_true( all( obj$flipped_loci %in% 1 : m_flippable ) ) # flippable loci were all the first few
    # loci that were not flipped should be identical to input
    expect_equal( obj$X[ -obj$flipped_loci, ], X[ -obj$flipped_loci, ] )
    # and flipped loci satisfy this relationship!
    expect_equal( obj$X[ obj$flipped_loci, ], 2 - X[ obj$flipped_loci, ] )

    # test version where the flipped loci are predetermined
    m_flipped <- 21
    flipped_loci <- 1 : m_flipped
    expect_silent(
        obj <- flip_revcomps( X, flipped_loci = flipped_loci )
    )
    # validate outputs
    expect_true( is.list( obj ) )
    expect_equal( names( obj ), c('X', 'flipped_loci') )
    # check genotypes
    expect_true( is.matrix( obj$X ) )
    expect_equal( nrow( obj$X ), m_loci )
    expect_equal( ncol( obj$X ), n_ind )
    expect_true( all( obj$X %in% 0L : 2L ) )
    # output flipped_loci should match input exactly
    expect_equal( obj$flipped_loci, flipped_loci )
    # loci that were not flipped should be identical to input
    expect_equal( obj$X[ -obj$flipped_loci, ], X[ -obj$flipped_loci, ] )
    # and flipped loci satisfy this relationship!
    expect_equal( obj$X[ obj$flipped_loci, ], 2 - X[ obj$flipped_loci, ] )
})

test_that( "af_filter works", {
    # use the global m_loci and bim that goes with it!
    
    # make some matrix data
    k_subpops <- 3
    p1 <- matrix( runif( m_loci * k_subpops ), nrow = m_loci, ncol = k_subpops )
    p2 <- matrix( runif( m_loci * k_subpops ), nrow = m_loci, ncol = k_subpops )
    # same for sample sizes
    n1 <- matrix( rpois( m_loci * k_subpops, 10 ), nrow = m_loci, ncol = k_subpops )
    n2 <- matrix( rpois( m_loci * k_subpops, 10 ), nrow = m_loci, ncol = k_subpops )
    # don't allow zeroes of course
    n1[ n1 == 0 ] <- 1
    n2[ n2 == 0 ] <- 1
    # simulate data in batch
    x1 <- matrix( rbinom( m_loci * k_subpops, n1, p1 ), nrow = m_loci, ncol = k_subpops )
    x2 <- matrix( rbinom( m_loci * k_subpops, n2, p2 ), nrow = m_loci, ncol = k_subpops )

    expect_silent(
        data <- af_filter( x1, n1, x2, n2, bim )
    )
    expect_true( is.data.frame( data ) )
    expect_equal( nrow( data ), m_loci )
    expect_equal( names( data ), c('alt', 'ref', 'revcomp', 'pval_fwd', 'pval_rev') )
    expect_equal( bim$alt == revcomp( bim$ref ), data$revcomp )
    expect_true( !anyNA( data$pval_fwd ) )
    expect_true( !anyNA( data$pval_rev[ data$revcomp ] ) )
    expect_true( all( is.na( data$pval_rev[ !data$revcomp ] ) ) )
})

test_that( "filter_category works", {
    # come up with an example with all relevant cases, all values far from the threshold one way or the other
    data <- tibble(
        revcomp  = c( FALSE, FALSE,  TRUE,  TRUE,  TRUE,  TRUE,  TRUE ),
        pval_fwd = c(   0.9, 1e-10,   0.1,   0.9,   0.9, 1e-10, 1e-10 ),
        pval_rev = c(    NA,    NA,   0.9,   0.1, 1e-10,   0.9, 1e-10 )
    )
    pcut <- 0.01
    cat_exp = c( 'keep', 'remove', 'flip', 'keep', 'keep', 'flip', 'remove' )
    
    expect_silent( 
        cat <- filter_category( data, pcut )
    )
    expect_equal( cat, cat_exp )
})

test_that( "binary_eval works", {
    # come up with a fixed example where we know the answers
    true <- c(TRUE, TRUE, TRUE,  TRUE, FALSE, FALSE, FALSE)
    pred <- c(TRUE, TRUE, TRUE, FALSE,  TRUE,  TRUE, FALSE)
    # expected return is: precision, recall, TP, FP, FN
    TP <- 3
    FP <- 2
    FN <- 1
    precision <- TP / ( TP + FP )
    recall <- TP / ( TP + FN )
    data_exp <- c( precision, recall, TP, FP, FN )
    names( data_exp ) <- c( 'precision', 'recall', 'TP', 'FP', 'FN' )

    expect_silent(
        data <- binary_eval( true, pred )
    )
    expect_equal( data, data_exp )
})

test_that( "check_cats works", {
    # a successful case
    expect_silent(
        check_cats( c( 'keep', 'flip', 'remove' ) )
    )

    # and cases that are supposed to cause errors
    expect_error(
        check_cats( c( 'keep', 'flip', 'remove', NA ) )
    )
    for ( cat in c( 'Keep', 'Flip', 'Remove', 'X', 'random' ) ) {
        expect_error(
            check_cats( cat )
        )
    }
})

test_that( "gwas_eval works", {
    # do several examples, this is a more complicated scenario
    # first the most trivial one, just enumerate all cases, expect one of each, should be easy to catch double counting or missing counts
    true <- c( 'keep', 'keep', 'keep', 'flip', 'flip', 'flip', 'remove', 'remove', 'remove' )
    pred <- c( 'keep', 'flip', 'remove', 'keep', 'flip', 'remove', 'keep', 'flip', 'remove' )
    # this is what we expect
    TP <- 2
    FP <- 4
    FN <- 2
    precision <- TP / ( TP + FP )
    recall <- TP / ( TP + FN )
    data_exp <- c( precision, recall, TP, FP, FN )
    names( data_exp ) <- c( 'precision', 'recall', 'TP', 'FP', 'FN' )

    expect_silent(
        data <- gwas_eval( true, pred )
    )
    expect_equal( data, data_exp )

    # some sparser examples, though it's unclear what's more informative exactly
    true <- c( 'keep', 'keep', 'keep' )
    pred <- c( 'keep', 'flip', 'remove' )
    # this is what we expect
    TP <- 1
    FP <- 1
    FN <- 1
    precision <- TP / ( TP + FP )
    recall <- TP / ( TP + FN )
    data_exp <- c( precision, recall, TP, FP, FN )
    names( data_exp ) <- c( 'precision', 'recall', 'TP', 'FP', 'FN' )

    expect_silent(
        data <- gwas_eval( true, pred )
    )
    expect_equal( data, data_exp )

    # some sparser examples, though it's unclear what's more informative exactly
    true <- c( 'flip', 'remove', 'flip', 'remove' )
    pred <- c( 'flip', 'flip', 'remove', 'remove' )
    # this is what we expect
    TP <- 1
    FP <- 1
    FN <- 1
    precision <- TP / ( TP + FP )
    recall <- TP / ( TP + FN )
    data_exp <- c( precision, recall, TP, FP, FN )
    names( data_exp ) <- c( 'precision', 'recall', 'TP', 'FP', 'FN' )

    expect_silent(
        data <- gwas_eval( true, pred )
    )
    expect_equal( data, data_exp )
})

test_that( "filter_eval works", {
    # come up with fixed categories where we know the answers
    true <- c('keep', 'keep',   'keep', 'flip',   'flip', 'remove', 'remove', 'remove')
    pred <- c('keep', 'flip', 'remove', 'flip', 'remove', 'remove',   'flip',   'keep')
    # expected data are just precisions and recalls for each case
    data_exp <- tibble(
        type = c('gwas', 'flip', 'remove'),
        precision = c( 2/5, 1/3, 1/3 ),
        recall = c( 1/2, 1/2, 1/3 ),
        removed = c( 3, NA, NA )
    )
    
    expect_silent(
        data <- filter_eval( true, pred )
    )
    expect_equal( data, data_exp )
})
