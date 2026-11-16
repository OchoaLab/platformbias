#' Calculate locus categories depending on threshold
#'
#' This takes the output of [af_filter()] or an equivalent function and classifies loci depending on the pattern of significance, which is threshold-dependent.
#' For loci that only have forward p-values, the are classified as "remove" if they are significant (below the threshold) or "keep" otherwise.
#' For loci with both forward and reversed p-values (which should have reversed complement ref/alt alleles, a warning is issued if this is not true), they are classified as "remove" if both p-values are significant, "flip" if the forward p-value is the smallest, and "keep" otherwise.
#' Loci with `NA` forward p-values are always "keep" (even if the reverse p-value is not `NA`).
#'
#' @param data A tibble/data.frame with at least three columns:
#' - `revcomp`: Logical, whether locus is identified as having reverse complement (flippale) alleles.
#' - `pval_fwd`: p-values in given orientation.
#' - `pval_rev`: p-values in reversed orientation for the second dataset, for loci with `revcomp == TRUE`.
#' @param pcut The p-value threshold
#'
#' @return A vector with length equal to the number of rows of `data`, with one of 3 possible values depending on classification outcome per locus: "remove", "flip", or "keep".
#'
#' @examples
#' library(tibble)
#' # an example with all cases
#' data <- tibble(
#'     revcomp  = c( FALSE, FALSE,  TRUE,  TRUE,  TRUE,  TRUE,  TRUE ),
#'     pval_fwd = c(   0.9, 1e-10,   0.1,   0.9,   0.9, 1e-10, 1e-10 ),
#'     pval_rev = c(    NA,    NA,   0.9,   0.1, 1e-10,   0.9, 1e-10 )
#' )
#' pcut <- 0.01
#'
#' # get categories!
#' cat <- filter_category( data, pcut )
#'
#' @seealso
#' [af_filter()]
#' 
#' @export
filter_category <- function( data, pcut ) {
    # when this is a string, bad things happen, force numeric now!
    pcut <- as.numeric( pcut )

    # make output, defaulting to keep when all p-values are NA
    x <- rep.int( 'keep', nrow( data ) )
    
    # do the two cases separately, accounting for NA p-values which force some cases to switch (if we went by `revcomp` status exclusively)
    
    # first the cases handled as being forward only
    # this can include `revcomp` loci that haven't had both p-values calculated (lmm-filter early phases), so there's no need to include `revcomp` in the decision
    indexes <- !is.na( data$pval_fwd ) & is.na( data$pval_rev )
    # copy down values
    pvals <- data$pval_fwd[ indexes ]
    x[ indexes ] <- ifelse( pvals < pcut , "remove", "keep")

    # last case requires both p-values to be non-NA
    indexes <- !is.na( data$pval_fwd ) & !is.na( data$pval_rev )
    # shouldn't happen, let's not make it fatal but let users know that it is unexpected
    if ( !all( data$revcomp[ indexes ] ) )
        warning( 'There are loci with both FWD and REV p-values that are not revcomp!' )
    # copy down values
    pval_fwd <- data$pval_fwd[ indexes ]
    pval_rev <- data$pval_rev[ indexes ]
    x[ indexes ] <- ifelse( pval_fwd < pcut & pval_rev < pcut, "remove",
                    ifelse( pval_fwd < pval_rev, "flip", "keep" ) )

    # return category vector
    return( x )
}
