library(platformbias) 
library(optparse)

# main function/loop
# initial p-value is string because we want folder name to be exactly this ("1e-02" instead of "0.01"), but after dir is made we can turn to numeric
lmm_filter <- function(
                       input_data,
                       platform_file,
                       pval = '1e-02',
                       dir_out = 'lmm-filter',
                       flip = TRUE,
                       iid = 'IID',
                       fid = 'FID',
                       platform_col = 'PLATFORM'
                       ) {
    # validate inputs
    if ( missing( input_data ) || is.na( input_data ) )
        stop( '`input_data` is required!' )
    if ( missing( platform_file ) || is.na( platform_file ) )
        stop( '`platform_file` is required!' )
    # require that inputs exist
    for ( ext in c('bed', 'bim', 'fam') ) {
        file_in <- paste0( input_data, '.', ext )
        if ( !file.exists( file_in ) )
            stop( 'Input does not exist: ', file_in )
    }
    if ( !file.exists( platform_file ) )
        stop( 'Input does not exist: ', platform_file )
    
    # force inputs to have absolute paths
    # for `input_data`, need to add extension (so file exists, otherwise this doesn't work), then take it back off
    input_data <- sub( '.bed$', '', normalizePath( paste0( input_data, '.bed' ) ) )
    platform_file <- normalizePath( platform_file )

    # load BIM file, used to track edits and write final output
    bim <- read_bim( input_data, flip )
    
    # all outputs should go in a subdirectory
    # this directory may exist already, if we're running for different thresholds
    if ( !dir.exists( dir_out ) )
        dir.create( dir_out )
    setwd( dir_out )

    ########################################################################

    phase <- 1
    iter <- 0

    # iteration 0 is the same for all p-values, put it in base directory to share automatically
    message( 'phase ', phase, ', iter ', iter )
    saige <- run_saige( platform_file, iid, platform_col, input_data )
    
    # if the directory exists, delete entirely! (to avoid overwriting existing files)
    if ( dir.exists( pval ) )
        unlink( pval, recursive = TRUE )
    dir.create( pval )
    setwd( pval )
    # from here on, p-value must be numeric to work
    pval <- as.numeric( pval )

    # create summary file and record 'remove' for each iteration
    iteration_summary <- NULL

    repeat {
        out <- identify_sig_snps( phase, iter, pval, bim, saige, iteration_summary )
        bim <- out$bim
        count <- out$count
        iteration_summary <- out$iteration_summary
        
        # stop loop when there are zero removals
        if ( count == 0 ) break
        
        # plink removes SNPs
        plink_apply_edits( bim, input_data, platform_two_id_file )
        
        # Increment for next iteration
        iter <- iter + 1
        
        # Run SAIGE
        message( 'phase ', phase, ', iter ', iter )
        saige <- run_saige( platform_file, iid, platform_col )
    }

    if ( flip ) {
        ########################################################################

        # start phase 2
        phase <- 2
        iter <- 0

        # one time, create file that assists flips with plink2
        platform_two_id_file <- make_platform_two_ids_file( platform_file, iid, fid, platform_col )

        # at the beginning of phase2 only, reclassify SNPs that are currently "remove" and flippable as "flip" (to be further removed in subsequent iterations)
        bim$category[ bim$category == 'remove' & bim$revcomp ] <- 'flip'

        # process flip and removal of SNPs, created new plink files
        plink_apply_edits( bim, input_data, platform_two_id_file )
        
        repeat { 
            message( 'phase ', phase, ', iter ', iter )
            saige <- run_saige( platform_file, iid, platform_col )

            out <- identify_sig_snps( phase, iter, pval, bim, saige, iteration_summary )
            bim <- out$bim
            count <- out$count
            iteration_summary <- out$iteration_summary
            
            # stop loop when there are zero removals
            if ( count == 0 ) break
            
            # remove snps with plink
            plink_apply_edits( bim, input_data, platform_two_id_file )

            # Increment for next iteration
            iter <- iter + 1
        }
        
        ########################################################################

        # start phase 3
        phase <- 3
        iter <- 0

        # at the beginning of phase 3, flip all loci we hadn't yet flipped otherwise
        # keep things we already determined to flip, but add to that flippable SNPs we kept but haven't flipped yet
        # also exclude SNPs with FWD p-values that were NA (unclear how to decide whether to flip those or not), we'll keep
        # (temporarily make a `bim` copy for this, because things don't look edited otherwise and p-values aren't updated correctly)
        bim2 <- bim
        bim2$category[ bim$revcomp & bim$category != 'remove' & !is.na( bim$pval_fwd ) ] <- 'flip'
        
        # created new plink files with remaining flippable SNPs flipped, to test both orientations
        plink_apply_edits( bim2, input_data, platform_two_id_file )
        
        # phase 3 is a single iteration, since nothing else is getting removed
        message( 'phase ', phase, ', iter ', iter )
        saige <- run_saige( platform_file, iid, platform_col )
        
        # perform the final classification!
        out <- identify_sig_snps( phase, iter, pval, bim, saige, iteration_summary )
        bim <- out$bim
        count <- out$count
        iteration_summary <- out$iteration_summary
        
        # done with this file!
        unlink( platform_two_id_file )
    }

    # write summary data, now that it is complete!
    write_tsv( iteration_summary, "iteration_summary.txt" ) 

    # produce final output of method, preds.txt.gz, which summarizes which loci are keep, remove, or flip.
    write_tsv( bim, gzfile( 'preds.txt.gz' ) )
}

# read BIM table, mark flippable SNPs, initialize output
read_bim <- function( input_data, flip = TRUE ) {
    bim <- read.table( paste0( input_data, '.bim' ), header = FALSE )
    colnames( bim ) <- c('chr', 'id', 'posg', 'pos', 'alt', 'ref')
    # subset to desired columns and reorder
    bim <- bim[ , c('chr', 'id', 'ref', 'alt') ]
    # a persistent issue is, at least in simulations, IDs are read as numeric, but writeLines requires character (and normally IDs are character)
    bim$id <- as.character( bim$id )
    # identify reverse complement cases, which are "flippable"
    if ( flip )
        bim$revcomp <- bim$ref == revcomp( bim$alt )
    # initialize category, which will get overwritten as we go
    bim$category <- 'keep'
    # other info for when edits occured
    # there's no phase 0, so 0,0 means unedited (to avoid NAs, which make some queries more difficult)
    bim$phase <- 0
    bim$iter <- 0
    # and the p-values, these are best initialized to NA
    bim$pval_fwd <- NA
    if ( flip )
        bim$pval_rev <- NA
    return( bim )
}

# global vars used here: seed, script_dir
# default `input_bfile` matches default `output_base` from `plink_apply_edits` (needed except first time here)
run_saige <- function( platform_file, iid = 'IID', platform_col = 'PLATFORM', input_bfile = 'edited' ) {
    output_prefix <- 'saige'
    output_file <- paste0( output_prefix, '_output.txt' )
    # first case can be compressed if it exists already
    output_file_gz <- paste0( output_file, '.gz' )
    
    # don't run if output already exists!
    # case we actually want to avoid repeating will be compressed
    if ( file.exists( output_file_gz ) ) {
        # don't run, set time as zero
        time <- 0
        # and tell other script this is the file to use
        output_file <- output_file_gz
    } else {
        # actually run
        # to pass to SAIGE
        seedopt <- if ( !is.null( seed ) ) paste0( '-s ', seed ) else ''
        
        # merged steps
        command <- paste0( 'Rscript ', script_dir, '/saige.R --bfile "', input_bfile, '" --platform "', platform_file, '" --iid "', iid, '" --platform_col "', platform_col, '" -o "', output_prefix, '" ', seedopt, ' &> ', output_prefix, '.log' )
        time <- system.time( ret <- system( command ) )[3]
        if ( ret != 0 ) stop( 'SAIGE failed with return value: ', ret, ' (see logs)' )
        
        # cleanup, not needed if run was successful
        unlink( paste0( output_prefix, c( '_output.txt.index', '.rda', '.varianceRatio.txt', '.log' ) ) )

        # compress first output only
        if ( input_bfile != 'edited' ) {
            ret <- system2( 'gzip', output_file )
            if ( ret != 0 ) stop( 'gzip failed with return value: ', ret )
            # and tell other script this is the file to use
            output_file <- output_file_gz
        } else
            # delete input bfile if "edited" (not original)
            unlink( paste0( input_bfile, '.', c( 'bed', 'bim', 'fam', 'log' ) ) )
    }
    
    # return file actually created (with full path, since we're going to move) and runtime for log
    return( list( file = normalizePath( output_file ), time = time ) )
}

identify_sig_snps <- function( phase, iter, pval, bim, saige, iteration_summary ) {
    # read SAIGE summary statistics
    saige_output_file <- saige$file
    data <- read.table( saige_output_file, header = TRUE )
    
    # we want to remember all p-values (not just significant ones)
    # most of these get overwritten at every iteration; removed loci keep their last p-value before removal
    # revcomp loci have p-values assigned to reverse orientation (not forward), but other loci have forward p-values overwritten
    # assume SAIGE is the subset (due to removals); all SAIGE SNPs have to be in BIM
    # this is a vector of indexes, length of SAIGE table
    indexes <- match( data$MarkerID, bim$id )
    bim2 <- bim[ indexes, ]
    if ( phase %in% 1:2 ) {
        # assignment depends on whether the SNP has been flipped already or not
        # (note current "remove" loci never have SAIGE data, because SAIGE is only run on loci that haven't been removed, so each locus is either "keep" or "flip")
        # this is a logical vector same length as indexes
        indexes2 <- bim2$category == 'keep'
        # when true, assign to FWD
        bim$pval_fwd[ indexes[ indexes2 ] ] <- data$p.value[ indexes2 ]
        # when false, assign to REV
        # (fails when we don't allow flips, since column `pval_rev` doesn't exist)
        if ( any( !indexes2 ) )
            bim$pval_rev[ indexes[ !indexes2 ] ] <- data$p.value[ !indexes2 ]
    } else {
        # phase 3
        # fill in only REV p-values for flippable loci that were previously keep, and which didn't have NA FWD p-values
        # all older p-values stay the same, they are more accurate that way
        indexes2 <- bim2$revcomp & bim2$category == 'keep' & !is.na( bim2$pval_fwd )
        ## # confirm that all those p-values are NA (this will not be true in iterated phase3)
        ## stopifnot( all( is.na( bim2$pval_rev[ indexes2 ] ) ) )
        bim$pval_rev[ indexes[ indexes2 ] ] <- data$p.value[ indexes2 ]
    }

    # this is the "reclassify" step
    if ( phase == 3 ) {
        # phase 3 doesn't have removals, and no need to set p-value thresholds
        
        # now that we have both p-values, make the final decision on which loci get flipped
        # this automatically focuses on flippable (revcomp) loci
        # need to additionally force that previous category was not "remove", otherwise removes get incorrectly rescued (into keep or flip) in this new comparison
        # this is a vector of indexes
        indexes <- which( !is.na( bim$pval_fwd ) & !is.na( bim$pval_rev ) & bim$category != 'remove' )
        bim2 <- bim[ indexes, ]
        category_new <- bim$category
        category_new[ indexes ] <- ifelse( bim2$pval_fwd < bim2$pval_rev, 'flip', 'keep' )
    } else {
        # recategorize using package function!
        category_new <- filter_category( bim, pval )
    }
    # mark changed cases as to which phase/iter the change happened in
    indexes <- bim$category != category_new
    bim$phase[ indexes ] <- phase
    bim$iter[ indexes ] <- iter
    # number of changes, for log and determine convergence
    count <- sum( indexes )
    # can overwrite now
    bim$category <- category_new
    
    # cleanup: don't need SAIGE file anymore, unless it's the first one
    if ( phase != 1 || iter != 0 )
        unlink( saige_output_file )

    # update log, even for stoping iterations because they took time (saige)
    iteration_summary <- rbind( iteration_summary, data.frame( phase = phase, iter = iter, edits = count, saige_runtime = saige$time ) )
    
    # return a few things
    list( bim = bim, count = count, iteration_summary = iteration_summary )
}

# make flip ID file (really platform 2 IDs)
make_platform_two_ids_file <- function( platform_file, iid, fid, platform_col ) {
    platform_two_id_file <- 'PLATFORM-TWO-IDs.txt'

    # read platform file
    data <- read.table( platform_file, header = TRUE )
    # keep individuals with PLATFORM==1 only (second platform, first one is 0)
    data <- data[ data[[ platform_col ]] == 1, ]
    # after this, `platform_col` isn't used by plink2, so it doesn't matter what it is
    # however, IID and FID must be what plink expects, lets rename columns
    colnames( data )[ colnames( data ) == iid ] <- 'IID'
    colnames( data )[ colnames( data ) == fid ] <- 'FID'
    # write to output
    write_tsv( data, platform_two_id_file )
    
    return( platform_two_id_file )
}

# apply current edits specified in `bim` to original `input_data`
plink_apply_edits <- function( bim, input_data, platform_two_id_file, output_base = 'edited' ) {
    # write these SNP IDs into files to perform edit in plink file
    exclude_file <- 'remove.txt'
    flip_file <- 'flip.txt'
    writeLines( bim$id[ bim$category == 'remove' ], exclude_file )
    # do extra work if we have flips (and avoid errors otherwise)
    have_flips <- any( bim$category == 'flip' )
    if ( have_flips )
        writeLines( bim$id[ bim$category == 'flip' ], flip_file )
    
    # run plink to remove loci and flip subset
    ret <- system2(
        'plink2',
        c(
            '--bfile', input_data,
            '--exclude', exclude_file,
            '--make-bed',
            '--out', output_base,
            '--silent',
            if ( have_flips )
                c(
                    '--flip', flip_file,
                    '--flip-subset', platform_two_id_file
                )
        )
    )
    if ( ret != 0 ) stop( 'plink2 failed with return value: ', ret, ' (see logs)' )
    
    # cleanup (ok if files were never created)
    unlink( c( exclude_file, flip_file ) )
}

write_tsv <- function( data, output_file )
    write.table( data, output_file, quote = FALSE, sep = "\t", row.names = FALSE )





### RUN ###

# terminal inputs
option_list = list(
    make_option("--bfile", type = "character", default = NA,
                help = "Input plink binary file without extensions (bed/bim/fam)", metavar = "character"),
    make_option("--platform", type = "character", default = NA,
                help = "Platform file that matches with the input data", metavar = "character"),
    make_option(c( "-d", "--dir_out"), type = "character", default = 'lmm-filter', 
                help = "Output directory (default lmm-filter)", metavar = "character"),
    make_option(c( "-s", "--seed"), type = "integer", default = NULL, 
                help = "Seed for random number generator", metavar = "integer"),
    make_option("--pval", type = "character", default = '1e-02',
                help = "P-value threshold for identifying significant snps, and exact name of output subdirectory (default 1e-02)", metavar = "character"),
    make_option("--noflip", action = "store_true", default = FALSE,
                help = "Run phase 1 only (removals only, no flips).  This is best for data where reverse-complement strand flips are not expected, such as whole-genome sequencing.  (Flips are often expected in genotyping array data.)"),
    make_option("--iid", type = "character", default = 'IID',
                help = "Name of individual ID column in platform file (default IID)", metavar = "character"),
    make_option("--fid", type = "character", default = 'FID',
                help = "Name of family ID column in platform file (default FID)", metavar = "character"),
    make_option("--platform_col", type = "character", default = 'PLATFORM',
                help = "Name of platform column (treated as binary trait) in platform file (default PLATFORM)", metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
# get values
input_data <- opt$bfile
platform_file <- opt$platform
pval <- opt$pval
dir_out <- opt$dir_out
seed <- opt$seed
flip <- !opt$noflip
iid <- opt$iid
fid <- opt$fid
platform_col <- opt$platform_col

# informative errors for the two required inputs
if ( is.na( input_data ) )
    stop( '`--bfile` is required!' )
if ( is.na( platform_file ) )
    stop( '`--platform` is required!' )

# annoying work to get current script location, to call other scripts within it as we navigate a directory structure elsewhere
initial_options <- commandArgs(trailingOnly = FALSE)
file_arg <- "--file="
script_path <- sub(file_arg, "", initial_options[grep(file_arg, initial_options)])
# Get the directory
script_dir <- normalizePath( dirname(script_path) )

# global vars: seed, script_dir
lmm_filter( input_data, platform_file, pval = pval, dir_out = dir_out, flip = flip, iid = iid, fid = fid, platform_col = platform_col )
