library(platformbias) 
#library(SAIGE)
library(optparse)

# big pic TODOs
# - make flip (phase2) optional (in WGS data it shouldn't be needed)
# - add phase3 experiment
# minor tasks
# - streamline tracking of removals and flips without writing files

## rlang::global_entrace()
## options(rlang_backtrace_on_error = "full")

# main function/loop
# initial p-value is string because we want folder name to be exactly this ("1e-02" instead of "0.01"), but after dir is made we can turn to numeric
lmm_filter <- function( input_data, platform_file, pval = '1e-02', dir_out = 'lmm-filter' ) {
    # validate inputs
    if ( missing( input_data ) )
        stop( '`input_data` is required!' )
    if ( missing( platform_file ) )
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

    # all outputs should go in a subdirectory
    # this directory may exist already, if we're running for different thresholds
    if ( !dir.exists( dir_out ) )
        dir.create( dir_out )
    setwd( dir_out )

    ########################################################################

    phase <- 1
    iter <- 0
    message( 'phase ', phase, ', iter ', iter )

    # iteration 0 is the same for all p-values, put it in base directory to share automatically
    runtime <- run_saige( phase, iter, platform_file, input_data )
    # compress the first output only, if there's need
    if ( file.exists( 'saige_phase1_0_output.txt' ) )
        system2( 'gzip', 'saige_phase1_0_output.txt' )

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
        out <- run_identify_sig_snps( phase, iter, pval )
        remove_file <- out$remove_file
        remove_count <- out$remove_count
        
        # update log, even for stoping iterations because they took time (saige)
        iteration_summary <- rbind( iteration_summary, data.frame( phase = phase, iter = iter, remove = remove_count, saige_runtime = runtime ) )
        
        # stop loop when there are zero removals
        if ( remove_count == 0 )
            break
        
        # plink removes SNPs
        # get previous iteration file as the input file for plink processing
        input_bfile <- if ( iter == 0 ) input_data else iter
        run_plink_remove( input_bfile, remove_file, iter + 1, input_data )
        
        # Increment for next iteration
        iter <- iter + 1
        
        message( 'phase ', phase, ', iter ', iter )
        
        # Run SAIGE
        runtime <- run_saige( phase, iter, platform_file )
    }

    # more cleanup
    # delete last bed/bim/fam files (not used in phase2)
    delete_plink_bed( iter )

    ########################################################################

    # start phase 2
    phase <- 2
    iter <- 0

    # identify snps to flip and remove
    run_phase2_flip_snps( input_data )

    # process flip and removal of SNPs, created new plink files
    run_plink_flip( input_data, 'phase2_init_remove.txt', iter, 'phase2_init_flip.txt', platform_file )
    
    repeat { 
        message( 'phase ', phase, ', iter ', iter )
        runtime <- run_saige( phase, iter, platform_file )

        # identify significant SNPs
        out <- run_identify_sig_snps( phase, iter, pval )
        remove_file <- out$remove_file
        remove_count <- out$remove_count

        # update log, even for stoping iterations because they took time (saige)
        iteration_summary <- rbind( iteration_summary, data.frame( phase = phase, iter = iter, remove = remove_count, saige_runtime = runtime ) )
        
        # stop loop when there are zero removals
        if ( remove_count == 0 )
            break
        
        # remove snps with plink
        run_plink_remove( iter, remove_file, iter + 1, input_data )

        # Increment for next iteration
        iter <- iter + 1
    }
    
    # write summary data, now that it is complete!
    write.table( iteration_summary, "iteration_summary.txt", quote = FALSE, sep = "\t", row.names = FALSE )
    
    # delete BED files again.  These are not what we may want in the end because it's controls only.  In practice we want to process full data separately, following `preds.txt.gz`.
    delete_plink_bed( iter )
    
    # produce final output of method, preds.txt.gz, which summarizes which loci are keep, remove, or flip.
    lmm_classify( input_data )
}

# global vars used here: seed, script_dir
run_saige <- function( phase, iter, platform_file, input_bfile = iter ) {
    output_prefix <- paste0( 'saige_phase', phase, '_', iter )
    # don't run if output already exists!
    output_file <- paste0( output_prefix, '_output.txt' )
    # case we actually want to avoid repeating will be compressed
    if ( file.exists( paste0( output_file, '.gz' ) ) )
        return( 0 )
    
    # to pass to SAIGE
    seedopt <- if ( !is.null( seed ) ) paste0( '-s ', seed ) else ''
    
    # merged steps
    command <- paste0( 'Rscript ', script_dir, '/saige.R -f "', input_bfile, '" -p "', platform_file, '" -o "', output_prefix, '" ', seedopt, ' > ', output_prefix, '.log' )
    time <- system.time( system( command ) )[3]
    
    # cleanup, not needed if run was successful
    unlink( paste0( output_prefix, c( '_output.txt.index', '.rda', '.varianceRatio.txt', '.log' ) ) )

    # add time to log
    return( time )
}

# TODO: this could all be handled virtually (not writing files)
run_identify_sig_snps <- function( phase, iter, pval ) {
    # organize file names across different phases and iterations
    # main file
    main_file <- paste0("remove_phase", phase, ".txt")

    # input file
    saige_output_file <- paste0("saige_phase", phase, '_', iter, "_output.txt")
    # only this one is compressed and one level down, because it's shared
    shared <- phase == 1 && iter == 0
    if ( shared )
        saige_output_file <- paste0( '../', saige_output_file, '.gz' )
    
    # output file for current iteration
    # plink2 will use this to remove SNPs from BED file, then it gets deleted
    current_file <- paste0("remove_phase", phase, "_", iter, ".txt")
    
    # read SAIGE summary statistics
    data <- read.table( saige_output_file, header = TRUE )
    sig_snps <- as.character( data$MarkerID[ data$p.value < pval ] )
    remove_count <- length( sig_snps )

    # writes separate file just for this iteration
    # don't bother writing an empty file
    if ( remove_count > 0 )
        writeLines( sig_snps, current_file )
    
    if ( iter == 0 ) {
        # this creates main file
        writeLines( sig_snps, main_file )
    } else if ( remove_count > 0 ) {
        # Load existing main SNPs and update it
        snps_main <- readLines( main_file )
        combined_snps <- union( snps_main, sig_snps )
        writeLines( combined_snps, main_file )
    }

    # cleanup: don't need SAIGE file anymore, unless it's the shared one
    if ( !shared )
        unlink( saige_output_file )
    
    # return a few things
    list( remove_file = current_file, remove_count = remove_count )
}

run_phase2_flip_snps <- function( input_data ) {
    # identify SNPs that can be flipped
    remove_snps <- readLines( "remove_phase1.txt" )

    bim <- read.table( paste0( input_data, '.bim' ), header = FALSE )
    colnames( bim ) <- c('chr', 'id', 'posg', 'pos', 'alt', 'ref')
    ids_revcomp <- bim$id[ bim$ref == revcomp( bim$alt ) ]

    sig_snps_flip <- intersect( remove_snps, ids_revcomp )
    sig_snps_remove <- setdiff( remove_snps, ids_revcomp )

    # write output and rerun saige
    writeLines( sig_snps_flip, "phase2_init_flip.txt" )
    writeLines( sig_snps_remove, "phase2_init_remove.txt" )
}

delete_plink_bed <- function( input_bfile )
    unlink( paste0( input_bfile, '.', c( 'bed', 'bim', 'fam', 'log' ) ) )

run_plink_remove <- function( input_bfile, exclude_file, output_prefix, input_data ) {
    # remove SNPs with plink2
    system2(
        'plink2',
        c(
            '--bfile', input_bfile,
            '--exclude', exclude_file,
            '--make-bed',
            '--out', output_prefix,
            '--silent'
        )
    )

    # cleanup: we don't need input anymore unless it's the original file!
    if ( input_bfile != input_data )
        delete_plink_bed( input_bfile )
    # we're also done with this file
    unlink( exclude_file )
}

run_plink_flip <- function( input_data, exclude_file, output_prefix, flip_file, platform_file ) {
    # make flip ID file (really platform 2 IDs)
    platform_two_id_file <- 'PLATFORM-TWO-IDs.txt'

    # read platform file
    data <- read.table( platform_file, header = TRUE )
    # keep individuals with PLATFORM==1 only (second platform, first one is 0)
    data <- data[ data$PLATFORM == 1, ]
    # write to output
    write.table( data, platform_two_id_file, quote = FALSE, sep = "\t", row.names = FALSE )

    # run plink to remove loci and flip subset
    system2(
        'plink2',
        c(
            '--bfile', input_data,
            '--exclude', exclude_file,
            '--make-bed',
            '--out', output_prefix,
            '--silent',
            '--flip', flip_file,
            '--flip-subset', platform_two_id_file
        )
    )

    # cleanup
    unlink( platform_two_id_file )
    # NOTE: this uses the original input data always, never delete it! (unlike run_plink_remove)
}

lmm_classify <- function( input_file ) {
    # takes predictions (disorganized remove and flip lists) for a certain run, and produces a nice classification vector and other useful info, similar to preds.txt.gz for the simpler methods

    # read BIM file to add data to
    bim <- read.table( paste0( input_file, '.bim' ), header = FALSE )
    colnames( bim ) <- c('chr', 'id', 'posg', 'pos', 'alt', 'ref')
    # subset to desired columns and reorder
    bim <- bim[ , c('chr', 'id', 'ref', 'alt') ]

    # identify reverse complement cases, which are "flippable" (not strictly necessary, though these files I've made usually have them and allow for internal checks)
    bim$revcomp <- bim$ref == revcomp( bim$alt )

    ### LMM DATA ###

    # first read final "phase2" edits, including permanently removed and temporarily flipped locus
    remove <- readLines( "phase2_init_remove.txt" )
    flip <- readLines( "phase2_init_flip.txt" )

    # check that remove and flip are unique and disjoint
    stopifnot( length( remove ) == length( unique( remove ) ) )
    stopifnot( length( flip ) == length( unique( flip ) ) )
    stopifnot( length( intersect( remove, flip ) ) == 0 )

    # phase 2 removals 
    phase2_snps <- readLines( "remove_phase2.txt" )
    total_remove <- union( remove, phase2_snps )

    # final text classification
    # this is most cases
    bim$category <- 'keep'
    # "flip" will have cases overwritten into "remove" if they were subsequently removed (in phase 3), happens in next step
    bim$category[ bim$id %in% flip ] <- 'flip'
    bim$category[ bim$id %in% total_remove ] <- 'remove'

    # confirm that all loci that were flipped are actually flippable
    stopifnot( all( bim$revcomp[ bim$category == 'flip' ] ) )

    # save updated/extended `bim`, the key calculation!
    write.table( bim, 'preds.txt.gz', quote = FALSE, sep = "\t", row.names = FALSE )
}






### RUN ###

# terminal inputs
option_list = list(
    make_option(c( "-f", "--file"), type = "character",
                help = "input plink binary file without extensions", metavar = "character"),
    make_option("--platform", type = "character",
                help = "Platform file that matches with the input data", metavar = "character"),
    make_option(c( "-d", "--dir_out"), type = "character", default = 'lmm-filter', 
                help = "Output directory", metavar = "character"),
    make_option(c( "-s", "--seed"), type = "integer", default = NULL, 
                help = "Seed for random number generator", metavar = "integer"),
    make_option("--pval", type = "character", default = '1e-02',
                help = "pvalue threshold for identifying significant snps", metavar = "numeric")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
# get values
input_data <- opt$file
platform_file <- opt$platform
pval <- opt$pval
dir_out <- opt$dir_out
seed <- opt$seed

# annoying work to get current script location, to call other scripts within it as we navigate a directory structure elsewhere
initial_options <- commandArgs(trailingOnly = FALSE)
file_arg <- "--file="
script_path <- sub(file_arg, "", initial_options[grep(file_arg, initial_options)])
# Get the directory
script_dir <- normalizePath( dirname(script_path) )

# global vars: seed, script_dir
lmm_filter( input_data, platform_file, pval = pval, dir_out = dir_out )
