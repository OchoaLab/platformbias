library(platformbias) 
#library(SAIGE)
library(optparse)

# big pic TODOs
# - make flip (phase2) optional (in WGS data it shouldn't be needed)
# - add phase3 experiment

## rlang::global_entrace()
## options(rlang_backtrace_on_error = "full")

# main function/loop
# initial p-value is string because we want folder name to be exactly this ("1e-02" instead of "0.01"), but after dir is made we can turn to numeric
lmm_filter <- function( input_data, platform_file, pval = '1e-02', dir_out = 'lmm-filter', script_dir = '.', seed = NULL ) {
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
    
    # to pass to SAIGE
    seedopt <- if ( !is.null( seed ) ) paste0( '-s ', seed ) else ''
    
    # force inputs to have absolute paths
    # for `input_data`, need to add extension (so file exists, otherwise this doesn't work), then take it back off
    input_data <- sub( '.bed$', '', normalizePath( paste0( input_data, '.bed' ) ) )
    platform_file <- normalizePath( platform_file )

    # all outputs should go in a subdirectory
    dir.create( dir_out )
    setwd( dir_out )
    
    ITER <- 0

    # iteration 0 is the same for all p-values, put it in base directory to share automatically
    run_saige( input_data, platform_file, paste0( 'saige_', ITER ), script_dir, seedopt )
    # compress the first output only
    system2( 'gzip', 'saige_0_output.txt' )

    # if the directory exists, delete entirely! (to avoid overwriting existing files)
    if ( dir.exists( pval ) )
        unlink( pval, recursive = TRUE )
    dir.create( pval )
    setwd( pval )
    # from here on, p-value must be numeric to work
    pval <- as.numeric( pval )

    ########################################################################

    message( "Starting Phase 1" )

    # create summary file and record 'remove' for each iteration
    iteration_summary <- NULL

    repeat {
        # Rscript identifies significant SNPs; writes remove_phase1_$ITER.txt
        out <- run_identify_sig_snps( "phase1", ITER, pval )
        remove_file <- out$remove_file
        remove_count <- out$remove_count
        
        # stop loop when there are zero removals
        if ( remove_count == 0 ) {
            message( 'No removals at ITER=', ITER, '. Exiting loop.' )
            break
        }
        # else update log
        iteration_summary <- rbind( iteration_summary, data.frame( phase = 1, ITER = ITER, remove = remove_count ) )

        # plink removes SNPs
        # get previous iteration file as the input file for plink processing
        input_bfile <- if ( ITER == 0 ) input_data else ITER - 1
        run_plink_remove( input_bfile, remove_file, ITER, input_data )
        
        # Run SAIGE
        run_saige( ITER, platform_file, paste0( 'saige_phase1_', ITER ), script_dir, seedopt )

        # Increment for next iteration
        ITER <- ITER + 1
    }

    # more cleanup
    # delete last bed/bim/fam files (not used in phase2)
    unlink( paste0( ITER - 1, '.', c( 'bed', 'bim', 'fam', 'log' ) ) )

    ########################################################################

    # start phase 2
    message( "Starting Phase 2" )

    # identify snps to flip and remove
    run_phase2_flip_snps( input_data )

    # process flip and removal of SNPs, created new plink files
    run_plink_flip( input_data, 'phase2_init_remove.txt', 0, 'phase2_init_flip.txt', platform_file )
    
    # output file written as 'control_0' in phase2 folder, next step is to run saige
    ITER <- 0

    repeat { 
        run_saige( ITER, platform_file, paste0( 'saige_phase2_', ITER ), script_dir, seedopt )

        # identify significant SNPs
        out <- run_identify_sig_snps( "phase2", ITER, pval )
        remove_file <- out$remove_file
        remove_count <- out$remove_count

        # stop loop when there are zero removals
        if ( remove_count == 0 ) {
            message( 'No removals at ITER=', ITER, '. Exiting loop.' )
            # write summary data, now that it is complete!
            write.table( iteration_summary, "iteration_summary.txt", quote = FALSE, sep = "\t", row.names = FALSE )
            break
        }
        # else update log
        iteration_summary <- rbind( iteration_summary, data.frame( phase = 2, ITER = ITER, remove = remove_count ) )
        
        # remove snps with plink
        run_plink_remove( ITER, remove_file, ITER + 1, input_data )

        # Increment for next iteration
        ITER <- ITER + 1
    }

    # delete BED files again.  These are not what we may want in the end because it's controls only.  In practice we want to process full data separately, following `preds.txt.gz`.
    unlink( paste0( ITER, '.', c( 'bed', 'bim', 'fam', 'log' ) ) )
    
    # produce final output of method, preds.txt.gz, which summarizes which loci are keep, remove, or flip.
    lmm_classify( input_data )
    
    # if we're here, presumably everything was good and we don't need those bulky logs anymore... (they're only for troubleshooting)
    #rm saige_phase?_*_step?.log
    # for now let's keep the full saige files, in case we wanted their p-values...
    # (otherwise consider deleting them)
    system( 'gzip saige_phase?_*_output.txt' )
}

run_saige <- function( input_bfile, platform_file, output_prefix, script_dir, seedopt ) {
    # don't run if output already exists!
    output_file <- paste0( output_prefix, '_output.txt' )
    # case we actually want to avoid repeating will be compressed
    if ( file.exists( paste0( output_file, '.gz' ) ) )
        return()

    message( 'SAIGE step 1' )
    system( paste0( 'time Rscript ', script_dir, '/saige_step1_nocovar.R -f "', input_bfile, '" -p "', platform_file, '" -o "', output_prefix, '" ', seedopt, ' > ', output_prefix, '_step1.log' ) )
    
    ## capture.output(
        ## fitNULLGLMM(
        ##     plinkFile = input_bfile,
        ##     phenoFile = platform_file,
        ##     phenoCol = 'PLATFORM',
        ##     sampleIDColinphenoFile = 'IID',
        ##     traitType = 'binary',
        ##     outputPrefix = output_prefix,
        ##     IsOverwriteVarianceRatioFile = TRUE,
        ##     LOCO = FALSE,
        ##     minMAFforGRM = 0,
        ##     maxMissingRateforGRM = 1
        ## ) #,
    ##     file = paste0( output_prefix, '_step1.log' ),
    ##     type = c("output", "message")
    ## )

    message( 'SAIGE step 2' )
    system( paste0( 'time Rscript ', script_dir, '/saige_step2.R -f "', input_bfile, '" -o "', output_prefix, '" ', seedopt, ' > ', output_prefix, '_step2.log' ) )
    ## capture.output(
        ## SPAGMMATtest(
        ##     bedFile = paste0( input_bfile, ".bed" ),
        ##     bimFile = paste0( input_bfile, ".bim" ),
        ##     famFile = paste0( input_bfile, ".fam" ),
        ##     AlleleOrder = 'alt-first',
        ##     is_imputed_data = TRUE,
        ##     GMMATmodelFile = paste0( output_prefix, '.rda' ),
        ##     varianceRatioFile = paste0( output_prefix, '.varianceRatio.txt' ),
        ##     SAIGEOutputFile = output_file,
        ##     is_output_moreDetails = TRUE,
        ##     is_overwrite_output = TRUE,
        ##     is_Firth_beta = TRUE,
        ##     LOCO = FALSE,
        ##     min_MAF = 0,
        ##     min_MAC = 0.5,
        ##     max_missing = 1,
        ##     dosage_zerod_cutoff = 0,
        ##     dosage_zerod_MAC_cutoff = 0
        ## )
## ,
##         file = paste0( output_prefix, '_step2.log' ),
##         type = c("output", "message")
##     )
    
    # cleanup, not needed if run was successful
    unlink( paste0( output_prefix, c( '_output.txt.index', '.rda', '.varianceRatio.txt', '_step1.log', '_step2.log' ) ) )
}

# TODO: this could all be handled virtually (not writing files)
run_identify_sig_snps <- function( phase, iter_num, pval ) {
    # organize file names across different phases and iterations
    # main file
    main_file <- paste0("remove_", phase, ".txt")

    # input file
    saige_output <- if (phase == "phase1" && iter_num == 0) {
                        "../saige_0_output.txt.gz" # only this one is already compressed
                    } else if (phase == "phase1") {
                        # read saige output file from previous iteration
                        paste0("saige_", phase, '_', iter_num -1, "_output.txt")
                    } else if (phase == "phase2") {
                        paste0("saige_", phase, '_', iter_num, "_output.txt")
                    }

    # output file for current iteration
    current_file <- paste0("remove_", phase, "_", iter_num, ".txt")

    # main script
    data <- read.table( saige_output, header = TRUE )
    sig_snps <- as.character( data$MarkerID[ data$p.value < pval ] )
    remove_count <- length( sig_snps )

    # writes separate file just for this iteration
    # this can be an empty file (loop outside decides to stop when that happens)
    writeLines( sig_snps, current_file )

    if (iter_num == 0) {
        # this creates main file
        writeLines( sig_snps, main_file )
    } else if ( remove_count > 0 ) {
        # Load existing main SNPs and update it
        snps_main <- readLines( main_file )
        combined_snps <- union( snps_main, sig_snps )
        writeLines( combined_snps, main_file )
    }

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
	unlink( paste0( input_bfile, '.', c( 'bed', 'bim', 'fam', 'log' ) ) )
    # TODO: do we still need ${input_bfile}_output.txt ???
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

    # assumes we're inside the phase2 dir

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

#set.seed( opt$seed )

# annoying work to get current script location, to call other scripts within it as we navigate a directory structure elsewhere
initial_options <- commandArgs(trailingOnly = FALSE)
file_arg <- "--file="
script_path <- sub(file_arg, "", initial_options[grep(file_arg, initial_options)])
# Get the directory
script_dir <- normalizePath( dirname(script_path) )

lmm_filter( input_data, platform_file, pval = pval, dir_out = dir_out, script_dir = script_dir, seed = opt$seed )
