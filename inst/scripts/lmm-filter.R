library(platformbias) 
library(optparse)

# big pic TODOs
# - make flip (phase2) optional (in WGS data it shouldn't be needed)
# - add phase3 experiment

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

    # load BIM file, used to track edits and write final output
    bim <- read_bim( input_data )
    
    # all outputs should go in a subdirectory
    # this directory may exist already, if we're running for different thresholds
    if ( !dir.exists( dir_out ) )
        dir.create( dir_out )
    setwd( dir_out )

    ########################################################################

    phase <- 1
    iter <- 0

    # iteration 0 is the same for all p-values, put it in base directory to share automatically
    saige <- run_saige( phase, iter, platform_file, input_data )

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
        out <- identify_sig_snps( phase, iter, pval, bim, saige$file )
        bim <- out$bim
        count <- out$count
        
        # update log, even for stoping iterations because they took time (saige)
        iteration_summary <- rbind( iteration_summary, data.frame( phase = phase, iter = iter, remove = count, saige_runtime = saige$time ) )
        
        # stop loop when there are zero removals
        if ( count == 0 ) break
        
        # plink removes SNPs
        run_plink_remove( phase, iter, bim, input_data )
        
        # Increment for next iteration
        iter <- iter + 1
        
        # Run SAIGE
        saige <- run_saige( phase, iter, platform_file )
    }

    # more cleanup
    # delete last bed/bim/fam files (not used in phase2)
    delete_plink_bed( iter )

    ########################################################################

    # start phase 2
    phase <- 2
    iter <- 0

    # process flip and removal of SNPs, created new plink files
    bim <- run_plink_flip( bim, input_data, iter, platform_file )
    
    repeat { 
        saige <- run_saige( phase, iter, platform_file )

        # identify significant SNPs
        out <- identify_sig_snps( phase, iter, pval, bim, saige$file )
        bim <- out$bim
        count <- out$count

        # update log, even for stoping iterations because they took time (saige)
        iteration_summary <- rbind( iteration_summary, data.frame( phase = phase, iter = iter, remove = count, saige_runtime = saige$time ) )
        
        # stop loop when there are zero removals
        if ( count == 0 ) break
        
        # remove snps with plink
        run_plink_remove( phase, iter, bim, input_data )

        # Increment for next iteration
        iter <- iter + 1
    }
    
    # write summary data, now that it is complete!
    write.table( iteration_summary, "iteration_summary.txt", quote = FALSE, sep = "\t", row.names = FALSE )
    
    # delete BED files again.  These are not what we may want in the end because it's controls only.  In practice we want to process full data separately, following `preds.txt.gz`.
    delete_plink_bed( iter )
    
    # produce final output of method, preds.txt.gz, which summarizes which loci are keep, remove, or flip.
    write_preds( bim )
}

# read BIM table, mark flippable SNPs, initialize output
read_bim <- function( input_data ) {
    bim <- read.table( paste0( input_data, '.bim' ), header = FALSE )
    colnames( bim ) <- c('chr', 'id', 'posg', 'pos', 'alt', 'ref')
    # subset to desired columns and reorder
    bim <- bim[ , c('chr', 'id', 'ref', 'alt') ]
    # a persistent issue is, at least in simulations, IDs are read as numeric, but writeLines requires character (and normally IDs are character)
    bim$id <- as.character( bim$id )
    # identify reverse complement cases, which are "flippable"
    bim$revcomp <- bim$ref == revcomp( bim$alt )
    # initialize category, which will get overwritten as we go
    bim$category <- 'keep'
    # other info for when edits occured
    # there's no phase 0, so 0,0 means unedited (to avoid NAs, which make some queries more difficult)
    bim$phase <- 0
    bim$iter <- 0
    return( bim )
}

# global vars used here: seed, script_dir
run_saige <- function( phase, iter, platform_file, input_bfile = iter ) {
    # every iteration starts with a SAIGE run, might as well put this here
    message( 'phase ', phase, ', iter ', iter )
    
    output_prefix <- paste0( 'saige_', phase, '_', iter )
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
        command <- paste0( 'Rscript ', script_dir, '/saige.R -f "', input_bfile, '" -p "', platform_file, '" -o "', output_prefix, '" ', seedopt, ' > ', output_prefix, '.log' )
        time <- system.time( system( command ) )[3]
        
        # cleanup, not needed if run was successful
        unlink( paste0( output_prefix, c( '_output.txt.index', '.rda', '.varianceRatio.txt', '.log' ) ) )

        # compress first output only
        if ( phase == 1 && iter == 0 ) {
            system2( 'gzip', output_file )
            # and tell other script this is the file to use
            output_file <- output_file_gz
        }
    }
    
    # return file actually created (with full path, since we're going to move) and runtime for log
    return( list( file = normalizePath( output_file ), time = time ) )
}

identify_sig_snps <- function( phase, iter, pval, bim, saige_output_file ) {
    # read SAIGE summary statistics
    data <- read.table( saige_output_file, header = TRUE )
    # get significant SNPs, which are the new removals
    sig_snps <- as.character( data$MarkerID[ data$p.value < pval ] )
    
    # mark new removals in BIM file, with detailed info
    indexes <- bim$id %in% sig_snps
    bim$category[ indexes ] <- 'remove'
    bim$phase[ indexes ] <- phase
    bim$iter[ indexes ] <- iter
    
    # cleanup: don't need SAIGE file anymore, unless it's the first one
    if ( phase != 1 || iter != 0 )
        unlink( saige_output_file )
    
    # return a few things
    list( bim = bim, count = length( sig_snps ) )
}

delete_plink_bed <- function( input_bfile )
    unlink( paste0( input_bfile, '.', c( 'bed', 'bim', 'fam', 'log' ) ) )

run_plink_remove <- function( phase, iter, bim, input_data ) {
    # special behavior for very first case only
    first <- phase == 1 && iter == 0
    input_bfile <- if ( first ) input_data else iter
    
    # make file with SNP IDs to remove using bim table
    # for maximum efficiency, fish out removals from this phase/iteration only
    exclude_file <- paste0("remove_phase", phase, "_", iter, ".txt")
    ids_rm <- bim$id[ bim$phase == phase & bim$iter == iter ]
    writeLines( ids_rm, exclude_file )
    
    # remove SNPs with plink2
    system2(
        'plink2',
        c(
            '--bfile', input_bfile,
            '--exclude', exclude_file,
            '--make-bed',
            '--out', iter + 1,
            '--silent'
        )
    )

    # cleanup: we don't need input anymore unless it's the original file!
    if ( !first )
        delete_plink_bed( input_bfile )
    # we're also done with this file
    unlink( exclude_file )
}

run_plink_flip <- function( bim, input_data, output_prefix, platform_file ) {
    # first, reclassify SNPs that are currently "remove" and flippable as "flip" (to be further removed in subsequent iterations)
    bim$category[ bim$category == 'remove' & bim$revcomp ] <- 'flip'
    # write these SNP IDs into files to perform edit in plink file
    exclude_file <- 'phase2_init_remove.txt'
    flip_file <- 'phase2_init_flip.txt'
    writeLines( bim$id[ bim$category == 'remove' ], exclude_file )
    writeLines( bim$id[ bim$category == 'flip' ], flip_file )
    
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
    unlink( c( platform_two_id_file, exclude_file, flip_file ) )
    # NOTE: this uses the original `input_data` always, never delete it! (unlike run_plink_remove)
    
    # this was edited, return!
    return( bim )
}

write_preds <- function( bim ) {
    # HACK TMP
    # remove new columns, to simplify comparisons to old outputs
    #bim <- bim[ , 1:6 ]
    
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
