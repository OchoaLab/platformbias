suppressWarnings(library(SAIGE))
library(optparse) 

# terminal inputs
option_list = list(
    make_option("--bfile", type = "character",
                help = "Input plink binary file without extensions (bed/bim/fam)", metavar = "character"),
    make_option("--platform", type = "character",
                help = "Platform file that matches with the input data", metavar = "character"),
    make_option("--iid", type = "character", default = 'IID',
                help = "Name of individual ID column in platform file", metavar = "character"),
    make_option("--platform_col", type = "character", default = 'PLATFORM',
                help = "Name of platform column (treated as binary trait) in platform file", metavar = "character"),
    make_option(c( "-o", "--out"), type = "character", default = NA, 
                help = "Output prefix", metavar = "character"),
    make_option(c( "-s", "--seed"), type = "integer", default = NULL, 
                help = "Seed for random number generator", metavar = "integer")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
# get values
plinkFile <- opt$bfile
phenoFile <- opt$platform
outputPrefix <- opt$out
iid <- opt$iid
platform_col <- opt$platform_col

set.seed( opt$seed )

message( 'SAIGE step 1' )
fitNULLGLMM(
    plinkFile = plinkFile,
    phenoFile = phenoFile,
    phenoCol = platform_col,
    sampleIDColinphenoFile = iid,
    traitType = 'binary',
    outputPrefix = outputPrefix,
    IsOverwriteVarianceRatioFile = TRUE,
    LOCO = FALSE,
    minMAFforGRM = 0,
    maxMissingRateforGRM = 1
)

message( 'SAIGE step 2' )
SPAGMMATtest(
    bedFile = paste0( plinkFile, ".bed" ),
    bimFile = paste0( plinkFile, ".bim" ),
    famFile = paste0( plinkFile, ".fam" ),
    AlleleOrder = 'alt-first',
    is_imputed_data = TRUE,
    GMMATmodelFile = paste0( outputPrefix, '.rda' ),
    varianceRatioFile = paste0( outputPrefix, '.varianceRatio.txt' ),
    SAIGEOutputFile = paste0( outputPrefix, '_output.txt' ),
    is_output_moreDetails = TRUE,
    is_overwrite_output = TRUE,
    is_Firth_beta = TRUE,
    LOCO = FALSE,
    min_MAF = 0,
    min_MAC = 0.5,
    max_missing = 1,
    dosage_zerod_cutoff = 0,
    dosage_zerod_MAC_cutoff = 0
)
