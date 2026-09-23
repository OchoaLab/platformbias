library(SAIGE)
library(optparse) 

# terminal inputs
option_list = list(
    make_option(c( "-f", "--file"), type = "character",
                help = "input file, after plink removed sig snps", metavar = "character"),
    make_option(c( "-p", "--platform"), type = "character",
                help = "Platform file that matches with the input data", metavar = "character"),
    make_option(c( "-o", "--out"), type = "character", default = NA, 
                help = "Output prefix", metavar = "character"),
    make_option(c( "-s", "--seed"), type = "integer", default = NULL, 
                help = "Seed for random number generator", metavar = "integer")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
# get values
plinkFile <- opt$file
phenoFile <- opt$platform
outputPrefix <- opt$out

set.seed( opt$seed )

fitNULLGLMM(
    plinkFile = plinkFile,
    phenoFile = phenoFile,
    phenoCol = 'PLATFORM',
    sampleIDColinphenoFile = 'IID',
    traitType = 'binary',
    outputPrefix = outputPrefix,
    IsOverwriteVarianceRatioFile = TRUE,
    LOCO = FALSE,
    minMAFforGRM = 0,
    maxMissingRateforGRM = 1
)
