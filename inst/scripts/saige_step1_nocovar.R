library(SAIGE)
library(optparse) 

# terminal inputs
option_list = list(
  make_option(c( "-f", "--file"), type = "character",
              help = "input file, after plink removed sig snps", metavar = "character"),
  make_option(c( "-c", "--covar"), type = "character",
              help = "covar file that matches with the input data", metavar = "character"),
  make_option(c( "-o", "--out"), type = "character", default = NA, 
              help = "Output prefix", metavar = "character"),
  make_option(c( "-s", "--seed"), type = "integer", default = NULL, 
              help = "Seed for random number generator", metavar = "integer")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
# get values
plinkFile <- opt$f
phenoFile <- opt$c
outputPrefix <- opt$o

set.seed( opt$seed )

phenoCol= 'pheno'
sampleIDColinphenoFile='iid' 
traitType='binary'        
IsOverwriteVarianceRatioFile=TRUE

fitNULLGLMM(plinkFile = plinkFile,
            phenoFile = phenoFile,
            phenoCol = phenoCol,
            sampleIDColinphenoFile = sampleIDColinphenoFile,
            traitType = traitType,
            outputPrefix = outputPrefix,
            #covarColList = c("famid", "sex"),
            #qCovarCol = c("famid", "sex"),
            IsOverwriteVarianceRatioFile = TRUE,
            LOCO = FALSE,
            #sexCol = "sex",
            minMAFforGRM = 0,
            maxMissingRateforGRM = 1
)
