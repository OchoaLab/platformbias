library(SAIGE)
library(optparse) 

# terminal inputs
option_list = list(
  make_option(c( "-f", "--file"), type = "character",
              help = "input file, after plink removed sig snps", metavar = "character"),
  make_option(c( "-o", "--out"), type = "character", default = NA, 
              help = "Output prefix", metavar = "character"),
  make_option(c( "-s", "--seed"), type = "integer", default = NULL, 
              help = "Seed for random number generator", metavar = "integer")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
# get values
plinkFile <- opt$file
outputPrefix <- opt$out

set.seed( opt$seed )

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
