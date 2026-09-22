library(SAIGE)
library(optparse) 
# terminal inputs
option_list = list(
  make_option(c( "-f", "--file"), type = "character",
              help = "input file, after plink removed sig snps", metavar = "character"),
  make_option(c( "-o", "--out"), type = "character", default = NA, 
              help = "Output prefix", metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
# get values
plinkFile <- opt$f
outputPrefix <- opt$o

GMMATmodelFile = paste0(outputPrefix, '.rda') 
varianceRatioFile = paste0(outputPrefix, '.varianceRatio.txt') 
SAIGEOutputFile = paste0(outputPrefix, '_output.txt') 

SPAGMMATtest(bedFile=paste0(plinkFile, ".bed"),
             bimFile=paste0(plinkFile, ".bim"),
             famFile=paste0(plinkFile, ".fam"),
             AlleleOrder= 'alt-first',
             is_imputed_data=TRUE,
             #impute_method = opt$impute_method,
             GMMATmodelFile=GMMATmodelFile,
             varianceRatioFile=varianceRatioFile,
             SAIGEOutputFile=SAIGEOutputFile,
             is_output_moreDetails =TRUE,
             is_overwrite_output = TRUE,
             #SPAcutoff = opt$SPAcutoff, default 2
             is_Firth_beta = TRUE, # for binary traits
             #pCutoffforFirth = opt$pCutoffforFirth, # default 0.01
             LOCO = FALSE,
             min_MAF=0,
             min_MAC=0.5,
             max_missing = 1,
             dosage_zerod_cutoff = 0,
             dosage_zerod_MAC_cutoff = 0
)
