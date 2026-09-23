library(optparse) 

# terminal inputs
option_list = list(
    make_option(c( "-p", "--platform"), type = "character",
                help = "Platform file that matches with the input data", metavar = "character"),
    make_option(c( "-o", "--out"), type = "character", default = NA, 
                help = "Output file", metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
# get values
phenoFile <- opt$platform
outputFile <- opt$out

# read covaraites file
data <- read.table( phenoFile, header = TRUE )
# keep individuals with PLATFORM==1 only (second platform, first one is 0)
data <- data[ data$PLATFORM == 1, ]
# write to output
write.table( data, outputFile, quote = FALSE, sep = "\t", row.names = FALSE )
