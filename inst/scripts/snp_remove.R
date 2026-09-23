library(optparse) 

# terminal inputs
option_list = list(
    make_option(c( "-l", "--loop"), type = "character", default = NA, 
                help = "phase 1 or 2", metavar = "character"),
    make_option(c( "-i", "--iteration"), type = "integer", default = NA, 
                help = "numeric number for iteration of saige control vs control", metavar = "integer"),
    make_option(c( "-p", "--pval"), type = "numeric", default = NA, 
                help = "pvalue threshold for identifying significant snps", metavar = "numeric")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
# get values
phase <- opt$l
iter_num <- opt$i
pval <- opt$p

###########################################
# organize file names across different phases and iterations
# main file
main_file = paste0("remove_", phase, ".txt")

# input file
saige_output = if (phase == "phase1" && iter_num == 0) {
                   "../saige_0_output.txt.gz" # only this one is already compressed
               } else if (phase == "phase1") {
                   # read saige output file from previous iteration
                   paste0("saige_", phase, '_', iter_num -1, "_output.txt")
               } else if (phase == "phase2") {
                   paste0("saige_", phase, '_', iter_num, "_output.txt")
               }

# output file for current iteration
current_file = paste0("remove_", phase, "_", iter_num, ".txt")

###########################################
# main script
data = read.table(saige_output, header = TRUE)
sig_snps = as.character( data$MarkerID[ data$p.value < pval ] )

# write output file depending on phase & iteration, also update main list

# writes separate file just for this iteration as well
# this can be an empty file (loop outside decides to stop when that happens)
writeLines(sig_snps, current_file)

if (iter_num == 0) {
    # this creates main file
    writeLines(sig_snps, main_file)
} else if (length(sig_snps) > 0) {
    # Load existing main SNPs and update list
    snps_main <- readLines(main_file)
    combined_snps <- union( snps_main, sig_snps )
    writeLines(combined_snps, main_file)
}
