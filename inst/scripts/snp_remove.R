library(readr)
library(dplyr)
library(genio)
library(optparse) 

# terminal inputs
option_list = list(
  make_option(c( "-l", "--loop"), type = "character", default = '1', 
              help = "phase 1 or 2", metavar = "character"),
  make_option(c( "-i", "--iteration"), type = "character", default = '1', 
              help = "numeric number for iteration of saige control vs control", metavar = "character"),
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
# master file
file_name_master = function(){
  paste0("remove_", phase, ".txt")
}

# input file
file_name_input <- function(iter_num, phase) {
  iter_num <- as.numeric(iter_num)
  if (phase == "phase1" && iter_num == 0) {
    return("../saige_0_output.txt.gz") # only this one is already compressed
  } else if (phase == "phase1") {
    # read saige output file from previous iteration
    return(paste0("saige_", phase, '_', iter_num -1, "_output.txt"))
  } else if (phase == "phase2") {
    return(paste0("saige_", phase, '_', iter_num, "_output.txt"))
  }
}
# output file for current iteration
file_name_output = function(iter_num){
  paste0("remove_", phase, "_", iter_num, ".txt")
}

master_file = file_name_master()
current_file = file_name_output(iter_num)
saige_output = file_name_input(iter_num, phase)
###########################################
# main script
sig_snps = read.table(saige_output, header = TRUE) %>% filter(p.value < pval) %>% pull(MarkerID)
print("remove SNPs:")
print(length(sig_snps))

# write output file depending on phase & iteration, also update master list

if (iter_num == 0) {
  write_lines(sig_snps, master_file)
  write_lines(sig_snps, current_file)
  
} else {
  if (length(sig_snps) == 0) {
    message("No SNPs pass the p-value threshold. Stopping.")
    stop("No significant SNPs found.")
  }
  
  # Load existing master SNPs and update list
  snps_master <- readLines(master_file)
  combined_snps <- unique(c(snps_master, sig_snps))
  write_lines(combined_snps, master_file)
  # writes separate file just for this iteration as well
  write_lines(sig_snps, current_file)
  
  print("remove SNPs accumulated count:")
  print(length(combined_snps))
}

