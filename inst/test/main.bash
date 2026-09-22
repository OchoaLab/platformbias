### TEST DATA ###

# test data for local tests

# input files
dccDir=$dcc:/hpc/group/ochoalab/tt207/platform_bias/sim1/rep1/
scp $dccDir/biased.{bed,bim,fam} .
scp $dccDir/covar.txt .
# get this for now, but a good version shouldn't need it
scp $dccDir/id_p2.txt .

# make a smaller example that will run faster
# thin input data to have fewer SNPs, but same number of individuals (don't want to lose power)
wc -l biased.bim
# 500000 biased.bim
#plink2 --bfile biased --thin-count 10000 --make-bed --out test
# go even smaller!
plink2 --bfile biased --thin-count 1000 --make-bed --out test


### RUN ###

# then I manually ran the main script at this location

# generate data for replicates
PVAL=1e-02
# these need to exist already
input_data=test
# in practice this will be deducible from BIM file
flip_id_file="id_p2.txt"
# this appears to be misnamed, it's saige's PhenoFile (i.e. platform indicator)
covar_file="covar.txt"
# more hacks to run saige from pixi environment
# NOTE: also need to install other R packages in pixi env: optparse, genio, tidyverse, pak, then platformbias from github (for now)
shopt -s expand_aliases
alias Rscript="pixi run -m ~/bin/src/github/SAIGE/ Rscript"

# run lmm-filter!
time . /home/viiia/docs/ochoalab/packages/platformbias/inst/scripts/lmm-filter.bash

