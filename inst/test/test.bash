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

# set a seed so saige (which is random) gives reproducible results in tests
# saige is still sort of random, but maybe less so?
SEEDOPT='-s 2026'

# run lmm-filter!
time . ../scripts/lmm-filter.bash

# compare new output to most recent one
#zdiff lmm-filter{,_EXPECTED}/1e-02/preds.txt.gz 
