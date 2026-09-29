# if next path is correct, adjusted, or removed (if not needed), run test this way:
# bash test.bash

# more hacks to run saige from pixi environment
# NOTE: also need to install other R packages in pixi env: optparse, genio, tidyverse, pak, then platformbias from github (for now)
shopt -s expand_aliases
alias Rscript="pixi run -m ~/bin/src/github/SAIGE/ Rscript"


## STANDARD test

time Rscript ../scripts/lmm-filter.R --bfile test --platform platform.txt -s 2026

# compare new output to most recent one, with tolerance for p-value precision and runtime variance
Rscript compare_outputs.R
Rscript validate_preds.R

# cleanup!
rm -r lmm-filter/



## non-default IID,FID,PLATFORM column names in platform file

# use exact file from simulation, which we want to work with
# since platform assignments are identical, the answer should match standard run

time Rscript ../scripts/lmm-filter.R --bfile test --platform covar.txt -s 2026 --iid iid --fid famid --platform_col pheno

# compare new output to most recent one, with tolerance for p-value precision and runtime variance
Rscript compare_outputs.R
Rscript validate_preds.R

# cleanup!
rm -r lmm-filter/



## NOFLIP version!

time Rscript ../scripts/lmm-filter.R --bfile test --platform platform.txt -s 2026 --noflip
Rscript compare_outputs.R --noflip
Rscript validate_preds.R --noflip
rm -r lmm-filter/
