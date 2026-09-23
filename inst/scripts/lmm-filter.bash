# determine location of this script (and its dependencies), absolute path!
SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &> /dev/null && pwd)

# TODO: check that $PVAL, $input_data and $platform_file are defined!
# add absolute paths to all files, so they keep working as we change paths
# this works even on strings without complete extensions like this one, and doesn't change absolute paths
input_data=$(readlink -f "$input_data")
# used by saige only
platform_file=$(readlink -f "$platform_file")

# load scripts
source $SCRIPT_DIR/functions.sh

# all outputs should go in a subdirectory
# NOTE all existing paths are now absolute, no need to adjust for this move
mkdir lmm-filter
cd lmm-filter

#############################################################################
# check if iteration 0 saige output exists, if not, run saige for iteration 0
ITER=0

# iteration 0 is the same for all p-values, put it in base directory to share automatically
file_path="saige_0_output.txt.gz"
if [ -f "$file_path" ]; then
  echo "SAIGE output for iteration 0 already exists. Skipping initial SAIGE execution"
else
  run_saige "$input_data" "$platform_file" "saige_$ITER"
fi
# HACK OUTPUTS: saige_0_output.txt
# DELETED IMMEDIATELY: saige_0_output.txt.index  saige_0.rda  saige_0.varianceRatio.txt

# more cleanup
# assuming we didn't have any errors, should be safe to delete logs
rm saige_0_step?.log
# will reuse this one a lot, save!
gzip saige_0_output.txt

if [ -d "$PVAL" ]; then
  # if the directory exists, delete (to avoid overwriting existing files)
  rm -r "$PVAL"
fi
mkdir $PVAL
cd $PVAL

#############################################################################
# start phase 1
echo "Starting Phase 1"

# create summary file and record 'remove' for each iteration
summary_file="iteration_summary.txt"
echo -e "phase\tITER\tremove" > "$summary_file"

while :; do
# Rscript identifies significant SNPs; writes remove_phase1_$ITER.txt
  run_identify_sig_snps "phase1" "$ITER" "$PVAL"

  remove_file="remove_phase1_$ITER.txt"
  remove_count=$(wc -l < "$remove_file")

  # stop loop when there are zero removals
  if (( remove_count == 0 )); then
    echo "No removals at ITER=$ITER. Exiting loop."
    break
  fi
  # else update log
  echo -e "1\t$ITER\t$remove_count" >> "$summary_file"

  # plink removes SNPs
  if [ "$ITER" -eq 0 ]; then
    input_bfile="$input_data"
  else
    # get previous iteration file as the input file for plink processing
    input_bfile=$((ITER - 1))
  fi
  run_plink_remove "$input_bfile" "$remove_file" "$ITER"

  # Run SAIGE
  run_saige "$ITER" "$platform_file" "saige_phase1_$ITER"

  # Increment for next iteration
  ITER=$((ITER + 1))
done

# more cleanup
# delete last bed/bim/fam files (not used in phase2)
rm $((ITER - 1)).{bed,bim,fam,log}

#############################################################################
# start phase 2
echo "Starting Phase 2"

# identify snps to flip and remove
run_phase2_flip_snps "$input_data"

# process flip and removal of SNPs, created new plink files

# TODO: flip-subset file should be created automatically, somehow...
run_plink_flip "$input_data" phase2_init_remove.txt 0 phase2_init_flip.txt "$platform_file"

# output file written as 'control_0' in phase2 folder, next step is to run saige
ITER=0

while :; do
  run_saige "$ITER" "$platform_file" "saige_phase2_$ITER"

  # identify significant SNPs
  run_identify_sig_snps "phase2" "$ITER" "$PVAL"

  # Count removed SNPs
  remove_file="remove_phase2_$ITER.txt"
  remove_count=$(wc -l < "$remove_file")
  
  # stop loop when there are zero removals
  if (( remove_count == 0 )); then
    echo "No removals at ITER=$ITER. Exiting loop."
    break
  fi
  # else update log
  echo -e "2\t$ITER\t$remove_count" >> "$summary_file"

  # remove snps with plink
  run_plink_remove "$ITER" "$remove_file" "$((ITER + 1))"

  # Increment for next iteration
  ITER=$((ITER + 1))
done

# delete BED files again.  These are not what we may want in the end because it's controls only.  In practice we want to process full data separately, following `preds.txt.gz`.
rm $ITER.{bed,bim,fam,log}

# produce final output of method, preds.txt.gz, which summarizes which loci are keep, remove, or flip.
Rscript $SCRIPT_DIR/lmm-classify.R -f "$input_data"

# if we're here, presumably everything was good and we don't need those bulky logs anymore... (they're only for troubleshooting)
rm saige_phase?_*_step?.log
# for now let's keep the full saige files, in case we wanted their p-values...
# (otherwise consider deleting them)
gzip saige_phase?_*_output.txt

# go back to starting directory, it's weird otherwise
cd ../..
