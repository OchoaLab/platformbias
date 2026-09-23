#!/bin/bash

run_saige() {
    local input_bfile=$1
    local covar=$2
    local output_prefix=$3

    echo "SAIGE step1..."
    time Rscript $SCRIPT_DIR/saige_step1_nocovar.R -f "$input_bfile" -c "$covar" -o "$output_prefix" > ${output_prefix}_step1.log $SEEDOPT
    echo "SAIGE step2..."
    time Rscript $SCRIPT_DIR/saige_step2.R -f "$input_bfile" -o "$output_prefix" > ${output_prefix}_step2.log $SEEDOPT
    
    # cleanup!
    rm ${output_prefix}{_output.txt.index,.rda,.varianceRatio.txt}
}

run_identify_sig_snps() {
    local phase=$1
    local index=$2
    local pval=$3

    Rscript $SCRIPT_DIR/snp_remove.R -l "$phase" -i "$index" -p "$pval"
}

run_plink_remove() {
    local input_bfile=$1
    local exclude_file=$2
    local output_prefix=$3
    
    plink2 --bfile "$input_bfile" \
        --exclude "$exclude_file" \
        --make-bed --out "$output_prefix" \
	--silent

    # cleanup: we don't need input anymore unless it's the original file!
    if [[ $input_bfile != "$input_data" ]]; then
	rm $input_bfile.{bed,bim,fam,log}
	# TODO: do we still need ${input_bfile}_output.txt ???
    fi
}

run_plink_flip() {
    local input_data=$1
    local exclude_file=$2
    local output_prefix=$3
    local flip_file=$4
    local flip_id_file=$5
    
    plink2 --bfile "$input_data"  \
	 --flip "$flip_file" \
	 --flip-subset "$flip_id_file" \
	 --exclude "$exclude_file" \
	 --silent \
	 --make-bed --out "$output_prefix"
}

run_phase2_flip_snps() {
    local input=$1
    Rscript $SCRIPT_DIR/snp_flip.R -f "$input"
}
