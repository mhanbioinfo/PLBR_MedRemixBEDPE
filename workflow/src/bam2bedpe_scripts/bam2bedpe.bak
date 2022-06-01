#!/bin/bash

#SBATCH --mail-type=ALL
#SBATCH --mail-user=ming.han@uhn.ca
#SBATCH -t 1-00:00:00
#SBATCH --mem=16G
#SBATCH -J bam2bedpe
#SBATCH -p himem
#SBATCH -c 1
#SBATCH -N 1
#SBATCH -o %j-%x.out
#SBATCH -e %j-%x.err

# getopts ###################################################
usage(){
    echo 
    echo "Usage: bash/sbatch bam2bedpe.sh -s [slurm|local] -c num_of_chunks -b bam_input_path -o output_dir -x SRC_DIR -p PICARD_DIR -a CONDA -e ENV -t"
    echo 
}
no_args="true"
KEEP_TMP=false

## Help 
Help()
{
    # Display Help
    echo 
    echo "Processes bam to bedpe in chunks, preserving FLAG, TLEN and CIGAR fields."
    echo
    echo "Usage: bash/sbatch bam2bedpe.sh -s [slurm|local] -c num_of_chunks -b bam_input_path -o output_dir -x SRC_DIR -p PICARD_DIR -a CONDA -e ENV -t"
    echo "options:"
    echo "-h   [HELP]      print help"
    echo "-s   [REQUIRED]  type either 'slurm' or 'local', local is with nohup"
    echo "-c   [REQUIRED]  number of chunks to process in parallel"
    echo "-b   [REQUIRED]  path to bam input (full path)"
    echo "-o   [REQUIRED]  output directory (full path)"
    echo "-x   [REQUIRED]  src directory with pipeline scripts (full path)"
    echo "-p   [REQUIRED]  full path to directory containing picard.jar"
    echo "-a   [REQUIRED]  full path to conda activate (e.g. /cluster/home/t110409uhn/bin/miniconda3/bin/activate)"
    echo "-e   [REQUIRED]  conda env with pysam (e.g. pipeline-medremixBEDPE)"
    echo "-t   [OPTIONAL]  keep tmp_dir"
    echo
}

## Get the options
while getopts ":hs:c:b:o:x:p:a:e:t" option; do
    case "${option}" in
        h) Help
           exit;;
        s) SLURMLOCAL=${OPTARG};;
        c) NUM_OF_CHUNKS=${OPTARG};;
        b) INPUT_BAM_PATH=${OPTARG};;
        o) OUT_DIR=${OPTARG};;
        x) SRC_DIR=${OPTARG};;
        p) PICARD_DIR=${OPTARG};;
        a) CONDA_ACTIVATE=${OPTARG};;
        e) CONDA_ENV=${OPTARG};;
        t) KEEP_TMP=true;;
       \?) echo "Error: Invalid option"
           exit;;
    esac
    no_args="false"
done

[[ "$no_args" == "true" ]] && { usage; exit 1; }

echo "Processing step8_bam2bedpe..."
echo "number of chunks: $NUM_OF_CHUNKS"
echo "input bam path:   $INPUT_BAM_PATH"
echo "output path:      $OUT_DIR"
echo "processing on:    $SLURMLOCAL"

# Main program ##############################################

echo "Job started at "$(date) 
time1=$(date +%s)

#NUM_OF_CORES=$( nproc )
#echo "${NUM_OF_CORES} cores available, processing ${NUM_OF_CORES} chunks in parallel."
#NUM_OF_CHUNKS=$NUM_OF_CORES

PY_SCRIPT_DIR=${SRC_DIR}
PY_SCRIPT_PATH="${PY_SCRIPT_DIR}/bam2bedpe_pysam_v5_wCIGAR.py"
INPUT_DIR="${INPUT_BAM_PATH%/*}"
INPUT_BAM="${INPUT_BAM_PATH##*/}"

TMP_DIR="${OUT_DIR}/tmp_dir/${INPUT_BAM%.*}"
mkdir -p $TMP_DIR

OUT_FRAG_NAMES="${INPUT_BAM%.*}.fragNames"
OUT_MERGED_SORTD_BEDPE="${INPUT_BAM%.*}_coordSortd.bedpe"

## -------------------------------------- ##
## get fragNames
samtools view ${INPUT_DIR}/${INPUT_BAM} \
    | cut -f1 \
    | awk '{ a[$1]++ } END { for (b in a) { print b } }' \
    > ${TMP_DIR}/${OUT_FRAG_NAMES}

## -------------------------------------- ##
## split bam into chunks
NUM_ROWS=$(cat ${TMP_DIR}/${OUT_FRAG_NAMES} | wc -l)
CHUNK_SIZE=$(( $NUM_ROWS / $NUM_OF_CHUNKS ))   ## won't have decimal, sometimes will be short, last chunk has to go to last row
echo ".bam has $NUM_ROWS fragments."
echo "Number of chunks was set to ${NUM_OF_CHUNKS}."
echo "Each chunk will be $CHUNK_SIZE rows."

for ((i=0; i<=$(( $NUM_OF_CHUNKS - 2 )); i++)); do
    CHUNK=$(( i + 1 ))
    ROW_START=$(( i*CHUNK_SIZE + 1))
    ROW_END=$(( CHUNK*CHUNK_SIZE ))
    echo "Processing CHUNK${CHUNK}... starting with row ${ROW_START}, ending in row ${ROW_END}."
    sed -n "$ROW_START,$ROW_END p" ${TMP_DIR}/${OUT_FRAG_NAMES} \
        > ${TMP_DIR}/${OUT_FRAG_NAMES}.CHUNK${CHUNK}

done

## since chunk might not be evenly divided, last chunk will just be to last row
ith=$(( $NUM_OF_CHUNKS - 1 ))
ROW_START=$(( ith*CHUNK_SIZE + 1 ))
echo "Processing CHUNK$NUM_OF_CHUNKS... starting with row ${ROW_START}, ending in row ${NUM_ROWS}."
sed -n "$ROW_START,$NUM_ROWS p" ${TMP_DIR}/${OUT_FRAG_NAMES} \
    > ${TMP_DIR}/${OUT_FRAG_NAMES}.CHUNK${NUM_OF_CHUNKS}

## make log dir
mkdir -p ${OUT_DIR}/logs_slurm

## -------------------------------------- ##
## picard FilterSamReads on chunks
for CHUNK in ${TMP_DIR}/${OUT_FRAG_NAMES}.CHUNK*; do
    echo "Picard FilterSamReads for ${CHUNK##*/}..."
    echo "Output is ${INPUT_BAM%.*}_${CHUNK##*.}.bam"
    
    echo "Writing out ${INPUT_BAM%.*}_bam2bedpe_${CHUNK##*.}.sh"
    ## write out sample sbatch script
    cat <<- EOF > "${TMP_DIR}/${INPUT_BAM%.*}_bam2bedpe_${CHUNK##*.}.sh" 
	#!/bin/bash
	#SBATCH -t 3-00:00:00
	#SBATCH -J bam2bedpe_chunks_${CHUNK##*.}
	#SBATCH -D ${OUT_DIR}/logs_slurm
	#SBATCH --mail-type=ALL
	#SBATCH --mail-user=ming.han@uhn.ca
	#SBATCH -p himem
	#SBATCH -c 1
	#SBATCH --mem=16G
	#SBATCH -o ./%j-%x.out
	#SBATCH -e ./%j-%x.err
	
	source ${CONDA_ACTIVATE} ${CONDA_ENV}
	PICARD_DIR=${PICARD_DIR}
	
	echo "Job started at "\$(date) 
	time1=\$(date +%s)
	
	java -jar ${PICARD_DIR}/picard.jar FilterSamReads \
	    I=${INPUT_DIR}/${INPUT_BAM} \
	    O=${TMP_DIR}/"${INPUT_BAM%.*}_${CHUNK##*.}.bam" \
	    READ_LIST_FILE=${CHUNK} \
	    FILTER=includeReadList \
	    WRITE_READS_FILES=false \
	    USE_JDK_DEFLATER=true \
	    USE_JDK_INFLATER=true \
	
	python ${PY_SCRIPT_PATH} \
	    --sort_bam_by_qname \
	    --bam_input ${TMP_DIR}/"${INPUT_BAM%.*}_${CHUNK##*.}.bam" \
	    --bedpe_output ${TMP_DIR}/"${INPUT_BAM%.*}_${CHUNK##*.}.bedpe"
	
	rm ${TMP_DIR}/"${INPUT_BAM%.*}_${CHUNK##*.}.bam"
	
	time2=\$(date +%s)
	echo "Job ended at "\$(date) 
	echo "Job took \$(((time2-time1)/3600)) hours \$((((time2-time1)%3600)/60)) minutes \$(((time2-time1)%60)) seconds"
	EOF
    
    if [ $SLURMLOCAL == "slurm" ]; then
        sbatch "${TMP_DIR}/${INPUT_BAM%.*}_bam2bedpe_${CHUNK##*.}.sh"
    elif [ $SLURMLOCAL == "local" ]; then
        nohup bash "${TMP_DIR}/${INPUT_BAM%.*}_bam2bedpe_${CHUNK##*.}.sh" &> "${TMP_DIR}/${INPUT_BAM%.*}_bam2bedpe_${CHUNK##*.}.log" &
    fi

done

## periodically check if chunks have been written to completely
while true; do
    if [ $(find ${TMP_DIR} -maxdepth 1 -mmin +3 -type f -regex ".*_CHUNK[1-9][0-9]*\.bedpe" | wc -l) -eq $NUM_OF_CHUNKS ] 
    then
        echo "All bedpe chunks have been written to, merging bedpe chunks..."
        BEDPE_CHUNKS=$(find ${TMP_DIR} -maxdepth 1 -mmin +3 -type f -regex ".*_CHUNK[1-9][0-9]*\.bedpe" | sort -V)
        cat $BEDPE_CHUNKS \
            | sort -k1,1V -k2,2n -k3,3n -k5,5n -k6,6n \
            | gzip -c \
            > ${OUT_DIR}/${OUT_MERGED_SORTD_BEDPE}.gz
        break
    else
        echo "Sleeping for 5 minutes, then recheck if bedpe chunks were written to completely."
        sleep 5m
    fi
done

## remove TMP_DIR
if "$KEEP_TMP"; then
    echo "Keeping temp dir"
else
    echo "Removed temp dir after completion"
    rm -rf ${TMP_DIR}
fi


time2=$(date +%s)
echo "Job ended at "$(date) 
echo "Job took $(((time2-time1)/3600)) hours $((((time2-time1)%3600)/60)) minutes $(((time2-time1)%60)) seconds"
echo ""

## EOF
