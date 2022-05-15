#!/bin/bash

#SBATCH --mail-type=ALL
#SBATCH --mail-user=ming.han@uhn.ca
#SBATCH -t 1-00:00:00
#SBATCH -D ./logs_slurm/
#SBATCH --mem=60G
#SBATCH -J get_clean_bedpe4medremix
#SBATCH -p himem
#SBATCH -c 4
#SBATCH -N 1
#SBATCH -o ./%j-%x.out
#SBATCH -e ./%j-%x.err

# getopts ###################################################
usage(){
    echo 
    echo "Usage: bash get_clean_bedpe4medremix.sh -i INPUT_PATH -s SAMPLE_NAME -p SPECIES -c CHR -f FRAGMENT_LENGTH_LIMIT -o OUT_DIR" 
    echo 
}
no_args="true"

## Help 
Help()
{ 
    # Display Help
    echo 
    echo "Get clean bedpe for input into MedRemixBEDPE."
    echo
    echo "Usage: bash get_clean_bedpe4medremix.sh -i INPUT_PATH -s SAMPLE_NAME -p SPECIES -c CHR -f FRAGMENT_LENGTH_LIMIT -o OUT_DIR"
    echo "options:"
    echo "-h   [HELP]      print help"
    echo "-i   [REQUIRED]  input file path for bedpe.gz (full path)"
    echo "-s   [REQUIRED]  sample name"
    echo "-p   [REQUIRED]  species (e.g. human, arabidopsis_F19K16, or arabidopsis_F24B22)"
    echo "-c   [REQUIRED]  chromosome to process (e.g. chr2)"
    echo "-f   [REQUIRED]  fragment length limit (e.g. 500)"
    echo "-o   [REQUIRED]  bedpe_bin_stats output directory (full path)"
    echo
}

## Get the options
while getopts ":hi:s:p:c:f:o:" option; do
    case "${option}" in
        h) Help
           exit;;
        i) INPUT_PATH=${OPTARG};;
        s) SAMPLE_NAME=${OPTARG};;
        p) SPECIES=${OPTARG};;
        c) CHR=${OPTARG};;
        f) FRAGMENT_LENGTH_LIMIT=${OPTARG};;
        o) OUT_DIR=${OPTARG};;
       \?) echo "Error: Invalid option"
           exit;;
    esac
    no_args="false"
done

[[ "$no_args" == "true" ]] && { usage; exit 1; }


# Main program ##############################################

echo "Processing get_clean_bedpe4medremix... " 
echo "Job started at "$(date) 
time1=$(date +%s)

INPUT_DIR=${INPUT_PATH%/*}
INPUT_BEDPE_GZ=${INPUT_PATH##*/}
BEDPE_CHR="${SAMPLE_NAME}_${SPECIES}_${CHR}.bedpe"
BEDPE_NOAMBI_NO113n177_mapQge20_isMate="${BEDPE_CHR%.*}.noAmbiguous.no113n177.mapQge20.isMate.bedpe"
BEDPE_lean="${BEDPE_NOAMBI_NO113n177_mapQge20_isMate%.*}.lean"
BEDPE_lean_withEnd="${BEDPE_lean}.withEnd"
BEDPE_lenFiltd="${BEDPE_lean_withEnd}.lenFiltd.bedpe4medremix"

## get specific chr whether as mate1 or mate2
zcat ${INPUT_DIR}/${INPUT_BEDPE_GZ} \
    | awk -v CHR=$CHR '($1 == CHR) || ($4 == CHR) {print}' \
    > ${OUT_DIR}/${BEDPE_CHR}

## 4 steps filtering:
## remove 'ambiguous reads' as defined by Rsamtools
  ## i.e. remove any .bedpe alignment if flag_mate1or2==SUPP && chr_mate1==chr_mate2
  ## ( chr_mate1!=chr_mate2 supp reads dealt with separately )
## remove 'FLAG 113,177' tandem reads
## remove fragments if either mate's MAPQ < 20
## remove unmated alignments 
  ## ( as defined by Rsamtools: read1-read2; both secondary, or none; no wrongly oriented alignments )
cat ${OUT_DIR}/${BEDPE_CHR} \
    | awk '((and($14,0x800) || and($15,0x800)) && $1==$4) {next} {print}' \
    | awk '(($14 == 113 && $15 == 177) || ($14 == 177 && $15 == 113)) {next} {print}' \
    | awk '($8>=20 && $9>=20) {print}' \
    | awk '(and($14,0x40) && and($15,0x80)) ||
           (and($14,0x80) && and($15,0x40)) {print}' \
    | awk '(!and($14,0x100) && !and($15, 0x100)) ||
           (and($14,0x100) && and($15, 0x100)) {print}' \
    | awk '(and($14,0x10) && and($15,0x20)) ||
           (and($14,0x20) && and($15,0x10)) {print}' \
    > ${OUT_DIR}/${BEDPE_NOAMBI_NO113n177_mapQge20_isMate}

## get lean bedpe for medremix and paste with qwidth for mate1 and mate2
## { qwidth calculated ignoring DHNPXB in CIGAR, as defined by Rsamtools }
paste \
<(cat ${OUT_DIR}/${BEDPE_NOAMBI_NO113n177_mapQge20_isMate} \
    | awk 'BEGIN{OFS="\t"}{print $1, $7, $2+1, $5+1, ($8+$9)/2, $10, $11}') \
<(cat ${OUT_DIR}/${BEDPE_NOAMBI_NO113n177_mapQge20_isMate} \
    | cut -f12 \
 	| awk '{gsub(/[1-9][0-9]*[D|H|N|P|X|B]/,"",$1); print $1}' \
    | awk 'BEGIN{OFS=""}{gsub(/[A-Z]/,"+",$1); print $1,0}' \
    | bc) \
<(cat ${OUT_DIR}/${BEDPE_NOAMBI_NO113n177_mapQge20_isMate} \
    | cut -f13 \
    | awk '{gsub(/[1-9][0-9]*[D|H|N|P|X|B]/,"",$1); print $1}' \
    | awk 'BEGIN{OFS=""}{gsub(/[A-Z]/,"+",$1); print $1,0}' \
    | bc) \
> ${OUT_DIR}/${BEDPE_lean}

## use qwidth to calculate 'end'
cat ${OUT_DIR}/${BEDPE_lean} \
    | awk 'BEGIN{OFS="\t"} ($6=="+") {print $1,$2, $3,$4, $5, $6,$7, $8,$9}
                           ($6=="-") {print $1,$2, $4,$3, $5, $7,$6, $9,$8}' \
    | awk 'BEGIN{OFS="\t"} {print $1,$2, $6":"$3"-"$3+$8","$7":"$4"-"$4+$9,
                                  $3, $4+$9, $4+$9-$3, $5}' \
    > ${OUT_DIR}/${BEDPE_lean_withEnd}

## filter by FRAGMENT_LENGTH_LIMIT
cat ${OUT_DIR}/${BEDPE_lean_withEnd} \
    | awk -v FRAGMENT_LENGTH_LIMIT=$FRAGMENT_LENGTH_LIMIT 'BEGIN{OFS="\t"} {if(($6 > 0) && ($6 < FRAGMENT_LENGTH_LIMIT)) {print}}' \
    > ${OUT_DIR}/${BEDPE_lenFiltd}

## remove intermediate files
rm ${OUT_DIR}/${BEDPE_CHR}
rm ${OUT_DIR}/${BEDPE_NOAMBI_NO113n177_mapQge20_isMate}
rm ${OUT_DIR}/${BEDPE_lean}
rm ${OUT_DIR}/${BEDPE_lean_withEnd}



time2=$(date +%s)
echo "Job ended at "$(date) 
echo "Job took $(((time2-time1)/3600)) hours $((((time2-time1)%3600)/60)) minutes $(((time2-time1)%60)) seconds"
echo ""
## EOF
