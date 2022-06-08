#!/bin/bash

PROJ_DIR="/cluster/projects/pughlab/projects/cfMeDIP_compare_pipelines/cfmedip_medremix_bedpe/analysis003_2toys_input_ftype"

## .bedgraph
echo "Writing to ${PROJ_DIR}/md5sum.bedgraph"
find ${PROJ_DIR} -type f -name "*.bedgraph" -exec md5sum {} \; | sort > ${PROJ_DIR}/md5sum.bedgraph

## .nbglm.tsv
echo "Writing to ${PROJ_DIR}/md5sum.nbglm.tsv"
find ${PROJ_DIR} -type f -name "*nbglm.tsv" -exec md5sum {} \; | sort > ${PROJ_DIR}/md5sum.nbglm.tsv

## qc_full
echo "Writing to ${PROJ_DIR}/md5sum.qc_full.txt"
find ${PROJ_DIR} -type f -name "*qc_full.txt" -exec sh -c 'echo "$(sed -n "5,\$p" $1 | md5sum), $1"' _ {} \; | sort > ${PROJ_DIR}/md5sum.qc_full.txt


## EOF
