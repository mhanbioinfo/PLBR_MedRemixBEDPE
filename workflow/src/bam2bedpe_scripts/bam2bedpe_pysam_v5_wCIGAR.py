#!/usr/bin/env python
## bam2bedpe with FLAG, NM, TLEN
## USAGE:
## python bam2bedpe_pysam_v4.py [--sort_bam_by_qname] --bam_input --bedpe_output

import os
import pysam
import argparse
from copy import deepcopy
from timeit import default_timer as timer
import time

parser = argparse.ArgumentParser()
parser.add_argument('--sort_bam_by_qname', help='Will start with `samtools sort -n ...`', action='store_true')
parser.add_argument('--bam_input', help='Specify full path to BAM file.', required=True)
parser.add_argument('--bedpe_output', help='Specify full path of bedpe file.', required=True)
args = parser.parse_args()

sort_bam = args.sort_bam_by_qname
bam_input_path = args.bam_input
bedpe_output_path = args.bedpe_output

output_dir = os.path.dirname(bedpe_output_path)

## sort bam by queryname
if sort_bam is True:
    print("Sorting bam by queryname... ")
    bam_namedSortd_path = os.path.join(output_dir, os.path.splitext(os.path.basename(bam_input_path))[0] + "_nameSortd.bam")
    # print(bam_namedSortd_path)
    pysam.sort(bam_input_path, "-n", "-o", bam_namedSortd_path)
    bam_input_path = bam_namedSortd_path


## create new bedpe file
open(bedpe_output_path, 'w').close()

## start timer
start = timer()

## ------------------------------------------------------------ ##
## function to write out bedpe from bam
def write_bedpe(frag_dict):
    for key, an_fragmnt in frag_dict.items():

        if len(an_fragmnt) == 1:
            print("Only 1 read for fragment: ", key)

        bedpe_file_out = open(bedpe_output_path, "a")
        # a_fragmnt = frags_dict['NB551051:196:H3GGWBGXH:1:12104:6290:3128']
        
        ## deal with flags 329_393_377_441 reads
        a_fragmnt_alignmts_lst = [] ## for removing duplicate lines later
        for read_a in an_fragmnt:
            
            ## deal with flags 329_393_377_441 reads
            if read_a.flag in [329,393,377,441]:
                a_fragmnt_alignmts_lst.append("".join(str(x) for x in [ \
                    read_a.reference_name, "\t", read_a.reference_start, "\t", read_a.reference_start, "\t", \
                    read_a.reference_name, "\t", read_a.reference_start, "\t", read_a.reference_start, "\t", \
                    read_a.query_name, "\t", read_a.mapping_quality, "\t", read_a.mapping_quality, "\t", \
                    ".", "\t", ".", read_a.cigarstring, "\t", read_a.cigarstring, "\t", \
                    read_a.flag, "\t", read_a.flag, "\t", \
                    read_a.template_length, "\t", read_a.template_length, "\t",
                    "NA", "\t", "NA", "\n"]))

            for read_b in an_fragmnt:
                if read_a != read_b:
                    if read_a.next_reference_name == read_b.reference_name and \
                    read_a.next_reference_start == read_b.reference_start:
                        
                        ## define chr, chr_end if unmapped
                        if read_a.is_unmapped:
                            read_a_ref_name = "."
                            read_a_ref_start = "-1"
                            read_a_ref_end = "-1"
                            read_a_cigar = "*"
                        else:
                            read_a_ref_name = read_a.reference_name
                            read_a_ref_start = read_a.reference_start
                            read_a_ref_end = read_a.reference_end
                            read_a_cigar = read_a.cigarstring
                        if read_b.is_unmapped:
                            read_b_ref_name = "."
                            read_b_ref_start = "-1"
                            read_b_ref_end = "-1"
                            read_b_cigar = "*"
                        else:
                            read_b_ref_name = read_b.reference_name
                            read_b_ref_start = read_b.reference_start
                            read_b_ref_end = read_b.reference_end
                            read_b_cigar = read_b.cigarstring
                        
                        ## define strand
                        if read_a.is_reverse: read_a_strand = "-"
                        else: read_a_strand = "+"
                        if read_a_ref_name == ".": read_a_strand = "."

                        if read_b.is_reverse: read_b_strand = "-"
                        else: read_b_strand = "+"
                        if read_b_ref_name == ".": read_b_strand = "."
                        
                        ## define NM tag
                        if 'NM' in dict(read_a.tags): read_a_NM_tag = dict(read_a.tags)['NM']
                        else: read_a_NM_tag = "NA"
                        if 'NM' in dict(read_b.tags): read_b_NM_tag = dict(read_b.tags)['NM']
                        else: read_b_NM_tag = "NA"
                        
                        ## read1 always in col1-3, same as bedtools bedpe -mate1
                        if read_a.is_read1:
                            # pass
                            a_fragmnt_alignmts_lst.append("".join(str(x) for x in [ \
                                read_a_ref_name, "\t", read_a_ref_start, "\t", read_a_ref_end, "\t", \
                                read_b_ref_name, "\t", read_b_ref_start, "\t", read_b_ref_end, "\t", \
                                read_a.query_name, "\t", read_a.mapping_quality, "\t", read_b.mapping_quality, "\t", \
                                read_a_strand, "\t", read_b_strand, "\t", \
                                read_a_cigar, "\t", read_b_cigar, "\t", \
                                read_a.flag, "\t", read_b.flag, "\t", \
                                read_a.template_length, "\t", read_b.template_length, "\t",
                                read_a_NM_tag, "\t", read_b_NM_tag, "\n"]))
                        elif read_b.is_read1:
                            # pass
                            a_fragmnt_alignmts_lst.append("".join(str(x) for x in [ \
                                read_b_ref_name, "\t", read_b_ref_start, "\t", read_b_ref_end, "\t", \
                                read_a_ref_name, "\t", read_a_ref_start, "\t", read_a_ref_end, "\t", \
                                read_b.query_name, "\t", read_b.mapping_quality, "\t", read_a.mapping_quality, "\t", \
                                read_b_strand, "\t", read_a_strand, "\t", \
                                read_b_cigar, "\t", read_a_cigar, "\t", \
                                read_b.flag, "\t", read_a.flag, "\t", \
                                read_b.template_length, "\t", read_a.template_length, "\t",
                                read_b_NM_tag, "\t", read_a_NM_tag, "\n"]))
                        
        for line in list(set(a_fragmnt_alignmts_lst)):
            # print(line)
            bedpe_file_out.write(line)
        bedpe_file_out.close()

## ------------------------------------------------------------ ##
## iterate over pysam generator

## Initial call to print 0% progress
#samfile1 = pysam.AlignmentFile(bam_input_path) ## just to get len of generator
#loop_len = len(list(samfile1))
#printProgressBar(0, loop_len, prefix = 'Progress:', suffix = 'Complete', length = 50)

## read in pysam generator
samfile = pysam.AlignmentFile(bam_input_path)

qname_previous = ""
frag_dict = {}
frag_next_dict = {}

## loop through generator
for i, line in enumerate(samfile):
    
    qname = line.query_name
    # print(qname)
    
    ## just to get it started
    if qname_previous == "":
        frag_dict[qname] = [line]
        qname_previous = deepcopy(qname)
        continue
    
    ## compare current line to previous line
    if qname == qname_previous:
        frag_dict[qname].append(line)
        qname_previous = deepcopy(qname)
    
    ## as soon as current line qname != previous line qname,
    ## save iter to next fragment dict, and
    ## execute write_bedpe()
    else:
        frag_next_dict[qname] = [line]
        # print(frag_next_dict)
        
        write_bedpe(frag_dict)
        
        qname_previous = deepcopy(qname)
        frag_dict = deepcopy(frag_next_dict)
        frag_next_dict = {}
    
## write final fragment
write_bedpe(frag_dict)

## ------------------------------------------------------------ ##
## sort bedpe by coordinates
#bedpe_sortd_path = os.path.splitext(bedpe_output_path)[0] + "_coordSortd.bedpe"
#os.system("sort -k1,1V -k2,2n -k3,3n -k5,5n -k6,6n %s > %s" % (bedpe_output_path, bedpe_sortd_path))

## ------------------------------------------------------------ ##

end = timer()
#print("Processing bam2bedpe took ", str(end - start), " seconds.")


## EOF
