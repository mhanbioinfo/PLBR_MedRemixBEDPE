# --------------------------- #
#  Align FASTQs to reference  #
# --------------------------- #

def get_read_group_from_fastq(fastq_file, sample_name):
    """Extracts the read group data from FASTQs automatically.
    This information is used in the shell script of rule bwa_mem.
    """
    with gzip.open(fastq_file, 'rt') as fastq:
        header = next(fastq)
        (instrument, run_number, flowcell, lane, tile, xpos, ypos) = header.split(' ')[0].split(':')
        lib_value = header.strip().split(' ')[1].split(':')[3]
        rg_line = r"@RG\tID:{flowcell}_{lane}\tSM:{sample}\tPL:Illumina\tPU:.\tLB:{lib_value}".format(
            flowcell = flowcell,
            lane = lane,
            sample = sample_name,
            lib_value = lib_value
        )
        return(rg_line)

## NO READ GROUP INFO (for SRA FASTQs)
# Run BWA mem on FASTQs after extracting barcodes.
rule bwa_mem:
    input:
        path_to_data + '/{cohort}/tmp/trim_fastq/{sample}_lib{lib}_extract_barcode_R1_val_1.fq',
        path_to_data + '/{cohort}/tmp/trim_fastq/{sample}_lib{lib}_extract_barcode_R2_val_2.fq',
    output:
        temp(path_to_data + '/{cohort}/tmp/bwa_mem/{sample}_lib{lib}.sam')
        #path_to_data + '/{cohort}/tmp/bwa_mem/{sample}_lib{lib}.sam'
    resources: cpus=4, mem_mb=16000, time_min='72:00:00'
    params:
        bwa_index = lambda wildcards: get_cohort_config(wildcards.cohort)['bwa_index']
    conda: '../conda_env/samtools.yml'
    shell:
        clean(r"""
        bwa mem -M -t4
        {params.bwa_index}
        {input} > {output}""")

# Converts SAM to BAM with sorting
rule sam_to_sorted_bam:
    input:
        path_to_data + '/{cohort}/tmp/bwa_mem/{sample}_lib{lib}.sam'
    output:
        bam = temp(path_to_data + '/{cohort}/tmp/bwa_mem/{sample}_lib{lib}.sorted.bam'),
        #bam = path_to_data + '/{cohort}/tmp/bwa_mem/{sample}_lib{lib}.sorted.bam',
        index = temp(path_to_data + '/{cohort}/tmp/bwa_mem/{sample}_lib{lib}.sorted.bam.bai'),
        #index = path_to_data + '/{cohort}/tmp/bwa_mem/{sample}_lib{lib}.sorted.bam.bai',
    resources: cpus=16, mem_mb=30000, time_min='72:00:00'
    conda: '../conda_env/samtools.yml'
    shell:
        # Try it without fixmate for replication purposes
        #"samtools view -buS -f 2 -F 4 -@4 {input} | samtools fixmate -m - - | samtools sort -@4 -o {output.bam} && samtools index {output.bam}"
        clean(r'''
        samtools view -buS -f 2 -F 4 -@32 {input} |
        samtools fixmate -m - - |
        samtools sort -@32 -o {output.bam} && samtools index {output.bam}
        ''')

def get_libraries_of_sample(sample):
    """Returns all library indices of a sample based on samplesheet."""
    filtered_table = get_all_samples()[get_all_samples().sample_name == sample]
    return(list(set(filtered_table.library_index.to_list())))

# If there are multiple libraries for a given sample, as specified in samplesheet,
# these libraries are automatically merged at this step into a single unified BAM.
rule merge_bam:
    input:
        lambda wildcards: expand(
                path_to_data + '/{{cohort}}/tmp/bwa_mem/' + wildcards.sample + '_lib{lib}.sorted.bam',
                lib=get_libraries_of_sample(wildcards.sample)
        )
    output:
        temp(path_to_data + '/{cohort}/tmp/merge_bam/{sample}.aligned.sorted.bam')
        #path_to_data + '/{cohort}/tmp/merge_bam/{sample}.aligned.sorted.bam'
    resources: cpus=1, mem_mb=8000, time_min='24:00:00'
    conda: '../conda_env/samtools.yml'
    shell:
        'samtools merge {output} {input} && samtools index {output}'

# Bam markdup and create index.
# This step finalizes the definitive BAM file.
rule bam_markdup:
    input:
        path_to_data + '/{cohort}/tmp/merge_bam/{sample}.aligned.sorted.bam'
    output:
        bam = path_to_data + '/{cohort}/results/bam_markdup/{sample}.aligned.sorted.markdup.bam',
        index = path_to_data + '/{cohort}/results/bam_markdup/{sample}.aligned.sorted.markdup.bam.bai'
    resources: cpus=1, mem_mb=8000, time_min='24:00:00'
    conda: '../conda_env/samtools.yml'
    shell:
        #"samtools markdup -r {input} {output.bam} && samtools index {output.bam}"
        "samtools markdup {input} {output.bam} && samtools index {output.bam}"

