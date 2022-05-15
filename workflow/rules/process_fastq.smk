# -------------------------- #
#  Pre-process input FASTQs  #
# -------------------------- #

def get_fastq_path(cohort, sample, library, read_in_pair=1):
    """Retrieves the path to the fastq file.

    Keyword arguments:
        cohort -- name of the cohort whose samplesheet should be accessed.
        sample -- identifier of the sample as specified in the samplesheet.
        library -- integer representing the library index as specified in the samplesheet.
        read_in_pair -- 1 or 2 - representing read 1 or read 2 in paired end data.
    """
    library = int(library)
    cohort_data = get_cohort_data(cohort)
    sample_line = cohort_data[
        (cohort_data.sample_name == sample) &
        (cohort_data.library_index == library) &
        (cohort_data.read_in_pair == read_in_pair)
    ]
    return(sample_line.path.to_list()[0])


rule gunzip_fastq:
    input:
        lambda wildcards: get_fastq_path(wildcards.cohort, wildcards.sample, int(wildcards.lib), int(wildcards.read)),
    output:
        temp(path_to_data + '/{cohort}/tmp/gunzip_fastq/{sample}_lib{lib}_R{read}.fastq')
        #path_to_data + '/{cohort}/tmp/gunzip_fastq/{sample}_lib{lib}_R{read}.fastq'
    resources: cpus=1, mem_mb=8000, time_min='24:00:00'
    shell:
        'gunzip -dc {input} > {output}'

# Extract Barcodes using ConsensusCruncher
# Pulls the path to extract_barcodes.py from config > paths > dependencies > extract_barcodes_path
rule extract_barcodes:
    input:
        R1_qc = path_to_data + '/{cohort}/results/qc_input/{sample}_lib{lib}_R1_fastqc.html',
        R2_qc = path_to_data + '/{cohort}/results/qc_input/{sample}_lib{lib}_R2_fastqc.html',
        R1 = path_to_data + '/{cohort}/tmp/gunzip_fastq/{sample}_lib{lib}_R1.fastq',
        R2 = path_to_data + '/{cohort}/tmp/gunzip_fastq/{sample}_lib{lib}_R2.fastq'
    output:
        R1 = temp(path_to_data + '/{cohort}/tmp/extract_barcodes/{sample}_lib{lib}_extract_barcode_R1.fastq'),
        #R1 = path_to_data + '/{cohort}/tmp/extract_barcodes/{sample}_lib{lib}_extract_barcode_R1.fastq',
        R2 = temp(path_to_data + '/{cohort}/tmp/extract_barcodes/{sample}_lib{lib}_extract_barcode_R2.fastq')
        #R2 = path_to_data + '/{cohort}/tmp/extract_barcodes/{sample}_lib{lib}_extract_barcode_R2.fastq'
    params:
        outprefix = lambda wildcards, output: output.R1.split('_barcode_')[0],
        barcodes = lambda wildcards: get_cohort_config(wildcards.cohort)['barcodes']
    resources: cpus=1, mem_mb=16000, time_min='5-00:00:00'
    conda: '../conda_env/biopython.yml'
    shell:
        clean(r'''
        python {extract_barcodes}
            --read1 {{input.R1}}
            --read2 {{input.R2}}
            --outfile {{params.outprefix}}
            {{params.barcodes}}
        '''.format(extract_barcodes = config['paths']['dependencies']['extract_barcodes_path']))

# Trims FASTQ using trimgalore to remove barcode sequences
# By default, trims 10 base pairs from the 5' end, which seems to be correct for OICR cfMeDIP-seq output.
# This can be configured in the config.yml under data > cohorts > settings > trimgalore.
rule trim_fastq:
    input:
        R1 = path_to_data + '/{cohort}/tmp/extract_barcodes/{sample}_lib{lib}_extract_barcode_R1.fastq',
        R2 = path_to_data + '/{cohort}/tmp/extract_barcodes/{sample}_lib{lib}_extract_barcode_R2.fastq'
    output:
        #trimmed_1 = temp(path_to_data + '/{cohort}/tmp/trim_fastq/{sample}_lib{lib}_extract_barcode_R1_val_1.fq'),
        trimmed_1 = path_to_data + '/{cohort}/tmp/trim_fastq/{sample}_lib{lib}_extract_barcode_R1_val_1.fq',
        #trimmed_2 = temp(path_to_data + '/{cohort}/tmp/trim_fastq/{sample}_lib{lib}_extract_barcode_R2_val_2.fq'),
        trimmed_2 = path_to_data + '/{cohort}/tmp/trim_fastq/{sample}_lib{lib}_extract_barcode_R2_val_2.fq',
        report_1 = path_to_data + '/{cohort}/results/qc/{sample}_lib{lib}_extract_barcode_R1.fastq_trimming_report.txt',
        report_2 = path_to_data + '/{cohort}/results/qc/{sample}_lib{lib}_extract_barcode_R2.fastq_trimming_report.txt'
    params:
        outdir = lambda wildcards, output: '/'.join(output.trimmed_1.split('/')[0:-1]),
        trimgalore_settings = lambda wildcards: get_cohort_config(wildcards.cohort)['trimgalore']
    resources: cpus=4, mem_mb=8000, time_min='24:00:00'
    conda: '../conda_env/trimgalore.yml'
    shell:
        'trim_galore --cores 4 --dont_gzip --paired {params.trimgalore_settings} --output_dir {params.outdir} {input.R1} {input.R2} && cp ' + path_to_data + '/{wildcards.cohort}/tmp/trim_fastq/{wildcards.sample}_lib{wildcards.lib}_extract_barcode_R*.fastq_trimming_report.txt ' + path_to_data + '/{wildcards.cohort}/results/qc'
