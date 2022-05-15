# ------------------------- #
# QC of FASTQ and BAM files #
# ------------------------- #

rule fastqc_fastq:
    input:
        path_to_data + '/{cohort}/tmp/gunzip_fastq/{sample}_lib{lib}_R{read}.fastq'
    output:
        html=path_to_data + '/{cohort}/results/qc_input/{sample}_lib{lib}_R{read}_fastqc.html',
        zipfile=path_to_data + '/{cohort}/results/qc_input/{sample}_lib{lib}_R{read}_fastqc.zip'
    resources: cpus=1, mem_mb=8000, time_min='24:00:00'
    params:
        outdir = lambda wildcards, output: '/'.join(output.html.split('/')[0:-1])
    conda: '../conda_env/fastqc.yml'
    shell:
        'fastqc --outdir {params.outdir} {input}'

# Run FASTQC on final BAM files.
rule fastqc_bam:
    input:
        path_to_data + '/{cohort}/results/bam_markdup/{sample}.aligned.sorted.markdup.bam',
    output:
        html = path_to_data + '/{cohort}/results/qc/{sample}.aligned.sorted.markdup_fastqc.html',
        zipfile = path_to_data + '/{cohort}/results/qc/{sample}.aligned.sorted.markdup_fastqc.zip',
    resources: cpus=1, mem_mb=8000, time_min='24:00:00'
    params:
        outdir = lambda wildcards, output: '/'.join(output.html.split('/')[0:-1])
    conda: '../conda_env/fastqc.yml'
    shell:
        'fastqc --outdir {params.outdir} {input}'

# Run Qualimap QC metrics on final BAM files.
#rule qualimap_bam:
#    input:
#        path_to_data + '/{cohort}/results/bam_markdup/{sample}.aligned.sorted.markdup.bam',
#    output:
#        html = path_to_data + '/{cohort}/results/qc/{sample}.aligned.sorted.markdup_fastqc.html',

# Run Flagstat to get basic stats on final BAM files.
rule bam_flagstat:
    input:
        path_to_data + '/{cohort}/results/bam_markdup/{sample}.aligned.sorted.markdup.bam',
    output:
        path_to_data + '/{cohort}/results/qc/{sample}.aligned.sorted.markdup.bam.flagstat',
    resources: cpus=1, mem_mb=8000, time_min='1:00:00'
    conda: '../conda_env/samtools.yml'
    shell:
        'samtools flagstat {input} > {output}'

# Unified QC rule which runs all of the above QCs
rule bam_qc:
    input:
        fastqc = path_to_data + '/{cohort}/results/qc/{sample}.aligned.sorted.markdup_fastqc.html',
        flagstat = path_to_data + '/{cohort}/results/qc/{sample}.aligned.sorted.markdup.bam.flagstat',
    output:
        path_to_data + '/{cohort}/results/qc/{sample}_qc_complete',
    resources: cpus=1, mem_mb=1000, time_min='00:01:00'
    conda: '../conda_env/samtools.yml'
    shell:
        'touch {output}'

