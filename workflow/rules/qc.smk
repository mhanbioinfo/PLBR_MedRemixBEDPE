# ------------------------- #
# QC of FASTQ and BAM files #
# ------------------------- #

# Run FastQC on input fastq files.
rule fastqc_fastq:
    input:
        path_to_data + '/{cohort}/tmp/gunzip_fastq/{sample}_lib{lib}_R{read}.fastq'
    output:
        html=path_to_data + '/{cohort}/qc/fastqc_input_fastq/{sample}_lib{lib}_R{read}_fastqc.html',
        zipfile=path_to_data + '/{cohort}/qc/fastqc_input_fastq/{sample}_lib{lib}_R{read}_fastqc.zip'
    resources: cpus=1, mem_mb=8000, time_min='24:00:00'
    params:
        outdir = lambda wildcards, output: '/'.join(output.html.split('/')[0:-1])
    conda: '../conda_env/fastqc.yml'
    shell:
        'fastqc --outdir {params.outdir} {input}'

## Get hg38 only bam without F19K16 and F24B22 reads for QC purposes
rule get_hg38_bam:
    input:
        path_to_data + '/{cohort}/results/bam_markdup/{sample}.aligned.sorted.markdup.bam',
    output:
        path_to_data + '/{cohort}/results/bam_markdup/{sample}.aligned.sorted.markdup_hg38only.bam',
    resources: cpus=1, mem_mb=12000, time_min='5:00:00'
    conda: '../conda_env/samtools.yml'
    shell:
        'samtools view {input} | grep -v "F19K16\|F24B22" | cat <(samtools view -H {input} | grep -v "F19K16\|F24B22") - | samtools view -b - > {output} && samtools index {output}'

# Run FastQC on final (aligned, sorted, dup marked) BAM files.
rule fastqc_bam:
    input:
        path_to_data + '/{cohort}/results/bam_markdup/{sample}.aligned.sorted.markdup_hg38only.bam',
    output:
        html = path_to_data + '/{cohort}/qc/fastqc_bam/{sample}.aligned.sorted.markdup_hg38only_fastqc.html',
        zipfile = path_to_data + '/{cohort}/qc/fastqc_bam/{sample}.aligned.sorted.markdup_hg38only_fastqc.zip',
    resources: cpus=1, mem_mb=8000, time_min='24:00:00'
    params:
        outdir = lambda wildcards, output: '/'.join(output.html.split('/')[0:-1])
    conda: '../conda_env/fastqc.yml'
    shell:
        'fastqc --outdir {params.outdir} {input}'

# Run Flagstat to get basic stats on final (aligned, sorted, dup marked) BAM files.
rule bam_flagstat:
    input:
        path_to_data + '/{cohort}/results/bam_markdup/{sample}.aligned.sorted.markdup_hg38only.bam',
    output:
        path_to_data + '/{cohort}/qc/flagstat/{sample}.aligned.sorted.markdup_hg38only.bam.flagstat',
    resources: cpus=1, mem_mb=8000, time_min='1:00:00'
    conda: '../conda_env/samtools.yml'
    shell:
        'samtools flagstat {input} > {output}'

# Run MEDIPS QC on final (aligned, sorted, dup marked) BAM files.
rule medips_qc:
    input:
        path_to_data + '/{cohort}/results/bam_markdup/{sample}.aligned.sorted.markdup_hg38only.bam',
    output:
        csv = temp(path_to_data + '/{cohort}/qc/medips_qc/{sample}_QC_MEDIPS.csv'),
    params:
        outdir = lambda wildcards, output: '/'.join(output.csv.split('/')[0:-1]),
        winsize = config['pipeline_params']['window_size'],
    resources: cpus=1, mem_mb=8000, time_min='2:00:00'
    conda: '../conda_env/cfmedip_r.yml'
    shell:
        'Rscript src/QC/QC_MEDIPS.R --bamFile {input} --outputDir {params.outdir} --windowSize {params.winsize}'


# Run Picard QC on final (aligned, sorted, dup marked) BAM files.

## CollectGcBiasMetrics
rule CollectGcBiasMetrics:
    input:
        path_to_data + '/{cohort}/results/bam_markdup/{sample}.aligned.sorted.markdup_hg38only.bam',
    output:
        metric = path_to_data + '/{cohort}/qc/picard_qc/{sample}.gc_bias_metrics.txt',
        chart = path_to_data + '/{cohort}/qc/picard_qc/{sample}.gc_bias_metrics.pdf',
        summary = temp(path_to_data + '/{cohort}/qc/picard_qc/{sample}.gc_bias_summary_metrics.txt'),
    params:
        fasta = config['data']['defaults']['hg38only_genome'],
    resources: cpus=1, mem_mb=32000, time_min='5:00:00'
    conda: '../conda_env/samtools.yml'
    shell:
        'picard CollectGcBiasMetrics -I {input} -O {output.metric} -CHART {output.chart} -S {output.summary} -R {params.fasta} && bash src/QC/parse_picard_QC.sh picard.CollectGcBiasMetrics {output.summary} {output.summary}.parsed'

## CollectInsertSizeMetrics
rule CollectInsertSizeMetrics:
    input:
        path_to_data + '/{cohort}/results/bam_markdup/{sample}.aligned.sorted.markdup_hg38only.bam',
    output:
        metric = path_to_data + '/{cohort}/qc/picard_qc/{sample}.picardInsertSize_metrics.txt',
        histo = path_to_data + '/{cohort}/qc/picard_qc/{sample}.picardInsertSize_metrics.pdf',
    resources: cpus=1, mem_mb=16000, time_min='5:00:00'
    conda: '../conda_env/samtools.yml'
    shell:
        'picard CollectInsertSizeMetrics -INCLUDE_DUPLICATES true -I {input} -O {output.metric} -H {output.histo} -M 0 -W 600 && bash src/QC/parse_picard_QC.sh picard.CollectInsertSizeMetrics {output.metric} {output.metric}.parsed'

## MarkDuplicates
rule MarkDuplicates_QC:
    input:
        path_to_data + '/{cohort}/results/bam_markdup/{sample}.aligned.sorted.markdup_hg38only.bam',
    output:
        markdup_bam = temp(path_to_data + '/{cohort}/qc/picard_qc/{sample}.picardMarkDup.bam'),
        metric = temp(path_to_data + '/{cohort}/qc/picard_qc/{sample}.marked_dup_metrics.txt'),
    resources: cpus=1, mem_mb=32000, time_min='5:00:00'
    conda: '../conda_env/samtools.yml'
    shell:
        'picard MarkDuplicates -I {input} -O {output.markdup_bam} -M {output.metric} && bash src/QC/parse_picard_QC.sh picard.MarkDuplicates {output.metric} {output.metric}.parsed'

## CollectAlignmentSummaryMetrics
rule CollectAlignmentSummaryMetrics:
    input:
        path_to_data + '/{cohort}/results/bam_markdup/{sample}.aligned.sorted.markdup_hg38only.bam',
    output:
        metric = temp(path_to_data + '/{cohort}/qc/picard_qc/{sample}.alignmentMetrics.txt'),
    params:
        fasta = config['data']['defaults']['hg38only_genome'],
    resources: cpus=1, mem_mb=32000, time_min='5:00:00'
    conda: '../conda_env/samtools.yml'
    shell:
        'picard CollectAlignmentSummaryMetrics -I {input} -O {output.metric} -R {params.fasta} && bash src/QC/parse_picard_QC.sh picard.CollectAlignmentSummaryMetrics {output.metric} {output.metric}.parsed'

## CollectQualityYieldMetrics
rule CollectQualityYieldMetrics:
    input:
        path_to_data + '/{cohort}/results/bam_markdup/{sample}.aligned.sorted.markdup_hg38only.bam',
    output:
        metric = temp(path_to_data + '/{cohort}/qc/picard_qc/{sample}.quality_yield_metrics.txt'),
    resources: cpus=1, mem_mb=32000, time_min='5:00:00'
    conda: '../conda_env/samtools.yml'
    shell:
        'picard CollectQualityYieldMetrics -I {input} -O {output.metric} && bash src/QC/parse_picard_QC.sh picard.CollectQualityYieldMetrics {output.metric} {output.metric}.parsed'

## { complains of missing files 'toy02.Rrbs completed successfully, but some output files are missing.' }
### CollectRrbsMetrics
#rule CollectRrbsMetrics:
#    input:
#        path_to_data + '/{cohort}/results/bam_markdup/{sample}.aligned.sorted.markdup_hg38only.bam',
#    output:
#        metric = path_to_data + '/{cohort}/qc/picard_qc/{sample}.Rrbs',
#        #picard_rrbs_detailed_metric = path_to_data + '/{cohort}/qc/picard_qc/{sample}.Rrbs.rrbs_detail_metrics',
#        #picard_rrbs_qc_pdf = path_to_data + '/{cohort}/qc/picard_qc/{sample}.Rrbs.rrbs_qc.pdf',
#        #picard_rrbs_summary = path_to_data + '/{cohort}/qc/picard_qc/{sample}.Rrbs.rrbs_summary_metrics',
#    params:
#        picard_dir = config['pipeline_params']['picard_dir'],
#        fasta = config['data']['defaults']['hg38only_genome'],
#    resources: cpus=1, mem_mb=32000, time_min='5:00:00'
#    conda: '../conda_env/samtools.yml'
#    shell:
#        'java -jar {params.picard_dir}/picard.jar CollectRrbsMetrics I={input} M={output.metric} R={params.fasta}' 


# F19K16 and F24B22 methylated filler QC
rule F19K16_F24B22_QC:
    input:
        path_to_data + '/{cohort}/results/bam_markdup/{sample}.aligned.sorted.markdup.bam',
    output:
        methyl_counts = temp(path_to_data + '/{cohort}/qc/methyl_qc/{sample}.methyl_counts'),
        methyl_summary = temp(path_to_data + '/{cohort}/qc/methyl_qc/{sample}.methyl_summary.txt'),
    resources: cpus=1, mem_mb=16000, time_min='2:00:00'
    conda: '../conda_env/samtools.yml'
    shell:
        'bash src/QC/methyl_QC.sh {input} {output.methyl_counts} {output.methyl_summary}'


# Unified QC rule which runs all of the above QCs
rule QC_out:
    input:
        #fastqc_input_fastq: had to be input to trim_galore in process_fastq.smk for wildcards to populate
        fastqc_bam = path_to_data + '/{cohort}/qc/fastqc_bam/{sample}.aligned.sorted.markdup_hg38only_fastqc.html',
        flagstat = path_to_data + '/{cohort}/qc/flagstat/{sample}.aligned.sorted.markdup_hg38only.bam.flagstat',
        medips_qc = path_to_data + '/{cohort}/qc/medips_qc/{sample}_QC_MEDIPS.csv',
        picard_gcbias_metric = path_to_data + '/{cohort}/qc/picard_qc/{sample}.gc_bias_metrics.txt',
        picard_gcbias_chart = path_to_data + '/{cohort}/qc/picard_qc/{sample}.gc_bias_metrics.pdf',
        picard_gcbias_summary = path_to_data + '/{cohort}/qc/picard_qc/{sample}.gc_bias_summary_metrics.txt',
        picard_insertsize_metric = path_to_data + '/{cohort}/qc/picard_qc/{sample}.picardInsertSize_metrics.txt',
        picard_markdup_metric = path_to_data + '/{cohort}/qc/picard_qc/{sample}.marked_dup_metrics.txt',
        picard_alignment_metric = path_to_data + '/{cohort}/qc/picard_qc/{sample}.alignmentMetrics.txt',
        picard_quality_metric = path_to_data + '/{cohort}/qc/picard_qc/{sample}.quality_yield_metrics.txt',
        #picard_rrbs_metric = path_to_data + '/{cohort}/qc/picard_qc/{sample}.Rrbs',
        #picard_rrbs_detailed_metric = path_to_data + '/{cohort}/qc/picard_qc/{sample}.Rrbs.rrbs_detail_metrics',
        #picard_rrbs_qc_pdf = path_to_data + '/{cohort}/qc/picard_qc/{sample}.Rrbs.rrbs_qc.pdf',
        #picard_rrbs_summary = path_to_data + '/{cohort}/qc/picard_qc/{sample}.Rrbs.rrbs_summary_metrics',
        F19K16_F24B22_methyl_summary = path_to_data + '/{cohort}/qc/methyl_qc/{sample}.methyl_summary.txt',
    output:
        full_qc = path_to_data + '/{cohort}/qc/{sample}_qc_full.txt',
    resources: cpus=1, mem_mb=1000, time_min='00:01:00'
    conda: '../conda_env/samtools.yml'
    shell:
        'cat {input.medips_qc} {input.F19K16_F24B22_methyl_summary}.parsed {input.picard_gcbias_summary}.parsed {input.picard_insertsize_metric}.parsed {input.picard_markdup_metric}.parsed {input.picard_alignment_metric}.parsed {input.picard_quality_metric}.parsed > {output.full_qc} && rm {input.F19K16_F24B22_methyl_summary}.parsed {input.picard_gcbias_summary}.parsed {input.picard_insertsize_metric}.parsed {input.picard_markdup_metric}.parsed {input.picard_alignment_metric}.parsed {input.picard_quality_metric}.parsed'

