# -------------------- #
# convert bam to bedpe #
# -------------------- #

rule bam2bedpe:
    input:
        path_to_data + '/{cohort}/results/bam_markdup/{sample}.aligned.sorted.markdup.bam',
    output:
        bedpe = path_to_data + '/{cohort}/results/bedpe_out/{sample}.aligned.sorted.markdup_coordSortd.bedpe.gz',
    params:
        slurm_local = config['pipeline_params']['slurm_local'],
        chunks = config['pipeline_params']['chunks'],
        out_dir = path_to_data + '/{cohort}/results/bedpe_out',
        src_dir = config['pipeline_params']['src_dir'],
        picard_dir = config['pipeline_params']['picard_dir'],
        conda_activate = config['pipeline_params']['conda_activate'],
        conda_env = config['pipeline_params']['conda_env'],
    resources: cpus=1, mem_mb=60000, time_min='3-00:00:00'
    conda: '../conda_env/samtools.yml'
    shell:
        'bash src/bam2bedpe_scripts/bam2bedpe.sh -s {params.slurm_local} -c {params.chunks} -b {input} -o {params.out_dir} -x {params.src_dir} -p {params.picard_dir} -a {params.conda_activate} -e {params.conda_env}'

rule get_clean_bedpe4medremix:
    input:
        bedpe = path_to_data + '/{cohort}/results/bedpe_out/{sample}.aligned.sorted.markdup_coordSortd.bedpe.gz',
    output:
        bedpe4medremix = temp(path_to_data + '/{cohort}/results/bedpe4medremix_out/{sample}_{species}_{chrom}.noAmbiguous.no113n177.mapQge20.isMate.lean.withEnd.lenFiltd.bedpe4medremix'),
    params:
        frag_len_limit = 500,
        out_dir = path_to_data + '/{cohort}/results/bedpe4medremix_out',
    resources: cpus=4, mem_mb=30000, time_min='1-00:00:00'
    shell:
        'bash src/bam2bedpe_scripts/get_clean_bedpe4medremix_v4_getopts.sh -i {input} -s {wildcards.sample} -p {wildcards.species} -c {wildcards.chrom} -f {params.frag_len_limit} -o {params.out_dir}'
