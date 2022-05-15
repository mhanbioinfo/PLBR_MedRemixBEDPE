# ---------------------- #
#  Compute BAM bin stats #
# ---------------------- #

# Generate bin stats for each chromosome individually.
rule bam_bin_stats:
    input:
        path_to_data + '/{cohort}/results/bam_markdup/{sample}.aligned.sorted.markdup.bam',
    output:
        #binstat = temp(path_to_data + '/{cohort}/tmp/bam_bin_stats/bam_bin_stats_{sample}_{species}_{chrom}.tsv'),
        binstat = path_to_data + '/{cohort}/tmp/bam_bin_stats/bam_bin_stats_{sample}_{species}_{chrom}.tsv',
        filtered = path_to_data + '/{cohort}/results/bam_bin_stats_filtered_out/bam_removed_bins_{sample}_{species}_{chrom}.tsv',
    params:
        bsgenome = lambda wildcards: get_cohort_config(wildcards.cohort)['bsgenome'][wildcards.species],
        winsize = config['pipeline_params']['window_size'],
    resources: cpus=1, mem_mb=30000, time_min='24:00:00'
    conda: '../conda_env/cfmedip_r.yml'
    shell:
        clean('''
        Rscript src/R/bin_stats_bam.R
            -i {input}
            -g {params.bsgenome}
            -c {wildcards.chrom}
            -o {output.binstat}
            --filtered {output.filtered}
            --winsize {params.winsize}
        ''')

# Preload chromosome values for below rule
chromosomes = {}
for c in config['data']['cohorts']:
    if config['data']['cohorts'][c]['active']:
        chr_data = get_cohort_config(c)['chromosomes']
        chromosome_tuples = [(species, chrom) for species in chr_data for chrom in chr_data[species].split(',')]
        chromosomes[c] = chromosome_tuples

# Merge bin stats across all chromosomes.
rule bam_merged_bin_stats:
    input:
        lambda wildcards: [path_to_data + '/{{cohort}}/tmp/bam_bin_stats/bam_bin_stats_{{sample}}_{species}_{chrom}.tsv'.format(species=a[0], chrom=a[1]) for a in chromosomes[wildcards.cohort]],
    output:
        path_to_data + '/{cohort}/results/bam_merged_bin_stats/bam_bin_stats_{sample}.feather'
    params:
        paths = lambda wildcards, input: ','.join(input)
    resources: cpus=1, mem_mb=8000, time_min='24:00:00'
    conda: '../conda_env/cfmedip_r.yml'
    shell:
        'Rscript src/R/row_bind_tables.R -p "{params.paths}" -o {output} --in-tsv --out-feather --omit-paths'

# ------------------------ #
#  Compute BEDPE bin stats #
# ------------------------ #

rule bedpe_bin_stats:
    input:
        bedpe4medremix = path_to_data + '/{cohort}/results/bedpe4medremix_out/{sample}_{species}_{chrom}.noAmbiguous.no113n177.mapQge20.isMate.lean.withEnd.lenFiltd.bedpe4medremix',
    output:
        #bedpe_binstat = temp(path_to_data + '/{cohort}/tmp/bedpe_bin_stats/bedpe_bin_stats_{sample}_{species}_{chrom}.tsv'),
        bedpe_binstat = path_to_data + '/{cohort}/tmp/bedpe_bin_stats/bedpe_bin_stats_{sample}_{species}_{chrom}.tsv',
    params:
        bsgenome = lambda wildcards: get_cohort_config(wildcards.cohort)['bsgenome'][wildcards.species],
        winsize = config['pipeline_params']['window_size'],
    resources: cpus=1, mem_mb=30000, time_min='24:00:00'
    conda: '../conda_env/cfmedip_r.yml'
    shell:
        clean('''
        Rscript src/R/bin_stats_bedpe.R
            -i {input}
            -g {params.bsgenome}
            -c {wildcards.chrom}
            -o {output.bedpe_binstat}
            --winsize {params.winsize}
        ''')

# Merge bedpe4medremix across all chromosomes.
rule bedpe_merged_bin_stats:
    input:
        lambda wildcards: [path_to_data + '/{{cohort}}/tmp/bedpe_bin_stats/bedpe_bin_stats_{{sample}}_{species}_{chrom}.tsv'.format(species=a[0], chrom=a[1]) for a in chromosomes[wildcards.cohort]],
    output:
        path_to_data + '/{cohort}/results/bedpe_merged_bin_stats/bedpe_bin_stats_{sample}.feather'
    params:
        paths = lambda wildcards, input: ','.join(input)
    resources: cpus=1, mem_mb=8000, time_min='24:00:00'
    conda: '../conda_env/cfmedip_r.yml'
    shell:
        'Rscript src/R/row_bind_tables.R -p "{params.paths}" -o {output} --in-tsv --out-feather --omit-paths'
