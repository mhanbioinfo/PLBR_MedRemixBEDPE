# ----------------------------------------------------------------------------- #
#  Fit cfMeDIP-seq coverage stats to infer absolute methylation using MedReMix  #
# ----------------------------------------------------------------------------- #

bam_or_bedpe = config['pipeline_params']['bam_or_bedpe']

def get_cfmedip_nbglm_input(bam_or_bedpe):
    if bam_or_bedpe == 'bam':
        return(path_to_data + '/{cohort}/results/bam_merged_bin_stats/bam_bin_stats_{sample}.feather')
    elif bam_or_bedpe == 'bedpe': 
        return(path_to_data + '/{cohort}/results/bedpe_merged_bin_stats/bedpe_bin_stats_{sample}.feather')

def get_cfmedip_nbglm_output(bam_or_bedpe):
    if bam_or_bedpe == 'bam':
        return([path_to_data + '/{cohort}/results/bam_cfmedip_nbglm/bam_{sample}_fit_nbglm.tsv',
                path_to_data + '/{cohort}/results/bam_cfmedip_nbglm/bam_{sample}_fit_nbglm_model.Rds'])
    elif bam_or_bedpe == 'bedpe':
        return([path_to_data + '/{cohort}/results/bedpe_cfmedip_nbglm/bedpe_{sample}_fit_nbglm.tsv',
                path_to_data + '/{cohort}/results/bedpe_cfmedip_nbglm/bedpe_{sample}_fit_nbglm_model.Rds'])

# This is the currently used negative binomial GLM approach to fitting
rule cfmedip_nbglm:
    input:
        get_cfmedip_nbglm_input(bam_or_bedpe)
    output:
        #fit=path_to_data + '/{cohort}/results/bedpe_cfmedip_nbglm/bedpe_{sample}_fit_nbglm.feather',
        fit=get_cfmedip_nbglm_output(bam_or_bedpe)[0],
        model=get_cfmedip_nbglm_output(bam_or_bedpe)[1]
    resources: cpus=1, time_min='1-00:00:00', mem_mb=lambda wildcards, attempt: 16000 if attempt == 1 else 30000
    conda: '../conda_env/cfmedip_r.yml'
    shell:
        #'Rscript src/R/cfmedip_nbglm.R -i {input} -o {output.fit}'
        'Rscript src/R/cfmedip_nbglm.R -i {input} -o {output.fit} --modelout {output.model}'
