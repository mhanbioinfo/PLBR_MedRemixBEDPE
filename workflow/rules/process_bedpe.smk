# ---------------------- #
#  Process input BEDPEs  #
# ---------------------- #

def get_bedpe_path(cohort, sample):
    """Retrieves the path to the bedpe.gz file.

    Keyword arguments:
        cohort -- name of the cohort whose samplesheet should be accessed.
        sample -- identifier of the sample as specified in the samplesheet.
    """
    cohort_data = get_cohort_data(cohort)
    sample_line = cohort_data[
        (cohort_data.sample_name == sample)
    ]
    return(sample_line.bedpe_path.to_list()[0])


rule softlink_bedpe:
    input:
        lambda wildcards: get_bedpe_path(wildcards.cohort, wildcards.sample),
    output:
        #temp(path_to_data + '/{cohort}/results/bedpe_out/{sample}.aligned.sorted.markdup_coordSortd.bedpe.gz')
        path_to_data + '/{cohort}/results/bedpe_out/{sample}.aligned.sorted.markdup_coordSortd.bedpe.gz'
    resources: cpus=1, mem_mb=1000, time_min='2:00:00'
    shell:
        'ln -s {input} {output}'
