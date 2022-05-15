# -------------------------- #
#  Pre-process input BAMs  #
# -------------------------- #

def get_bam_path(cohort, sample):
    """Retrieves the path to the bam file.

    Keyword arguments:
        cohort -- name of the cohort whose samplesheet should be accessed.
        sample -- identifier of the sample as specified in the samplesheet.
    """
    cohort_data = get_cohort_data(cohort)
    sample_line = cohort_data[
        (cohort_data.sample_name == sample)
    ]
    return(sample_line.bam_path.to_list()[0])

rule softlink_bam:
    input:
        lambda wildcards: get_bam_path(wildcards.cohort, wildcards.sample),
    output:
        #temp(path_to_data + '/{cohort}/results/bam_markdup/{sample}.aligned.sorted.markdup.bam')
        path_to_data + '/{cohort}/results/bam_markdup/{sample}.aligned.sorted.markdup.bam'
    resources: cpus=1, mem_mb=1000, time_min='2:00:00'
    conda: '../conda_env/samtools.yml'
    shell:
        'ln -s {input} {output} && samtools index {output}'

