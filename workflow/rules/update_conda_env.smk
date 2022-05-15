# ------------------------------------------ #
#  Install conda dependencies automatically  #
# ------------------------------------------ #
# This section is only necessary to support air-gapped environments,
# where the cluster nodes have no direct access to internet.
# In this setting, you must first load all necessary packages from
# a build node with internet access. This can be done simply
# by running `bash update_conda.sh`
# from a node with internet access, which loads install_dependencies
# rule below to install all of the necessary conda dependencies
# in advance. Any time the git repo is pulled or if you modify
# any of the conda environments in ../conda_env, this needs to be
# re-run.

rule install_dependencies:
    input:
        'conda_env/biopython_installed',
        'conda_env/cfmedip_r_installed',
        'conda_env/fastqc_installed',
        'conda_env/samtools_installed',
        'conda_env/trimgalore_installed'
    output:
        'conda_env/dependencies'
    shell:
        'touch {output}'

rule biopython:
    output:
        'conda_env/biopython_installed'
    conda: '../conda_env/biopython.yml'
    shell:
        'touch {output}'

rule cfmedip_r:
    output:
        'conda_env/cfmedip_r_installed'
    conda: '../conda_env/cfmedip_r.yml'
    shell:
        'touch {output}'

rule fastqc:
    output:
        'conda_env/fastqc_installed'
    conda: '../conda_env/fastqc.yml'
    shell:
        'touch {output}'

rule samtools:
    output:
        'conda_env/samtools_installed'
    conda: '../conda_env/samtools.yml'
    shell:
        'touch {output}'

rule trimgalore:
    output:
        'conda_env/trimgalore_installed'
    conda: '../conda_env/trimgalore.yml'
    shell:
        'touch {output}'

