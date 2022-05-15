#!/bin/bash -login
#SBATCH -J cfMeDIP_MedRemix_bedpe
#SBATCH -t 7-00:00:00
#SBATCH -p himem
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ming.han@uhn.ca

. env.sh

echo 'Running cfMeDIP_MedRemix_bedpe on H4H cluster'
snakemake -p --profile $HOME/.config/snakemake/slurm --configfile ./config.yml
