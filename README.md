# PLBR_MedRemixBEDPE

## Summary

- Run MedRemix methylation profiler with bedpe input as part of PLBR database workflow.
- Pipeline is an extension of the original MedRemix pipeline (https://github.com/pughlab/cfMeDIP-seq-analysis-pipeline)
- Pipeline can take FASTQs, BAM or BEDPE as input.
- Pipeline is written in Snakemake, designed to run on SLURM cluster, but can run locally as well.
- Note: 
  - .bam filtering is slightly different from original MedRemix pipeline, therefore results will differ slightly
  - however this pipeline will output same results whether using .bam or .bedpe as input.

## Workflow overview

![MedRemixBEDPE](https://user-images.githubusercontent.com/98410560/168604501-ce838357-b058-4df8-b6f7-cd4c1293e22d.png)

Please see [wiki](https://github.com/mhanbioinfo/PLBR_MedRemixBEDPE/wiki) for tutorial on settings up and running pipeline.
