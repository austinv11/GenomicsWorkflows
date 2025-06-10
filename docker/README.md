# Docker environments
This directory contains docker images that can be used for various genomics analyses.

Additionally, it can generate docker images from the conda environments in the `conda` directory. 

To build all these images, run `deploy.sh` script.

## Available Images
- `10xrangers`: Contains installations for cellranger, spaceranger, and xeniumranger.
- `multiqc`: Fastp, multiqc for quality control of sequencing data.
- `samtools`: Contains samtools, bcftools, and bedtools for manipulating sequencing data.
- `bowtie2`: Contains bowtie2 for alignment of sequencing reads.
- `bwa`: Contains bwa-mem2 for alignment of sequencing reads.
- `minimap2`: Contains minimap2 for alignment of sequencing reads.
- `gatk`: Contains GATK for variant calling and other genomic analyses.
- `vep`: Contains VEP for annotating variants.
