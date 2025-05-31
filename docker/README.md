# Docker environments
This directory contains docker images that can be used for various genomics analyses.

Additionally, it can generate docker images from the conda environments in the `conda` directory. 

To build all these images, run `deploy.sh` script.

## Available Images
- `10xrangers`: Contains installations for cellranger, spaceranger, and xeniumranger.
- `multiqc`: Fastp, multiqc for quality control of sequencing data.
- `samtools`: Contains samtools, bcftools, and bedtools for manipulating sequencing data.
