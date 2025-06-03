# 10X Rangers
This image contains installations for cellranger, spaceranger, and xeniumranger.

**Note**: You must have downloaded cellranger-version.tar.gz, spaceranger-version.tar.gz, and xeniumranger-version.tar.gz in the current directory for this to build correctly.


## Usage:
Simply enter a shell to have access to all three tools:
```bash
docker run -it --rm austinv11/10xrangers:latest
```

Additionally, there are convenience scripts for downloading and extracting references that you can directly pass to the
tools:
- `cellranger_human_reference`: Downloads and extracts the human reference for cellranger.
- `spaceranger_human_reference`: Downloads and extracts the human reference for spaceranger.
- `cellranger_mouse_reference`: Downloads and extracts the mouse reference for cellranger.
- `spaceranger_mouse_reference`: Downloads and extracts the mouse reference for spaceranger.
- `cellranger_human_and_mouse_reference`: Downloads and extracts the combined human and mouse reference for cellranger.
- `cellranger_human_flex_probe_set`: Downloads and extracts the human flex probe set for cellranger.
- `cellranger_mouse_flex_probe_set`: Downloads and extracts the mouse flex probe set for spaceranger.
- `spaceranger_human_probe_set`: Downloads and extracts the human probe set for spaceranger.
- `spaceranger_mouse_probe_set`: Downloads and extracts the mouse probe set for spaceranger.

Example of running cellranger with the human reference:
```bash
docker run -it --rm \
  -v /path/to/your/data:/data \
  austinv11/10xrangers:latest \
  cellranger count --id=sample_id --fastqs=/data/fastq_files \
    --transcriptome="$(cellranger_human_reference)" \
    --sample=sample_name
```

## Tool Versions
- cellranger: `9.0.1`
- spaceranger: `3.1.3`
- xeniumranger: `3.1.1`
