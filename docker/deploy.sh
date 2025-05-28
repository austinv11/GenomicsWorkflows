#!/usr/bin/env sh

# If docker1 is an executable, we are in biohpc
IS_BIOHPC=$(command -v docker1)
username=$(whoami)

if [ -z "$IS_BIOHPC" ]; then
    docker_exec="docker"
    docker_user=$username
else # biohpc
    docker_exec="docker1"
    docker_user="biohpc_$username"
fi

# Delete old sif files
rm -f $(pwd)/*/*.sif
$docker_exec build -t 10xrangers "$(pwd)/10xrangers/"
$docker_exec save -o "$(pwd)/10xrangers/10x_ranger_images.tar" "$docker_user/10xrangers:latest"
apptainer build "$(pwd)/10xrangers/10x_ranger_image.sif" "docker-archive://$(pwd)/10xrangers/10x_ranger_images.tar"
rm -f "$(pwd)/10xrangers/10x_ranger_images.tar"
$docker_exec build -t multiqc "$(pwd)/multiqc/"
$docker_exec save -o "$(pwd)/multiqc/multiqc_images.tar" "$docker_user/multiqc:latest"
apptainer build "$(pwd)/multiqc/multiqc_image.sif" "docker-archive://$(pwd)/multiqc/multiqc_images.tar"
rm -f "$(pwd)/multiqc/multiqc_images.tar"
$docker_exec build -t samtools "$(pwd)/samtools/"
$docker_exec save -o "$(pwd)/samtools/samtools_images.tar" "$docker_user/samtools:latest"
apptainer build "$(pwd)/samtools/samtools_image.sif" "docker-archive://$(pwd)/samtools/samtools_images.tar"
rm -f "$(pwd)/samtools/samtools_images.tar"

# Generate conda environments on-the-fly using conda-lock
rm -rf "$(pwd)/generated"
mkdir -p "$(pwd)/generated"
conda-lock -f "$(pwd)/../conda/scanpy-environment.yml" \
  --platform linux-64 \
  --kind explicit \
  --filename-template "$(pwd)/generated/scanpy-environment.lock"

# Generate the Dockerfile for the scanpy environment
echo "FROM continuumio/miniconda3:latest AS builder" > "$(pwd)/generated/Dockerfile"
echo "COPY scanpy-environment.lock /tmp/" >> "$(pwd)/generated/Dockerfile"
echo "RUN conda create -p /opt/env --copy --file /tmp/scanpy-environment.lock" >> "$(pwd)/generated/Dockerfile"
echo "FROM gcr.io/distroless/base-debian10" >> "$(pwd)/generated/Dockerfile"
echo "COPY --from=builder /opt/env /opt/env" >> "$(pwd)/generated/Dockerfile"
echo "ENV PATH=/opt/env/bin:\$PATH" >> "$(pwd)/generated/Dockerfile"
# Build the scanpy environment Docker image
$docker_exec build -t scanpy-environment "$(pwd)/generated/"
$docker_exec save -o "$(pwd)/generated/scanpy-environment.tar" "$docker_user/scanpy-environment:latest"
# Convert the Docker image to a Singularity image file
apptainer build "$(pwd)/generated/scanpy-environment.sif" "docker-archive://$(pwd)/generated/scanpy-environment.tar"
rm -f "$(pwd)/generated/scanpy-environment.tar"
# Clean up the generated directory
rm -f "$(pwd)/generated/Dockerfile"
