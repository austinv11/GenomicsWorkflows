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
