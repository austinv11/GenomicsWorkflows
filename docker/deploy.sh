#!/usr/bin/env sh

# If docker1 is an executable, we are in biohpc
IS_BIOHPC=$(command -v docker1)

if [ -z "$IS_BIOHPC" ]; then
    docker_exec="docker"
else # biohpc
    docker_exec="docker1"
fi

$docker_exec build -t 10xrangers 10xrangers/
$docker_exec build -t multiqc multiqc/
$docker_exec build -t samtools samtools/
