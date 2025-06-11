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

build_docker_image() {
    local image_name="$1"
    $docker_exec build -t "$image_name" "$(pwd)/$image_name/"
    $docker_exec save "$docker_user/$image_name:latest" > "$(pwd)/$image_name/$image_name.tar"
    apptainer build "$(pwd)/$image_name/$image_name.sif" "docker-archive://$(pwd)/$image_name/$image_name.tar"
    rm -f "$(pwd)/$image_name/$image_name.tar"
}

build_docker_image "10xrangers"
build_docker_image "multiqc"
build_docker_image "samtools"
build_docker_image "bowtie2"
build_docker_image "bwa"
build_docker_image "minimap2"
build_docker_image "gatk"
build_docker_image "vep"

# Generate conda environments on-the-fly using conda-lock
generate_conda_docker_images() {
  local conda_dir="${1:-$(pwd)/../conda}"
  local generated_dir="$(pwd)/generated"
  rm -rf "$generated_dir"
  mkdir -p "$generated_dir"

  for env_file in "$conda_dir"/*.yml; do
    [ -e "$env_file" ] || continue
    env_name=$(basename "$env_file")
    env_name="${env_name%.yml}" # Remove the .yml extension
    lock_file="$generated_dir/${env_name}.lock"
    dockerfile="$generated_dir/Dockerfile"

    conda-lock -f "$env_file" \
      --platform linux-64 \
      --kind explicit \
      --filename-template "$lock_file"

    echo "FROM continuumio/miniconda3:latest AS builder" > "$dockerfile"
    echo "COPY ${env_name}.lock /tmp/" >> "$dockerfile"
    echo "ENV BIOC_MIRROR=https://bioconductor.posit.com" >> "$dockerfile"
    echo "RUN conda create -p /opt/env --copy --file /tmp/${env_name}.lock" >> "$dockerfile"
    # FIXME: Temporary fix for some packages not being updated in conda
    # If scdblfinder in the environment name, manually install with biocmanager in R
    if echo "$env_name" | grep -q "scdblfinder"; then
      # Install the development version of plger/scDblFinder
      # Overwrite domain resolution for Bioconductor packages in the hosts file
      echo "RUN conda run -p /opt/env R -e \"options(BioC_mirror=Sys.getenv('BIOC_MIRROR')); BiocManager::install('plger/scDblFinder', force = TRUE, ask = FALSE, update = FALSE)\"" >> "$dockerfile"
    fi
    echo "FROM debian:bookworm-slim" >> "$dockerfile"
    echo "COPY --from=builder /opt/env /opt/env" >> "$dockerfile"
    echo "ENV PATH=/opt/env/bin:\$PATH" >> "$dockerfile"

    $docker_exec build -t "${env_name}" "$generated_dir/"
    $docker_exec save "$docker_user/${env_name}:latest" > "$generated_dir/${env_name}.tar"
    apptainer build "$generated_dir/${env_name}.sif" "docker-archive://$generated_dir/${env_name}.tar"
    rm -f "$generated_dir/${env_name}.tar" "$dockerfile" "$lock_file"
  done
}

generate_conda_docker_images "$(pwd)/../conda"
