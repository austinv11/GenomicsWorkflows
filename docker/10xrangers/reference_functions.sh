#!/usr/bin/env bash

# Function to download from a URL and save to a file
# If the file ends in tar.gz, it will be extracted
download_and_extract() {
    local url="$1"
    local output_dir=$(pwd)  #"$2"

    local filename=$(basename "$url")
    local filepath="$output_dir/$filename"
    local extract_dir="$output_dir/${filename%.tar.gz}"

    # Create output directory if it doesn't exist
    mkdir -p "$output_dir"
    # Download the file
    if ! curl -L -o "$filepath" "$url"; then
        echo "Error downloading $url"
        return 1
    fi
    echo "Downloaded $filename to $output_dir"
    # Check if the file is a tar.gz file
    if [[ "$filename" == *.tar.gz ]]; then
        # Extract the tar.gz file
        if ! tar -xzf "$filepath" -C "$output_dir"; then
            echo "Error extracting $filepath"
            return 1
        fi
        echo "Extracted $filename to $extract_dir"
        # Return the path to the extracted directory
        return "$extract_dir"
    else
        # Return the path to the downloaded file
        return "$extract_dir"
    fi
}


# Download and extract the 10x Genomics reference data, returning the final path
cellranger_human_reference() {
    local path=$(download_and_extract $CELLRANGER_HUMAN_REFERENCE_URL)
    if [ $? -ne 0 ]; then
        echo "Failed to download or extract the reference data."
        return 1
    fi
    echo "Reference data downloaded and extracted to: $path"
    return "$path"
}

cellranger_mouse_reference() {
    local path=$(download_and_extract $CELLRANGER_MOUSE_REFERENCE_URL)
    if [ $? -ne 0 ]; then
        echo "Failed to download or extract the reference data."
        return 1
    fi
    echo "Reference data downloaded and extracted to: $path"
    return "$path"
}

cellranger_human_and_mouse_reference() {
    local path=$(download_and_extract $CELLRANGER_HUMAN_AND_MOUSE_REFERENCE_URL)
    if [ $? -ne 0 ]; then
        echo "Failed to download or extract the reference data."
        return 1
    fi
    echo "Reference data downloaded and extracted to: $path"
    return "$path"
}

cellranger_human_flex_probe_set() {
    local path=$(download_and_extract $CELLRANGER_HUMAN_FLEX_PROBE_SET_URL)
    if [ $? -ne 0 ]; then
        echo "Failed to download or extract the flex probe set."
        return 1
    fi
    echo "Flex probe set downloaded and extracted to: $path"
    return "$path"
}

cellranger_mouse_flex_probe_set() {
    local path=$(download_and_extract $CELLRANGER_MOUSE_FLEX_PROBE_SET_URL)
    if [ $? -ne 0 ]; then
        echo "Failed to download or extract the flex probe set."
        return 1
    fi
    echo "Flex probe set downloaded and extracted to: $path"
    return "$path"
}

spaceranger_human_reference() {
    local path=$(download_and_extract $SPACERANGER_HUMAN_REFERENCE_URL)
    if [ $? -ne 0 ]; then
        echo "Failed to download or extract the reference data."
        return 1
    fi
    echo "Reference data downloaded and extracted to: $path"
    return "$path"
}

spaceranger_mouse_reference() {
    local path=$(download_and_extract $SPACERANGER_MOUSE_REFERENCE_URL)
    if [ $? -ne 0 ]; then
        echo "Failed to download or extract the reference data."
        return 1
    fi
    echo "Reference data downloaded and extracted to: $path"
    return "$path"
}

spaceranger_human_probe_set() {
    local path=$(download_and_extract $SPACERANGER_HUMAN_PROBE_SET_URL)
    if [ $? -ne 0 ]; then
        echo "Failed to download or extract the probe set."
        return 1
    fi
    echo "Probe set downloaded and extracted to: $path"
    return "$path"
}

spaceranger_mouse_probe_set() {
    local path=$(download_and_extract $SPACERANGER_MOUSE_PROBE_SET_URL)
    if [ $? -ne 0 ]; then
        echo "Failed to download or extract the probe set."
        return 1
    fi
    echo "Probe set downloaded and extracted to: $path"
    return "$path"
}
