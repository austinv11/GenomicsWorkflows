FROM ubuntu:latest
LABEL authors="austinv11"

# Disable 10X Telemetry
ENV TENX_DISABLE_TELEMETRY=1

WORKDIR /build

# Install tree command
RUN apt-get update \
    && apt-get install -y tree curl \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Set the URLS
## CellRanger
ENV CELLRANGER_DOWNLOAD_VER="9.0.1"
ENV CELLRANGER_HUMAN_REFERENCE_URL="https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2024-A.tar.gz"
ENV CELLRANGER_MOUSE_REFERENCE_URL="https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCm39-2024-A.tar.gz"
#ENV CELLRANGER_RAT_REFERENCE_URL="https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-mRatBN7-2-2024-A.tar.gz"
ENV CELLRANGER_HUMAN_AND_MOUSE_REFERENCE_URL="https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38_and_GRCm39-2024-A.tar.gz"
#ENV CELLRANGER_HUMAN_VDJ_REFERENCE_URL="https://cf.10xgenomics.com/supp/cell-vdj/refdata-cellranger-vdj-GRCh38-alts-ensembl-7.1.0.tar.gz"
#ENV CELLRANGER_MOUSE_VDJ_REFERENCE_URL="https://cf.10xgenomics.com/supp/cell-vdj/refdata-cellranger-vdj-GRCm38-alts-ensembl-7.0.0.tar.gz"
#probe/antibody probe barcodes: https://cf.10xgenomics.com/supp/cell-exp/probeset/probe-barcodes-fixed-rna-profiling.txt
ENV CELLRANGER_FLEX_HUMAN_PROBE_SET_URL="https://cf.10xgenomics.com/supp/cell-exp/probeset/Chromium_Human_Transcriptome_Probe_Set_v1.1.0_GRCh38-2024-A.csv"
ENV CELLRANGER_FLEX_MOUSE_PROBE_SET_URL="https://cf.10xgenomics.com/supp/cell-exp/probeset/Chromium_Mouse_Transcriptome_Probe_Set_v1.1.1_GRCm39-2024-A.csv"
## Spaceranger
ENV SPACERANGER_DOWNLOAD_VER="4.0.1"
ENV SPACERANGER_HUMAN_REFERENCE_URL="https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2024-A.tar.gz"
ENV SPACERANGER_MOUSE_REFERENCE_URL="https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCm39-2024-A.tar.gz"
ENV SPACERANGER_HUMAN_PROBE_SET_URL="https://cf.10xgenomics.com/supp/spatial-exp/probeset/Visium_Human_Transcriptome_Probe_Set_v2.1.0_GRCh38-2024-A.csv"
ENV SPACERANGER_MOUSE_PROBE_SET_URL="https://cf.10xgenomics.com/supp/spatial-exp/probeset/Visium_Mouse_Transcriptome_Probe_Set_v2.1.0_GRCm39-2024-A.csv"
## Xeniumranger
ENV XENIUMRANGER_DOWNLOAD_VER="3.1.1"
ENV XENIUMRANGER_DOWNLOAD_VER_TRUNCATED="3.1"

# Copy rangers...
ADD "cellranger-${CELLRANGER_DOWNLOAD_VER}.tar.gz" /bin/
ADD "spaceranger-${SPACERANGER_DOWNLOAD_VER}.tar.gz" /bin/
ADD "xeniumranger-${XENIUMRANGER_DOWNLOAD_VER}.tar.gz" /bin/
ENV CELLRANGER_BIN_PATH="/bin/cellranger-${CELLRANGER_DOWNLOAD_VER}/bin/"
ENV SPACERANGER_BIN_PATH="/bin/spaceranger-${SPACERANGER_DOWNLOAD_VER}/bin/"
ENV XENIUMRANGER_BIN_PATH="/bin/xeniumranger-xenium${XENIUMRANGER_DOWNLOAD_VER_TRUNCATED}/bin/"
ENV PATH="$PATH:${CELLRANGER_BIN_PATH}:${SPACERANGER_BIN_PATH}:${XENIUMRANGER_BIN_PATH}"

# Wrapper script to get the reference data
ADD reference_functions.sh /bin
# Ensure the bash script is run to expose the functions
RUN chmod a+x /bin/reference_functions.sh && \
    . /bin/reference_functions.sh

## Download references
#ADD $CELLRANGER_HUMAN_REFERENCE_URL /ref/cellranger/
#ADD $CELLRANGER_MOUSE_REFERENCE_URL /ref/cellranger/
##ADD $CELLRANGER_RAT_REFERENCE_URL /ref/cellranger/
#ADD $CELLRANGER_HUMAN_AND_MOUSE_REFERENCE_URL /ref/cellranger/
##ADD $CELLRANGER_HUMAN_VDJ_REFERENCE_URL /ref/cellranger/
##ADD $CELLRANGER_MOUSE_VDJ_REFERENCE_URL /ref/cellranger/
#ADD $CELLRANGER_FLEX_HUMAN_PROBE_SET_URL /ref/cellranger/
#ADD $CELLRANGER_FLEX_MOUSE_PROBE_SET_URL /ref/cellranger/
#
#ADD $SPACERANGER_HUMAN_REFERENCE_URL /ref/spaceranger/
#ADD $SPACERANGER_MOUSE_REFERENCE_URL /ref/spaceranger/
#ADD $SPACERANGER_HUMAN_PROBE_SET_URL /ref/spaceranger/
#ADD $SPACERANGER_MOUSE_PROBE_SET_URL /ref/spaceranger/
#
## Extract the references and delete the tarballs but iterating over all tar.gz files
#RUN for file in /ref/cellranger/*.tar.gz; do \
#        tar -xzf "$file" -C /ref/cellranger/ && \
#        rm -f "$file"; \
#    done && \
#    for file in /ref/spaceranger/*.tar.gz; do \
#        tar -xzf "$file" -C /ref/spaceranger/ && \
#        rm -f "$file"; \
#    done

WORKDIR /work
RUN rm -rf /build

COPY wrapper.sh /bin/

RUN cellranger telemetry disable all \
    && spaceranger telemetry disable all
    #\
    #&& xeniumranger telemetry disable all

# Final mounting
VOLUME /work

CMD []
ENTRYPOINT ["/bin/wrapper.sh"]
