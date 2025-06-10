#FROM alpine:latest
FROM adoptopenjdk/openjdk8-openj9:debianslim
LABEL authors="austinv11"

# Install Python
RUN apt-get update && \
    apt-get install -y python3 python3-pip aria2 && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

ENV GATK_VERSION="4.6.2.0"

ENV GATK_URL="https://github.com/broadinstitute/gatk/releases/download/${GATK_VERSION}/gatk-${GATK_VERSION}.zip"

WORKDIR /build
ADD $GATK_URL /build/

RUN apt-get update && \
    apt-get install -y unzip && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* && \
    unzip gatk-${GATK_VERSION}.zip && \
    rm gatk-${GATK_VERSION}.zip && \
    chmod +x gatk-${GATK_VERSION}/gatk

ENV PATH="/build/gatk-${GATK_VERSION}/:$PATH"

WORKDIR /work

VOLUME /work

CMD []
ENTRYPOINT ["bash", "-c", "$@"]
