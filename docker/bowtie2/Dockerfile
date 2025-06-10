#FROM alpine:latest
FROM debian:bookworm-slim
LABEL authors="austinv11"

ENV BOWTIE2_VERSION="2.5.4"

ENV BOWTIE2_URL="https://github.com/BenLangmead/bowtie2/releases/download/v${BOWTIE2_VERSION}/bowtie2-${BOWTIE2_VERSION}-linux-x86_64.zip"

WORKDIR /build
ADD $BOWTIE2_URL /build/

RUN apt-get update && \
    apt-get install -y unzip && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* && \
    unzip bowtie2-${BOWTIE2_VERSION}-linux-x86_64.zip && \
    rm bowtie2-${BOWTIE2_VERSION}-linux-x86_64.zip && \
    chmod +x bowtie2-${BOWTIE2_VERSION}-linux-x86_64/bowtie2* && \

ENV PATH="/build/bowtie2-${BOWTIE2_VERSION}-linux-x86_64/:$PATH"

WORKDIR /work

VOLUME /work

CMD []
ENTRYPOINT ["bash", "-c", "$@"]
