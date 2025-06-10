#FROM alpine:latest
FROM debian:bookworm-slim
LABEL authors="austinv11"

ENV BWAMEM_VERSION="2.2.1"

ENV BWAMEM_URL="https://github.com/bwa-mem2/bwa-mem2/releases/download/v${BWAMEM_VERSION}/bwa-mem2-${BWAMEM_VERSION}_x64-linux.tar.bz2"

WORKDIR /build
ADD $BWAMEM_URL /build/

RUN apt-get update && \
    apt-get install -y bzip2 && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* && \
    tar -xvjf bwa-mem2-${BWAMEM_VERSION}_x64-linux.tar.bz2 && \
    rm bwa-mem2-${BWAMEM_VERSION}_x64-linux.tar.bz2 && \
    chmod +x bwa-mem2-${BWAMEM_VERSION}_x64-linux/bwa-mem2*

ENV PATH="/build/bwa-mem2-${BWAMEM_VERSION}_x64-linux/:$PATH"

WORKDIR /work

VOLUME /work

CMD []
ENTRYPOINT ["bash", "-c", "$@"]
