#FROM alpine:latest
FROM debian:bookworm-slim
LABEL authors="austinv11"

ENV MINIMAP2_VERSION="2.29"

ENV MINIMAP2_URL="https://github.com/lh3/minimap2/releases/download/v${MINIMAP2_VERSION}/minimap2-${MINIMAP2_VERSION}_x64-linux.tar.bz2"

WORKDIR /build
ADD $MINIMAP2_URL /build/

RUN apt-get update && \
    apt-get install -y bzip2 && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* && \
    tar -xvjf minimap2-${MINIMAP2_VERSION}_x64-linux.tar.bz2 && \
    rm minimap2-${MINIMAP2_VERSION}_x64-linux.tar.bz2 && \
    chmod +x minimap2-${MINIMAP2_VERSION}_x64-linux/minimap2

ENV PATH="/build/minimap2-${MINIMAP2_VERSION}_x64-linux/:$PATH"

WORKDIR /work

VOLUME /work

CMD []
ENTRYPOINT ["bash", "-c", "$@"]
