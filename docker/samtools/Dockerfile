FROM alpine:latest
#FROM debian:slim-bookworm
LABEL authors="austinv11"

RUN apk add --no-cache \
    autoconf \
    automake \
    make \
    gcc \
    musl-dev \
    perl \
    bash \
    zlib-dev \
    bzip2-dev \
    xz-dev \
    curl-dev \
    openssl-dev \
    ncurses-dev \
    libdeflate-dev

ENV SAMTOOLS_VER="1.21"
ENV BCFTOOLS_VER="1.21"

ENV SAMTOOLS_URL="https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VER}/samtools-${SAMTOOLS_VER}.tar.bz2"
ENV BCFTOOLS_URL="https://github.com/samtools/bcftools/releases/download/${BCFTOOLS_VER}/bcftools-${BCFTOOLS_VER}.tar.bz2"

WORKDIR /build
ADD $SAMTOOLS_URL /build/
ADD $BCFTOOLS_URL /build/

RUN tar xjf samtools-${SAMTOOLS_VER}.tar.bz2 && \
    cd samtools-${SAMTOOLS_VER}/ \
    && ./configure --prefix=/ \
    && make \
    && make install

RUN tar xjf bcftools-${BCFTOOLS_VER}.tar.bz2 && \
    cd bcftools-${BCFTOOLS_VER}/ \
    && ./configure --prefix=/ \
    && make \
    && make install

ENV PATH="/bin/:$PATH"

WORKDIR /work

VOLUME /work

CMD []
ENTRYPOINT ["/bin/bash", "-c", "$@"]
