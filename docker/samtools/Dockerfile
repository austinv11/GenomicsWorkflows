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
    libdeflate-dev \
    g++ \
    python3 \
    bash

ENV SAMTOOLS_VER="1.21"
ENV BCFTOOLS_VER="1.21"
#ENV BEDTOOLS_VER="2.31.1"

ENV SAMTOOLS_URL="https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VER}/samtools-${SAMTOOLS_VER}.tar.bz2"
ENV BCFTOOLS_URL="https://github.com/samtools/bcftools/releases/download/${BCFTOOLS_VER}/bcftools-${BCFTOOLS_VER}.tar.bz2"
#ENV BEDTOOLS_URL="https://github.com/arq5x/bedtools2/releases/download/v${BEDTOOLS_VER}/bedtools-${BEDTOOLS_VER}.tar.gz"

WORKDIR /build
ADD $SAMTOOLS_URL /build/
ADD $BCFTOOLS_URL /build/
#ADD $BEDTOOLS_URL /build/

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

#RUN tar xzf bedtools-${BEDTOOLS_VER}.tar.gz && \
#    cd bedtools2*/ \
#    && make \
#    && mv bin/* /bin/


ENV PATH="/bin/:$PATH"

WORKDIR /work

VOLUME /work

CMD []
ENTRYPOINT ["bash", "-c", "$@"]
