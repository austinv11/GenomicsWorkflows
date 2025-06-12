#FROM alpine:latest
FROM debian:bookworm-slim
LABEL authors="austinv11"

# Requires HTSlib >= 1.7 and BOOST >= 1.65
RUN apt-get update && \
    apt-get install -y \
        build-essential libbz2-dev libcurl4-openssl-dev autoconf libssl-dev wget zlib1g-dev liblzma-dev libdeflate-dev && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*


ENV GLIMPSE2_VERSION="f72ff67"

ENV GLIMPSE2_URL="https://github.com/odelaneau/GLIMPSE/archive/${GLIMPSE2_VERSION}.zip"

WORKDIR /build

# From: https://github.com/odelaneau/GLIMPSE/blob/master/Dockerfile

ADD https://archives.boost.io/release/1.78.0/source/boost_1_78_0.tar.gz /build/
RUN tar -xf boost_1_78_0.tar.gz && \
    rm boost_1_78_0.tar.gz && \
    cd boost_1_78_0/ && \
    ./bootstrap.sh --with-libraries=iostreams,program_options,serialization --prefix=../boost && \
    ./b2 install && \
    cd .. && \
    cp boost/lib/libboost_iostreams.a boost/lib/libboost_program_options.a boost/lib/libboost_serialization.a /usr/local/lib/ && \
    cp -r boost/include/boost/ /usr/include/ && \
    rm -r boost_1_78_0 boost


ADD https://github.com/samtools/htslib/releases/download/1.16/htslib-1.16.tar.bz2 /build/

RUN tar -xf htslib-1.16.tar.bz2 && \
    rm htslib-1.16.tar.bz2 && \
    cd htslib-1.16 && \
    autoheader && \
    autoconf && \
    ./configure --enable-libcurl && \
    make install && \
    cd .. && \
    rm -r htslib-1.16

ADD $GLIMPSE2_URL /build/

RUN unzip GLIMPSE*.zip && \
    rm GLIMPSE*.zip && \
    mv GLIMPSE* GLIMPSE

WORKDIR /build/GLIMPSE
RUN make clean && make COMPILATION_ENV=docker
RUN mv GLIMPSE/chunk/bin/GLIMPSE2_chunk GLIMPSE/split_reference/bin/GLIMPSE2_split_reference GLIMPSE/phase/bin/GLIMPSE2_phase GLIMPSE/ligate/bin/GLIMPSE2_ligate GLIMPSE/concordance/bin/GLIMPSE2_concordance /bin && \
    chmod +x /bin/GLIMPSE2* && \
    rm -rf GLIMPSE

ENV PATH="/bin/:$PATH"

WORKDIR /work

VOLUME /work

CMD []
ENTRYPOINT ["bash", "-c", "$@"]
