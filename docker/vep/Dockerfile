#FROM alpine:latest
FROM perl:5-slim
LABEL authors="austinv11"

# Install Aria2, gcc, g++, make
RUN apt-get update && \
    apt-get install -y aria2 gcc g++ make build-essential libdbd-mysql-perl libmariadb-dev default-libmysqlclient-dev libexpat1-dev libdb-dev zlib1g-dev unzip gzip tar && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Install Perl libraries: Archive::Zip and DBI
RUN cpanm Archive::Zip DBI Archive::Extract

# Install extra dependencies
RUN cpanm --notest Set::IntervalTree JSON PerlIO::gzip

ENV VEP_VERSION="114.1"

ENV VEP_URL="https://github.com/Ensembl/ensembl-vep/archive/refs/tags/release/${VEP_VERSION}.tar.gz"

WORKDIR /build
ADD $VEP_URL /build/

RUN tar -xzf ${VEP_VERSION}.tar.gz && \
    rm ${VEP_VERSION}.tar.gz && \
    mv ensembl-vep-release-${VEP_VERSION} ensembl-vep && \
    cd ensembl-vep && \
    perl INSTALL.pl --AUTO acfp --ASSEMBLY GRCh38 --PLUGINS all --SPECIES homo_sapiens --NO_UPDATE

ENV PATH="/build/ensembl-vep:${PATH}"

WORKDIR /work

VOLUME /work

CMD []
ENTRYPOINT ["bash", "-c", "$@"]
