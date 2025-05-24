FROM python:3.11-slim-bookworm
LABEL authors="austinv11"

ENV FASTP_VER="0.24.0"
ENV FASTPLONG_VER="0.2.2"

ADD https://opengene.org/fastp/fastp.${FASTP_VER} /bin/fastp
ADD https://opengene.org/fastplong/fastplong.${FASTPLONG_VER} /bin/fastplong

RUN chmod a+x /bin/fastp && chmod a+x /bin/fastplong

ENV PATH="/bin/:$PATH"

WORKDIR /work

RUN pip install multiqc==1.29

VOLUME /work

CMD []
ENTRYPOINT ["/bin/bash", "-c", "$@"]
