FROM condaforge/mambaforge

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
    git g++ cmake autoconf libtool liblzma-dev zlib1g-dev libbz2-dev libcurl3-dev libssl-dev \
    ca-certificates cpp make libltdl-dev wget unzip \
    && apt-get clean && rm -rf /var/lib/apt/lists/

WORKDIR /tmp

# Install MuSE
RUN wget https://github.com/wwylab/MuSE/archive/refs/tags/v2.0.1.tar.gz && \
    tar -vxf v2.0.1.tar.gz && \
    cd MuSE-2.0.1 && \
    ./install_muse.sh && \
    cp MuSE /usr/local/bin/

ENV PATH=$PATH:/usr/local/bin/MuSE

# Install VEP
COPY --from=ensemblorg/ensembl-vep /opt/vep/ /opt/vep/
COPY --from=ensemblorg/ensembl-vep /etc/perl/ /etc/perl/
COPY --from=ensemblorg/ensembl-vep /usr/local/bin/ /usr/local/bin/

ENV OPT /opt/vep
ENV OPT_SRC $OPT/src
ENV PATH $OPT_SRC/ensembl-vep:$OPT_SRC/var_c_code:$PATH
RUN echo >> $OPT/.profile && \
    echo PATH=$PATH:\$PATH >> $OPT/.profile && \
    echo export PATH >> $OPT/.profile


# Install COSAP
COPY requirements.txt /tmp/requirements.txt
RUN --mount=type=cache,target=/opt/conda/pkgs mamba install -c conda-forge -c bioconda --yes --name base --file /tmp/requirements.txt
RUN mkdir /app
COPY . /app/.

WORKDIR /app
RUN mamba run --no-capture-output -n base pip install celery redis
RUN mamba run --no-capture-output -n base pip install -e .

ENV COSAP /app
ENV COSAP_LIBRARY_PATH /cosap_data

RUN mkdir /workdir
WORKDIR /workdir
