FROM condaforge/mambaforge

# Set ensembl-vep
COPY --from=ensemblorg/ensembl-vep /opt/vep/ /opt/vep/
COPY --from=ensemblorg/ensembl-vep /etc/perl/ /etc/perl/
COPY --from=ensemblorg/ensembl-vep /usr/local/bin/ /usr/local/bin/

ENV OPT /opt/vep
ENV OPT_SRC $OPT/src
ENV PATH $OPT_SRC/ensembl-vep:$OPT_SRC/var_c_code:$PATH
RUN echo >> $OPT/.profile && \
    echo PATH=$PATH:\$PATH >> $OPT/.profile && \
    echo export PATH >> $OPT/.profile

RUN mkdir /app
COPY . /app/.

WORKDIR /app
RUN --mount=type=cache,target=/opt/conda/pkgs mamba install -c bioconda -c conda-forge --yes --name base --file requirements.txt
RUN mamba run --no-capture-output -n base pip install celery redis
RUN mamba run --no-capture-output -n base pip install -e .

ENV COSAP /app
ENV COSAP_LIBRARY_PATH /cosap_data

RUN mkdir /workdir
WORKDIR /workdir
