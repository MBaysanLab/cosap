FROM mambaorg/micromamba:1.5.10

USER root

ENV PATH=/opt/conda/bin:$PATH

RUN micromamba install -y -n base -c conda-forge \
        git \
        r-base=4.3.* \
        r-optparse \
        r-robustbase \
        r-data.table \
        r-foreach \
        r-doparallel \
        r-mclust \
        r-r.utils \
        r-rcolorbrewer \
        r-party \
        r-dbscan \
        r-umap \
        r-rcpp \
    && micromamba clean -a -y

RUN git clone --depth 1 https://github.com/imgag/ClinCNV.git /opt/ClinCNV

RUN perl -0pi -e 's/custom\.settings\$metric = "pearson"/custom.settings\$metric = "pearson"\n  custom.settings\$n_neighbors = max\(2, min\(custom.settings\$n_neighbors, ncol\(coverageForClustering\) - 1\)\)/g' \
        /opt/ClinCNV/germline/helpersGermline.R

WORKDIR /work
