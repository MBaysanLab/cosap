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

# Install COSAP
COPY requirements.txt /tmp/requirements.txt
RUN --mount=type=cache,target=/opt/conda/pkgs mamba create --name cosap -c conda-forge -c bioconda --file /tmp/requirements.txt
RUN mkdir /app
COPY . /app/.

WORKDIR /app
RUN mamba run --no-capture-output -n cosap pip install celery redis docker
RUN mamba run --no-capture-output -n cosap pip install -e .

# Set environment variables
ENV COSAP /app
ENV COSAP_LIBRARY_PATH /cosap_data

# Activate conda environment
RUN conda init
RUN echo "conda activate cosap" >> ~/.bashrc
ENV PATH /opt/conda/envs/cosap/bin:$PATH

# Create working directory
RUN mkdir /workdir
WORKDIR /workdir