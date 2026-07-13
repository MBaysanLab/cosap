FROM condaforge/mambaforge

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
    git g++ cmake autoconf libtool liblzma-dev zlib1g-dev libbz2-dev libcurl3-dev libssl-dev \
    ca-certificates cpp make libltdl-dev wget unzip graphviz software-properties-common \
    curl cpio rpm2cpio \
    && apt-get clean && rm -rf /var/lib/apt/lists/

RUN curl -fsSL https://raw.githubusercontent.com/apptainer/apptainer/main/tools/install-unprivileged.sh \
    | bash -s - /usr/local/apptainer \
    && ln -s /usr/local/apptainer/bin/apptainer /usr/local/bin/apptainer \
    && test -x /usr/local/apptainer/bin/apptainer \
    && /usr/local/apptainer/bin/apptainer --version

WORKDIR /tmp

# Install COSAP
COPY requirements.txt /tmp/requirements.txt
RUN --mount=type=cache,target=/opt/conda/pkgs mamba create --name cosap -c conda-forge -c bioconda --file /tmp/requirements.txt
RUN mkdir /app
COPY . /app/.

WORKDIR /app
RUN mamba run --no-capture-output -n cosap uv pip install --no-cache .

# Set environment variables
ENV COSAP=/app
ENV COSAP_LIBRARY_PATH=/cosap_data

# Activate conda environment
RUN conda init
RUN echo "conda activate cosap" >> ~/.bashrc
ENV PATH=/opt/conda/envs/cosap/bin:$PATH

# Create non-root user and working directory
RUN groupadd -r cosap && useradd -r -g cosap -m -s /bin/bash cosap
RUN mkdir /workdir && chown cosap:cosap /workdir

USER cosap
WORKDIR /workdir
