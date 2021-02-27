FROM continuumio/miniconda3

RUN apt-get update

MKDIR /app
COPY . /app/.

WORKDIR /app

RUN conda install -c bioconda --file requirements.txt

