FROM continuumio/miniconda3

RUN apt-get update

RUN mkdir /app
COPY . /app/.

WORKDIR /app

RUN conda env create -f environment.yml
