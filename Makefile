
up: build start

build:
	docker build . 

start:
	docker run -it cosap /bin/bash

install:
	conda install -c conda-forge mamba
	mamba env create -f environments/default_environment.yml
	conda activate cosap
	pip install .
	conda deactivate

develop:
	conda install -c conda-forge mamba
	mamba env create -f environments/default_environment.yml
	conda activate cosap
	pip install -e .
	conda deactivate