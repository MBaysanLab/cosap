build:
	DOCKER_BUILDKIT=1 docker build -t itubioinformatics/cosap .

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

cosap_data_path = $(HOME)/cosap_data
download_files:
	@echo The files will be downloaded to $(cosap_data_path)\; continue? [Y/n]
	@read line; if [ $$line = "n" ]; then echo aborting; exit 1 ; fi
	wget -i ./required_files.txt -P $(cosap_data_path)