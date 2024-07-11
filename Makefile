PWD = $(shell pwd)

docker_build:
	DOCKER_BUILDKIT=1 docker build --platform linux/amd64 -t itubioinformatics/cosap:$(version) -t itubioinformatics/cosap:latest . 

install:
	# Check if miniforge3 is installed and install if not
	@if [ ! -d "$(HOME)/miniforge3" ]; then \
		wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh; \
		bash Miniforge3-Linux-x86_64.sh -b -p $(HOME)/miniforge3; \
		rm Miniforge3-Linux-x86_64.sh; \
	fi;

	mamba create --name cosap -c conda-forge -c bioconda --yes --file requirements.txt

	# Set environment variables in conda environment
	conda env config vars set COSAP=$(PWD) -n cosap

develop:
	# Check if miniforge3 is installed and install if not
	@if [ ! -d "$(HOME)/miniforge3" ]; then \
		wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh; \
		bash Miniforge3-Linux-x86_64.sh -b -p $(HOME)/miniforge3; \
		rm Miniforge3-Linux-x86_64.sh; \
	fi;
	mamba create --name cosap -c conda-forge -c bioconda --yes --file requirements.txt
	mamba run --no-capture-output -n cosap pip install -e . celery redis docker

download:
	@if [ -n $(COSAP_LIBRARY_PATH) ] && [ "$(COSAP_LIBRARY_PATH)" != "" ] ; \
		then echo Reading COSAP_LIBRARY_PATH: $$COSAP_LIBRARY_PATH; fi
	@COSAP_DATA=$${COSAP_LIBRARY_PATH:-$(HOME)/cosap_data/}; \
	\
	echo The files will be downloaded to \` $$COSAP_DATA \`.; \
	echo The directory will be created if it does not exist.; \
	read -p "Continue? [Y/n]: " line; \
	if [ -z $$line ] || [ $$line = "y" ] || [ $$line = "Y" ]; \
	then \
		mkdir -p $$HOME/cosap_data; \
		wget -nc -i ./required_files_test.txt -P $$COSAP_DATA; \
	else \
		echo No files downloaded.; \
	fi
