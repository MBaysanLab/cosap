
up: build start

build:
	docker build . 

start:
	docker run -it cosap /bin/bash

