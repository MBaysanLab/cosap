
up: build start

build:
	docker-compose build 

start:
	docker-compose start -d --name=cosap
	