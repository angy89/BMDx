docker_run: docker_build
	docker run --rm -it -v$(shell pwd):/home/docker/shiny_app -p8080:8080 bmdx

docker_build: Dockerfile 
	docker build -t bmdx .

run_shiny: install_dependencies
	Rscript shiny_start.r 

install_dependencies: 
	Rscript install_dependencies.R
