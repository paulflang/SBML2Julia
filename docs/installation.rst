Installation
============

`DisFit` depends on several Python and Julia packages. If you have Docker installed on your machine, the easiest way of installing these dependencies is to pull the latest `DisFit docker image <https://hub.docker.com/repository/docker/paulflang/disfit>`_ from dockerhub and build a container.:

	docker pull paulflang/disfit:latest
	docker run -it paulflang/disfit:latest

To install the latest `DisFit` version in the docker container, run::

	git clone https://github.com/paulflang/DisFit.git
	python3 -m pip install -e DisFit

To check if the installation was succesful, run::

	DisFit -h


Alternatively, the `DisFit` dependencies can be installed as indicated in the `Dockerfile <https://github.com/paulflang/DisFit/blob/master/Dockerfile>`_ in the `DisFit` Git repository. Once these dependencie are installed, `DisFit` can be installed as above::

	git clone https://github.com/paulflang/DisFit.git
	python3 -m pip install -e DisFit
	DisFit -h