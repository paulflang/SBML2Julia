Installation
============

`DisFit` depends on several python and julia packages. If you have docker installed on your machine the easiest way of installing these dependencies is to pull the latest DisFit docker image from dockerhub and build a container.::

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
