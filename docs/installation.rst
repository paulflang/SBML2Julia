Installation
============

`SBML2Julia` depends on several Python and Julia packages. If you have Docker installed on your machine, the easiest way of installing these dependencies is to pull the latest `SBML2Julia docker image <https://hub.docker.com/repository/docker/paulflang/sbml2julia>`_ from dockerhub and build a container.::

	docker pull paulflang/sbml2julia:latest
	docker run -it paulflang/sbml2julia:latest

To install the latest `SBML2Julia` version in the docker container, run::

	git clone https://github.com/paulflang/sbml2julia.git
	python3 -m pip install -e sbml2julia

To check if the installation was succesful, run::

	sbml2julia -h


Alternatively, the `SBML2Julia` dependencies can be installed on Ubuntu machines as indicated in the `Dockerfile <https://github.com/paulflang/sbml2julia/blob/master/Dockerfile>`_ in the `SBML2Julia` GitHub repository. Once these dependencie are installed, `SBML2Julia` can be installed as above::

	git clone https://github.com/paulflang/sbml2julia.git
	python3 -m pip install -e sbml2julia
	sbml2julia -h

Optional installation of further linear solvers
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To solve more involved problems, we recommend the linear solvers provided in HSL, for which academics can `request a free license <http://www.hsl.rl.ac.uk/>`_.