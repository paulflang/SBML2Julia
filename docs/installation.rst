Installation
============

`SBML2JuliaMP` depends on several Python and Julia packages. If you have Docker installed on your machine, the easiest way of installing these dependencies is to pull the latest `SBML2JuliaMP docker image <https://hub.docker.com/repository/docker/paulflang/sbml2juliamp>`_ from dockerhub and build a container.::

	docker pull paulflang/sbml2juliamp:latest
	docker run -it paulflang/sbml2juliamp:latest

To install the latest `SBML2JuliaMP` version in the docker container, run::

	git clone https://github.com/paulflang/SBML2JuliaMP.git
	python3 -m pip install -e SBML2JuliaMP

To check if the installation was succesful, run::

	SBML2JuliaMP -h


Alternatively, the `SBML2JuliaMP` dependencies can be installed on Ubuntu machines as indicated in the `Dockerfile <https://github.com/paulflang/SBML2JuliaMP/blob/master/Dockerfile>`_ in the `SBML2JuliaMP` GitHub repository. Once these dependencie are installed, `SBML2JuliaMP` can be installed as above::

	git clone https://github.com/paulflang/SBML2JuliaMP.git
	python3 -m pip install -e SBML2JuliaMP
	SBML2JuliaMP -h

Optional installation of further linear solvers
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To solve more involved problems, we recommend the linear solvers provided in HSL, for which academics can `request a free license <http://www.hsl.rl.ac.uk/>`_.