Installation
============

`SBML2JuliaMP` depends on several Python and Julia packages. If you have Docker installed on your machine, the easiest way of installing these dependencies is to pull the latest `SBML2JuliaMP docker image <https://hub.docker.com/repository/docker/paulflang/sbml2juliamp>`_ from dockerhub and build a container.:

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

`SBML2JuliaMP` uses nonlinear optimization solver `IPOPT` as a core optimization engine. The performance of optimization solver `IPOPT` critically relies on the efficiency of lienar solver that is used within the solver. In order to use efficient HSL linear solver libraries, custom installation of `IPOPT` is necessary. This option is recommended if the estimation problem faces intractability. For custom installation of `IPOPT` with HSL libraries, follow the instructions in `Ipopt.jl` repository (https://github.com/JuliaOpt/Ipopt.jl) and `IPOPT` documentation (https://coin-or.github.io/Ipopt/INSTALL.html). After custom installation of IPOPT, HSL solvers can be used, for example, by::

	problem = SBML2JuliaMP.SBML2JuliaMPProblem('tests/fixtures/G2M_copasi.xml', 'tests/fixtures/G2M_copasi.csv',optimizer_options={"linear_solver":"ma57"})
