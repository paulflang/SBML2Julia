[![Documentation](https://readthedocs.org/projects/SBML2JuliaMP/badge/?version=latest)](https://sbml2juliamp.readthedocs.io/en/latest/)
[![Test results](https://circleci.com/gh/paulflang/SBML2JuliaMP.svg?style=shield)](https://app.circleci.com/pipelines/github/paulflang/SBML2JuliaMP)
[![License](https://img.shields.io/github/license/paulflang/SBML2JuliaMP.svg)](LICENSE)

# `SBML2JuliaMP`

`SBML2JuliaMP` is a tool to for optimizing parameters of an ordinary differential equation (ODE) model in sbml format.

## Optimization method

`SBML2JuliaMP` uses the optimization method presented in [Scalable nonlinear programming framework for parameter estimation in dynamic biological system models](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1006828). In brief, contrary to typical parameter optimization methods for ODE systems, `SBML2JuliaMP` does not rely on simulation of the ODE system. Instead `SBML2JuliaMP` uses an implicit Euler scheme to time-discretize an ODE system of n equations into m time steps. This transforms the ODE system into a system of n * (m - 1) algebraic equations with n * m variables. These n * m variables (or a subset thereof) can then be cast into an objective function. Per default, `SBML2JuliaMP` uses a least square objective. `SBML2JuliaMP` then uses interior-point optimization implemented in the Julia language to minimize the objective function constraint to the n * (m - 1) algebraic equations.

## Installation

`SBML2JuliaMP` depends on several Python and Julia packages. If you have Docker installed on your machine, the easiest way of installing these dependencies is to pull the latest [SBML2JuliaMP docker image](https://hub.docker.com/repository/docker/paulflang/sbml2juliamp) from Docker Hub and build a container.
  ```
  user@bash:/$ docker pull paulflang/sbml2juliamp:latest
  user@bash:/$ docker run -it paulflang/sbml2juliamp:latest
  ```
To install the latest `SBML2JuliaMP` version in the Docker container, run:
  ```
  user@bash:/$ git clone https://github.com/paulflang/SBML2JuliaMP.git
  user@bash:/$ python3 -m pip install -e SBML2JuliaMP
  ```
To check if the installation was succesful, run:
  ```
  user@bash:/$ SBML2JuliaMP -h
  ```

Alternatively, the `SBML2JuliaMP` dependencies can be installed as indicated in the [Dockerfile](https://github.com/paulflang/SBML2JuliaMP/blob/master/Dockerfile) in the `SBML2JuliaMP` Git repository. Once these dependencie are installed, `SBML2JuliaMP` can be installed as above:
  ```
  user@bash:/$ git clone https://github.com/paulflang/SBML2JuliaMP.git
  user@bash:/$ python3 -m pip install -e SBML2JuliaMP
  user@bash:/$ SBML2JuliaMP -h
  ```

## Interfaces

Optimization tasks can be performed from a Python API or a command line interface.

## Tutorial, and documentation
Please see the [documentation](https://sbml2juliamp.readthedocs.io/en/latest/index.html) for a description of how to use `SBML2JuliaMP`. 

## License
The package is released under the [MIT license](LICENSE).

## Development team
This package was developed by the [Paul F. Lang](https://www.linkedin.com/in/paul-lang-7b54a81a3/) at the University of Oxford, UK.


## Questions and comments
Please contact [Paul F. Lang](mailto:paul.lang@wolfson.ox.ac.uk) with any questions or comments.
