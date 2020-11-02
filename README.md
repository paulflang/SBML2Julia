[![Documentation Status](https://readthedocs.org/projects/sbml2julia/badge/?version=latest)](https://sbml2julia.readthedocs.io/en/latest/?badge=latest)
[![Test results](https://circleci.com/gh/paulflang/sbml2julia.svg?style=shield)](https://app.circleci.com/pipelines/github/paulflang/sbml2julia)
[![License](https://img.shields.io/github/license/paulflang/sbml2julia.svg)](LICENSE)

# `SBML2Julia`

`SBML2Julia` is a tool to for optimizing parameters of ordinary differential equation (ODE) models. `SBML2Julia` translates a model from SBML/[PEtab](https://petab.readthedocs.io/en/stable/) format into Julia for Mathematical Programming ([JuMP](https://jump.dev/JuMP.jl/stable/)), performes the optimization task and returns the results.

## Optimization method

`SBML2Julia` uses the optimization method presented in [Scalable nonlinear programming framework for parameter estimation in dynamic biological system models](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1006828). In brief, contrary to typical parameter optimization methods for ODE systems, `SBML2Julia` does not rely on simulation of the ODE system. Instead `SBML2Julia` uses an implicit Euler scheme to time-discretize an ODE system of n equations into m time steps. This transforms the ODE system into a system of n * (m - 1) algebraic equations with n * m variables. These n * m variables (or a subset thereof) can then be cast into an objective function. `SBML2Julia` then uses interior-point optimization implemented in the Julia language to minimize the objective function constraint to the n * (m - 1) algebraic equations.

## Installation

`SBML2Julia` depends on several Python and Julia packages. If you have Docker installed on your machine, the easiest way of installing these dependencies is to pull the latest [SBML2Julia docker image](https://hub.docker.com/repository/docker/paulflang/sbml2julia) from Docker Hub and build a container.
  ```
  user@bash:/$ docker pull paulflang/sbml2julia:latest
  user@bash:/$ docker run -it paulflang/sbml2julia:latest
  ```
To install the latest `SBML2Julia` version in the Docker container, run:
  ```
  user@bash:/$ git clone https://github.com/paulflang/sbml2julia.git
  user@bash:/$ python3 -m pip install -e sbml2julia
  ```
To check if the installation was succesful, run:
  ```
  user@bash:/$ sbml2julia -h
  ```

Alternatively, the `SBML2Julia` dependencies can be installed as indicated in the [Dockerfile](https://github.com/paulflang/sbml2julia/blob/master/Dockerfile) in the `SBML2Julia` GitHub repository. Once these dependencie are installed, `SBML2Julia` can be installed as above:
  ```
  user@bash:/$ git clone https://github.com/paulflang/sbml2julia.git
  user@bash:/$ python3 -m pip install -e sbml2julia
  user@bash:/$ sbml2julia -h
  ```

## Interfaces

Optimization tasks can be performed from a Python API or a command line interface.

## Tutorial, and documentation
Please see the [documentation](https://sbml2julia.readthedocs.io/en/latest/index.html) for a description of how to use `SBML2Julia`. 

## License
The package is released under the [MIT license](LICENSE).

## Development team
This package was developed by [Paul F. Lang](https://www.linkedin.com/in/paul-lang-7b54a81a3/) at the University of Oxford, UK and [Sungho Shin](https://www.sunghoshin.com/) at the University of Wisconsin-Madison, USA..


## Questions and comments
Please contact [Paul F. Lang](mailto:paul.lang@wolfson.ox.ac.uk) with any questions or comments.
