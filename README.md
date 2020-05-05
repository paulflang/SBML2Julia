[![Documentation](https://readthedocs.org/projects/disfit/badge/?version=latest)](https://disfit.readthedocs.io/en/documentation/)
[![Test results](https://circleci.com/gh/paulflang/disfit.svg?style=shield)](https://app.circleci.com/pipelines/github/paulflang/DisFit)
[![License](https://img.shields.io/github/license/paulflang/disfit.svg)](LICENSE)

# `DisFit`

`DisFit` is a tool to for optimizing parameters of an ordinary differential equation (ODE) model in sbml format.

## Optimization method

`DisFit` uses the optimization method presented in [Scalable nonlinear programming framework for parameter estimation in dynamic biological system models](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1006828). In brief, contrary to typical parameter optimization methods for ODE systems, `DisFit` does not rely on simulation of the ODE system. Instead `DisFit` uses an implicit Euler scheme to time-discretize an ODE system of n equations into m time steps. This transforms the ODE system into a system of n * (m - 1) algebraic equations with n * m variables. These n * m variables (or a subset thereof) can then be cast into an objective function. Per default, `DisFit` uses a least square objective. `DisFit` then uses interior-point optimization implemented in the Julia language to minimize the objective function constraint to the n * (m - 1) algebraic equations.

## Installation

`DisFit` depends on several Python and Julia packages. If you have Docker installed on your machine, the easiest way of installing these dependencies is to pull the latest [DisFit docker image](https://hub.docker.com/repository/docker/paulflang/disfit) from Docker Hub and build a container.
  ```
  user@bash:/$ docker pull paulflang/disfit:latest
  user@bash:/$ docker run -it paulflang/disfit:latest
  ```
To install the latest `DisFit` version in the Docker container, run:
  ```
  user@bash:/$ git clone https://github.com/paulflang/DisFit.git
  user@bash:/$ python3 -m pip install -e DisFit
  ```
To check if the installation was succesful, run:
  ```
  user@bash:/$ DisFit -h
  ```

Alternatively, the `DisFit` dependencies can be installed as indicated in the [Dockerfile](https://github.com/paulflang/DisFit/blob/master/Dockerfile) in the `DisFit` Git repository. Once these dependencie are installed, `DisFit` can be installed as above:
  ```
  user@bash:/$ git clone https://github.com/paulflang/DisFit.git
  user@bash:/$ python3 -m pip install -e DisFit
  user@bash:/$ DisFit -h
  ```

## Interfaces

Optimization tasks can be performed from a Python API or a command line interface.

## Tutorial, and documentation
Please see the [documentation](https://disfit.readthedocs.io/en/documentation/index.html) for a description of how to use `DisFit`. 

## License
The package is released under the [MIT license](LICENSE).

## Development team
This package was developed by the [Paul F. Lang](https://www.linkedin.com/in/paul-lang-7b54a81a3/) at the University of Oxford, UK.


## Questions and comments
Please contact the [Paul F. Lang](mailto:paul.lang@wolfson.ox.ac.uk) with any questions or comments.
