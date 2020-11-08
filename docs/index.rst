.. image:: https://circleci.com/gh/paulflang/SBML2Julia.svg?style=shield
   :target: https://app.circleci.com/pipelines/github/paulflang/SBML2Julia
   :alt: Build status
.. image:: https://readthedocs.org/projects/sbml2julia/badge/?version=latest
   :target: https://sbml2julia.readthedocs.io/en/latest/?badge=latest
   :alt: Documentation Status
.. image:: https://img.shields.io/github/license/paulflang/sbml2julia.svg
   :target: LICENSE
   :alt: License

Welcome to `SBML2Julia`'s documentation!
==========================================

`SBML2Julia` is a tool to for optimizing parameters of ordinary differential equation (ODE) models. `SBML2Julia` translates a model from SBML/`PEtab <https://petab.readthedocs.io/en/stable/>`_ format into Julia for Mathematical Programming (`JuMP <https://jump.dev/JuMP.jl/stable/>`_), performs the optimization task and returns the results.

| Source code: https://github.com/paulflang/SBML2Julia

Optimization method
-------------------

`SBML2Julia` uses the optimization method presented in `Scalable nonlinear programming framework for parameter estimation in dynamic biological system models <https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1006828>`_. In brief, contrary to typical parameter optimization methods for ODE systems, `SBML2Julia` does not rely on simulation of the ODE system. Instead `SBML2Julia` uses an implicit Euler scheme to time-discretize an ODE system of n equations into m time steps. This transforms the ODE system into a system of n * (m - 1) algebraic equations with n * m variables. These n * m variables (or a subset thereof) can then be cast into an objective function. Per default, `SBML2Julia` uses a least square objective. `SBML2Julia` then uses interior-point optimization implemented in the Julia language to minimize the objective function constraint to the n * (m - 1) algebraic equations.

Interfaces
----------

Optimization tasks can be performed from a Python API or a command line interface.

Contents
--------

.. toctree::
   :maxdepth: 2
   :numbered:

   installation.rst
   tutorial.rst
   examples.rst
   known_limitations.rst
   API_reference.rst
   contributing.rst
   about.rst
