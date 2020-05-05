Welcome to `DisFit`'s documentation!
====================================

`DisFit` is a tool to for optimizing parameters of an ordinary differential equation (ODE) model in sbml format.

Optimization method
-------------------

`DisFit` uses the optimization method presented in `Scalable nonlinear programming framework for parameter estimation in dynamic biological system models <https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1006828>`_. In brief, contrary to typical parameter optimization methods for ODE systems, `DisFit` does not rely on simulation of the ODE system. Instead `DisFit` uses an implicit Euler scheme to time-discretize an ODE system of n equations into m time steps. This transforms the ODE system into a system of n * (m - 1) algebraic equations with n * m variables. These n * m variables (or a subset thereof) can then be cast into an objective function. Per default, `DisFit` uses a least square objective. `DisFit` then uses interior-point optimization implemented in the Julia language to minimize the objective function constraint to the n * (m - 1) algebraic equations.

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
   example.rst
   known_limitations.rst
   contributing.rst
   about.rst