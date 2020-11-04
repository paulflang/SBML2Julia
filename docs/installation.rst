Installation
============

`SBML2Julia` depends on several Python and Julia packages. If you have Docker installed on your machine, the easiest way of installing these dependencies is to pull the latest `SBML2Julia docker image <https://hub.docker.com/repository/docker/paulflang/sbml2julia>`_ from dockerhub and build a container.::

    user@bash:/$ docker pull paulflang/sbml2julia:latest
    user@bash:/$ docker run -it --mount type=bind,source=<my_host_dir>,target=/media paulflang/sbml2julia:latest

To install the latest `SBML2Julia` release in the docker container, run::

    user@bash:/$ python3 -m pip install sbml2julia

Alternatively, to install the latest `SBML2Julia` version from GitHub, run::

    user@bash:/$ git clone https://github.com/paulflang/sbml2julia.git
    user@bash:/$ python3 -m pip install sbml2julia

To check if the installation was succesful, run::

    user@bash:/$ sbml2julia -h

If you do not want to use Docker, the `SBML2Julia` dependencies can be installed on Ubuntu machines as indicated in the `Dockerfile <https://github.com/paulflang/sbml2julia/blob/master/Dockerfile>`_. Once these dependencie are installed, `SBML2Julia` can be installed as above.

Optional installation of efficient HSL linear solvers
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

`SBML2Julia` uses the nonlinear optimization solver `Ipopt` as core optimization engine. Its performance critically relies on the efficiency of the linear solver used within `Ipopt`. If the estimation problem faces intractability, we recommend custom installation of efficient HSL linear solvers. Since HSL linear solvers run under a different license that SBML2Julia, we cannot distribute them with `SBML2Julia`. However, academics can request a `free license for HSL linear solvers <http://www.hsl.rl.ac.uk/ipopt/>`_. Using these HSL linear solvers within `SBML2Julia` requires custom `installation of Ipopt <https://coin-or.github.io/Ipopt/INSTALL.html>`_ and its `Julia interface <https://github.com/JuliaOpt/Ipopt.jl>`_.
