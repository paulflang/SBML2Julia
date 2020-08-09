.. _python_api:

Python API
----------

The following tutorial illustrates how to use the `DisFit` Python API.

Importing `DisFit`
^^^^^^^^^^^^^^^^^^^

Run this command to import `DisFit`::

    >>> import DisFit


Specifying an optimization problem
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

DisFit uses the PEtab format for specifying biological parameter estimation problems. PEtab is built around SBML and based on tab-separated values (TSV) files. Please visit the PEtab docs and have a look at the PEtab examples for detailed instructions on how to specify an optimization problem in PEtab.

DisFit also contains the following optimization hyperparameters:

* **t_ratio**: ratio between experimental observation intervals and simulation time-discretization intervals. Default ``2``.
* **n_starts**: number of multistarts. Default ``1``.
* **infer_ic_from_sbml**: if missing initial conditions shall be infered from SBML model. Default ``False``.

The problem is then specified as::

    >>> problem = DisFit.DisFitProblem(petab_promlem.yaml, t_ratio=2, n_starts=1, infer_ic_from_sbml=False)

Once the problem is specified, `DisFit` has transformed the problem to a julia JuMP model. The code for this model can be accessed via::

    >>> code = problem.julia_code

or written to a file via::

    >>> problem.write_jl_file(path='path_to_jl_file.jl')

If you want to change the optimization problem in a way that is not yet supported by `DisFit`, you can manually modify the julia code and run the optimization in julia yourself.

Running the optimization
^^^^^^^^^^^^^^^^^^^^^^^^

The optimization can be run with::

    >>> problem.optimize()

Please note that this may take a while.

Accessing the results
^^^^^^^^^^^^^^^^^^^^^

The results can be accessed via::

    >>> results = problem.results

written to an excel file via::

    >>> problem.wirte_results(path='path_to_results.xlsx')

Time courses for the optimal solution of condition `cond` and corresponding experimental datapoints can be plotted by::

    >>> problem.plot_results(cond, path='path_to_plot.pdf', observables=[], size=(6, 5))

where the optional ``observable`` argument accepts a list of observables that shall be plotted. The optional ``size`` argument specifies the size of the figure.