.. _python_api:

Python API
----------

The following tutorial illustrates how to use the `SBML2JuliaMP` Python API.

Importing `SBML2JuliaMP`
^^^^^^^^^^^^^^^^^^^^^^^^

Run this command to import `SBML2JuliaMP`::

    >>> import SBML2JuliaMP


Specifying an optimization problem
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

SBML2JuliaMP uses the PEtab format for specifying biological parameter estimation problems. PEtab is built around SBML and based on tab-separated values (TSV) files. Please visit the `PEtab documentation <https://petab.readthedocs.io/en/stable/documentation_data_format.html>`_ and have a look at the `PEtab examples <https://github.com/PEtab-dev/petab_test_suite/tree/master/cases>`_ for detailed instructions on how to specify an optimization problem in PEtab. If you also want to customise upper and lower boundaries for the model species, you can provide an additional species table (see `species_Vinod_FEBS2015.tsv <https://github.com/paulflang/SBML2JuliaMP/blob/main/examples/Vinod_FEBS2015/species_Vinod_FEBS2015.tsv>`_ as an example).

SBML2JuliaMP also contains the following optimization hyperparameters:

* **t_steps**: number of time-discretization steps. Default ``None``.
* **n_starts**: number of multistarts. Default ``1``.
* **infer_ic_from_sbml**: infer missing initial conditions from SBML. Default ``False``.
* **optimizer_options**: optimization solver options. Default ``{}``.
* **custom_code_dict**: dict with replaced code as keys and replacement code as values. Default ``{}``.

The problem is then specified as::

    >>> problem = SBML2JuliaMP.SBML2JuliaMPProblem(my_petab_promlem.yaml, t_steps=100, n_starts=1, infer_ic_from_sbml=False, optimizer_options={}, custom_code_dict={})

Once the problem is specified, `SBML2JuliaMP` has transformed the problem to a julia JuMP model. The code for this model can be accessed via::

    >>> code = problem.julia_code

or written to a file via::

    >>> problem.write_jl_file(path='path_to_jl_file.jl')

If you want to change the optimization problem in a way that is not yet supported by `SBML2JuliaMP`, you can manually modify the julia code and run the optimization in julia yourself. Alternatively, you can change problem.julia_code via::

    >>> problem.insert_custom_code({'<replaced lines>': '<replacement lines>'})

Running the optimization
^^^^^^^^^^^^^^^^^^^^^^^^

The optimization can be run with::

    >>> problem.optimize()

Please note that this may take a while.

Accessing the results
^^^^^^^^^^^^^^^^^^^^^

The results can be accessed via::

    >>> results = problem.results

and written to an excel file via::

    >>> problem.write_results(path='path_to_results.xlsx')

Time courses for the optimal solution of condition ``cond`` and corresponding experimental datapoints can be plotted by::

    >>> problem.plot_results(cond, path='path_to_plot.pdf', observables=[], size=(6, 5))

where the optional ``observable`` argument accepts a list of observables that shall be plotted. The optional ``size`` argument specifies the size of the figure.

If you want to update the nominal parameter values in your PEtab problem parameter table with the fitted values, run::

    >>> problem.write_optimized_Parameter_table()

This will create a `post_fit_parameters.tsv` file in your PEtab problem directory. This can be useful to perform sensitivity analysis in other PEtab compatible optimization toolboxes such as `PyPesto <https://pypesto.readthedocs.io/en/latest/>`_.
