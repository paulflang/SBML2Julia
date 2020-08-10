Example
=======

The `DisFit GitHub repository <https://github.com/paulflang/DisFit/tree/master/tests/fixtures>`_ contains a version of the `Vinod et Novak model <https://www.sciencedirect.com/science/article/pii/S0014579315000873>`_ of the G2/M cell cycle transition, along with simulated experimental data.

Using the Python API
--------------------

The `DisFit` problem can be created using the Python API (and assuming that the current working directory is the `DisFit` root directory) via::

	>>> import DisFit

	>>> problem = DisFit.DisFitProblem('tests/fixtures/G2M_copasi/G2M_copasi.yaml')

and solved by::

	>>> problem.optimize()

The results are then available under ``problem.results``, which returns a dictionary containing ``'par_best'`` (the best found parameter set), ``'species'``, ``'observables'``, ``'fval'`` (the negative log-likelihood) and ``'chi2'`` (chi2 values of residuals). For example, to access the best parameter set, type::

>>> problem.results['x_best']                                                                        
                Name    par_0   par_best  par_best_to_par_0
0        fB55_Cb_low   1.1000   1.375000           1.250000
1          fB55_iWee   0.9000   1.125000           1.250000
2      fB55_pGw_weak   1.0000   1.250000           1.250000
3            fB55_wt   1.0000   1.250000           1.250000
4                fCb   2.0000   2.000040           1.000020
5              jiWee   0.1000   0.114439           1.144389
6           kAspEB55  57.0000  52.057389           0.913288
7           kCdc25_1   0.1000   0.103031           1.030312
8           kCdc25_2   0.9000   0.913392           1.014880
9           kDipEB55   0.0068   0.003419           0.502858
10          kDpCdc25  10.0000   9.841409           0.984141
11           kDpEnsa   0.0500   0.049844           0.996879
12            kDpGw1   0.2500   0.248459           0.993836
13            kDpGw2  10.0000   9.993553           0.999355
14            kDpWee  10.0000  10.052517           1.005252
15          kPhCdc25   1.0000   0.995569           0.995569
16    kPhEnsa_Cb_low   0.1000   0.125000           1.250000
17      kPhEnsa_iWee   0.1000   0.125000           1.250000
18  kPhEnsa_pGw_weak   0.0900   0.112500           1.250000
19        kPhEnsa_wt   0.1000   0.125000           1.250000
20             kPhGw   1.0000   0.995924           0.995924
21            kPhWee   1.0000   1.013151           1.013151
22             kWee1   0.0100   0.009879           0.987921
23             kWee2   0.9900   1.016440           1.026707


Using the command line interface
--------------------------------

Similarly, the same example problem can be solved from the comand line interface::

	user@bash:/DisFit$ DisFit optimize 'tests/fixtures/G2M_copasi/G2M_copasi.yaml' -o './DisFit_results'

The results can be found in the output directory given to the ``-o`` argument.