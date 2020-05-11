Example
=======

The `DisFit GitHub repository <https://github.com/paulflang/DisFit/tree/master/tests/fixtures>`_ contains a version of the `Vinod et Novak model <https://www.sciencedirect.com/science/article/pii/S0014579315000873>`_ of the G2/M cell cycle transition, along with simulated experimental data.

Using the Python API
--------------------

The `DisFit` problem can be created using the Python API (and assuming that the current working directory is the `DisFit` root directory) via::

	>>> import DisFit

	>>> problem = DisFit.DisFitProblem('tests/fixtures/G2M_copasi.xml', 'tests/fixtures/G2M_copasi.csv')

and solved by::

	>>> problem.optimize()

The results are then available under ``problem.results``, which returns a dictionary containing the ``'states'``, all found parameter sets ``'x'`` and the best parameter set  ``'x_best'``. For example, to access the best parameter set, type::

	>>> problem.results['x_best']                                                                        
	        Name      x_0     x_best  x_best_to_x_0
	0   kAspEB55  57.0000  47.529257       0.833847
	1   kCdc25_1   0.1000   0.121924       1.219243
	2   kCdc25_2   0.9000   1.024421       1.138246
	3   kDipEB55   0.0068   0.003400       0.500004
	4   kDpCdc25  10.0000  10.282131       1.028213
	5    kDpEnsa   0.0500   0.049753       0.995059
	6     kDpGw1   0.2500   0.260846       1.043383
	7     kDpGw2  10.0000  10.450302       1.045030
	8     kDpWee  10.0000  10.422349       1.042235
	9   kPhCdc25   1.0000   1.051618       1.051618
	10   kPhEnsa   0.1000   0.099802       0.998017
	11     kPhGw   1.0000   1.047174       1.047174
	12    kPhWee   1.0000   1.061641       1.061641
	13     kWee1   0.0100   0.010209       1.020938
	14     kWee2   0.9900   1.188540       1.200546

Using the command line interface
--------------------------------

Similarly, the same example problem can be solved from the comand line interface::

	user@bash:/DisFit$ DisFit optimize 'tests/fixtures/G2M_copasi.xml' 'tests/fixtures/G2M_copasi.csv' -o './DisFit_results'

The results can be found in the output directory given to the ``-o`` argument.