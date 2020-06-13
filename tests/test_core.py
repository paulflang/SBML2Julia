""" Test of DisFit.core
:Author: Paul F Lang <paul.lang@wolfson.ox.ac.uk>
:Date: 2020-04-16
:Copyright: 2020, Paul F Lang
:License: MIT
"""

import os
import pandas as pd
import pickle
import pkg_resources
import re
import shutil
import tempfile
import unittest
from DisFit import core
from numpy.testing import assert_allclose
from pandas.testing import assert_frame_equal

#Todo: write resimulations test

FIXTURES = pkg_resources.resource_filename('tests', 'fixtures')
PETAB_YAML = os.path.join(FIXTURES, 'G2M_copasi', 'G2M_copasi.yaml')
jl_file_gold = os.path.join(FIXTURES, 'jl_file_gold.jl')
with open(jl_file_gold, 'r') as f:
    JL_CODE_GOLD = f.read()
JL_CODE_GOLD = re.sub('/media/sf_DPhil_Project/Project07_Parameter Fitting/df_software/DisFit/tests/fixtures',
    FIXTURES, JL_CODE_GOLD)

class DisFitProblemTestCase(unittest.TestCase):
    
    def setUp(self):
        self.dirname = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.dirname)

    def test_constructor(self):
        problem = core.DisFitProblem(PETAB_YAML)

    def test_n_starts_setter(self):
        problem = core.DisFitProblem(PETAB_YAML)
        self.assertEqual(problem.n_starts, 1)
        with self.assertRaises(ValueError):
            problem.n_starts = 1.1
        with self.assertRaises(ValueError):
            problem.n_starts = 'a'
        self.assertEqual(problem.julia_code, JL_CODE_GOLD)
        problem.n_starts = 2
        self.assertNotEqual(problem.julia_code, JL_CODE_GOLD)

    def test_julia_code(self):
        problem = core.DisFitProblem(PETAB_YAML)
        self.assertEqual(problem.julia_code, JL_CODE_GOLD)

    def test_write_jl_file(self):
        problem = core.DisFitProblem(PETAB_YAML)
        problem.write_jl_file(path=os.path.join(self.dirname, 'jl_code.jl'))
        with open(os.path.join(self.dirname, 'jl_code.jl'), 'r') as f:
            jl_code = f.read()
        self.assertEqual(jl_code, JL_CODE_GOLD)

    def test_optimize_results_plot(self):
        # test_optimize()
        problem = core.DisFitProblem(PETAB_YAML)
        problem.optimize()

        with open(os.path.join(FIXTURES, 'results_gold.pickle'), 'rb') as f:
            results_gold = pickle.load(f)

        self.assertEqual(problem.results, results_gold) # Todo: for some reason the output is none even if it shouldn't be)
        self.assertEqual(set(results.keys()), set(['x', 'x_best', 'states', 'observables']))
        self.assertEqual(results['x'].keys(), results_gold['x'].keys())
        for i_iter in results_gold['x'].keys():
            self.assertEqual(results['x'][i_iter].keys(), results_gold['x'][i_iter].keys())
            for param in results_gold['x'][i_iter].keys():
                self.assertAlmostEqual(results['x'][i_iter][param], results_gold['x'][i_iter][param])

        self.assertEqual(results['states'].keys(), results_gold['states'].keys())
        for i_iter in results_gold['states'].keys():
            self.assertEqual(results['states'][i_iter].keys(), results_gold['states'][i_iter].keys())
            for state in results_gold['states'][i_iter].keys():
                assert_allclose(results['states'][i_iter][state], results_gold['states'][i_iter][state])

        self.assertEqual(results['observables'].keys(), results_gold['observables'].keys())
        for i_iter in results_gold['observables'].keys():
            self.assertEqual(results['observables'][i_iter].keys(), results_gold['observables'][i_iter].keys())
            for observable in results_gold['observables'][i_iter].keys():
                assert_allclose(results['observables'][i_iter][observable], results_gold['observables'][i_iter][observable])
        
        assert_frame_equal(results['x_best'], results_gold['x_best'])

        problem.write_results(path=os.path.join(self.dirname, 'results.xlsx'))
        self.assertTrue(os.path.isfile(os.path.join(self.dirname, 'results.xlsx')))

        # test_plot_results()
        problem.plot_results('wt', path=os.path.join(self.dirname, 'plot.pdf'))
        problem.plot_results('wt', path=os.path.join(self.dirname, 'plot.pdf'), observables=['obs_Ensa', 'obs_pEnsa'])
        self.assertTrue(os.path.isfile(os.path.join(self.dirname, 'plot.pdf')))


# import os
# import pandas as pd
# import pickle
# import pkg_resources
# import re
# import shutil
# import tempfile
# import unittest
# from DisFit import core
# importlib.reload(core)
# from numpy.testing import assert_allclose
# from pandas.testing import assert_frame_equal

# #Todo: write resimulations test

# FIXTURES = pkg_resources.resource_filename('tests', 'fixtures')
# PETAB_YAML = os.path.join(FIXTURES, 'G2M_copasi', 'G2M_copasi.yaml')
# jl_file_gold = os.path.join(FIXTURES, 'jl_file_gold.jl')
# with open(jl_file_gold, 'r') as f:
#     JL_CODE_GOLD = f.read()
# JL_CODE_GOLD = re.sub('/media/sf_DPhil_Project/Project07_Parameter Fitting/df_software/DisFit/tests/fixtures',
#     FIXTURES, JL_CODE_GOLD)

# problem = core.DisFitProblem(PETAB_YAML)
# problem.write_jl_file(path='jl_code_20200611.jl')
# problem.optimize()
# problem.plot_results('wt', path='plot.pdf')