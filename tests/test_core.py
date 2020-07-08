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

    def test_t_ratio_setter(self):
        problem = core.DisFitProblem(PETAB_YAML)
        self.assertEqual(problem.t_ratio, 2)
        with self.assertRaises(ValueError):
            problem.t_ratio = -1
        with self.assertRaises(ValueError):
            problem.t_ratio = 'a'
        self.assertEqual(problem.julia_code, JL_CODE_GOLD)
        problem.t_ratio = 1
        self.assertNotEqual(problem.julia_code, JL_CODE_GOLD)

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
        problem = core.DisFitProblem(PETAB_YAML, t_ratio=1.99999)
        problem.optimize()
        results = problem.results

        with open(os.path.join(FIXTURES, 'results_gold.pickle'), 'rb') as f:
            results_gold = pickle.load(f)

        # self.assertEqual(problem.results, results_gold) # Todo: for some reason the output is none even if it shouldn't be)
        self.assertEqual(set(results.keys()), set(['parameters', 'par_best', 'species', 'observables']))
        self.assertEqual(results['parameters'].keys(), results_gold['parameters'].keys())
        for i_iter in results_gold['parameters'].keys():
            self.assertEqual(results['parameters'][i_iter].keys(), results_gold['parameters'][i_iter].keys())
            for param in results_gold['parameters'][i_iter].keys():
                assert_allclose(results['parameters'][i_iter][param], results_gold['parameters'][i_iter][param])

        self.assertEqual(results['species'].keys(), results_gold['species'].keys())
        for i_iter in results_gold['species'].keys():
            self.assertEqual(results['species'][i_iter].keys(), results_gold['species'][i_iter].keys())
            for specie in results_gold['species'][i_iter].keys():
                assert_allclose(results['species'][i_iter][specie],
                    results_gold['species'][i_iter][specie], rtol=1e-05, atol=1e-08)

        self.assertEqual(results['observables'].keys(), results_gold['observables'].keys())
        for i_iter in results_gold['observables'].keys():
            self.assertEqual(results['observables'][i_iter].keys(), results_gold['observables'][i_iter].keys())
            for observable in results_gold['observables'][i_iter].keys():
                assert_allclose(results['observables'][i_iter][observable],
                    results_gold['observables'][i_iter][observable], rtol=1e-05, atol=1e-08)
        
        assert_frame_equal(results['par_best'], results_gold['par_best'])

        problem.write_results(path=os.path.join(self.dirname, 'results.xlsx'))
        self.assertTrue(os.path.isfile(os.path.join(self.dirname, 'results.xlsx')))

        # test_plot_results()
        problem.plot_results('wt', path=os.path.join(self.dirname, 'plot.pdf'))
        problem.plot_results('wt', path=os.path.join(self.dirname, 'plot.pdf'), observables=['obs_Ensa', 'obs_pEnsa'])
        self.assertTrue(os.path.isfile(os.path.join(self.dirname, 'plot.pdf')))


import os
import pandas as pd
import pickle
import pkg_resources
import re
import shutil
import tempfile
import unittest
from DisFit import core
importlib.reload(core)
from numpy.testing import assert_allclose
from pandas.testing import assert_frame_equal

#Todo: write resimulations test

FIXTURES = pkg_resources.resource_filename('tests', 'fixtures')
PETAB_YAML = os.path.join(FIXTURES, 'petab_suite_0001', '_0001.yaml')
jl_file_gold = os.path.join(FIXTURES, 'jl_file_gold.jl')
with open(jl_file_gold, 'r') as f:
    JL_CODE_GOLD = f.read()
JL_CODE_GOLD = re.sub('/media/sf_DPhil_Project/Project07_Parameter Fitting/df_software/DisFit/tests/fixtures',
    FIXTURES, JL_CODE_GOLD)

problem = core.DisFitProblem(PETAB_YAML)
problem.write_jl_file(path='jl_code_20200708.jl')
problem.optimize()
problem.plot_results('wt', path='plot.pdf')
problem.write_results()
# problem.results['par_best']