""" Test of DisFit.core
:Author: Paul F Lang <paul.lang@wolfson.ox.ac.uk>
:Date: 2020-04-16
:Copyright: 2020, Paul F Lang
:License: MIT
"""

import filecmp
import os
import pandas as pd
import pickle
import pkg_resources
import shutil
import tempfile
import unittest
from DisFit import core
from numpy.testing import assert_allclose
from pandas.testing import assert_frame_equal

#Todo: write resimulations test

FIXTURES = pkg_resources.resource_filename('tests', 'fixtures')
SBML_PATH = os.path.join(FIXTURES, 'G2M_copasi.xml')
DATA_PATH = os.path.join(FIXTURES, 'G2M_copasi.csv')
JL_FILE_GOLD = os.path.join(FIXTURES, 'jl_file_gold.jl')
with open(JL_FILE_GOLD, 'r') as f:
    JL_CODE_GOLD = f.read()

class DisFitProblemTestCase(unittest.TestCase):
    
    def setUp(self):
        self.dirname = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.dirname)

    def test_constructor(self):
        problem = core.DisFitProblem(SBML_PATH, DATA_PATH)

    def test_sbml_path_setter(self):
        problem = core.DisFitProblem(SBML_PATH, DATA_PATH)
        problem.sbml_path = SBML_PATH
        self.assertEqual(problem.sbml_path, SBML_PATH)
        with self.assertRaises(ValueError):
            problem.sbml_path = 0
        with self.assertRaises(ValueError):
            problem.sbml_path = DATA_PATH
        with self.assertRaises(ValueError):
            problem.sbml_path = FIXTURES

    def test_data_path_setter(self):
        problem = core.DisFitProblem(SBML_PATH, DATA_PATH)
        problem.data_path = DATA_PATH
        self.assertEqual(problem.data_path, DATA_PATH)
        with self.assertRaises(ValueError): 
            problem.data_path = 0
        with self.assertRaises(ValueError):
            problem.data_path = SBML_PATH
        with self.assertRaises(ValueError):
            problem.data_path = FIXTURES
    
    def test_t_ratio_setter(self):
        problem = core.DisFitProblem(SBML_PATH, DATA_PATH)
        self.assertEqual(problem.t_ratio, 2)
        with self.assertRaises(ValueError):
            problem.t_ratio = 0
        with self.assertRaises(ValueError):
            problem.t_ratio = 'a'
        self.assertEqual(problem.julia_code, JL_CODE_GOLD)
        problem.t_ratio = 10
        self.assertNotEqual(problem.julia_code, JL_CODE_GOLD)

    def test_fold_change_setter(self):
        problem = core.DisFitProblem(SBML_PATH, DATA_PATH)
        self.assertEqual(problem.fold_change, 2)
        with self.assertRaises(ValueError):
            problem.fold_change = 1.0
        with self.assertRaises(ValueError):
            problem.fold_change = 'a'
        self.assertEqual(problem.julia_code, JL_CODE_GOLD)
        problem.fold_change = 10
        self.assertNotEqual(problem.julia_code, JL_CODE_GOLD)

    def test_n_starts_setter(self):
        problem = core.DisFitProblem(SBML_PATH, DATA_PATH)
        self.assertEqual(problem.n_starts, 1)
        with self.assertRaises(ValueError):
            problem.n_starts = 1.1
        with self.assertRaises(ValueError):
            problem.n_starts = 'a'
        self.assertEqual(problem.julia_code, JL_CODE_GOLD)
        problem.n_starts = 2
        self.assertNotEqual(problem.julia_code, JL_CODE_GOLD)

    def test_julia_code(self):
        problem = core.DisFitProblem(SBML_PATH, DATA_PATH)
        self.assertEqual(problem.julia_code, JL_CODE_GOLD)

    def test_write_jl_file(self):
        problem = core.DisFitProblem(SBML_PATH, DATA_PATH)
        problem.write_jl_file(path=os.path.join(self.dirname, 'jl_code.jl'))
        self.assertTrue(filecmp.cmp(os.path.join(self.dirname, 'jl_code.jl'), JL_FILE_GOLD))

    def test_optimize_results_plot(self):

        # test_optimize()
        problem = core.DisFitProblem(SBML_PATH, DATA_PATH)
        problem.optimize()

        results = problem.results
        with open(os.path.join(FIXTURES, 'results_gold.pickle'), 'rb') as f:
            results_gold = pickle.load(f)
        # self.assertEqual(problem.results, results_gold) # Todo: for some reason the output is none even if it shouldn't be)
        self.assertEqual(set(results.keys()), set(['x', 'x_best', 'states']))
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
        
        assert_frame_equal(results['x_best'], results_gold['x_best'])

        problem.write_results(path=os.path.join(self.dirname, 'results.xlsx'))
        self.assertTrue(os.path.isfile(os.path.join(self.dirname, 'results.xlsx')))

        # test_plot_results()
        problem.plot_results(path=os.path.join(self.dirname, 'plot.pdf'))
        problem.plot_results(path=os.path.join(self.dirname, 'plot.pdf'), variables=['Ensa', 'pEnsa'])
        self.assertTrue(os.path.isfile(os.path.join(self.dirname, 'plot.pdf')))