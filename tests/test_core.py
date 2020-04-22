""" Test of DisFit.core

:Author: Paul F Lang <paul.lang@wolfson.ox.ac.uk>
:Date: 2020-04-16
"""

import os
import pandas as pd
import pkg_resources
import shutil
import tempfile
import unittest
from DisFit import core

#Todo: write resimulations test

FIXTURES = pkg_resources.resource_filename('tests', 'fixtures')
SBML_PATH = os.path.join(FIXTURES, 'G2M_copasi.xml')
DATA_PATH = os.path.join(FIXTURES, 'G2M_copasi.csv')
OPTIMIZER_PATH = os.path.join(FIXTURES, 'optimizer.jl')
with open(OPTIMIZER_PATH) as f:
    OPTIMIZER_GOLD = f.read()
X_BEST_PATH = os.path.join(FIXTURES, 'x_best.jl')
with open(X_BEST_PATH) as f:
    X_BEST_GOLD = f.read()
VALUES_PATH = os.path.join(FIXTURES, 'values.jl')
with open(VALUES_PATH) as f:
    VALUES_GOLD = f.read()

class DisFitProblemTestCase(unittest.TestCase):
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
        self.assertEqual(problem.julia_code['optimizer'], OPTIMIZER_GOLD)
        problem.t_ratio = 10
        self.assertNotEqual(problem.julia_code['optimizer'], OPTIMIZER_GOLD)

    def test_fold_change_setter(self):
        problem = core.DisFitProblem(SBML_PATH, DATA_PATH)
        self.assertEqual(problem.fold_change, 2)
        with self.assertRaises(ValueError):
            problem.fold_change = 1.0
        with self.assertRaises(ValueError):
            problem.fold_change = 'a'
        self.assertEqual(problem.julia_code['optimizer'], OPTIMIZER_GOLD)
        problem.fold_change = 10
        self.assertNotEqual(problem.julia_code['optimizer'], OPTIMIZER_GOLD)

    def test_julia_code(self):
        problem = core.DisFitProblem(SBML_PATH, DATA_PATH)
        self.assertEqual(problem.julia_code['optimizer'], OPTIMIZER_GOLD)
        self.assertEqual(problem.julia_code['x_best'], X_BEST_GOLD)
        self.assertEqual(problem.julia_code['values'], VALUES_GOLD)

    def test_write_jl_file(self):
        def setUp(self):
            self.filename = tempfile.mkstemp()

        def tearDown(self):
            shutil.rmtree(self.filename)

        def test_write_jl_file(self):
            jl_file_gold = os.path.join(FIXTURES, jl_file_gold.jl)
            problem = core.DisFitProblem(SBML_PATH, DATA_PATH)
            problem.write_jl_file(filename)
            self.assertTrue(filecmp(filename, jl_file_gold))

    def test_optimize_results_plot(self):
        def setUp(self):
            self.dirname = tempfile.mkdtemp()

        def tearDown(self):
            shutil.rmtree(self.dirname)

        # test_optimize()
        problem = core.DisFitProblem(SBML_PATH, DATA_PATH)
        problem.optimize()
        with open(os.path.join(FIXTURES, 'results_gold.pickle'), 'rb') as f:
            results_gold = pickle.load(f)

        self.assertEqual(problem.results['optimizer_log'], results_gold['optimizer_log']) # Todo: for some reason the output is none even if it shouldn't be)
        self.assertEqual(problem.results['parameters'], results_gold['parameters'])
        self.assertEqual(problem.results['values'], results_gold['values'])

        # test_plot_results()
        problem.plot_results(path=os.path.join(self.dirname, 'plot.pdf'))
        problem.plot_results(path=os.path.join(self.dirname, 'plot.pdf'), variables=['Ensa', 'pENSA'])