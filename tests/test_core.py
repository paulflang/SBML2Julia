""" Test of DisFit.core

:Author: Paul F Lang <paul.lang@wolfson.ox.ac.uk>
:Date: 2020-04-16
"""

import os
import pkg_resources
import shutil
import tempfile
import unittest
from DisFit import core

FIXTURES = pkg_resources.resource_filename('tests', 'fixtures')
SBML_PATH = os.path.join(FIXTURES, 'G2M_copasi.xml')
DATA_PATH = os.path.join(FIXTURES, 'G2M_copasi.csv')
OPTIMIZER_PATH = os.path.join(FIXTURES, 'Optimizer.jl')
with open(OPTIMIZER_PATH) as f:
    OPTIMIZER_GOLD = f.read()
PARAMETERS_PATH = os.path.join(FIXTURES, 'Parameters.jl')
with open(PARAMETERS_PATH) as f:
    PARAMETERS_GOLD = f.read()
VALUES_PATH = os.path.join(FIXTURES, 'Values.jl')
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
            problem.sbml_path = FIXTURES
    
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
        self.assertEqual(problem.julia_code['parameters'], PARAMETERS_GOLD)
        self.assertEqual(problem.julia_code['values'], VALUES_GOLD)

    def test_optimize(self):
        problem = core.DisFitProblem(SBML_PATH, DATA_PATH)
        results = problem.optimize()
        self.assertEqual(results['parameters'], problem.results['parameters']) # Todo: this test says very little
        self.assertEqual(results['values'], problem.results['values'])