""" Test of DisFit.core
:Author: Paul F Lang <paul.lang@wolfson.ox.ac.uk>
:Date: 2020-04-16
:Copyright: 2020, Paul F Lang
:License: MIT
"""

import copy
import numpy as np
import os
import pandas as pd
import petab
import pickle
import pkg_resources
import re
import shutil
import tempfile
import unittest
from DisFit import core
import importlib
importlib.reload(core)
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
with open(os.path.join(FIXTURES, 'results_gold.pickle'), 'rb') as f:
    RESULTS_GOLD = pickle.load(f)


class Mock(core.DisFitProblem):
    def __init__(self):
        self._initialization = True


class DisFitProblemTestCase(unittest.TestCase):
    
    def setUp(self):
        self.dirname = tempfile.mkdtemp()


    def tearDown(self):
        shutil.rmtree(self.dirname)


    def test_constructor(self):
        problem = core.DisFitProblem(PETAB_YAML)


    def test_petab_yaml_dict(self):
        problem = Mock()
        problem._petab_yaml_dict = dict()
        self.assertTrue(isinstance(problem.petab_yaml_dict, dict))


    def test_petab_problem(self):
        problem = Mock()
        problem._petab_problem = petab.problem.Problem()
        self.assertTrue(isinstance(problem.petab_problem, petab.problem.Problem))


    def test_t_ratio_setter(self):
        problem = Mock()
        problem.t_ratio = 2.1
        self.assertEqual(problem.t_ratio, 2.1)
        with self.assertRaises(ValueError):
            problem.t_ratio = -1
        with self.assertRaises(ValueError):
            problem.t_ratio = 'a'

        # Test resetting
        problem = core.DisFitProblem(PETAB_YAML)
        self.assertEqual(problem.julia_code, JL_CODE_GOLD)
        problem.t_ratio = 1
        self.assertNotEqual(problem.julia_code, JL_CODE_GOLD)


    def test_n_starts_setter(self):
        problem = Mock()
        problem.n_starts = 2
        self.assertEqual(problem.n_starts, 2)
        with self.assertRaises(ValueError):
            problem.n_starts = 1.1
        with self.assertRaises(ValueError):
            problem.n_starts = 'a'

        # Test resetting
        problem = core.DisFitProblem(PETAB_YAML)
        self.assertEqual(problem.julia_code, JL_CODE_GOLD)
        problem.n_starts = 2
        self.assertNotEqual(problem.julia_code, JL_CODE_GOLD)


    def test_infer_ic_from_sbml_setter(self):
        problem = Mock()
        problem.infer_ic_from_sbml = True
        self.assertTrue(problem.infer_ic_from_sbml)
        with self.assertRaises(ValueError):
            problem.infer_ic_from_sbml = 1
        
        # Test resetting
        problem = core.DisFitProblem(PETAB_YAML)
        self.assertEqual(problem.julia_code, JL_CODE_GOLD)
        problem.n_starts = 2
        self.assertNotEqual(problem.julia_code, JL_CODE_GOLD)


    def test_julia_code(self):
        problem = Mock()
        problem._julia_code = JL_CODE_GOLD
        self.assertEqual(problem.julia_code, JL_CODE_GOLD)


    def test_results(self):
        problem = Mock()
        problem._results = RESULTS_GOLD
        self.assertTrue(isinstance(RESULTS_GOLD, dict))


    def test_set_petab_problem(self):
        problem = Mock()
        problem._set_petab_problem(PETAB_YAML)
        self.assertTrue(isinstance(problem._petab_problem, petab.problem.Problem))
        self.assertTrue(isinstance(problem.petab_yaml_dict, dict))
        self.assertTrue(len(problem._condition2index), 4)
        self.assertEqual(problem._n_conditions, 4)
        self.assertEqual(len(problem._condition_specific_pars), 4)
        self.assertEqual(len(problem._condition_specific_pars['iWee_0']), 4)
        self.assertEqual(len(problem._global_pars), 24)


    def test_check_for_not_implemented_features(self):
        problem = Mock()

        petab_problem = petab.problem.Problem()                                                                                                                                                                          
        petab_problem = petab_problem.from_yaml(PETAB_YAML)
        problem._check_for_not_implemented_features(petab_problem)

        petab_problem_wrong = copy.deepcopy(petab_problem)
        petab_problem_wrong.measurement_df['preequilibrationConditionId'] = np.empty(len(petab_problem_wrong.measurement_df.index))
        with self.assertRaises(NotImplementedError) as c:
            problem._check_for_not_implemented_features(petab_problem_wrong)
        self.assertTrue('Preequilibration' in str(c.exception))

        petab_problem_wrong = copy.deepcopy(petab_problem)
        petab_problem_wrong.measurement_df['time'] = np.inf
        with self.assertRaises(NotImplementedError) as c:
            problem._check_for_not_implemented_features(petab_problem_wrong)
        print(c.exception)
        self.assertTrue('steady state' in str(c.exception))

        petab_problem_wrong = copy.deepcopy(petab_problem)
        petab_problem_wrong.measurement_df.loc[2, 'time'] = 1
        with self.assertRaises(NotImplementedError) as c:
            problem._check_for_not_implemented_features(petab_problem_wrong)
        self.assertTrue('Measurement time points' in str(c.exception))


    def test_sort_condition_df_problem(self):
        problem = Mock()

        petab_problem = petab.problem.Problem()                                                                                                                                                                          
        petab_problem = petab_problem.from_yaml(PETAB_YAML)
        petab_problem.measurement_df = petab_problem.measurement_df.iloc[np.arange(-1, len(petab_problem.measurement_df.index)-1)]
        petab_problem = problem._sort_condition_df_problem(petab_problem)
        print(petab_problem.condition_df)
        self.assertEqual(petab_problem.condition_df.index[0], 'pGw_weak')


    def test_get_translation_vars(self):
        problem = Mock()

        petab_problem = petab.problem.Problem()                                                                                                                                                                          
        petab_problem = petab_problem.from_yaml(PETAB_YAML)
        yaml_dict, condition2index, n_conditions, condition_specific_pars, global_pars =\
            problem._get_translation_vars(PETAB_YAML, petab_problem)
        self.assertTrue(isinstance(yaml_dict, dict))
        self.assertTrue(len(condition2index), 4)
        self.assertEqual(n_conditions, 4)
        self.assertEqual(len(condition_specific_pars), 4)
        self.assertEqual(len(condition_specific_pars['iWee_0']), 4)
        self.assertEqual(len(global_pars), 24)


    def test_write_jl_file(self):



        problem = core.DisFitProblem(PETAB_YAML)
        problem.write_jl_file(path=os.path.join(self.dirname, 'jl_code.jl'))
        with open(os.path.join(self.dirname, 'jl_code.jl'), 'r') as f:
            jl_code = f.read()
        self.assertEqual(jl_code, JL_CODE_GOLD)






    def test_plot(self):
        problem = core.DisFitProblem(PETAB_YAML)
        with open(os.path.join(FIXTURES, 'results_gold.pickle'), 'rb') as f:
            results_gold = pickle.load(f)
        problem._results = results_gold
        problem.plot_results('wt')

    # def test_optimize_results_plot(self):
    #     # test_optimize()
    #     problem = core.DisFitProblem(PETAB_YAML, t_ratio=1.99999)
    #     problem.optimize()
    #     results = problem.results

    #     with open(os.path.join(FIXTURES, 'results_gold.pickle'), 'rb') as f:
    #         results_gold = pickle.load(f)

    #     # self.assertEqual(problem.results, results_gold) # Todo: for some reason the output is none even if it shouldn't be)
    #     self.assertEqual(set(results.keys()), set(['parameters', 'par_best', 'species', 'observables']))
    #     self.assertEqual(results['parameters'].keys(), results_gold['parameters'].keys())
    #     for i_iter in results_gold['parameters'].keys():
    #         self.assertEqual(results['parameters'][i_iter].keys(), results_gold['parameters'][i_iter].keys())
    #         for param in results_gold['parameters'][i_iter].keys():
    #             assert_allclose(results['parameters'][i_iter][param], results_gold['parameters'][i_iter][param])

    #     self.assertEqual(results['species'].keys(), results_gold['species'].keys())
    #     for i_iter in results_gold['species'].keys():
    #         self.assertEqual(results['species'][i_iter].keys(), results_gold['species'][i_iter].keys())
    #         for specie in results_gold['species'][i_iter].keys():
    #             assert_allclose(results['species'][i_iter][specie],
    #                 results_gold['species'][i_iter][specie], rtol=1e-05, atol=1e-08)

    #     self.assertEqual(results['observables'].keys(), results_gold['observables'].keys())
    #     for i_iter in results_gold['observables'].keys():
    #         self.assertEqual(results['observables'][i_iter].keys(), results_gold['observables'][i_iter].keys())
    #         for observable in results_gold['observables'][i_iter].keys():
    #             assert_allclose(results['observables'][i_iter][observable],
    #                 results_gold['observables'][i_iter][observable], rtol=1e-05, atol=1e-08)
        
    #     assert_frame_equal(results['par_best'], results_gold['par_best'])

    #     problem.write_results(path=os.path.join(self.dirname, 'results.xlsx'))
    #     self.assertTrue(os.path.isfile(os.path.join(self.dirname, 'results.xlsx')))

    #     # test_plot_results()
    #     problem.plot_results('wt', path=os.path.join(self.dirname, 'plot.pdf'))
    #     problem.plot_results('wt', path=os.path.join(self.dirname, 'plot.pdf'), observables=['obs_Ensa', 'obs_pEnsa'])
    #     self.assertTrue(os.path.isfile(os.path.join(self.dirname, 'plot.pdf')))


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
# # FIXTURES = os.path.join('/media/sf_DPhil_Project/Project07_Parameter Fitting/df_software',
# #     'petab_test_suite', 'cases')
# # PETAB_YAML = os.path.join(FIXTURES, '0003', '_0003.yaml')
# # FIXTURES = os.path.join('/media/sf_DPhil_Project/Project07_Parameter Fitting/df_software',
# #     'Benchmark-Models-PEtab', 'Benchmark-Models')
# # PETAB_YAML = os.path.join(FIXTURES, 'Borghans_BiophysChem1997', 'Borghans_BiophysChem1997.yaml')

# # jl_file_gold = os.path.join(FIXTURES, 'jl_file_gold.jl')
# # with open(jl_file_gold, 'r') as f:
# #     JL_CODE_GOLD = f.read()
# # JL_CODE_GOLD = re.sub('/media/sf_DPhil_Project/Project07_Parameter Fitting/df_software/DisFit/tests/fixtures',
# #     FIXTURES, JL_CODE_GOLD)

# problem = core.DisFitProblem(PETAB_YAML)
# problem.write_jl_file()
# problem.optimize()
# problem.plot_results('wt', path='plot.pdf')
# problem.write_results()
# # problem.results['par_best']


# import pickle
# with open('results_gold.pickle', 'wb') as f:
#     pickle.dump(problem.results, f)