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
from julia.api import Julia
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
        self._jl = Julia(compiled_modules=False)


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
        problem = Mock()
        problem._julia_code = JL_CODE_GOLD
        problem.write_jl_file(path=os.path.join(self.dirname, 'jl_code.jl'))
        with open(os.path.join(self.dirname, 'jl_code.jl'), 'r') as f:
            jl_code = f.read()
        self.assertEqual(jl_code, JL_CODE_GOLD)


    def test_optimize(self):
        problem = core.DisFitProblem(PETAB_YAML)
        results = problem.optimize()

        print(problem._best_iter)
        self.assertEqual(problem._best_iter, '1')

        self.assertEqual(set(results.keys()), set(['par_best', 'species', 'observables', 'fval', 'chi2']))
        
        assert_frame_equal(results['par_best'], RESULTS_GOLD['par_best'])
        assert_frame_equal(results['species'], RESULTS_GOLD['species'])
        assert_frame_equal(results['observables'], RESULTS_GOLD['observables'])
        self.assertAlmostEqual(results['fval'], RESULTS_GOLD['fval'])
        self.assertAlmostEqual(results['chi2'], RESULTS_GOLD['chi2'])

        self.assertTrue(isinstance(problem.petab_problem.simulation_df, pd.DataFrame))


    def test_get_param_ratios(self):
        problem = Mock()
        problem._best_iter = '1'
        problem._n_conditions = 4
        problem._petab_problem = petab.problem.Problem()                                                                                                                                                                          
        problem._petab_problem = problem._petab_problem.from_yaml(PETAB_YAML)

        par_dict = {'1':
            {'fB55_pGw_weak': 1.250000005, 'fB55_wt': 1.250000005, 'fB55_iWee': 1.1250000039999999,
            'kPhCdc25': 0.9955692429176523, 'kDpEnsa': 0.04984394161435494, 'kDpGw2': 9.993552511531249,
            'kCdc25_1': 0.10303120213630562, 'kPhEnsa_iWee': 0.125, 'kDpGw1': 0.24845896113749727,
            'kDpCdc25': 9.841409443285976, 'kPhEnsa_pGw_weak': 0.11249999999999999,
            'kPhWee': 1.013150760597791, 'kDpWee': 10.052516838989787, 'kPhGw': 0.995923588764344,
            'kPhEnsa_wt': 0.125, 'kPhEnsa_Cb_low': 0.125, 'kWee1': 0.009879213854455046,
            'fCb': 2.0000399837133416, 'fB55_Cb_low': 1.375000006, 'kAspEB55': 52.05738881531509,
            'jiWee': 0.11443890616434232, 'kWee2': 1.0164397373870508, 'kCdc25_2': 0.9133917145270964,
            'kDipEB55': 0.0034194365965829255}
            }

        par_best = problem._get_param_ratios(par_dict)
        assert_frame_equal(par_best, RESULTS_GOLD['par_best'])


    def test_results_to_frame(self):
        problem = Mock()
        problem._petab_problem = petab.problem.Problem()                                                                                                                                                                          
        problem._petab_problem = problem._petab_problem.from_yaml(PETAB_YAML)
        problem._best_iter = '1'
        problem._t_ratio = 1
        res_dict_gold = {'speciesId': ['a', 'a', 'b', 'b', 'a', 'a', 'b', 'b'],
            'simulationConditionId': ['c0', 'c0', 'c0', 'c0', 'c1', 'c1', 'c1', 'c1'],
            'time': [0., 1., 0., 1., 0., 1., 0., 1.], 'simulation': [1,2,3,4,5,6,7,8]}
        simulation_dict = {'1': {'a': [[1,2], [5,6]], 'b': [[3,4], [7,8]]}}
        problem._condition2index = {'c0': 0, 'c1': 1}
        problem._petab_problem.measurement_df = pd.DataFrame(res_dict_gold).sort_values(['simulationConditionId', 'speciesId', 'time'])
        df = problem._results_to_frame(simulation_dict, variable_type='speciesId')

        assert_frame_equal(df, problem._petab_problem.measurement_df)


    def test_set_simulation_df(self):
        problem = Mock()
        problem._petab_problem = petab.problem.Problem()                                                                                                                                                                          
        problem._petab_problem = problem._petab_problem.from_yaml(PETAB_YAML)

        problem._petab_problem.measurement_df = pd.DataFrame({'observableId': ['obs_a', 'obs_a'],
            'simulationConditionId': ['c0', 'c0'], 'time': [0, 10], 'measurement': [0.7, 0.1]})
        df = pd.DataFrame({'observableId': ['obs_a'], 'simulationConditionId': ['c0'],
            'time': [10.1], 'measurement': [1]})
        problem._results = {'observables': problem._petab_problem.measurement_df.append(df).reset_index(drop=True)}
        print(problem.results)

        problem._set_simulation_df()

        simulation_df_gold = copy.deepcopy(problem._petab_problem.measurement_df)
        assert_frame_equal(problem._petab_problem.simulation_df, simulation_df_gold)


    def test_write_results(self):
        problem = Mock()
        problem._results = RESULTS_GOLD
        problem._condition2index = {'wt': 0, 'iWee': 1, 'Cb_low': 2, 'pGw_weak': 3}

        problem.write_results(path=os.path.join(self.dirname, 'results_1.xlsx'))
        self.assertTrue(os.path.isfile(os.path.join(self.dirname, 'results_1.xlsx')))

        problem.write_results(path=os.path.join(self.dirname, 'results_2.xlsx'), df_format='long')
        self.assertTrue(os.path.isfile(os.path.join(self.dirname, 'results_2.xlsx')))


    def test_plot(self):
        problem = Mock()
        problem._t_ratio = 2
        problem._petab_problem = petab.problem.Problem()                                                                                                                                                                          
        problem._petab_problem = problem._petab_problem.from_yaml(PETAB_YAML)
        problem._results = RESULTS_GOLD

        problem.plot_results('wt', path=os.path.join(self.dirname, 'plot_1.pdf'))
        self.assertTrue(os.path.isfile(os.path.join(self.dirname, 'plot_1.pdf')))

        problem.plot_results('wt', path=os.path.join(self.dirname, 'plot_2.pdf'), observables=['obs_Ensa', 'obs_pEnsa'])
        self.assertTrue(os.path.isfile(os.path.join(self.dirname, 'plot_2.pdf')))


    def test_set_julia_code(self):
        problem = core.DisFitProblem(PETAB_YAML)
        self.assertEqual(problem.julia_code, JL_CODE_GOLD)


    def test_write_overrides(self):
        problem = Mock()
        problem._petab_problem = petab.problem.Problem()                                                                                                                                                                          
        problem._petab_problem = problem._petab_problem.from_yaml(PETAB_YAML)

        code, set_of_params = problem._write_overrides('', 'observable')
        self.assertTrue(code in JL_CODE_GOLD)
        self.assertEqual(set_of_params, set())

        code, set_of_params = problem._write_overrides('', 'noise')
        self.assertTrue(code in JL_CODE_GOLD)
        self.assertEqual(set_of_params, set())
        
        sdfs # Todo: add a problem with noise parameters here.


    # def test_plot(self):
    #     problem = core.DisFitProblem(PETAB_YAML)
    #     with open(os.path.join(FIXTURES, 'results_gold.pickle'), 'rb') as f:
    #         results_gold = pickle.load(f)
    #     problem._results = results_gold
    #     problem.plot_results('wt')

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