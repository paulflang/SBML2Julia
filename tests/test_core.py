""" Test of sbml2julia.core
:Author: Paul F Lang <paul.lang@wolfson.ox.ac.uk>
:Date: 2020-04-16
:Copyright: 2020, Paul F Lang
:License: MIT
"""

import copy
import importlib
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
from julia.api import Julia
from pandas.testing import assert_frame_equal
from sbml2julia import core
importlib.reload(core)


FIXTURES = pkg_resources.resource_filename('tests', 'fixtures')
PETAB_YAML = os.path.join(FIXTURES, '0015_objectivePrior', '_0015_objectivePrior.yaml')
JL_FILE_GOLD = os.path.join(FIXTURES, 'jl_file_gold.jl')
with open(JL_FILE_GOLD, 'r') as f:
    JL_CODE_GOLD = f.read()
JL_CODE_GOLD = re.sub('/media/sf_DPhil_Project/Project07_Parameter Fitting/'
                      'df_software/sbml2julia/tests/fixtures', FIXTURES, JL_CODE_GOLD)
with open(os.path.join('.', 'substituted_code.jl'), 'w') as f:
    f.write(JL_CODE_GOLD)
with open(os.path.join(FIXTURES, 'results_gold.pickle'), 'rb') as f:
    RESULTS_GOLD = pickle.load(f)


class Mock(core.SBML2JuliaProblem):
    def __init__(self):
        self._initialization = True
        self._jl = Julia(compiled_modules=False)


class SBML2JuliaProblemTestCase(unittest.TestCase):

    def setUp(self):
        self.dirname = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.dirname)

    def test_constructor(self):
        _ = core.SBML2JuliaProblem(PETAB_YAML)

    def test_petab_yaml_dict(self):
        problem = Mock()
        problem._petab_yaml_dict = dict()
        self.assertTrue(isinstance(problem.petab_yaml_dict, dict))

    def test_petab_problem(self):
        problem = Mock()
        problem._petab_problem = petab.problem.Problem()
        self.assertTrue(isinstance(problem.petab_problem, petab.problem.Problem))

    def test_t_steps_setter(self):
        problem = Mock()
        problem._set_petab_problem(PETAB_YAML)
        problem.t_steps = None
        self.assertEqual(problem.t_steps, 101)
        problem.t_steps = 3
        self.assertEqual(problem.t_steps, 3)
        with self.assertRaises(ValueError):
            problem.t_steps = -1
        with self.assertRaises(ValueError):
            problem.t_steps = 1.1
        with self.assertRaises(ValueError):
            problem.t_steps = 'a'

    def test_n_starts_setter(self):
        problem = Mock()
        problem.n_starts = 2
        self.assertEqual(problem.n_starts, 2)
        with self.assertRaises(ValueError):
            problem.n_starts = 1.1
        with self.assertRaises(ValueError):
            problem.n_starts = 'a'

    def test_infer_ic_from_sbml_setter(self):
        problem = Mock()
        problem.infer_ic_from_sbml = True
        self.assertTrue(problem.infer_ic_from_sbml)
        with self.assertRaises(ValueError):
            problem.infer_ic_from_sbml = 1

        # Test resetting
        problem = core.SBML2JuliaProblem(PETAB_YAML)
        problem.write_jl_file(path=os.path.join('.', 'test_infer_ic_from_sbml_setter.jl'))
        self.assertEqual(problem.julia_code, JL_CODE_GOLD)
        # problem.infer_ic_from_sbml = True # Todo: add species to SBML for this test case
        # self.assertNotEqual(problem.julia_code, JL_CODE_GOLD)

    def test_optimizer_options_setter(self):
        problem = Mock()
        problem.optimizer_options = {'linear_solver': 'MA27'}
        self.assertEqual(problem.optimizer_options['linear_solver'], 'MA27')
        with self.assertRaises(ValueError):
            problem.optimizer_options = 1

        # Test resetting
        problem = core.SBML2JuliaProblem(PETAB_YAML)
        self.assertEqual(problem.julia_code, JL_CODE_GOLD)
        problem.optimizer_options = {'linear_solver': 'MA27'}
        self.assertNotEqual(problem.julia_code, JL_CODE_GOLD)

    def test_custom_code_dict_setter(self):
        problem = Mock()
        problem._julia_code = 'abc\n123'
        problem.custom_code_dict = {'abc': 'abcd'}
        self.assertEqual(problem.custom_code_dict['abc'], 'abcd')
        with self.assertRaises(ValueError):
            problem.custom_code_dict = 1
        with self.assertRaises(ValueError):
            problem.custom_code_dict = {1: '1234'}
        with self.assertRaises(ValueError):
            problem.custom_code_dict = {'abc': 1}

        # Test resetting
        problem = core.SBML2JuliaProblem(PETAB_YAML)
        self.assertEqual(problem.julia_code, JL_CODE_GOLD)
        problem.optimizer_options = {'# Write global parameters': '# Write global parameters1'}
        self.assertNotEqual(problem.julia_code, JL_CODE_GOLD)

    def test_julia_code(self):
        problem = Mock()
        problem._julia_code = JL_CODE_GOLD
        self.assertEqual(problem.julia_code, JL_CODE_GOLD)

    def test_results(self):
        problem = Mock()
        problem._results = RESULTS_GOLD
        self.assertTrue(isinstance(RESULTS_GOLD, dict))

    def test_import_julia_code(self):
        problem = Mock()
        problem.import_julia_code(JL_FILE_GOLD)
        problem._julia_code = re.sub('/media/sf_DPhil_Project/Project07_Parameter Fitting/'
                                     'df_software/sbml2julia/tests/fixtures',
                                     FIXTURES, problem._julia_code)
        self.assertEqual(problem.julia_code, JL_CODE_GOLD)

    def test_set_petab_problem(self):
        problem = Mock()
        problem._set_petab_problem(PETAB_YAML)
        self.assertTrue(isinstance(problem._petab_problem, petab.problem.Problem))
        self.assertTrue(isinstance(problem._petab_yaml_dict, dict))
        self.assertTrue(len(problem._condition2index), 3)
        self.assertEqual(problem._n_conditions, 3)
        self.assertEqual(len(problem._condition_specific_pars), 1)
        self.assertEqual(len(problem._condition_specific_pars['k1']), 3)
        self.assertEqual(len(problem._global_pars), 9)

    def test_check_for_not_implemented_features(self):
        problem = Mock()

        petab_problem = petab.problem.Problem()
        petab_problem = petab_problem.from_yaml(PETAB_YAML)
        problem._check_for_not_implemented_features(petab_problem)

        petab_problem_wrong = copy.deepcopy(petab_problem)
        petab_problem_wrong.measurement_df['time'] = np.inf
        with self.assertRaises(NotImplementedError) as c:
            problem._check_for_not_implemented_features(petab_problem_wrong)
        print(c.exception)
        self.assertTrue('steady state' in str(c.exception))

        petab_problem_wrong = copy.deepcopy(petab_problem)
        petab_problem_wrong.measurement_df['preequilibrationConditionId'] = 1*['p1']+43*['p2']
        with self.assertRaises(NotImplementedError) as c:
            problem._check_for_not_implemented_features(petab_problem_wrong)
        self.assertTrue('with <= 1 preequilibrationConditionIds' in str(c.exception))

    def test_sort_condition_df_problem(self):
        problem = Mock()

        petab_problem = petab.problem.Problem()
        petab_problem = petab_problem.from_yaml(PETAB_YAML)
        print(petab_problem.condition_df)
        petab_problem.measurement_df = petab_problem.measurement_df.iloc[
            np.arange(-1, len(petab_problem.measurement_df.index)-1)]
        print(petab_problem.measurement_df['simulationConditionId'])
        petab_problem = problem._sort_condition_df_problem(petab_problem)
        print(petab_problem.condition_df)
        self.assertEqual(petab_problem.condition_df.index[0], 'c1')

    def test_get_translation_vars(self):
        problem = Mock()

        petab_problem = petab.problem.Problem()
        petab_problem = petab_problem.from_yaml(PETAB_YAML)
        problem._petab_problem = petab_problem
        yaml_dict, condition2index, j_to_parameters, n_conditions, condition_specific_pars,\
            global_pars = problem._get_translation_vars(PETAB_YAML, petab_problem)
        self.assertTrue(isinstance(yaml_dict, dict))
        print(condition2index)
        self.assertTrue(len(condition2index), 3)

        self.assertEqual(j_to_parameters[0], [1, 2])
        self.assertEqual(j_to_parameters[1], ['', 3])

        self.assertEqual(n_conditions, 3)
        self.assertEqual(len(condition_specific_pars), 1)
        self.assertEqual(len(global_pars), 9)

    def test_write_jl_file(self):
        problem = Mock()
        problem._julia_code = JL_CODE_GOLD
        problem.write_jl_file(path=os.path.join(self.dirname, 'jl_code.jl'))
        with open(os.path.join(self.dirname, 'jl_code.jl'), 'r') as f:
            jl_code = f.read()
        self.assertEqual(jl_code, JL_CODE_GOLD)

    def test_insert_custom_code(self):
        problem = Mock()
        problem._julia_code = 'abc\n123'
        problem.insert_custom_code({'abc': 'abcd', '123': '01234'})
        self.assertEqual(problem._julia_code, 'abcd\n01234')

    def test_optimize(self):

        def _assert_frame_almost_equal(df1, df2):
            return (np.allclose(df1.select_dtypes(exclude=[object]), df2.select_dtypes(
                exclude=[object])) & df1.select_dtypes(include=[object]).equals(
                    df2.select_dtypes(include=[object])))

        problem = core.SBML2JuliaProblem(PETAB_YAML, n_starts=3)
        results = problem.optimize()

        self.assertTrue(problem._best_iter in ['1', '2', '3'])
        self.assertEqual(set(results.keys()),
                         set(['par_est', 'species', 'observables', 'fval', 'chi2']))

        print(results['par_est'])
        _assert_frame_almost_equal(results['par_est'], RESULTS_GOLD['par_est'])
        _assert_frame_almost_equal(results['species'], RESULTS_GOLD['species'])
        _assert_frame_almost_equal(results['observables'], RESULTS_GOLD['observables'])
        self.assertAlmostEqual(results['fval'], RESULTS_GOLD['fval'], delta=1.)
        self.assertAlmostEqual(results['chi2'], RESULTS_GOLD['chi2'], delta=1.)

        self.assertTrue(isinstance(problem.petab_problem.simulation_df, pd.DataFrame))

    def test_prior_code(self):

        problem = core.SBML2JuliaProblem(PETAB_YAML)
        problem._global_pars = {k: (0 if k != 'a0' else 1) for k in problem._global_pars.keys()}
        problem.petab_problem.parameter_df['objectivePriorParameters'].iloc[0] = '2; 0.1'
        problem._set_julia_code()
        results_shifted = problem.optimize()
        self.assertTrue(results_shifted['par_est'].loc['a0', 'estimatedValue']
                        > RESULTS_GOLD['par_est'].loc['a0', 'estimatedValue'])

        problem = core.SBML2JuliaProblem(PETAB_YAML)
        problem._global_pars = {k: (0 if k != 'b0' else 1) for k in problem._global_pars.keys()}
        problem.petab_problem.parameter_df['objectivePriorParameters'].iloc[1] = '0.05; 1'
        problem._set_julia_code()
        results_shifted = problem.optimize()
        self.assertTrue(results_shifted['par_est'].loc['b0', 'estimatedValue']
                        > RESULTS_GOLD['par_est'].loc['b0', 'estimatedValue'])

        problem = core.SBML2JuliaProblem(PETAB_YAML)
        problem._global_pars = {k: (0 if k != 'k1_free' else 1)
                                for k in problem._global_pars.keys()}
        problem.petab_problem.parameter_df['objectivePriorParameters'].iloc[2] = '1; 0.1'
        problem._set_julia_code()
        results_shifted = problem.optimize()
        self.assertTrue(results_shifted['par_est'].loc['k1_free', 'estimatedValue']
                        > RESULTS_GOLD['par_est'].loc['k1_free', 'estimatedValue'])

        problem = core.SBML2JuliaProblem(PETAB_YAML)
        problem._global_pars = {k: (0 if k != 'k2' else 1) for k in problem._global_pars.keys()}
        problem._petab_problem.parameter_df['objectivePriorParameters'].iloc[3] = '-0.4; 0.01'
        problem._set_julia_code()
        results_shifted = problem.optimize()
        self.assertTrue(results_shifted['par_est'].loc['k2', 'estimatedValue']
                        > RESULTS_GOLD['par_est'].loc['k2', 'estimatedValue'])

    def test_get_param_ratios(self):
        problem = Mock()
        problem._best_iter = '1'
        problem._n_conditions = 3
        problem._petab_problem = petab.problem.Problem()
        problem._petab_problem = problem._petab_problem.from_yaml(PETAB_YAML)

        par_dict = {'1':
                    {'a0': 1, 'b0': 0, 'k1_free': 0.8, 'k2': 0.6, 'noise_A1': 0.05,
                     'noise_A2': 0.1, 'noise_B': 0.1, 'offset_B': 0.5, 'scaling_B': 1}
                    }

        par_est = problem._get_param_ratios(par_dict)
        par_gold = problem._petab_problem.parameter_df.copy()
        par_gold['estimatedValue'] = [1, 0, 0.8, 0.6, 0.05, 0.1, 0.1, 0.5, 1]
        par_gold['estimate_to_nominal'] = [1., 'NA (diff=0)', 1., 1., 1., 1., 1., 1., 0.5]
        assert_frame_equal(par_est, par_gold)

    def test_results_to_frame(self):
        problem = Mock()
        problem._petab_problem = petab.problem.Problem()
        problem._petab_problem = problem._petab_problem.from_yaml(PETAB_YAML)
        problem._best_iter = '1'
        problem._t_steps = 11
        problem._j_to_parameters = ([1, 2], ['', 3])
        res_dict_gold = {'speciesId': 11*['A']+11*['B']+11*['A']+11*['B'],
                         'simulationConditionId': 22*['c0'] + 22*['c1'],
                         'time': 4*[float(i) for i in range(0, 11)],
                         'simulation': list(range(1, 45))}

        simulation_dict = {'A': {1: list(range(1, 12)), 2: list(range(23, 34))},
                           'B': {1: list(range(12, 23)), 2: list(range(34, 45))}}

        problem._condition2index = {'c0': 0, 'c1': 1}
        problem._petab_problem.measurement_df = pd.DataFrame(res_dict_gold).sort_values(
            ['simulationConditionId', 'speciesId', 'time'])
        df = problem._results_to_frame(simulation_dict, variable_type='speciesId')

        assert_frame_equal(df, problem._petab_problem.measurement_df)

        with self.assertRaises(ValueError):
            df = problem._results_to_frame(simulation_dict, variable_type=1)

    def test_write_results(self):
        problem = Mock()
        problem._results = RESULTS_GOLD
        problem._condition2index = {'c0': 0, 'c1': 1, 'p1': 2}
        problem._j_to_parameters = ([1, 2], ['', 3])

        problem.write_results(path=os.path.join(self.dirname, 'results_1.xlsx'))
        self.assertTrue(os.path.isfile(os.path.join(self.dirname, 'results_1.xlsx')))

        problem.write_results(path=os.path.join(self.dirname, 'results_2.xlsx'), df_format='wide')
        self.assertTrue(os.path.isfile(os.path.join(self.dirname, 'results_2.xlsx')))

        problem.write_results(path=os.path.join(self.dirname, 'results'))
        self.assertTrue(os.path.isfile(os.path.join(self.dirname, 'results', 'species.tsv')))
        self.assertTrue(os.path.isfile(os.path.join(self.dirname, 'results', 'observables.tsv')))
        self.assertTrue(os.path.isfile(os.path.join(self.dirname, 'results', 'parameters.tsv')))

        with self.assertRaises(ValueError):
            problem.write_results(path=os.path.join(self.dirname, 'results.docx'))

    def test_write_optimized_parameter_table(self):
        problem = Mock()
        problem._petab_problem = petab.problem.Problem()
        problem._petab_problem = problem._petab_problem.from_yaml(PETAB_YAML)
        problem._petab_dirname = os.path.dirname(PETAB_YAML)
        problem._results = RESULTS_GOLD

        problem.write_optimized_parameter_table()
        out_path = os.path.join(FIXTURES, '0015_objectivePrior', 'post_fit_parameters.tsv')
        self.assertTrue(os.path.isfile(out_path))
        os.remove(out_path)

    def test_plot(self):
        problem = Mock()
        problem._t_steps = 2
        problem._petab_problem = petab.problem.Problem()
        problem._petab_problem = problem._petab_problem.from_yaml(PETAB_YAML)
        problem._results = RESULTS_GOLD
        problem._condition2index = {'c0': 0, 'c1': 1}
        problem._obs_to_conditions = {'obs_a': [0, 1], 'obs_b': [0, 1]}

        problem.plot_results('c0', path=os.path.join(self.dirname, 'plot_1.pdf'))
        self.assertTrue(os.path.isfile(os.path.join(self.dirname, 'plot_1.pdf')))

        problem.plot_results('c0', path=os.path.join(self.dirname, 'plot_2.pdf'),
                             observables=['obs_a'])
        self.assertTrue(os.path.isfile(os.path.join(self.dirname, 'plot_2.pdf')))

    def test_set_julia_code(self):
        problem = core.SBML2JuliaProblem(PETAB_YAML)
        problem.write_jl_file(path=os.path.join('.', 'test_set_julia_code.jl'))
        self.assertEqual(problem.julia_code, JL_CODE_GOLD)

    def test_write_overrides(self):  # Todo: add a problem with noise parameters here.
        problem = Mock()
        problem._petab_problem = petab.problem.Problem()
        problem._petab_problem = problem._petab_problem.from_yaml(PETAB_YAML)
        problem._obs_to_conditions = {obs: [1, 2, 3, 4]
                                      for obs in problem.petab_problem.measurement_df.index}

        code, set_of_params = problem._write_overrides('observable')
        print('code')
        print(code)
        self.assertTrue(code in JL_CODE_GOLD)
        self.assertEqual(set_of_params, set())

        code, set_of_params = problem._write_overrides('noise')
        print('code')
        print(code)
        self.assertTrue(code in JL_CODE_GOLD)
        self.assertEqual(set_of_params, set())

    def test_write_prior_code(self):
        problem = Mock()
        problem._petab_problem = petab.problem.Problem()
        problem._petab_problem = problem._petab_problem.from_yaml(PETAB_YAML)

        code = problem._write_prior_code()
        print('code')
        print(code)
        self.assertTrue(code in JL_CODE_GOLD)

    def test_resetting(self):
        problem = core.SBML2JuliaProblem(PETAB_YAML)
        problem.write_jl_file(path=os.path.join('.', 'test_resetting.jl'))
        self.assertEqual(problem.julia_code, JL_CODE_GOLD)
        problem.t_steps = 3
        self.assertNotEqual(problem.julia_code, JL_CODE_GOLD)

        problem = core.SBML2JuliaProblem(PETAB_YAML)
        self.assertEqual(problem.julia_code, JL_CODE_GOLD)
        problem.n_starts = 2
        self.assertNotEqual(problem.julia_code, JL_CODE_GOLD)

        problem = core.SBML2JuliaProblem(PETAB_YAML)
        self.assertEqual(problem.julia_code, JL_CODE_GOLD)
