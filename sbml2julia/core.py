"""Class to fit PEtab problem via time discretization in Julia
:Author: Paul Lang <paul.lang@wolfson.ox.ac.uk>
:Date: 2020-04-15
:Copyright: 2020, Paul F Lang
:License: MIT
"""

import importlib
import libsbml
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import petab
import re
import sys
import warnings
import yaml
from julia.api import Julia
importlib.reload(libsbml)


class SBML2JuliaProblem(object):

    """Class to create and solve an optimization and retreive the results
    """
    
    def __init__(self, petab_yaml, t_steps=None, n_starts=1, infer_ic_from_sbml=False,
                 optimizer_options={}, custom_code_dict={}):
        """
        Args:
            petab_yaml (:obj:`str`): path to petab yaml file
            t_steps (:obj:`int`, optional): number of time-discretiation steps
            n_starts (:obj:`int`): number of multistarts
            infer_ic_from_sbml (:obj:`bool`): infer missing initial conditions from SBML
            optimizer_options (:obj:`dict`): optimization solver options
            custom_code_dict (:obj:`dict`): dict with replaced code as keys
            and replacement code as values
        """
        print('Initialising problem...')
        self._calling_function = sys._getframe(1).f_code.co_name
        self._initialization = True
        self._optimized = False
        self._files_written = False
        self._plotted = False

        self._jl = Julia(compiled_modules=False)
        # self._results = {}
        self._petab_dirname = os.path.dirname(petab_yaml)
        self._set_petab_problem(petab_yaml)
        self.t_steps = t_steps
        self.n_starts = n_starts
        self.optimizer_options = optimizer_options
        self.infer_ic_from_sbml = infer_ic_from_sbml
        self._set_julia_code()
        self.custom_code_dict = custom_code_dict

        self._initialization = False
        print('Problem initialized.')

    @property
    def petab_yaml_dict(self):
        """Get petab_yaml_dict
        
        Returns:
            :obj:`dict`: petab_yaml_dict
        """
        return self._petab_yaml_dict

    @property
    def petab_problem(self):
        """Get petab_problem
        
        Returns:
            :obj:`petab.problem.Problem`: petab problem
        """
        return self._petab_problem

    @property
    def t_steps(self):
        """Get t_steps
        
        Returns:
            t_steps (:obj:`int`, optional): number of time-discretiation steps
        """
        return self._t_steps

    @t_steps.setter
    def t_steps(self, value):
        """Set t_steps
        
        Args:
            value (:obj:`int`, optional): number of time-discretiation steps
        
        Raises:
            ValueError: if t_steps is not a positive integer.
        """
        if value is None:
            n_exp = len(set(self.petab_problem.measurement_df['time']))
            if n_exp == 1:
                value = 101
            else:
                value = int(np.ceil(100/(n_exp-1))*(n_exp-1) + 1)
        if not (isinstance(value, int)) or (value <= 0):
            raise ValueError('`t_steps` must be a positive integer.')
        self._t_steps = value
        if not self._initialization:
            self._set_julia_code()

    @property
    def n_starts(self):
        """Get n_starts
        
        Returns:
            :obj:`int`: number of multistarts
        """
        return self._n_starts

    @n_starts.setter
    def n_starts(self, value):
        """Set n_starts
        
        Args:
            value (:obj:`int`): number of multistarts
        
        Raises:
            ValueError: if n_starts is not a positive integer
        """
        if not isinstance(value, int) or not (value > 0):
            raise ValueError('`n_starts` must be a positive integer')
        self._n_starts = value
        if not self._initialization:
            self._set_julia_code()

    @property
    def infer_ic_from_sbml(self):
        """Get infer_ic_from_sbml
        
        Returns:
            :obj:`bool`: if missing initial conditions shall be infered from SBML model
        """
        return self._infer_ic_from_sbml

    @infer_ic_from_sbml.setter
    def infer_ic_from_sbml(self, value):
        """Set infer_ic_from_sbml
        
        Args:
            value (:obj:`bool): if missing initial conditions shall be infered from SBML model
        
        Raises:
            ValueError: if infer_ic_from_sbml is not boolean
        """
        if not isinstance(value, bool):
            raise ValueError('`infer_ic_from_sbml` must be boolean')
        self._infer_ic_from_sbml = value
        if not self._initialization:
            self._set_julia_code()

    @property
    def optimizer_options(self):
        """Get optimizer_options
        
        Returns:
            :obj:`dict`: optimization solver options
        """
        return self._optimizer_options

    @optimizer_options.setter
    def optimizer_options(self, value):
        """Set optimizer_options
        
        Args:
            value (:obj:`dict`): optimization solver options
        
        Raises:
            ValueError: if optimizer_options is not a dict
        """
        if not isinstance(value, dict):
            raise ValueError('`optimizer_options` must be a dictionary')
        self._optimizer_options = value
        if not self._initialization:
            self._set_julia_code()

    @property
    def custom_code_dict(self):
        """Get custom_code_dict
        
        Returns:
            :obj:`dict`: custom code dict
        """
        return self._custom_code_dict

    @custom_code_dict.setter
    def custom_code_dict(self, value):
        """Set custom_code_dict
        
        Args:
            value (:obj:`dict`): custom code dict
        
        Raises:
            ValueError: if custom_code_dict is not a dict
            ValueError: if dict values are not strings
            ValueError: if dict keys are not string
        """
        if not isinstance(value, dict):
            raise ValueError('`custom_code_dict` must be a dictionary')
        for k in value.keys():
            if not isinstance(k, str):
                raise ValueError('Keys of `custom_code_dict` must be strings '
                                 f'but `{k}` is `{type(k)}`.')
        for v in value.values():
            if not isinstance(v, str):
                raise ValueError('Values of `custom_code_dict` must be strings '
                                 f'but `{v}` is `{type(v)}`.')

        self._custom_code_dict = value
        if not self._initialization:
            self._set_julia_code()
        if value:
            self.insert_custom_code(value)

    @property
    def julia_code(self):
        """Get julia_code
        
        Returns:
            :obj:`str`: julia code for optimization
        """
        return self._julia_code

    @property
    def results(self):
        """Get results
        
        Returns:
            :obj:`dict`: optimization results
        """
        return self._results

    def import_julia_code(self, file):
        """Summary
        
        Args:
            file (TYPE): Description
        """
        with open(file, 'r') as f:
            self._julia_code = f.read()

    def _set_petab_problem(self, petab_yaml):
        """Converts petab yaml to dict and creates petab.problem.Problem object
        
        Args:
            petab_yaml (:obj:`str`): path to petab yaml file
        """
        petab_problem = petab.problem.Problem()
        petab_problem = petab_problem.from_yaml(petab_yaml)

        petab.lint.lint_problem(petab_problem)  # Returns `False` if no error occured. Raises exception otherwise.
        self._check_for_not_implemented_features(petab_problem)
        petab_problem = self._sort_condition_df_problem(petab_problem)

        self._petab_problem = petab_problem

        self._petab_yaml_dict, self._condition2index, self._j_to_parameters,\
            self._n_conditions, self._condition_specific_pars, self._global_pars =\
            self._get_translation_vars(petab_yaml, petab_problem)

    def _check_for_not_implemented_features(self, petab_problem):
        """Checks if the petab_problem contains not implemented features
        
        Args:
            petab_problem (:obj:`petab.problem.Problem`): PEtab problem
        
        Raises:
            NotImplementedError: if `time` column in measurement table contains `Inf`
            NotImplementedError: if simulation conditions are associated with > 1 preequilibration
            NotImplementedError: if conditions are used for both, preequilibration and simulation
            NotImplementedError: if priors for noise parameters are provided
            NotImplementedError: if priors for observable parameters are provided

        """
        if np.inf in list(petab_problem.measurement_df['time']):
            raise NotImplementedError('Fitting steady state problems is not implemented.')  # Todo: consider implementing it.

        t_conds = []
        for obs, data_1 in petab_problem.measurement_df.groupby('observableId'):
            t_conds.append(tuple(data_1['time']))
        if 'preequilibrationConditionId' in petab_problem.measurement_df.columns \
                and not petab_problem.measurement_df['preequilibrationConditionId'].empty:
            p2c = petab_problem.measurement_df.loc[:, ['preequilibrationConditionId',
                                                   'simulationConditionId']].drop_duplicates()
            for sCId, data in p2c.groupby('simulationConditionId'):
                if len(data.index) > 1:
                    raise NotImplementedError(f'{sCId} must be assiciated with <= 1 '
                                              'preequilibrationConditionIds. Please modify PEtab'
                                              'problem accordingly.')
            condition_conflicts = set(p2c['preequilibrationConditionId'])\
                .intersection(set(p2c['simulationConditionId']))
            if len(condition_conflicts) > 0:
                raise NotImplementedError('The following conditions are used for both, simulation '
                                          'and preequilibration. Please modify '
                                          f'PEtab problem. {condition_conflicts}.')

        noiseParameter_names = set()
        if 'noiseParameters' in petab_problem.measurement_df.columns:
            for vals in petab_problem.measurement_df['noiseParameters']:
                pars = str(vals).rstrip(';').split(';')
                for par in pars:
                    noiseParameter_names.add(par.strip())
        observableParameter_names = set()
        if 'observableParameters' in petab_problem.measurement_df.columns:
            for vals in petab_problem.measurement_df['observableParameters']:
                pars = str(vals).rstrip(';').split(';')
                for par in pars:
                    observableParameter_names.add(par.strip())
        if 'objectivePriorType' in petab_problem.parameter_df.columns:
            for l, par in enumerate(petab_problem.parameter_df.index):
                if par in noiseParameter_names and\
                        isinstance(petab_problem.parameter_df['objectivePriorType'][l], str):
                    raise NotImplementedError('Priors for noiseParameter '
                                              'overrides are not implemented.')
                if par in observableParameter_names and\
                        isinstance(petab_problem.parameter_df['objectivePriorType'][l], str):
                    raise NotImplementedError('Priors for observableParameter overrides are '
                                              'not implemented.')

    def _sort_condition_df_problem(self, petab_problem):
        """Sorts the rows of the contition table based on the first
        occurence of the respective condition in the measurement table
        
        Args:
            petab_problem (:obj:`petab.problem.Problem`): PEtab problem
        
        Returns:
            :obj:`petab.problem.Problem`: PEtab problem
        """
        idx = 1e6*np.ones(len(petab_problem.condition_df.index))
        for i, cond in enumerate(petab_problem.measurement_df['simulationConditionId']
                                 .drop_duplicates()):
            for j, c in enumerate(petab_problem.condition_df.index):
                if c == cond:
                    idx[j] = i
        petab_problem.condition_df['sorting'] = idx
        petab_problem.condition_df = petab_problem.condition_df\
            .sort_values(by='sorting').drop(columns=['sorting'])

        return petab_problem

    def _get_translation_vars(self, petab_yaml, petab_problem):
        """Gets variables required for translation from PEtab to JuMP
        
        Args:
            petab_yaml (:obj:`string`): path to petab yaml file
            petab_problem (:obj:`petab.problem.Problem`): PEtab problem

        Raises:
            SystemExit: if `yaml.YAMLError` occured
        
        Returns:
            :obj:`tuple`: (yaml_dict, condition2index, j_to_parameters, n_conditions,
                           condition_specific_pars, global_pars)
        """
        with open(petab_yaml, 'r') as f:
            try:
                yaml_dict = yaml.safe_load(f)
            except yaml.YAMLError as error:
                raise SystemExit('Error occured: {}'.format(str(error)))

        condition2index = {petab_problem.condition_df.index[i]:
                           i for i in range(len(petab_problem.condition_df.index))}

        simulationConditionIdx = list(np.arange(len(condition2index))+1)
        preequilibrationConditionIdx = ['']*len(condition2index)
        if 'preequilibrationConditionId' in petab_problem.measurement_df.columns\
                and not self.petab_problem.measurement_df['preequilibrationConditionId'].empty:
            p2c = self.petab_problem.measurement_df\
                .loc[:, ['preequilibrationConditionId', 'simulationConditionId']].drop_duplicates()
            simulationConditionIdx = [condition2index[c]+1 for c in p2c['simulationConditionId']]
            preequilibrationConditionIdx = [condition2index[c]+1 if isinstance(c, str)
                                            else '' for c in p2c['preequilibrationConditionId']]
        j_to_parameters = (simulationConditionIdx, preequilibrationConditionIdx)

        n_conditions = len(petab_problem.condition_df.index)

        condition_specific_pars = {}
        for parameter in petab_problem.condition_df.columns:
            if parameter != 'conditionName':
                condition_specific_pars[parameter] = [val for val in
                                                      petab_problem.condition_df[parameter]]

        global_pars = {}
        for parameter in petab_problem.parameter_df.index:
            global_pars[parameter] = petab_problem.parameter_df.loc[parameter, 'estimate']

        return (yaml_dict, condition2index, j_to_parameters, n_conditions,
                condition_specific_pars, global_pars)

    def insert_custom_code(self, custom_code_dict):
        """Inserts custom code into Julia code
        
        Args:
            custom_code_dict (:obj:`dict`): dict with replaced code as keys
                                            and replacement code as values
        """
        positions = custom_code_dict.keys()
        code = self.julia_code
        for pos in positions:
            code = re.sub(pos, custom_code_dict[pos], code)
        self._julia_code = code

    def write_jl_file(self, path=os.path.join('.', 'julia_code.jl')):
        """Write code to julia file
        
        Args:
            path (:obj:`str`, optional): path to output Julia file
        """
        with open(path, 'w') as f:
            f.write(self.julia_code)
            self._julia_file = path
            self._files_written = True

    def optimize(self):
        """Optimize SBML2JuliaProblem
        
        Returns:
            :obj:`dict`: Results in a dict with keys 'species',
                         'observables', 'parameters' and 'par_est'
        """
        print('Entering Julia for optimization...')
        self._results_all = self._jl.eval(self.julia_code)
        print('Results transferred. Exited Julia.')

        self._best_iter = min(self._results_all['objective_value'],
                              key=self._results_all['objective_value'].get)

        self._results = {}
        self._results['par_est'] = self._get_param_ratios(self._results_all['parameters'])
        self._results['species'] = self._results_to_frame(
            self._results_all['species'][self._best_iter], variable_type='speciesId')
        self._results['observables'] = self._results_to_frame(
            self._results_all['observables'][self._best_iter], variable_type='observableId')

        self.petab_problem.simulation_df = self._results['observables']
        # Todo: remove the removal of the `observableParamters` column once the bug it causes in petab.calculate_llh is fixed.
        cols = [not b for b in self.petab_problem.measurement_df.columns
                .isin(['observableParameters'])]  # , 'noiseParameters'])]
        ndf = pd.DataFrame()
        if 'noiseParameters' in self.petab_problem.measurement_df.columns:
            ndf = self.petab_problem.measurement_df['noiseParameters']
        try:
            self._results['fval'] = -petab.calculate_llh(
                self.petab_problem.measurement_df.loc[:, cols],
                pd.concat([self.petab_problem.simulation_df.rename(
                    columns={'measurement': 'simulation'}), ndf], axis=1),
                self.petab_problem.observable_df,
                self.petab_problem.parameter_df)
            if not ('objectivePriorType' in self.petab_problem.parameter_df.columns
                    or 'objectivePriorParameters' in self.petab_problem.parameter_df.columns):
                if not np.isclose(self._results['fval'],
                                  self._results_all['objective_value'][self._best_iter]):
                    warnings.warn(
                        'Optimization algorithm may not have used correct objective'
                        f'(Julia llh: {self._results_all["objective_value"][self._best_iter]};'
                        f'PEtab llh: {self._results["fval"]}).')
            self._results['chi2'] = petab.calculate_chi2(
                self.petab_problem.measurement_df.loc[:, cols],
                pd.concat([self.petab_problem.simulation_df.rename(
                    columns={'measurement': 'simulation'}), ndf], axis=1),
                self.petab_problem.observable_df, self.petab_problem.parameter_df)
        except TypeError:
            warnings.warn('Could not calculate llh and/or chi2 using PEtab. '
                          'Using llh from Julia instead.')
            self._results['fval'] = self._results_all['objective_value'][self._best_iter]
            self._results['chi2'] = np.nan

            self._optimized = True
        return self.results

    def _get_param_ratios(self, par_dict):
        """Get ratios between optimized and nominal parameters
        
        Args:
            par_dict (:obj:`dict`): dict with parameter values
        
        Returns:
            :obj:`pandas.DataFrame`: Dataframe with optimized and nominal parameters
                                     and the ratio between them
        """
        local_par_names = {}
        for par_type in self.petab_problem.condition_df.columns:
            if str(self.petab_problem.condition_df[par_type].dtype) == 'object':
                for i in range(self._n_conditions):
                    local_par_names[
                        self.petab_problem.condition_df.iloc[i][par_type]
                        ] = (par_type, i)

        par_est = par_dict[str(self._best_iter)]
        par_est_col = []
        par_est_to_nom_col = []
        for parameterId in list(self.petab_problem.parameter_df.index):
            if self.petab_problem.parameter_df.loc[parameterId, 'estimate'] == 0:
                par_est_col.append('')
                par_est_to_nom_col.append('')
                continue
            nominal = self.petab_problem.parameter_df.loc[parameterId, 'nominalValue']
            if parameterId in par_est.keys():
                estimated = par_est[parameterId]
            else:
                par_type, i = local_par_names[parameterId]
                estimated = par_est[par_type][i]
            par_est_col.append(estimated)
            if nominal == 0:
                par_est_to_nom_col.append(f'NA (diff={estimated})')
            else:
                par_est_to_nom_col.append(estimated/nominal)

        df = self.petab_problem.parameter_df.copy()
        df['estimatedValue'] = par_est_col
        df['estimate_to_nominal'] = par_est_to_nom_col

        return df

    def _results_to_frame(self, simulation_dict, variable_type='observableId'):
        """Summary
        
        Args:
            simulation_dict (:obj:`dict`): dictionary of species or observable
                                           time-courses
            variable_type (:ob:`str`, optional): `speciesId` or `observableId`

        Raises:
            ValueError: if `variable_type` is not `speciesId` or `observableId`
        
        Returns:
            :obj:`pandas.DataFrame`: Dataframe of time-courses.
        """
        if variable_type not in ('speciesId', 'observableId'):
            raise ValueError('variable_type must be `speciesId` or `observableId`.')
        t_max = self.petab_problem.measurement_df['time'].max()
        time = np.linspace(start=0, stop=t_max, num=self.t_steps)
        index2condition = {v: k for k, v in self._condition2index.items()}
        res_dict = {variable_type: [], 'simulationConditionId': [], 'time': [], 'simulation': []}
        for variableId, dictionary in simulation_dict.items():
            for cond_idx, simulation in dictionary.items():
                simulationConditionId = index2condition[cond_idx-1]  # This is wrong
                if variable_type == 'observableId':
                    cond_bool = list(self.petab_problem.measurement_df['simulationConditionId']
                                     == simulationConditionId)
                    obs_bool = list(self.petab_problem.measurement_df[variable_type] == variableId)
                    time = self.petab_problem.measurement_df['time'].iloc[np.logical_and(
                        cond_bool, obs_bool)]
                res_dict[variable_type] = res_dict[variable_type] + [variableId]*len(time)
                res_dict['simulationConditionId'] = res_dict['simulationConditionId']\
                    + [simulationConditionId]*len(time)
                res_dict['time'] = res_dict['time'] + list(time)
                res_dict['simulation'] = res_dict['simulation'] + list(simulation)

        return pd.DataFrame(res_dict).sort_values(
            ['simulationConditionId', variable_type, 'time']).reset_index(drop=True)

    def write_results(self, path=os.path.join('.', 'results'), df_format='long'):
        """Write results to excel file
        
        Args:
            path (:obj:`str`, optional): path of excel file to write results to
            df_format (:obj:`str`, optional): long or wide table format
        
        Raises:
            ValueError: if `path` is not a directory or Excel file
        """
        if df_format not in ['long', 'wide']:
            warnings.warn('`df_format` must be `long` or `wide` but is {}. Defaulting to `long`.')
            df_format = 'long'

        pd.set_option("display.max_rows", None)
        sheets = {'parameters': self.results['par_est'],
                  'fval_chi2': pd.DataFrame({'value': [self._results['fval'],
                                                       self._results['chi2']]},
                                            index=['fval', 'chi2'])}
        if df_format == 'wide':
            for var_type, Id in [('species', 'speciesId'), ('observables', 'observableId')]:
                df = self.results[var_type].groupby('simulationConditionId')
                for condition, i in self._condition2index.items():
                    if i+1 in self._j_to_parameters[0]:
                        dfg = df.get_group(condition)
                        dfg = dfg.set_index(['time', Id])  # .drop_duplicates()
                        dfg = dfg.unstack()
                        dfg = dfg.loc[:, 'simulation']
                        sheets[var_type+'_'+condition] = dfg
        elif df_format == 'long':
            for var_type in ['species', 'observables']:
                df = self.results[var_type]
                sheets[var_type] = df

        name, ext = os.path.splitext(path)
        if ext == '.xlsx':
            with pd.ExcelWriter(path) as writer:
                for sheet_name, df in sheets.items():
                    write_idx = True
                    if (df_format == 'long') and (sheet_name not in ('parameters', 'fval_chi2')):
                        write_idx = False
                    df.to_excel(writer, sheet_name=sheet_name, index=write_idx)
        elif ext == '':
            if not os.path.isdir(path):
                os.mkdir(path)
            for sheet_name, df in sheets.items():
                write_idx = True
                if (df_format == 'long') and (sheet_name not in ('parameters', 'fval_chi2')):
                    write_idx = False
                df.to_csv(os.path.join(path, sheet_name+'.tsv'), sep='\t', index=write_idx)
        else:
            raise ValueError('`path` must be a directory or excel file.')

    def write_optimized_parameter_table(self):
        """Writes a new parameter table were nominal values are replaced with optimized values
        """
        df = self.petab_problem.parameter_df
        df['nominalValue'] = list(self.results['par_est']['estimatedValue'])
        out_file = os.path.join(self._petab_dirname, 'post_fit_parameters.tsv')
        df.to_csv(out_file, sep='\t')
        warnings.warn('Wrote post_fit_parameters.tsv. Please edit `yaml` file accordingly '
                      'if it shall be added to the petab problem.')

    def plot_results(self, condition, path=os.path.join('.', 'plot.pdf'),
                     observables=[], size=(6, 5)):
        """Plot results
        
        Args:
            condition (:obj:`str`): experimental condition to plot
            path (:obj:`str`, optional): path to output plot
            observables (:obj:`list`, optional): list of observables to be plotted
            size (:obj:`tuple`, optional): size of image
        
        Raises:
            ValueError: if `observables` is not a list
        """
        # Options
        x_label = 'time'
        y_label = 'Abundance'
        cols = [not b for b in self.petab_problem.measurement_df.columns.isin(
            ['observableParameters', 'noiseParameters', 'preequilibrationConditionId'])]
        measurement_df = self.petab_problem.measurement_df\
            .loc[:, cols].set_index(['simulationConditionId', 'time', 'observableId'])
        if sum(measurement_df.index.duplicated(keep='first')) > 0:
            warnings.warn('Some time points contain replicate measurements. '
                          'Only the first measurement is plotted.')

        measurement_df = measurement_df[~measurement_df.index.duplicated(keep='first')]
        measurement_df = measurement_df.unstack().loc[str(condition), :]  # measurement_df.drop_duplicates().unstack().loc[str(condition), :]
        measurement_df.columns = measurement_df.columns.droplevel()
        t = [measurement_df.index[i] for i in range(len(measurement_df.index))]
        t_sim = np.linspace(start=0, stop=t[-1], num=self.t_steps)
        if not isinstance(observables, list):
            raise ValueError('`observables` must be a list of observables.')
        if not observables:
            observables = self.petab_problem.observable_df.index
            observables = [obs for obs in observables if self._condition2index[condition]+1
                           in self._obs_to_conditions[obs]]

        df = self.results['observables'].groupby('simulationConditionId')
        df = df.get_group(condition)
        df = df.set_index(['time', 'observableId']).drop_duplicates()
        df = df.unstack()
        t_sim = [df.index[i] for i in range(len(df.index))]
        values = df.loc[:, 'simulation'].reset_index(drop=True, inplace=False)[observables]
        exp_data = measurement_df[observables]

        # Determine the size of the figure
        plt.figure(figsize=size)
        axes = plt.axes([0.1, 0.1, 0.8, 0.8])

        sim_lines = axes.plot(t_sim, np.array(values), linewidth=3)
        axes.set_prop_cycle(None)  # reset the color cycle
        exp_points = axes.plot(t, np.array(exp_data), 'x')
        sim_legend = axes.legend(sim_lines, observables, frameon=True,
                                 title='Simulation', loc='upper right')
        # axes.set_prop_cycle(None)
        axes.legend(exp_points, observables, frameon=True, title='Experiment', loc='upper left')
        axes.add_artist(sim_legend)

        plt.xlim(np.min(t), np.max(t_sim))
        plt.ylim(0, 1.05 * max(values.max().max(), exp_data.max().max()))
        plt.xlabel(x_label, fontsize=18)
        plt.ylabel(y_label, fontsize=18)
        plt.title('sbml2julia time course')

        plt.savefig(path)
        plt.close()

        self._plotted = True
        self._plot_file = path

    def _set_julia_code(self):
        """Transform petab.problem.Problem to Julia JuMP model.
        
        Raises:
            NotImplementedError: if SBML model contains assignment rules
            ValueError: if column `estimate` in parameter table contains values other than `0` or `1`.
        """
        # ---------------------------------------------------------------------- #
        """
        `_set_julia_code` is adapted from Frank T. Bergman
        Date: 2019
        Availability: https://groups.google.com/forum/#!topic/sbml-discuss/inS4Lzp3Ri8 or
        https://www.dropbox.com/s/2bfpiausejp0gd0/convert_reactions.py?dl=0
        and based on the methods published by Sungho Shin et al. in "Scalable Nonlinear
        Programming Framework for Parameter Estimation in Dynamic Biological System Models"
        Date: 2019
        Availability: https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1006828
        """
        # ---------------------------------------------------------------------- #

        if self._calling_function == '_execute_case':
            warnings.warn('Problem is called from PEtab test suite. Simulating with nominal parameter values.')

        # read the SBML from file
        sbml_filename = os.path.join(self._petab_dirname, self.petab_yaml_dict['problems'][0]['sbml_files'][0])
        doc = libsbml.readSBMLFromFile(sbml_filename)
        if doc.getNumErrors(libsbml.LIBSBML_SEV_FATAL):
            print('Encountered serious errors while reading file')
            print(doc.getErrorLog().toString())
            sys.exit(1)

        # clear errors
        doc.getErrorLog().clearLog()

        # perform conversions
        props = libsbml.ConversionProperties()
        props.addOption("promoteLocalParameters", True)

        if doc.convert(props) != libsbml.LIBSBML_OPERATION_SUCCESS:
            print('The document could not be converted')
            print(doc.getErrorLog().toString())

        props = libsbml.ConversionProperties()
        props.addOption("expandFunctionDefinitions", True)  # Todo: ask PEtab developers set this to `True` when creating `petab.problem.Problem()`

        if doc.convert(props) != libsbml.LIBSBML_OPERATION_SUCCESS:
            print('The document could not be converted')
            print(doc.getErrorLog().toString())

        mod = doc.getModel()

        assignments = {}
        for a in mod.getListOfRules():
            if (a.getMath().getName() is None) or (a.getMath().getName() == 'piecewise'):
                raise NotImplementedError('Assignment rules are not implemented.')
            assignments[a.getId()] = a.getMath().getName()

        initial_assignments = {}
        for a in mod.getListOfInitialAssignments():
            initial_assignments[a.getId()] = a.getMath().getName()

        reactions = {}  # dict of reaction and kinetic formula in JuMP format
        for i in range(mod.getNumReactions()):
            reaction = mod.getReaction(i)
            kinetics = reaction.getKineticLaw()
            kinetic_components = kinetics.getFormula()
            reactions[reaction.getId()] = kinetic_components  # jump_formula

        species = {}  # dict of species and stoichiometry-reactionId tuple they are involved in
        for i in range(mod.getNumSpecies()):
            specie = mod.getSpecies(i)
            if specie.getBoundaryCondition() is True or (specie.getId() in species):
                continue
            species[specie.getId()] = []

        for i in range(mod.getNumReactions()):
            reaction = mod.getReaction(i)
            kinetics = reaction.getKineticLaw()
            for j in range(reaction.getNumReactants()):
                ref = reaction.getReactant(j)
                specie = mod.getSpecies(ref.getSpecies())
                products = [r.getSpecies() for r in reaction.getListOfProducts()]
                if (specie.getBoundaryCondition() is True):
                    print('continueing...')
                    continue
                stoich = -ref.getStoichiometry()

                if specie.getName() in products:
                    stoich_p = reaction.getProduct(specie.getId()).getStoichiometry()
                    stoich = stoich_p + stoich
                if stoich < 0:
                    stoich = str(stoich)
                elif stoich > 0:
                    stoich = '+'+str(stoich)
                else:
                    stoich = ''
                if stoich:
                    species[specie.getId()].append((stoich, reaction.getId()))
                # print('added reaction {} to specie {}'.format(reaction.getID(), specie))
            for j in range(reaction.getNumProducts()):
                ref = reaction.getProduct(j)
                specie = mod.getSpecies(ref.getSpecies())
                reactants = [r.getSpecies() for r in reaction.getListOfReactants()]
                if (specie.getBoundaryCondition() is True) or (specie.getName() in reactants):
                    print('continueing')
                    continue
                species[specie.getId()].append((('+'+str(ref.getStoichiometry()), reaction.getId())))

        parameters = {}
        for i in range(mod.getNumParameters()):
            par = mod.getParameter(i)
            parameters[par.getId()] = par.getValue()

        try:
            species_file = os.path.join(self._petab_dirname, self.petab_yaml_dict['problems'][0]['species_files'][0])
            self.petab_problem.species_df = pd.read_csv(species_file, sep='\t', index_col='speciesId')
        except KeyError:
            ub = 11 * self.petab_problem.measurement_df.loc[:, 'measurement'].max()
            warnings.warn('Could not find `species_files` that specify lower and upper species boundaries in {}. Setting lower species boundaries to zero and upper species boundaries to {}.'.format(self._petab_dirname, ub))
            self.petab_problem.species_df = pd.DataFrame(
                {'speciesId': list(species.keys()), 'lowerBound': np.zeros(len(species)),
                 'upperBound': ub * np.ones(len(species))}).set_index('speciesId')


# -------start generating the code by appending to bytearray------- #
        generated_code = bytearray('', 'utf8')
        generated_code.extend(bytes('using CSV\n', 'utf8'))
        generated_code.extend(bytes('using DataFrames\n', 'utf8'))
        generated_code.extend(bytes('using Ipopt\n', 'utf8'))
        generated_code.extend(bytes('using JuMP\n\n', 'utf8'))

        generated_code.extend(bytes('n_steps = {} # Setting number of ODE discretisation steps\n\n'.format(self.t_steps), 'utf8'))

        generated_code.extend(bytes('# Data\n', 'utf8'))
        generated_code.extend(bytes('println("Reading measurement data...")\n', 'utf8'))
        generated_code.extend(bytes('data_path = "{}"\n'.format(os.path.join(self._petab_dirname, self.petab_yaml_dict['problems'][0]['measurement_files'][0])), 'utf8'))
        generated_code.extend(bytes('df = CSV.read(data_path)\n', 'utf8'))
        generated_code.extend(bytes('insert!(df, 1, (1:length(df[:,1])), :id)\n\n', 'utf8'))

        generated_code.extend(bytes('df_by_o = groupby(df, :observableId)\n', 'utf8'))
        generated_code.extend(bytes('df_by_c = groupby(df, :simulationConditionId)\n\n', 'utf8'))

        generated_code.extend(bytes('# Setting variables\n', 'utf8'))
        generated_code.extend(bytes('t_sim = range(0, stop=maximum(df[!, :time]), length=n_steps)\n', 'utf8'))
        generated_code.extend(bytes('t_exp = Dict()\n', 'utf8'))
        generated_code.extend(bytes('m_exp = Dict()\n', 'utf8'))
        generated_code.extend(bytes('t_sim_to_exp = Dict()\n', 'utf8'))
        generated_code.extend(bytes('deduplicated_time_idx = Dict()\n\n', 'utf8'))

        generated_code.extend(bytes('cond2idx = Dict()\n', 'utf8'))
        generated_code.extend(bytes('for i_cond in 1:length(keys(df_by_c))\n', 'utf8'))
        generated_code.extend(bytes('    cond_name = values(keys(df_by_c)[i_cond])[1]\n', 'utf8'))
        generated_code.extend(bytes('    cond2idx[cond_name] = i_cond\n', 'utf8'))
        generated_code.extend(bytes('end\n\n', 'utf8'))

        generated_code.extend(bytes('for i_obs in 1:length(keys(df_by_o))\n', 'utf8'))
        generated_code.extend(bytes('    obs_name = values(keys(df_by_o)[i_obs])[1]\n', 'utf8'))
        generated_code.extend(bytes('    t_exp[obs_name] = Dict()\n', 'utf8'))
        generated_code.extend(bytes('    m_exp[obs_name] = Dict()\n', 'utf8'))
        generated_code.extend(bytes('    t_sim_to_exp[obs_name] = Dict()\n', 'utf8'))
        generated_code.extend(bytes('    deduplicated_time_idx[obs_name] = Dict()\n', 'utf8'))
        generated_code.extend(bytes('    df_by_o_c = groupby(DataFrame(df_by_o[i_obs]), :simulationConditionId)\n', 'utf8'))
        generated_code.extend(bytes('    for i_cond in 1:length(keys(df_by_o_c))\n', 'utf8'))
        generated_code.extend(bytes('        cond_name = values(keys(df_by_o_c)[i_cond])[1]\n', 'utf8'))
        generated_code.extend(bytes('        cond_idx = cond2idx[cond_name]\n', 'utf8'))
        generated_code.extend(bytes('        t_exp[obs_name][cond_idx] = df_by_o_c[i_cond][!, :time]\n', 'utf8'))
        generated_code.extend(bytes('        m_exp[obs_name][cond_idx] = df_by_o_c[i_cond][!, :measurement]\n\n', 'utf8'))

        generated_code.extend(bytes('        t = df_by_o_c[i_cond][!, :time]\n', 'utf8'))
        generated_code.extend(bytes('        tmp = []\n', 'utf8'))
        generated_code.extend(bytes('        for i in 1:length(t)\n', 'utf8'))
        generated_code.extend(bytes('            idx = argmin(abs.(t[i] .- t_sim))\n', 'utf8'))
        generated_code.extend(bytes('            append!(tmp, idx)\n', 'utf8'))
        generated_code.extend(bytes('        end\n', 'utf8'))
        generated_code.extend(bytes('        t_sim_to_exp[obs_name][cond_idx] = tmp\n\n', 'utf8'))

        generated_code.extend(bytes('        deduplicated_time_idx[obs_name][cond_idx] = []\n', 'utf8'))
        generated_code.extend(bytes('        idx = 0\n', 'utf8'))
        generated_code.extend(bytes('        prev = "a"\n', 'utf8'))
        generated_code.extend(bytes('        for t in df_by_o_c[i_cond][!, :time]\n', 'utf8'))
        generated_code.extend(bytes('            if t != prev\n', 'utf8'))
        generated_code.extend(bytes('                idx = idx + 1\n', 'utf8'))
        generated_code.extend(bytes('            end\n', 'utf8'))
        generated_code.extend(bytes('            push!(deduplicated_time_idx[obs_name][cond_idx], idx)\n', 'utf8'))
        generated_code.extend(bytes('            prev = t\n', 'utf8'))
        generated_code.extend(bytes('        end\n', 'utf8'))
        generated_code.extend(bytes('    end\n', 'utf8'))
        generated_code.extend(bytes('end\n\n', 'utf8'))

        generated_code.extend(bytes('obs2conds = Dict()\n', 'utf8'))
        generated_code.extend(bytes('for obs_dict in t_exp\n', 'utf8'))
        generated_code.extend(bytes('    obs_name = obs_dict[1]\n', 'utf8'))
        generated_code.extend(bytes('    obs2conds[obs_name] = []\n', 'utf8'))
        generated_code.extend(bytes('    for cond_dict in obs_dict[2]\n', 'utf8'))
        generated_code.extend(bytes('        cond_idx = cond_dict[1]\n', 'utf8'))
        generated_code.extend(bytes('        push!(obs2conds[obs_name], cond_idx)\n', 'utf8'))
        generated_code.extend(bytes('    end\n', 'utf8'))
        generated_code.extend(bytes('    sort!(obs2conds[obs_name])\n', 'utf8'))
        generated_code.extend(bytes('end\n\n', 'utf8'))

        generated_code.extend(bytes('results = Dict()\n', 'utf8'))
        generated_code.extend(bytes('results["objective_value"] = Dict()\n', 'utf8'))
        generated_code.extend(bytes('results["parameters"] = Dict()\n', 'utf8'))
        generated_code.extend(bytes('results["species"] = Dict()\n', 'utf8'))
        generated_code.extend(bytes('results["observables"] = Dict()\n\n', 'utf8'))

        generated_code.extend(bytes(f'j_to_cond_par = {self._j_to_parameters[0]}\n', 'utf8'))

        preequ_bool = [False if isinstance(item, str) else True for item in self._j_to_parameters[1]]
        cond_without_preequ = [item for item, b in zip(self._j_to_parameters[0], preequ_bool) if not b]
        cond_with_preequ = [item for item, b in zip(self._j_to_parameters[0], preequ_bool) if b]
        preequ = [item for item, b in zip(self._j_to_parameters[1], preequ_bool) if b]

        generated_code.extend(bytes(f'cond_without_preequ = {cond_without_preequ}\n', 'utf8'))
        generated_code.extend(bytes(f'cond_with_preequ = {cond_with_preequ}\n', 'utf8'))
        generated_code.extend(bytes(f'preequ = {preequ}\n\n', 'utf8'))

        if self.n_starts > 1:
            generated_code.extend(bytes('for i_start in 1:{}\n'.format(self._n_starts), 'utf8'))
        else:
            generated_code.extend(bytes('i_start = 1\n\n', 'utf8'))
        generated_code.extend(bytes('    species_dict = Dict()\n', 'utf8'))
        generated_code.extend(bytes('    obs_dict = Dict()\n\n', 'utf8'))
        generated_code.extend(bytes('    m = Model(with_optimizer(Ipopt.Optimizer))\n', 'utf8'))
        for k, v in self.optimizer_options.items():
            if isinstance(v, str):
                generated_code.extend(bytes('    set_optimizer_attribute(m,"{}","{}")\n'.format(k, v), 'utf8'))
            else:
                generated_code.extend(bytes('    set_optimizer_attribute(m,"{}",{})\n'.format(k, v), 'utf8'))
        generated_code.extend(bytes('\n', 'utf8'))

        # Write global parameters
        generated_code.extend(bytes('    # Define global parameters\n', 'utf8'))
        generated_code.extend(bytes('    println("Defining global parameters...")\n', 'utf8'))
        for parameter, estimate in self._global_pars.items():
            lb = self.petab_problem.parameter_df.loc[parameter, 'lowerBound']
            ub = self.petab_problem.parameter_df.loc[parameter, 'upperBound']
            nominal = self.petab_problem.parameter_df.loc[parameter, 'nominalValue']
            if self._calling_function == '_execute_case':
                estimate = 0
            if estimate == 1:
                generated_code.extend(bytes('    @variable(m, {0} <= {1} <= {2}, start={0}+({2}-({0}))*rand(Float64))\n'.format(lb, parameter, ub), 'utf8'))
            elif estimate == 0:
                generated_code.extend(bytes('    @variable(m, {0} == {1}, start={1})\n'.format(parameter, nominal), 'utf8'))
            else:
                raise ValueError('Column `estimate` in parameter table must contain only `0` or `1`.')
        generated_code.extend(bytes('\n', 'utf8'))

        # Write condition-local parameters
        generated_code.extend(bytes('    # Define condition-specific parameters\n', 'utf8'))
        generated_code.extend(bytes('    println("Defining condition-specific parameters...")\n', 'utf8'))
        species_interpreted_as_ic = []
        for k, v in self._condition_specific_pars.items():
            if k in species.keys():
                species_interpreted_as_ic.append(k)
                k = k+'_0'
            generated_code.extend(bytes('    @variable(m, {0}[1:{1}])\n'.format(k, self._n_conditions), 'utf8'))
            for i, par in enumerate(v):
                if str(par).replace('.', '', 1).replace('e-', '', 1).replace('e', '', 1).isdigit():
                    generated_code.extend(bytes('    @constraint(m, {}[{}] == {})\n'.format(k, i+1, par), 'utf8'))
                elif not (isinstance(par, float) and np.isnan(par)):
                    lb = self.petab_problem.parameter_df.loc[par, 'lowerBound']
                    ub = self.petab_problem.parameter_df.loc[par, 'upperBound']
                    nominal = self.petab_problem.parameter_df.loc[par, 'nominalValue']
                    estimate = self.petab_problem.parameter_df.loc[par, 'estimate']
                    if self._calling_function == '_execute_case':  # The test cases always simulate from the nominal value
                        estimate = 0
                    if estimate == 1:
                        generated_code.extend(bytes('    @constraint(m, {}[{}] == {})\n'.format(k, i+1, par), 'utf8'))
                    elif estimate == 0:
                        generated_code.extend(bytes('    @constraint(m, {}[{}] == {})\n'.format(k, i+1, nominal), 'utf8'))
            generated_code.extend(bytes('\n', 'utf8'))
        generated_code.extend(bytes('\n', 'utf8'))

        # Write overrides:
        generated_code.extend(bytes('    # Define overrides\n', 'utf8'))
        self._obs_to_conditions = {}
        for observable in self.petab_problem.observable_df.index:
            obs_in_condition = [j+1 for c, j in self._condition2index.items() if j+1 in self._j_to_parameters[0] and
                                c in list(self.petab_problem.measurement_df
                                .loc[self.petab_problem.measurement_df['observableId']
                                     == observable, 'simulationConditionId'])]
            self._obs_to_conditions[observable] = obs_in_condition
        override_code, set_of_observable_params = self._write_overrides('observable')
        generated_code.extend(bytes(override_code, 'utf8'))
        override_code, set_of_noise_params = self._write_overrides('noise')
        generated_code.extend(bytes(override_code, 'utf8'))

        # Write out compartment values
        generated_code.extend(bytes('    # Model compartments\n', 'utf8'))
        generated_code.extend(bytes('    println("Defining compartments...")\n', 'utf8'))
        for i in range(mod.getNumCompartments()):
            element = mod.getCompartment(i)
            if element.getId() not in self.petab_problem.condition_df.columns:
                generated_code.extend(bytes('    @variable(m, {0} == {1}, start={1})\n'.format(element.getId(), element.getSize()), 'utf8'))
        generated_code.extend(bytes('\n', 'utf8'))

        # Write species
        n_j = len(self._j_to_parameters[0])
        generated_code.extend(bytes('    # Model species\n', 'utf8'))
        generated_code.extend(bytes('    println("Defining species...")\n', 'utf8'))
        for specie in species.keys():
            if species[specie]:
                lb = self.petab_problem.species_df.loc[specie, 'lowerBound']  # Todo: write somwehere a method that check that the set of sbml model species == self.petab_problem.species_df.index
                ub = self.petab_problem.species_df.loc[specie, 'upperBound']
                generated_code.extend(bytes('    species_dict["{}"] = @variable(m, {} <= {}[j in 1:{}, k in 1:(length(t_sim)+1)] <= {})\n'.format(specie, lb, specie, n_j, ub), 'utf8'))
            else:
                generated_code.extend(bytes('    species_dict["{}"] = @variable(m, {}[j in 1:{}, k in 1:(length(t_sim)+1)])\n'.format(specie, specie, n_j), 'utf8'))
        generated_code.extend(bytes('\n', 'utf8'))

        # Write initial assignments
        generated_code.extend(bytes('    # Model initial assignments\n', 'utf8'))
        write_ic = False
        if ('preequilibrationConditionId' not in self.petab_problem.measurement_df.columns) or (np.nan in list(self.petab_problem.measurement_df['preequilibrationConditionId'])):
            write_ic = True
        generated_code.extend(bytes('    println("Defining initial assignments...")\n', 'utf8'))
        for specie, par in initial_assignments.items():
            if specie not in species_interpreted_as_ic and write_ic:
                if par in self._global_pars:
                    generated_code.extend(bytes('    @constraint(m, [j in 1:{}], {}[cond_without_preequ[j],1] == {})\n'.format(len(cond_without_preequ), specie, initial_assignments[specie]), 'utf8'))
                elif par in self._condition_specific_pars.keys():
                    generated_code.extend(bytes('    @constraint(m, [j in 1:{}], {}[cond_without_preequ[j],1] == {}[cond_without_preequ[j]])\n'.format(len(cond_without_preequ), specie, initial_assignments[specie]), 'utf8'))
                elif self.infer_ic_from_sbml:
                    formula = str(par).split()  # par can sometimes be a float, which would cause an error when splitting
                    for i in range(len(formula)):
                        if (formula[i] in parameters.keys()) and\
                            (formula[i] not in list(self._condition_specific_pars)
                             + list(self._global_pars)):
                            generated_code.extend(bytes('    @variable(m, {0} == {1}, start={1})\n'.format(formula[i], parameters[formula[i]]), 'utf8'))
                    generated_code.extend(bytes('    @constraint(m, [j in 1:{}], {}[cond_without_preequ[j],1] == {})\n'.format(len(cond_without_preequ), specie, initial_assignments[specie]), 'utf8'))
        for specie in species_interpreted_as_ic:
            generated_code.extend(bytes('    @constraint(m, [j in 1:{}], {}[j,1] == {}[j_to_cond_par[j]])\n'.format(n_j, specie, specie+'_0'), 'utf8'))
        generated_code.extend(bytes('\n', 'utf8'))

        # Write ODEs
        generated_code.extend(bytes('    # Model ODEs\n', 'utf8'))
        generated_code.extend(bytes('    println("Defining ODEs...")\n', 'utf8'))
        patterns = [par for par in self.petab_problem.condition_df.columns if par not in species.keys()]
        for specie in species:  # For every species
            if species[specie]:
                generated_code.extend(bytes('    println("{}")\n'.format(specie), 'utf8'))
                generated_code.extend(bytes('    @NLconstraint(m, [j in 1:{}, k in 1:length(t_sim)-1],\n'.format(n_j), 'utf8'))
                generated_code.extend(bytes('        {}[j, k+1] == {}[j, k] + ('.format(specie, specie), 'utf8'))
                for (coef, reaction_name) in species[specie]:  # For every reaction
                    reaction_formula = ' {}*( {} )'.format(coef, reactions[reaction_name])
                    for pattern in patterns:  # Add iterator `j` to condition-defined and local parameters
                        reaction_formula = re.sub('[( ]'+pattern+'[, )]', lambda matchobj: matchobj.group(0)[:-1]+'[j_to_cond_par[j]]'+matchobj.group(0)[-1:], reaction_formula)  # The matchobject starts with a `(` or ` ` and ends with a `,` ` ` or `)`. I insert `[j]` just before the ending of the matchobject.
                    for spec in species.keys():
                        tmp_iterator = '[j]'
                        if species[spec]:
                            tmp_iterator = '[j, k+1]'
                        reaction_formula = re.sub('[^a-zA-Z0-9_]'+spec+'[^a-zA-Z0-9_]', lambda matchobj: matchobj.group(0)[:-1]+tmp_iterator+matchobj.group(0)[-1:], reaction_formula)
                    reaction_formula = re.sub('pow', '^', reaction_formula)
                    if 'piecewise(' in reaction_formula:
                        raise NotImplementedError('`libsbml` model contains `piecewise` conditition in rate law. This is not implemented.')
                    generated_code.extend(bytes(reaction_formula, 'utf8'))
                generated_code.extend(bytes('     ) * ( t_sim[k+1] - t_sim[k] ) )\n', 'utf8'))
            else:  # Todo: think about handling if initial_assignments[specie] == None
                generated_code.extend(bytes('    @constraint(m, [j in 1:{}, k in 1:length(t_sim)], {}[j, k] == {}[j])\n'.format(n_j, specie, initial_assignments[specie]), 'utf8'))
        generated_code.extend(bytes('\n', 'utf8'))

        # Write pre-equilibration constraints
        generated_code.extend(bytes('    # Pre-equilibration constraints\n', 'utf8'))
        condition_iterator = f'1:{len(cond_with_preequ)}'
        if cond_with_preequ or preequ:
            generated_code.extend(bytes('    println("Defining pre-equilibration constraints...")\n', 'utf8'))
            if self._calling_function == '_execute_case':
                generated_code.extend(bytes('    @constraint(m, [j in {}], A[j, length(t_sim)+1] + B[j, length(t_sim)+1] == 1)\n       '.format(condition_iterator), 'utf8'))

            for specie in species:  # For every species
                # Non-preequilibrated conditions
                if cond_without_preequ:
                    generated_code.extend(bytes('    @constraint(m, [j in cond_without_preequ], {}[j, length(t_sim)+1] == 0) # Dummy preequilibration for these conditions\n'.format(specie), 'utf8'))  # Dummy preequilibration

                # Preequilibrated conditions
                if species[specie]:  # Species has reaction
                    generated_code.extend(bytes('    println("{}")\n'.format(specie), 'utf8'))
                    generated_code.extend(bytes('    @NLconstraint(m, [j in {}],\n       '.format(condition_iterator), 'utf8'))
                    for (coef, reaction_name) in species[specie]:  # For every reaction
                        reaction_formula = ' {}*( {} )'.format(coef, reactions[reaction_name])
                        for pattern in patterns:  # Add iterator `[preequilibration_idx[j]` to condition-defined and local parameters
                            reaction_formula = re.sub('[( ]'+pattern+'[, )]', lambda matchobj: matchobj.group(0)[:-1]+'[preequ[j]]'+matchobj.group(0)[-1:], reaction_formula)  # The matchobject starts with a `(` or ` ` and ends with a `,` ` ` or `)`. I insert `[j]` just before the ending of the matchobject.
                        for spec in species.keys():
                            tmp_iterator = '[cond_with_preequ[j]]'
                            if species[spec]:
                                tmp_iterator = '[cond_with_preequ[j], length(t_sim)+1]'
                            reaction_formula = re.sub('[^a-zA-Z0-9_]'+spec+'[^a-zA-Z0-9_]', lambda matchobj: matchobj.group(0)[:-1]+tmp_iterator+matchobj.group(0)[-1:], reaction_formula)
                        reaction_formula = re.sub('pow', '^', reaction_formula)
                        if 'piecewise(' in reaction_formula:
                            raise NotImplementedError('`libsbml` model contains `piecewise` conditition in rate law. This is not implemented.')
                        generated_code.extend(bytes(reaction_formula, 'utf8'))
                    generated_code.extend(bytes('  == 0 )\n', 'utf8'))
                    if specie not in species_interpreted_as_ic:  # Link preequilibration to ic.
                        generated_code.extend(bytes('    @constraint(m, [j in {0}], {1}[cond_with_preequ[j], length(t_sim)+1] == {1}[cond_with_preequ[j], 1])\n'.format(condition_iterator, specie), 'utf8'))
                else:  # Species has no reaction
                    generated_code.extend(bytes('    @constraint(m, [j in {}], {}[j, length(t_sim)+1] == 0) # Dummy preequilibration\n'.format(condition_iterator, specie), 'utf8'))  # Dummy preequilibration
        else:  # Just a dummy to fill the last time value that will never be used.
            for specie in species:
                generated_code.extend(bytes('    @constraint(m, [j in cond_without_preequ], {}[j, length(t_sim)+1] == 0) # Dummy preequilibration\n'.format(specie), 'utf8'))  # Todo: I think if else block is superfluous

        generated_code.extend(bytes('\n', 'utf8'))

        # Write observables
        generated_code.extend(bytes('    # Define observables\n', 'utf8'))
        generated_code.extend(bytes('    println("Defining observables...")\n', 'utf8'))
        for observable in self.petab_problem.observable_df.index:
            if observable not in list(self.petab_problem.measurement_df['observableId']):
                warnings.warn(f'Observable {observable} is not in the PEtab measurement table and will be omitted during optimization.')
                continue
            min_exp_val = np.min(self.petab_problem.measurement_df.loc[self.petab_problem.measurement_df.loc[:, 'observableId'] == observable, 'measurement'])
            max_exp_val = np.max(self.petab_problem.measurement_df.loc[self.petab_problem.measurement_df.loc[:, 'observableId'] == observable, 'measurement'])
            diff = max_exp_val - min_exp_val
            lb = min_exp_val - 1*diff
            ub = max_exp_val + 1*diff
            if self._calling_function == '_execute_case':
                lb = min(0, lb)
                ub = max(16, ub)
            if np.isnan(lb):
                lb = 'NaN'
            if np.isnan(ub):
                ub = 'NaN'

            generated_code.extend(bytes(f'    obs_dict["{observable}"] = @variable(m, {lb} <= {observable}[j in obs2conds["{observable}"], k in 1:length(t_exp["{observable}"][j])] <= {ub}, start=1.)\n', 'utf8'))
            formula = self.petab_problem.observable_df.loc[observable, 'observableFormula'].split()
            formula = ' '+''.join(formula)+' '
            for spec in species.keys():
                formula = re.sub('[^a-zA-Z0-9_]'+spec+'[^a-zA-Z0-9_]',
                                 lambda matchobj: matchobj.group(0)[:-1]+f'[j, t_sim_to_exp["{observable}"][j][k]]'+matchobj.group(0)[-1:],
                                 formula)
            for pat in patterns:
                formula = re.sub('[^a-zA-Z0-9_]'+pat+'[^a-zA-Z0-9_]',
                                 lambda matchobj: matchobj.group(0)[:-1]+'[j]'+matchobj.group(0)[-1:],
                                 formula)
            for obs_par_name in set_of_observable_params:
                formula = re.sub('[^a-zA-Z0-9_]'+obs_par_name[0]+'[^a-zA-Z0-9_]',
                                 lambda matchobj: matchobj.group(0)[:-1]+f'[j{obs_par_name[1]}]'+matchobj.group(0)[-1:],  # Todo: not sure if this is correct
                                 formula)
            formula = formula.strip()

            generated_code.extend(bytes(f'    @NLconstraint(m, [j in obs2conds["{observable}"], k in 1:length(t_exp["{observable}"][j])], {observable}[j, k] == {formula})\n', 'utf8'))
        generated_code.extend(bytes('\n', 'utf8'))

        # Write priors
        nlp = ''
        if 'objectivePriorType' in self.petab_problem.parameter_df.columns \
                or 'objectivePriorParameters' in self.petab_problem.parameter_df.columns:
            prior_code = self._write_prior_code()
            generated_code.extend(bytes(prior_code, 'utf8'))
            if prior_code:
                nlp = '+ neg_log_priors'

        # Write objective
        generated_code.extend(bytes('    # Define objective\n', 'utf8'))
        generated_code.extend(bytes('    println("Defining objective...")\n', 'utf8'))
        generated_code.extend(bytes('    @NLobjective(m, Min, ', 'utf8'))
        sums_of_nllhs = []
        for observable in self.petab_problem.observable_df.index:

            if observable not in list(self.petab_problem.measurement_df['observableId']):
                continue

            sigma = '( '+str(self.petab_problem.observable_df.loc[observable, 'noiseFormula'])+' )'
            for noise_par_name in set_of_noise_params:
                sigma = re.sub('[^a-zA-Z0-9_]'+noise_par_name[0]+'[^a-zA-Z0-9_]',
                    lambda matchobj: matchobj.group(0)[:-1]+f'[j{noise_par_name[1]}]'+matchobj.group(0)[-1:],
                    sigma)

            for spec in species.keys():
                sigma = re.sub('[^a-zA-Z0-9_]'+spec+'[^a-zA-Z0-9_]',
                    lambda matchobj: matchobj.group(0)[:-1]+f'[t_sim_to_exp["{observable}"][j][k]]'+matchobj.group(0)[-1:],
                    sigma)
            for obs in self.petab_problem.observable_df.index:
                sigma = re.sub('[^a-zA-Z0-9_"]'+obs+'[^a-zA-Z0-9_"]',
                    lambda matchobj: matchobj.group(0)[:-1]+'[j, k]'+matchobj.group(0)[-1:],
                    sigma)
            scale = 'lin'
            if 'observableTransformation' in self.petab_problem.observable_df.columns:
                scale = self.petab_problem.observable_df.loc[observable, 'observableTransformation']
            if isinstance(scale, float) and np.isnan(scale):
                scale = 'lin'
            if scale not in ['lin', 'log', 'log10']:
                raise ValueError(f'`scale` must be `lin`, `log`, or `log10` but is {scale}.')

            noise_distribution = 'normal'
            if 'noiseDistribution' in self.petab_problem.observable_df.columns:
                noise_distribution = self.petab_problem.observable_df.loc[observable, 'noiseDistribution']
            if isinstance(noise_distribution, float) and np.isnan(noise_distribution):
                noise_distribution = 'normal'
            if noise_distribution not in ['normal', 'laplace']:
                raise ValueError(f'`noiseDistribution` must be `normal` or `laplace` but is {noise_distribution}.')

            if noise_distribution == 'normal' and scale == 'lin':
                sums_of_nllhs.append('sum(0.5 * log(2*pi*({1})^2) + 0.5*(({0}[j, deduplicated_time_idx["{0}"][j][k]]-m_exp["{0}"][j][k])/({1}))^2 for j in obs2conds["{0}"] for k in 1:length(t_exp["{0}"][j]))\n'.format(observable, sigma))  # 1:length(dfg[j][:, :time])
            elif noise_distribution == 'normal' and scale == 'log':
                sums_of_nllhs.append('sum(0.5*log(2*pi*({1})^2*(m_exp["{0}"][j][k])^2) + 0.5*((log({0}[j, deduplicated_time_idx["{0}"][j][k]])-log(m_exp["{0}"][j][k]))/({1}))^2 for j in obs2conds["{0}"] for k in 1:length(t_exp["{0}"][j]))\n'.format(observable, sigma))
            elif noise_distribution == 'normal' and scale == 'log10':
                sums_of_nllhs.append('sum(0.5*log(2*pi*({1})^2*(m_exp["{0}"][j][k])^2*log(10)^2) + 0.5*((log10({0}[j, deduplicated_time_idx["{0}"][j][k]])-log10(m_exp["{0}"][j][k]))/({1}))^2 for j in obs2conds["{0}"] for k in 1:length(t_exp["{0}"][j]))\n'.format(observable, sigma))
            elif noise_distribution == 'laplace' and scale == 'lin':
                sums_of_nllhs.append('sum(log(2*{1}) + abs({0}[j, deduplicated_time_idx["{0}"][j][k]]-m_exp["{0}"][j][k])/({1})) for j in obs2conds["{0}"] for k in 1:length(t_exp["{0}"][j]))\n'.format(observable, sigma))
            elif noise_distribution == 'laplace' and scale == 'log':
                sums_of_nllhs.append('sum(log(2*{1}*m_exp["{0}"][j][k]) + abs(log({0}[j, deduplicated_time_idx["{0}"][j][k]])-log(m_exp["{0}"][j][k]))/({1})) for j in obs2conds["{0}"] for k in 1:length(t_exp["{0}"][j]))\n'.format(observable, sigma))
            elif noise_distribution == 'laplace' and scale == 'log10':
                sums_of_nllhs.append('sum(log(2*{1}*m_exp["{0}"][j][k]*log(10)) + abs(log10({0}[j, deduplicated_time_idx["{0}"][j][k]])-log10(m_exp["{0}"][j][k]))/({1})) for j in obs2conds["{0}"] for k in 1:length(t_exp["{0}"][j]))\n'.format(observable, sigma))

        generated_code.extend(bytes('        + '.join(sums_of_nllhs), 'utf8'))
        generated_code.extend(bytes(f'        {nlp})\n\n', 'utf8'))

        generated_code.extend(bytes('    println("Optimizing:")\n', 'utf8'))
        generated_code.extend(bytes('    optimize!(m)\n\n', 'utf8'))

        # Write code to get the solution
        julia_pars = list(self._global_pars.keys())
        generated_code.extend(bytes('    println("Transfering results to Python...")\n', 'utf8'))
        generated_code.extend(bytes('    parameter_names = ' + str(julia_pars).replace('\'', '') + '\n', 'utf8'))
        generated_code.extend(bytes('    parameter_values = Dict()\n', 'utf8'))
        generated_code.extend(bytes('    for p in parameter_names\n', 'utf8'))
        generated_code.extend(bytes('        if occursin("[", string(p))\n', 'utf8'))
        generated_code.extend(bytes('            parameter_values[split(string(p[1]), "[")[1]] = JuMP.value.(p)\n', 'utf8'))
        generated_code.extend(bytes('        else\n', 'utf8'))
        generated_code.extend(bytes('            parameter_values[string(p)] = JuMP.value.(p)\n', 'utf8'))
        generated_code.extend(bytes('        end\n', 'utf8'))
        generated_code.extend(bytes('    end\n\n', 'utf8'))

        generated_code.extend(bytes('    species_names = [', 'utf8'))
        for specie in species:
            generated_code.extend(bytes('"'+specie+'", ', 'utf8'))
        generated_code.extend(bytes(']\n', 'utf8'))
        generated_code.extend(bytes('    species_values = Dict(string(spec)=>Dict(cond_idx=>[value(species_dict[spec][cond_idx, val_idx]) for val_idx=1:length(t_sim)]\n', 'utf8'))
        generated_code.extend(bytes('        for cond_idx in 1:length(keys(df_by_c))) for spec in species_names)\n\n', 'utf8'))

        generated_code.extend(bytes('    observable_values = Dict(obs=>Dict(cond_idx=>[value(obs_dict[obs][cond_idx, val_idx]) for val_idx=1:length(list)]\n', 'utf8'))
        generated_code.extend(bytes('        for (cond_idx, list) in dict) for (obs, dict) in m_exp)\n\n', 'utf8'))

        generated_code.extend(bytes('    objective_val = objective_value(m)\n\n', 'utf8'))
        generated_code.extend(bytes('    results["objective_value"][string(i_start)] = objective_val\n', 'utf8'))
        generated_code.extend(bytes('    results["parameters"][string(i_start)] = parameter_values\n', 'utf8'))
        generated_code.extend(bytes('    results["species"][string(i_start)] = species_values\n', 'utf8'))
        generated_code.extend(bytes('    results["observables"][string(i_start)] = observable_values\n\n', 'utf8'))
        if self.n_starts > 1:
            generated_code.extend(bytes('    println("--------------------Finished $i_start. iteration.--------------------")\n', 'utf8'))
            generated_code.extend(bytes('end\n\n', 'utf8'))

        generated_code.extend(bytes('results', 'utf8'))


        # Updating self and files
        code = generated_code.decode()
        self._julia_code = code        
        # if self._optimized == True:
        #     self.optimize()
        if self._files_written == True:
            self.write_jl_file(self._julia_file)
        if self._plotted == True:
            self.plot_results(self._plot_file)


    def _write_overrides(self, var_type):
        """Write observableParameter or noiseParameter overrides
        
        Args:
            var_type (:obj:`str`): 'observable' or 'noise'
        
        Returns:
            :obj:`tuple`: (override_code.decode(), set_of_params)
        """
        set_of_params = set()
        override_code = bytearray('', 'utf8')
        override_code.extend(bytes('    # Define '+var_type+' overrides\n', 'utf8'))
        override_code.extend(bytes('    println("Defining '+var_type+'Parameter overrides...")\n', 'utf8'))
        if var_type+'Parameters' in self.petab_problem.measurement_df.columns and not self.petab_problem.measurement_df[var_type+'Parameters'].empty:
            params = {}
            for obs, data_1 in self.petab_problem.measurement_df.groupby('observableId'):
                params[obs] = {}
                same_par_for_all_t = True
                for cond, data_2 in data_1.groupby('simulationConditionId'):
                    if len(set(data_2[var_type+'Parameters'])) > 1:
                        same_par_for_all_t = False
                for cond, data_2 in data_1.groupby('simulationConditionId'):
                    if same_par_for_all_t:
                        params[obs][cond] = [data_2[var_type+'Parameters'].values[0]]
                    else:
                        params[obs][cond] = data_2[var_type+'Parameters'].values

            str_1 = ''
            str_2 = ''
            str_3 = ''
            for data_1 in iter(params.values()):
                for element in data_1.values():
                    if len(set(element)) > 1:
                        str_3 = True
            for obs, obs_in_condition in self._obs_to_conditions.items():
                if not self._obs_to_conditions[obs]:
                    continue
                data_1 = params[obs]
                n_par = len(str(next(iter(data_1.values()))[0]).rstrip(';').split(';'))

                if str_3:
                    str_1 = f', k in 1:length(t_exp["{obs}"][j])'
                    str_2 = ', k'
                for i in range(n_par):
                    override_code.extend(bytes(f'    @variable(m, {var_type}Parameter{i+1}_{obs}[j in obs2conds["{obs}"]{str_1}], start=1.)\n', 'utf8'))
                    set_of_params.add((f'{var_type}Parameter{i+1}_{obs}', str_2))

                j = 0
                for cond, arr in data_1.items():
                    j += 1
                    j = self._condition2index[cond]+1
                    for k, pars in enumerate(arr):
                        if str_3:
                            str_3 = f', {k+1}'
                        for element in iter(data_1.values()):
                            if len(element) > 1:
                                str_3 = f', {k+1}'
                        p = [par.strip() for par in str(pars).rstrip(';').split(';')]
                        for n, par in enumerate(p):
                            if par != 'nan':
                                override_code.extend(bytes(f'    @constraint(m, {var_type}Parameter{n+1}_{obs}[{j}{str_3}] == {par})\n', 'utf8'))
                override_code.extend(bytes('\n', 'utf8'))
        override_code.extend(bytes('\n', 'utf8'))

        return (override_code.decode(), set_of_params)


    def _write_prior_code(self):
        """Write code for objectivePriors
        
        Raises:
            Exception: if either `objectivePriorType` or `objectivePriorParameters` but not both are specicied
            NotImplementedError: if objectivePriorType in parameter table is not `normal`, `laplace`, `logNormal` or `logLaplace`
        
        Returns:
            :obj:`str`: prior code
        """
        prior_code = bytearray('', 'utf8')
        idx_with_prior_1 = list(~self.petab_problem.parameter_df['objectivePriorType'].isna())
        idx_with_prior_2 = list(~self.petab_problem.parameter_df['objectivePriorParameters'].isna())
        if sum(np.equal(idx_with_prior_1, idx_with_prior_2)) != len (idx_with_prior_1):
            raise Exception('All rows and only those rows containing an `objectivePriorType` must also contain `objectivePriorParameters`.')
        parameter_df = self.petab_problem.parameter_df.iloc[idx_with_prior_1, :]
        n_par_with_prior = sum(idx_with_prior_1)
        if sum(idx_with_prior_1) == 0:
            return prior_code.decode()

        prior_code.extend(bytes('    # Defining objectivePriors\n', 'utf8'))
        prior_code.extend(bytes('    println("Defining objectivePriors")\n', 'utf8'))
        prior_code.extend(bytes(f'    @variable(m, prior_mean[l in 1:{n_par_with_prior}], start=1.)\n', 'utf8'))
        for l, vals in enumerate(parameter_df['objectivePriorParameters']):
            val = vals.rstrip(';').split(';')[0].strip()
            prior_code.extend(bytes(f'    @constraint(m, prior_mean[{l+1}] == {val})\n', 'utf8'))
        prior_code.extend(bytes('\n', 'utf8'))

        prior_code.extend(bytes(f'    @variable(m, 0 <= prior_std[l in 1:{n_par_with_prior}], start=1.)\n', 'utf8'))
        for l, vals in enumerate(parameter_df['objectivePriorParameters']):
            val = vals.rstrip(';').split(';')[1].strip()
            prior_code.extend(bytes(f'    @constraint(m, prior_std[{l+1}] == {val})\n', 'utf8'))
        prior_code.extend(bytes('\n', 'utf8'))

        prior_code.extend(bytes(f'    @variable(m, par_est[l in 1:{n_par_with_prior}], start=1.)\n', 'utf8'))
        for l, par in enumerate(parameter_df.index):
            prior_code.extend(bytes(f'    @constraint(m, par_est[{l+1}] == {par})\n', 'utf8'))
        prior_code.extend(bytes('\n', 'utf8'))

        prior_types = {'normal': [], 'laplace': [], 'logNormal': [], 'logLaplace': []}
        for l, prior_type in enumerate(parameter_df['objectivePriorType']):
            if prior_type == 'normal':
                prior_types['normal'].append(l+1)
            elif prior_type == 'laplace':
                prior_types['laplace'].append(l+1)
            elif prior_type == 'logNormal':
                vals = parameter_df.loc[parameter_df.index[l], 'objectivePriorParameters']
                vals = vals.rstrip(';').split(';')
                mode = np.exp(float(vals[0])-float(vals[1])**2)
                if mode < 2e-8:
                    warnings.warn(f'Mode of {parameter_df.index[l]} is {mode}. A small parameter value may cause an error in JuMP when taking its log.')
                prior_types['logNormal'].append(l+1)
            elif prior_type == 'logLaplace':
                prior_types['logLaplace'].append(l+1)
            else:
                raise NotImplementedError('Only `objectivePriorType`s `normal`, `laplace`, `logNormal` and `logLaplace` are implemented.')

        for k, v in prior_types.items():
            if v:
                prior_code.extend(bytes(f'    {k}_priors = {v}\n', 'utf8'))
        if 'laplace' in list(parameter_df['objectivePriorType']):
            prior_code.extend(bytes(f'    @variable(m, delta_par[l in {prior_types["laplace"]}])\n', 'utf8'))
            prior_code.extend(bytes(f'    @constraint(m, [l in {prior_types["laplace"]}], delta_par[l] >= par_est[l] - prior_mean[l])\n', 'utf8'))
            prior_code.extend(bytes(f'    @constraint(m, [l in {prior_types["laplace"]}], delta_par[l] >= prior_mean[l] - par_est[l])\n', 'utf8'))
            prior_code.extend(bytes('\n', 'utf8'))

        if 'logLaplace' in list(parameter_df['objectivePriorType']):
            prior_code.extend(bytes(f'    @variable(m, delta_log_par[l in {prior_types["logLaplace"]}])\n', 'utf8'))
            prior_code.extend(bytes(f'    @NLconstraint(m, [l in {prior_types["logLaplace"]}], delta_log_par[l] >= log(par_est[l]) - prior_mean[l])\n', 'utf8'))
            prior_code.extend(bytes(f'    @NLconstraint(m, [l in {prior_types["logLaplace"]}], delta_log_par[l] >= prior_mean[l] - log(par_est[l]))\n', 'utf8'))
            prior_code.extend(bytes('\n', 'utf8'))

        for k, v in prior_types.items():
            if v:
                prior_code.extend(bytes(f'    @variable(m, nlp_{k})\n', 'utf8'))
                if k == 'normal':
                    prior_code.extend(bytes(f'    @NLconstraint(m, nlp_{k} == sum( 0.5 * log(2*pi*prior_std[l]^2) + 0.5*((par_est[l]-prior_mean[l])/prior_std[l])^2 for l in {k}_priors))\n', 'utf8'))
                elif k == 'laplace':
                    prior_code.extend(bytes(f'    @NLconstraint(m, nlp_{k} == sum( log(2*prior_std[l]) + delta_par[l]/prior_std[l] for l in {k}_priors))\n', 'utf8'))
                elif k == 'logNormal':
                    prior_code.extend(bytes(f'    @NLconstraint(m, nlp_{k} == sum( log(par_est[l]*prior_std[l]*sqrt(2*pi)) + 0.5*((log(par_est[l])-prior_mean[l])/prior_std[l])^2 for l in {k}_priors))\n', 'utf8'))
                elif k == 'logLaplace':
                    prior_code.extend(bytes(f'    @NLconstraint(m, nlp_{k} == sum( log(2*prior_std[l]*par_est[l]) + delta_log_par[l]/prior_std[l] for l in {k}_priors))\n', 'utf8'))
        prior_code.extend(bytes('\n', 'utf8'))

        prior_code.extend(bytes('    @variable(m, neg_log_priors)\n', 'utf8'))
        neg_log_priors = list(set(['nlp_'+v for v in parameter_df['objectivePriorType']]))
        prior_code.extend(bytes(f"    @constraint(m, neg_log_priors == {' + '.join(np.sort(neg_log_priors))})\n\n", 'utf8'))

        return prior_code.decode()
