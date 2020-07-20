""" Class to fit PEtab problem via time discretization in Julia
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
import scipy as sp
import sys
import tempfile
import warnings
import yaml
from julia.api import Julia
importlib.reload(libsbml)


class DisFitProblem(object):

    def __init__(self, petab_yaml, t_ratio=2, n_starts=1, infer_ic_from_sbml=False):
        """        
        Args:
            petab_yaml (:obj:`str`): path petab yaml file
            t_ratio (:obj:`int` or `float`, optional): number of time discretiation steps per time unit
            n_starts (:obj:`int`): number of multistarts
        """

        print('Initialising problem...')
        self._calling_function = sys._getframe(1).f_code.co_name
        self._initialization = True
        self._optimized = False
        self._files_written = False
        self._pickled = False
        self._plotted = False
        self._jl = Julia(compiled_modules=False)
        self._initialization = True
        self._results = {}
        self._petab_dirname = os.path.dirname(petab_yaml)
        self._set_petab_problem(petab_yaml)
        self.t_ratio = t_ratio
        self.n_starts = n_starts
        self.infer_ic_from_sbml = infer_ic_from_sbml
        self._set_julia_code()
        self._initialization = False

    @property
    def petab_yaml_dict(self):
        """Get petab_yaml_dict
        
        Returns:
            :obj:`dict`: petab_yaml_dict 
        """
        return self._petab_yaml_dict

    @property
    def t_ratio(self):
        """Get t_ratio
        
        Returns:
            :obj:`int` or `float`: number of time discretiation steps per time unit
        """
        return self._t_ratio

    @t_ratio.setter
    def t_ratio(self, value):
        """Set t_ratio
        
        Args:
            value (:obj:`int` or `float`): number of time discretiation steps per time unit
        
        Raises:
            ValueError: if t_ratio is not a positive real number.
        """
        if not (isinstance(value, int) or isinstance(value, float)) or (value <= 0):
            raise ValueError('`t_ratio` must be a positive real number.')
        self._t_ratio = value
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

    @property
    def petab_problem(self):
        """Get petab_problem
        
        Returns:
            :obj:`petab.problem.Problem`: petab problem
        """
        return self._petab_problem

    def write_jl_file(self, path=os.path.join('.', 'julia_code.jl')):
        """Write code to julia file
        
        Args:
            path (:obj:`str`, optional): path to output julia file
        """
        with open(path, 'w') as f:
            f.write(self.julia_code)
            self._julia_file = path
            self._files_written = True

    def optimize(self):
        """Optimize DisFitProblem
        
        Returns:
            :obj:`dict`: Results in a dict with keys 'species', 'observables', 'parameters' and 'par_best'
        """
        print('Entering Julia for optimization...')
        self._results_all = self._jl.eval(self.julia_code)
        print('Results transferred. Exited Julia.')

        self._best_iter = min(self._results_all['objective_value'], key=self._results_all['objective_value'].get)

        self._results = {}
        self._results['par_best'] = self._get_param_ratios(self._results_all['parameters'])
        df = self._results_to_frame(self._results_all['species'], variable_type='speciesId')
        self._results['species'] = df.sort_values(['speciesId', 'simulationConditionId', 'time'])
        df = self._results_to_frame(self._results_all['observables'], variable_type='observableId')
        self._results['observables'] = df.sort_values(['observableId', 'simulationConditionId', 'time'])
        self._set_simulation_df()
        # Todo: remove the removal of the `observableParamters` column once the bug it causes in petab.calculate_llh is fixed.
        cols = [not b for b in self.petab_problem.measurement_df.columns.isin(['observableParameters'])] #, 'noiseParameters'])]
        try:
            self._results['fval'] = -petab.calculate_llh(self.petab_problem.measurement_df.loc[:, cols],
                self.petab_problem.simulation_df.rename(columns={'measurement': 'simulation'}),
                self.petab_problem.observable_df, self.petab_problem.parameter_df) # self._results_all['objective_value'][self._best_iter]
            if not np.isclose(self._results['fval'], self._results_all['objective_value'][self._best_iter]):
                warnings.warn('Optimization algorithm may not have used correct objective (Julia llh: {}; PEtab llh: {}).'.format(self._results[self._results_all['objective_value'][self._best_iter], 'fval']))
            self._results['chi2'] = petab.calculate_chi2(self.petab_problem.measurement_df.loc[:, cols],
                self.petab_problem.simulation_df.rename(columns={'measurement': 'simulation'}),
                self.petab_problem.observable_df, self.petab_problem.parameter_df)
        except:
            warnings.warn('Could not calculate llh and/or chi2 using PEtab. Using llh from Julia instead.')
            self._results['fval'] = self._results_all['objective_value'][self._best_iter]
            self._results['chi2'] = np.nan

        self._optimized = True
        return self.results

    def plot_results(self, condition, path=os.path.join('.', 'plot.pdf'), observables=[], size=(6, 5)):
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
        cols = [not b for b in self.petab_problem.measurement_df.columns.isin(['observableParameters', 'noiseParameters'])]
        measurement_df = self.petab_problem.measurement_df\
            .loc[:, cols].set_index(['simulationConditionId', 'time', 'observableId'])\
            .unstack().loc[str(condition), :]
        measurement_df.columns = measurement_df.columns.droplevel()
        t = [measurement_df.index[i] for i in range(len(measurement_df.index))]
        t_sim = np.linspace(start=0, stop=t[-1], num=np.int(np.ceil(t[-1]*self.t_ratio+1)))
        if not isinstance(observables, list):
            raise ValueError('`observables` must be a list of observables.')
        if not observables:
            observables = self.petab_problem.observable_df.index
        
        df = self.results['observables'].groupby('simulationConditionId')
        df = df.get_group(condition)
        df = df.set_index(['time', 'observableId'])
        df = df.unstack()
        t_sim = [df.index[i] for i in range(len(df.index))]
        values = df.loc[:, 'simulation'].reset_index(drop=True, inplace=False)[observables]
        print('measureent_df')
        print(measurement_df)
        print(self.petab_problem.measurement_df)
        exp_data = measurement_df[observables]
        
        # Determine the size of the figure
        plt.figure(figsize=size)
        axes = plt.axes([0.1, 0.1, 0.8, 0.8])
        
        sim_lines = axes.plot(t_sim, values, linewidth=3)
        axes.set_prop_cycle(None) # reset the color cycle
        exp_points = axes.plot(t, exp_data, 'x')
        sim_legend = axes.legend(sim_lines, observables, frameon=True, title='Simulation', loc='upper right')
        # axes.set_prop_cycle(None)
        axes.legend(exp_points, observables, frameon=True, title='Experiment', loc='upper left')
        axes.add_artist(sim_legend)
        
        plt.xlim(np.min(t), np.max(t_sim))
        plt.ylim(0, 1.05 * max(values.max().max(), exp_data.max().max()))
        plt.xlabel(x_label, fontsize=18)
        plt.ylabel(y_label, fontsize=18)
        plt.title('DisFit time course')

        plt.savefig(path)
        plt.close()

        self._plotted = True
        self._plot_file = path

    def write_results(self, path=os.path.join('.', 'results.xlsx'), df_format='wide'):
        """Write results to excel file
        
        Args:
            path (:obj:`str`, optional): path of excel file to write results to.
        """
        if df_format not in ['long', 'wide']:
            warnings.warn('`df_format` must be `long` or `wide` but is {}. Defaulting to `wide`.')
            df_format = 'wide'

        with pd.ExcelWriter(path) as writer:
            self.results['par_best'].to_excel(writer, sheet_name='par_best', index=False)

            if df_format == 'wide':
                for var_type, Id in [('species', 'speciesId'), ('observables', 'observableId')]:
                    df = self.results[var_type].groupby('simulationConditionId')
                    for condition in self._condition2index.keys():
                        dfg = df.get_group(condition)
                        dfg = dfg.set_index(['time', Id])
                        dfg = dfg.unstack()
                        dfg = dfg.loc[:, 'simulation']
                        dfg.to_excel(writer, sheet_name=var_type+'_'+condition, index=True)
            else:
                for var_type in ['species', 'observables']:
                    df = self.results[var_type]
                    df.to_excel(writer, sheet_name=var_type, index=True)

    def _get_param_ratios(self, par_dict):
        
        par_best = par_dict[str(self._best_iter)]
        par_0 = dict(zip(list(self.petab_problem.parameter_df.index), self.petab_problem.parameter_df.loc[:, 'nominalValue']))

        local_par_names = {}
        for par_type in self.petab_problem.condition_df.columns:
            if str(self.petab_problem.condition_df[par_type].dtype) == 'object':
                for i in range(self._n_conditions):
                    local_par_names[self.petab_problem.condition_df.iloc[i][par_type]] = (par_type, i)

        par_best_to_par_0_col = []
        par_best_col = []
        for key in par_0.keys():
            if key in par_best.keys():
                if par_0[key] != 0:
                    par_best_to_par_0_col.append(par_best[key] / par_0[key])
                else:
                    par_best_to_par_0_col.append('NA (diff={})'.format(par_best[key]-par_0[key]))
                par_best_col.append(par_best[key])
            else:
                par_type, i = local_par_names[key]
                par_best_to_par_0_col.append(par_best[par_type][i] / par_0[key])
                par_best_col.append(par_best[par_type][i])

        name_col = [str(key) for key in par_0.keys()]
        par_0_col = [par_0[str(key)] for key in par_0.keys()]

        df = pd.DataFrame(list(zip(name_col, par_0_col, par_best_col, par_best_to_par_0_col)),
            columns=['Name', 'par_0', 'par_best', 'par_best_to_par_0'])
        df = df.sort_values(by=['Name']).reset_index(drop=True)

        return df

    def _set_simulation_df(self):

        simulation_df = self.results['observables']
        t_sim_to_gt_sim = []
        idx = -1
        for i in range(len(self.petab_problem.measurement_df.index)): # Todo: expand this for all dataframes in gt_simulation_dfs
            idx = np.argmin(abs(simulation_df.loc[(idx+1):, 'time']
                - self.petab_problem.measurement_df.loc[:, 'time'].iloc[i])) + idx+1
            t_sim_to_gt_sim.append(idx)
        simulation_df = simulation_df.iloc[t_sim_to_gt_sim, :].reset_index(drop=True)
        simulation_df[petab.TIME] = simulation_df[petab.TIME].astype(int)
        self.petab_problem.simulation_df = simulation_df.rename(columns={'simulation': 'measurement'})

    def _results_to_frame(self, simulation_dict, variable_type='observableId'):

        t_max = self.petab_problem.measurement_df['time'].max()
        t_sim = np.linspace(start=0, stop=t_max, num=np.int(np.ceil(t_max*self.t_ratio+1)))
        t_exp = sorted(set(self.petab_problem.measurement_df['time']))
        res_dict = {variable_type: [], 'simulationConditionId': [], 'time': [], 'simulation': []}
        for variable in simulation_dict[self._best_iter].keys():
            for c in self._condition2index.keys():
                value = simulation_dict[self._best_iter][variable][self._condition2index[c]]
                res_dict['simulationConditionId'] = res_dict['simulationConditionId'] + [c]*len(value)
                res_dict[variable_type] = res_dict[variable_type] + [variable]*len(value)
                if variable_type == 'speciesId':
                    res_dict['time'] = res_dict['time'] + list(t_sim)
                elif variable_type == 'observableId':
                    res_dict['time'] = res_dict['time'] + list(t_exp)
                res_dict['simulation'] = res_dict['simulation'] + list(value)

        print(t_exp)
        print(value)
                
        return pd.DataFrame(res_dict)

    def _set_petab_problem(self, petab_yaml):
        """Converts petab yaml to dict and creates petab.problem.Problem object
        
        Args:
            petab_yaml (:obj:`str`): path petab yaml file
        
        Raises:
            SystemExit: if petab yaml file cannot be loaded.
        """
        problem = petab.problem.Problem()                                                                                                                                                                          
        problem = problem.from_yaml(petab_yaml)
        petab.lint.lint_problem(problem) # Returns `False` if no error occured and raises exception otherwise.
        self._petab_problem = problem
        with open(petab_yaml, 'r') as f:
            try:
                self._petab_yaml_dict = yaml.safe_load(f)
            except yaml.YAMLError as error:
                raise SystemExit('Error occured: {}'.format(str(error)))
        if 'preequilibrationConditionId' in self._petab_problem.measurement_df.columns:
            raise NotImplementedError('Preequilibration is not implemented (DisFit does not simulate ODEs. Therefore it cannot determine the time until equilibration).')

        self._condition2index = {self.petab_problem.condition_df.index[i]: i for i in range(len(self.petab_problem.condition_df.index))}

        self._condition_defined_pars = {}
        self._local_pars = {}
        self._n_conditions = self.petab_problem.condition_df.shape[0]
        for parameter in self.petab_problem.condition_df.columns:
            if parameter != 'conditionName':
                if str(self.petab_problem.condition_df[parameter].dtype) in ('float64', 'int16', 'int64'):
                    self._condition_defined_pars[parameter] = [val for val in self.petab_problem.condition_df[parameter]]
                elif str(self.petab_problem.condition_df[parameter].dtype) == 'object':
                    self._local_pars[parameter] = [par for par in self.petab_problem.condition_df[parameter]]

        self._global_pars = {}
        for parameter in self.petab_problem.parameter_df.index:
            tmp = []
            for v in self._local_pars.values():
                tmp = tmp + v
            if parameter not in tmp:
                self._global_pars[parameter] = self.petab_problem.parameter_df.loc[parameter, 'estimate']

    def _set_julia_code(self):
        """Transform petab.problem.Problem to Julia JuMP model.
        """
        #----------------------------------------------------------------------#
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
        #----------------------------------------------------------------------#

        NEG_LLH_DICT = {'normal': '0.5 * (log(2*pi*{2}^2) + {2}^(-2) * ({0}{3} - {1})^2)',
            'laplace': 'log(2*{2}) + abs({0}{3} - {1})/{2}',
            'logNormal': 'log({0}{3}*{2}*sqrt(2*pi)) + (log({0}{3}) - {1})^2/(2*{2}^2)',
            # 'logLaplace': '{1} + log(sqrt(2)*{2}) - min((sqrt(2)/{2} - {1})*(log({0}{3}) - {1}) , (sqrt(2)/{2} + {1})*(-log({0}{3}) + {1}))'
            } # Todo: ask Sungho to get logLaplace working

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
        props.addOption("expandFunctionDefinitions", True) # Todo: ask PEtab developers set this to `True` when creating `petab.problem.Problem()`

        if doc.convert(props) != libsbml.LIBSBML_OPERATION_SUCCESS: 
            print('The document could not be converted')
            print(doc.getErrorLog().toString())

        mod = doc.getModel()

        initial_assignments = {}
        for a in mod.getListOfInitialAssignments():
            initial_assignments[a.getId()] = a.getMath().getName()

        reactions = {} # dict of reaction and kinetic formula in JuMP format
        for i in range(mod.getNumReactions()):
            reaction = mod.getReaction(i)
            kinetics = reaction.getKineticLaw()
            kinetic_components = kinetics.getFormula() #.split(' * ')[1:]
            # kinetic_components = re.sub('compartment \* ', '', kinetic_components)
            reactions[reaction.getId()] = kinetic_components #jump_formula

        species = {} # dict of species and stoichiometry-reactionId tuple they are involved in
        for i in range(mod.getNumSpecies()): 
            specie = mod.getSpecies(i)
            if specie.getBoundaryCondition() == True or (specie.getId() in species):
                continue
            species[specie.getId()] = []

        # if self.infer_ic_from_sbml: # This is to import ini
        #     props = libsbml.ConversionProperties()
        #     props.addOption("expandInitialAssignments", True)

        #     if doc.convert(props) != libsbml.LIBSBML_OPERATION_SUCCESS: 
        #         print('The document could not be converted')
        #         print(doc.getErrorLog().toString())

        #     mod_numeric_ic = doc.getModel()
        #     for key in species.keys():
        #         if key not in self.petab_problem.condition_df.columns:
        #             if key not in initial_assignments.keys() or initial_assignments[key] == None:
        #                 if key in initial_assignments.keys() and initial_assignments[key] == None:
        #                     warnings.warn(f'Could not parse initial assignment for {key}. Overwriting with its numerical value. This may overwrite parameters that should have been fitted.')
        #                 # get initialValue
        #                 element = mod_numeric_ic.getElementBySId(key)
        #                 if element.getTypeCode() == libsbml.SBML_SPECIES: 
        #                     if element.isSetInitialConcentration(): 
        #                         initial_assignments[key] = element.getInitialConcentration()
        #                     else: 
        #                         initial_assignments[key] = element.getInitialAmount()
        #                 else: 
        #                     initial_assignments[key] = element.getSize()


        for i in range(mod.getNumReactions()): 
            reaction = mod.getReaction(i)
            kinetics = reaction.getKineticLaw()   
            for j in range(reaction.getNumReactants()): 
                ref = reaction.getReactant(j)
                specie = mod.getSpecies(ref.getSpecies())
                products = [r.getSpecies() for r in reaction.getListOfProducts()]
                if (specie.getBoundaryCondition() == True) or (specie.getName() in products):
                    # print('continueing...')
                    continue
                species[specie.getId()].append(('-'+str(ref.getStoichiometry()), reaction.getId()))
                # print('added reaction {} to specie {}'.format(reaction.getID(), specie))
            for j in range(reaction.getNumProducts()): 
                ref = reaction.getProduct(j)
                specie = mod.getSpecies(ref.getSpecies())
                reactants = [r.getSpecies() for r in reaction.getListOfReactants()]
                if (specie.getBoundaryCondition() == True) or (specie.getName() in reactants): 
                    continue
                species[specie.getId()].append((('+'+str(ref.getStoichiometry()), reaction.getId())))

        parameters = {}
        for i in range(mod.getNumParameters()):
            par = mod.getParameter(i)
            parameters[par.getId()] = par.getValue()

        try:
            species_file = os.path.join(self._petab_dirname, self.petab_yaml_dict['problems'][0]['species_files'][0])
            self.petab_problem.species_df = pd.read_csv(species_file, sep='\t', index_col='speciesId')
        except:
            ub = 11 * self.petab_problem.measurement_df.loc[:, 'measurement'].max()
            warnings.warn('Could not find `species_files` that specify lower and upper species boundaries in {}. Setting lower species boundaries to zero and upper species boundaries to {}.'.format(self._petab_dirname, ub))
            self.petab_problem.species_df = pd.DataFrame({'speciesId': list(species.keys()), 'lowerBound': np.zeros(len(species)),
                'upperBound': ub * np.ones(len(species))}).set_index('speciesId')


#-------start generating the code by appending to bytearray-------#
        generated_code = bytearray('', 'utf8')
        generated_code.extend(bytes('using CSV\n', 'utf8'))
        generated_code.extend(bytes('using DataFrames\n', 'utf8'))
        generated_code.extend(bytes('using Ipopt\n', 'utf8'))
        generated_code.extend(bytes('using JuMP\n\n', 'utf8'))

        generated_code.extend(bytes('t_ratio = {} # Setting number of ODE discretisation steps\n\n'.format(self.t_ratio), 'utf8'))
 
        generated_code.extend(bytes('# Data\n', 'utf8'))
        generated_code.extend(bytes('data_path = "{}"\n'.format(os.path.join(self._petab_dirname, self.petab_yaml_dict['problems'][0]['measurement_files'][0])), 'utf8'))
        generated_code.extend(bytes('df = CSV.read(data_path)\n', 'utf8'))
        generated_code.extend(bytes('dfg = groupby(df, :simulationConditionId)\n', 'utf8'))
        generated_code.extend(bytes('data = []\n', 'utf8'))
        generated_code.extend(bytes('for condition in keys(dfg)\n', 'utf8'))
        generated_code.extend(bytes('    push!(data,unstack(dfg[condition], :time, :observableId, :measurement))\n', 'utf8'))
        generated_code.extend(bytes('end\n\n', 'utf8'))

        generated_code.extend(bytes('t_exp = Vector(DataFrame(groupby(dfg[1], :observableId)[1])[!, :time])\n', 'utf8'))
        generated_code.extend(bytes('t_sim = range(0, stop=t_exp[end], length=Int64(ceil(t_exp[end]*t_ratio+1)))\n', 'utf8'))
        generated_code.extend(bytes('t_sim_to_exp = []\n', 'utf8'))
        generated_code.extend(bytes('for i in 1:length(t_exp)\n', 'utf8'))
        generated_code.extend(bytes('    idx = argmin(abs.(t_exp[i] .- t_sim))\n', 'utf8'))
        generated_code.extend(bytes('    append!(t_sim_to_exp, idx)\n', 'utf8'))
        generated_code.extend(bytes('end\n\n', 'utf8'))

        generated_code.extend(bytes('results = Dict()\n', 'utf8'))
        generated_code.extend(bytes('results["objective_value"] = Dict()\n', 'utf8'))
        generated_code.extend(bytes('results["parameters"] = Dict()\n', 'utf8'))
        generated_code.extend(bytes('results["species"] = Dict()\n', 'utf8'))
        generated_code.extend(bytes('results["observables"] = Dict()\n', 'utf8'))
        generated_code.extend(bytes('for i_start in 1:{}\n'.format(self._n_starts), 'utf8'))  
        generated_code.extend(bytes('    m = Model(with_optimizer(Ipopt.Optimizer, tol=1e-6))\n\n', 'utf8'))
        
        # Write condition-defined parameters
        generated_code.extend(bytes('    # Define condition-defined parameters\n', 'utf8'))
        species_interpreted_as_ic = []
        for k, v in self._condition_defined_pars.items():
            if k in species.keys():
                species_interpreted_as_ic.append(k)
                k = k+'_0'
            generated_code.extend(bytes('    @variable(m, {0}[1:{1}])\n'.format(k, self._n_conditions), 'utf8'))
            for i, val in enumerate(v):
                generated_code.extend(bytes('    @constraint(m, {0}[{1}] == {2})\n'.format(k, i+1, val), 'utf8'))
            generated_code.extend(bytes('\n', 'utf8'))

        # Write condition-local parameters
        generated_code.extend(bytes('    # Define condition-local parameters\n', 'utf8'))
        for k, v in self._local_pars.items():
            if k in species.keys():
                species_interpreted_as_ic.append(k)
                k = k+'_0'
            generated_code.extend(bytes('    @variable(m, {0}[1:{1}])\n'.format(k, self._n_conditions), 'utf8'))
            for i, par in enumerate(v):
                lb = self.petab_problem.parameter_df.loc[par, 'lowerBound']
                ub = self.petab_problem.parameter_df.loc[par, 'upperBound']
                nominal = self.petab_problem.parameter_df.loc[par, 'nominalValue']
                estimate = self.petab_problem.parameter_df.loc[par, 'estimate']
                if self._calling_function == '_execute_case': # The test cases always simulate from the nominal value
                    estimate = 0
                if estimate == 1:
                    generated_code.extend(bytes('    @constraint(m, {} <= {}[{}] <= {})\n'.format(lb, k, i+1, ub), 'utf8'))
                elif estimate == 0:
                    generated_code.extend(bytes('    @constraint(m, {}[{}] == {})\n'.format(k, i+1, nominal), 'utf8'))
                else:
                    raise ValueError('Column `estimate` in parameter table must contain only `0` or `1`.')
            generated_code.extend(bytes('\n', 'utf8'))

        # Write global parameters
        generated_code.extend(bytes('    # Define global parameters\n', 'utf8'))
        for parameter, estimate in self._global_pars.items():
            lb = self.petab_problem.parameter_df.loc[parameter, 'lowerBound']
            ub = self.petab_problem.parameter_df.loc[parameter, 'upperBound']
            nominal = self.petab_problem.parameter_df.loc[parameter, 'nominalValue']
            if self._calling_function == '_execute_case':
                estimate = 0
            if estimate == 1:
                    generated_code.extend(bytes('    @variable(m, {0} <= {1} <= {2}, start={0}+({2}-{0})*rand(Float64))\n'.format(lb, parameter, ub), 'utf8'))
            elif estimate == 0:
                generated_code.extend(bytes('    @variable(m, {} == {})\n'.format(parameter, nominal), 'utf8'))
            else:
                raise ValueError('Column `estimate` in parameter table must contain only `0` or `1`.')
        generated_code.extend(bytes('\n', 'utf8'))

        # Write observables overrides:
        set_of_observable_params = set()
        if 'observableParameters' in self.petab_problem.measurement_df.columns:
            observable_params = {}
            for obs, data_1 in self.petab_problem.measurement_df.groupby('observableId'):
                observable_params[obs] = {}
                for cond, data_2 in data_1.groupby('simulationConditionId'):
                    if len(set(data_2['observableParameters'])) <= 1:
                        print('one')
                        observable_params[obs][cond] = [data_2['observableParameters'].values[0]]
                    else:
                        print('two')
                        observable_params[obs][cond] = data_2['observableParameters'].values
            print(observable_params)

            generated_code.extend(bytes('    # Define observable overrides\n', 'utf8'))
            for obs, data_1 in observable_params.items():
                n_par = len(str(next(iter(data_1.values()))[0]).split(';'))
                str_1 = ', k in 1:length(t_exp)'
                str_2 = ', k'
                if len(next(iter(data_1.values()))) <= 1:
                    str_1 = ''
                    str_2 = ''
                for i in range(n_par):
                    generated_code.extend(bytes('    @variable(m, observableParameter{}_{}[j in 1:{}{}])\n'.format(i+1, obs, self._n_conditions, str_1), 'utf8'))
                    set_of_observable_params.add((f'observableParameter{i+1}_{obs}', str_2))
                print(observable_params[obs][cond])
                for cond, values in data_1.items():
                    for k, params in enumerate(values):
                        str_3 = f', {k+1}'
                        if len(next(iter(data_1.values()))) <= 1:
                            str_3 = ''
                        pars = [par.strip() for par in str(params).split(';')]
                        for n, par in enumerate(pars):
                            generated_code.extend(bytes(f'    @constraint(m, observableParameter{n+1}_{obs}[{self._condition2index[cond]+1}{str_3}] == {par})\n', 'utf8'))
            generated_code.extend(bytes('\n', 'utf8'))
        
        # Write noise overrides:
        set_of_noise_params = set()
        if 'noiseParameters' in self.petab_problem.measurement_df.columns:
            noise_params = {}
            for obs, data_1 in self.petab_problem.measurement_df.groupby('observableId'):
                noise_params[obs] = {}
                for cond, data_2 in data_1.groupby('simulationConditionId'):
                    if len(set(data_2['noiseParameters'])) <= 1:
                        noise_params[obs][cond] = [data_2['noiseParameters'].values[0]]
                    else:
                        noise_params[obs][cond] = data_2['noiseParameters'].values
            print(noise_params)

            generated_code.extend(bytes('    # Define noise overrides\n', 'utf8'))
            for obs, data_1 in noise_params.items():
                n_par = len(str(next(iter(data_1.values()))[0]).split(';'))
                str_1 = ', k in 1:length(t_exp)'
                str_2 = ', k'
                if len(data_1.values()) <= 1:
                    str_1 = ''
                    str_2 = ''
                for i in range(n_par):
                    generated_code.extend(bytes('    @variable(m, noiseParameter{}_{}[j in 1:{}{}])\n'.format(i+1, obs, self._n_conditions, str_1), 'utf8'))
                    set_of_noise_params.add((f'noiseParameter{i+1}_{obs}', str_2))
                print(noise_params[obs][cond])
                for cond, values in data_1.items():
                    for k, params in enumerate(values):
                        str_3 = f', {k+1}'
                        if len(next(iter(data_1.values()))) <= 1:
                            str_3 = ''
                        pars = [par.strip() for par in str(params).split(';')]
                        for n, par in enumerate(pars):
                            generated_code.extend(bytes(f'    @constraint(m, noiseParameter{n+1}_{obs}[{self._condition2index[cond]+1}{str_3}] == {par})\n', 'utf8'))
            generated_code.extend(bytes('\n', 'utf8'))


        # Write out compartment values 
        generated_code.extend(bytes('    # Model compartments\n', 'utf8'))
        generated_code.extend(bytes('    println("Defining compartments...")\n', 'utf8'))
        for i in range(mod.getNumCompartments()):
            element = mod.getCompartment(i)
            if element.getId() not in self.petab_problem.condition_df.columns:
                generated_code.extend(bytes('    @variable(m, {} == {})\n'.format(element.getId(), element.getSize()), 'utf8'))
        generated_code.extend(bytes('\n', 'utf8'))

        # Write species
        generated_code.extend(bytes('    # Model species\n', 'utf8'))
        generated_code.extend(bytes('    println("Defining species...")\n', 'utf8'))
        for specie in species.keys():
            if species[specie]:
                lb = self.petab_problem.species_df.loc[specie, 'lowerBound'] #Todo: write somhere a linter that check that the set of sbml model species == self.petab_problem.species_df.index
                ub = self.petab_problem.species_df.loc[specie, 'upperBound']
                generated_code.extend(bytes('    @variable(m, {} <= {}[j in 1:{}, k in 1:length(t_sim)] <= {})\n'.format(lb, specie, self._n_conditions, ub), 'utf8'))
            else:
                generated_code.extend(bytes('    @variable(m, {}[j in 1:{}, k in 1:length(t_sim)])\n'.format(specie, self._n_conditions), 'utf8'))
        generated_code.extend(bytes('\n', 'utf8'))

        # Write initial assignments
        generated_code.extend(bytes('    # Model initial assignments\n', 'utf8'))
        generated_code.extend(bytes('    println("Defining initial assignments...")\n', 'utf8'))
        for specie, par in initial_assignments.items():
            if specie not in species_interpreted_as_ic:
                if par in self._global_pars:
                    generated_code.extend(bytes('    @constraint(m, [j in 1:{}], {}[j,1] == {})\n'.format(self._n_conditions, specie, initial_assignments[specie]), 'utf8'))
                elif par in list(self._local_pars.keys())+list(self._condition_defined_pars.keys()):
                    generated_code.extend(bytes('    @constraint(m, [j in 1:{}], {}[j,1] == {}[j])\n'.format(self._n_conditions, specie, initial_assignments[specie]), 'utf8'))
                elif self.infer_ic_from_sbml:
                    formula = str(par).split() # par can sometimes be a float, which would cause an error when splitting
                    for i in range(len(formula)):
                        if (formula[i] in parameters.keys()) and (formula[i] not in
                            list(self._local_pars)+list(self._condition_defined_pars)+list(self._global_pars)):
                            generated_code.extend(bytes('    @variable(m, {} == {})\n'.format(formula[i], self._n_conditions, parameters[formula[i]]), 'utf8'))
                    generated_code.extend(bytes('    @constraint(m, [j in 1:{}], {}[j,1] == {})\n'.format(self._n_conditions, specie, initial_assignments[specie]), 'utf8'))
        for specie in species_interpreted_as_ic:
            generated_code.extend(bytes('    @constraint(m, [j in 1:{}], {}[j,1] == {}[j])\n'.format(self._n_conditions, specie, specie+'_0'), 'utf8'))
        generated_code.extend(bytes('\n', 'utf8'))

        # Write ODEs
        generated_code.extend(bytes('    # Model ODEs\n', 'utf8'))
        generated_code.extend(bytes('    println("Defining ODEs...")\n', 'utf8'))
        patterns = [par for par in self.petab_problem.condition_df.columns if par not in species.keys()]       
        for specie in species: # For every species
            if species[specie]:
                generated_code.extend(bytes('    println("{}")\n'.format(specie), 'utf8'))
                generated_code.extend(bytes('    @NLconstraint(m, [j in 1:{}, k in 1:length(t_sim)-1],\n'.format(self._n_conditions), 'utf8'))
                generated_code.extend(bytes('        {}[j, k+1] == {}[j, k] + ('.format(specie, specie), 'utf8'))
                for (coef, reaction_name) in species[specie]: # For every reaction
                    reaction_formula = ' {}*( {} )'.format(coef, reactions[reaction_name])
                    for pattern in patterns: # Add iterator `j` to condition-defined and local parameters
                        reaction_formula = re.sub('[( ]'+pattern+'[, )]', lambda matchobj: matchobj.group(0)[:-1]+'[j]'+matchobj.group(0)[-1:], reaction_formula) # The matchobject starts with a `(` or ` ` and ends with a `,` ` ` or `)`. I insert `[j]` just before the ending of the matchobject. 
                        # reaction_formula = re.sub(' '+pattern+' ', ' '+pattern.rstrip()+'[j] ', reaction_formula) # Todo: not sure if the tailing whitespace is always in the pattern.
                        # reaction_formula = re.sub(' '+pattern+',', ' '+pattern.rstrip(',')+'[j], ', reaction_formula) # Todo: I think the rstrip is superfluous in the last two lines.
                    for spec in species.keys():
                        tmp_iterator = '[j]'
                        if species[spec]:
                            tmp_iterator = '[j, k+1]'
                        reaction_formula = re.sub('[^a-zA-Z0-9_]'+spec+'[^a-zA-Z0-9_]', lambda matchobj: matchobj.group(0)[:-1]+tmp_iterator+matchobj.group(0)[-1:], reaction_formula)
                    reaction_formula = re.sub('pow', '^', reaction_formula)
                    generated_code.extend(bytes(reaction_formula, 'utf8'))
                generated_code.extend(bytes('     ) * ( t_sim[k+1] - t_sim[k] ) )\n', 'utf8'))
            else:
                generated_code.extend(bytes('    @constraint(m, [j in 1:{}, k in 1:length(t_sim)-1], {}[j, k] == {}[j])\n'.format(self._n_conditions, specie, initial_assignments[specie]), 'utf8'))
        generated_code.extend(bytes('\n', 'utf8'))

        # Write observables
        generated_code.extend(bytes('    # Define observables\n', 'utf8'))
        generated_code.extend(bytes('    println("Defining observables...")\n', 'utf8'))
        for observable in self.petab_problem.observable_df.index:
            min_exp_val = np.min(self.petab_problem.measurement_df.loc[self.petab_problem.measurement_df.loc[:, 'observableId'] == observable, 'measurement'])
            max_exp_val = np.max(self.petab_problem.measurement_df.loc[self.petab_problem.measurement_df.loc[:, 'observableId'] == observable, 'measurement'])
            diff = max_exp_val - min_exp_val
            lb = min_exp_val - 1*diff
            ub = max_exp_val + 1*diff
            ub = 5 # Todo: delete this again
            if self._calling_function == '_execute_case':
                lb = min(0, lb)
                ub = max(16, ub)
            generated_code.extend(bytes('    @variable(m, {} <= {}[j in 1:{}, k in 1:length(t_exp)] <= {})\n'.format(lb, observable, self._n_conditions, ub), 'utf8'))
            formula = self.petab_problem.observable_df.loc[observable, 'observableFormula'].split()
            for i in range(len(formula)):
                if formula[i] in species.keys():
                    formula[i] = formula[i]+'[j, t_sim_to_exp[k]]'
                elif formula[i] in patterns:
                    formula[i] = formula[i]+'[j]'
            formula = ''.join(formula)
            for obs_par_name in set_of_observable_params:
                formula = re.sub(obs_par_name[0], obs_par_name[0]+f'[j{obs_par_name[1]}]', formula)
            generated_code.extend(bytes('    @NLconstraint(m, [j in 1:{}, k in 1:length(t_exp)], {}[j, k] == {})\n'.format(self._n_conditions,observable, formula), 'utf8'))
        generated_code.extend(bytes('\n', 'utf8'))

        # Write objective
        generated_code.extend(bytes('    # Define objective\n', 'utf8'))
        generated_code.extend(bytes('    println("Defining objective...")\n', 'utf8'))
        generated_code.extend(bytes('    @NLobjective(m, Min, ', 'utf8'))
        sums_of_nllhs = []
        for observable in self.petab_problem.observable_df.index:
            
            sigma = self.petab_problem.observable_df.loc[observable, 'noiseFormula']
            for noise_par_name in set_of_noise_params:
                sigma = '( '+re.sub(noise_par_name[0], noise_par_name[0]+f'[j{noise_par_name[1]}]', sigma)+' )'

            scale = 'lin'
            if 'observableTransformation' in self.petab_problem.observable_df.columns:
                scale = self.petab_problem.observable_df.loc[observable, 'observableTransformation']
                        
            noise_distribution = 'normal'
            if 'noiseDistribution' in self.petab_problem.observable_df.columns:
                noise_distribution = self.petab_problem.observable_df.loc[observable, 'noiseDistribution']

            if noise_distribution not in ['normal', 'laplace']:
                warnings.warn('`noiseDistribution is {}, but must be either `laplace` or `normal`.'.format(noise_distribution))
            if scale not in ['lin', 'log', 'log10']:
                warnings.warn('`scale is {}, but must be either `lin`, `log` or `log10`.'.format(scale_distribution))
            
            # if noise_distribution == 'laplace':
            #     sums_of_nllhs.append('sum(log(2*{0}) + abs({1}[j, k]-data[j][k, :{1}])/{0}) for j in 1:{2} for k in 1:length(t_exp))\n'.format(sigma, observable, self._n_conditions))
            # else:
            #     sums_of_nllhs.append('sum(0.5 * (log(2*pi) + log({0}^2) + {0}^(-2) * ({1}[j, k]-data[j][k, :{1}])^2) for j in 1:{2} for k in 1:length(t_exp))\n'.format(sigma, observable, self._n_conditions))

            if noise_distribution == 'normal' and scale == 'lin':
                # nllh = 0.5*log(2*pi*sigma**2) + 0.5*((s-m)/sigma)**2
                sums_of_nllhs.append('sum(0.5 * log(2*pi*({1})^2) + 0.5*(({0}[j, k]-data[j][k, :{0}])/({1}))^2 for j in 1:{2} for k in 1:length(t_exp))\n'.format(observable, sigma, self._n_conditions))
            elif noise_distribution == 'normal' and scale == 'log':
                # nllh = 0.5*log(2*pi*sigma**2*m**2) + 0.5*((log(s)-log(m))/sigma)**2
                sums_of_nllhs.append('sum(0.5*log(2*pi*({1})^2*(data[j][k, :{0}])^2) + 0.5*((log({0}[j, k])-log(data[j][k, :{0}]))/({1}))^2 for j in 1:{2} for k in 1:length(t_exp))\n'.format(observable, sigma, self._n_conditions))
            elif noise_distribution == 'normal' and scale == 'log10':
                # nllh = 0.5*log(2*pi*sigma**2*m**2*log(10)**2) + \
                    # 0.5*((log10(s)-log10(m))/sigma)**2
                sums_of_nllhs.append('sum(0.5*log(2*pi*({1})^2*(data[j][k, :{0}])^2*log(10)^2) + 0.5*((log10({0}[j, k])-log10(data[j][k, :{0}]))/({1}))^2 for j in 1:{2} for k in 1:length(t_exp))\n'.format(observable, sigma, self._n_conditions))
            elif noise_distribution == 'laplace' and scale == 'lin':
                # nllh = log(2*sigma) + abs((s-m)/sigma)
                sums_of_nllhs.append('sum(log(2*{1}) + abs({0}[j, k]-data[j][k, :{0}])/({1})) for j in 1:{2} for k in 1:length(t_exp))\n'.format(observable, sigma, self._n_conditions))
            elif noise_distribution == 'laplace' and scale == 'log':
                # nllh = log(2*sigma*m) + abs((log(s)-log(m))/sigma)
                sums_of_nllhs.append('sum(log(2*{1}*data[j][k, :{0}]) + abs(log({0}[j, k])-log(data[j][k, :{0}]))/({1})) for j in 1:{2} for k in 1:length(t_exp))\n'.format(observable, sigma, self._n_conditions))
            elif noise_distribution == 'laplace' and scale == 'log10':
                # nllh = log(2*sigma*m*log(10)) + abs((log10(s)-log10(m))/sigma)
                sums_of_nllhs.append('sum(log(2*{1}*data[j][k, :{0}]*log(10)) + abs(log10({0}[j, k])-log10(data[j][k, :{0}]))/({1})) for j in 1:{2} for k in 1:length(t_exp))\n'.format(observable, sigma, self._n_conditions))                




        priors = []
        if 'objectivePriorType' in self.petab_problem.parameter_df.columns \
            or 'objectivePriorParameters' in self.petab_problem.parameter_df.columns:
            for index, row in self.petab_problem.parameter_df.iterrows():
                prior_type = row['objectivePriorType']
                prior_parameters = [par.strip() for par in row['objectivePriorParameters'].split(';')]
                if index in self._global_pars.keys():
                    priors.append(str(self._n_conditions)+' * '+NEG_LLH_DICT[prior_type].format(index, prior_parameters[0], prior_parameters[1], '')+'\n')
                else:
                    priors.append('sum('+NEG_LLH_DICT[prior_type].format(index, prior_parameters[0], prior_parameters[1], '[j]')+' for j in 1:{})\n'.format(self._n_conditions))


        generated_code.extend(bytes('        + '.join(sums_of_nllhs+priors), 'utf8'))
        generated_code.extend(bytes('        )\n\n', 'utf8'))

        generated_code.extend(bytes('    println("Optimizing:")\n', 'utf8'))
        generated_code.extend(bytes('    optimize!(m)\n\n', 'utf8'))

        # Write code to get the solution
        julia_pars = list(self._global_pars.keys()) + list(self._local_pars.keys())
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
            generated_code.extend(bytes(specie+', ', 'utf8'))
        generated_code.extend(bytes(']\n', 'utf8'))
        generated_code.extend(bytes('    species_values = Dict()\n', 'utf8'))
        generated_code.extend(bytes('    for s in species_names\n', 'utf8'))
        generated_code.extend(bytes('        species_values[split(string(s[1]), "[")[1]] = JuMP.value.(s)\n', 'utf8')) # Todo: maybe replace `Vector` with `Array`
        generated_code.extend(bytes('    end\n\n', 'utf8'))

        generated_code.extend(bytes('    observable_names = [', 'utf8'))
        for observable in self.petab_problem.observable_df.index:
            generated_code.extend(bytes(observable+', ', 'utf8'))
        generated_code.extend(bytes(']\n', 'utf8'))
        generated_code.extend(bytes('    observable_values = Dict()\n', 'utf8'))
        generated_code.extend(bytes('    for o in observable_names\n', 'utf8'))
        generated_code.extend(bytes('        observable_values[split(string(o[1]), "[")[1]] = Array(JuMP.value.(o))\n', 'utf8'))
        generated_code.extend(bytes('    end\n\n', 'utf8'))

        generated_code.extend(bytes('    objective_val = objective_value(m)\n\n', 'utf8'))
        generated_code.extend(bytes('    results["objective_value"][string(i_start)] = objective_val\n', 'utf8'))
        generated_code.extend(bytes('    results["parameters"][string(i_start)] = parameter_values\n', 'utf8'))
        generated_code.extend(bytes('    results["species"][string(i_start)] = species_values\n', 'utf8'))
        generated_code.extend(bytes('    results["observables"][string(i_start)] = observable_values\n\n', 'utf8'))
        generated_code.extend(bytes('end\n\n', 'utf8'))

        generated_code.extend(bytes('results', 'utf8'))


        # Updating self and files
        code = generated_code.decode()
        self._julia_code = code        
        if self._optimized == True:
            self.optimize()
        if self._files_written == True:
            self.write_jl_file(self._julia_file)
        if self._plotted == True:
            self.plot_results(self._plot_file)