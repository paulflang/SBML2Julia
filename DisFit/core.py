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
import yaml
from julia.api import Julia
importlib.reload(libsbml)


class DisFitProblem(object):

    def __init__(self, petab_yaml, t_ratio=2, n_starts=1):
        """        
        Args:
            petab_yaml (:obj:`str`): path petab .yaml file
            t_ratio (:obj:`int`, optional): number of time discretiation steps per time unit.
            fold_change (:obj:`float`, optional): fold change window of parameter search range wrt sbml parameters
            n_starts (:obj:`int`): number of multistarts
        """
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
            :obj:`int`: number of time discretiation steps per time unit.
        """
        return self._t_ratio

    @t_ratio.setter
    def t_ratio(self, value):
        """Set t_ratio
        
        Args:
            value (:obj:`int`): number of time discretiation steps per time unit.
        
        Raises:
            ValueError: if t_ratio is not an integer >= 1
        """
        if not isinstance(value, int) or (value < 1):
            raise ValueError('`t_ratio` must be an integer >= 1.')
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
            :obj:`dict`: Results in a dict with keys 'states', 'observables', 'x' and 'x_best'
        """
        print('Running optimization problem in julia...')
        out = self._jl.eval(self.julia_code)
        self._results['states'] = out['states']
        self._results['observables'] = out['observables']
        print('Finished optimization in julia.')
        self._best_iter = min(out['objective_val'], key=out['objective_val'].get)
        self._results['x'] = {}
        for i_iter in range(1, self._n_starts+1):
            self._results['x'][str(i_iter)] = {str(k).split()[1].rstrip('>'): v for k, v in out['x'][str(i_iter)].items()}
        x_best = self.results['x'][self._best_iter]

        x_0 = dict(zip(list(self.petab_problem.parameter_df.index), self.petab_problem.parameter_df.loc[:, 'nominalValue']))
        x_best_to_x_0_col = [x_best[key] / x_0[str(key)] for key in x_best.keys()]
        name_col = [str(key) for key in x_best.keys()]
        x_0_col = [x_0[str(key)] for key in x_best.keys()]
        x_best_col = [x_best[str(key)] for key in x_best.keys()]
        self._results['x_best'] = pd.DataFrame(list(zip(name_col, x_0_col, x_best_col,
            x_best_to_x_0_col)), columns = ['Name', 'x_0', 'x_best', 'x_best_to_x_0'])
        self._results['x_best'] = self._results['x_best'].sort_values(by=['Name']).reset_index(drop=True)

        self._optimized = True
        return self.results

    def plot_results(self, path=os.path.join('.', 'plot.pdf'), observables=[], size=(6, 5)):
        """Plot results
        
        Args:
            path (:obj:`str`, optional): path to output plot
            observables (:obj:`list`, optional): list of observables to be plotted
            size (:obj:`tuple`, optional): size of image
        
        Raises:
            ValueError: if `observables` is not a list
        """
        # Options
        x_label = 'time'
        y_label = 'Abundance'
        measurement_df = self.petab_problem.measurement_df.pivot(index='time', columns='observableId')
        t = list(measurement_df.index)
        t_sim = np.linspace(start=0, stop=t[-1], num=t[-1]*self.t_ratio+1)
        if not isinstance(observables, list):
            raise ValueError('`observables` must be a list of observables.')
        if not observables:
            values = pd.DataFrame(self.results['observables'][self._best_iter])
            exp_data = measurement_df.drop(['simulationConditionId'], axis=1)
        else:
            values = pd.DataFrame(self.results['observables'][self._best_iter]).loc[:, observables]
            exp_data = measurement_df.loc[:, variables]

        # Determine the size of the figure
        plt.figure(figsize=size)

        a = plt.axes([0.1, 0.1, 0.8, 0.8])
        a.plot(t_sim, values, linewidth=3)
        a.legend(tuple(values.columns))
        plt.plot(t, exp_data, 'x')
        a.legend(tuple(values.columns)) # Todo: create a legend for experimental data.
        plt.xlim(np.min(t), np.max(t))
        plt.ylim(0, 1.1 * max(values.max()))
        plt.xlabel(x_label, fontsize=18)
        plt.ylabel(y_label, fontsize=18)
        plt.title('DisFit time course')

        plt.savefig(path)
        plt.close()

        self._plotted = True
        self._plot_file = path

    def write_results(self, path=os.path.join('.', 'results.xlsx')):
        """Write results to excel file
        
        Args:
            path (:obj:`str`, optional): path of excel file to write results to.
        """
        with pd.ExcelWriter(path) as writer:  
            self.results['x_best'].to_excel(writer, sheet_name='x_best')

            df = pd.DataFrame(self.results['states'][self._best_iter])
            df.axes[0].names = ['time']
            df.to_excel(writer, sheet_name='states')
            
            df = pd.DataFrame(self.results['observables'][self._best_iter])
            df.axes[0].names = ['time']
            df.to_excel(writer, sheet_name='observables')

    def _set_petab_problem(self, petab_yaml):
        """Converts petab yaml to dict and creates petab.problem.Problem object
        
        Args:
            petab_yaml (:obj:`str`): path petab .yaml file
        
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

        mod = self.petab_problem.sbml_model
        parameter_df = self.petab_problem.parameter_df
        condition_df = self.petab_problem.condition_df
        observable_df = self.petab_problem.observable_df
        measurement_df = self.petab_problem.measurement_df

        species_file = os.path.join(self._petab_dirname, self.petab_yaml_dict['problems'][0]['species_files'][0])
        species_df = pd.read_csv(species_file, sep='\t', index_col='speciesId')

        observableIds = list(observable_df.index)

        n_params = mod.getNumParameters()
        
        variables = {}
        for i in range(mod.getNumSpecies()): 
            species = mod.getSpecies(i)
            if species.getBoundaryCondition() == True or (species.getId() in variables):
                continue
            variables[species.getId()] = []
        self._var_names = list(variables.keys())

        par_values = []
        par_names = []
        for i in range(mod.getNumParameters()):
            element = mod.getParameter(i)
            par_values.append(element.getValue())
            par_names.append(element.getId())
        par_values_string = str(par_values)
        par_names_string = str(par_names)
        self._par_values = par_values
        self._par_names = par_names
        
        x_0 = []
        for variable in variables:
            # get initialValue 
            element = mod.getElementBySId(variable)
            if element.getTypeCode() == libsbml.SBML_PARAMETER: 
                x_0.append(element.getValue())
            elif element.getTypeCode() == libsbml.SBML_SPECIES: 
                if element.isSetInitialConcentration(): 
                    x_0.append(element.getInitialConcentration())
                else: 
                    x_0.append(element.getInitialAmount())
            else: 
                x_0.append(element.getSize())
        n_x_0 = len(x_0)
        x_0_string = str(x_0)

        # start generating the code by appending to bytearray
        generated_code = bytearray('', 'utf8')
        
        
        generated_code.extend(bytes('using CSV\n', 'utf8'))
        generated_code.extend(bytes('using DataFrames\n', 'utf8'))
        generated_code.extend(bytes('using Ipopt\n', 'utf8'))
        generated_code.extend(bytes('using JuMP\n', 'utf8'))

        generated_code.extend(bytes('\n', 'utf8'))
        # generated_code.extend(bytes('fc = {} # Setting parameter search span\n'.format(self.fold_change), 'utf8'))
        generated_code.extend(bytes('t_ratio = {} # Setting number of ODE discretisation steps\n'.format(self.t_ratio), 'utf8'))

        generated_code.extend(bytes('\n', 'utf8'))  
        generated_code.extend(bytes('# Data\n', 'utf8'))
        generated_code.extend(bytes('data_path = "{}"\n'.format(os.path.join(self._petab_dirname, self.petab_yaml_dict['problems'][0]['measurement_files'][0])), 'utf8'))
        generated_code.extend(bytes('df = CSV.read(data_path)\n', 'utf8'))

        # Todo Sungho: I did not manage to read in the data correctly.
        # Could you please make sure that experimental data is read into Julia in a way that is compatible with the objective function definition?
        # Note that `df` is in molten/long form (see DisFit/tests/fixtures/G2M_copasi/measurementData_G2M_copasi.tsv).
        # As of now, I would have thought of accessing observable_x of condition j and time point k as `data[j][k, :observable_x]` (see a few lines below).
        generated_code.extend(bytes('dfg = groupby(df, :simulationConditionId)\n', 'utf8'))
        generated_code.extend(bytes('data = [DataFrame() for i in length(dfg)]\n', 'utf8'))
        generated_code.extend(bytes('data[i] = unstack(dfg[condition], :time, :observableId, :measurement)\n', 'utf8'))
        generated_code.extend(bytes('i = 1\n', 'utf8'))
        generated_code.extend(bytes('for condition in keys(dfg)\n', 'utf8'))
        generated_code.extend(bytes('    data[i] = unstack(dfg[condition], :time, :observableId, :measurement)\n', 'utf8'))
        generated_code.extend(bytes('    i = i+1\n', 'utf8'))
        generated_code.extend(bytes('end\n', 'utf8'))
        generated_code.extend(bytes('\n', 'utf8'))
        generated_code.extend(bytes('t_exp = Vector(DataFrame(groupby(dfg[1], :observableId)[1])[!, :time])\n', 'utf8'))


        generated_code.extend(bytes('t_sim = range(0, stop=t_exp[end], length=t_exp[end]*t_ratio+1)\n\n', 'utf8'))

        generated_code.extend(bytes('results = Dict()\n', 'utf8'))
        generated_code.extend(bytes('results["objective_val"] = Dict()\n', 'utf8'))
        generated_code.extend(bytes('results["x"] = Dict()\n', 'utf8'))
        generated_code.extend(bytes('results["states"] = Dict()\n', 'utf8'))
        generated_code.extend(bytes('results["observables"] = Dict()\n', 'utf8'))
        generated_code.extend(bytes('for i_start in 1:{}\n'.format(self._n_starts), 'utf8'))  
        generated_code.extend(bytes('    m = Model(with_optimizer(Ipopt.Optimizer))\n\n', 'utf8'))
        # i = 0
        
        generated_code.extend(bytes('    # Define condition-defined parameters\n', 'utf8'))
        n_conditions = condition_df.shape[0]
        for parameter in condition_df.columns:
            if str(condition_df[parameter].dtype) in ('float64', 'int16'):
                generated_code.extend(bytes('    @variable(m, {0}[1:{1}])\n'.format(parameter, n_conditions), 'utf8'))
                for i in range(1, n_conditions+1):
                    generated_code.extend(bytes('    @constraint(m, {0}[{1}] == {2})\n'.format(parameter, i, condition_df.iloc[i-1][parameter]), 'utf8'))
        generated_code.extend(bytes('\n', 'utf8'))

        local_pars = []
        generated_code.extend(bytes('    # Define condition-local parameters\n', 'utf8'))
        for parameter in condition_df.columns:
            if str(condition_df[parameter].dtype) == 'object':
                generated_code.extend(bytes('    @variable(m, {0}[1:{1}])\n'.format(parameter, n_conditions), 'utf8'))
                for i in range(1, n_conditions+1):
                    local_pars.append(condition_df.iloc[i-1][parameter])
                    correct_row = (parameter_df.index == condition_df.iloc[i-1][parameter])
                    lb = parameter_df.loc[correct_row, 'lowerBound']
                    ub = parameter_df.loc[correct_row, 'upperBound']
                    generated_code.extend(bytes('    @constraint(m, {} <= {}[{}] <= {})\n'.format(lb, parameter, i, ub), 'utf8'))
        generated_code.extend(bytes('\n', 'utf8'))

        generated_code.extend(bytes('    # Define global parameters\n', 'utf8'))
        for element in parameter_df.index:
            if element not in local_pars:
                lb = parameter_df.loc[element, 'lowerBound']
                ub = parameter_df.loc[element, 'upperBound']
                nominal = parameter_df.loc[element, 'nominalValue']
                generated_code.extend(bytes('    @variable(m, {0} <= {1} <= {2}, start={0}+({2}-{0})*rand(Float64))\n'.format(lb, str(element), ub), 'utf8'))

        generated_code.extend(bytes('\n', 'utf8'))
        generated_code.extend(bytes('    # Model states\n', 'utf8'))
        generated_code.extend(bytes('    println("Defining states ...")\n', 'utf8'))
        for variable in variables.keys():
            lb = species_df.loc[variable, 'lowerBound'] #Todo: write somhere a linter that check that the set of sbml model species == species_df.index
            ub = species_df.loc[variable, 'upperBound']
            generated_code.extend(bytes('    @variable(m, {} <= {}[j in 1:{}, k in 1:length(t_sim)] <= {})\n'.format(lb, variable, n_conditions, ub), 'utf8'))

        generated_code.extend(bytes('\n', 'utf8'))
        generated_code.extend(bytes('    # Model ODEs\n', 'utf8'))
        generated_code.extend(bytes('    println("Defining ODEs ...")\n', 'utf8'))
        
        reactions = {}
        for i in range(mod.getNumReactions()):
            reaction = mod.getReaction(i)
            kinetics = reaction.getKineticLaw()
            kinetic_components = kinetics.getFormula().split(' * ')[1:]
            for i in range(len(kinetic_components[1:])):
                kinetic_components[i+1] = kinetic_components[i+1] + '[k+1]'
            jump_formula = ' * '.join(kinetic_components)
            reactions[reaction.getId()] = jump_formula
        
        for i in range(mod.getNumReactions()): 
            reaction = mod.getReaction(i)
            kinetics = reaction.getKineticLaw()   
            for j in range(reaction.getNumReactants()): 
                ref = reaction.getReactant(j)
                species = mod.getSpecies(ref.getSpecies())
                products = [r.getSpecies() for r in reaction.getListOfProducts()]
                if (species.getBoundaryCondition() == True) or (species.getName() in products):
                    # print('continueing...')
                    continue
                variables[species.getId()].append(('-'+str(ref.getStoichiometry()), reaction.getId()))
                # print('added reaction {} to species {}'.format(reaction.getID(), species))
            for j in range(reaction.getNumProducts()): 
                ref = reaction.getProduct(j)
                species = mod.getSpecies(ref.getSpecies())
                reactants = [r.getSpecies() for r in reaction.getListOfReactants()]
                if (species.getBoundaryCondition() == True) or (species.getName() in reactants): 
                    continue
                variables[species.getId()].append((('+'+str(ref.getStoichiometry()), reaction.getId())))

        patterns = [par+' ' for par in condition_df.columns]
        for variable in variables:
            generated_code.extend(bytes('    @NLconstraint(m, [j in 1:{}, k in 1:length(t_sim)-1],\n'.format(n_conditions), 'utf8'))
            generated_code.extend(bytes('        {}[j, k+1] == {}[j, k] + ('.format(variable, variable), 'utf8'))
            for (coef, reaction_name) in variables[variable]:
                reaction_formula = ' {}*({})'.format(coef, reactions[reaction_name])
                for pattern in patterns:
                    reaction_formula = re.sub(pattern, pattern.rstrip()+'[j] ', reaction_formula) # Todo: not sure if the tailing whitespace is always in the pattern.
                generated_code.extend(bytes(reaction_formula, 'utf8'))
            generated_code.extend(bytes('     ) * ( t_sim[k+1] - t_sim[k] ) )\n', 'utf8'))
        generated_code.extend(bytes('\n', 'utf8'))

        # Define observables
        generated_code.extend(bytes('    # Define observables\n', 'utf8'))
        generated_code.extend(bytes('    println("Defining observables ...")\n', 'utf8'))
        for observable in observableIds:
            min_exp_val = np.min(measurement_df.loc[measurement_df.loc[:, 'observableId'] == observable, 'measurement'])
            max_exp_val = np.max(measurement_df.loc[measurement_df.loc[:, 'observableId'] == observable, 'measurement'])
            diff = max_exp_val - min_exp_val
            lb = min_exp_val - 0.2*diff
            ub = max_exp_val + 0.2*diff
            generated_code.extend(bytes('    @variable(m, {} <= {}[j in 1:{}, k in 1:length(t_sim)] <= {})\n'.format(lb, observable, n_conditions, ub), 'utf8'))
            formula = observable_df.loc[observable, 'observableFormula'].split()
            for i in range(len(formula)):
                if formula[i] in self._var_names:
                    formula[i] = formula[i]+'[k]'
                elif formula[i] in patterns:
                    formula[i] = formula[i]+'[j]'
            formula = ''.join(formula)
            generated_code.extend(bytes('    @NLconstraint(m, [j in 1:{}, k in 1:length(t_sim)], {}[k] == {})\n'.format(n_conditions,observable, formula), 'utf8'))
        generated_code.extend(bytes('\n', 'utf8'))

        # Define objective
        generated_code.extend(bytes('    # Define objective\n', 'utf8'))
        generated_code.extend(bytes('    println("Defining objective ...")\n', 'utf8'))
        generated_code.extend(bytes('    @NLobjective(m, Min,', 'utf8'))
        sums_of_squares = []
        for observable in observableIds:
            sums_of_squares.append('sum(({0}[j, (k-1)*t_ratio+1]-data[j][k, :{0}])^2 for j in 1:{1} for k in 1:length(t_exp))\n'.format(observable, n_conditions))
        generated_code.extend(bytes('        + '.join(sums_of_squares), 'utf8'))
        generated_code.extend(bytes('        )\n\n', 'utf8'))

        generated_code.extend(bytes('    println("Optimizing...")\n', 'utf8'))
        generated_code.extend(bytes('    optimize!(m)\n\n', 'utf8'))

        # Retreiving the solution
        generated_code.extend(bytes('    println("Retreiving solution...")\n', 'utf8'))
        # generated_code.extend(bytes('    species_to_plot = {}\n'.format(species_to_plot), 'utf8'))
        generated_code.extend(bytes('    params = ' + str(list(parameter_df.index)).replace('\'', '') + '\n', 'utf8'))
        generated_code.extend(bytes('    paramvalues = Dict()\n', 'utf8'))
        generated_code.extend(bytes('    for param in params\n', 'utf8'))
        generated_code.extend(bytes('        paramvalues[param] = JuMP.value.(param)\n', 'utf8'))
        generated_code.extend(bytes('    end\n\n', 'utf8'))

        generated_code.extend(bytes('    variables = [', 'utf8'))
        for variable in variables:
            generated_code.extend(bytes(variable+', ', 'utf8'))
        generated_code.extend(bytes(']\n', 'utf8'))
        generated_code.extend(bytes('    variablevalues = Dict()\n', 'utf8'))
        generated_code.extend(bytes('    for v in variables\n', 'utf8'))
        generated_code.extend(bytes('        variablevalues[string(v[1])[1:end-3]] = Vector(JuMP.value.(v))\n', 'utf8')) # Todo: maybe replace `Vector` with `Array`
        generated_code.extend(bytes('    end\n\n', 'utf8'))

        generated_code.extend(bytes('    observables = [', 'utf8'))
        for observable in observableIds:
            generated_code.extend(bytes(observable+', ', 'utf8'))
        generated_code.extend(bytes(']\n', 'utf8'))
        generated_code.extend(bytes('    observablevalues = Dict()\n', 'utf8'))
        generated_code.extend(bytes('    for o in observables\n', 'utf8'))
        generated_code.extend(bytes('        observablevalues[string(o[1])[1:end-3]] = Array(JuMP.value.(o))\n', 'utf8'))
        generated_code.extend(bytes('    end\n\n', 'utf8'))

        generated_code.extend(bytes('    v = objective_value(m)\n\n', 'utf8'))
        generated_code.extend(bytes('    results["objective_val"][string(i_start)] = v\n', 'utf8'))
        generated_code.extend(bytes('    results["x"][string(i_start)] = paramvalues\n', 'utf8'))
        generated_code.extend(bytes('    results["states"][string(i_start)] = variablevalues\n', 'utf8'))
        generated_code.extend(bytes('    results["observables"][string(i_start)] = observablevalues\n\n', 'utf8'))
        generated_code.extend(bytes('end\n\n', 'utf8'))

        generated_code.extend(bytes('results', 'utf8'))

        code = generated_code.decode()
        self._julia_code = code
        
        # Updating self and files if needed
        if self._optimized == True:
            self.optimize()
        if self._files_written == True:
            self.write_jl_file(self._julia_file)
        if self._plotted == True:
            self.plot_results(self._plot_file)