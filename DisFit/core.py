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
import scipy as sp
import sys
import tempfile
from julia.api import Julia
importlib.reload(libsbml)


class DisFitProblem(object):

    def __init__(self, sbml_path, data_path, t_ratio=2, fold_change=2, n_starts=1):
        """        
        Args:
            sbml_path (:obj:`str`): path to sbml file
            data_path (:obj:`str`): path to data file
            t_ratio (:obj:`int`, optional): number of time discretiaation steps
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
        self.sbml_path = sbml_path
        self.data_path = data_path
        self.t_ratio = t_ratio
        self.fold_change = fold_change
        self.n_starts = n_starts
        self._set_julia_code()
        self._initialization = False

    @property
    def sbml_path(self):
        """Get sbml path
        
        Returns:
            :obj:`str`: path to `sbml` file
        """
        return self._sbml_path

    @sbml_path.setter
    def sbml_path(self, value):
        """Set sbml path
        
        Args:
            value (:obj:`str`): path to `sbml` file
        
        Raises:
            ValueError: if sbml_path is not a sbml file
        """
        if not isinstance(value, str) or not (value.endswith('.sbml') or value.endswith('.xml')):
            raise ValueError('`sbml_path` must be a path to a sbml file')
        if not os.path.isfile(value):
            raise ValueError('Cannot find file `{}`'.format(value))
        self._sbml_path = value
        if not self._initialization:
            self._set_julia_code()

    @property
    def data_path(self):
        """Get data path
        
        Returns:
            :obj:`str`: path to data file
        """
        return self._data_path

    @data_path.setter
    def data_path(self, value):
        """Set data path
        
        Args:
            value (:obj:`str`): path to data file
        
        Raises:
            ValueError: if data_path is not a csv file
        """
        if not isinstance(value, str) or not (value.endswith('.csv')):
            raise ValueError('`data_path` must be a path to a csv file')
        if not os.path.isfile(value):
            raise ValueError('Cannot find file `{}`'.format(value))
        df = pd.read_csv(value)
        self._data_path = value
        self._exp_data = df
        if not self._initialization:
            self._set_julia_code()

    @property
    def t_ratio(self):
        """Get t_ratio
        
        Returns:
            :obj:`int`: ratio between experimental observation intervals and time
            discretization intervals
        """
        return self._t_ratio

    @t_ratio.setter
    def t_ratio(self, value):
        """Set t_ratio
        
        Args:
            value (:obj:`int`): number of time discretization steps
        
        Raises:
            ValueError: if t_ratio is not an integer >= 1
        """
        if not isinstance(value, int) or (value < 1):
            raise ValueError('`t_ratio` must be an integer >= 1.')
        self._t_ratio = value
        if not self._initialization:
            self._set_julia_code()

    @property
    def fold_change(self):
        """Get fold_change
        
        Returns:
            :obj:`float`: fold change window of parameter search range wrt sbml parameters
        """
        return self._fold_change

    @fold_change.setter
    def fold_change(self, value):
        """Set fold_change
        
        Args:
            value (:obj:`float`): fold change window of parameter search range wrt sbml parameters
        
        Raises:
            ValueError: if fold_change <= 1
        """
        if not (isinstance(value, float) or isinstance(value, int)) or (value <= 1):
            raise ValueError('`fold_change` must be a float > 1')
        self._fold_change = float(value)
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
    def exp_data(self):
        """Get exp_data
        
        Returns:
            :obj:`pandas.DataFrame`: experimental data
        """
        return self._exp_data

    @property
    def results(self):
        """Get results
        
        Returns:
            :obj:`dict`: optimization results
        """
        return self._results

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
            :obj:`dict`: Results in a dict with keys 'states', 'x' and 'x_best'
        """
        print('Running optimization problem in julia...')
        out = self._jl.eval(self.julia_code)
        self._results['states'] = out['states']
        print('Finished optimization in julia.')
        self._best_iter = min(out['objective_val'], key=out['objective_val'].get)
        self._results['x'] = {}
        for i_iter in range(1, self._n_starts+1):
            self._results['x'][str(i_iter)] = {str(k).split()[1].rstrip('>'): v for k, v in out['x'][str(i_iter)].items()}
        x_best = self.results['x'][self._best_iter]

        x_0 = dict(zip(self._par_names, self._par_values))
        x_best_to_x_0_col = [x_best[key] / x_0[str(key)] for key in x_best.keys()]
        name_col = [str(key) for key in x_best.keys()]
        x_0_col = [x_0[str(key)] for key in x_best.keys()]
        x_best_col = [x_best[str(key)] for key in x_best.keys()]
        self._results['x_best'] = pd.DataFrame(list(zip(name_col, x_0_col, x_best_col,
            x_best_to_x_0_col)), columns = ['Name', 'x_0', 'x_best', 'x_best_to_x_0'])
        self._results['x_best'] = self._results['x_best'].sort_values(by=['Name']).reset_index(drop=True)

        self._optimized = True
        return self.results

    def plot_results(self, path=os.path.join('.', 'plot.pdf'), variables=[], size=(6, 5)):
        """Plot results
        
        Args:
            path (:obj:`str`, optional): path to output plot
            variables (:obj:`list`, optional): list of variables to be plotted
            size (:obj:`tuple`, optional): size of image
        
        Raises:
            ValueError: if `variables` is not a list
        """
        # Options
        x_label = 'time'
        y_label = 'Abundance'
        t = self.exp_data.loc[:, 't'].values
        t_sim = np.linspace(start=0, stop=t[-1], num=t[-1]*self.t_ratio+1)
        if not isinstance(variables, list):
            raise ValueError('`variables` must be a list of variables.')
        if not variables:
            values = pd.DataFrame(self.results['states'][self._best_iter])
            var_names = self._var_names
            exp_data = self._exp_data
        else:
            values = pd.DataFrame(self.results['states'][self._best_iter]).loc[:, variables]
            var_names = pd.DataFrame(self.results['states'][self._best_iter]).loc[:, variables].index
            exp_data = self.exp_data.loc[:, variables]

        # Determine the size of the figure
        plt.figure(figsize=size)

        a = plt.axes([0.1, 0.1, 0.8, 0.8])
        a.plot(t_sim, values, linewidth=3)
        a.legend(tuple(values.columns))
        plt.plot(t, exp_data, 'x')
        a.legend(tuple(values.columns))
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
            pd.DataFrame(self.results['states'][self._best_iter]).to_excel(writer, sheet_name='states')

    def _set_julia_code(self):
        """Transform sbml file to Julia JuMP model.
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

        doc = libsbml.readSBMLFromFile(self.sbml_path)
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
        props.addOption("expandInitialAssignments", True)
         
        if doc.convert(props) != libsbml.LIBSBML_OPERATION_SUCCESS: 
            print('The document could not be converted')
            print(doc.getErrorLog().toString())
           
        props = libsbml.ConversionProperties()
        props.addOption("expandFunctionDefinitions", True)
         
        if doc.convert(props) != libsbml.LIBSBML_OPERATION_SUCCESS: 
            print('The document could not be converted')
            print(doc.getErrorLog().toString())

        # figure out which species are variable 
        mod = doc.getModel()

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
        generated_code.extend(bytes('using JuMP\n', 'utf8'))
        generated_code.extend(bytes('using Ipopt\n', 'utf8'))
        generated_code.extend(bytes('using CSV\n', 'utf8'))

        generated_code.extend(bytes('\n', 'utf8'))
        generated_code.extend(bytes('fc = {} # Setting parameter search span\n'.format(self.fold_change), 'utf8'))
        generated_code.extend(bytes('t_ratio = {} # Setting number of ODE discretisation steps\n'.format(self.t_ratio), 'utf8'))

        generated_code.extend(bytes('\n', 'utf8'))  
        generated_code.extend(bytes('# Data\n', 'utf8'))
        generated_code.extend(bytes('data_path = "{}"\n'.format(self.data_path), 'utf8'))
        generated_code.extend(bytes('df = CSV.read(data_path)\n', 'utf8'))
        generated_code.extend(bytes('t_exp = Vector(df[!, :t]) # set of simulation times.)\n', 'utf8'))
        generated_code.extend(bytes('t_sim = range(0, stop=t_exp[end], length=t_exp[end]*t_ratio+1)\n\n', 'utf8'))

        generated_code.extend(bytes('results = Dict()\n', 'utf8'))
        generated_code.extend(bytes('results["objective_val"] = Dict()\n', 'utf8'))
        generated_code.extend(bytes('results["x"] = Dict()\n', 'utf8'))
        generated_code.extend(bytes('results["states"] = Dict()\n', 'utf8'))
        generated_code.extend(bytes('for i_start in 1:{}\n'.format(self._n_starts), 'utf8'))  
        generated_code.extend(bytes('m = Model(with_optimizer(Ipopt.Optimizer))\n\n', 'utf8'))
        i = 0
        for i in range(mod.getNumParameters()):
            element = mod.getParameter(i)
            generated_code.extend(bytes('    @variable(m, {0}/fc <= {1} <= {0}*fc, start={0}/fc+({0}*fc-{0}/fc)*rand(Float64))\n'.format(par_values[i], element.getId()), 'utf8'))
            i =+ 1

        generated_code.extend(bytes('\n', 'utf8'))
        generated_code.extend(bytes('    # Model states\n', 'utf8'))
        generated_code.extend(bytes('    println("Defining states ...")\n', 'utf8'))
        for variable in variables.keys():
            generated_code.extend(bytes('    @variable(m, 0 <= {0}[k in 1:length(t_sim)] <= 1)\n'.format(variable), 'utf8'))

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

        for variable in variables:
            generated_code.extend(bytes('    @NLconstraint(m, [k in 1:length(t_sim)-1],\n', 'utf8'))
            generated_code.extend(bytes('        {}[k+1] == {}[k] + ('.format(variable, variable), 'utf8'))
            for (coef, reaction_name) in variables[variable]:
                reaction_formula = ' {}*({})'.format(coef, reactions[reaction_name])
                generated_code.extend(bytes(reaction_formula, 'utf8'))
            generated_code.extend(bytes('     ) * ( t_sim[k+1] - t_sim[k] ) )\n', 'utf8'))
        generated_code.extend(bytes('\n', 'utf8'))

        generated_code.extend(bytes('    # Define objective\n', 'utf8'))
        generated_code.extend(bytes('    println("Defining objective ...")\n', 'utf8'))
        generated_code.extend(bytes('    @NLobjective(m, Min,', 'utf8'))
        sums_of_squares = []
        for variable in variables:
            sums_of_squares.append('sum(({}[(k-1)*t_ratio+1]-df[k, :{}])^2 for k in 1:length(t_exp))\n'.format(variable, variable))
        generated_code.extend(bytes('        + '.join(sums_of_squares), 'utf8'))
        generated_code.extend(bytes('        )\n\n', 'utf8'))

        generated_code.extend(bytes('    println("Optimizing...")\n', 'utf8'))
        generated_code.extend(bytes('    optimize!(m)\n\n', 'utf8'))

        # Retreiving the solution
        generated_code.extend(bytes('    println("Retreiving solution...")\n', 'utf8'))
        # generated_code.extend(bytes('    species_to_plot = {}\n'.format(species_to_plot), 'utf8'))
        generated_code.extend(bytes('    params = ' + str(par_names).replace('\'', '') + '\n', 'utf8'))
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
        generated_code.extend(bytes('        variablevalues[string(v[1])[1:end-3]] = Vector(JuMP.value.(v))\n', 'utf8'))
        generated_code.extend(bytes('    end\n\n', 'utf8'))

        generated_code.extend(bytes('    v = objective_value(m)\n\n', 'utf8'))
        generated_code.extend(bytes('    results["objective_val"][string(i_start)] = v\n', 'utf8'))
        generated_code.extend(bytes('    results["x"][string(i_start)] = paramvalues\n', 'utf8'))
        generated_code.extend(bytes('    results["states"][string(i_start)] = variablevalues\n\n', 'utf8'))
        generated_code.extend(bytes('end\n\n', 'utf8'))

        generated_code.extend(bytes('results\n\n', 'utf8'))

        code = generated_code.decode()
        self._julia_code = code
        
        # Updating self and files if needed
        if self._optimized == True:
            self.optimize()
        if self._files_written == True:
            self.write_jl_file(self._julia_file)
        if self._plotted == True:
            self.plot_results(self._plot_file)