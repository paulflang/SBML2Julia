""" Class to fit PEtab problem via time discretization in Julia
:Author: Paul Lang <paul.lang@wolfson.ox.ac.uk>
:Date: 2020-04-15
"""

import importlib
import libsbml
import numpy as np
import os
import scipy as sp
import shutil
import sys
import subprocess
import tempfile
import time as timer
import warnings
from scipy.signal import find_peaks
from julia.api import Julia
importlib.reload(libsbml)


class DisFitProblem(object):

    def __init__(self, sbml_path, data_path, t_ratio=2, fold_change=2):
        """        
        Args:
            sbml_path (:obj:`str`): path to sbml file
            data_path (:obj:`str`): path to data file
            t_ratio (:obj:`int`, optional): number of time discretiaation steps
            fold_change (:obj:`float`, optional): fold change window of parameter search range wrt sbml parameters
        """
        self._initialization = True
        self._jl = Julia(compiled_modules=False)
        self._initialization = True
        self.sbml_path = sbml_path
        self.data_path = data_path
        self.t_ratio = t_ratio
        self.fold_change = fold_change
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
        if not isinstance(value, str) or not value.endswith('.csv'):
            raise ValueError('`data_path` must be a path to a csv file')
        if not os.path.isfile(value):
            raise ValueError('Cannot find file `{}`'.format(value))
        self._data_path = value
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
        self._fold_change = value
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

    def julia_code_to_file(self):
        pass

    def optimize(self):
        self._results = {}
        print('Running optimization problem in julia')
        optimizer = self._jl.eval(self.julia_code['optimizer'])
        print('Retreiving parameters')
        self._results['parameters'] = self._jl.eval(self.julia_code['parameters'])
        print('Retreiving values')
        self._results['values'] = self._jl.eval(self.julia_code['values'])
        return(self.results)

    def plot_results(self):
        pass

    def _set_julia_code(self):
        #----------------------------------------------------------------------#
        """ `_set_julia_code` is adapted from Frank T. Bergman
        Date: 2019
        Availability: https://groups.google.com/forum/#!topic/sbml-discuss/inS4Lzp3Ri8 or
        https://www.dropbox.com/s/2bfpiausejp0gd0/convert_reactions.py?dl=0 """
        #----------------------------------------------------------------------#

        self._julia_code = {}
        doc = libsbml.readSBMLFromFile(self.sbml_path)
        if doc.getNumErrors(libsbml.LIBSBML_SEV_FATAL):
            print('Encountered serious errors while reading file')
            print(doc.getErrorLog().toString())
            sys.exit(1)
           
        # clear errors
        doc.getErrorLog().clearLog()
         
        #
        # perform conversions
        #
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
           
        #
        # figure out which species are variable 
        #
        mod = doc.getModel()

        n_params = mod.getNumParameters()
        
        variables = {}
        for i in range(mod.getNumSpecies()): 
            species = mod.getSpecies(i)
            if species.getBoundaryCondition() == True or (species.getId() in variables):
                continue
            variables[species.getId()] = []

        par_values = []
        par_names = []
        for i in range(mod.getNumParameters()):
            element = mod.getParameter(i)
            par_values.append(element.getValue())
            par_names.append(element.getId())
        par_values_string = str(par_values)
        par_names_string = str(par_names)
        
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

        #
        # start generating the code by appending to bytearray
        #
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
        generated_code.extend(bytes('t_sim = range(0, stop=t_exp[end], length=t_exp[end]*t_ratio+1)\n', 'utf8'))
        
        generated_code.extend(bytes('\n', 'utf8'))
        generated_code.extend(bytes('# create JuMP model object")\n', 'utf8'))
        generated_code.extend(bytes('# m = Model(with_optimizer(KNITRO.Optimizer, ms_enable=1, ms_maxtime_real=10))\n', 'utf8'))
        generated_code.extend(bytes('m = Model(with_optimizer(Ipopt.Optimizer))\n', 'utf8'))

        generated_code.extend(bytes('\n', 'utf8'))
        generated_code.extend(bytes('# Model parameters\n', 'utf8'))
        generated_code.extend(bytes('println("Defining parameters ...")\n', 'utf8'))  
        i = 0
        for i in range(mod.getNumParameters()):
            element = mod.getParameter(i)
            generated_code.extend(bytes('@variable(m, {0}/fc <= {1} <= {0}*fc)\n'.format(par_values[i], element.getId()), 'utf8'))
            i =+ 1

        generated_code.extend(bytes('\n', 'utf8'))
        generated_code.extend(bytes('# Model states\n', 'utf8'))
        generated_code.extend(bytes('println("Defining states ...")\n', 'utf8'))
        for variable in variables.keys():
            generated_code.extend(bytes('@variable(m, 0 <= {0}[k in 1:length(t_sim)] <= 1)\n'.format(variable), 'utf8'))

        generated_code.extend(bytes('\n', 'utf8'))
        generated_code.extend(bytes('# Model ODEs\n', 'utf8'))
        generated_code.extend(bytes('println("Defining ODEs ...")\n', 'utf8'))
        
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
            generated_code.extend(bytes('@NLconstraint(m, [k in 1:length(t_sim)-1],\n', 'utf8'))
            generated_code.extend(bytes('    {}[k+1] == {}[k] + ('.format(variable, variable), 'utf8'))
            for (coef, reaction_name) in variables[variable]:
                reaction_formula = ' {}*({})'.format(coef, reactions[reaction_name])
                generated_code.extend(bytes(reaction_formula, 'utf8'))
            generated_code.extend(bytes(' ) * ( t_sim[k+1] - t_sim[k] ) )\n', 'utf8'))
        generated_code.extend(bytes('\n', 'utf8'))

        generated_code.extend(bytes('# Define objective\n', 'utf8'))
        generated_code.extend(bytes('println("Defining objective ...")\n', 'utf8'))
        generated_code.extend(bytes('@NLobjective(m, Min,', 'utf8'))
        sums_of_squares = []
        for variable in variables:
            sums_of_squares.append('sum(({}[(k-1)*t_ratio+1]-df[k, :{}])^2 for k in 1:length(t_exp))\n'.format(variable, variable))
        generated_code.extend(bytes('    + '.join(sums_of_squares), 'utf8'))
        generated_code.extend(bytes('    )\n\n', 'utf8'))

        generated_code.extend(bytes('println("Optimizing...")\noptimize!(m)\n\n\n', 'utf8'))
        code = generated_code.decode()
        self._julia_code['optimizer'] = code

        # Retreiving the solution
        generated_code = bytearray('', 'utf8')
        generated_code.extend(bytes('println("# Obtain the solution")\n', 'utf8'))
        generated_code.extend(bytes('println("Retreiving solution...")\n', 'utf8'))
        # generated_code.extend(bytes('species_to_plot = {}\n'.format(species_to_plot), 'utf8'))
        generated_code.extend(bytes('params = ' + str(par_names).replace('\'', '') + '\n', 'utf8'))
        generated_code.extend(bytes('paramvalues = Dict()\n', 'utf8'))
        generated_code.extend(bytes('for param in params\n', 'utf8'))
        generated_code.extend(bytes('    paramvalues[param] = JuMP.value.(param)\n', 'utf8'))
        generated_code.extend(bytes('end\n\n', 'utf8'))
        generated_code.extend(bytes('paramvalues\n\n\n', 'utf8'))
        code = generated_code.decode()
        self._julia_code['parameters'] = code

        generated_code = bytearray('', 'utf8')
        generated_code.extend(bytes('variables = [', 'utf8'))
        for variable in variables:
            generated_code.extend(bytes(' :'+variable+'\n', 'utf8'))
        generated_code.extend(bytes(']\n', 'utf8'))
        generated_code.extend(bytes('variablevalues = Dict()\n', 'utf8'))
        generated_code.extend(bytes('for v in variables\n', 'utf8'))
        generated_code.extend(bytes('    variablevalues[string(v)] = Vector(JuMP.value.(eval(v)))\n', 'utf8'))
        generated_code.extend(bytes('end\n\n', 'utf8'))

        # generated_code.extend(bytes('data_matrix = zeros(length(t_sim), length(species_to_plot))\n', 'utf8'))
        # generated_code.extend(bytes('i = 1\n', 'utf8'))
        # generated_code.extend(bytes('for s in species_to_plot\n', 'utf8'))
        # generated_code.extend(bytes('    data_matrix[:, i] = variablevalues[string(s)]\n', 'utf8'))
        # generated_code.extend(bytes('i = i+1\n', 'utf8'))
        # generated_code.extend(bytes('end\n\n\n', 'utf8'))
        generated_code.extend(bytes('variablevalues\n\n\n', 'utf8'))
        code = generated_code.decode()
        self._julia_code['values'] = code

        # # Ploting
        # generated_code.extend(bytes('using Plots\n', 'utf8'))
        # generated_code.extend(bytes('using DataFrames\n', 'utf8'))
        # generated_code.extend(bytes('using StatPlots\n', 'utf8'))
        # generated_code.extend(bytes('p = plot(t_sim, data_matrix, xlabel="Time (hr)", ylabel="Abundance", label=species_to_plot\n, legend=:topleft)\n\n', 'utf8'))
        # generated_code.extend(bytes('@df df plot!(p, t_exp, seriestype=:scatter, markershape = :x,\n', 'utf8'))
        # generated_code.extend(bytes('    markersize = 2,\n    markerstrokewidth = 0.1,\n    markerstrokealpha = 1)', 'utf8'))

        # code = generated_code.decode()
        # self._julia_code = code