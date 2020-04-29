""" DisFit command line interface
:Author: Paul Lang <paul.lang@wolfson.ox.ac.uk>
:Date: 2020-04-26
:Copyright: 2020, Paul F Lang
:License: MIT
"""

from .core import DisFitProblem
# import DisFit
import cement
import os
import re
import sys


class BaseController(cement.Controller):
    """ Base controller for command line application """

    class Meta:
        label = 'base'
        description = "DisFit"
        help = "DisFit"
        arguments = [
            # (['-v', '--version'], dict(action='version', version=DisFit.__version__)),
        ]

    # @cement.ex(hide=False)
    def _default(self):
        raise SystemExit(self._parser.print_help())


class OptimizeController(cement.Controller):
    """ Optimize fitting problem """

    class Meta:
        label = 'optimize'
        description = 'Optimize a fitting problem'
        help = 'Optimize a fitting problem'
        stacked_on = 'base'
        stacked_type = 'nested'
        arguments = [
            (['sbml_file'], dict(type=str, help='sbml file')),
            (['data_file'], dict(type=str, help='data file in csv format')),
            (['-t', '--t_ratio'], dict(default=2, type=int,
                            help='ratio between experimental observation intervals and time discretization intervals')),
            (['-f', '--fold_change'], dict(default=2, type=float,
                            help='fold change window of parameter search range wrt sbml parameters')),
            (['-n', '--n_starts'], dict(default=1, type=int,
                            help='number of multistarts')),
            (['-o', '--out_dir'], dict(default='./DisFit_results', type=str,
                            help='output directory for julia_code, results and plot')),
            (['-p', '--plot_vars'], dict(default='[]', type=str,
                            help='list of variables to be plotted'))
        ]

    @cement.ex(hide=False)
    def _default(self):
        args = self.app.pargs

        try:
            print('--- Generating optimization problem ---')
            problem = DisFitProblem(args.sbml_file, args.data_file, t_ratio=args.t_ratio, 
                fold_change=args.fold_change, n_starts=args.n_starts)
        except Exception as error:
            raise SystemExit('Error occured: {}'.format(str(error)))

        try:
            print('\n--- Writing problem to julia file ---')
            if not os.path.isdir(args.out_dir):
                print('Creating {}'.format(args.out_dir))
                os.makedirs(args.out_dir)
            problem.write_jl_file(path=os.path.join(args.out_dir, 'julia_code.jl'))
        except Exception as error:
            print('Error occured: {}'.format(error), file=sys.stderr)
        
        try:
            print('\n--- Optimizing ---')
            problem.optimize()
        except Exception as error:
            raise SystemExit('Error occured: {}'.format(str(error)))

        try:
            print('\n--- Plotting results ---')
            variables = re.split(', |,', args.plot_vars.strip('[]'))
            problem.plot_results(path=os.path.join(args.out_dir, 'plot.pdf'),
                variables=variables)
        except Exception as error:
            print('Error occured: {}'.format(error), file=sys.stderr)

        try:
            print('\n--- Writing results to excel ---')
            problem.write_results(path=os.path.join(args.out_dir, 'results.xlsx'))
        except Exception as error:
            print('Error occured: {}'.format(error), file=sys.stderr)


class App(cement.App):
    """ Command line application """
    class Meta:
        label = 'DisFit'
        base_controller = 'base'
        handlers = [
            BaseController,
            OptimizeController
        ]


# if __name__ == '__main__':
def main():
    with App() as app:
        app.run()