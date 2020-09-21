""" DisFit command line interface
:Author: Paul Lang <paul.lang@wolfson.ox.ac.uk>
:Date: 2020-04-26
:Copyright: 2020, Paul F Lang
:License: MIT
"""

from .core import DisFitProblem
import DisFit
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
            (['-v', '--version'], dict(action='version', version=DisFit.__version__)),
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
            (['petab_yaml'], dict(type=str, help='PEtab yaml problem specification')),
            (['-t', '--t_steps'], dict(default='', type=str,
                            help='number of time-discretization steps')),
            (['-n', '--n_starts'], dict(default=1, type=int,
                            help='number of multistarts')),
            (['-i', '--infer_ic_from_sbml'], dict(default=False, type=bool,
                            help='if missing initial conditions shall be infered from SBML model')),
            (['-o', '--optimizer_options'], dict(default='{}', type=str,
                            help='optimization solver options')),
            (['-c', '--custom_code_dict'], dict(default='{}', type=str,
                            help='dict with replaced code as keys and replacement code as values')),
            (['-d', '--out_dir'], dict(default='./DisFit_results', type=str,
                            help='output directory for julia_code, results and plot')),
            (['-p', '--plot_obs'], dict(default='[]', type=str,
                            help='list of observables to be plotted'))
        ]

    @cement.ex(hide=False)
    def _default(self):
        args = self.app.pargs
        try:
            print('--- Generating optimization problem ---')
            if args.t_steps == '':
                t_steps = None
            else:
                t_steps = int(args.t_steps)
            items = re.split(',', args.optimizer_options.strip('{}'))
            if items == ['']:
                optimizer_options = {}
            else:
                optimizer_options = {re.split(':', item)[0].strip(): re.split(':', item)[1].strip() for item in items}
            items = re.split(',', args.custom_code_dict.strip('{}'))
            if items == ['']:
                custom_code_dict = {}
            else:
                custom_code_dict = {re.split(':', item)[0]: re.split(':', item)[1].lstrip() for item in items}
            problem = DisFitProblem(args.petab_yaml, t_steps=t_steps,
                n_starts=args.n_starts, infer_ic_from_sbml=args.infer_ic_from_sbml,
                optimizer_options=optimizer_options,
                custom_code_dict=custom_code_dict)
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
            # if args.plot_obs != 'None':
            if not os.path.isdir(os.path.join(args.out_dir, 'plots')):
                print('Creating {}'.format(os.path.join(args.out_dir, 'plots')))
                os.makedirs(os.path.join(args.out_dir, 'plots'))
            observables = re.split(',', args.plot_obs.strip('[]'))
            observables = [o.strip() for o in observables]
            if observables == ['']:
                observables = []
            for c in [c for c, c_ind in problem._condition2index.items() if c_ind+1 in problem._j_to_parameters[0]]:
                problem.plot_results(c, path=os.path.join(args.out_dir, 'plots', 'plot_'+c+'.pdf'),
                    observables=observables)
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