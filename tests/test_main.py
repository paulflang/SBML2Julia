""" Tests of SBML2JuliaMP command line interface (SBML2JuliaMP.__main__)
:Author: Paul Lang <paul.lang@wolfson.ox.ac.uk>
:Date: 2020-04-26
:Copyright: 2020, Paul F Lang
:License: MIT
"""

import capturer
import mock
import os
import pkg_resources
import re
import SBML2JuliaMP
import shutil
import tempfile
import unittest
from SBML2JuliaMP import __main__

FIXTURES = pkg_resources.resource_filename('tests', 'fixtures')
YAML_PATH = os.path.join(FIXTURES, '0015_objectivePrior', '_0015_objectivePrior.yaml')
jl_file_gold = os.path.join(FIXTURES, 'jl_file_gold.jl')
with open(jl_file_gold, 'r') as f:
    JL_CODE_GOLD = f.read()
JL_CODE_GOLD = re.sub('/media/sf_DPhil_Project/Project07_Parameter Fitting/'
                      'df_software/SBML2JuliaMP/tests/fixtures', FIXTURES, JL_CODE_GOLD)


class CliVersionTestCase(unittest.TestCase):
    def test_version(self):
        with __main__.App(argv=['-v']) as app:
            with capturer.CaptureOutput(merged=False, relay=False) as captured:
                with self.assertRaises(SystemExit):
                    app.run()
                self.assertEqual(captured.stdout.get_text(), SBML2JuliaMP.__version__)
                self.assertEqual(captured.stderr.get_text(), '')

        with __main__.App(argv=['--version']) as app:
            with capturer.CaptureOutput(merged=False, relay=False) as captured:
                with self.assertRaises(SystemExit):
                    app.run()
                self.assertEqual(captured.stdout.get_text(), SBML2JuliaMP.__version__)
                self.assertEqual(captured.stderr.get_text(), '')


class CliTestCase(unittest.TestCase):
    def setUp(self):
        self.tempdir_1 = tempfile.mkdtemp()
        self.tempdir_2 = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tempdir_1)
        shutil.rmtree(self.tempdir_2)

    def test_cli(self):
        with mock.patch('sys.argv', ['SBML2JuliaMP', '--help']):
            with self.assertRaises(SystemExit) as context:
                __main__.main()
                self.assertRegex(context.Exception, 'usage: SBML2JuliaMP')

    def test_help(self):
        with self.assertRaises(SystemExit):
            with __main__.App(argv=[]) as app:
                app.run()

        with self.assertRaises(SystemExit):
            with __main__.App(argv=['--help']) as app:
                app.run()

    def test_optimize(self):
        with __main__.App(argv=['optimize', YAML_PATH, '-n', '1', '-i', 'False', '-o', '{}',
                                '-c', '{}', '-d', self.tempdir_1, '-p', '[obs_a, obs_b]']) as app:
            app.run()

            # test that the CLI produced the correct output
            with open(os.path.join(self.tempdir_1, 'julia_code.jl'), 'r') as f:
                jl_code = f.read()
            self.assertEqual(jl_code, JL_CODE_GOLD)
            self.assertTrue(os.path.exists(os.path.join(self.tempdir_1, 'plots', 'plot_c1.pdf')))
            self.assertTrue(os.path.exists(os.path.join(self.tempdir_1, 'parameters.tsv')))
            self.assertTrue(os.path.exists(os.path.join(self.tempdir_1, 'fval_chi2.tsv')))
            self.assertTrue(os.path.exists(os.path.join(self.tempdir_1, 'observables.tsv')))
            self.assertTrue(os.path.exists(os.path.join(self.tempdir_1, 'species.tsv')))
            self.assertTrue(os.path.exists(os.path.join(self.tempdir_1, 'julia_code.jl')))

        with __main__.App(argv=['optimize', YAML_PATH, '-t', '101', '-n', '2',
                                '-i', 'True', '-o', '{linear_solver: MA27}',
                                '-c', '{# Write global parameters: # Write global parameters1}',
                                '-d', self.tempdir_2, '-p', '[obs_a, obs_b]']) as app:
            app.run()
            self.assertTrue(os.path.exists(os.path.join(self.tempdir_2, 'plots', 'plot_c1.pdf')))
            self.assertTrue(os.path.exists(os.path.join(self.tempdir_2, 'parameters.tsv')))
            self.assertTrue(os.path.exists(os.path.join(self.tempdir_2, 'fval_chi2.tsv')))
            self.assertTrue(os.path.exists(os.path.join(self.tempdir_2, 'observables.tsv')))
            self.assertTrue(os.path.exists(os.path.join(self.tempdir_2, 'species.tsv')))
            self.assertTrue(os.path.exists(os.path.join(self.tempdir_2, 'julia_code.jl')))

        with __main__.App(argv=['optimize', 'a', '-t', '3', '-n', '2', '-i', 'True',
                                '-o', '{linear_solver: MA27}',
                                '-c', '{# Write global parameters: # Write global parameters1}',
                                '-d', self.tempdir_1, '-p', '[obs_Cb, obs_pCb]']) as app:
            with self.assertRaises(SystemExit):
                app.run()
