""" Tests of DisFit command line interface (DisFit.__main__)
:Author: Paul Lang <paul.lang@wolfson.ox.ac.uk>
:Date: 2020-04-26
:Copyright: 2020, Paul F Lang
:License: MIT
"""

import capturer
import DisFit
import filecmp
import mock
import os
import pkg_resources
import shutil
import tempfile
import unittest
from DisFit import __main__

FIXTURES = pkg_resources.resource_filename('tests', 'fixtures')
SBML_PATH = os.path.join(FIXTURES, 'G2M_copasi.xml')
DATA_PATH = os.path.join(FIXTURES, 'G2M_copasi.csv')
JL_FILE_GOLD = os.path.join(FIXTURES, 'jl_file_gold.jl')

class CliTestCase(unittest.TestCase):
    def setUp(self):
        self.tempdir_1 = tempfile.mkdtemp()
        self.tempdir_2 = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tempdir_1)
        shutil.rmtree(self.tempdir_2)

    def test_cli(self):
        with mock.patch('sys.argv', ['DisFit', '--help']):
            with self.assertRaises(SystemExit) as context:
                __main__.main()
                self.assertRegex(context.Exception, 'usage: DisFit')

    def test_help(self):
        with self.assertRaises(SystemExit):
            with __main__.App(argv=[]) as app:
                app.run()

        with self.assertRaises(SystemExit):
            with __main__.App(argv=['--help']) as app:
                app.run()

    '''
    def test_version(self):
        with __main__.App(argv=['-v']) as app:
            with capturer.CaptureOutput(merged=False, relay=False) as captured:
                with self.assertRaises(SystemExit):
                    app.run()
                self.assertEqual(captured.stdout.get_text(), DisFit.__version__)
                self.assertEqual(captured.stderr.get_text(), '')

        with __main__.App(argv=['--version']) as app:
            with capturer.CaptureOutput(merged=False, relay=False) as captured:
                with self.assertRaises(SystemExit):
                    app.run()
                self.assertEqual(captured.stdout.get_text(), DisFit.__version__)
                self.assertEqual(captured.stderr.get_text(), '')
                '''

    def test_optimize(self):
        with __main__.App(argv=['optimize', SBML_PATH, DATA_PATH,
            '-t', '2', '-f', '2', '-n', '1', '-o',
            self.tempdir_1, '-p', '[Cb, pCb]']) as app:
            # run app
            app.run()

            # test that the CLI produced the correct output
            self.assertTrue(filecmp.cmp(os.path.join(self.tempdir_1, 'julia_code.jl'), JL_FILE_GOLD))
            self.assertTrue(os.path.exists(os.path.join(self.tempdir_1, 'plot.pdf')))
            self.assertTrue(os.path.exists(os.path.join(self.tempdir_1, 'results.xlsx')))

        with __main__.App(argv=['optimize', 'a', DATA_PATH,
            '-t', '3', '-f', '3.2', '-n', '2', '-o', self.tempdir_1, '-p', '[Cb, pCb]']) as app:
            with self.assertRaises(SystemExit):
                app.run()

        with __main__.App(argv=['optimize', SBML_PATH, DATA_PATH,
            '-t', '3', '-f', '3.2', '-n', '2', '-o', self.tempdir_2, '-p', '[a, b]']) as app:
            app.run()
            self.assertTrue(os.path.exists(os.path.join(self.tempdir_2, 'results.xlsx')))