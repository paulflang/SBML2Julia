""" Tests of calc command line interface (calc.__main__)

:Author: Jonathan Karr <karr@mssm.edu>
:Date: 2019-01-31
:Copyright: 2019, Karr Lab
:License: MIT
"""

from calc import __main__
import calc
import capturer
import mock
import os
import shutil
import tempfile
import unittest


class CliTestCase(unittest.TestCase):
    def setUp(self): # What does this even do. Never seems to be used in bpforms
        self.tempdir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tempdir)

    def test_cli(self):
        with mock.patch('sys.argv', ['calc', '--help']):
            with self.assertRaises(SystemExit) as context:
                __main__.main()
                self.assertRegex(context.Exception, 'usage: calc')

    def test_help(self):
        with self.assertRaises(SystemExit):
            with __main__.App(argv=[]) as app:
                app.run()

        with self.assertRaises(SystemExit):
            with __main__.App(argv=['--help']) as app:
                app.run()

    # def test_version(self):
    #     with __main__.App(argv=['-v']) as app:
    #         with capturer.CaptureOutput(merged=False, relay=False) as captured:
    #             with self.assertRaises(SystemExit):
    #                 app.run()
    #             self.assertEqual(captured.stdout.get_text(), calc.__version__)
    #             self.assertEqual(captured.stderr.get_text(), '')

    #     with __main__.App(argv=['--version']) as app:
    #         with capturer.CaptureOutput(merged=False, relay=False) as captured:
    #             with self.assertRaises(SystemExit):
    #                 app.run()
    #             self.assertEqual(captured.stdout.get_text(), calc.__version__)
    #             self.assertEqual(captured.stderr.get_text(), '')

    def test_isprime(self):
        with capturer.CaptureOutput(merged=False, relay=False) as captured:
            with __main__.App(argv=['isprime', '5']) as app:
                # run app
                app.run()

                # test that the CLI produced the correct output
                self.assertEqual(captured.stdout.get_text(), 'True') # note that you need .get_text(), cause whatever is printed to stdout is text (as apposed to integer, logical, etc.)
                self.assertEqual(captured.stderr.get_text(), '')

    def test_isprime(self):
        with capturer.CaptureOutput(merged=False, relay=False) as captured:
            with __main__.App(argv=['isprime', '4']) as app:
                # run app
                app.run()

                # test that the CLI produced the correct output
                self.assertEqual(captured.stdout.get_text(), 'False')
                self.assertEqual(captured.stderr.get_text(), '')

        # with self.assertRaisesRegex(SystemExit, '^Form is invalid'):
        #     with __main__.App(argv=['isprime', 'five']) as app:
        #         # run app
        #         app.run()

if __name__ == "__main__":
    unittest.main()