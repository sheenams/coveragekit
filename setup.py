#!/usr/bin/env python

from setuptools import setup, find_packages
from setuptools.command.test import test as TestCommand
import os
import sys

import coveragekit.version

here = os.path.abspath(os.path.dirname(__file__))

class PyTest(TestCommand):
    def finalize_options(self):
        TestCommand.finalize_options(self)
        self.test_args = []
        self.test_suite = True

    def run_tests(self):
        import pytest
        errcode = pytest.main(self.test_args)
        sys.exit(errcode)

setup(
    name='coveragekit',
    version=coveragekit.version.__version__,
    url='',
    license='DBAD',
    author='Christopher Hale',
    tests_require=['pytest'],
    install_requires=['pysam==0.8.4','pytest==2.7.2'],
    cmdclass={'test': PyTest},
    author_email='chris.joel.hale@gmail.com',
    description='NGS coverage analysis package.',
    scripts = ['coveragekit.py'],
    packages=['coveragekit'],
    include_package_data=True,
    platforms='Linux',
    test_suite='coveragekit.test.test_coveragekit',
    classifiers = [
        'Programming Language :: Python',
        'Natural Language :: English',
        'Environment :: Linux Command Line',
        'Intended Audience :: Bioinformaticists',
        'Operating System :: Linux'
        ],
    extras_require={
        'testing': ['pytest'],
    }
)