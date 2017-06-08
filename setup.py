#!/usr/bin/env python
# -*- coding: utf-8 -*-

from setuptools import setup, find_packages
import versioneer

with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

requirements = [i.strip() for i in open('requirements.txt').readlines()]

setup(
    name='lcdblib',
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    description="A set of helper functions for bioinformatics analysis with snakemake.",
    long_description=readme + '\n\n' + history,
    author="Ryan Dale",
    author_email='dalerr@niddk.nih.gov',
    url='https://github.com/lcdb/lcdblib',
    packages=find_packages(),
    include_package_data=True,
    install_requires=requirements,
    license="MIT license",
    entry_points={
        'console_scripts':
        [
            'chrom_convert = lcdblib.utils.chrom_convert:main',
        ],
    },
    setup_requires=['pytest-runner'],
    tests_require=['pytest'],
)
