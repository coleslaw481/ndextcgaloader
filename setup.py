#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""The setup script."""
import os
import re
from setuptools import setup, find_packages


with open(os.path.join('ndextcgaloader', '__init__.py')) as ver_file:
    for line in ver_file:
        if line.startswith('__version__'):
            version=re.sub("'", "", line[line.index("'"):])

with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

requirements = ['ndex2>=3.1.0a1,<=4.0.0',
                'ndexutil>=0.1.0a3,<=1.0.0']

setup_requirements = [ ]

test_requirements = [ ]

setup(
    author="Chris Churas",
    author_email='churas.camera@gmail.com',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: BSD License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
    ],
    description="Loads TCGA data into NDEx",
    install_requires=requirements,
    license="BSD license",
    long_description=readme + '\n\n' + history,
    include_package_data=True,
    keywords='ndextcgaloader',
    name='ndextcgaloader',
    packages=find_packages(include=['ndextcgaloader']),
    scripts=[ 'ndextcgaloader/ndexloadtcga.py'],
    setup_requires=setup_requirements,
    test_suite='tests',
    tests_require=test_requirements,
    url='https://github.com/coleslaw481/ndextcgaloader',
    version=version,
    zip_safe=False,
)
