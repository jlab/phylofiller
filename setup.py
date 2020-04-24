#!/usr/bin/env python

import re
import ast
from setuptools import find_packages, setup


# version parsing from __init__ pulled from Flask's setup.py
# https://github.com/mitsuhiko/flask/blob/master/setup.py
_version_re = re.compile(r'__version__\s+=\s+(.*)')

with open('phylofiller/__init__.py', 'rb') as f:
    hit = _version_re.search(f.read().decode('utf-8')).group(1)
    version = str(ast.literal_eval(hit))

classifiers = [
    'Development Status :: 1 - Pre-Alpha',
    'License :: OSI Approved :: BSD License',
    'Environment :: Console',
    'Topic :: Software Development :: Libraries :: Application Frameworks',
    'Topic :: Scientific/Engineering',
    'Topic :: Scientific/Engineering :: Bio-Informatics',
    'Programming Language :: Python',
    'Programming Language :: Python :: 3',
    'Programming Language :: Python :: 3.7',
    'Operating System :: Unix',
    'Operating System :: POSIX',
    'Operating System :: MacOS :: MacOS X',
    'Operating System :: Microsoft :: Windows']


description = 'PhyloPhiller: '
with open('README.md') as f:
    long_description = f.read()

keywords = '',

setup(name='phylophiller',
      version=version,
      license='BSD',
      description=description,
      long_description=long_description,
      keywords=keywords,
      classifiers=classifiers,
      author="Stefan Janssen",
      author_email="stefan.janssen@computational.bio.uni-giessen.de",
      maintainer="jlab development team",
      maintainer_email="stefan.janssen@computational.bio.uni-giessen.de",
      url='https://github.com/jlab/phylofiller',
      test_suite='nose.collector',
      packages=find_packages(),
      install_requires=[
      ],
      extras_require={'test': ["nose", "pep8", "flake8"],
                      'coverage': ["coverage"]})
