"""
This is a Python reproduction of the Mathematica package that provides the GWR
function, ``NumericalLaplaceInversion.m``.

https://library.wolfram.com/infocenter/MathSource/4738/

This package provides only one function: ``gwr``. The function calculates the
value of the inverse of a Laplace transform at a specified time value,
``Sequence`` of time values, or numpy array of time values.

The Laplace transform should be provided as a function that uses the ``mpmath``
library for a scalar value of the Laplace parameter.  The ``math`` library and
``numpy`` functions do not support multiprecision math and will return invalid
results if they are used.

The method is described in: Valkó, P.P., and Abate J. 2002. Comparison of
Sequence Accelerators for the Gaver Method of Numerical Laplace Transform
Inversion. *Computers and Mathematics with Application* **48** (3): 629–636.
https://doi.org/10.1016/j.camwa.2002.10.017.

More information on multi-precision inversion can be found in: Valkó, P.P.and
Vajda, S. 2002. Inversion of Noise-Free Laplace Transforms: Towards a
Standardized Set of Test Problems. *Inverse Problems in Engineering* **10** (5):
467-483. https://doi.org/10.1080/10682760290004294.

Author
------
Peter P. Valko
Joe Abate

Python version by
-----------------
David S. Fulford

Notes
-----
Created on June 24, 2020
"""

import os
import sys
import re

try:
    from setuptools import setup  # type: ignore
except ImportError:
    from distutils.core import setup


__version__ = '1.0.1'


def get_long_description() -> str:
    # Fix display issues on PyPI caused by RST markup
    with open('README.rst', 'r') as f:
        return f.read()

if sys.argv[-1] == 'build':
    print(f'\nBuilding version {__version__}...\n')
    os.system('rm -r dist\\')  # clean out dist/
    os.system('python setup.py sdist bdist_wheel')
    sys.exit()


setup(
    name='gwr-inversion',
    version=__version__,
    description='GWR Algorithm Numerical Laplace Inversion',
    long_description=get_long_description(),
    long_description_content_type="text/x-rst",
    url='https://github.com/petbox-dev/gwr',
    author='David S. Fulford',
    author_email='petbox.dev@gmail.com',
    install_requires=['numpy>=1.17', 'mpmath>=1.1.0'],
    zip_safe=False,
    py_modules=['gwr_inversion'],
    package_data={
        '': ['py.typed']
    },
    python_requires='>=3.7',
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Science/Research',
        'Intended Audience :: Education',
        'Intended Audience :: Developers',
        'Natural Language :: English',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: Implementation :: CPython',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Mathematics',
        'Topic :: Software Development :: Libraries',
        'Typing :: Typed'
    ],
    keywords=[
        'laplace', 'inversion', 'transform', 'gwr'
    ],
)
