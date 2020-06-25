"""
This is a Python reproduction of the Mathematica package that provides the GWR
function, ``NumericalLaplaceInversion.m``.

https://library.wolfram.com/infocenter/MathSource/4738/

This package provides only one function: ``GWR``. The function calculates the
value of the inverse of a Laplace transform at a specified time point. The
Laplace transform should be provided as a function ready for multiple-precision
evaluation. In other words, approximate numbers (with decimal point) are not
allowed. ``Sympy`` can be used to rationalize a function before passing to GWR.

The method is described in: Valkó, P.P.and Abate J.: Comparison of Sequence
Accelerators for the Gaver Method of Numerical Laplace Transform Inversion,
Computers & Mathematics with Applications, (2002), accepted for publication
(CAM 5307)

More information on multi-precision inversion can be found in: Valkó, P.P.and
Vajda, S : Inversion of noise-free Laplace transforms: Towards a standardized
set of test problems, Inverse Problems in Engineering, (2002) vol .10.No.5,
pp 467-483.

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


__version__ = '1.0.0'


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
    name='gwr',
    version=__version__,
    description='Numerical Laplace Inversion',
    long_description=get_long_description(),
    long_description_content_type="text/x-rst",
    url='https://github.com/petbox-dev/gwr',
    author='David S. Fulford',
    author_email='petbox.dev@gmail.com',
    install_requires=['numpy>=1.17', 'mpmath>=1.1.0'],
    zip_safe=False,
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
        'laplace', 'inversion', 'transform'
    ],
)
