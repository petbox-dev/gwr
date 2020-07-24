Gaver-Wynn-Rho Algorithm
------------------------

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

Simple Example
--------------

.. code-block:: python

    >>> import math
    >>> from gwr_inversion import gwr
    >>> from mpmath import mp

    >>> def lap_ln_fn(s: float):
    ...     # log function
    ...     return -mp.log(s) / s - 0.577216 / s

    >>> gwr(lap_log_fn, time=5.0, M=32)
    mpf('1.6094375773356333')

    >>> math.log(5.0)
    1.6094379124341003


See the notebooks in ``test\`` for other use examples.

The method is described in: Valkó, P.P., and Abate J. 2002. Comparison of
Sequence Accelerators for the Gaver Method of Numerical Laplace Transform
Inversion. *Computers and Mathematics with Application* **48** (3): 629–636.
https://doi.org/10.1016/j.camwa.2002.10.017.

More information on multi-precision inversion can be found in: Valkó, P.P.and
Vajda, S. 2002. Inversion of Noise-Free Laplace Transforms: Towards a
Standardized Set of Test Problems. *Inverse Problems in Engineering* **10** (5):
467-483. https://doi.org/10.1080/10682760290004294.
