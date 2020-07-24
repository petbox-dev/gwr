This is a Python reproduction of the Mathematica package that provides the GWR
function, ``NumericalLaplaceInversion.m``.

https://library.wolfram.com/infocenter/MathSource/4738/

This package provides only one function: ``GWR``. The function calculates the
value of the inverse of a Laplace transform at a specified time value,
``Sequence`` of time values, or numpy array of time values.

The Laplace transform should be provided as a function that uses the ``mpmath``
library for a scalar value of the Laplace parameter.  The ``math`` library and
``numpy`` functions do not support multiprecision math and will return invalid
results if they are used.

The method is described in: Valkó, P.P.and Abate J.: Comparison of Sequence
Accelerators for the Gaver Method of Numerical Laplace Transform Inversion,
Computers & Mathematics with Applications, (2002), accepted for publication
(CAM 5307)

More information on multi-precision inversion can be found in: Valkó, P.P.and
Vajda, S : Inversion of noise-free Laplace transforms: Towards a standardized
set of test problems, Inverse Problems in Engineering, (2002) vol .10.No.5,
pp 467-483.
