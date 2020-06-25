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
