(* :Name: User`NumericalLaplaceInversion` *)

(* :Title: Numerical Inversion of Laplace Transform with Multiple Precision*)

(* :Author: Peter P. Valko and Joe Abate*)

(* :Summary:
This package provides only one function: GWR. The function calculates the value of the inverse of a Laplace transform at a specified time point. The Laplace transform should be provided as a function ready for multiple-precision evaluation. In other words, approximate numbers (with decimal point) or Mathematica functions starting with the letter 'N' are not allowed.
*)

(* :Context: User`NumericalLaplaceInversion` *)

(* :Package Version: 1.0 *)

(* :Copyright: Copyright 2002, Peter P. Valko and Joe Abate*)

(* :History: Originally written by Peter P. Valko, Dec 01, 2002. *)

(* :Keywords:
Laplace transform, Numerical inversion, Multiple-precision, Gaver functional, Wynn-Rho algorithm
*)

(* :Source:
Valk�, P. P. and Abate, J.:  Comparison of Sequence Accelerators for the Gaver Method of Numerical Laplace Transform Inversion,  Computers & Mathematics with Applications, (2002), accepted for publication (CAM 5307)
*)

(* :Warnings: None. *)

(* :Mathematica Version: 4.1 *)

(* :Limitations:
The Laplace transform function should be given in a form, that permits multiple precision evaluation.
*)

(* :Discussion:
The basic difference between this package and all other Mathematica packages aimed to do the same job is the systematic use of multiple precision arithmetics. The algorithm is also new.
*)

BeginPackage["User`NumericalLaplaceInversion`"]

GWR::usage =

"GWR[F, t, M, precin ] gives the inverse of the Laplace transform function named 'F' for a given time point 't'. The method involves the calculation of 'M' terms of the Gaver-functional. The obtained series is accelerated using the Wynn-Rho convergence acceleration scheme. The precision of internal calculations is set to 'precin'.
\n
\n
GWR[F, t, M ] does the same, but the precision of the internal calculations is selected automatically:  precin = 2.1 M).
\n
\n
GWR[F, t] uses M = 32 terms and precin = 67 as defaults. It should give reasonable results for many problems.
\n
\n
Important note: The Laplace transform should be degfined as a function of one argument. It can involve anything from a simple Mathematica expression to a sophisticated Module. Since the Laplace transform will be evaluated with non-standard (multiple) precision, approximate numbers (with decimal point) or Mathematica functions starting with the letter 'N' are not allowed in the function definition.
\n
\n
Example usage:
\n
\n
fun[s_]=(1/s) Exp[-Sqrt[ (s^2 + (37/100) s + 1)/(s^2 + s + Pi)]]
\n
t0=100
\n
GWR[fun,t0]
"
Unprotect[GWR];

Begin["User`NumericalLaplaceTransformInversion`Private`"]
GWR[F_, t_, M_:32, precin_:0] := Module[
    {M1, G0, Gm, Gp, best, expr, tau = Log[2]/t, Fi, broken, prec},
    If[precin <= 0, prec = 21 M/10, prec = precin];
    If[prec <= $MachinePrecision, prec = $MachinePrecision];
    broken = False;
    If[Precision[tau]<prec,tau=SetPrecision[tau,prec]];
    Do[Fi[i] = N[F[i tau], prec], {i, 1, 2 M}];
    M1 = M;
    Do[
      G0[n - 1] = tau(2n)!/(n!(n - 1)!)Sum[
            Binomial[n, i](-1)^i  Fi[n + i], {i, 0, n}];
      If[Not[NumberQ[G0[n - 1]]], M1 = n - 1; G0[n - 1] =.; Break[]];
      , {n, 1, M}];
    Do[Gm[n] = 0, {n, 0, M1}];
    best = G0[M1 - 1];

    Do[
      Do[
        expr = G0[n + 1] - G0[n];
        If[Or[Not[NumberQ[expr]], expr == 0], broken = True; Break[]];
        expr = Gm[n + 1] + (k + 1)/expr;
        Gp[n] = expr;
        If[OddQ[k],
          If[n == M1 - 2 - k, best = expr]
          ];
        , {n, M1 - 2 - k, 0, -1}];
      If[broken, Break[]];
      Do[Gm[n] = G0[n]; G0[n] = Gp[n], {n, 0, M1 - k}];
      , {k, 0, M1 - 2}];
    best
    ]

End[]  (* User`NumericalLaplaceTransformInversion` *)

Protect[GWR];

EndPackage[] (* User`NumericalLaplaceInversion` *)
