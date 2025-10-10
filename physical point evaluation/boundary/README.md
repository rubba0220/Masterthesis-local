### Ancillary files for "One-loop QCD helicity amplitudes for $pp\to t\tb j$ to $O(\eps^2)$"

This directory contains a set of files for the boundary constants for all the master integrals  in an analytic form as well as their numerical values. The numerical values files have been used to evaluate the master integrals using differential equations and DiffExp.

The base directory contains 5 sub-directories:

* T1 - analytic expressions of boundary constants in the files f's and the numerical value of these constants in bound_Num_T1.m 
* T2 - analytic expressions of boundary constants in the files f's and the numerical value of these constants in bound_Num_T2.m 
* T3 - analytic expressions of boundary constants in the files f's and the numerical value of these constants in bound_Num_T3.m 
* T4 - analytic expressions of boundary constants in the files f's and the numerical value of these constants in bound_Num_T4.m 
* pentagonT2-  instructions for evaluating slightly more complicated pentagon integrals can be evaluated numerically.

# Files inside Ti:
* **boundary_work_T1.wl**: file for the numerical evaluation of the boundary values with Ginac;
* **boundary_work_T2.wl**: file for the numerical evaluation of the boundary values with Ginac;
* **boundary_work_T3_Analytic.wl**: file for the numerical evaluation of the analytic boundary conditions with Ginac;
* **boundary_work_T3_diffexp.wl**: file for the numerical evaluation of the non-analytic boundary conditions with DiffExp (it has to be run after
                                   boundary_work_T3_Analytic.wl);
* **boundary_work_T4_Analytic.wl**: file for the numerical evaluation of the analytic boundary conditions with Ginac;
* **boundary_work_T4_diffexp.wl**: file for the numerical evaluation of the non-analytic boundary conditions with DiffExp (it has to be run after
                                   boundary_work_T3_Analytic.wl);

# Permutation.m file:
  It contains commands for the evaluation of the starting topologies.
