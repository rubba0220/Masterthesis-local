### Ancillary files for "One-loop QCD helicity amplitudes for $pp\to t\tb j$ to $O(\eps^2)$"

This directory contains a set of ancillary files with analytic expressions for the mass renormalised amplitudes in terms of master integrals and a set of routines to evaluate the master intgerals using differential equations and DiffExp (Hidding 2020, https://arxiv.org/abs/2006.05510).

The base directory contains two Mathematica script

* polechecks.wl - all independent helicity configurations for pp > ttj are evaluated using a high precision evaluation of the master integrals and cross checked against the universal pole structure.
* testpoint.wl - all independent helicity configurations for pp > ttj are evaluated up to O(eps^2)

both of these evaluations are performed at the level of the sub-amplitudes in the spin projection. The definitions are listed in the README.md file in the ttj_amplitudes directory.

Two directories contain computer readable files containing the results of the paper

ttj_amplitudes/README.md - details of the file format for the amplitudes
ttj_amplitudes/ttggg.tar.xz - tree and 1L sub-amplitudes for ttggg (must be unzipped)
ttj_amplitudes/ttqqg.tar.xz - tree and 1L sub-amplitudes for ttqqg (must be unzipped)

MIs_evaluation/README.md - details of the differential equations for the 4 independent pentagon topologies
MIs_evaluation/T*" - files for the differential equations
MIs_evaluation/T1T2T3T4_all_Permutations_eval_mt_not1_High_prec.m - high precision values for the MIs at a single point

boundary/README.md - details of the boundary constants for all the master integrals in all topologies
boundary/T*" - files for the boundary constants
boundary/pentagonT2 - files for the one-parametric boundary constants for pentagon for T2
