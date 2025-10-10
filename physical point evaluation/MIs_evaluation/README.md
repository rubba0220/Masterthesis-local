This folder contains the ancillary material for the numerical evaluation of the MIs with DiffExp.  
# Directories
## Ti:
Contains files for the numerical evaluation of MIs of topology Ti: 
### files description:
* **DEQs_Ti**: contains mathematica files for the differential equations matrices, as described in section 4 of the paper. 
               The file ddij_1.m is the differential equations matrix with respect to the kinematic invariant dij, 
               the file dmt2_1.m is the differential equations matrix with respect to the top-mass squared.
* **masters_Ti.m**: contains the list of MIs for topology Ti.
* **utdefsTi.m**: contains the defintion of the canonical MIs for topology Ti.
* **masters_Ti_norm.m**: contains additional normalization for the canonical MIs.
* **Analytic_Continuation_Ti.m**: contains polynomial lists for the analytic continuation within DiffExp.
* **DIFFEXP_Ti.wl**: mathematica .wl for the numerical evaluation of the MIs for topology Ti.
# Other files:
* **MIsReplacementRules.wl**: file to create a substitution list for the MIs into their numerical values.
* **T1T2T3T4_all_Permutations_eval_mt_not1_High_prec.m**: contains the results for the MIs with mt!=1,         
                                                          at the point specified inside the file.
