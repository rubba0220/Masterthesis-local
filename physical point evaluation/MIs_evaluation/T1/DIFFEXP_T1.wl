(* ::Package:: *)

(* ::Section:: *)
(*Load*)


$dir = If[$FrontEnd=!=Null,NotebookDirectory[],DirectoryName[$InputFileName]];
SetDirectory[$dir]


(*Load DiffExp*)


PathTo = "C:/Users/Laptop von Ruben/AppData/Roaming/Wolfram/Applications/DiffExp/";
Get[PathTo<>"DiffExp.m"];


(*Load the list of masters*)


masters=Get["masters_T1.m"];


(*Load the substitution list that writes explicitly masters[[k]] in terms of j[]: masters[[k]]//.mis*)


mis=Get["utdefsT1.m"];


(*Load the normalization for the masters*)


misnorm = Import["masters_T1_norm.m"];


(*Set $MaxExtraPrecision=500 allows avoiding warning messaged when generating boundaries in DiffExp*)


$MaxExtraPrecision=500


(*Loading numerical boundary values*)


boundary=Get["../../boundary/T1/bound_Num_T1.m"];


(*Loading polynomial set for the analytic continuation of topology T1*)


AC=Get["Analytic_Continuation_T1.m"];


(*Loading commands to evaluate topology permutations in DiffExp*)


Get["../../boundary/Permutation.m"];


(* ::Section:: *)
(*DIFFEXP*)


(*Load DiffExp configuration:
- Variables: set of variables for the differential equations;
- DeltaPrescriptions: polynomials for the analytic continuation;
- MatrixDirectory: path to the directory with the differential equations matrices;
- UseMobius: allows for using Mobius transformation in parametrization of the integration path. Setting it to false speed up the run;
- UsePade: allows for pad\[EGrave] approximands in the numerical evaluation: Setting it to false speed up the run;
- EpsilonOrder: set max order in the epsilon expansion for the solution;
- AccuracyGoal: set accuracy asked for the numerical evaluation;
- ExpansionOrder: set max order in the expansion of the differential equations around a given point;
- ChopPrecision: same as in mathematica*)


OneLoopConfiguration={
Variables->{d12,d23,d34,d45,d15,mt2},
DeltaPrescriptions->AC,
MatrixDirectory->"DEQs_T1/",
UseMobius->False,
UsePade->False,
EpsilonOrder->4,
AccuracyGoal->16,ExpansionOrder->50,ChopPrecision->200
};

LoadConfiguration[OneLoopConfiguration];


(*Point in which the boundaries are computed with respect to the ordering (1,2,3,4,5)*)


point=<|d12->-2,d15->-2,d23->-2,d34->-2,d45->-2,mt2->1|>;


(*Generate boundary values in DiffExp format*)


DiffExpBC=PrepareBoundaryConditions[boundary,point];


(*point where we want to evaluate the amplitude*)


point1=<|d12->0.1030127934867479,d23->-0.07518278423934056,d34->0.5,d45->-0.319871570505173,d15->0.31832575235540994,mt2->0.029821836099999898|>;


(*TransportTo: main command of DiffExp, it evolves the boundaries DiffExpBC to the value=point*)


ResultsP1=TransportTo[DiffExpBC,point1];


(*Code for creating a substitution list "masterrulesP1" for the MIs of permutation 1 into their numerical values.*)


perm1=permutations[[1]][[1]];


mastersnorm = Get["masters_T1.m"];
mastersnorm = Get["masters_T1_norm.m"] /. ResultsP1[[1]] /. Rule[a_,b_]:>Sqrt[b];
masterrulesP1 =Rule@@@Transpose[{masters /. i_Integer:>p[i] /. p[i_Integer]:>perm1[[i]],(({1,eps,eps^2,eps^3,eps^4} . #)&/@ResultsP1[[2]])/mastersnorm}] /. Rule[a_,b_]:>Rule[a,Collect[b,eps,Factor]];


Export["Numerics/resT1_p12345.m",masterrulesP1];


(* ::Section:: *)
(*P2*)


(*Numerical evaluation for permutation P2*)


(*point2 is the permuted kinematic point of point1*)


inv0=InvariantsPermutations[[1,2]]/.mmp[p4,p3]->d34/.mmp[p5,p4]->d45;
inv1=inv0[[All,1]];
inv2=inv0[[All,2]]/.sijreplace/.point1;
point2=Inner[Rule,inv1,inv2,List];


ResultsP2=TransportTo[DiffExpBC,point2];


perm2=permutations[[1]][[2]];


mastersnorm = Get["masters_T1.m"];
mastersnorm = Get["masters_T1_norm.m"] /. ResultsP2[[1]] /. Rule[a_,b_]:>Sqrt[b];
masterrulesP2 =Rule@@@Transpose[{masters /. i_Integer:>p[i] /. p[i_Integer]:>perm2[[i]],(({1,eps,eps^2,eps^3,eps^4} . #)&/@ResultsP2[[2]])/mastersnorm}] /. Rule[a_,b_]:>Rule[a,Collect[b,eps,Factor]];


Export["Numerics/resT1_p12435.m",masterrulesP2];


(* ::Section:: *)
(*P3*)


(*Numerical evaluation for permutation P3*)


(*point3 is the permuted kinematic point of point1*)


inv0=InvariantsPermutations[[1,3]]/.mmp[p5,p3]->d35/.mmp[p5,p4]->d45;
inv1=inv0[[All,1]];
inv2=inv0[[All,2]]/.sijreplace/.point1;
point3=Inner[Rule,inv1,inv2,List];


ResultsP3=TransportTo[DiffExpBC,point3];


perm3=permutations[[1]][[3]];


mastersnorm = Get["masters_T1.m"];
mastersnorm = Get["masters_T1_norm.m"] /. ResultsP3[[1]] /. Rule[a_,b_]:>Sqrt[b];
masterrulesP3 =Rule@@@Transpose[{masters /. i_Integer:>p[i] /. p[i_Integer]:>perm3[[i]],(({1,eps,eps^2,eps^3,eps^4} . #)&/@ResultsP3[[2]])/mastersnorm}] /. Rule[a_,b_]:>Rule[a,Collect[b,eps,Factor]];


Export["Numerics/resT1_p12453.m",masterrulesP3];


(* ::Section:: *)
(*P4*)


(*Numerical evaluation for permutation 4*)


(*point4 is the permuted kinematic point of point1*)


inv0=InvariantsPermutations[[1,4]]/.mmp[p4,p3]->d34/.mmp[p5,p4]->d45;
inv1=inv0[[All,1]];
inv2=inv0[[All,2]]/.sijreplace/.point1;
point4=Inner[Rule,inv1,inv2,List];


ResultsP4=TransportTo[DiffExpBC,point4];


perm4=permutations[[1]][[4]];


mastersnorm = Get["masters_T1.m"];
mastersnorm = Get["masters_T1_norm.m"] /. ResultsP4[[1]] /. Rule[a_,b_]:>Sqrt[b];
masterrulesP4 =Rule@@@Transpose[{masters /. i_Integer:>p[i] /. p[i_Integer]:>perm4[[i]],(({1,eps,eps^2,eps^3,eps^4} . #)&/@ResultsP4[[2]])/mastersnorm}] /. Rule[a_,b_]:>Rule[a,Collect[b,eps,Factor]];


Export["Numerics/resT1_p12543.m",masterrulesP4];


(* ::Section:: *)
(*P5*)


(*Numerical evaluation for permutation 5*)


(*point5 is the permuted kinematic point of point1*)


inv0=InvariantsPermutations[[1,5]]/.mmp[p4,p3]->d34/.mmp[p5,p4]->d45;
inv1=inv0[[All,1]];
inv2=inv0[[All,2]]/.sijreplace/.point1;
point5=Inner[Rule,inv1,inv2,List];


ResultsP5=TransportTo[DiffExpBC,point5];


perm5=permutations[[1]][[5]];


mastersnorm = Get["masters_T1.m"];
mastersnorm = Get["masters_T1_norm.m"] /. ResultsP5[[1]] /. Rule[a_,b_]:>Sqrt[b];
masterrulesP5 =Rule@@@Transpose[{masters /. i_Integer:>p[i] /. p[i_Integer]:>perm5[[i]],(({1,eps,eps^2,eps^3,eps^4} . #)&/@ResultsP5[[2]])/mastersnorm}] /. Rule[a_,b_]:>Rule[a,Collect[b,eps,Factor]];


Export["Numerics/resT1_p12354.m",masterrulesP5];


(* ::Section:: *)
(*P6*)


(*Numerical evaluation for permutation 6*)


(*point6 is the permuted kinematic point of point1*)


inv0=InvariantsPermutations[[1,6]]/.mmp[p5,p3]->d35/.mmp[p5,p4]->d45;
inv1=inv0[[All,1]];
inv2=inv0[[All,2]]/.sijreplace/.point1;
point6=Inner[Rule,inv1,inv2,List];


ResultsP6=TransportTo[DiffExpBC,point6];


perm6=permutations[[1]][[6]];


mastersnorm = Get["masters_T1.m"];
mastersnorm = Get["masters_T1_norm.m"] /. ResultsP6[[1]] /. Rule[a_,b_]:>Sqrt[b];
masterrulesP6 =Rule@@@Transpose[{masters /. i_Integer:>p[i] /. p[i_Integer]:>perm6[[i]],(({1,eps,eps^2,eps^3,eps^4} . #)&/@ResultsP6[[2]])/mastersnorm}] /. Rule[a_,b_]:>Rule[a,Collect[b,eps,Factor]];


Export["Numerics/resT1_p12534.m",masterrulesP6];
