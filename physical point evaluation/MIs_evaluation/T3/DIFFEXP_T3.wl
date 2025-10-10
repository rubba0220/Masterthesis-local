(* ::Package:: *)

(* ::Section:: *)
(*Load*)


$dir = If[$FrontEnd=!=Null,NotebookDirectory[],DirectoryName[$InputFileName]];
SetDirectory[$dir]


(*Load DiffExp*)


PathTo = "C:/Users/Laptop von Ruben/AppData/Roaming/Wolfram/Applications/DiffExp/";
Get[PathTo<>"DiffExp.m"];


(*Load the list of masters*)


masters=Get["masters_T3.m"];


(*Load the substitution list that writes explicitly masters[[k]] in terms of j[]: masters[[k]]//.mis*)


mis=Get["utdefsT3.m"];


(*Load the normalization for the masters*)


misnorm = Import["masters_T3_norm.m"];


(*Set $MaxExtraPrecision=500 allows avoiding warning messaged when generating boundaries in DiffExp*)


$MaxExtraPrecision=500


(*Loading numerical boundary values*)


boundary=Get["../../boundary/T3/bound_Num_T3.m"];


(*Loading commands to evaluate topology permutations in DiffExp*)


Get["../../boundary/Permutation.m"];


(* ::Section:: *)
(*DIFFEXP Permutation P2,P4,P5,P6*)


AC=Get["Analytic_Continuation_T3.m"][[1]];


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
Variables->{d12,d15,d23,d34,d45,mt2},
DeltaPrescriptions->AC,
MatrixDirectory->"DEQs_T3/",
Verbosity->1,
UseMobius->False,
UsePade->False,
AccuracyGoal->16,
ExpansionOrder->80,
EpsilonOrder->4,WorkingPrecision->300,ChopPrecision->200
};

LoadConfiguration[OneLoopConfiguration];


(*Point in which the boundaries are computed with respect to the ordering (1,3,2,4,5)*)


point1=<|d12->0.1030127934867479,d23->-0.07518278423934056,d34->0.5,d45->-0.319871570505173,d15->0.31832575235540994,mt2->0.029821836099999898|>;


(*Generate boundary values in DiffExp format*)


DiffExpBC=PrepareBoundaryConditions[boundary,point];


(*point where we want to evaluate the amplitude*)


point5=<|d12->-(3/7),d23->-(7/11),d34->-(5/13),d45->-(7/5),d15->-(13/17),mt2->1|>;


(*TransportTo: main command of DiffExp, it evolves the boundaries DiffExpBC to the value=point*)


ResultsP5=TransportTo[DiffExpBC,point5];


(*Code for creating a substitution list "masterrulesP5" for the MIs of permutation 5 into their numerical values.*)


mastersnorm = Get["masters_T3.m"];
mastersnorm = Get["masters_T3_norm.m"] /. ResultsP5[[1]] /. Rule[a_,b_]:>Sqrt[b];
masterrulesP5 =Rule@@@Transpose[{masters /.Inner[Rule,{1,3,2,4,5},permutations[[3,5]],List],(({1,eps,eps^2,eps^3,eps^4} . #)&/@ResultsP5[[2]])/mastersnorm}] /. Rule[a_,b_]:>Rule[a,Collect[b,eps,Factor]];


Export["Numerics/resT3_p13245.m",masterrulesP5];


(* ::Section:: *)
(*P2*)


(*Numerical evaluation for permutation P2*)


(*point2 is the permuted kinematic point of point5*)


inv0=InvariantsPermutations[[3,2]]/.mmp[p4,p3]->d34/.mmp[p5,p4]->d45;
inv1=inv0[[All,1]];
inv2=inv0[[All,2]]/.sijreplace/.point5;
point2=Inner[Rule,inv1,inv2,List]


ResultsP2=TransportTo[DiffExpBC,point2];


mastersnorm = Get["masters_T3.m"];
mastersnorm = Get["masters_T3_norm.m"] /. ResultsP2[[1]] /. Rule[a_,b_]:>Sqrt[b];
masterrulesP2 =Rule@@@Transpose[{masters /.Inner[Rule,{1,3,2,4,5},permutations[[3,2]],List],(({1,eps,eps^2,eps^3,eps^4} . #)&/@ResultsP2[[2]])/mastersnorm}] /. Rule[a_,b_]:>Rule[a,Collect[b,eps,Factor]];


Export["Numerics/resT3_p15243.m",masterrulesP2];


(* ::Section:: *)
(*P3*)


(*Numerical evaluation for permutation 3*)


(*point3 is the permuted kinematic point of point5*)


InvariantsPermutations[[3,3]]/.mmp[p4,p3]->d34/.mmp[p5,p4]->d45;
%[[All,1]];
%%[[All,2]]/.sijreplace/.point5;
point3=Inner[Rule,%%,%,List]


ResultsP3=TransportTo[DiffExpBC,point3];


mastersnorm = Get["masters_T3.m"];
mastersnorm = Get["masters_T3_norm.m"] /. ResultsP3[[1]] /. Rule[a_,b_]:>Sqrt[b];
masterrulesP3 =Rule@@@Transpose[{masters /.Inner[Rule,{1,3,2,4,5},permutations[[3,3]],List],(({1,eps,eps^2,eps^3,eps^4} . #)&/@ResultsP3[[2]])/mastersnorm}] /. Rule[a_,b_]:>Rule[a,Collect[b,eps,Factor]];


Export[PathTo<>"/resT3_p14235.m",masterrulesP3];


(* ::Section:: *)
(*P4*)


(*Numerical evaluation for permutation P4*)


(*point4 is the permuted kinematic point of point5*)


inv0=InvariantsPermutations[[3,4]]/.mmp[p5,p3]->d35/.mmp[p5,p4]->d45;
inv1=inv0[[All,1]];
inv2=inv0[[All,2]]/.sijreplace/.point5;
point4=Inner[Rule,inv1,inv2,List]


ResultsP4=TransportTo[DiffExpBC,point4];


mastersnorm = Get["masters_T3.m"];
mastersnorm = Get["masters_T3_norm.m"] /. ResultsP4[[1]] /. Rule[a_,b_]:>Sqrt[b];
masterrulesP4 =Rule@@@Transpose[{masters /.Inner[Rule,{1,3,2,4,5},permutations[[3,4]],List],(({1,eps,eps^2,eps^3,eps^4} . #)&/@ResultsP4[[2]])/mastersnorm}] /. Rule[a_,b_]:>Rule[a,Collect[b,eps,Factor]];


Export["Numerics/resT3_p15234.m",masterrulesP4];


(* ::Section:: *)
(*P6*)


(*Numerical evaluation for permutation P6*)


(*point6 is the permuted kinematic point of point5*)


inv0=InvariantsPermutations[[3,6]]/.mmp[p5,p3]->d35/.mmp[p5,p4]->d45;
inv1=inv0[[All,1]];
inv2=inv0[[All,2]]/.sijreplace/.point5;
point6=Inner[Rule,inv1,inv2,List]


ResultsP6=TransportTo[DiffExpBC,point6];


mastersnorm = Get["masters_T3.m"];
mastersnorm = Get["masters_T3_norm.m"] /. ResultsP6[[1]] /. Rule[a_,b_]:>Sqrt[b];
masterrulesP6 =Rule@@@Transpose[{masters /.Inner[Rule,{1,3,2,4,5},permutations[[3,6]],List],(({1,eps,eps^2,eps^3,eps^4} . #)&/@ResultsP6[[2]])/mastersnorm}] /. Rule[a_,b_]:>Rule[a,Collect[b,eps,Factor]];


Export["Numerics/resT3_p13254.m",masterrulesP6];


(* ::Section:: *)
(*DIFFEXP Permutation P1,P3*)


ACP=Get["Analytic_Continuation_T3.m"][[2]];


OneLoopConfiguration={
Variables->{d12,d15,d23,d34,d45,mt2},
DeltaPrescriptions->ACP,
MatrixDirectory->"DEQs_T3/",
Verbosity->1,
UseMobius->False,
UsePade->False,
AccuracyGoal->16,
ExpansionOrder->80,
EpsilonOrder->4,WorkingPrecision->300,ChopPrecision->200
};

LoadConfiguration[OneLoopConfiguration];


(*Point in which the boundaries are computed with respect to the ordering (1,3,2,4,5)*)


point=<|d12->-2,d15->-2,d23->-2,d34->-2,d45->-2,mt2->1|>;


(*Generate boundary values in DiffExp format*)


DiffExpBC=PrepareBoundaryConditions[boundary,point];


(*point where we want to evaluate the amplitude*)


point5=<|d12->-(3/7),d23->-(7/11),d34->-(5/13),d45->-(7/5),d15->-(13/17),mt2->1|>;


(* ::Section:: *)
(*P1*)


(*Numerical evaluation for permutation P1*)


(*point1 is the permuted kinematic point of point5*)


inv0=InvariantsPermutations[[3,1]]/.mmp[p5,p3]->d35/.mmp[p5,p4]->d45;
inv1=inv0[[All,1]];
inv2=inv0[[All,2]]/.sijreplace/.point5;
point1=Inner[Rule,inv1,inv2,List]


ResultsP1=TransportTo[DiffExpBC,point1];


mastersnorm = Get["masters_T3.m"];
mastersnorm = Get["masters_T3_norm.m"] /. ResultsP1[[1]] /. Rule[a_,b_]:>Sqrt[b];
masterrulesP1 =Rule@@@Transpose[{masters /.Inner[Rule,{1,3,2,4,5},permutations[[3,1]],List],(({1,eps,eps^2,eps^3,eps^4} . #)&/@ResultsP1[[2]])/mastersnorm}] /. Rule[a_,b_]:>Rule[a,Collect[b,eps,Factor]];


Export["Numerics/resT3_p14253.m",masterrulesP1];


(* ::Section:: *)
(*P3*)


(*Numerical evaluation for permutation P3*)


(*point3 is the permuted kinematic point of point5*)


inv0=InvariantsPermutations[[3,3]]/.mmp[p4,p3]->d34/.mmp[p5,p4]->d45;
inv1=inv0[[All,1]];
inv2=inv0[[All,2]]/.sijreplace/.point5;
point3=Inner[Rule,inv1,inv2,List]


ResultsP3=TransportTo[DiffExpBC,point3];


mastersnorm = Get["masters_T3.m"];
mastersnorm = Get["masters_T3_norm.m"] /. ResultsP3[[1]] /. Rule[a_,b_]:>Sqrt[b];
masterrulesP3 =Rule@@@Transpose[{masters /.Inner[Rule,{1,3,2,4,5},permutations[[3,3]],List],(({1,eps,eps^2,eps^3,eps^4} . #)&/@ResultsP3[[2]])/mastersnorm}] /. Rule[a_,b_]:>Rule[a,Collect[b,eps,Factor]];


Export["Numerics/resT3_p14235.m",masterrulesP3];
