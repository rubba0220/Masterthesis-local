(* ::Package:: *)

(* ::Section:: *)
(*Load*)


Quit[];


$dir = If[$FrontEnd=!=Null,NotebookDirectory[],DirectoryName[$InputFileName]];
SetDirectory[$dir]


Get[PathTo<>"DiffExp.m"];


(*Set $MaxExtraPrecision=500 allows avoiding warning messaged when generating boundaries in DiffExp*)


$MaxExtraPrecision=500


(*Loading commands to evaluate topology permutations in DiffExp*)


Get["../Permutation.m"];


(*Loading numerical evaluation of analytically computed boundary values*)


boundAnNum=Import["boundAnNum.m"];


(* ::Subsection:: *)
(*T1*)


(*Loading polynomial set for the analytic continuation of topology T1*)


AC=Get["../../MIs_Evaluation/T1/Analytic_Continuation_T1.m"];


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
MatrixDirectory->"../../MIs_Evaluation/T1/DEQs_T1/",
Verbosity->1,
UseMobius->False,
UsePade->False,
AccuracyGoal->16,
ExpansionOrder->80,
EpsilonOrder->4,WorkingPrecision->300,ChopPrecision->200
};

LoadConfiguration[OneLoopConfiguration];


(*Loading numerical boundary values*)


boundary=Get["../T1/bound_Num_T1.m"];


(*Point in which the boundaries are computed with respect to the ordering (1,2,3,4,5)*)


point=<|d12->-2,d15->-2,d23->-2,d34->-2,d45->-2,mt2->1|>;


(*Generate boundary values in DiffExp format*)


DiffExpBC=PrepareBoundaryConditions[boundary,point];


(*Generate the point "pointP6" in which the boundary are computed for permutation P6 of topology T1*)


inv0=InvariantsPermutations[[1,6]]/.mmp[p5,p3]->d35/.mmp[p5,p4]->d45;
inv1=inv0[[All,1]];
inv2=inv0[[All,2]]/.sijreplace/.point;
pointP6=Inner[Rule,inv1,inv2,List]


(*Transport the values of the boundary conditions from the point "point" to "pointP6"*)


ResultsP6=TransportTo[DiffExpBC,pointP6];


(*The boundary values for the MIs f15 and f16 in topology T4 are obtained from the values of the MIs f11 and f14 in permutation P6 of topology T1*)


b[16]=ResultsP6[[2]][[14]];
b[15]=ResultsP6[[2]][[11]];


(* ::Subsection:: *)
(*T2*)


(*Loading polynomial set for the analytic continuation of topology T2*)


AC=Get["../../MIs_Evaluation/T2/Analytic_Continuation_T2.m"];


OneLoopConfiguration={
Variables->{d12,d23,d34,d45,d15,mt2},
DeltaPrescriptions->AC,
MatrixDirectory->"../../MIs_Evaluation/T2/DEQs_T2/",
UseMobius->False,
UsePade->False,
EpsilonOrder->4,
AccuracyGoal->16,ExpansionOrder->50,ChopPrecision->200
};

LoadConfiguration[OneLoopConfiguration];


(*Loading numerical boundary values*)


boundary=Get["../T2/bound_Num_T2.m"];


(*Point in which the boundaries are computed with respect to the ordering (1,2,3,4,5)*)


point=<|d12->-2,d15->-2,d23->-2,d34->-2,d45->-2,mt2->1|>;


(*Generate boundary values in DiffExp format*)


DiffExpBC=PrepareBoundaryConditions[boundary,point];


(*Generate the point "pointP6" in which the boundary are computed for permutation P6 of topology T6*)


inv0=InvariantsPermutations[[2,6]]/.mmp[p5,p3]->d35/.mmp[p5,p4]->d45;
inv1=inv0[[All,1]];
inv2=inv0[[All,2]]/.sijreplace/.point;
pointP6=Inner[Rule,inv1,inv2,List]


(*In order to avoid singularities along the integration path we first transport the boundary values to a different point, "pointP6AUX"*)


pointP6AUX={d12->-2,d23->1,d34->4,d45->-2,d15->2,mt2->1};


ResultsP6AUX=TransportTo[DiffExpBC,pointP6AUX];


(*Transport the values of the boundary conditions to the point "pointP6"*)


ResultsP6=TransportTo[ResultsP6AUX,pointP6];


(*The boundary values for the MI f13,f12,f7 and f2 in topology T4 are obtained from the values of the MIs f5,f10,f11 and f4 in permutation P6 of topology T2*)


b[13]=-ResultsP6[[2]][[15]];
b[12]=ResultsP6[[2]][[10]];
b[7]=ResultsP6[[2]][[11]];
b[2]=-ResultsP6[[2]][[4]];


(* ::Subsection:: *)
(*T3*)


(*Loading polynomial set for the analytic continuation of topology T3*)


AC=Get["../../MIs_Evaluation/T3/Analytic_Continuation_T3.m"][[1]];


OneLoopConfiguration={
Variables->{d12,d15,d23,d34,d45,mt2},
DeltaPrescriptions->AC,
MatrixDirectory->"../../MIs_Evaluation/T3/DEQs_T3/",
Verbosity->1,
UseMobius->False,
UsePade->False,
AccuracyGoal->16,
ExpansionOrder->80,
EpsilonOrder->4,WorkingPrecision->300,ChopPrecision->200
};

LoadConfiguration[OneLoopConfiguration];


boundary=Get["../T3/bound_Num_T3.m"];


point=<|d12->-2,d15->-2,d23->-2,d34->-2,d45->-2,mt2->1|>;


DiffExpBC=PrepareBoundaryConditions[boundary,point];


inv0=InvariantsPermutations[[3,6]]/.mmp[p5,p3]->d35/.mmp[p5,p4]->d45;
inv1=inv0[[All,1]];
inv2=inv0[[All,2]]/.sijreplace/.point;
pointP6=Inner[Rule,inv1,inv2,List]


ResultsP6=TransportTo[DiffExpBC,pointP6];


b[4]=ResultsP6[[2]][[2]];


(* ::Subsection:: *)
(*Export*)


(*Create complete list of new boundary values for topology T4*)


listT4Prec={boundAnNum[[1]],b[2],boundAnNum[[2]],b[4],boundAnNum[[3]],boundAnNum[[4]],b[7],boundAnNum[[5]],boundAnNum[[6]],boundAnNum[[7]],boundAnNum[[8]],b[12],b[13],boundAnNum[[9]],b[15],b[16],boundAnNum[[10]],boundAnNum[[11]],boundAnNum[[12]]};


(*Loading boundary values evaluated with 100 digits accuracy*)


boundOG=Get["bound_Num_T4.m"];


(*Comparison between the new evaluation and the 100 digits one*)


Chop[listT4Prec-boundOG,10^(-4)]


Export["listT4Prec.m",listT4Prec];
