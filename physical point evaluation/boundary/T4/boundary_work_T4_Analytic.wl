(* ::Package:: *)

(* ::Section:: *)
(*Load*)


$dir = If[$FrontEnd=!=Null,NotebookDirectory[],DirectoryName[$InputFileName]];
SetDirectory[$dir]


<<PolyLogTools`


(*Loading analytic expressions for boundary values of topology T4*)


Do[
f[ToExpression[i]]=Get["f"<>i<>".m"];
,{i,{"1","3","5","6","8","10","11","14","17","18","19"}}]


(* ::Section:: *)
(*Boundary Evaluation*)


(* ::Subsection:: *)
(*Direct Integration*)


(*The Prec variable sets the precision of the evaluation of the boundaries*)


Prec=16;


Do[
b[i]=Table[Ginsh[f[i][[1]][[k]]/.f[i][[2]],{},PrecisionGoal->Prec],{k,5}];
,{i,{1,3,5,6,8,10,11,14,17,18,19}}]


b[9]={0,0,0,0,0};


boundAnNum={b[1],b[3],b[5],b[6],b[8],b[9],b[10],b[11],b[14],b[17],b[18],b[19]};


(*Store numerical evaluation of the analytically computed boundary values*)


Export["boundAnNum.m",boundAnNum];
