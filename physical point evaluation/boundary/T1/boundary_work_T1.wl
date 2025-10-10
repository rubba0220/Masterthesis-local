(* ::Package:: *)

(* ::Section:: *)
(*Load*)


$dir = If[$FrontEnd=!=Null,NotebookDirectory[],DirectoryName[$InputFileName]];
SetDirectory[$dir]


<<PolyLogTools`


(*Loading analytic expressions for boundary values of topology T1*)


Do[
f[ToExpression[i]]=Get["f"<>i<>".m"];
,{i,{"1","2","3","4","5","6","7","8","9","10","11","12","13","14","15"}}]


(* ::Section:: *)
(*Boundary Evaluation*)


(*The Prec variable sets the precision of the evaluation of the boundaries*)


Prec=16;


(*Evaluate numerically the boundaries with Ginac*)


Do[
b[i]=Table[Ginsh[f[i][[1]][[k]]/.f[i][[2]],{},PrecisionGoal->Prec],{k,5}];
,{i,1,15}]


listT1Prec={b[1],b[2],b[3],b[4],b[5],b[6],b[7],b[8],b[9],b[10],b[11],b[12],b[13],b[14],b[15]};


(*Loading boundary values evaluated with 100 digits accuracy*)


boundOG=Get["bound_Num_T1.m"];


(*Comparison between the new evaluation and the 100 digits one*)


Chop[listT1Prec-boundOG,10^(-15)]


Export["listT1Prec.m",listT1Prec];
