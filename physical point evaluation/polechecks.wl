(* ::Package:: *)

$dir = If[$FrontEnd=!=Null,NotebookDirectory[],DirectoryName[$InputFileName]];
SetDirectory[$dir];


(* include integral evaluations *)
integraleval = Get["MIs_evaluation/T1T2T3T4_all_Permutations_eval_mt_not1_High_prec.m"] /. minorm[__]:>1 /. \[Epsilon]->eps;


replacemi = integraleval[[2]];
dijvalues = integraleval[[1]];
dij2mt = {
d12 -> -(s34*t12*(-t45 + 2*t23*t51 - 2*t23*x5123 + 2*t12*t23*x5123 + 2*t45*x5123 + 2*t12*t51*x5123 - 2*t45*t51*x5123 - 2*t12*x5123^2 + 2*t12^2*x5123^2 - 2*t12*t45*x5123^2))/(2*t45),
d15 -> (s34*t51)/2,
d23 -> (s34*t23)/2,
d34 -> s34/2,
d45 -> (s34*t45)/2,
mt2 -> (s34*t12*(t23*t51 - t23*x5123 + t12*t23*x5123 + t45*x5123 + t12*t51*x5123 - t45*t51*x5123 - t12*x5123^2 + t12^2*x5123^2 - t12*t45*x5123^2))/t45
};
mt2dij = Solve[dij2mt /. Rule->Equal,{s34,t45,t51,t12,t23,x5123}] // Simplify;
momentumtwistorvalues = mt2dij[[1]] /. dijvalues // FullSimplify;


dijrules = {
d41 -> d14, d31->d13, d42->d24, d32->d23, d53->d35,
d13 -> d45-d12-d23-mt2,
d14 -> d23-d45-d15,
d24 -> d15-d23-d34,
d25 -> d34-d15-d12-mt2,
d35 -> d12-d34-d45+mt2
};


(* ::Section::Closed:: *)
(*ttggg*)


(* ::Text:: *)
(*In all ttggg one-loop amplitudes an additional overall factor of -8 is included which originates from the choice to polarisation vector normalisation. Since this normalisation differs from the one used at tree level we apply a correction of -1/8 to all expression for the poles obtained from Catani, Dittmaier and Trocsanyi.*)


(* ::Subsection::Closed:: *)
(*Nc T23451*)


(* ::Subsubsection::Closed:: *)
(*ds=2*)


poles = (-1/8)*A0[1, 2, 3, 4, 5]*(
- 3/eps^2
- (Log[-mt2/(2*d15)] + Log[-MuR2/(2*d15)])/(2*eps)
- (Log[-mt2/(2*d23)] + Log[-MuR2/(2*d23)])/(2*eps)
- Log[-MuR2/(2*d34)]/eps
- Log[-MuR2/(2*d45)]/eps
);


poleeval = poles /. A0[1,2,3,4,5]->Table[AA[i]["+++++"],{i,1,4}] 


Do[
Print[hel];
amp0L = Get["ttj_amplitudes/ttggg/tt3g_0L_T23451_"<>hel<>".m"];
amp1L = Get["ttj_amplitudes/ttggg/tt3g_1L_NcT23451_"<>hel<>"_dsm2p0.m"];

amp0Lvalues = amp0L[[1]] /. amp0L[[2]] /. amp0L[[3]] /. momentumtwistorvalues;
loopeval = amp1L[[1]] /. Rule[a_,b_]:>Rule[a,Normal@Series[b /. amp1L[[2]] /. amp1L[[3]] /. momentumtwistorvalues /. replacemi,{eps,0,-1}]];
poleeval = poles /. A0[1,2,3,4,5]->Table[AA[i][hel],{i,1,4}] /. amp0Lvalues /. dijvalues /. MuR2->1 // N[#,200]&;

Print[Collect[loopeval[[;;,2]]-poleeval,eps,Chop[#,10^(-130)]&]];
,{hel,{"+++++","++++-","+++-+"}}]


(* ::Subsubsection::Closed:: *)
(*ds-2*)


poles = (-1/8)*A0[1, 2, 3, 4, 5]*(
  1/4/eps
);


Do[
Print[hel];
amp0L = Get["ttj_amplitudes/ttggg/tt3g_0L_T23451_"<>hel<>".m"];
amp1L = Get["ttj_amplitudes/ttggg/tt3g_1L_NcT23451_"<>hel<>"_dsm2p1.m"];

amp0Lvalues = amp0L[[1]] /. amp0L[[2]] /. amp0L[[3]] /. momentumtwistorvalues;
loopeval = amp1L[[1]] /. Rule[a_,b_]:>Rule[a,Normal@Series[b /. amp1L[[2]] /. amp1L[[3]] /. momentumtwistorvalues /. replacemi,{eps,0,-1}]];
poleeval = poles /. A0[1,2,3,4,5]->Table[AA[i][hel],{i,1,4}] /. amp0Lvalues /. dijvalues /. MuR2->1 // N[#,200]&;

Print[Collect[loopeval[[;;,2]]-poleeval,eps,Chop[#,10^(-130)]&]];
,{hel,{"+++++","++++-","+++-+"}}]


(* ::Subsection::Closed:: *)
(*Nf T23451*)


(* ::Subsubsection::Closed:: *)
(*ds=2*)


Do[
Print[hel];

amp1L = Get["ttj_amplitudes/ttggg/tt3g_1L_NfT23451_"<>hel<>"_dsm2p0.m"];
loopeval = amp1L[[1]] /. Rule[a_,b_]:>Rule[a,Normal@Series[b /. amp1L[[2]] /. amp1L[[3]] /. momentumtwistorvalues /. replacemi,{eps,0,-1}]];

Print[Collect[loopeval[[;;,2]],eps,Chop[#,10^(-130)]&]];
,{hel,{"+++++","++++-","+++-+"}}]


(* ::Subsection::Closed:: *)
(*NhT23451*)


(* ::Subsubsection::Closed:: *)
(*ds=2*)


Do[
Print[hel];

amp1L = Get["ttj_amplitudes/ttggg/tt3g_1L_NhT23451_"<>hel<>"_dsm2p0.m"];
loopeval = amp1L[[1]] /. Rule[a_,b_]:>Rule[a,Normal@Series[b /. amp1L[[2]] /. amp1L[[3]] /. momentumtwistorvalues /. replacemi,{eps,0,-1}]];

Print[Collect[loopeval[[;;,2]],eps,Chop[#,10^(-130)]&]];
,{hel,{"+++++","++++-","+++-+"}}]


(* ::Subsection::Closed:: *)
(*1/Nc T23451*)


(* ::Subsubsection::Closed:: *)
(*ds=2*)


poles = (-1/8)*A0[1, 2, 3, 4, 5]*(
  d12/(d12+mt2)/\[Beta][1,2]/eps*Log[-(1-\[Beta][1,2])/(1+\[Beta][1,2])]
) /. \[Beta][1,2]->Sqrt[1-4*mt2/2/(d12+mt2)];


Do[
Print[hel];
amp0L = Get["ttj_amplitudes/ttggg/tt3g_0L_T23451_"<>hel<>".m"];
amp1L = Get["ttj_amplitudes/ttggg/tt3g_1L_Ncpm1T23451_"<>hel<>"_dsm2p0.m"];

amp0Lvalues = amp0L[[1]] /. amp0L[[2]] /. amp0L[[3]] /. momentumtwistorvalues;
loopeval = amp1L[[1]] /. Rule[a_,b_]:>Rule[a,Normal@Series[b /. amp1L[[2]] /. amp1L[[3]] /. momentumtwistorvalues /. replacemi,{eps,0,-1}]];
poleeval = poles /. A0[1,2,3,4,5]->Table[AA[i][hel],{i,1,4}] /. amp0Lvalues /. dijvalues /. MuR2->1 // N[#,200]&;

Print[Collect[loopeval[[;;,2]]-poleeval,eps,Chop[#,10^(-110)]&]];
,{hel,{"+++++","++++-","+++-+"}}]


(* ::Subsubsection::Closed:: *)
(*ds-2*)


poles = (-1/8)*A0[1, 2, 3, 4, 5]*(
  -1/4/eps
);


Do[
Print[hel];
amp0L = Get["ttj_amplitudes/ttggg/tt3g_0L_T23451_"<>hel<>".m"];
amp1L = Get["ttj_amplitudes/ttggg/tt3g_1L_Ncpm1T23451_"<>hel<>"_dsm2p1.m"];

amp0Lvalues = amp0L[[1]] /. amp0L[[2]] /. amp0L[[3]] /. momentumtwistorvalues;
loopeval = amp1L[[1]] /. Rule[a_,b_]:>Rule[a,Normal@Series[b /. amp1L[[2]] /. amp1L[[3]] /. momentumtwistorvalues /. replacemi,{eps,0,-1}]];
poleeval = poles /. A0[1,2,3,4,5]->Table[AA[i][hel],{i,1,4}] /. amp0Lvalues /. dijvalues /. MuR2->1 // N[#,200]&;

Print[Collect[loopeval[[;;,2]]-poleeval,eps,Chop[#,10^(-110)]&]];
,{hel,{"+++++","++++-","+++-+"}}]


(* ::Subsection::Closed:: *)
(*\[Delta]45 T231*)


(* ::Subsubsection::Closed:: *)
(*ds=2*)


poles = (
-(A0[1, 2, 3, 5, 4]*(
  (Log[-mt2/(2*d13)] + Log[-MuR2/(2*d13)])/(2*eps) - (Log[-mt2/(2*d15)] + Log[-MuR2/(2*d15)])/(2*eps)
  - Log[-MuR2/(2*d34)]/eps + Log[-MuR2/(2*d45)]/eps))
- A0[1, 2, 4, 5, 3]*(
  (Log[-mt2/(2*d23)] + Log[-MuR2/(2*d23)])/(2*eps) - (Log[-mt2/(2*d25)] + Log[-MuR2/(2*d25)])/(2*eps)
  - Log[-MuR2/(2*d34)]/eps + Log[-MuR2/(2*d45)]/eps)
- A0[1, 2, 3, 4, 5]*(
   (Log[-mt2/(2*d13)] + Log[-MuR2/(2*d13)])/(2*eps) - (Log[-mt2/(2*d14)] + Log[-MuR2/(2*d14)])/(2*eps)
   - Log[-MuR2/(2*d35)]/eps + Log[-MuR2/(2*d45)]/eps)
- A0[1, 2, 5, 4, 3]*(
  (Log[-mt2/(2*d23)] + Log[-MuR2/(2*d23)])/(2*eps) - (Log[-mt2/(2*d24)] + Log[-MuR2/(2*d24)])/(2*eps)
  - Log[-MuR2/(2*d35)]/eps + Log[-MuR2/(2*d45)]/eps)
)


(* The sub-leading colour poles contain permutation of the tree level amplitudes. As a result some
   manipulation is required to extract the pole structure of the sub-amplitudes in the spin projection.
   We provide the result of this manipulation to simplify the cross check *)


{subpole["+++++"],subpole["++++-"],subpole["+++-+"],treesub} = Get["ttj_amplitudes/ttggg/ttggg_d45T231_poles.m"];


Do[
Print[hel];
amp1L = Get["ttj_amplitudes/ttggg/tt3g_1L_d45T231_"<>hel<>"_dsm2p0.m"];

loopeval = amp1L[[1]] /. Rule[a_,b_]:>Rule[a,Normal@Series[b /. amp1L[[2]] /. amp1L[[3]] /. momentumtwistorvalues /. replacemi,{eps,0,-1}]];
poleeval = (-1/8)*subpole[hel][[;;,2]] /. treesub /. dijrules /. dijvalues /. momentumtwistorvalues /. MuR2->1 // N[#,150]&;

Print[Collect[loopeval[[;;,2]]-poleeval,eps,Chop[#,10^(-100)]&]];
,{hel,{"+++++","++++-","+++-+"}}]


(* ::Subsubsection::Closed:: *)
(*ds-2*)


Do[
Print[hel];
amp1L = Get["ttj_amplitudes/ttggg/tt3g_1L_d45T231_"<>hel<>"_dsm2p1.m"];
loopeval = amp1L[[1]] /. Rule[a_,b_]:>Rule[a,Normal@Series[b /. amp1L[[2]] /. amp1L[[3]] /. momentumtwistorvalues /. replacemi,{eps,0,-1}]];

Print[Collect[loopeval[[;;,2]],eps,Chop[#,10^(-105)]&]];
,{hel,{"+++++","++++-","+++-+"}}]


(* ::Subsection::Closed:: *)
(*\[Delta]12 Tr345*)


(* ::Subsubsection::Closed:: *)
(*ds=2*)


poles = (
 A0[1,2,4,5,3] ((Log[-(mt2/(2 d14))]+Log[-(MuR2/(2 d14))])/(2 eps)+(Log[-(mt2/(2 d23))]+Log[-(MuR2/(2 d23))])/(2 eps)-Log[-(MuR2/(2 d34))]/eps-(d12 Log[(-1+\[Beta][1,2])/(1+\[Beta][1,2])])/(eps (d12+mt2) \[Beta][1,2]))
+A0[1,2,3,4,5] ((Log[-(mt2/(2 d13))]+Log[-(MuR2/(2 d13))])/(2 eps)+(Log[-(mt2/(2 d25))]+Log[-(MuR2/(2 d25))])/(2 eps)-Log[-(MuR2/(2 d35))]/eps-(d12 Log[(-1+\[Beta][1,2])/(1+\[Beta][1,2])])/(eps (d12+mt2) \[Beta][1,2]))
+A0[1,2,5,3,4] ((Log[-(mt2/(2 d15))]+Log[-(MuR2/(2 d15))])/(2 eps)+(Log[-(mt2/(2 d24))]+Log[-(MuR2/(2 d24))])/(2 eps)-Log[-(MuR2/(2 d45))]/eps-(d12 Log[(-1+\[Beta][1,2])/(1+\[Beta][1,2])])/(eps (d12+mt2) \[Beta][1,2]))
)


(*poles = poles /. \[Beta][1,2]->Sqrt[1-4*mt2/2/(d12+mt2)];*)


(* The sub-leading colour poles contain permutation of the tree level amplitudes. As a result some
   manipulation is required to extract the pole structure of the sub-amplitudes in the spin projection.
   We provide the result of this manipulation to simplify the cross check *)


{subpole["+++++"],subpole["++++-"],subpole["+++-+"],treesub} = Get["ttj_amplitudes/ttggg/ttggg_d21TR345_poles.m"];


Do[
Print[hel];
amp1L = Get["ttj_amplitudes/ttggg/tt3g_1L_d21TR345_"<>hel<>"_dsm2p0.m"];

loopeval = amp1L[[1]] /. Rule[a_,b_]:>Rule[a,Normal@Series[b /. amp1L[[2]] /. amp1L[[3]] /. momentumtwistorvalues /. replacemi,{eps,0,-1}]];
poleeval = (-1/8)*subpole[hel][[;;,2]] /. treesub /. dijrules /. dijvalues /. momentumtwistorvalues /. MuR2->1 // N[#,150]&;

Print[Collect[loopeval[[;;,2]]-poleeval,eps,Chop[#,10^(-100)]&]];
,{hel,{"+++++","++++-","+++-+"}}]


(* ::Subsubsection::Closed:: *)
(*ds-2*)


Do[
Print[hel];
amp1L = Get["ttj_amplitudes/ttggg/tt3g_1L_d21TR345_"<>hel<>"_dsm2p1.m"];
loopeval = amp1L[[1]] /. Rule[a_,b_]:>Rule[a,Normal@Series[b /. amp1L[[2]] /. amp1L[[3]] /. momentumtwistorvalues /. replacemi,{eps,0,-1}]];

Print[Collect[loopeval[[;;,2]],eps,Chop[#,10^(-105)]&]];
,{hel,{"+++++","++++-","+++-+"}}]


(* ::Section::Closed:: *)
(*ttqqg*)


(* ::Text:: *)
(*In all ttqqg one-loop amplitudes an additional overall factor of -1/2 is included which originates from the choice to polarisation vector normalisation. Since this normalisation differs from the one used at tree level we apply a correction of -2 to all expression for the poles obtained from Catani, Dittmaier and Trocsanyi.*)


(* ::Subsection:: *)
(*Nc d14T253*)


(* ::Subsubsection::Closed:: *)
(*ds=2*)


poles = (2 A0["b",1,2,3,4,5])/eps+A0["b",1,2,3,4,5] (-(2/eps^2)-(Log[-(mt2/(2 d25))]+Log[-(MuR2/(2 d25))])/(2 eps)-Log[-(MuR2/(2 d35))]/eps-(Log[-(mt2/(2 d41))]+Log[-(MuR2/(2 d41))])/(2 eps))


Do[
Print[hel];
amp0L = Get["ttj_amplitudes/ttqqg/ttqqg_0L_d14T253_"<>hel<>".m"];
amp1L = Get["ttj_amplitudes/ttqqg/ttqqg_1L_d14T253_Ncp1_"<>hel<>"_dsm2p0.m"];

amp0Lvalues = amp0L[[1]] /. amp0L[[2]] /. amp0L[[3]] /. momentumtwistorvalues;
loopeval = amp1L[[1]] /. Rule[a_,b_]:>Rule[a,Normal@Series[b /. amp1L[[2]] /. amp1L[[3]] /. momentumtwistorvalues /. replacemi,{eps,0,-1}]];
poleeval = (-2)*poles /. A0["b",1,2,3,4,5]->Table[AA[i][hel],{i,1,4}] /. amp0Lvalues //. dijrules /. dijvalues /. MuR2->1 // N[#,200]&;

Print[Collect[loopeval[[;;,2]]-poleeval,eps,Chop[#,10^(-100)]&]];
,{hel,{"+++-+","++-++"}}]


(* ::Subsubsection::Closed:: *)
(*ds-2*)


poles = A0["b",1, 2, 3, 4, 5]*(
  1/3/eps
);


Do[
Print[hel];
amp0L = Get["ttj_amplitudes/ttqqg/ttqqg_0L_d14T253_"<>hel<>".m"];
amp1L = Get["ttj_amplitudes/ttqqg/ttqqg_1L_d14T253_Ncp1_"<>hel<>"_dsm2p1.m"];

amp0Lvalues = amp0L[[1]] /. amp0L[[2]] /. amp0L[[3]] /. momentumtwistorvalues;
loopeval = amp1L[[1]] /. Rule[a_,b_]:>Rule[a,Normal@Series[b /. amp1L[[2]] /. amp1L[[3]] /. momentumtwistorvalues /. replacemi,{eps,0,-1}]];
poleeval = (-2)*poles /. A0["b",1,2,3,4,5]->Table[AA[i][hel],{i,1,4}] /. amp0Lvalues //. dijrules /. dijvalues /. MuR2->1 // N[#,200]&;

Print[Collect[loopeval[[;;,2]]-poleeval,eps,Chop[#,10^(-105)]&]];
,{hel,{"+++-+","++-++"}}]


(* ::Subsection:: *)
(*1/Nc d14T253*)


(* ::Subsubsection::Closed:: *)
(*ds=2*)


poles = ((2 A0["b",1,2,3,4,5])/eps
+A0["c",1,2,3,4,5] ((Log[-(mt2/(2 d15))]+Log[-(MuR2/(2 d15))])/(2 eps)-(Log[-(mt2/(2 d25))]+Log[-(MuR2/(2 d25))])/(2 eps)-(Log[-(mt2/(2 d41))]+Log[-(MuR2/(2 d41))])/(2 eps)+(Log[-(mt2/(2 d42))]+Log[-(MuR2/(2 d42))])/(2 eps))
+A0["d",1,2,3,4,5] ((Log[-(mt2/(2 d31))]+Log[-(MuR2/(2 d31))])/(2 eps)-Log[-(MuR2/(2 d35))]/eps-(Log[-(mt2/(2 d41))]+Log[-(MuR2/(2 d41))])/(2 eps)+Log[-(MuR2/(2 d45))]/eps)
+A0["b",1,2,3,4,5] (1/eps^2-(Log[-(mt2/(2 d31))]+Log[-(MuR2/(2 d31))])/(2 eps)+(Log[-(mt2/(2 d32))]+Log[-(MuR2/(2 d32))])/(2 eps)+Log[-(MuR2/(2 d34))]/eps+(Log[-(mt2/(2 d41))]+Log[-(MuR2/(2 d41))])/(2 eps)-(Log[-(mt2/(2 d42))]+Log[-(MuR2/(2 d42))])/(2 eps)+(d12 Log[(-1+\[Beta][1,2])/(1+\[Beta][1,2])])/(eps (d12+mt2) \[Beta][1,2]))
)/. \[Beta][1,2]->Sqrt[1-4*mt2/2/(d12+mt2)];


Do[
Print[hel];
amp0La = Get["ttj_amplitudes/ttqqg/ttqqg_0L_d23T451_"<>hel<>".m"]/. AA[i_][h_]:>AA[i]["a",h];
amp0Lb = Get["ttj_amplitudes/ttqqg/ttqqg_0L_d14T253_"<>hel<>".m"]/. AA[i_][h_]:>AA[i]["b",h];
amp0Lc = Get["ttj_amplitudes/ttqqg/ttqqg_0L_d12T453_"<>hel<>".m"]/. AA[i_][h_]:>AA[i]["c",h];
amp0Ld = Get["ttj_amplitudes/ttqqg/ttqqg_0L_d34T251_"<>hel<>".m"]/. AA[i_][h_]:>AA[i]["d",h];

amp0Lavalues = amp0La[[1]] /. amp0La[[2]] /. amp0La[[3]] /. momentumtwistorvalues;
amp0Lbvalues = amp0Lb[[1]] /. amp0Lb[[2]] /. amp0Lb[[3]] /. momentumtwistorvalues;
amp0Lcvalues = amp0Lc[[1]] /. amp0Lc[[2]] /. amp0Lc[[3]] /. momentumtwistorvalues;
amp0Ldvalues = amp0Ld[[1]] /. amp0Ld[[2]] /. amp0Ld[[3]] /. momentumtwistorvalues;

amp0Lvalues = Join[amp0Lavalues,amp0Lbvalues,amp0Lcvalues,amp0Ldvalues];

amp1L = Get["ttj_amplitudes/ttqqg/ttqqg_1L_d14T253_Ncp-1_"<>hel<>"_dsm2p0.m"];
loopeval = amp1L[[1]] /. Rule[a_,b_]:>Rule[a,Normal@Series[b /. amp1L[[2]] /. amp1L[[3]] /. momentumtwistorvalues /. replacemi,{eps,0,-1}]];

poleeval = (-2)*poles /. A0[x_,1,2,3,4,5]:>Table[AA[i][x,hel],{i,1,4}] /. amp0Lvalues //. dijrules /. dijvalues /. MuR2->1 // N[#,200]&;

Print[Collect[loopeval[[;;,2]]-poleeval,eps,Chop[#,10^(-100)]&]];
,{hel,{"+++-+","++-++"}}]


(* ::Subsubsection::Closed:: *)
(*ds-2*)


poles = -1/2*A0["b",1,2,3,4,5]/eps


Do[
Print[hel];
amp0La = Get["ttj_amplitudes/ttqqg/ttqqg_0L_d23T451_"<>hel<>".m"]/. AA[i_][h_]:>AA[i]["a",h];
amp0Lb = Get["ttj_amplitudes/ttqqg/ttqqg_0L_d14T253_"<>hel<>".m"]/. AA[i_][h_]:>AA[i]["b",h];
amp0Lc = Get["ttj_amplitudes/ttqqg/ttqqg_0L_d12T453_"<>hel<>".m"]/. AA[i_][h_]:>AA[i]["c",h];
amp0Ld = Get["ttj_amplitudes/ttqqg/ttqqg_0L_d34T251_"<>hel<>".m"]/. AA[i_][h_]:>AA[i]["d",h];

amp0Lavalues = amp0La[[1]] /. amp0La[[2]] /. amp0La[[3]] /. momentumtwistorvalues;
amp0Lbvalues = amp0Lb[[1]] /. amp0Lb[[2]] /. amp0Lb[[3]] /. momentumtwistorvalues;
amp0Lcvalues = amp0Lc[[1]] /. amp0Lc[[2]] /. amp0Lc[[3]] /. momentumtwistorvalues;
amp0Ldvalues = amp0Ld[[1]] /. amp0Ld[[2]] /. amp0Ld[[3]] /. momentumtwistorvalues;

amp0Lvalues = Join[amp0Lavalues,amp0Lbvalues,amp0Lcvalues,amp0Ldvalues];

amp1L = Get["ttj_amplitudes/ttqqg/ttqqg_1L_d14T253_Ncp-1_"<>hel<>"_dsm2p1.m"];
loopeval = amp1L[[1]] /. Rule[a_,b_]:>Rule[a,Normal@Series[b /. amp1L[[2]] /. amp1L[[3]] /. momentumtwistorvalues /. replacemi,{eps,0,-1}]];

poleeval = (-2)*poles /. A0[x_,1,2,3,4,5]:>Table[AA[i][x,hel],{i,1,4}] /. amp0Lvalues //. dijrules /. dijvalues /. MuR2->1 // N[#,200]&;

Print[Collect[loopeval[[;;,2]]-poleeval,eps,Chop[#,10^(-100)]&]];
,{hel,{"+++-+","++-++"}}]


(* ::Subsection:: *)
(*Nc d23T451*)


(* ::Subsubsection::Closed:: *)
(*ds=2*)


poles = (2 A0["a",1,2,3,4,5])/eps+A0["a",1,2,3,4,5] (-(2/eps^2)-(Log[-(mt2/(2 d15))]+Log[-(MuR2/(2 d15))])/(2 eps)-(Log[-(mt2/(2 d32))]+Log[-(MuR2/(2 d32))])/(2 eps)-Log[-(MuR2/(2 d45))]/eps);


Do[
Print[hel];
amp0L = Get["ttj_amplitudes/ttqqg/ttqqg_0L_d23T451_"<>hel<>".m"];
amp1L = Get["ttj_amplitudes/ttqqg/ttqqg_1L_d23T451_Ncp1_"<>hel<>"_dsm2p0.m"];

amp0Lvalues = amp0L[[1]] /. amp0L[[2]] /. amp0L[[3]] /. momentumtwistorvalues;
loopeval = amp1L[[1]] /. Rule[a_,b_]:>Rule[a,Normal@Series[b /. amp1L[[2]] /. amp1L[[3]] /. momentumtwistorvalues /. replacemi,{eps,0,-1}]];
poleeval = (-2)*poles /. A0["a",1,2,3,4,5]->Table[AA[i][hel],{i,1,4}] /. amp0Lvalues //. dijrules /. dijvalues /. MuR2->1 // N[#,200]&;

Print[Collect[loopeval[[;;,2]]-poleeval,eps,Chop[#,10^(-130)]&]];
,{hel,{"+++-+","++-++"}}]


(* ::Subsubsection::Closed:: *)
(*ds-2*)


poles = 1/3*A0["a",1,2,3,4,5]/eps


Do[
Print[hel];
amp0L = Get["ttj_amplitudes/ttqqg/ttqqg_0L_d23T451_"<>hel<>".m"];
amp1L = Get["ttj_amplitudes/ttqqg/ttqqg_1L_d23T451_Ncp1_"<>hel<>"_dsm2p1.m"] /. mi["bubmT2",{{1},{2,3,4,5}}]->mi["bubmt2T1",{{1},{2,3,4,5}}];

amp0Lvalues = amp0L[[1]] /. amp0L[[2]] /. amp0L[[3]] /. momentumtwistorvalues;
loopeval = amp1L[[1]] /. Rule[a_,b_]:>Rule[a,Normal@Series[b /. amp1L[[2]] /. amp1L[[3]] /. momentumtwistorvalues /. replacemi,{eps,0,-1}]];
poleeval = (-2)*poles /. A0["a",1,2,3,4,5]->Table[AA[i][hel],{i,1,4}] /. amp0Lvalues //. dijrules /. dijvalues /. MuR2->1 // N[#,200]&;

Print[Collect[loopeval[[;;,2]]-poleeval,eps,Chop[#,10^(-130)]&]];
,{hel,{"+++-+","++-++"}}]


(* ::Subsection:: *)
(*1/Nc d23T451*)


(* ::Subsubsection::Closed:: *)
(*ds=2*)


poles = ((2 A0["a",1,2,3,4,5])/eps
+A0["c",1,2,3,4,5] (-((Log[-(mt2/(2 d15))]+Log[-(MuR2/(2 d15))])/(2 eps))+(Log[-(mt2/(2 d25))]+Log[-(MuR2/(2 d25))])/(2 eps)+(Log[-(mt2/(2 d31))]+Log[-(MuR2/(2 d31))])/(2 eps)-(Log[-(mt2/(2 d32))]+Log[-(MuR2/(2 d32))])/(2 eps))
+A0["d",1,2,3,4,5] (-((Log[-(mt2/(2 d32))]+Log[-(MuR2/(2 d32))])/(2 eps))+Log[-(MuR2/(2 d35))]/eps+(Log[-(mt2/(2 d42))]+Log[-(MuR2/(2 d42))])/(2 eps)-Log[-(MuR2/(2 d45))]/eps)
+A0["a",1,2,3,4,5] (1/eps^2-(Log[-(mt2/(2 d31))]+Log[-(MuR2/(2 d31))])/(2 eps)+(Log[-(mt2/(2 d32))]+Log[-(MuR2/(2 d32))])/(2 eps)+Log[-(MuR2/(2 d34))]/eps+(Log[-(mt2/(2 d41))]+Log[-(MuR2/(2 d41))])/(2 eps)-(Log[-(mt2/(2 d42))]+Log[-(MuR2/(2 d42))])/(2 eps)+(d12 Log[(-1+\[Beta][1,2])/(1+\[Beta][1,2])])/(eps (d12+mt2) \[Beta][1,2]))
)/. \[Beta][1,2]->Sqrt[1-4*mt2/2/(d12+mt2)];


Do[
Print[hel];
amp0La = Get["ttj_amplitudes/ttqqg/ttqqg_0L_d23T451_"<>hel<>".m"]/. AA[i_][h_]:>AA[i]["a",h];
amp0Lb = Get["ttj_amplitudes/ttqqg/ttqqg_0L_d14T253_"<>hel<>".m"]/. AA[i_][h_]:>AA[i]["b",h];
amp0Lc = Get["ttj_amplitudes/ttqqg/ttqqg_0L_d12T453_"<>hel<>".m"]/. AA[i_][h_]:>AA[i]["c",h];
amp0Ld = Get["ttj_amplitudes/ttqqg/ttqqg_0L_d34T251_"<>hel<>".m"]/. AA[i_][h_]:>AA[i]["d",h];

amp0Lavalues = amp0La[[1]] /. amp0La[[2]] /. amp0La[[3]] /. momentumtwistorvalues;
amp0Lbvalues = amp0Lb[[1]] /. amp0Lb[[2]] /. amp0Lb[[3]] /. momentumtwistorvalues;
amp0Lcvalues = amp0Lc[[1]] /. amp0Lc[[2]] /. amp0Lc[[3]] /. momentumtwistorvalues;
amp0Ldvalues = amp0Ld[[1]] /. amp0Ld[[2]] /. amp0Ld[[3]] /. momentumtwistorvalues;

amp0Lvalues = Join[amp0Lavalues,amp0Lbvalues,amp0Lcvalues,amp0Ldvalues];

amp1L = Get["ttj_amplitudes/ttqqg/ttqqg_1L_d23T451_Ncp-1_"<>hel<>"_dsm2p0.m"];
loopeval = amp1L[[1]] /. Rule[a_,b_]:>Rule[a,Normal@Series[b /. amp1L[[2]] /. amp1L[[3]] /. momentumtwistorvalues /. replacemi,{eps,0,-1}]];

poleeval = (-2)*poles /. A0[x_,1,2,3,4,5]:>Table[AA[i][x,hel],{i,1,4}] /. amp0Lvalues //. dijrules /. dijvalues /. MuR2->1 // N[#,200]&;

Print[Collect[loopeval[[;;,2]]-poleeval,eps,Chop[#,10^(-100)]&]];
,{hel,{"+++-+","++-++"}}]


(* ::Subsubsection::Closed:: *)
(*ds-2*)


poles = -1/2*A0["a",1,2,3,4,5]/eps


Do[
Print[hel];
amp0La = Get["ttj_amplitudes/ttqqg/ttqqg_0L_d23T451_"<>hel<>".m"]/. AA[i_][h_]:>AA[i]["a",h];
amp0Lb = Get["ttj_amplitudes/ttqqg/ttqqg_0L_d14T253_"<>hel<>".m"]/. AA[i_][h_]:>AA[i]["b",h];
amp0Lc = Get["ttj_amplitudes/ttqqg/ttqqg_0L_d12T453_"<>hel<>".m"]/. AA[i_][h_]:>AA[i]["c",h];
amp0Ld = Get["ttj_amplitudes/ttqqg/ttqqg_0L_d34T251_"<>hel<>".m"]/. AA[i_][h_]:>AA[i]["d",h];

amp0Lavalues = amp0La[[1]] /. amp0La[[2]] /. amp0La[[3]] /. momentumtwistorvalues;
amp0Lbvalues = amp0Lb[[1]] /. amp0Lb[[2]] /. amp0Lb[[3]] /. momentumtwistorvalues;
amp0Lcvalues = amp0Lc[[1]] /. amp0Lc[[2]] /. amp0Lc[[3]] /. momentumtwistorvalues;
amp0Ldvalues = amp0Ld[[1]] /. amp0Ld[[2]] /. amp0Ld[[3]] /. momentumtwistorvalues;

amp0Lvalues = Join[amp0Lavalues,amp0Lbvalues,amp0Lcvalues,amp0Ldvalues];

amp1L = Get["ttj_amplitudes/ttqqg/ttqqg_1L_d23T451_Ncp-1_"<>hel<>"_dsm2p1.m"];
loopeval = amp1L[[1]] /. Rule[a_,b_]:>Rule[a,Normal@Series[b /. amp1L[[2]] /. amp1L[[3]] /. momentumtwistorvalues /. replacemi,{eps,0,-1}]];

poleeval = (-2)*poles /. A0[x_,1,2,3,4,5]:>Table[AA[i][x,hel],{i,1,4}] /. amp0Lvalues //. dijrules /. dijvalues /. MuR2->1 // N[#,200]&;

Print[Collect[loopeval[[;;,2]]-poleeval,eps,Chop[#,10^(-100)]&]];
,{hel,{"+++-+","++-++"}}]


(* ::Subsection:: *)
(*d12T453*)


(* ::Subsubsection::Closed:: *)
(*ds=2*)


poles = ((2 A0["c",1,2,3,4,5])/eps
+A0["a",1,2,3,4,5] ((Log[-(mt2/(2 d25))]+Log[-(MuR2/(2 d25))])/(2 eps)+(Log[-(mt2/(2 d31))]+Log[-(MuR2/(2 d31))])/(2 eps)-Log[-(MuR2/(2 d35))]/eps-(d12 Log[(-1+\[Beta][1,2])/(1+\[Beta][1,2])])/(eps (d12+mt2) \[Beta][1,2]))
+A0["c",1,2,3,4,5] (-(2/eps^2)-Log[-(MuR2/(2 d35))]/eps-Log[-(MuR2/(2 d45))]/eps-(d12 Log[(-1+\[Beta][1,2])/(1+\[Beta][1,2])])/(eps (d12+mt2) \[Beta][1,2]))
+A0["b",1,2,3,4,5] ((Log[-(mt2/(2 d15))]+Log[-(MuR2/(2 d15))])/(2 eps)+(Log[-(mt2/(2 d42))]+Log[-(MuR2/(2 d42))])/(2 eps)-Log[-(MuR2/(2 d45))]/eps-(d12 Log[(-1+\[Beta][1,2])/(1+\[Beta][1,2])])/(eps (d12+mt2) \[Beta][1,2]))
) /. \[Beta][1,2]->Sqrt[1-4*mt2/2/(d12+mt2)];


Do[
Print[hel];
amp0La = Get["ttj_amplitudes/ttqqg/ttqqg_0L_d23T451_"<>hel<>".m"]/. AA[i_][h_]:>AA[i]["a",h];
amp0Lb = Get["ttj_amplitudes/ttqqg/ttqqg_0L_d14T253_"<>hel<>".m"]/. AA[i_][h_]:>AA[i]["b",h];
amp0Lc = Get["ttj_amplitudes/ttqqg/ttqqg_0L_d12T453_"<>hel<>".m"]/. AA[i_][h_]:>AA[i]["c",h];
amp0Ld = Get["ttj_amplitudes/ttqqg/ttqqg_0L_d34T251_"<>hel<>".m"]/. AA[i_][h_]:>AA[i]["d",h];

amp1L = Get["ttj_amplitudes/ttqqg/ttqqg_1L_d12T453_Ncp0_"<>hel<>"_dsm2p0.m"];

amp0Lavalues = amp0La[[1]] /. amp0La[[2]] /. amp0La[[3]] /. momentumtwistorvalues;
amp0Lbvalues = amp0Lb[[1]] /. amp0Lb[[2]] /. amp0Lb[[3]] /. momentumtwistorvalues;
amp0Lcvalues = amp0Lc[[1]] /. amp0Lc[[2]] /. amp0Lc[[3]] /. momentumtwistorvalues;
amp0Ldvalues = amp0Ld[[1]] /. amp0Ld[[2]] /. amp0Ld[[3]] /. momentumtwistorvalues;
amp0Lvalues = Join[amp0Lavalues,amp0Lbvalues,amp0Lcvalues,amp0Ldvalues];

loopeval = amp1L[[1]] /. Rule[a_,b_]:>Rule[a,Normal@Series[b /. amp1L[[2]] /. amp1L[[3]] /. momentumtwistorvalues /. replacemi,{eps,0,-1}]];
poleeval = (-2)*poles /. A0[x_,1,2,3,4,5]:>Table[AA[i][x,hel],{i,1,4}] /. amp0Lvalues //. dijrules /. dijvalues /. MuR2->1 // N[#,200]&;

Print[Collect[loopeval[[;;,2]]-poleeval,eps,Chop[#,10^(-100)]&]];
,{hel,{"+++-+","++-++"}}]


(* ::Subsubsection::Closed:: *)
(*ds-2*)


poles = 1/3*A0["c",1,2,3,4,5]/eps;


Do[
Print[hel];
amp0La = Get["ttj_amplitudes/ttqqg/ttqqg_0L_d23T451_"<>hel<>".m"]/. AA[i_][h_]:>AA[i]["a",h];
amp0Lb = Get["ttj_amplitudes/ttqqg/ttqqg_0L_d14T253_"<>hel<>".m"]/. AA[i_][h_]:>AA[i]["b",h];
amp0Lc = Get["ttj_amplitudes/ttqqg/ttqqg_0L_d12T453_"<>hel<>".m"]/. AA[i_][h_]:>AA[i]["c",h];
amp0Ld = Get["ttj_amplitudes/ttqqg/ttqqg_0L_d34T251_"<>hel<>".m"]/. AA[i_][h_]:>AA[i]["d",h];

amp1L = Get["ttj_amplitudes/ttqqg/ttqqg_1L_d12T453_Ncp0_"<>hel<>"_dsm2p1.m"];

amp0Lavalues = amp0La[[1]] /. amp0La[[2]] /. amp0La[[3]] /. momentumtwistorvalues;
amp0Lbvalues = amp0Lb[[1]] /. amp0Lb[[2]] /. amp0Lb[[3]] /. momentumtwistorvalues;
amp0Lcvalues = amp0Lc[[1]] /. amp0Lc[[2]] /. amp0Lc[[3]] /. momentumtwistorvalues;
amp0Ldvalues = amp0Ld[[1]] /. amp0Ld[[2]] /. amp0Ld[[3]] /. momentumtwistorvalues;
amp0Lvalues = Join[amp0Lavalues,amp0Lbvalues,amp0Lcvalues,amp0Ldvalues];

loopeval = amp1L[[1]] /. Rule[a_,b_]:>Rule[a,Normal@Series[b /. amp1L[[2]] /. amp1L[[3]] /. momentumtwistorvalues /. replacemi,{eps,0,-1}]];
poleeval = (-2)*poles /. A0[x_,1,2,3,4,5]:>Table[AA[i][x,hel],{i,1,4}] /. amp0Lvalues //. dijrules /. dijvalues /. MuR2->1 // N[#,200]&;

Print[Collect[loopeval[[;;,2]]-poleeval,eps,Chop[#,10^(-100)]&]];
,{hel,{"+++-+","++-++"}}]


(* ::Subsection:: *)
(*1/Nc^2 d12T453*)


(* ::Subsubsection::Closed:: *)
(*ds=2*)


poles = ((2 A0["c",1,2,3,4,5])/eps
+A0["c",1,2,3,4,5] (1/eps^2-(Log[-(mt2/(2 d31))]+Log[-(MuR2/(2 d31))])/(2 eps)+(Log[-(mt2/(2 d32))]+Log[-(MuR2/(2 d32))])/(2 eps)+Log[-(MuR2/(2 d34))]/eps+(Log[-(mt2/(2 d41))]+Log[-(MuR2/(2 d41))])/(2 eps)-(Log[-(mt2/(2 d42))]+Log[-(MuR2/(2 d42))])/(2 eps)+(d12 Log[(-1+\[Beta][1,2])/(1+\[Beta][1,2])])/(eps (d12+mt2) \[Beta][1,2]))
)/. \[Beta][1,2]->Sqrt[1-4*mt2/2/(d12+mt2)];


Do[
Print[hel];
amp0La = Get["ttj_amplitudes/ttqqg/ttqqg_0L_d23T451_"<>hel<>".m"]/. AA[i_][h_]:>AA[i]["a",h];
amp0Lb = Get["ttj_amplitudes/ttqqg/ttqqg_0L_d14T253_"<>hel<>".m"]/. AA[i_][h_]:>AA[i]["b",h];
amp0Lc = Get["ttj_amplitudes/ttqqg/ttqqg_0L_d12T453_"<>hel<>".m"]/. AA[i_][h_]:>AA[i]["c",h];
amp0Ld = Get["ttj_amplitudes/ttqqg/ttqqg_0L_d34T251_"<>hel<>".m"]/. AA[i_][h_]:>AA[i]["d",h];

amp1L = Get["ttj_amplitudes/ttqqg/ttqqg_1L_d12T453_Ncp-2_"<>hel<>"_dsm2p0.m"];

amp0Lavalues = amp0La[[1]] /. amp0La[[2]] /. amp0La[[3]] /. momentumtwistorvalues;
amp0Lbvalues = amp0Lb[[1]] /. amp0Lb[[2]] /. amp0Lb[[3]] /. momentumtwistorvalues;
amp0Lcvalues = amp0Lc[[1]] /. amp0Lc[[2]] /. amp0Lc[[3]] /. momentumtwistorvalues;
amp0Ldvalues = amp0Ld[[1]] /. amp0Ld[[2]] /. amp0Ld[[3]] /. momentumtwistorvalues;
amp0Lvalues = Join[amp0Lavalues,amp0Lbvalues,amp0Lcvalues,amp0Ldvalues];

loopeval = amp1L[[1]] /. Rule[a_,b_]:>Rule[a,Normal@Series[b /. amp1L[[2]] /. amp1L[[3]] /. momentumtwistorvalues /. replacemi,{eps,0,-1}]];
poleeval = (-2)*poles /. A0[x_,1,2,3,4,5]:>Table[AA[i][x,hel],{i,1,4}] /. amp0Lvalues //. dijrules /. dijvalues /. MuR2->1 // N[#,200]&;

Print[Collect[loopeval[[;;,2]]-poleeval,eps,Chop[#,10^(-100)]&]];
,{hel,{"+++-+","++-++"}}]


(* ::Subsubsection::Closed:: *)
(*ds-2*)


poles = 1/3*A0["c",1,2,3,4,5]/eps;


Do[
Print[hel];
amp0La = Get["ttj_amplitudes/ttqqg/ttqqg_0L_d23T451_"<>hel<>".m"]/. AA[i_][h_]:>AA[i]["a",h];
amp0Lb = Get["ttj_amplitudes/ttqqg/ttqqg_0L_d14T253_"<>hel<>".m"]/. AA[i_][h_]:>AA[i]["b",h];
amp0Lc = Get["ttj_amplitudes/ttqqg/ttqqg_0L_d12T453_"<>hel<>".m"]/. AA[i_][h_]:>AA[i]["c",h];
amp0Ld = Get["ttj_amplitudes/ttqqg/ttqqg_0L_d34T251_"<>hel<>".m"]/. AA[i_][h_]:>AA[i]["d",h];

amp1L = Get["ttj_amplitudes/ttqqg/ttqqg_1L_d12T453_Ncp0_"<>hel<>"_dsm2p1.m"];

amp0Lavalues = amp0La[[1]] /. amp0La[[2]] /. amp0La[[3]] /. momentumtwistorvalues;
amp0Lbvalues = amp0Lb[[1]] /. amp0Lb[[2]] /. amp0Lb[[3]] /. momentumtwistorvalues;
amp0Lcvalues = amp0Lc[[1]] /. amp0Lc[[2]] /. amp0Lc[[3]] /. momentumtwistorvalues;
amp0Ldvalues = amp0Ld[[1]] /. amp0Ld[[2]] /. amp0Ld[[3]] /. momentumtwistorvalues;
amp0Lvalues = Join[amp0Lavalues,amp0Lbvalues,amp0Lcvalues,amp0Ldvalues];

loopeval = amp1L[[1]] /. Rule[a_,b_]:>Rule[a,Normal@Series[b /. amp1L[[2]] /. amp1L[[3]] /. momentumtwistorvalues /. replacemi,{eps,0,-1}]];
poleeval = (-2)*poles /. A0[x_,1,2,3,4,5]:>Table[AA[i][x,hel],{i,1,4}] /. amp0Lvalues //. dijrules /. dijvalues /. MuR2->1 // N[#,200]&;

Print[Collect[loopeval[[;;,2]]-poleeval,eps,Chop[#,10^(-100)]&]];
,{hel,{"+++-+","++-++"}}]


(* ::Subsection:: *)
(*d34T251*)


(* ::Subsubsection::Closed:: *)
(*ds=2*)


poles = ((2 A0["d",1,2,3,4,5])/eps
+A0["d",1,2,3,4,5] (-(2/eps^2)-(Log[-(mt2/(2 d15))]+Log[-(MuR2/(2 d15))])/(2 eps)-(Log[-(mt2/(2 d25))]+Log[-(MuR2/(2 d25))])/(2 eps)-Log[-(MuR2/(2 d34))]/eps)
+A0["a",1,2,3,4,5] (-((Log[-(mt2/(2 d25))]+Log[-(MuR2/(2 d25))])/(2 eps))-Log[-(MuR2/(2 d34))]/eps+Log[-(MuR2/(2 d35))]/eps+(Log[-(mt2/(2 d42))]+Log[-(MuR2/(2 d42))])/(2 eps))
+A0["b",1,2,3,4,5] (-((Log[-(mt2/(2 d15))]+Log[-(MuR2/(2 d15))])/(2 eps))+(Log[-(mt2/(2 d31))]+Log[-(MuR2/(2 d31))])/(2 eps)-Log[-(MuR2/(2 d34))]/eps+Log[-(MuR2/(2 d45))]/eps)
) /. \[Beta][1,2]->Sqrt[1-4*mt2/2/(d12+mt2)];


Do[
Print[hel];
amp0La = Get["ttj_amplitudes/ttqqg/ttqqg_0L_d23T451_"<>hel<>".m"]/. AA[i_][h_]:>AA[i]["a",h];
amp0Lb = Get["ttj_amplitudes/ttqqg/ttqqg_0L_d14T253_"<>hel<>".m"]/. AA[i_][h_]:>AA[i]["b",h];
amp0Lc = Get["ttj_amplitudes/ttqqg/ttqqg_0L_d12T453_"<>hel<>".m"]/. AA[i_][h_]:>AA[i]["c",h];
amp0Ld = Get["ttj_amplitudes/ttqqg/ttqqg_0L_d34T251_"<>hel<>".m"]/. AA[i_][h_]:>AA[i]["d",h];

amp1L = Get["ttj_amplitudes/ttqqg/ttqqg_1L_d34T251_Ncp0_"<>hel<>"_dsm2p0.m"];

amp0Lavalues = amp0La[[1]] /. amp0La[[2]] /. amp0La[[3]] /. momentumtwistorvalues;
amp0Lbvalues = amp0Lb[[1]] /. amp0Lb[[2]] /. amp0Lb[[3]] /. momentumtwistorvalues;
amp0Lcvalues = amp0Lc[[1]] /. amp0Lc[[2]] /. amp0Lc[[3]] /. momentumtwistorvalues;
amp0Ldvalues = amp0Ld[[1]] /. amp0Ld[[2]] /. amp0Ld[[3]] /. momentumtwistorvalues;
amp0Lvalues = Join[amp0Lavalues,amp0Lbvalues,amp0Lcvalues,amp0Ldvalues];

loopeval = amp1L[[1]] /. Rule[a_,b_]:>Rule[a,Normal@Series[b /. amp1L[[2]] /. amp1L[[3]] /. momentumtwistorvalues /. replacemi,{eps,0,-1}]];
poleeval = (-2)*poles /. A0[x_,1,2,3,4,5]:>Table[AA[i][x,hel],{i,1,4}] /. amp0Lvalues //. dijrules /. dijvalues /. MuR2->1 // N[#,200]&;

Print[Collect[loopeval[[;;,2]]-poleeval,eps,Chop[#,10^(-100)]&]];
,{hel,{"+++-+","++-++"}}]


(* ::Subsubsection::Closed:: *)
(*ds-2*)


poles = 1/3*A0["d",1,2,3,4,5]/eps;


Do[
Print[hel];
amp0La = Get["ttj_amplitudes/ttqqg/ttqqg_0L_d23T451_"<>hel<>".m"]/. AA[i_][h_]:>AA[i]["a",h];
amp0Lb = Get["ttj_amplitudes/ttqqg/ttqqg_0L_d14T253_"<>hel<>".m"]/. AA[i_][h_]:>AA[i]["b",h];
amp0Lc = Get["ttj_amplitudes/ttqqg/ttqqg_0L_d12T453_"<>hel<>".m"]/. AA[i_][h_]:>AA[i]["c",h];
amp0Ld = Get["ttj_amplitudes/ttqqg/ttqqg_0L_d34T251_"<>hel<>".m"]/. AA[i_][h_]:>AA[i]["d",h];

amp1L = Get["ttj_amplitudes/ttqqg/ttqqg_1L_d34T251_Ncp0_"<>hel<>"_dsm2p1.m"];

amp0Lavalues = amp0La[[1]] /. amp0La[[2]] /. amp0La[[3]] /. momentumtwistorvalues;
amp0Lbvalues = amp0Lb[[1]] /. amp0Lb[[2]] /. amp0Lb[[3]] /. momentumtwistorvalues;
amp0Lcvalues = amp0Lc[[1]] /. amp0Lc[[2]] /. amp0Lc[[3]] /. momentumtwistorvalues;
amp0Ldvalues = amp0Ld[[1]] /. amp0Ld[[2]] /. amp0Ld[[3]] /. momentumtwistorvalues;
amp0Lvalues = Join[amp0Lavalues,amp0Lbvalues,amp0Lcvalues,amp0Ldvalues];

loopeval = amp1L[[1]] /. Rule[a_,b_]:>Rule[a,Normal@Series[b /. amp1L[[2]] /. amp1L[[3]] /. momentumtwistorvalues /. replacemi,{eps,0,-1}]];
poleeval = (-2)*poles /. A0[x_,1,2,3,4,5]:>Table[AA[i][x,hel],{i,1,4}] /. amp0Lvalues //. dijrules /. dijvalues /. MuR2->1 // N[#,200]&;

Print[Collect[loopeval[[;;,2]]-poleeval,eps,Chop[#,10^(-100)]&]];
,{hel,{"+++-+","++-++"}}]


(* ::Subsection:: *)
(*1/Nc^2 d34T251*)


(* ::Subsubsection::Closed:: *)
(*ds=2*)


poles = ((2 A0["d",1,2,3,4,5])/eps
+A0["d",1,2,3,4,5] (1/eps^2-(Log[-(mt2/(2 d31))]+Log[-(MuR2/(2 d31))])/(2 eps)+(Log[-(mt2/(2 d32))]+Log[-(MuR2/(2 d32))])/(2 eps)+Log[-(MuR2/(2 d34))]/eps+(Log[-(mt2/(2 d41))]+Log[-(MuR2/(2 d41))])/(2 eps)-(Log[-(mt2/(2 d42))]+Log[-(MuR2/(2 d42))])/(2 eps)+(d12 Log[(-1+\[Beta][1,2])/(1+\[Beta][1,2])])/(eps (d12+mt2) \[Beta][1,2]))
) /. \[Beta][1,2]->Sqrt[1-4*mt2/2/(d12+mt2)];


Do[
Print[hel];
amp0La = Get["ttj_amplitudes/ttqqg/ttqqg_0L_d23T451_"<>hel<>".m"]/. AA[i_][h_]:>AA[i]["a",h];
amp0Lb = Get["ttj_amplitudes/ttqqg/ttqqg_0L_d14T253_"<>hel<>".m"]/. AA[i_][h_]:>AA[i]["b",h];
amp0Lc = Get["ttj_amplitudes/ttqqg/ttqqg_0L_d12T453_"<>hel<>".m"]/. AA[i_][h_]:>AA[i]["c",h];
amp0Ld = Get["ttj_amplitudes/ttqqg/ttqqg_0L_d34T251_"<>hel<>".m"]/. AA[i_][h_]:>AA[i]["d",h];

amp1L = Get["ttj_amplitudes/ttqqg/ttqqg_1L_d34T251_Ncp-2_"<>hel<>"_dsm2p0.m"];

amp0Lavalues = amp0La[[1]] /. amp0La[[2]] /. amp0La[[3]] /. momentumtwistorvalues;
amp0Lbvalues = amp0Lb[[1]] /. amp0Lb[[2]] /. amp0Lb[[3]] /. momentumtwistorvalues;
amp0Lcvalues = amp0Lc[[1]] /. amp0Lc[[2]] /. amp0Lc[[3]] /. momentumtwistorvalues;
amp0Ldvalues = amp0Ld[[1]] /. amp0Ld[[2]] /. amp0Ld[[3]] /. momentumtwistorvalues;
amp0Lvalues = Join[amp0Lavalues,amp0Lbvalues,amp0Lcvalues,amp0Ldvalues];

loopeval = amp1L[[1]] /. Rule[a_,b_]:>Rule[a,Normal@Series[b /. amp1L[[2]] /. amp1L[[3]] /. momentumtwistorvalues /. replacemi,{eps,0,-1}]];
poleeval = (-2)*poles /. A0[x_,1,2,3,4,5]:>Table[AA[i][x,hel],{i,1,4}] /. amp0Lvalues //. dijrules /. dijvalues /. MuR2->1 // N[#,200]&;

Print[Collect[loopeval[[;;,2]]-poleeval,eps,Chop[#,10^(-100)]&]];
,{hel,{"+++-+","++-++"}}]


(* ::Subsubsection::Closed:: *)
(*ds-2*)


poles = -1/2*A0["d",1,2,3,4,5]/eps;


Do[
Print[hel];
amp0La = Get["ttj_amplitudes/ttqqg/ttqqg_0L_d23T451_"<>hel<>".m"]/. AA[i_][h_]:>AA[i]["a",h];
amp0Lb = Get["ttj_amplitudes/ttqqg/ttqqg_0L_d14T253_"<>hel<>".m"]/. AA[i_][h_]:>AA[i]["b",h];
amp0Lc = Get["ttj_amplitudes/ttqqg/ttqqg_0L_d12T453_"<>hel<>".m"]/. AA[i_][h_]:>AA[i]["c",h];
amp0Ld = Get["ttj_amplitudes/ttqqg/ttqqg_0L_d34T251_"<>hel<>".m"]/. AA[i_][h_]:>AA[i]["d",h];

amp1L = Get["ttj_amplitudes/ttqqg/ttqqg_1L_d34T251_Ncp-2_"<>hel<>"_dsm2p1.m"];

amp0Lavalues = amp0La[[1]] /. amp0La[[2]] /. amp0La[[3]] /. momentumtwistorvalues;
amp0Lbvalues = amp0Lb[[1]] /. amp0Lb[[2]] /. amp0Lb[[3]] /. momentumtwistorvalues;
amp0Lcvalues = amp0Lc[[1]] /. amp0Lc[[2]] /. amp0Lc[[3]] /. momentumtwistorvalues;
amp0Ldvalues = amp0Ld[[1]] /. amp0Ld[[2]] /. amp0Ld[[3]] /. momentumtwistorvalues;
amp0Lvalues = Join[amp0Lavalues,amp0Lbvalues,amp0Lcvalues,amp0Ldvalues];

loopeval = amp1L[[1]] /. Rule[a_,b_]:>Rule[a,Normal@Series[b /. amp1L[[2]] /. amp1L[[3]] /. momentumtwistorvalues /. replacemi,{eps,0,-1}]];
poleeval = (-2)*poles /. A0[x_,1,2,3,4,5]:>Table[AA[i][x,hel],{i,1,4}] /. amp0Lvalues //. dijrules /. dijvalues /. MuR2->1 // N[#,200]&;

Print[Collect[loopeval[[;;,2]]-poleeval,eps,Chop[#,10^(-100)]&]];
,{hel,{"+++-+","++-++"}}]


(* ::Subsection:: *)
(*Nf d14T253*)


(* ::Subsubsection::Closed:: *)
(*ds=2*)


poles = A0["b",1, 2, 3, 4, 5]*(
  -2/3/eps
);


Do[
Print[hel];
amp0L = Get["ttj_amplitudes/ttqqg/ttqqg_0L_d14T253_"<>hel<>".m"];
amp1L = Get["ttj_amplitudes/ttqqg/ttqqg_1L_Nfd14T253_Ncp0_"<>hel<>"_dsm2p0.m"];

amp0Lvalues = amp0L[[1]] /. amp0L[[2]] /. amp0L[[3]] /. momentumtwistorvalues;
loopeval = amp1L[[1]] /. Rule[a_,b_]:>Rule[a,Normal@Series[b /. amp1L[[2]] /. amp1L[[3]] /. momentumtwistorvalues /. replacemi,{eps,0,-1}]];
poleeval = (-2)*poles /. A0["b",1,2,3,4,5]->Table[AA[i][hel],{i,1,4}] /. amp0Lvalues //. dijrules /. dijvalues /. MuR2->1 // N[#,200]&;

Print[Collect[loopeval[[;;,2]]-poleeval,eps,Chop[#,10^(-105)]&]];
,{hel,{"+++-+","++-++"}}]


(* ::Subsection:: *)
(*Nf d23T451*)


(* ::Subsubsection::Closed:: *)
(*ds=2*)


poles = A0["a",1, 2, 3, 4, 5]*(
  -2/3/eps
);


Do[
Print[hel];
amp0L = Get["ttj_amplitudes/ttqqg/ttqqg_0L_d23T451_"<>hel<>".m"];
amp1L = Get["ttj_amplitudes/ttqqg/ttqqg_1L_Nfd23T451_Ncp0_"<>hel<>"_dsm2p0.m"];

amp0Lvalues = amp0L[[1]] /. amp0L[[2]] /. amp0L[[3]] /. momentumtwistorvalues;
loopeval = amp1L[[1]] /. Rule[a_,b_]:>Rule[a,Normal@Series[b /. amp1L[[2]] /. amp1L[[3]] /. momentumtwistorvalues /. replacemi,{eps,0,-1}]];
poleeval = (-2)*poles /. A0["a",1,2,3,4,5]->Table[AA[i][hel],{i,1,4}] /. amp0Lvalues //. dijrules /. dijvalues /. MuR2->1 // N[#,200]&;

Print[Collect[loopeval[[;;,2]]-poleeval,eps,Chop[#,10^(-105)]&]];
,{hel,{"+++-+","++-++"}}]


(* ::Subsection:: *)
(*Nf/Nc d12T453*)


(* ::Subsubsection::Closed:: *)
(*ds=2*)


poles = A0["c",1, 2, 3, 4, 5]*(
  -2/3/eps
);


Do[
Print[hel];
amp0L = Get["ttj_amplitudes/ttqqg/ttqqg_0L_d12T453_"<>hel<>".m"];
amp1L = Get["ttj_amplitudes/ttqqg/ttqqg_1L_Nfd12T453_Ncp-1_"<>hel<>"_dsm2p0.m"];

amp0Lvalues = amp0L[[1]] /. amp0L[[2]] /. amp0L[[3]] /. momentumtwistorvalues;
loopeval = amp1L[[1]] /. Rule[a_,b_]:>Rule[a,Normal@Series[b /. amp1L[[2]] /. amp1L[[3]] /. momentumtwistorvalues /. replacemi,{eps,0,-1}]];
poleeval = (-2)*poles /. A0["c",1,2,3,4,5]->Table[AA[i][hel],{i,1,4}] /. amp0Lvalues //. dijrules /. dijvalues /. MuR2->1 // N[#,200]&;

Print[Collect[loopeval[[;;,2]]-poleeval,eps,Chop[#,10^(-105)]&]];
,{hel,{"+++-+","++-++"}}]


(* ::Subsection:: *)
(*Nf/Nc d34T251*)


(* ::Subsubsection::Closed:: *)
(*ds=2*)


poles = A0["d",1, 2, 3, 4, 5]*(
  -2/3/eps
);


Do[
Print[hel];
amp0L = Get["ttj_amplitudes/ttqqg/ttqqg_0L_d34T251_"<>hel<>".m"];
amp1L = Get["ttj_amplitudes/ttqqg/ttqqg_1L_Nfd34T251_Ncp-1_"<>hel<>"_dsm2p0.m"];

amp0Lvalues = amp0L[[1]] /. amp0L[[2]] /. amp0L[[3]] /. momentumtwistorvalues;
loopeval = amp1L[[1]] /. Rule[a_,b_]:>Rule[a,Normal@Series[b /. amp1L[[2]] /. amp1L[[3]] /. momentumtwistorvalues /. replacemi,{eps,0,-1}]];
poleeval = (-2)*poles /. A0["d",1,2,3,4,5]->Table[AA[i][hel],{i,1,4}] /. amp0Lvalues //. dijrules /. dijvalues /. MuR2->1 // N[#,200]&;

Print[Collect[loopeval[[;;,2]]-poleeval,eps,Chop[#,10^(-105)]&]];
,{hel,{"+++-+","++-++"}}]


(* ::Subsection:: *)
(*Nh d14T253*)


(* ::Subsubsection::Closed:: *)
(*ds=2*)


poles = A0["b",1, 2, 3, 4, 5]*(
  -2/3/eps
);


Do[
Print[hel];
amp0L = Get["ttj_amplitudes/ttqqg/ttqqg_0L_d14T253_"<>hel<>".m"];
amp1L = Get["ttj_amplitudes/ttqqg/ttqqg_1L_Nhd14T253_Ncp0_"<>hel<>"_dsm2p0.m"];

amp0Lvalues = amp0L[[1]] /. amp0L[[2]] /. amp0L[[3]] /. momentumtwistorvalues;
loopeval = amp1L[[1]] /. Rule[a_,b_]:>Rule[a,Normal@Series[b /. amp1L[[2]] /. amp1L[[3]] /. momentumtwistorvalues /. replacemi,{eps,0,-1}]];
poleeval = (-2)*poles /. A0["b",1,2,3,4,5]->Table[AA[i][hel],{i,1,4}] /. amp0Lvalues //. dijrules /. dijvalues /. MuR2->1 // N[#,200]&;

Print[Collect[loopeval[[;;,2]]-poleeval,eps,Chop[#,10^(-105)]&]];
,{hel,{"+++-+","++-++"}}]


(* ::Subsection:: *)
(*Nh d23T451*)


(* ::Subsubsection::Closed:: *)
(*ds=2*)


poles = A0["a",1, 2, 3, 4, 5]*(
  -2/3/eps
);


Do[
Print[hel];
amp0L = Get["ttj_amplitudes/ttqqg/ttqqg_0L_d23T451_"<>hel<>".m"];
amp1L = Get["ttj_amplitudes/ttqqg/ttqqg_1L_Nhd23T451_Ncp0_"<>hel<>"_dsm2p0.m"];

amp0Lvalues = amp0L[[1]] /. amp0L[[2]] /. amp0L[[3]] /. momentumtwistorvalues;
loopeval = amp1L[[1]] /. Rule[a_,b_]:>Rule[a,Normal@Series[b /. amp1L[[2]] /. amp1L[[3]] /. momentumtwistorvalues /. replacemi,{eps,0,-1}]];
poleeval = (-2)*poles /. A0["a",1,2,3,4,5]->Table[AA[i][hel],{i,1,4}] /. amp0Lvalues //. dijrules /. dijvalues /. MuR2->1 // N[#,200]&;

Print[Collect[loopeval[[;;,2]]-poleeval,eps,Chop[#,10^(-105)]&]];
,{hel,{"+++-+","++-++"}}]


(* ::Subsection:: *)
(*Nh/Nc d12T453*)


(* ::Subsubsection::Closed:: *)
(*ds=2*)


poles = A0["c",1, 2, 3, 4, 5]*(
  -2/3/eps
);


Do[
Print[hel];
amp0L = Get["ttj_amplitudes/ttqqg/ttqqg_0L_d12T453_"<>hel<>".m"];
amp1L = Get["ttj_amplitudes/ttqqg/ttqqg_1L_Nhd12T453_Ncp-1_"<>hel<>"_dsm2p0.m"];

amp0Lvalues = amp0L[[1]] /. amp0L[[2]] /. amp0L[[3]] /. momentumtwistorvalues;
loopeval = amp1L[[1]] /. Rule[a_,b_]:>Rule[a,Normal@Series[b /. amp1L[[2]] /. amp1L[[3]] /. momentumtwistorvalues /. replacemi,{eps,0,-1}]];
poleeval = (-2)*poles /. A0["c",1,2,3,4,5]->Table[AA[i][hel],{i,1,4}] /. amp0Lvalues //. dijrules /. dijvalues /. MuR2->1 // N[#,200]&;

Print[Collect[loopeval[[;;,2]]-poleeval,eps,Chop[#,10^(-105)]&]];
,{hel,{"+++-+","++-++"}}]


(* ::Subsection:: *)
(*Nh/Nc d34T251*)


(* ::Subsubsection::Closed:: *)
(*ds=2*)


poles = A0["d",1, 2, 3, 4, 5]*(
  -2/3/eps
);


Do[
Print[hel];
amp0L = Get["ttj_amplitudes/ttqqg/ttqqg_0L_d34T251_"<>hel<>".m"];
amp1L = Get["ttj_amplitudes/ttqqg/ttqqg_1L_Nhd34T251_Ncp-1_"<>hel<>"_dsm2p0.m"];

amp0Lvalues = amp0L[[1]] /. amp0L[[2]] /. amp0L[[3]] /. momentumtwistorvalues;
loopeval = amp1L[[1]] /. Rule[a_,b_]:>Rule[a,Normal@Series[b /. amp1L[[2]] /. amp1L[[3]] /. momentumtwistorvalues /. replacemi,{eps,0,-1}]];
poleeval = (-2)*poles /. A0["d",1,2,3,4,5]->Table[AA[i][hel],{i,1,4}] /. amp0Lvalues //. dijrules /. dijvalues /. MuR2->1 // N[#,200]&;

Print[Collect[loopeval[[;;,2]]-poleeval,eps,Chop[#,10^(-105)]&]];
,{hel,{"+++-+","++-++"}}]
