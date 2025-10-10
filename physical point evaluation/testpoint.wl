(* ::Package:: *)

$dir = If[$FrontEnd=!=Null,NotebookDirectory[],DirectoryName[$InputFileName]];
SetDirectory[$dir];


(* include integral evaluations *)
integraleval = Get["MIs_evaluation/T1T2T3T4_all_Permutations.m"] /. minorm[__]:>1 /. \[Epsilon]->eps;


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


(* ::Section:: *)
(*ttggg*)


(* ::Text:: *)
(*In all ttggg one-loop amplitudes an additional overall factor of -8 is included which originates from the choice to polarisation vector normalisation.*)


dspart["NcT23451"] = {"dsm2p0","dsm2p1"};
dspart["NfT23451"] = {"dsm2p0"};
dspart["NhT23451"] = {"dsm2p0"};
dspart["Ncpm1T23451"] = {"dsm2p0","dsm2p1"};
dspart["d45T231"] = {"dsm2p0","dsm2p1"};
dspart["d21TR345"] = {"dsm2p0","dsm2p1"};


Do[
Print[" === ",part," ",hel," ",dspow," === "];
amp1L = Get["ttj_amplitudes/ttggg/tt3g_1L_"<>part<>"_"<>hel<>"_"<>dspow<>".m"];
loopeval = -8*amp1L[[1]] /. Rule[a_,b_]:>Rule[a,Normal@Series[b /. amp1L[[2]] /. amp1L[[3]] /. momentumtwistorvalues /. replacemi,{eps,0,2}]];

Do[
Print[i," ",Collect[loopeval[[i,2]],eps,Chop[N[Chop[#],5]]&] /. eps->\[Epsilon]];
,{i,1,4}];

,{part,{"NcT23451","NfT23451","NhT23451","Ncpm1T23451","d45T231","d21TR345"}}
,{hel,{"+++++","++++-","+++-+"}}
,{dspow,dspart[part]}]


(* ::Section:: *)
(*ttqqg*)


(* ::Text:: *)
(*In all ttqqg one-loop amplitudes an additional overall factor of -1/2 is included which originates from the choice to polarisation vector normalisation.*)


dspowlist["d23T451"] = {"dsm2p0","dsm2p1"};
dspowlist["d14T253"] = {"dsm2p0","dsm2p1"};
dspowlist["d12T453"] = {"dsm2p0","dsm2p1"};
dspowlist["d34T251"] = {"dsm2p0","dsm2p1"};

dspowlist["Nfd23T451"] = {"dsm2p0"};
dspowlist["Nfd14T253"] = {"dsm2p0"};
dspowlist["Nfd12T453"] = {"dsm2p0"};
dspowlist["Nfd34T251"] = {"dsm2p0"};

dspowlist["Nhd23T451"] = {"dsm2p0"};
dspowlist["Nhd14T253"] = {"dsm2p0"};
dspowlist["Nhd12T453"] = {"dsm2p0"};
dspowlist["Nhd34T251"] = {"dsm2p0"};

ncpowlist["d23T451"] = {"Ncp1","Ncp-1"};
ncpowlist["d14T253"] = {"Ncp1","Ncp-1"};
ncpowlist["d12T453"] = {"Ncp0","Ncp-2"};
ncpowlist["d34T251"] = {"Ncp0","Ncp-2"};

ncpowlist["Nfd23T451"] = {"Ncp0"};
ncpowlist["Nfd14T253"] = {"Ncp0"};
ncpowlist["Nfd12T453"] = {"Ncp-1"};
ncpowlist["Nfd34T251"] = {"Ncp-1"};

ncpowlist["Nhd23T451"] = {"Ncp0"};
ncpowlist["Nhd14T253"] = {"Ncp0"};
ncpowlist["Nhd12T453"] = {"Ncp-1"};
ncpowlist["Nhd34T251"] = {"Ncp-1"};


Do[
Print[" === ",part," ",ncpow," ",hel," ",dspow," === "];
amp1L = Get["ttj_amplitudes/ttqqg/ttqqg_1L_"<>part<>"_"<>ncpow<>"_"<>hel<>"_"<>dspow<>".m"];
loopeval = (-1/2)*amp1L[[1]] /. Rule[a_,b_]:>Rule[a,Normal@Series[b /. amp1L[[2]] /. amp1L[[3]] /. momentumtwistorvalues /. replacemi,{eps,0,2}]];

Do[
Print[i," ",Collect[loopeval[[i,2]],eps,Chop[N[Chop[#],5]]&] /. eps->\[Epsilon]];
,{i,1,4}];

,{part,{"d23T451","d14T253","d12T453","d34T251","Nfd23T451","Nfd14T253","Nfd12T453","Nfd34T251","Nhd23T451","Nhd14T253","Nhd12T453","Nhd34T251"}}
,{ncpow,ncpowlist[part]}
,{hel,{"++-++","+++-+"}}
,{dspow,dspowlist[part]}]
