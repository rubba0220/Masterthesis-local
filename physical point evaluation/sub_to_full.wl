(* ::Package:: *)

ClearAll;

$dir = If[$FrontEnd=!=Null,NotebookDirectory[],DirectoryName[$InputFileName]];
SetDirectory[$dir];

Get[ToFileName[{"C:/Users/Laptop von Ruben/AppData/Roaming/Wolfram/Applications/Spinors-1.0"},"init.windows" ]];
<<Spinors`

Get["C:/Users/Laptop von Ruben/AppData/Roaming/Wolfram/Applications/SU3basics.m"];


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


mdot[a_, b_] := a[[1]]b[[1]]-a[[2]]b[[2]]-a[[3]]b[[3]]-a[[4]]b[[4]];

p1v = {451.160381942158, -254.955478242156, -127.963213069012, -303.886449763003};
p2v = {181.67424764459, 46.8336036978071, -3.17432367453475, 31.3086791659089};
p3v = -{500., 0, 0, 500.};
p4v = -{500., 0, 0, -500.};
p5v = {367.165370413252, 208.121874544348, 131.137536743547, 272.577770597094};
mtv = Sqrt[mdot[p5v,p5v]];

n1v = {100,100,0,0};
n2v = {100,0,100,0};
f1v = p1v - mdot[p1v, p1v]/(2*mdot[p1v, n1v]) * n1v;
f2v = p2v - mdot[p2v, p2v]/(2*mdot[p2v, n2v]) * n2v;

s34v = mdot[p3v+p4v, p3v+p4v];

Clear[f1,f2,n1,n2,p3,p4,p5];
DeclareSpinor[f1,f2,n1,n2,p3,p4,p5];
spinv = {f1v,f2v,n1v,n2v,p3v,p4v,p5v}/Sqrt[s34v];
spin = {f1,f2,n1,n2,p3,p4,p5};
Do[
DeclareSpinorMomentum[spin[[i]], spinv[[i]]];
,{i,{1,2,3,4,5,6,7}}];

spA[a_, b_] := Spaa[a, b];                  
spB[a_, b_] := Spbb[a, b];               
spBB[a_, b_, c_, d_] := Spbb[a, b, c, d];    
spAA[a_, b_, c_, d_] := Spaa[a, b, c, d];

phasepm = mtv/Sqrt[s34v]*spA[p3,p4]/(spA[p3,p5]^2)//N;
phasemp =  mtv/Sqrt[s34v]*spA[p3,p4]/(spA[p4,p5]^2)//N;

theta1 = spA[n1,n2]*1*spA[f1,n1]/spA[f2,n2]//N;
theta2 = spA[n1,p3]*spA[n2,p4]*spB[p3,p4]/spA[f1,n1]/spA[f2,n2]//N;
theta3 = spA[n1,p3]*spA[n2,p3]*spBB[p3,p4,p5,p3]/1/spA[f1,n1]/spA[f2,n2]//N;
theta4 = spA[n1,p4]*spA[n2,p4]*spBB[p4,p5,p3,p4]/1/spA[f1,n1]/spA[f2,n2]//N;

theta = {theta1, theta2, theta3, theta4};
Print[theta . Transpose[{theta}]]
Print[theta]


Clear[CT];
CT[1][i1_, i2bar_, i3_, i4bar_, a5_] := KroneckerDelta[i1, i4bar]*T[a5][[i3, i2bar]];
CT[2][i1_, i2bar_, i3_, i4bar_, a5_] := KroneckerDelta[i3, i2bar]*T[a5][[i1, i4bar]];
CT[3][i1_, i2bar_, i3_, i4bar_, a5_] := -1/Nc * KroneckerDelta[i1, i2bar]*T[a5][[i3, i4bar]];
CT[4][i1_, i2bar_, i3_, i4bar_, a5_] := -1/Nc * KroneckerDelta[i3, i4bar]*T[a5][[i1, i2bar]];

Gij[i_,j_] := Sum[Conjugate[CT[i][i1, i2bar, i3, i4bar, a5]]*CT[j][i1, i2bar, i3, i4bar, a5]
,{i1,1,3},{i2bar,1,3},{i3,1,3},{i4bar,1,3},{a5,1,8}];



cScheme["tHV"] = 2 - 2 * eps;
cScheme["FDH"] = 2;
scheme = "tHV"; (* ??? *)

alphas = 0.1;

eval[amp_] := amp[[1]] /. Rule[a_, b_] :> Rule[a,(-1/2)*Normal@Series[b /. amp[[2]] /. amp[[3]] /. momentumtwistorvalues /. replacemi,{eps, 0, 2}]];


(*Tree level*)
part={"d23T451","d14T253","d12T453","d34T251"};
helicity = {"+++-+", "++-++"};

a1pm = Get["ttj_amplitudes/ttqqg/ttqqg_0L_"<>part[[1]]<>"_"<>helicity[[1]]<>".m"];
a2pm = Get["ttj_amplitudes/ttqqg/ttqqg_0L_"<>part[[2]]<>"_"<>helicity[[1]]<>".m"];
a3pm = Get["ttj_amplitudes/ttqqg/ttqqg_0L_"<>part[[3]]<>"_"<>helicity[[1]]<>".m"];
a4pm = Get["ttj_amplitudes/ttqqg/ttqqg_0L_"<>part[[4]]<>"_"<>helicity[[1]]<>".m"];

e1pm = eval[a1pm];
e2pm = eval[a2pm];
e3pm = eval[a3pm];
e4pm = eval[a4pm];

a1mp = Get["ttj_amplitudes/ttqqg/ttqqg_0L_"<>part[[1]]<>"_"<>helicity[[2]]<>".m"];
a2mp = Get["ttj_amplitudes/ttqqg/ttqqg_0L_"<>part[[2]]<>"_"<>helicity[[2]]<>".m"];
a3mp = Get["ttj_amplitudes/ttqqg/ttqqg_0L_"<>part[[3]]<>"_"<>helicity[[2]]<>".m"];
a4mp = Get["ttj_amplitudes/ttqqg/ttqqg_0L_"<>part[[4]]<>"_"<>helicity[[2]]<>".m"];

e1mp = eval[a1mp];
e2mp = eval[a2mp];
e3mp = eval[a3mp];
e4mp = eval[a4mp];

e1pm = Rule[AA["+++-+"], phasepm*Normal@Series[Total[theta . (e1mp // Normal)[[All,2]]], {eps, 0, 2}]];
e2pm = Rule[AA["+++-+"], phasepm*Normal@Series[Total[theta . (e2mp // Normal)[[All,2]]], {eps, 0, 2}]];
e3pm = Rule[AA["+++-+"], phasepm*Normal@Series[Total[theta . (e3mp // Normal)[[All,2]]], {eps, 0, 2}]];
e4pm = Rule[AA["+++-+"], phasepm*Normal@Series[Total[theta . (e4mp // Normal)[[All,2]]], {eps, 0, 2}]];

e1mp = Rule[AA["++-++"], phasemp*Normal@Series[Total[theta . (e1mp // Normal)[[All,2]]], {eps, 0, 2}]];
e2mp = Rule[AA["++-++"], phasemp*Normal@Series[Total[theta . (e2mp // Normal)[[All,2]]], {eps, 0, 2}]];
e3mp = Rule[AA["++-\.08++"], phasemp*Normal@Series[Total[theta . (e3mp // Normal)[[All,2]]], {eps, 0, 2}]];
e4mp = Rule[AA["++-++"], phasemp*Normal@Series[Total[theta . (e4mp // Normal)[[All,2]]], {eps, 0, 2}]];

ampspm = {e1pm[[2]], e2pm[[2]], e3pm[[2]], e4pm[[2]]};
ampsmp = {e1mp[[2]], e2mp[[2]], e3mp[[2]], e4mp[[2]]};

epmsquare = Rule[AA2["+++-+"], Sum[Conjugate[ampspm[[i]]]*Gij[i, j]*ampspm[[j]], {i,1,4}, {j,1,4}]]
empsquare = Rule[AA2["++-++"], Sum[Conjugate[ampsmp[[i]]]*Gij[i, j]*ampsmp[[j]], {i,1,4}, {j,1,4}]]

ampsquare = epmsquare[[2]]+empsquare[[2]]





combds[partil_, ncpow_, hel_] := Module[{p=partil, n=ncpow, h=hel},

If[ Length[dspowlist[p]] != 2,
    amp1L = Get["ttj_amplitudes/ttqqg/ttqqg_1L_"<>p<>"_"<>n<>"_"<>h<>"_"<>"dsm2p0"<>".m"];
    eamp1L = eval[amp1L];];
If[ Length[dspowlist[p]] == 2,
	amp1L0 = Get["ttj_amplitudes/ttqqg/ttqqg_1L_"<>p<>"_"<>n<>"_"<>h<>"_"<>"dsm2p0"<>".m"];
    amp1L1 = Get["ttj_amplitudes/ttqqg/ttqqg_1L_"<>p<>"_"<>n<>"_"<>h<>"_"<>"dsm2p1"<>".m"];
    eamp1L0 = eval[amp1L0];
    eamp1L1 = eval[amp1L1];
    eamp1L = Table[ Rule[ eamp1L0[[i, 1]], Normal@Series[eamp1L0[[i, 2]] + cScheme[scheme] * eamp1L1[[i, 2]],{eps,0,2}]],{i, 1, Length[eamp1L0]}];];
    
    eamp1L];



ncfac[part_, ncpow_] := Module[{p=part, n=ncpow},

fac = 1;
If[ n == "Ncp0",
	fac = fac*1;];
If [n == "Ncp-1",
	fac = fac*1/3;];
If [n == "Ncp1",
	fac = fac*3;];
If [n == "Ncp-2",
	fac = fac*1/9;];
If [n == "Ncp2",
	fac = fac*9;];
If[ StringTake[p, {1,2}] == "Nf",
	fac = fac*5;];
If[ StringTake[p, {1,2}] == "Nh",
	fac = fac*1;];
	
fac];



combnc[partil_, hel_] := Module[{p=partil, h=hel},

ampa = combds[p, First[ncpowlist[p]], h];
ampb = combds[p, Last[ncpowlist[p]], h];

faca = ncfac[p, First[ncpowlist[p]]];
facb = ncfac[p, Last[ncpowlist[p]]];
If[ Length[ncpowlist[p]] == 1,
	facb = 0;];

eamp1L = Table[Rule[ ampa[[i, 1]], Normal@Series[faca*ampa[[i, 2]] + facb*ampb[[i, 2]], {eps, 0, 2}]],{i, 1, Length[ampa]}];

eamp1L];



part = "Nhd23T451";

eee = combnc[part, "++-++"];

Do[
Print[i," ",Collect[eee[[i,2]],eps,Chop[N[Chop[#],5]]&] /. eps->\[Epsilon]];
,{i,1,4}];















(* numeric 4-vectors (massless) *)
p3v = {500., 0., 0.,  500.};
p4v = {500., 0., 0., -500.};

Clear[p3, p4];

DeclareSpinor[p3, p4];
DeclareSpinorMomentum[p3, p3v];
DeclareSpinorMomentum[p4, p4v];

(* try again \[Dash] these should now be numbers (possibly complex) *)
N@Spaa[p3, p4]
N@Spbb[p3, p4]

(* sanity check for massless: \:27e8ij\:27e9[ji] = 2 p_i\[CenterDot]p_j *)
mdot[a_, b_] := a[[1]] b[[1]] - a[[2]] b[[2]] - a[[3]] b[[3]] - a[[4]] b[[4]];
{ N[Spaa[p3,p4] Spbb[p4,p3]], N[2 mdot[p3v,p4v]] }



