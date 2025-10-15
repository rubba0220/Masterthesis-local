(* ::Package:: *)

ClearAll;

$dir = If[$FrontEnd=!=Null,NotebookDirectory[],DirectoryName[$InputFileName]];
SetDirectory[$dir];

Get[ToFileName[{"C:/Users/Laptop von Ruben/AppData/Roaming/Wolfram/Applications/Spinors-1.0"},"init.windows" ]];
<<Spinors`

Get["C:/Users/Laptop von Ruben/AppData/Roaming/Wolfram/Applications/SU3basics.m"];


(* include integral evaluations *)
integraleval = Get["MIs_evaluation/T1T2T3T4_all_Permutations.m"] /. minorm[__]:>1 /. \[Epsilon]->eps;


mdot[a_, b_] := a[[1]]b[[1]]-a[[2]]b[[2]]-a[[3]]b[[3]]-a[[4]]b[[4]];

p3v = -{243.958874340707, 0, 0, 243.958874340707};
p4v = -{243.274422280898, 0, 0, -243.274422280898};
p5v = SetPrecision[{72.4210128964652, 50.6623055703701, 44.0522222052368, 27.1576070747199},50];
p1v = {202.733585822752, -81.9445526000469, 60.951575089163, 30.2432242639469};
p2v = {212.078697902388, 31.2822470296768, -105.0037972944, -56.7163792788578};

(*mdot[p5v,p5v]*)
(*p1v = {451.16038194215792601, -254.95547824215555011, -127.96321306901207038, -303.88644976300315648};
p2v = {181.67424764458991149, 46.833603697807077992, -3.1743236745347509498, 31.308679165908877451};
p3v = -{500., 0, 0, 500.};
p4v = -{500., 0, 0, -500.};
p5v = {367.16537041325221935, 208.12187454434831579, 131.13753674354686041, 272.57777059709434297};*)

mtv = Sqrt[mdot[p1v,p1v]];

ni = SetPrecision[{100.,100.,0,0},50];
Rot[a_, b_, c_] := {{1, 0, 0, 0}, {0, Cos[b]Cos[c], Sin[a]Sin[b]Cos[c]-Cos[a]Sin[c], Cos[a]Sin[b]Cos[c]+Sin[a]Sin[c]},{0, Cos[b]Sin[c], Sin[a]Sin[b]Sin[c]+Cos[a]Cos[c], Cos[a]Sin[b]Sin[c]-Sin[a]Cos[c]},{0, -Sin[b], Sin[a]Cos[b], Cos[a]Cos[b]}};
rot1v = N[Rot[SetPrecision[0.,50],SetPrecision[0.4,50],SetPrecision[0.,50]],50];
rot2v = N[Rot[SetPrecision[0.2,50],SetPrecision[0.9,50],SetPrecision[0.5,50]],50];
(*rot2v = rot1v;*)

n1v = rot1v . ni;
n2v = rot2v . ni;

f1v = SetPrecision[p1v - mdot[p1v, p1v]/(2*mdot[p1v, n1v]) * n1v,50];
f2v = SetPrecision[p2v - mdot[p2v, p2v]/(2*mdot[p2v, n2v]) * n2v,50];

projMassless[v_] := With[{sp = v[[2;;4]], E = Norm[v[[2;;4]]]}, Flatten[{E, sp}]];

(*f1v = projMassless[f1v];
f2v = projMassless[f2v];
p5v = projMassless[p5v];*)


(*result is n1, n2 dependent ???*)
(*aaaa = 70.7106781186547524400844362104849039284835937688474036588339868995366239231;
n1vt = {100, aaaa, aaaa, 0};
n2vt = {100, aaaa, -aaaa, 0};
f1vt = p1v - mdot[p1v, p1v]/(2*mdot[p1v, n1vt]) * n1vt;
f2vt = p2v - mdot[p2v, p2v]/(2*mdot[p2v, n2vt]) * n2vt;*)

s34v = mdot[p3v+p4v, p3v+p4v];
Print["s34 = ", s34v];

Print["Sum p = ", p1v+p2v+p3v+p4v+p5v];

Do[
Print["sqrt pi^2 = ", Sqrt[mdot[p,p]]];
,{p,{p1v,p2v,p3v,p4v,p5v}}];

Do[
Print["sqrt ni/fi^2 = ",Sqrt[mdot[p,p]]];
,{p,{n1v,n2v,f1v,f2v}}]


replacemi = integraleval[[2]];
(*dijvalues2 = integraleval[[1]];*)

dijvalues = {d12->mdot[p1v,p2v], d23->mdot[p2v,p3v], d34->mdot[p3v,p4v], d45->mdot[p4v,p5v], d15->mdot[p5v,p1v], mt2->mdot[p1v,p1v]};

(*Print[dijvalues2];
Print[dijvalues];*)

Clear[f1,f2,n1,n2, p3, p4, p5, m3, m4, m5];
spin = {f1,f2,n1,n2,p3,p4,p5};
spinv = {f1v,f2v,n1v,n2v,p3v,p4v,p5v};
(*spinv = {f1vt,f2vt,n1vt,n2vt,p3v,p4v,p5v};*)

(*final ampsquare scales with norm2 --> scales with Mass^-2 (expected for 2->5)*)
norm2 = 1.;
norm = Sqrt[norm2];

spinvn = spinv/norm;
mtvn = mtv/norm;
s34vn = s34v/norm2;
Print["normalized s34 = ", s34vn];
dijvaluesn = dijvalues /. Rule[a_, b_] :> Rule[a,b/norm2];
Print["normalized dij = ", dijvalues];

dij2mt = {
d12 -> -(s34*t12*(-t45 + 2*t23*t51 - 2*t23*x5123 + 2*t12*t23*x5123 + 2*t45*x5123 + 2*t12*t51*x5123 - 2*t45*t51*x5123 - 2*t12*x5123^2 + 2*t12^2*x5123^2 - 2*t12*t45*x5123^2))/(2*t45),
d15 -> (s34*t51)/2,
d23 -> (s34*t23)/2,
d34 -> s34/2,
d45 -> (s34*t45)/2,
mt2 -> (s34*t12*(t23*t51 - t23*x5123 + t12*t23*x5123 + t45*x5123 + t12*t51*x5123 - t45*t51*x5123 - t12*x5123^2 + t12^2*x5123^2 - t12*t45*x5123^2))/t45
};
mt2dij = Solve[dij2mt /. Rule->Equal,{s34,t45,t51,t12,t23,x5123}] // Simplify;
momentumtwistorvalues = mt2dij[[1]] /. dijvaluesn // FullSimplify;


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


DeclareSpinor[f1,f2,n1,n2,p3,p4,p5]

Do[
DeclareSpinorMomentum[spin[[i]], spinvn[[i]]];
,{i,{1,2,3,4,5,6,7}}];

phasemp = mtvn*Spaa[p3,p4]/(Spaa[p4,p5]^2)
phasepm = mtvn*Spaa[p3,p4]/(Spaa[p3,p5]^2)

Print["Phi = ", phasemp//N];
Print["Spinorprod = ", Spaa[f1,p5]//N];
Print["Spinorprod- = ", Spaa[p4,p3]//N];
Print["Angleprod = ", Spbb[p3,p4]//N];
Print[Spaa[f1,p5]*Spbb[p5,f1]/(2*mdot[f1v,p5v])//N];
Print[Spaa[p3,p4]*Spbb[p4,p3]/(2*mdot[p3v,p4v])//N];
Print["Spinorprod pp = ", Spbb[p3,p4,p5,p3]//N];

theta1 = Spaa[n1,n2]*s34vn/Spaa[f1,n1]/Spaa[f2,n2]
theta2 = Spaa[n1,p3]*Spaa[n2,p4]*Spbb[p3,p4]/Spaa[f1,n1]/Spaa[f2,n2]
theta3 = Spaa[n1,p3]*Spaa[n2,p3]*Spbb[p3,p4,p5,p3]/s34vn/Spaa[f1,n1]/Spaa[f2,n2]
theta4 = Spaa[n1,p4]*Spaa[n2,p4]*Spbb[p4,p5,p3,p4]/s34vn/Spaa[f1,n1]/Spaa[f2,n2]

theta = {theta1, theta2, theta3, theta4};

(*test: is declaring Spinor momentum sufficient*)
(*DeclareLVector[m3,m4,m5];
DeclareLVectorMomentum[m4,p4v/Sqrt[s34v]];
DeclareLVectorMomentum[m5,p5v/Sqrt[s34v]];
theta33 = Spaa[n1,p3]*Spaa[n2,p3]*Spbb[p3,m4,m5,p3]/1/Spaa[f1,n1]/Spaa[f2,n2]//N;

Print[theta33];
Print[theta3];*)

(* Gram to test n1,n2 dependence *)
 (*Print[Outer[Times, Conjugate[theta], theta]//N];*)


Clear[CT];
CT[1][i1_, i2_, i3_, i4_, a5_] := KroneckerDelta[i1, i4]*T[a5][[i3, i2]];
CT[2][i1_, i2_, i3_, i4_, a5_] := KroneckerDelta[i2, i3]*T[a5][[i1, i4]];
CT[3][i1_, i2_, i3_, i4_, a5_] := -1/Nc * KroneckerDelta[i1, i2]*T[a5][[i3, i4]];
CT[4][i1_, i2_, i3_, i4_, a5_] := -1/Nc * KroneckerDelta[i3, i4]*T[a5][[i1, i2]];

Gij[i_,j_] := Sum[Conjugate[CT[i][i1, i2, i3, i4, a5]]*CT[j][i1, i2, i3, i4, a5]
,{i1,1,3},{i2,1,3},{i3,1,3},{i4,1,3},{a5,1,8}];

Gijmatrix = {{0,0,0,0}, {0,0,0,0}, {0,0,0,0}, {0,0,0,0}};
Do[
Gijmatrix[[a]][[b]] = Gij[a,b];
,{a,{1,2,3,4}},{b,{1,2,3,4}}];
Print["Color metric = ",Gijmatrix]


alphas = 1/4/Pi;
(* mu *)

gs = Sqrt[4 Pi alphas];
Neps = Exp[EulerGamma*eps] * (Gamma[1-eps]^2) * Gamma[1+eps] / ((4*Pi)^(2-eps)) / Gamma[1-2eps];
(*Print[Normal@Series[Neps,{eps,0,2}] /. eps->\[Epsilon]];*)


cScheme["tHV"] = 2 - 2 * eps;
cScheme["FDH"] = 2;
scheme = "tHV"; (* ??? *)

eval[amp_] := amp[[1]] /. Rule[a_, b_] :> Rule[a,(-1/2)*Normal@Series[b /. amp[[2]] /. amp[[3]] /. momentumtwistorvalues /.replacemi,{eps, 0, 2}]];


(*Tree level*)
part={"d14T253","d23T451","d12T453","d34T251"}; (*d for delta, T for generator*)
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
Print["partial amp = ", (e1mp//Normal)[[All,2]], e1mp];

e1pm = phasepm*theta . (e1pm // Normal)[[All,2]];
e2pm = phasepm*theta . (e2pm // Normal)[[All,2]];
e3pm = phasepm*theta . (e3pm // Normal)[[All,2]];
e4pm = phasepm*theta . (e4pm // Normal)[[All,2]];

e1mp = phasemp*theta . (e1mp // Normal)[[All,2]];
e2mp = phasemp*theta . (e2mp // Normal)[[All,2]];
e3mp = phasemp*theta . (e3mp // Normal)[[All,2]];
e4mp = phasemp*theta . (e4mp // Normal)[[All,2]];

ampspm = {e1pm, e2pm, e3pm, e4pm}//N;
ampsmp = {e1mp, e2mp, e3mp, e4mp}//N;
Print["amps up to color ", ampspm, ampsmp];

epm2 = Rule[M["+++-+"], Sum[Conjugate[ampspm[[i]]]*Gij[i, j]*ampspm[[j]], {i,1,4}, {j,1,4}]]
emp2 = Rule[M["++-++"], Sum[Conjugate[ampsmp[[i]]]*Gij[i, j]*ampsmp[[j]], {i,1,4}, {j,1,4}]]

avampsquare = 1/2 * 1/2 * 1/3 * 1/3 * (gs^(3+0)*Neps^0)^2 (epm2[[2]]+emp2[[2]])
ampsquare = (gs^(3+0)*Neps^0)^2 (epm2[[2]]+emp2[[2]])
(* q,q_bar are incoming in openloops: 2 possible helicities, 3 possible colors *)



(*(* Loop contributions *)

combds[partil_, ncpow_, hel_] := Module[{p=partil, n=ncpow, h=hel},

If[ Length[dSpowlist[p]] != 2,
    amp1L = Get["ttj_amplitudes/ttqqg/ttqqg_1L_"<>p<>"_"<>n<>"_"<>h<>"_"<>"dsm2p0"<>".m"];
    eamp1L = eval[amp1L];];
If[ Length[dSpowlist[p]] == 2,
	amp1L0 = Get["ttj_amplitudes/ttqqg/ttqqg_1L_"<>p<>"_"<>n<>"_"<>h<>"_"<>"dsm2p0"<>".m"];
    amp1L1 = Get["ttj_amplitudes/ttqqg/ttqqg_1L_"<>p<>"_"<>n<>"_"<>h<>"_"<>"dsm2p1"<>".m"];
    eamp1L0 = eval[amp1L0];
    eamp1L1 = eval[amp1L1];
    eamp1L = Table[ Rule[ eamp1L0[[i, 1]], Normal@Series[eamp1L0[[i, 2]] + cScheme[scheme] * eamp1L1[[i, 2]],{eps,0,2}]],{i, 1, Length[eamp1L0]}];];
    
    eamp1L];*)



(*ncfac[part_, ncpow_] := Module[{p=part, n=ncpow},

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
*)


(*combnc[partil_, hel_] := Module[{p=partil, h=hel},

ampa = combds[p, First[ncpowlist[p]], h];
ampb = combds[p, Last[ncpowlist[p]], h];

faca = ncfac[p, First[ncpowlist[p]]];
facb = ncfac[p, Last[ncpowlist[p]]];
If[ Length[ncpowlist[p]] == 1,
	facb = 0;];

eamp1L = Table[Rule[ ampa[[i, 1]], Normal@Series[faca*ampa[[i, 2]] + facb*ampb[[i, 2]], {eps, 0, 2}]],{i, 1, Length[ampa]}];

eamp1L];*)



(*part = "Nhd23T451";

eee = combnc[part, "++-++"];

Do[
Print[i," ",Collect[eee[[i,2]],eps,Chop[N[Chop[#],5]]&] /. eps->\[Epsilon]];
,{i,1,4}];


*)












(*(* numeric 4-vectors (massless) *)
p3v = {500., 0., 0.,  500.};
p4v = {500., 0., 0., -500.};

Clear[p3, p4];

DeclareSpinor[p3, p4];
DeclareSpinorMomentum[p3, p3v];
DeclareSpinorMomentum[p4, p4v];

(* try again \[Dash] these should now be numbers (possibly complex) *)
N@Spaa[p3, p4]
N@Spbbb[p3, p4]

(* sanity check for massless: \:27e8ij\:27e9[ji] = 2 p_i\[CenterDot]p_j *)
mdot[a_, b_] := a[[1]] b[[1]] - a[[2]] b[[2]] - a[[3]] b[[3]] - a[[4]] b[[4]];
{ N[Spaa[p3,p4] Spbbb[p4,p3]], N[2 mdot[p3v,p4v]] }
*)


