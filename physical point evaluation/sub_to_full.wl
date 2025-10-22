(* ::Package:: *)

$dir = If[$FrontEnd=!=Null,NotebookDirectory[],DirectoryName[$InputFileName]];
SetDirectory[$dir];

Get[ToFileName[{"C:/Users/Laptop von Ruben/AppData/Roaming/Wolfram/Applications/Spinors-1.0"},"init.windows" ]];
<<Spinors`

Get["C:/Users/Laptop von Ruben/AppData/Roaming/Wolfram/Applications/SU3basics.m"];


(* include integral evaluations *)
(*integraleval = Get["MIs_evaluation/T1T2T3T4_all_Permutations.m"] /. minorm[__]:>1 /. \[Epsilon]->eps;*)

(* momentum related functions *)
mdot[a_, b_] := a[[1]]b[[1]]-a[[2]]b[[2]]-a[[3]]b[[3]]-a[[4]]b[[4]];
projMassless[v_] := With[{sp = v[[2;;4]], E = Norm[v[[2;;4]]]}, Flatten[{E, sp}]];


(*2-to-3 momenta converted to 0-to-5*)
sv = Sqrt[mdot[{243.274422280898, 0, 0, -243.274422280898}, {243.958874340707, 0, 0, 243.958874340707}]*2]
p4v = -{243.958874340707, 0, 0, 243.958874340707};
p3v = -{243.274422280898, 0, 0, -243.274422280898};
p5v = {72.4210128964652, 50.6623055703701, 44.0522222052368, 27.1576070747199};
p1v = {202.733585822752, -81.9445526000469, 60.951575089163, 30.2432242639469};
p2v = {212.078697902388, 31.2822470296768, -105.0037972944, -56.7163792788578};

(*p1v = {451.16038194215792601, -254.95547824215555011, -127.96321306901207038, -303.88644976300315648};
p2v = {181.67424764458991149, 46.833603697807077992, -3.1743236745347509498, 31.308679165908877451};
p3v = -{500., 0, 0, 500.};
p4v = -{500., 0, 0, -500.};
p5v = {367.16537041325221935, 208.12187454434831579, 131.13753674354686041, 272.57777059709434297};*)

ni = SetPrecision[{100.,100.,0,0},50];
Rot[a_, b_, c_] := {{1, 0, 0, 0}, {0, Cos[b]Cos[c], Sin[a]Sin[b]Cos[c]-Cos[a]Sin[c], Cos[a]Sin[b]Cos[c]+Sin[a]Sin[c]},{0, Cos[b]Sin[c], Sin[a]Sin[b]Sin[c]+Cos[a]Cos[c], Cos[a]Sin[b]Sin[c]-Sin[a]Cos[c]},{0, -Sin[b], Sin[a]Cos[b], Cos[a]Cos[b]}};
rot1v = N[Rot[SetPrecision[0.5,50],SetPrecision[0.1,50],SetPrecision[0.9,50]],50];
rot2v = N[Rot[SetPrecision[0.6,50],SetPrecision[0.9,50],SetPrecision[0.2,50]],50];

n1v = rot1v . ni;
n2v = rot2v . ni;
(*n1v = projMassless[p1v];
n2v = projMassless[p2v];*)

n1v = mdot[p1v, p1v]/(2*mdot[p1v, n1v]) * n1v;
n2v = mdot[p2v, p2v]/(2*mdot[p2v, n2v]) * n2v;

f1v = p1v - mdot[p1v, p1v]/(2*mdot[p1v, n1v]) * n1v;
f2v = p2v - mdot[p2v, p2v]/(2*mdot[p2v, n2v]) * n2v;

(*(* with higher precision *)
f1v = projMassless[f1v];
f2v = projMassless[f2v];
p5v = projMassless[p5v];*)

(*Print["Sum p = ", p1v+p2v+p3v+p4v+p5v];
*)
Do[
Print["sqrt pi^2 = ", Sqrt[mdot[p,p]]];
,{p,{p1v,p2v,p3v,p4v,p5v}}];
(*
Do[
Print["sqrt ni or fi^2 = ",Sqrt[mdot[p,p]]];
,{p,{n1v,n2v,f1v,f2v}}]*)


(*replacemi = integraleval[[2]];*)
(*dijvalues2 = integraleval[[1]];*)
dijvalues = {d12->mdot[p1v,p2v], d23->mdot[p2v,p3v], d34->mdot[p3v,p4v], d45->mdot[p4v,p5v], d15->mdot[p5v,p1v], mt2->mdot[p1v,p1v]}
Export["dij.m",dijvalues];

dij2mt = {
d12 -> -(s34*t12*(-t45 + 2*t23*t51 - 2*t23*x5123 + 2*t12*t23*x5123 + 2*t45*x5123 + 2*t12*t51*x5123 - 2*t45*t51*x5123 - 2*t12*x5123^2 + 2*t12^2*x5123^2 - 2*t12*t45*x5123^2))/(2*t45),
d15 -> (s34*t51)/2,
d23 -> (s34*t23)/2,
d34 -> s34/2,
d45 -> (s34*t45)/2,
mt2 -> (s34*t12*(t23*t51 - t23*x5123 + t12*t23*x5123 + t45*x5123 + t12*t51*x5123 - t45*t51*x5123 - t12*x5123^2 + t12^2*x5123^2 - t12*t45*x5123^2))/t45
};
mt2dij = Solve[dij2mt /. Rule->Equal,{s34,t45,t51,t12,t23,x5123}] // Simplify; (* first solution spikey brackets, second solution square brackets *)
momentumtwistorvalues = mt2dij[[1]] /. dijvalues // FullSimplify;

DeclareSpinor[p5,p3];
DeclareSpinorMomentum[p5,p5v];
DeclareSpinorMomentum[p3,p3v];
DeclareLVector[l1, l45];
DeclareLVectorMomentum[l1, p1v];
DeclareLVectorMomentum[l45, p4v+p5v];
mtvaluesspike = {s34->mdot[p3v+p4v,p3v+p4v], t12->mdot[p1v+p2v,p1v+p2v]/mdot[p3v+p4v,p3v+p4v], t23->(mdot[p3v+p2v,p3v+p2v] - mdot[p1v,p1v])/mdot[p3v+p4v,p3v+p4v], t45->mdot[p5v+p4v,p5v+p4v]/mdot[p3v+p4v,p3v+p4v], t15->(mdot[p1v+p5v,p1v+p5v] - mdot[p1v,p1v])/mdot[p3v+p4v,p3v+p4v], x5123->-Spaa[p5,l1,l45,p3]/Spaa[p5,p3]/mdot[p1v+p2v,p1v+p2v]} //N
mtvaluescube = {s34->mdot[p3v+p4v,p3v+p4v], t12->mdot[p1v+p2v,p1v+p2v]/mdot[p3v+p4v,p3v+p4v], t23->(mdot[p3v+p2v,p3v+p2v] - mdot[p1v,p1v])/mdot[p3v+p4v,p3v+p4v], t45->mdot[p5v+p4v,p5v+p4v]/mdot[p3v+p4v,p3v+p4v], t15->(mdot[p1v+p5v,p1v+p5v] - mdot[p1v,p1v])/mdot[p3v+p4v,p3v+p4v], x5123->-Spbb[p5,l1,l45,p3]/Spbb[p5,p3]/mdot[p1v+p2v,p1v+p2v]} //N
Print[mt2dij /.dijvalues]
Sqrt[MP[l1,l1]]//N


Clear[f1,f2,n1,n2, p3, p4, p5, tmp1, tmp2, tmp];
spin = {f1,f2,n1,n2,p3,p4,p5};
spinv = {f1v,f2v,n1v,n2v,p3v,p4v,p5v};

(*(* final ampsquare scales with norm2 --> scales with Mass^-2 (expected for 2->5) *)
norm2 = 1.;
norm = Sqrt[norm2];

spinvn = spinv/norm;
mtvn = mtv/norm;
s34vn = s34v/norm2;
Print["normalized s34 = ", s34vn];
dijvaluesn = dijvalues /. Rule[a_, b_] :> Rule[a,b/norm2];
Print["normalized dij = ", dijvalues]; *)

DeclareSpinor[f1,f2,n1,n2,p3,p4,p5,tmp1, tmp2, tmp];

Do[
DeclareSpinorMomentum[spin[[i]], spinv[[i]]];
,{i,{1,2,3,4,5,6,7}}];

Clear[global, phase];
global = Sqrt[mt2];

(* signs do not change anything here *)
phase["-++"] = Spaa[p3,p4]/(Spaa[p4,p5]^2);
phase["+-+"] = Spaa[p3,p4]/(Spaa[p3,p5]^2);
phase["-+-"] = Spbb[p3,p4]/(Spbb[p3,p5]^2);
phase["+--"] = Spbb[p3,p4]/(Spbb[p4,p5]^2);
(*phase["+-+"] = Spaa[p3,p4]/(Spaa[p4,p5]^2);
phase["-++"] = Spaa[p3,p4]/(Spaa[p3,p5]^2);
phase["+--"] = Spbb[p3,p4]/(Spbb[p3,p5]^2); 
phase["-+-"] = Spbb[p3,p4]/(Spbb[p4,p5]^2);*)

(*Print["Phi = ", phasemp//N];
Print["Spinorprod = ", Spaa[f1,p5]//N];
Print["Spinorprod- = ", Spaa[p4,p3]//N];
Print["Angleprod = ", Spbb[p3,p4]//N];
Print[Spaa[f1,p5]*Spbb[p5,f1]/(2*mdot[f1v,p5v])//N];
Print[Spaa[p3,p4]*Spbb[p4,p3]/(2*mdot[p3v,p4v])//N];
Print["Spinorprod pp = ", Spbb[p3,p4,p5,p3]//N];*)

Clear[theta];
theta[1] = Spaa[n1,n2]*s34/Spaa[f1,n1]/Spaa[f2,n2];
theta[2] = Spaa[n1,p3]*Spaa[n2,p4]*Spbb[p3,p4]/Spaa[f1,n1]/Spaa[f2,n2];
theta[3] = Spaa[n1,p3]*Spaa[n2,p3]*Spbb[p3,p4,p5,p3]/s34/Spaa[f1,n1]/Spaa[f2,n2];
theta[4] = Spaa[n1,p4]*Spaa[n2,p4]*Spbb[p4,p5,p3,p4]/s34/Spaa[f1,n1]/Spaa[f2,n2];

theta[3] * Spaa[n1,f1]/Sqrt[mt2] * Spaa[n2,f2]/Sqrt[mt2]/. {f1->tmp1, n1->f1, f2->tmp2, n2->f2} /. {tmp1->n1, tmp2->n2} /. mt2dij[[1]] /.dijvalues //N
theta[3]/. mt2dij[[1]] /.dijvalues //N
Spbb[n1,p3]*Spbb[n2,p3]*Spaa[p3,p4,p5,p3]/s34/Spbb[f1,n1]/Spbb[f2,n2]/. mt2dij[[1]] /.dijvalues //N
Spbb[n1,p3]*Spbb[n2,p3]*Spaa[p3,p4,p5,p3]/s34/Spbb[f1,n1]/Spbb[f2,n2]* Spbb[n1,f1]/Sqrt[mt2] * Spbb[n2,f2]/Sqrt[mt2]/. {f1->tmp1, n1->f1, f2->tmp2, n2->f2} /. {tmp1->n1, tmp2->n2}/. mt2dij[[1]] /.dijvalues //N


Abs[Spaa[n1, f1]/Sqrt[mt2]/.dijvalues //N]^2
mdot[f1v, n1v]*2/mt2 /.dijvalues

(*Spaa[n1,n2]*s34/Spaa[f1,n1]/Spaa[f2,n2]* Spaa[n1,f1]/Sqrt[mt2] * Spaa[n2,f2]/Sqrt[mt2] /. {f1->tmp1, n1->f1, f2->tmp2, n2->f2} /. {tmp1->n1, tmp2->n2}
Spaa[n1,p3]*Spaa[n2,p4]*Spbb[p3,p4]/Spaa[f1,n1]/Spaa[f2,n2]* Spaa[n1,f1]/Sqrt[mt2] * Spaa[n2,f2]/Sqrt[mt2] /. {f1->tmp1, n1->f1, f2->tmp2, n2->f2} /. {tmp1->n1, tmp2->n2}/.momentumtwistorvalues /.dijvalues//N
Spaa[n1,p3]*Spaa[n2,p3]*Spbb[p3,p4,p5,p3]/s34/Spaa[f1,n1]/Spaa[f2,n2]* Spaa[n1,f1]/Sqrt[mt2] * Spaa[n2,f2]/Sqrt[mt2] /. {f1->tmp1, n1->f1, f2->tmp2, n2->f2} /. {tmp1->n1, tmp2->n2}/.momentumtwistorvalues /.dijvalues//N
Spaa[n1,p4]*Spaa[n2,p4]*Spbb[p4,p5,p3,p4]/s34/Spaa[f1,n1]/Spaa[f2,n2]* Spaa[n1,f1]/Sqrt[mt2] * Spaa[n2,f2]/Sqrt[mt2] /. {f1->tmp1, n1->f1, f2->tmp2, n2->f2} /. {tmp1->n1, tmp2->n2}/.momentumtwistorvalues /.dijvalues//N

Spaa[n1,n2]*s34/Spaa[f1,n1]/Spaa[f2,n2]* Spaa[n1,f1]/Sqrt[mt2]  /. {f1->tmp1, n1->f1} /. {tmp1->n1}
Spaa[n1,p3]*Spaa[n2,p4]*Spbb[p3,p4]/Spaa[f1,n1]/Spaa[f2,n2]* Spaa[n1,f1]/Sqrt[mt2]  /. {f1->tmp1, n1->f1} /. {tmp1->n1}/.momentumtwistorvalues /.dijvalues//N
Spaa[n1,p3]*Spaa[n2,p3]*Spbb[p3,p4,p5,p3]/s34/Spaa[f1,n1]/Spaa[f2,n2]* Spaa[n1,f1]/Sqrt[mt2]/. {f1->tmp1, n1->f1} /. {tmp1->n1}/.momentumtwistorvalues /.dijvalues//N
Spaa[n1,p4]*Spaa[n2,p4]*Spbb[p4,p5,p3,p4]/s34/Spaa[f1,n1]/Spaa[f2,n2]* Spaa[n1,f1]/Sqrt[mt2] /. {f1->tmp1, n1->f1} /. {tmp1->n1}/.momentumtwistorvalues/.dijvalues //N

Spaa[n1,n2]*s34/Spaa[f1,n1]/Spaa[f2,n2] * Spaa[n2,f2]/Sqrt[mt2] /. {f2->tmp2, n2->f2} /. {tmp2->n2}
Spaa[n1,p3]*Spaa[n2,p4]*Spbb[p3,p4]/Spaa[f1,n1]/Spaa[f2,n2] * Spaa[n2,f2]/Sqrt[mt2] /. {f2->tmp2, n2->f2} /. {tmp2->n2}/.momentumtwistorvalues/.dijvalues //N
Spaa[n1,p3]*Spaa[n2,p3]*Spbb[p3,p4,p5,p3]/s34/Spaa[f1,n1]/Spaa[f2,n2] * Spaa[n2,f2]/Sqrt[mt2] /. {f2->tmp2, n2->f2} /. {tmp2->n2}/.momentumtwistorvalues /.dijvalues//N
Spaa[n1,p4]*Spaa[n2,p4]*Spbb[p4,p5,p3,p4]/s34/Spaa[f1,n1]/Spaa[f2,n2] * Spaa[n2,f2]/Sqrt[mt2] /. {f2->tmp2, n2->f2} /. {tmp2->n2}/.momentumtwistorvalues /.dijvalues//N
*)
(*test: is declaring Spinor momentum sufficient*)
(*DeclareLVector[m3,m4,m5];
DeclareLVectorMomentum[m4,p4v/Sqrt[s34v]];
DeclareLVectorMomentum[m5,p5v/Sqrt[s34v]];
theta33 = Spaa[n1,p3]*Spaa[n2,p3]*Spbb[p3,m4,m5,p3]/1/Spaa[f1,n1]/Spaa[f2,n2]//N;

Print[theta33];
Print[theta3];*)

(* Gram to test n1,n2 dependence --> kompletter Unfug*)
 (*Print[Outer[Times, Conjugate[theta], theta]//N];*)


Clear[CT];
CT[1][i1_, i2_, i3_, i4_, a5_] := KroneckerDelta[i1, i4]*T[a5][[i3, i2]];
CT[2][i1_, i2_, i3_, i4_, a5_] := KroneckerDelta[i2, i3]*T[a5][[i1, i4]];
CT[3][i1_, i2_, i3_, i4_, a5_] := -1/Nc * KroneckerDelta[i1, i2]*T[a5][[i3, i4]];
CT[4][i1_, i2_, i3_, i4_, a5_] := -1/Nc * KroneckerDelta[i3, i4]*T[a5][[i1, i2]];

Gij[i_,j_] := Sum[Conjugate[CT[i][i1, i2, i3, i4, a5]]*CT[j][i1, i2, i3, i4, a5]
,{i1,1,3},{i2,1,3},{i3,1,3},{i4,1,3},{a5,1,8}];

(*CTl[1][i1_, i2_, i3_, i4_, a5_] := KroneckerDelta[i1, i4]*T[a5][[i3, i2]];
CTl[2][i1_, i2_, i3_, i4_, a5_] := KroneckerDelta[i2, i3]*T[a5][[i1, i4]];
CTl[3][i1_, i2_, i3_, i4_, a5_] := -KroneckerDelta[i1, i2]*T[a5][[i3, i4]];
CTl[4][i1_, i2_, i3_, i4_, a5_] := -KroneckerDelta[i3, i4]*T[a5][[i1, i2]];

Gijl[i_,j_] := Sum[Conjugate[CTl[i][i1, i2, i3, i4, a5]]*CTl[j][i1, i2, i3, i4, a5]
,{i1,1,3},{i2,1,3},{i3,1,3},{i4,1,3},{a5,1,8}];*)


Gijmatrix = {{0,0,0,0}, {0,0,0,0}, {0,0,0,0}, {0,0,0,0}};
Do[
Gijmatrix[[a]][[b]] = Gij[a,b];
,{a,{1,2,3,4}},{b,{1,2,3,4}}];
Print["Color metric = ",Gijmatrix]


alphas = 1/4/Pi;
(* mu *)

gs = Sqrt[4 Pi alphas]
Neps = Exp[EulerGamma*eps] * (Gamma[1-eps]^2) * Gamma[1+eps] / ((4*Pi)^(2-eps)) / Gamma[1-2eps]; (*todo: fix convention with EulerGamma and log(4Pi)*)
(*Print[Normal@Series[Neps,{eps,0,2}] /. eps->\[Epsilon]];*)


cScheme["tHV"] = 2 - 2 * eps;
cScheme["FDH"] = 2;
scheme = "tHV"; (* ??? *)

evaltree[amp_] := amp[[1]] /. Rule[a_, b_] :> Rule[a, b /. amp[[2]] /. amp[[3]]];
evalloop[amp_] := amp[[1]] /. Rule[a_, b_] :> (-1./2.)*Normal@Series[b /. amp[[2]] /. amp[[3]] /. momentumtwistorvalues /.replacemi,{eps, 0, 2}];


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

helicity = {"+++-+", "++-++"};
hel1qqg  = {"+-+", "-++"};
hel2qqg = {"-+-", "+--"}; (*parity*)
heltt = {"++", "+-", "-+", "--"}; (*paper [30]*)


(*Tree level*)
part["0L"]={"d14T253","d23T451","d12T453","d34T251"}; (*d for delta, T for generator*)
level1 = "0L";

Do[
amp[level1, part[level1][[i]], helicity[[j]]] = Get["ttj_amplitudes/ttqqg/ttqqg_0L_"<>part[level1][[i]]<>"_"<>helicity[[j]]<>".m"];
eamp[level1, part[level1][[i]], helicity[[j]]] = evaltree[amp[level1, part[level1][[i]], helicity[[j]]]];

camp[level1, part[level1][[i]], heltt[[1]], hel1qqg[[j]]] = global * phase[hel1qqg[[j]]] * Sum[theta[l] * eamp[level1, part[level1][[i]], helicity[[j]]][[l,2]],{l, 1, 4}]; 
camp[level1, part[level1][[i]], heltt[[1]], hel2qqg[[j]]] = global * phase[hel2qqg[[j]]] * Sum[theta[l] * eamp[level1, part[level1][[i]], helicity[[j]]][[l,2]],{l, 1, 4}];

camp[level1, part[level1][[i]], heltt[[2]], hel1qqg[[j]]] = camp[level1, part[level1][[i]], heltt[[1]], hel1qqg[[j]]] * Spaa[n2,f2]/Sqrt[mt2] /. f2->tmp /. n2->f2 /. tmp->n2;
camp[level1, part[level1][[i]], heltt[[2]], hel2qqg[[j]]] = camp[level1, part[level1][[i]], heltt[[1]], hel2qqg[[j]]] * Spaa[n2,f2]/Sqrt[mt2] /. f2->tmp /. n2->f2 /. tmp->n2;

camp[level1, part[level1][[i]], heltt[[3]], hel1qqg[[j]]] = camp[level1, part[level1][[i]], heltt[[1]], hel1qqg[[j]]] * Spaa[n1,f1]/Sqrt[mt2] /. f1->tmp /. n1->f1 /. tmp->n1;
camp[level1, part[level1][[i]], heltt[[3]], hel2qqg[[j]]] = camp[level1, part[level1][[i]], heltt[[1]], hel2qqg[[j]]] * Spaa[n1,f1]/Sqrt[mt2] /. f1->tmp /. n1->f1 /. tmp->n1;

camp[level1, part[level1][[i]], heltt[[4]], hel1qqg[[j]]] = camp[level1, part[level1][[i]], heltt[[1]], hel1qqg[[j]]] * Spaa[n1,f1]/Sqrt[mt2] * Spaa[n2,f2]/Sqrt[mt2]/. {f1->tmp1, n1->f1, f2->tmp2, n2->f2} /. {tmp1->n1, tmp2->n2};
camp[level1, part[level1][[i]], heltt[[4]], hel2qqg[[j]]] = camp[level1, part[level1][[i]], heltt[[1]], hel2qqg[[j]]] * Spaa[n1,f1]/Sqrt[mt2] * Spaa[n2,f2]/Sqrt[mt2]/. {f1->tmp1, n1->f1, f2->tmp2, n2->f2} /. {tmp1->n1, tmp2->n2};

,{i,{1,2,3,4}},{j,{1,2}}];

(* mt2dij[1] or [2] makes a difference *)
(* Print[camp[level1, part[[1]], heltt[[4]], hel1qqg[[1]]] /. mt2dij[[1]] /.dijvalues //N]
Print[camp[level1, part[[1]], heltt[[3]], hel1qqg[[1]]] * Spaa[n2,f2]/Sqrt[mt2] /. f2->tmp /. n2->f2 /. tmp->n2 /. mt2dij[[1]] /.dijvalues //N] *)


Amp2tree[level_, heltt_, helqqg_] := Rule[M[level, heltt, helqqg], Sum[Conjugate[camp[level, part["0L"][[l]], heltt, helqqg]]* Gij[l,m] * camp[level, part["0L"][[m]], heltt, helqqg],{l,1,4},{m,1,4}] ];

(*Amp2["0L", "--", "+--"]*)

(*Do[

Print[Amp2tree["0L", heltt[[i]], hel1qqg[[l]]]/. mt2dij[[1]] /.dijvalues //N]; 
,{i,{1,2,3,4}}
,{l,{1,2}}

];

Do[

Print[Amp2tree["0L", heltt[[i]], hel2qqg[[l]]]/. mt2dij[[1]] /.dijvalues //N]; 
,{i,{1,2,3,4}}
,{l,{1,2}}

];
*)

Amp2gp = Sum[Amp2tree["0L", heltt[[i]], hel1qqg[[j]]][[2]],{i,1,4},{j,1,2}]/. mt2dij[[1]] /.dijvalues //N
Amp2gm = Sum[Amp2tree["0L", heltt[[i]], hel2qqg[[j]]][[2]],{i,1,4},{j,1,2}]/. mt2dij[[1]] /.dijvalues //N
fullamp2tree = (Amp2gp + Amp2gm)
(* q,q_bar are incoming in openloops: 2 possible helicities, 3 possible colors --> 36 for avage*)


(*(* Loop contributions *)

(*part["1L"]={"d14T253","d23T451","d12T453","d34T251",  "Nfd14T253","Nfd23T451","Nfd12T453","Nfd34T251",  "Nhd14T253","Nhd23T451","Nhd12T453","Nhd34T251"};
level2 = "1L";
dsfactor["dsm2p0"] = 1.;
dsfactor["dsm2p1"] = cScheme[scheme];
ncfactor["Ncp0"] = 1.;
ncfactor["Ncp1"] = Nc;
ncfactor["Ncp-1"] = 1/Nc;
ncfactor["Ncp2"] = Nc^2;
ncfactor["Ncp-2"] = 1/Nc^2;*)

Do[
Do[
Do[
amp[level2, part[level2][[i]], helicity[[j]], ncpow, dspow] = Get["ttj_amplitudes/ttqqg/ttqqg_1L_"<>part[level2][[i]]<>"_"<>ncpow<>"_"<>helicity[[j]]<>"_"<>dspow<>".m"];
samp[level2, part[level2][[i]], helicity[[j]], ncpow, dspow] = ncfactor[ncpow]*dsfactor[dspow]*evalloop[amp[level2, part[level2][[i]], helicity[[j]], ncpow, dspow]];
If[StringStartsQ[part[level2][[i]], "Nf"],
	samp[level2, part[level2][[i]], helicity[[j]], ncpow, dspow] = 5.*samp[level2, part[level2][[i]], helicity[[j]], ncpow, dspow],
	samp[level2, part[level2][[i]], helicity[[j]], ncpow, dspow] = 1.*samp[level2, part[level2][[i]], helicity[[j]], ncpow, dspow]];	
,{dspow,dspowlist[part[level2][[i]]]}];

chnamp[level2, part[level2][[i]], helicity[[j]], ncpow] = Sum[samp[level2, part[level2][[i]], helicity[[j]], ncpow, ds],{ds, dspowlist[part[level2][[i]]]}];
,{ncpow,ncpowlist[part[level2][[i]]]}];

eamp[level2, part[level2][[i]], helicity[[j]]] = Sum[chnamp[level2, part[level2][[i]], helicity[[j]], nc], {nc, ncpowlist[part[level2][[i]]]}];

camp[level2, part[level2][[i]], heltt[[1]], hel1qqg[[j]]] = global * phase[hel1qqg[[j]]] * Sum[theta[l] * eamp[level2, part[level2][[i]], helicity[[j]]][[l]],{l, 1, 4}];
camp[level2, part[level2][[i]], heltt[[1]], hel2qqg[[j]]] = global * phase[hel2qqg[[j]]] * Sum[theta[l] * eamp[level2, part[level2][[i]], helicity[[j]]][[l]],{l, 1, 4}];

camp[level2, part[level2][[i]], heltt[[2]], hel1qqg[[j]]] = camp[level2, part[level2][[i]], heltt[[1]], hel1qqg[[j]]] * Spaa[n2,f2]/Sqrt[mt2] /. f2->tmp /. n2->f2 /. tmp->n2;
camp[level2, part[level2][[i]], heltt[[2]], hel2qqg[[j]]] = camp[level2, part[level2][[i]], heltt[[1]], hel2qqg[[j]]] * Spaa[n2,f2]/Sqrt[mt2] /. f2->tmp /. n2->f2 /. tmp->n2;

camp[level2, part[level2][[i]], heltt[[3]], hel1qqg[[j]]] = camp[level2, part[level2][[i]], heltt[[1]], hel1qqg[[j]]] * Spaa[n1,f1]/Sqrt[mt2] /. f1->tmp /. n1->f1 /. tmp->n1;
camp[level2, part[level2][[i]], heltt[[3]], hel2qqg[[j]]] = camp[level2, part[level2][[i]], heltt[[1]], hel2qqg[[j]]] * Spaa[n1,f1]/Sqrt[mt2] /. f1->tmp /. n1->f1 /. tmp->n1;

camp[level2, part[level2][[i]], heltt[[4]], hel1qqg[[j]]] = camp[level2, part[level2][[i]], heltt[[2]], hel1qqg[[j]]] * Spaa[n1,f1]/Sqrt[mt2] /. f1->tmp /. n1->f1 /. tmp->n1;
camp[level2, part[level2][[i]], heltt[[4]], hel2qqg[[j]]] = camp[level2, part[level2][[i]], heltt[[2]], hel2qqg[[j]]] * Spaa[n1,f1]/Sqrt[mt2] /. f1->tmp /. n1->f1 /. tmp->n1;

,{i,{1,2,3,4}},{j,{1,2}}];*)



(*Amp2loop[level_, heltt_, helqqg_] := Rule[M[level, heltt, helqqg], Normal@Series[Sum[Conjugate[camp[level, part["1L"][[l]], heltt, helqqg]]* Gijl[l,m] * camp[level, part["1L"][[m]], heltt, helqqg],{l,1,4},{m,1,4}],{eps,0,2}] /. mt2dij[[1]] /.dijvalues //N ];

fullamp2loop = Normal@Series[Assuming[Element[eps,Reals],Simplify[Sum[Amp2loop["1L", heltt[[i]], hel1qqg[[j]]][[2]],{i,1,4},{j,1,2}] + Sum[Amp2loop["1L", heltt[[i]], hel2qqg[[j]]][[2]],{i,1,4},{j,1,2}]]],{eps,0,2}]
*)



Quit
