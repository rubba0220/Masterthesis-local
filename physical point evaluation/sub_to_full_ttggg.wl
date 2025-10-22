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
Do[
Print["sqrt pi^2 = ", Sqrt[mdot[p,p]]];
,{p,{p1v,p2v,p3v,p4v,p5v}}];
Do[
Print["sqrt ni or fi^2 = ",Sqrt[mdot[p,p]]];
,{p,{n1v,n2v,f1v,f2v}}]*)


(*replacemi = integraleval[[2]];*)
(*dijvalues2 = integraleval[[1]];*)
dijvalues = {d12->mdot[p1v,p2v], d23->mdot[p2v,p3v], d34->mdot[p3v,p4v], d45->mdot[p4v,p5v], d15->mdot[p5v,p1v], mt2->mdot[p1v,p1v]
,d13->mdot[p1v,p3v], d14->mdot[p1v,p4v], d24->mdot[p2v,p4v], d25->mdot[p2v,p5v], d35->mdot[p3v,p5v]};
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
phase["+++"] = Spaa[p3,p5]/(Spaa[p3,p4]*Spaa[p4,p5]);
phase["++-"] = Spaa[p5,p3,p4,p5]/(Spaa[p3,p4]^2);
phase["+-+"] = Spaa[p4,p5,p3,p4]/(Spaa[p3,p5]^2);
phase["-++"] = phase["+-+"] /.p3->tmp/.p4->p3/.tmp->p4;

Clear[theta];
theta[1] = Spaa[n1,n2]*s34/Spaa[f1,n1]/Spaa[f2,n2];
theta[2] = Spaa[n1,p3]*Spaa[n2,p4]*Spbb[p3,p4]/Spaa[f1,n1]/Spaa[f2,n2];
theta[3] = Spaa[n1,p3]*Spaa[n2,p3]*Spbb[p3,p4,p5,p3]/s34/Spaa[f1,n1]/Spaa[f2,n2];
theta[4] = Spaa[n1,p4]*Spaa[n2,p4]*Spbb[p4,p5,p3,p4]/s34/Spaa[f1,n1]/Spaa[f2,n2];


Clear[CT];
CT[1,1][i1_, i2_, a3_, a4_, a5_] := Sum[T[a3][[i2,l]]T[a4][[l,m]]T[a5][[m,i1]],{l,1,3},{m,1,3}];
CT[1,2][i1_, i2_, a3_, a4_, a5_] := CT[1,1][i1, i2, a3, a5, a4]; (*4,5*)
CT[1,3][i1_, i2_, a3_, a4_, a5_] := CT[1,1][i1, i2, a4, a5, a3]; (*4,3 3,5*)
CT[1,4][i1_, i2_, a3_, a4_, a5_] := CT[1,1][i1, i2, a4, a3, a5]; (*3,4*)
CT[1,5][i1_, i2_, a3_, a4_, a5_] := CT[1,1][i1, i2, a5, a3, a4]; (*3,4 4,5*)
CT[1,6][i1_, i2_, a3_, a4_, a5_] := CT[1,1][i1, i2, a5, a4, a3]; (*3,5*)

CT[2,1][i1_, i2_, a3_, a4_, a5_] := KroneckerDelta[a3,a4]*T[a5][[i2, i1]];
CT[2,2][i1_, i2_, a3_, a4_, a5_] := CT[2,1][i1, i2, a4, a5, a3]; (*3,4 3,5*)
CT[2,3][i1_, i2_, a3_, a4_, a5_] := CT[2,1][i1, i2, a5, a3, a4]; (*3,4 4,5*)

CT[3,1][i1_, i2_, a3_, a4_, a5_] := KroneckerDelta[i2, i1]*Sum[T[a3][[l, m]]*T[a4][[m,n]]*T[a5][[n,l]], {m,1,3},{n,1,3},{l,1,3}];
CT[3,2][i1_, i2_, a3_, a4_, a5_] := CT[3,2][i1, i2, a3, a5, a4]; (*4,5*)

Gij[x1_,x2_,y1_,y2_] := Sum[Conjugate[CT[x1, x2][i1, i2, a3, a4, a5]]*CT[y1,y2][i1, i2, a3, a4, a5]
,{i1,1,3},{i2,1,3},{a3,1,8},{a4,1,8},{a5,1,8}];


GijTreematrix = {{0,0,0,0,0,0}, {0,0,0,0,0,0}, {0,0,0,0,0,0}, {0,0,0,0,0,0}, {0,0,0,0,0,0}, {0,0,0,0,0,0}};
Do[
GijTreematrix[[a]][[b]] = Gij[1,a,1,b];
,{a,{1,2,3,4,5,6}},{b,{1,2,3,4,5,6}}];
Print["Color metric = ",GijTreematrix]


{1,0,0,0,0,0} . (GijTreematrix . {1,1,1,1,1,1})


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
(*ttggg*)


(*Tree level*)
part["0L"]={"T23451"};
level1 = "0L";
helicity = {"+++++","++++-","+++-+","++-++"};
hel1qqg  = {"+++","++-","+-+","-++"};
hel2qqg = {"---","--+","-+-","+--"}; (*parity*)
heltt = {"++", "+-", "-+", "--"}; (*paper [30]*)

Do[
amp[level1, part[level1][[i]], helicity[[j]]] = Get["ttj_amplitudes/ttggg/tt3g_0L_"<>part[level1][[i]]<>"_"<>helicity[[j]]<>".m"];
eamp[level1, part[level1][[i]], helicity[[j]]] = evaltree[amp[level1, part[level1][[i]], helicity[[j]]]];

camp[level1, part[level1][[i]], heltt[[1]], hel1qqg[[j]]] = global * phase[hel1qqg[[j]]] * Sum[theta[l] * eamp[level1, part[level1][[i]], helicity[[j]]][[l,2]],{l, 1, 4}]; 
(*camp[level1, part[level1][[i]], heltt[[1]], hel2qqg[[j]]] = global * phase[hel2qqg[[j]]] * Sum[theta[l] * eamp[level1, part[level1][[i]], helicity[[j]]][[l,2]],{l, 1, 4}];*)

camp[level1, part[level1][[i]], heltt[[2]], hel1qqg[[j]]] = camp[level1, part[level1][[i]], heltt[[1]], hel1qqg[[j]]] * Spaa[n2,f2]/Sqrt[mt2] /. f2->tmp /. n2->f2 /. tmp->n2;
(*camp[level1, part[level1][[i]], heltt[[2]], hel2qqg[[j]]] = camp[level1, part[level1][[i]], heltt[[1]], hel2qqg[[j]]] * Spaa[n2,f2]/Sqrt[mt2] /. f2->tmp /. n2->f2 /. tmp->n2;*)

camp[level1, part[level1][[i]], heltt[[3]], hel1qqg[[j]]] = camp[level1, part[level1][[i]], heltt[[1]], hel1qqg[[j]]] * Spaa[n1,f1]/Sqrt[mt2] /. f1->tmp /. n1->f1 /. tmp->n1;
(*camp[level1, part[level1][[i]], heltt[[3]], hel2qqg[[j]]] = camp[level1, part[level1][[i]], heltt[[1]], hel2qqg[[j]]] * Spaa[n1,f1]/Sqrt[mt2] /. f1->tmp /. n1->f1 /. tmp->n1;*)

camp[level1, part[level1][[i]], heltt[[4]], hel1qqg[[j]]] = camp[level1, part[level1][[i]], heltt[[1]], hel1qqg[[j]]] * Spaa[n1,f1]/Sqrt[mt2] * Spaa[n2,f2]/Sqrt[mt2]/. {f1->tmp1, n1->f1, f2->tmp2, n2->f2} /. {tmp1->n1, tmp2->n2};
(*camp[level1, part[level1][[i]], heltt[[4]], hel2qqg[[j]]] = camp[level1, part[level1][[i]], heltt[[1]], hel2qqg[[j]]] * Spaa[n1,f1]/Sqrt[mt2] * Spaa[n2,f2]/Sqrt[mt2]/. {f1->tmp1, n1->f1, f2->tmp2, n2->f2} /. {tmp1->n1, tmp2->n2};*)

,{i,{1}},{j,{1,2,3,4}}];


rule[1,1] = {};
rule[1,2] = {d34->d35, d15->d14, p4->p5, p5->p4}; (*todo spinor adjustments*)
rule[1,3] = {d23->d24, d34->d45, d45->d35, d15->d13, p3->p4, p4->p5, p5->p3};
rule[1,4] = {d23->d24, d45->d35, p3->p4, p4->p3};
rule[1,5] = {d23->d25, d34->d35, d45->d34, d15->d14, p3->p5, p4->p3, p5->p4};
rule[1,6] = {d23->d25, d34->d45, d45->d34, d15->d13, p3->p5, p5->p3};

rule[2,1] = rule[1,1];
rule[2,2] = rule[1,3];
rule[2,3] = rule[1,5];

rule[3,1] = rule[1,1];
rule[3,2] = rule[1,2];

x2s[1] = 6;
x2s[2] = 3;
x2s[3] = 2;

Amp2tree[level_, heltt_, helqqg_] := Rule[M[level, heltt, helqqg], Sum[Conjugate[camp[level, part["0L"][[1]], heltt, helqqg]/.rule[1,x2]]* Gij[1,x2, 1,y2] * (camp[level, part["0L"][[1]], heltt, helqqg]/.rule[1,y2]),{x2, 1, x2s[1]},{y2,1,x2s[1]}]];

Amp2gp = Sum[Amp2tree["0L", heltt[[i]], hel1qqg[[j]]][[2]],{i,1,4},{j,1,2}]/. mt2dij[[1]] /.dijvalues //N
(*Amp2gm = Sum[Amp2tree["0L", heltt[[i]], hel2qqg[[j]]][[2]],{i,1,4},{j,1,2}]/. mt2dij[[1]] /.dijvalues //N
fullamp2tree = (Amp2gp + Amp2gm)*)
(* q,q_bar are incoming in openloops: 2 possible helicities, 3 possible colors --> 36 for avage*)



