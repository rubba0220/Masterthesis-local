(* ::Package:: *)

MomentumConservation = {mm[p5]->-mm[p1]-mm[p2]-mm[p3]-mm[p4]};
Replacements = {
  mmp2[p1] -> mt2, 
  mmp[p1, p2] -> d12, 
  mmp[p1, p3] -> d45 - d12 - d23 - mt2, 
  mmp[p1, p4] -> -d15 + d23 - d45, 
  mmp[p1, p5] -> d15, 
  mmp2[p2] -> mt2, 
  mmp[p2, p3] -> d23, 
  mmp[p2, p4] -> d15 - d23 - d34, 
  mmp[p2, p5] -> d34 - d12 - d15 - mt2, 
  mmp2[p3] -> 0, 
  mmp[p3, p4] -> d34, 
  mmp[p3, p5] -> d12 - d34 - d45 + mt2, 
  mmp[p4, p5] -> d45, 
  mmp2[p4] -> 0
};


Invariants = {
  d12 -> mmp[p1,p2],
  d23 -> mmp[p2,p3],
  d34 -> mmp[p3,p4],
  d45 -> mmp[p4,p5],
  d15 -> mmp[p1,p5],
  mt2 -> mmp2[p1]
}/.MomentumConservation//Expand;


sijreplace = {
 
 s12 -> 2*(d12+mt2),
 s13 -> 2*d45-2*d12-2*d23-mt2,
 s14 -> 2*d23-2*d45-2*d15+mt2,
 s15 -> 2*d15+mt2,
 s23 -> 2*d23+mt2,
 s24 -> 2*d15-2*d23-2*d34+mt2,
 s25 -> 2*d34-2*d15-2*d12-mt2,
 s34 -> 2*d34,
 s35 -> 2*d12-2*d34-2*d45+2*mt2,
 s45 -> 2*d45,
 
 d13 -> d45-d12-d23-mt2,
 d14 -> d23-d45-d15,
 d24 -> d15-d23-d34,
 d25 -> d34-d15-d12-mt2,
 d35 -> d12-d34-d45+mt2

};


permutations=<|
T1->{{1,2,3,4,5}, {1,2,4,3,5}, {1,2,4,5,3}, {1,2,5,4,3}, {1,2,3,5,4}, {1,2,5,3,4}},
T2->{{1,2,3,4,5}, {1,2,3,5,4}, {1,2,4,3,5}, {1,2,4,5,3}, {1,2,5,4,3}, {1,2,5,3,4}},
T3->{{1,4,2,5,3}, {1,5,2,4,3}, {1,4,2,3,5}, {1,5,2,3,4}, {1,3,2,4,5}, {1,3,2,5,4}},
T4->{{1,3,2,5,4}, {1,3,2,4,5}, {1,5,2,4,3}, {1,5,2,3,4}, {1,4,2,3,5}, {1,4,2,5,3}}
|>;

families = {};
topofamilies = {};

TVec={T1,T2,T3,T4};

FamilyPermutations=Table[permutations[i], {i,TVec}];

Masses={
  {0,mt2,0,0,0},
  {mt2,0,mt2,mt2,mt2},
  {0,mt2,mt2,0,0},
  {mt2,0,0,mt2,mt2}
};

masterTopos = {
  {p1, p2, p3, p4, p5},
  {p1, p2, p3, p4, p5}, 
  {p1, p3, p2, p4, p5},
  {p1, p3, p2, p4, p5}
};
MomentaPermutations = Table[
  Inner[Rule,masterTopos[[i]],{p1,p2,p3,p4,p5}[[FamilyPermutations[[i,ii]]]],List]
,{i, Length[TVec]},{ii,Length[FamilyPermutations[[i]]]}];
InvariantsPermutations = Table[
  Inner[Rule,{d12,d23,d34,d45,d15, mt2},{d12,d23,d34,d45,d15, mt2}/.Invariants/.MomentaPermutations[[i,ii]]/.MomentumConservation/.Replacements//Expand,List]
,{i, Length[TVec]},{ii,Length[FamilyPermutations[[i]]]}];
lenVc=Prepend[Table[Length[FamilyPermutations[[k]]], {k, Length[FamilyPermutations]}],0];
