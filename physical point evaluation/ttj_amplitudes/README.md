### pp -> tt+j helicity amplitudes

Analytic helicity amplitudes for pp->ttj are presented using the spin decomposition for the massive quarks detailed in the article:

There are two compressed tar archives containing the results for each channel. These must be unzipped before the checks in the base directory can be executed.

tar -xJf 

Explicitly we use:

AmpSpinProj[f1_,n1_,f2_,n2_][hel_] := mt*phase[hel]*(
spA[n1,n2]*s[3,4]/spA[f1,n1]/spA[f2,n2]*AA[1][hel] + 
spA[n1,p[3]]*spA[n2,p[4]]*spB[p[3],p[4]]/spA[f1,n1]/spA[f2,n2]*AA[2][hel] + 
spA[n1,p[3]]*spA[n2,p[3]]*spBB[p[3],p[4],p[5],p[3]]/s[3,4]/spA[f1,n1]/spA[f2,n2]*AA[3][hel] + 
spA[n1,p[4]]*spA[n2,p[4]]*spBB[p[4],p[5],p[3],p[4]]/s[3,4]/spA[f1,n1]/spA[f2,n2]*AA[4][hel]);

with phases for ttggg defined by 

phase["+++++"] = spB[p[3],p[5]]/(spA[p[3],p[4]]*spA[p[4],p[5]]);
phase["++++-"] = spAA[p[5],p[3],p[4],p[5]]/(spA[p[3],p[4]]^2);
phase["+++-+"] = spAA[p[4],p[5],p[3],p[4]]/(spA[p[5],p[3]]^2);
phase["++-++"] = spAA[p[3],p[4],p[5],p[3]]/(spA[p[4],p[5]]^2);

and phases for ttqqg defined by

phase["+++-+"] = spA[p[3],p[4]]/spA[p[3],p[5]]^2;
phase["++-++"] = spA[p[3],p[4]]/spA[p[4],p[5]]^2;

Each partial colour amplitude is given as a separte file labelled by the loop order, colour structure, helicity and component in ds-2.
ds=4-2*eps for the `t Hooft-Veltman scheme and ds=4 for the four dimensional helicity scheme.

filename format: <process>_<loop order>_<colour structure>_<helicity>_<ds-2 component>.m

"dsm2p0" = (ds-2)^0
"dsm2p1" = (ds-2)^1

Each file contains the following information:

{subamplitudes, coefficients, factors}

where subamplitudes contains four rules for the four subamplitudes AA[i] written in terms of linearly independent coefficients and master integrals. The coefficients are rational functions written in terms of a set of independent factors which are polynomials in eps and the momentum twistor variables {t12,t23,s34,t45,t51,x5123} defined in the article.


