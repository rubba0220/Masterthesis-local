### Ancillary files for "One-loop QCD helicity amplitudes for $pp\to t\tb j$ to $O(\eps^2)$"

This folder contains files and instructions for obtaining the boundary constant for the pentagon integral for T2 

There are 3 machine-readable files in this directory:
* Boundary0pentagonT2.wl - boundary constant for the pentagon integral at eps order 0
* Boundary1pentagonT2.wl - boundary constant for the pentagon integral at eps order 1
* NumbersPentagonT2.wl - the numerical value for this integral at eps order 0 and 1

Following are the instructions to obtain the number in the NumbersPentagonT2.wl file from the boundary files.

We perform the numerical evaluation using Ginac and Mathematica. A few special features of PolyLogTools are also useful.

1) Start by analyzing the singular points of the analytic expression.
This can be done for instance by naively analyzing the singular points of all the denominators:
       den=Cases[Collect[GatherTranscendentals[expression], _G, h@*Denominator@*Together], h[e_] :> e, Infinity] // DeleteDuplicates;
       list=Table[Solve[den[[i]] == 0, {x1}], {i, 1, Length[den]}]//Flatten // DeleteDuplicates;
list contains the points where the numerical evaluation has to be done extra carefully, pick out the points within the region of integration, i.e. within 0 and 1.

2) Since the expression is not integrated in the last variable, it contains integrations over G's, where G is the symbol for Goncharov Polylogarithms. We perform the numerical evaluation in two steps: first we evaluate the integrand at all the points (within 0 and 1) and then numerically evaluate them. One way to perform this is as follows:
 a)f[x_? NumericQ,p_]:= (expression/.G[a__,b_]:>Ginsh[G[a,b],{x1->x},PrecisionGoal->p])/.x1->x; (* Here Ginsh is a feature of PolyLogTools which invokes GiNac but this can be skipped completely *)

 b)NIntegrate[f[t,25],{t,0,1/2 (2-Sqrt[2]),1/2,3/4, 1- 10^(-20)},PrecisionGoal->20,WorkingPrecision->25]   

Here we have done the numerical integration choosing the boudaries of the regions in consideration of the list obtained in step 1. For more precision the PrecisionGoal and WorkingPrecision has to be increased accodingly. 


3) The numerical integration for order 0 in epsilon can be done in much simpler ways, for instance for precision 100 we can simply do

   NIntegrate[Ex//Re,{x1,0,1/2 (2-Sqrt[2]),1/2,3/4,1},PrecisionGoal->100,WorkingPrecision->100]


