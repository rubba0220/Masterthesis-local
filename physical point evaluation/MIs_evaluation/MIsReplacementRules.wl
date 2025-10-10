(* ::Package:: *)

$dir = If[$FrontEnd=!=Null,NotebookDirectory[],DirectoryName[$InputFileName]];
SetDirectory[$dir]


l1T2=Get["T2/Numerics/resT2_p12345.m"];
l2T2=Get["T2/Numerics/resT2_p12354.m"];
l3T2=Get["T2/Numerics/resT2_p12435.m"];
l4T2=Get["T2/Numerics/resT2_p12453.m"];
l5T2=Get["T2/Numerics/resT2_p12543.m"];
l6T2=Get["T2/Numerics/resT2_p12534.m"];
l1T1=Get["T1/Numerics/resT1_p12345.m"];
l2T1=Get["T1/Numerics/resT1_p12354.m"];
l3T1=Get["T1/Numerics/resT1_p12435.m"];
l4T1=Get["T1/Numerics/resT1_p12453.m"];
l5T1=Get["T1/Numerics/resT1_p12543.m"];
l6T1=Get["T1/Numerics/resT1_p12534.m"];
l1T3=Get["T3/Numerics/resT3_p13245.m"];
l2T3=Get["T3/Numerics/resT3_p15243.m"];
l3T3=Get["T3/Numerics/resT3_p14235.m"];
l4T3=Get["T3/Numerics/resT3_p15234.m"];
l5T3=Get["T3/Numerics/resT3_p14253.m"];
l6T3=Get["T3/Numerics/resT3_p13254.m"];
l1T4=Get["T4/Numerics/resT4_p13254.m"];
l2T4=Get["T4/Numerics//resT4_p13245.m"];
l3T4=Get["T4/Numerics//resT4_p15243.m"];
l4T4=Get["T4/Numerics//resT4_p15234.m"];
l5T4=Get["T4/Numerics//resT4_p14235.m"];
l6T4=Get["T4/Numerics//resT4_p14253.m"];


list=Join[l1T1,l1T2,l2T1,l2T2,l3T1,l3T2,l4T1,l4T2,l5T1,l5T2,l6T1,l6T2,l1T3,l2T3,l3T3,l4T3,l5T3,l6T3,l1T4,l2T4,l3T4,l4T4,l5T4,l6T4];


listNew=Join[{{d12->0.1030127934867479,d23->-0.07518278423934056,d34->0.5,d45->-0.319871570505173,d15->0.31832575235540994,mt2->0.029821836099999898}},{list}];


Export["T1T2T3T4_all_Permutations.m",listNew];
