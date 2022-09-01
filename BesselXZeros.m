(* :Title: BesselXZeros *)

(* :Author: Jeff Kershaw, National Institute of Radiological Sciences, QST, Chiba, Japan *)

(* :Context: BesselXZeros` *)

(* :Version: 2.0 *)

(* :Date: 2022-08-30 *)

(* :Description: This package is a Mathematica-based implementation of some techniques
    developed to find the complex z-roots of the Bessel cross-product functions for
    general complex order. The following manuscript contains more detail:

     Zeros distribution for cross-product combinations of the Bessel functions
     with complex order
     J. Kershaw, T. Obata (submitted to Const Approx, Sept 2022)
*)

BeginPackage["BesselXZeros`"];

BesselX00Zeros::usage=
"\!\(BesselX00Zeros[\*StyleBox[\"\[Nu]\",\"TI\"], \*StyleBox[\"\[Lambda]\",\"TI\"]]\) \
returns the \*StyleBox[\"z\",\"TI\"]-zeros of the Bessel cross-product function \
\!\(\*SubsuperscriptBox[X,StyleBox[\"\[Nu]\",\"TI\"],\"00\"]\)(z, \[Lambda]) \
for complex order \*StyleBox[\"\[Nu]\",\"TI\"] and non-negative real number \
\*StyleBox[\"\[Lambda]\",\"TI\"]."

BesselX01Zeros::usage=
"\!\(BesselX01Zeros[\*StyleBox[\"\[Nu]\",\"TI\"], \*StyleBox[\"\[Lambda]\",\"TI\"]]\) \
returns the \*StyleBox[\"z\",\"TI\"]-zeros of the Bessel cross-product function \
\!\(\*SubsuperscriptBox[X,StyleBox[\"\[Nu]\",\"TI\"],\"01\"]\)(z, \[Lambda]) \
for complex order \*StyleBox[\"\[Nu]\",\"TI\"] and non-negative real number \
\*StyleBox[\"\[Lambda]\",\"TI\"]."

BesselX10Zeros::usage=
"\!\(BesselX10Zeros[\*StyleBox[\"\[Nu]\",\"TI\"], \*StyleBox[\"\[Lambda]\",\"TI\"]]\) \
returns the \*StyleBox[\"z\",\"TI\"]-zeros of the Bessel cross-product function \
\!\(\*SubsuperscriptBox[X,StyleBox[\"\[Nu]\",\"TI\"],\"10\"]\)(z, \[Lambda]) \
for complex order \*StyleBox[\"\[Nu]\",\"TI\"] and non-negative real number \
\*StyleBox[\"\[Lambda]\",\"TI\"]."

BesselX11Zeros::usage=
"\!\(BesselX11Zeros[\*StyleBox[\"\[Nu]\",\"TI\"], \*StyleBox[\"\[Lambda]\",\"TI\"]]\) \
returns the \*StyleBox[\"z\",\"TI\"]-zeros of the Bessel cross-product function \
\!\(\*SubsuperscriptBox[X,StyleBox[\"\[Nu]\",\"TI\"],\"11\"]\)(z, \[Lambda]) \
for complex order \*StyleBox[\"\[Nu]\",\"TI\"] and non-negative real number \
\*StyleBox[\"\[Lambda]\",\"TI\"]."

BesselX00::usage=
"\!\(BesselX00[\*StyleBox[\"z\",\"TI\"], \*StyleBox[\"\[Nu]\",\"TI\"], \*StyleBox[\"\[Lambda]\",\"TI\"]]\) \
gives the Bessel cross-product function \
\!\(\*SubsuperscriptBox[X,StyleBox[\"\[Nu]\",\"TI\"],\"00\"]\)(z, \[Lambda]) \
= \!\(\*SubscriptBox[J,StyleBox[\"\[Nu]\",\"TI\"]](z) \
\*SubscriptBox[Y,StyleBox[\"\[Nu]\",\"TI\"]](\*StyleBox[\"\[Lambda] z\"]) \
- \*SubscriptBox[Y,StyleBox[\"\[Nu]\",\"TI\"]](z) \
\*SubscriptBox[J,StyleBox[\"\[Nu]\",\"TI\"]](\*StyleBox[\"\[Lambda] z\"]) \
for complex order \*StyleBox[\"\[Nu]\",\"TI\"] and non-negative real number \
\*StyleBox[\"\[Lambda]\",\"TI\"]\)."

BesselX01::usage=
"\!\(BesselX01[\*StyleBox[\"z\",\"TI\"], \*StyleBox[\"\[Nu]\",\"TI\"], \*StyleBox[\"\[Lambda]\",\"TI\"]]\) \
gives the Bessel cross-product function \
\!\(\*SubsuperscriptBox[X,StyleBox[\"\[Nu]\",\"TI\"],\"01\"]\)(z, \[Lambda]) \
= \!\(\*SubscriptBox[J,StyleBox[\"\[Nu]\",\"TI\"]](z) \
\*SubscriptBox[Y,StyleBox[\"\[Nu]\",\"TI\"]]'(\*StyleBox[\"\[Lambda] z\"]) \
- \*SubscriptBox[Y,StyleBox[\"\[Nu]\",\"TI\"]](z) \
\*SubscriptBox[J,StyleBox[\"\[Nu]\",\"TI\"]]'(\*StyleBox[\"\[Lambda] z\"]) \
for complex order \*StyleBox[\"\[Nu]\",\"TI\"] and non-negative real number \
\*StyleBox[\"\[Lambda]\",\"TI\"]\)."

BesselX10::usage=
"\!\(BesselX10[\*StyleBox[\"z\",\"TI\"], \*StyleBox[\"\[Nu]\",\"TI\"], \*StyleBox[\"\[Lambda]\",\"TI\"]]\) \
gives the Bessel cross-product function \
\!\(\*SubsuperscriptBox[X,StyleBox[\"\[Nu]\",\"TI\"],\"10\"]\)(z, \[Lambda]) \
= \!\(\*SubscriptBox[J,StyleBox[\"\[Nu]\",\"TI\"]]'(z) \
\*SubscriptBox[Y,StyleBox[\"\[Nu]\",\"TI\"]](\*StyleBox[\"\[Lambda] z\"]) \
- \*SubscriptBox[Y,StyleBox[\"\[Nu]\",\"TI\"]]'(z) \
\*SubscriptBox[J,StyleBox[\"\[Nu]\",\"TI\"]](\*StyleBox[\"\[Lambda] z\"]) \
for complex order \*StyleBox[\"\[Nu]\",\"TI\"] and non-negative real number \
\*StyleBox[\"\[Lambda]\",\"TI\"]\)."

BesselX11::usage=
"\!\(BesselX11[\*StyleBox[\"z\",\"TI\"], \*StyleBox[\"\[Nu]\",\"TI\"], \*StyleBox[\"\[Lambda]\",\"TI\"]]\) \
gives the Bessel cross-product function \
\!\(\*SubsuperscriptBox[X,StyleBox[\"\[Nu]\",\"TI\"],\"11\"]\)(z, \[Lambda]) \
= \!\(\*SubscriptBox[J,StyleBox[\"\[Nu]\",\"TI\"]]'(z) \
\*SubscriptBox[Y,StyleBox[\"\[Nu]\",\"TI\"]]'(\*StyleBox[\"\[Lambda] z\"]) \
- \*SubscriptBox[Y,StyleBox[\"\[Nu]\",\"TI\"]]'(z) \
\*SubscriptBox[J,StyleBox[\"\[Nu]\",\"TI\"]]'(\*StyleBox[\"\[Lambda] z\"]) \
for complex order \*StyleBox[\"\[Nu]\",\"TI\"] and non-negative real number \
\*StyleBox[\"\[Lambda]\",\"TI\"]\)."

Unprotect[BesselX00Zeros];
Unprotect[BesselX01Zeros];
Unprotect[BesselX10Zeros];
Unprotect[BesselX11Zeros];
Unprotect[BesselX00];
Unprotect[BesselX01];
Unprotect[BesselX10];
Unprotect[BesselX11];


Begin["`Private`"]

(* ::Subsubsection:: *)
(*Definitions independent of a & b*)

\[Sigma][z_]:=Module[{argz=Arg[z],absz=Abs[z],sigp=ArcSech[z]},
  Which[argz==0&&absz>1,Conjugate[sigp],
  True,sigp]
]

\[Rho][z_]:=\[Sigma][z]-Tanh[\[Sigma][z]]

zi[z_]:=Block[{rho,zp,zip,th=Arg[z],abszip,argrho},
  If[0<th<=Pi,zp=z,zp=Conjugate[z]];
  rho=\[Rho][zp];
  If[Re[zp]>1&&Im[zp]==0,
   zip=(3/2Abs[rho])^(2/3)Exp[I \[Pi]],
   zip=(3/2rho)^(2/3);
   abszip=Abs[zip];
   argrho=Arg[rho]-If[-Pi<=Arg[zip]<=0,0,2Pi];
   zip=abszip Exp[I 2/3 argrho];
  ];
  If[0<th<=Pi,zip,Conjugate[zip]]
]

\[Eta][Z_,\[Lambda]_]:=\[Rho][Z]-\[Rho][\[Lambda] Z]

rho[Z_]:=Log[(1+Sqrt[1-Z^2])/Z]-Sqrt[1-Z^2]

zed[zi_]:=Module[{zip,rho,sig,s},
  If[Arg[zi]<=0,zip=zi,zip=Conjugate[zi]];
  Which[zi==0,sig=0,zi==0.,sig=0.,
   Re[zi]<0&&Im[zi]==0,rho=(2/3)Abs[zi]^(3/2)I;
   sig=(Re@s I)/.FindRoot[s I -Tanh[s I]-rho,{s,-\[Pi]/4,-\[Pi]/2,0}],
   True,rho=(2/3)zip^(3/2);sig=s/.FindRoot[s-Tanh[s]-rho,{s,10.^-10-I Pi/2}]
  ];
  Which[Im[zi]==0,Re@Sech[sig],Arg[zi]<=0,Sech[sig],True,Conjugate[Sech[sig]]]
]

ptO\[Alpha][\[Alpha]_,\[Lambda]_]:=Block[{t\[Alpha],\[Alpha]p,zO,t\[Alpha]0},
  If[\[Alpha]>0,\[Alpha]p=-\[Alpha],\[Alpha]p=\[Alpha]];
  If[Log10[-\[Alpha]p]<-5,If[Abs[\[Alpha]p]<Abs[-10^-5-\[Alpha]p],\[Alpha]p=0,\[Alpha]p=-10^-5]];
  t\[Alpha]0=If[\[Alpha]p<-0.1,10^-5.,-zi[\[Lambda]]/2];
  zO=Switch[\[Alpha]p,0,1,0.,1,-\[Pi]/2,1/\[Lambda],-\[Pi]/2.,1/\[Lambda],
    _,(zed[Re[t\[Alpha]/.#] Exp[I \[Pi]]/Exp[I \[Alpha]p]^(2/3)]/\[Lambda])&@
      FindRoot[Arg[Exp[2I \[Alpha]p/3]zi[zed[t\[Alpha] Exp[\[Pi] I]/Exp[2I \[Alpha]p/3]]/\[Lambda]]]+\[Pi]/3,
       {t\[Alpha],t\[Alpha]0,0,-zi[\[Lambda]]},Evaluated->False]];
  If[\[Alpha]>0,Conjugate[zO],zO]
]/;ReQ[\[Lambda]]&&(1<\[Lambda])&&(Abs[\[Alpha]]<=\[Pi]/2)

ReQ[x_]:=NumericQ[x]&&(Im[x]==0)

CmplxQ[x_]:=NumericQ[x]&&(Im[x]!=0)


(* ::Subsubsection:: *)
(*Definitions for X11*)

Xi11[Z_,\[Nu]_,\[Lambda]_]:=AiryAiPrime[\[Nu]^(2/3)zi[Z]]AiryBiPrime[\[Nu]^(2/3)zi[\[Lambda] Z]]-
  AiryBiPrime[\[Nu]^(2/3)zi[Z]]AiryAiPrime[\[Nu]^(2/3)zi[\[Lambda] Z]]

Xi11Z0m1[m1_(*Integer*),\[Nu]_(*Complex*),\[Lambda]_]:=
  zed[(3/2Abs[m1-3/4]\[Pi]/Abs[\[Nu]])^(2/3)Exp[(-\[Pi]+2Abs[Arg[\[Nu]]]/3)I]]/\[Lambda]/;
   ((*(1\[LessEqual]m1\[LessEqual]M11m1[\[Nu],\[Lambda]])&&*)ReQ[\[Lambda]]&&\[Lambda]>1)

Xi11Z0m2[m2_(*Integer*),\[Nu]_(*Complex*),\[Lambda]_]:=
  zed[(3/2Abs[m2-3/4]\[Pi]/Abs[\[Nu]])^(2/3)Exp[(-\[Pi]/3+2Abs[Arg[\[Nu]]]/3)I]]/;
   ((*(1\[LessEqual]m2\[LessEqual]M11m2[\[Nu],\[Lambda]])&&*)ReQ[\[Lambda]]&&\[Lambda]>1)

Xi11Z0m3[m3_(*Integer*), \[Nu]_(*Complex*),\[Lambda]_]:=
  Block[{Z,\[Nu]p,Zp},
   If[Arg[\[Nu]]<=0,\[Nu]p=\[Nu],\[Nu]p=Conjugate[\[Nu]]];
   Z=Zp/.FindRoot[\[Nu]p(\[Rho][Zp]-\[Rho][\[Lambda] Zp])+m3 \[Pi] I,{Zp,1+I},Evaluated->False];
   If[Arg[\[Nu]]<=0,Z,Conjugate[Z]]
  ]/;((*(M11m3[\[Nu],\[Lambda]]\[LessEqual]m3)&&*)CmplxQ[\[Nu]]&&ReQ[\[Lambda]]&&\[Lambda]>1)

Xi11Z0m3[m3_(*Integer*), \[Nu]_(*Real*),\[Lambda]_]:=
  Block[{Z},
   Re@Z/.FindRoot[\[Nu](\[Rho][Z]-\[Rho][\[Lambda] Z])+m3 \[Pi] I,{Z,1},Evaluated->False]
  ]/;((*(M11m3[\[Nu],\[Lambda]]\[LessEqual]m3)&&*)ReQ[\[Nu]]&&ReQ[\[Lambda]]&&\[Lambda]>1)

M11m1[\[Nu]_,\[Lambda]_]:=
  (Floor@Re@(-I # \[Rho][\[Lambda] ptO\[Alpha][Arg[#],\[Lambda]]]/\[Pi]+3/4))&@If[Arg[\[Nu]]>0,Conjugate[\[Nu]],\[Nu]]

M11m2[\[Nu]_,\[Lambda]_]:=
  (Floor@Re@(I # \[Rho][ptO\[Alpha][Arg[#],\[Lambda]]]/\[Pi]+3/4))&@If[Arg[\[Nu]]>0,Conjugate[\[Nu]],\[Nu]]

M11m3[\[Nu]_,\[Lambda]_]:=
  Block[{\[Nu]p=If[Arg[\[Nu]]>0,Conjugate[\[Nu]],\[Nu]]},
   ((Ceiling@Re@(I \[Nu]p (\[Rho][#]-\[Rho][\[Lambda] #])/\[Pi]))&@ptO\[Alpha][Arg[\[Nu]p],\[Lambda]])]

Xi11Zm1[m1_,\[Nu]_,\[Lambda]_,N_Integer:1,showAllZn_Symbol:False]:=
  Block[{z,zz,z0rule,znrules,p,prules,n,lhs,eqns,soln,\[Epsilon],N\[Epsilon]=If[N==0,0,Floor[(N-1)/2]],z0=Xi11Z0m1[m1,\[Nu],\[Lambda]]},
   If[N==0,Return[z0],z0rule=z[0]->z0];
   prules=Join[{p[z[0]]->\[Rho][z0],p[\[Lambda] z[0]]->\[Rho][\[Lambda] z0]},
    Table[Derivative[n][p][z[0]]->Derivative[n][rho][z0],{n,1,N}],
     Table[Derivative[n][p][\[Lambda] z[0]]->Derivative[n][rho][\[Lambda] z0],{n,1,N}]];
   zz=Sum[\[Epsilon]^n z[n],{n,0,N}]+O[\[Epsilon]]^(N+1);
   lhs=Series[2Exp[I \[Pi]/4+\[Nu] p[zz]]Cos[\[Pi]/4-I \[Nu] p[\[Lambda] zz]]-\[Epsilon] Exp[-\[Nu](p[zz]-p[\[Lambda] zz])],{\[Epsilon],0,N}];
   eqns=(Map[#==0&,CoefficientList[lhs,\[Epsilon]]][[2;;N+1]])/.prules//Simplify;
   soln=Solve[eqns,Table[z[n],{n,1,N}]];
   znrules=Flatten[Append[soln,{z0rule}]];
   If[N\[Epsilon]==0,
    If[showAllZn,
     {Normal[zz/.znrules]/.{\[Epsilon]->1},Table[z[n],{n,0,N}]/.znrules},
     Normal[zz/.znrules]/.{\[Epsilon]->1}],
    If[showAllZn,
     {PadeApproximant[(zz/.znrules),{\[Epsilon],0,N\[Epsilon]}]/.{\[Epsilon]->1},Table[z[n],{n,0,N}]/.znrules},
     PadeApproximant[(zz/.znrules),{\[Epsilon],0,N\[Epsilon]}]/.{\[Epsilon]->1}]
   ]
  ]/;((*(1\[LessEqual]m1\[LessEqual]M11m1[\[Nu],\[Lambda]])&&*)N>=0&&ReQ[\[Lambda]]&&\[Lambda]>1)

Xi11Zm2[m2_,\[Nu]_,\[Lambda]_,N_Integer:1,showAllZn_Symbol:False]:=
  Block[{z,zz,z0rule,znrules,p,prules,n,lhs,eqns,soln,\[Epsilon],N\[Epsilon]=If[N==0,0,Floor[(N-1)/2]],z0=Xi11Z0m2[m2,\[Nu],\[Lambda]]},
   If[N==0,Return[z0],z0rule=z[0]->z0];
   prules=Join[{p[z[0]]->\[Rho][z0],p[\[Lambda] z[0]]->\[Rho][\[Lambda] z0]},
    Table[Derivative[n][p][z[0]]->Derivative[n][rho][z0],{n,1,N}],
     Table[Derivative[n][p][\[Lambda] z[0]]->Derivative[n][rho][\[Lambda] z0],{n,1,N}]];
   zz=Sum[\[Epsilon]^n z[n],{n,0,N}]+O[\[Epsilon]]^(N+1);
   lhs=Series[2Exp[I 3\[Pi]/4+\[Nu] p[\[Lambda] zz]]Cos[\[Pi]/4+I \[Nu] p[zz]]+\[Epsilon] Exp[\[Nu](p[zz]-p[\[Lambda] zz])],{\[Epsilon],0,N}];
   eqns=(Map[#==0&,CoefficientList[lhs,\[Epsilon]]][[2;;N+1]])/.prules//Simplify;
   soln=Solve[eqns,Table[z[n],{n,1,N}]];
   znrules=Flatten[Append[soln,{z0rule}]];
   If[N\[Epsilon]==0,
    If[showAllZn,
     {Normal[zz/.znrules]/.{\[Epsilon]->1},Table[z[n],{n,0,N}]/.znrules},
     Normal[zz/.znrules]/.{\[Epsilon]->1}],
    If[showAllZn,
     {PadeApproximant[(zz/.znrules),{\[Epsilon],0,N\[Epsilon]}]/.{\[Epsilon]->1},Table[z[n],{n,0,N}]/.znrules},
     PadeApproximant[(zz/.znrules),{\[Epsilon],0,N\[Epsilon]}]/.{\[Epsilon]->1}]
   ]
  ]/;((*(1\[LessEqual]m2\[LessEqual]M11m2[\[Nu],\[Lambda]])&&*) N>=0&&ReQ[\[Lambda]]&&\[Lambda]>1)

Xi11Zm3[m3_,\[Nu]_,\[Lambda]_,N_Integer:1,showAllZn_Symbol:False]:=
  Block[{z,zz,z0rule,znrules,p,prules,n,lhs,eqns,soln,\[Epsilon],N\[Epsilon]=If[N==0,0,Floor[(N-1)/2]],z0=Xi11Z0m3[m3,\[Nu],\[Lambda]]},
   If[N==0,Return[z0],z0rule=z[0]->z0];
   prules=Join[{p[z[0]]->\[Rho][z0],p[\[Lambda] z[0]]->\[Rho][\[Lambda] z0]},
    Table[Derivative[n][p][z[0]]->Derivative[n][rho][z0],{n,1,N}],
     Table[Derivative[n][p][\[Lambda] z[0]]->Derivative[n][rho][\[Lambda] z0],{n,1,N}]];
   zz=Sum[\[Epsilon]^n z[n],{n,0,N}]+O[\[Epsilon]]^(N+1);
   lhs=Series[-2I Sin[I \[Nu](p[zz]- p[\[Lambda] zz])]+\[Epsilon] I Exp[\[Nu](p[zz]+p[\[Lambda] zz])],{\[Epsilon],0,N}];
   eqns=(Map[#==0&,CoefficientList[lhs,\[Epsilon]]][[2;;N+1]])/.prules//Simplify;
   soln=Solve[eqns,Table[z[n],{n,1,N}]];
   znrules=Flatten[Append[soln,{z0rule}]];
   If[N\[Epsilon]==0,
    If[showAllZn,
     {Normal[zz/.znrules]/.{\[Epsilon]->1},Table[z[n],{n,0,N}]/.znrules},
     Normal[zz/.znrules]/.{\[Epsilon]->1}],
    If[showAllZn,
     {PadeApproximant[(zz/.znrules),{\[Epsilon],0,N\[Epsilon]}]/.{\[Epsilon]->1},Table[z[n],{n,0,N}]/.znrules},
     PadeApproximant[(zz/.znrules),{\[Epsilon],0,N\[Epsilon]}]/.{\[Epsilon]->1}]
   ]
  ]/;((*(M11m3[\[Nu],\[Lambda]]\[LessEqual]m3)&&*)N>=0&&ReQ[\[Lambda]]&&\[Lambda]>1)

findRootXi11m1[seeds_,\[Nu]_,\[Lambda]_]:=
  Block[{Z,Z0,Z1,Z2,R,root,i,j,Nseeds=Length[seeds],m1,mm,debug=False},
   Table[
    Z0=seeds[[i]];
    root=Z/.FindRoot[N[Xi11[SetPrecision[Z,\[Infinity]],\[Nu],\[Lambda]],50],{Z,SetPrecision[Z0,\[Infinity]]},Evaluated->False];
    If[i==1,Z1=1/\[Lambda],Z1=seeds[[i-1]]];
    If[i==Nseeds,Z2=Xi11Z0m1[Nseeds+1,\[Nu],\[Lambda]],Z2=seeds[[i+1]]];
    R=Min@Abs[Z0-{Z1,Z2}];
    If[debug,Print[i,"    ",{Z1,Z0,Z2},"     ",R];Print[root]];
    j=0;m1=i;mm=m1;
    While[(Abs[Z0-root]>=0.95R)&&mm>If[m1==1,3/4,m1-1],
     j=j+1;mm=m1-j*.001;Z0=Xi11Zm1[mm,\[Nu],\[Lambda],1];
     root=Z/.FindRoot[N[Xi11[SetPrecision[Z,\[Infinity]],\[Nu],\[Lambda]],50],{Z,SetPrecision[Z0,\[Infinity]]},Evaluated->False,MaxIterations->1000];
     If[debug,Print["   ",j,"   ",mm,"    ",Z0,"    ",root]]
    ];
    root,{i,1,Nseeds}]
   ]

findRootXi11m2[seeds_,\[Nu]_,\[Lambda]_,m0_,roots1_,roots3_]:=
  Block[{Z,Z0,Z1,Z2,R,root,i,j,Nseeds=Length[seeds],M1=Length[roots1],m2,mm,sa,test,roots13,r1,r3,debug=False},
   sa[Z_,a_:10]:=SetAccuracy[Z,a];
   Which[
    (m0==-1)||(m0==0),
     roots13=Join[roots1,roots3];
     test[rts_,rt_]:=ContainsAny[rts,{rt},SameTest->((sa[#]==sa[#2])&)];,
    m0==1,
     If[M1==0,r1=Null,r1=Last[roots1]];r3=First[roots3];
     If[TrueQ[sa[r1]==sa[r3]],
      roots13=Join[roots1,Rest[roots3]];
      test[rts_,rt_]:=ContainsAny[rts,{rt},SameTest->((sa[#]==sa[#2])&)];,
      roots13=Join[Most[roots1],Rest[roots3]];
      test[rts_,rt_]:=(ContainsAny[rts,{rt},SameTest->((sa[#]==sa[#2])&)]|| !(TrueQ[sa[rt]==sa[r1]]||(sa[rt]==sa[r3])));
   ];];
   Table[
    Z0=seeds[[i]];
    root=Z/.FindRoot[N[Xi11[SetPrecision[Z,\[Infinity]],\[Nu],\[Lambda]],50],{Z,SetPrecision[Z0,\[Infinity]]},Evaluated->False];
    If[i==1,Z1=1+0.I,Z1=seeds[[i-1]]];
    If[i==Nseeds,Z2=ptO\[Alpha][Arg[\[Nu]],\[Lambda]],Z2=seeds[[i+1]]];
    R=Min@Abs[Z0-{Z1,Z2}];
    If[debug,Print[i,"    ",{Z2,Z0,Z1},"     ",R];Print[root]];
     j=0;m2=i;mm=m2;
     While[(Abs[Z0-root]>=0.95R)&&mm>If[m2==1,3/4,m2-1]&& test[roots13,root],
      j=j+1;mm=m2-j*.001;Z0=Xi11Zm2[mm,\[Nu],\[Lambda],1];
      root=Z/.FindRoot[N[Xi11[SetPrecision[Z,\[Infinity]],\[Nu],\[Lambda]],50],{Z,SetPrecision[Z0,\[Infinity]]},Evaluated->False,MaxIterations->1000];
      If[debug,Print["   ",j,"   ",mm,"    ",Z0,"    ",root]];
     ];
     root,{i,1,Nseeds}]
  ]

Xi11Roots[\[Nu]_(*Complex*),\[Lambda]_(*Real*),Nz_Integer:10]:=
  Block[{\[Nu]p,\[Lambda]p,seeds1,seeds2,seeds3,roots1,roots2,roots3,Z,m1,m2,m3,m0,M1,M2,M3,s1,s2,s3,r1,r2,r3,dd,pp,roots,sa,sortFail=False,O\[Alpha],s0,r0,eq7bits,debug=False},
   sa[Z_,a_:10]:=SetAccuracy[Z,a];
   \[Nu]p=If[Arg[\[Nu]]>0,Conjugate[\[Nu]],\[Nu]];
   \[Lambda]p=If[0<\[Lambda]<1,1/\[Lambda],\[Lambda]];
   If[debug,Print["Xi11Roots:0"];Print["\[Nu]: ",\[Nu]," \[Equivalent] ",\[Nu]p,", \[Lambda]: ",\[Lambda]," \[Equivalent] ",\[Lambda]p]];
   M1=M11m1[\[Nu]p,\[Lambda]p];
   seeds1=Table[Xi11Zm1[m1,\[Nu]p,\[Lambda]p,1],{m1,1,M1}];
   roots1=findRootXi11m1[Chop@seeds1,\[Nu]p,\[Lambda]p];
   If[M1!=#,Print["Unexpected number of roots on seedline 1. Expected: ",M1,", Found: ",#,"."]; Print[{\[Nu],\[Lambda]}];
    Print[seeds1]; Print[roots1]; Abort[];]&@(Length@Union[roots1,SameTest->((sa[#]==sa[#2])&)]);
   If[debug,Print["Xi11Roots:1"];Print[seeds1];Print[roots1]];
   M2=M11m2[\[Nu]p,\[Lambda]p];
   seeds2=Table[Xi11Zm2[m2,\[Nu]p,\[Lambda]p,5],{m2,1,M2}];
   M3=M11m3[\[Nu]p,\[Lambda]p];
   If[debug,Print["Xi11Roots:3.0"];Print[{M3,M3+(Nz-1),Nz}]];
   If[Abs@Im[\[Nu]p]<0.1,
    seeds3=Table[Xi11Z0m3[m3,\[Nu]p,\[Lambda]p],{m3,M3,M3+(Nz-1)}],
    seeds3=Table[Xi11Zm3[m3,\[Nu]p,\[Lambda]p,1],{m3,M3,M3+(Nz-1)}]];
   roots3=Z/.FindRoot[N[Xi11[SetPrecision[Z,\[Infinity]],\[Nu]p,\[Lambda]p],50],{Z,SetPrecision[#,\[Infinity]]},Evaluated->False]&/@seeds3;
   If[Nz!=#,Print["Unexpected number of roots on seedline 3. Expected: ",Nz,", Found: ",#,"."]; Print[{\[Nu],\[Lambda]}];
    Print[seeds3]; Print[roots3]; Abort[];]&@(Length@Union[roots3,SameTest->((sa[#]==sa[#2])&)]);
   If[debug,Print["Xi11Roots:3.1"];Print[seeds3];Print[roots3]];
   m0=M1+M2-M3;If[debug,Print[{M1,M2,M3}," m0 = ",m0]];
   roots2=findRootXi11m2[seeds2,\[Nu]p,\[Lambda]p,m0,roots1,roots3];
   If[M2!=#,Print["Unexpected number of roots on seedline 2. Expected: ",M2,", Found: ",#,"."];  Print[{\[Nu],\[Lambda]}];
    Print[seeds2]; Print[roots2]; Abort[];]&@ (Length@Union[roots2,SameTest->((sa[#]==sa[#2])&)]);
   If[debug,Print["Xi11Roots:2"];Print[seeds2];Print[roots2]];
   roots=Union[roots1,roots2,roots3,SameTest->((sa[#]==sa[#2])&)];
   Switch[m0,
    1, (* one extra seed *)
    If[M1+M2+Nz==Length[roots]+1,
     If[M1==0,s1=Null;r1=Null,s1=Last[seeds1];r1=Last[roots1]];
     If[M2==0,s2=Null;r2=Null,s2=Last[seeds2];r2=Last[roots2]];
     s3=First[seeds3];r3=First[roots3];
     dd=Abs@{r1-s1,r2-s2,r3-s3};
     Which[
      TrueQ[sa[r1]==sa[r2]],pp=Position[#,Max@#]&@ReplacePart[dd,3->0],
      TrueQ[sa[r1]==sa[r3]],pp=Position[#,Max@#]&@ReplacePart[dd,2->0],
      TrueQ[sa[r2]==sa[r3]],pp=Position[#,Max@#]&@ReplacePart[dd,1->0],
      True,Print["Unexpected sorting problem (m0 = ",m0,")."]; Print[{s1,s2,s3}]; Print[{r1,r2,r3}];Print[dd];Abort[];];
     If[debug,Print[{s1,s2,s3}];Print[{r1,r2,r3}];Print[dd];Print["pp = ",pp]];
     Switch[First@@pp,
      1,roots1=Most@roots1,
      2,roots2=Most@roots2 ,
      3,roots3=Rest@roots3
     ],
     Print["Unexpected number of roots. Expected: ",M1+M2+Nz-1,", Found: ",Length[roots]," (m0 = ",m0,")"];
     If[(M1+1)+(M2-1)+Nz==Length[roots],(* could be a roots = seeds case? *)
     Print["Could be a misclassified m0 = 2 case."]];sortFail=True;
    ],
    0, (* #roots == #seeds *)
    If[M1+M2+Nz!= Length[roots],
     Print["Unexpected number of roots. Expected: ",M1+M2+Nz,", Found: ",Length[roots]," (m0 = ",m0,")"];sortFail=True;],
    -1,  (* one extra root *)
    If[M1+M2+Nz==Length[roots],
     (* Manufacture seed for extra root by expanding Eq.(7) about O\[Alpha] *)
     If[debug,Print["Xi11Roots:4.0"]; Print["seeds1 = ",seeds1]; Print["roots1 = ",roots1]; Print["seeds2 = ",seeds2];
      Print["roots2 = ",roots2]; Print["seeds3 = ",seeds3]; Print["roots3 = ",roots3];];
     eq7bits[Z_,nu_,lambda_]:={I Exp[nu (rho[Z]+rho[lambda Z])],Exp[nu (rho[Z]-rho[lambda Z])],-Exp[nu (-rho[Z]+rho[lambda Z])]};
     O\[Alpha]=ptO\[Alpha][Arg[\[Nu]p],\[Lambda]p];If[debug,Print["O\[Alpha] = ",O\[Alpha]]];
     s0=(Z/.Solve[Normal[Quiet@Series[Plus@@eq7bits[Z,\[Nu]p,\[Lambda]p],{Z,O\[Alpha],1},Assumptions->Im[Z]>0]]==0,Z][[1]]);
     r0=(Z/.FindRoot[Xi11[Z,\[Nu]p,\[Lambda]p],{Z,s0},Evaluated->False]);
     If[debug,Print["Xi11Roots:4.1"];Print[{N@s0,r0,Abs[r0-s0]}];
      Print[{#,Abs@#,Position[Abs@#,Min@Abs@#]}&@(N@eq7bits[s0,\[Nu]p,\[Lambda]p])];
       Print[{#,Abs@#,Position[Abs@#,Min@Abs@#]}&@(N@eq7bits[r0,\[Nu]p,\[Lambda]p])];];
     If[Length[roots]==Length@Union[{r0},roots,SameTest->((sa[#]==sa[#2])&)],Print["Extra seed did not converge to central root."];
      Print[{\[Nu],\[Lambda]}]; Print[s0]; Print[r0];  Print[roots]; Abort[];,
      (* else sort central root into one of the seedlines *)
      If[Re[\[Nu]p]==0,s1=Null;,s1=If[Im[\[Nu]p]==0,Xi11Z0m1[M1+1,\[Nu]p,\[Lambda]p],Xi11Zm1[M1+1,\[Nu]p,\[Lambda]p]]];
      If[Im[\[Nu]p]==0,s2=Null;,s2=Xi11Zm2[M2+1,\[Nu]p,\[Lambda]p]];
      If[(M3-1)==0,s3=Null;,s3=If[Im[\[Nu]p]==0,Xi11Z0m3[M3-1,\[Nu]p,\[Lambda]p],Xi11Zm3[M3-1,\[Nu]p,\[Lambda]p,1]]];
      dd=Abs[r0-{s1,s2,s3}];
      If[debug,Print[dd]; Print[{s1,s2,s3}] ;Print[{TrueQ[(s1==Null)&&(s3==Null)],TrueQ[(s2==Null)&&(s3==Null)],TrueQ[s1==Null],TrueQ[s2==Null],TrueQ[s3==Null]}]];
      Which[
       TrueQ[(s1==Null)&&(s3==Null)],pp=2,
       TrueQ[(s2==Null)&&(s3==Null)],pp=1,
       TrueQ[s1==Null],pp=Position[#,Min@#]&@ReplacePart[dd,1->10^10],
       TrueQ[s2==Null],pp=Position[#,Min@#]&@ReplacePart[dd,2->10^10],
       TrueQ[s3==Null],pp=Position[#,Min@#]&@ReplacePart[dd,3->10^10],
       True,pp=Position[dd,Min@dd]
      ];
      If[debug,Print[{s1,s2,s3}];Print[dd];Print["pp = ",pp]];
      Switch[First@@pp,
       1,roots1=Append[roots1,r0],
       2,roots2=Append[roots2,r0],
       3,roots3=Prepend[roots3,r0]
      ]
     ],
     Print["Unexpected number of roots. Expected: ",M1+M2+Nz,", Found: ",Length[roots]," (m0 = ",m0,")"]; sortFail=True;
    ],
    _,Print["Unexpected value for m0: ",{m0,{\[Nu]p,\[Lambda]p,M1,M2,M3}}, ". Aborting."];Abort[];
   ];
   If[sortFail,Print["Check whether the seeds have converged to the expected root. (m0 = ",m0,")"]; Print[{\[Nu]p,\[Lambda]p}];
    Print["seeds1 = ",seeds1]; Print["roots1 = ",roots1]; Print["seeds2 = ",seeds2]; Print["roots2 = ",roots2];
     Print["seeds3 = ",seeds3]; Print["roots3 = ",roots3]; Print["filtered roots = ",roots]; Abort[];];
   If[debug,Print["Xi11Roots:5"]];
   roots={roots1,roots2,roots3};
   Which[Arg[\[Nu]]==0,roots=Re@roots,Arg[\[Nu]]>0,roots=Conjugate[roots],True,roots];
   If[0<\[Lambda]<1,roots=roots/\[Lambda],roots]
  ]/;ReQ[\[Lambda]]&&(\[Lambda]>1)

BesselX11[z_,\[Nu]_,\[Lambda]_]:=
  ((BesselY[\[Nu]-1,\[Lambda] z]-BesselY[\[Nu]+1,\[Lambda] z])(BesselJ[\[Nu]-1,z]-BesselJ[\[Nu]+1,z])
   -(BesselJ[\[Nu]-1,\[Lambda] z]-BesselJ[\[Nu]+1,\[Lambda] z])(BesselY[\[Nu]-1,z]-BesselY[\[Nu]+1,z]))/4/;
    ((*ReQ[\[Nu]]&&*)ReQ[\[Lambda]](*&&\[Lambda]\[NotEqual]1*))

BesselX11Zeros[\[Nu]_,\[Lambda]_,Nz_:10]:=
  Block[{\[Nu]p,\[Lambda]p,seeds,roots,z,i,debug=False},
   \[Lambda]p=SetPrecision[If[0<\[Lambda]<1,1/\[Lambda],\[Lambda]],\[Infinity]];
   \[Nu]p=SetPrecision[#,\[Infinity]]&@Which[-\[Pi]<Arg[\[Nu]]<-\[Pi]/2,Conjugate[-\[Nu]],\[Pi]/2<Arg[\[Nu]]<=\[Pi],-\[Nu],
    0<Arg[\[Nu]]<=\[Pi]/2,Conjugate[\[Nu]],-\[Pi]/2<=Arg[\[Nu]]<=0,\[Nu],True,Print["Unexpected value for Arg[\[Nu]]."];Abort[]];
   If[debug,Print["BesselX11Zeros:1"];Print[\[Nu]," \[Equivalent] ",\[Nu]p,", ",\[Lambda]," \[Equivalent] ",\[Lambda]p,", ",Nz]];
   seeds=Xi11Roots[\[Nu]p,\[Lambda]p,Nz];
   If[debug,Print["BesselX11Zeros:2"];Print[seeds];Print[\[Nu]p seeds]];
   roots=Table[z/.FindRoot[N[BesselX11[SetPrecision[z,\[Infinity]],\[Nu]p,\[Lambda]p],50],{z,SetPrecision[#,\[Infinity]]},
    Evaluated->False]&/@(\[Nu]p seeds[[i]]),{i,1,3}];
   If[debug,Print["BesselX11Zeros:3"];Print[roots]];
   Which[Arg[\[Nu]]==0,roots=Re@roots,-\[Pi]<Arg[\[Nu]]<-\[Pi]/2,roots=-Conjugate[roots],\[Pi]/2<Arg[\[Nu]]<=\[Pi],roots=-roots,
    0<Arg[\[Nu]]<=\[Pi]/2,roots=Conjugate[roots],-\[Pi]/2<=Arg[\[Nu]]<=0,roots];
   If[0<\[Lambda]<1,roots=roots/\[Lambda],roots]
  ]/;(ReQ[\[Lambda]]&&\[Lambda]>0 &&\[Lambda]!=1)


(* ::Subsubsection:: *)
(*Definitions for X00*)

Xi00[Z_,\[Nu]_,\[Lambda]_]:=AiryAi[\[Nu]^(2/3)zi[Z]]AiryBi[\[Nu]^(2/3)zi[\[Lambda] Z]]-
  AiryBi[\[Nu]^(2/3)zi[Z]]AiryAi[\[Nu]^(2/3)zi[\[Lambda] Z]]

Xi00Z0m1[m1_(*Integer*),\[Nu]_(*Complex*),\[Lambda]_Real]:=
  zed[(3/2Abs[m1-1/4]\[Pi]/Abs[\[Nu]])^(2/3)Exp[(-\[Pi]+2Abs[Arg[\[Nu]]]/3)I]]/\[Lambda]/;
   ((*(1\[LessEqual]m1\[LessEqual]M00m1[\[Nu],\[Lambda]]+1)&&*)\[Lambda]>1)

Xi00Z0m2[m2_(*Integer*),\[Nu]_(*Complex*),\[Lambda]_Real]:=
  zed[(3/2Abs[m2-1/4]\[Pi]/Abs[\[Nu]])^(2/3)Exp[(-\[Pi]/3+2Abs[Arg[\[Nu]]]/3)I]]/;
   ((*(1\[LessEqual]m2\[LessEqual]M00m2[\[Nu],\[Lambda]])&&*)\[Lambda]>1)

Xi00Z0m3[m3_(*Integer*), \[Nu]_Complex,\[Lambda]_Real]:=
Block[{Z,\[Nu]p,Zp},
 If[Arg[\[Nu]]<=0,\[Nu]p=\[Nu],\[Nu]p=Conjugate[\[Nu]]];
 Z=Zp/.FindRoot[\[Nu]p(\[Rho][Zp]-\[Rho][\[Lambda] Zp])+m3 \[Pi] I,{Zp,1+I},Evaluated->False];
 If[Arg[\[Nu]]<=0,Z,Conjugate[Z]]
]/;((*(M00m3[\[Nu],\[Lambda]]\[LessEqual]m3)&&*)\[Lambda]>1)

Xi00Z0m3[m3_(*Integer*), \[Nu]_Real,\[Lambda]_Real]:=
Block[{Z},
 Re@Z/.FindRoot[\[Nu](\[Rho][Z]-\[Rho][\[Lambda] Z])+m3 \[Pi] I,{Z,1},Evaluated->False]
]/;((*(M00m3[\[Nu],\[Lambda]]\[LessEqual]m3)&&*)\[Lambda]>1)

M00m1[\[Nu]_,\[Lambda]_]:=
 (Floor@Re@(-I # \[Rho][\[Lambda] ptO\[Alpha][Arg[#],\[Lambda]]]/\[Pi]+1/4))&@If[Arg[\[Nu]]>0,Conjugate[\[Nu]],\[Nu]]

M00m2[\[Nu]_,\[Lambda]_]:=
 (Floor@Re@(I # \[Rho][ptO\[Alpha][Arg[#],\[Lambda]]]/\[Pi]+1/4))&@If[Arg[\[Nu]]>0,Conjugate[\[Nu]],\[Nu]]

M00m3[\[Nu]_,\[Lambda]_]:=
Block[{\[Nu]p=If[Arg[\[Nu]]>0,Conjugate[\[Nu]],\[Nu]]},
 ((Ceiling@Re@(I \[Nu]p (\[Rho][#]-\[Rho][\[Lambda] #])/\[Pi]))&@ptO\[Alpha][Arg[\[Nu]p],\[Lambda]])
]

Xi00Zm1[m1_,\[Nu]_,\[Lambda]_Real,N_Integer:1,showAllZn_Symbol:False]:=
Block[{z,zz,z0rule,znrules,p,prules,n,lhs,eqns,soln,\[Epsilon],N\[Epsilon]=If[N==0,0,Floor[(N-1)/2]],z0=Xi00Z0m1[m1,\[Nu],\[Lambda]]},
 If[N==0,Return[z0],z0rule=z[0]->z0];
 prules=Join[{p[z[0]]->\[Rho][z0],p[\[Lambda] z[0]]->\[Rho][\[Lambda] z0]},
  Table[Derivative[n][p][z[0]]->Derivative[n][rho][z0],{n,1,N}],
   Table[Derivative[n][p][\[Lambda] z[0]]->Derivative[n][rho][\[Lambda] z0],{n,1,N}]];
 zz=Sum[\[Epsilon]^n z[n],{n,0,N}]+O[\[Epsilon]]^(N+1);
 lhs=Series[2I Exp[I \[Pi]/4+\[Nu] p[zz]]Sin[\[Pi]/4-I \[Nu] p[\[Lambda] zz]]+\[Epsilon] Exp[-\[Nu](p[zz]-p[\[Lambda] zz])],{\[Epsilon],0,N}];
 eqns=(Map[#==0&,CoefficientList[lhs,\[Epsilon]]][[2;;N+1]])/.prules//Simplify;
 soln=Solve[eqns,Table[z[n],{n,1,N}]];
 znrules=Flatten[Append[soln,{z0rule}]];
 If[N\[Epsilon]==0,
  If[showAllZn,
   {Normal[zz/.znrules]/.{\[Epsilon]->1},Table[z[n],{n,0,N}]/.znrules},
   Normal[zz/.znrules]/.{\[Epsilon]->1}],
  If[showAllZn,
   {PadeApproximant[(zz/.znrules),{\[Epsilon],0,N\[Epsilon]}]/.{\[Epsilon]->1},Table[z[n],{n,0,N}]/.znrules},
   PadeApproximant[(zz/.znrules),{\[Epsilon],0,N\[Epsilon]}]/.{\[Epsilon]->1}]
 ]
]/;((*(1\[LessEqual]m1\[LessEqual]M00m1[\[Nu],\[Lambda]])&&*)N>=0&&\[Lambda]>1)

Xi00Zm2[m2_,\[Nu]_,\[Lambda]_Real,N_Integer:1,showAllZn_Symbol:False]:=
Block[{z,zz,z0rule,znrules,p,prules,n,lhs,eqns,soln,\[Epsilon],N\[Epsilon]=If[N==0,0,Floor[(N-1)/2]],z0=Xi00Z0m2[m2,\[Nu],\[Lambda]]},
 If[N==0,Return[z0],z0rule=z[0]->z0];
 prules=Join[{p[z[0]]->\[Rho][z0],p[\[Lambda] z[0]]->\[Rho][\[Lambda] z0]},
  Table[Derivative[n][p][z[0]]->Derivative[n][rho][z0],{n,1,N}],
   Table[Derivative[n][p][\[Lambda] z[0]]->Derivative[n][rho][\[Lambda] z0],{n,1,N}]];
 zz=Sum[\[Epsilon]^n z[n],{n,0,N}]+O[\[Epsilon]]^(N+1);
 lhs=Series[2Exp[I \[Pi]/4+\[Nu] p[\[Lambda] zz]]Cos[\[Pi]/4-I \[Nu] p[zz]]-\[Epsilon] Exp[\[Nu](p[zz]-p[\[Lambda] zz])],{\[Epsilon],0,N}];
 eqns=(Map[#==0&,CoefficientList[lhs,\[Epsilon]]][[2;;N+1]])/.prules//Simplify;
 soln=Solve[eqns,Table[z[n],{n,1,N}]];
 znrules=Flatten[Append[soln,{z0rule}]];
 If[N\[Epsilon]==0,
  If[showAllZn,
   {Normal[zz/.znrules]/.{\[Epsilon]->1},Table[z[n],{n,0,N}]/.znrules},
   Normal[zz/.znrules]/.{\[Epsilon]->1}],
  If[showAllZn,
   {PadeApproximant[(zz/.znrules),{\[Epsilon],0,N\[Epsilon]}]/.{\[Epsilon]->1},Table[z[n],{n,0,N}]/.znrules},
   PadeApproximant[(zz/.znrules),{\[Epsilon],0,N\[Epsilon]}]/.{\[Epsilon]->1}]
 ]
]/;((*(1\[LessEqual]m2\[LessEqual]M00m2[\[Nu],\[Lambda]])&&*) N>=0&&\[Lambda]>1)

Xi00Zm3[m3_,\[Nu]_,\[Lambda]_Real,N_Integer:1,showAllZn_Symbol:False]:=
Block[{z,zz,z0rule,znrules,p,prules,n,lhs,eqns,soln,\[Epsilon],N\[Epsilon]=If[N==0,0,Floor[(N-1)/2]],z0=Xi00Z0m3[m3,\[Nu],\[Lambda]]},
 If[N==0,Return[z0],z0rule=z[0]->z0];
 prules=Join[{p[z[0]]->\[Rho][z0],p[\[Lambda] z[0]]->\[Rho][\[Lambda] z0]},
  Table[Derivative[n][p][z[0]]->Derivative[n][rho][z0],{n,1,N}],
   Table[Derivative[n][p][\[Lambda] z[0]]->Derivative[n][rho][\[Lambda] z0],{n,1,N}]];
 zz=Sum[\[Epsilon]^n z[n],{n,0,N}]+O[\[Epsilon]]^(N+1);
 lhs=Series[-2I Sin[-I \[Nu](p[zz]- p[\[Lambda] zz])]+\[Epsilon] I Exp[\[Nu](p[zz]+p[\[Lambda] zz])],{\[Epsilon],0,N}];
 eqns=(Map[#==0&,CoefficientList[lhs,\[Epsilon]]][[2;;N+1]])/.prules//Simplify;
 soln=Solve[eqns,Table[z[n],{n,1,N}]];
 znrules=Flatten[Append[soln,{z0rule}]];
 If[N\[Epsilon]==0,
  If[showAllZn,
   {Normal[zz/.znrules]/.{\[Epsilon]->1},Table[z[n],{n,0,N}]/.znrules},
   Normal[zz/.znrules]/.{\[Epsilon]->1}],
  If[showAllZn,
   {PadeApproximant[(zz/.znrules),{\[Epsilon],0,N\[Epsilon]}]/.{\[Epsilon]->1},Table[z[n],{n,0,N}]/.znrules},
   PadeApproximant[(zz/.znrules),{\[Epsilon],0,N\[Epsilon]}]/.{\[Epsilon]->1}]
 ]
]/;((*(M00m3[\[Nu],\[Lambda]]\[LessEqual]m3)&&*)N>=0&&\[Lambda]>1)

findRootXi00m1[seeds_,\[Nu]_,\[Lambda]_]:=
Block[{Z,Z0,Z1,Z2,R,root,i,j,Nseeds=Length[seeds],m1,mm,debug=False},
 Table[
  Z0=seeds[[i]];
  root=Z/.FindRoot[Xi00[Z,\[Nu],\[Lambda]],{Z,Z0},Evaluated->False];
  If[i==1,Z1=1/\[Lambda],Z1=seeds[[i-1]]];
  If[i==Nseeds,Z2=Xi00Z0m1[Nseeds+1,\[Nu],\[Lambda]],Z2=seeds[[i+1]]];
  R=Min@Abs[Z0-{Z1,Z2}];
  If[debug,Print[i,"    ",{Z1,Z0,Z2},"     ",R];Print[root]];
  j=0;m1=i;mm=m1;
  While[(Abs[Z0-root]>=0.95R)&&mm>If[m1==1,1/4,m1-1],
   j=j+1;mm=m1-j*.001;Z0=Xi00Zm1[mm,\[Nu],\[Lambda],1];
   root=Z/.FindRoot[Xi00[Z,\[Nu],\[Lambda]],{Z,Z0},Evaluated->False,MaxIterations->1000];
   If[debug,Print["   ",j,"   ",mm,"    ",Z0,"    ",root]]
  ];
  root,{i,1,Nseeds}]
]

findRootXi00m2[seeds_,\[Nu]_,\[Lambda]_,m0_,roots1_,roots3_]:=
Block[{Z,Z0,Z1,Z2,R,root,i,j,Nseeds=Length[seeds],M1=Length[roots1],m2,mm,sa,test,roots13,r1,r3,debug=False},
 sa[Z_,a_:10]:=SetAccuracy[Z,a];
 Which[
  (m0==-1)||(m0==0),
   roots13=Join[roots1,roots3];test[rts_,rt_]:=ContainsAny[rts,{rt},SameTest->((sa[#]==sa[#2])&)];,
  m0==1,
   If[M1==0,r1=Null,r1=Last[roots1]];r3=First[roots3];
   If[TrueQ[sa[r1]==sa[r3]],
   roots13=Join[roots1,Rest[roots3]];
   test[rts_,rt_]:=ContainsAny[rts,{rt},SameTest->((sa[#]==sa[#2])&)];,
   roots13=Join[Most[roots1],Rest[roots3]];
   test[rts_,rt_]:=(ContainsAny[rts,{rt},SameTest->((sa[#]==sa[#2])&)]|| !(TrueQ[sa[rt]==sa[r1]]||(sa[rt]==sa[r3])));
 ];];
 Table[
  Z0=seeds[[i]];
  root=Z/.FindRoot[Xi00[Z,\[Nu],\[Lambda]],{Z,Z0},Evaluated->False];
  If[i==1,Z1=1+0.I,Z1=seeds[[i-1]]];
  If[i==Nseeds,Z2=ptO\[Alpha][Arg[\[Nu]],\[Lambda]],Z2=seeds[[i+1]]];
  R=Min@Abs[Z0-{Z1,Z2}];
  If[debug,Print[i,"    ",{Z2,Z0,Z1},"     ",R];Print[root]];
  j=0;m2=i;mm=m2;
  While[(Abs[Z0-root]>=0.95R)&&mm>If[m2==1,1/4,m2-1]&& test[roots13,root],
   j=j+1;mm=m2-j*.001;Z0=Xi00Zm2[mm,\[Nu],\[Lambda],1];
   root=Z/.FindRoot[Xi00[Z,\[Nu],\[Lambda]],{Z,Z0},Evaluated->False,MaxIterations->1000];
   If[debug,Print["   ",j,"   ",mm,"    ",Z0,"    ",root]];
  ];
  root,{i,1,Nseeds}]
]

Xi00Roots[\[Nu]_(*Complex*),\[Lambda]_(*Real*),Nz_Integer:10]:=
Block[{\[Nu]p,\[Lambda]p,seeds1,seeds2,seeds3,roots1,roots2,roots3,Z,m1,m2,m3,m0,M1,M2,M3,
      s1,s2,s3,r1,r2,r3,dd,pp,roots,sa,sortFail=False,O\[Alpha],s0,r0,eq22bits,debug=False},
 sa[Z_,a_:10]:=SetAccuracy[Z,a];
 \[Nu]p=If[Arg[\[Nu]]>0,Conjugate[\[Nu]],\[Nu]];
 \[Lambda]p=If[0<\[Lambda]<1,1/\[Lambda],\[Lambda]];
 If[debug,Print["Xi00Roots:0"];Print["\[Nu]: ",\[Nu]," \[Equivalent] ",\[Nu]p,", \[Lambda]: ",\[Lambda]," \[Equivalent] ",\[Lambda]p]];
 M1=M00m1[\[Nu]p,\[Lambda]p];
 seeds1=Table[Xi00Zm1[m1,\[Nu]p,\[Lambda]p,1],{m1,1,M1}];
 roots1=findRootXi00m1[Chop@seeds1,\[Nu]p,\[Lambda]p];
 If[M1!=#,Print["Unexpected number of roots on seedline 1. Expected: ",M1,", Found: ",#,"."]; Print[{\[Nu],\[Lambda]}];
  Print[seeds1]; Print[roots1]; Abort[];]&@(Length@Union[roots1,SameTest->((sa[#]==sa[#2])&)]);
 If[debug,Print["Xi00Roots:1"];Print[seeds1];Print[roots1]];
 M2=M00m2[\[Nu]p,\[Lambda]p];
 seeds2=Table[Xi00Zm2[m2,\[Nu]p,\[Lambda]p,5],{m2,1,M2}];
 M3=M00m3[\[Nu]p,\[Lambda]p];
 If[debug,Print["Xi00Roots:3.0"];Print[{M3,M3+(Nz-1),Nz}]];
 If[Abs@Im[\[Nu]p]<0.1,
  seeds3=Table[Xi00Z0m3[m3,\[Nu]p,\[Lambda]p],{m3,M3,M3+(Nz-1)}],
  seeds3=Table[Xi00Zm3[m3,\[Nu]p,\[Lambda]p,1],{m3,M3,M3+(Nz-1)}]];
 roots3=Z/.FindRoot[Xi00[Z,\[Nu]p,\[Lambda]p],{Z,#},Evaluated->False]&/@seeds3;
 If[Nz!=#,Print["Unexpected number of roots on seedline 3. Expected: ",Nz,", Found: ",#,"."]; Print[{\[Nu],\[Lambda]}];
  Print[seeds3]; Print[roots3]; Abort[];]&@(Length@Union[roots3,SameTest->((sa[#]==sa[#2])&)]);
 If[debug,Print["Xi00Roots:3.1"];Print[seeds3];Print[roots3]];
 m0=M1+M2-M3+1;If[debug,Print[{M1,M2,M3}," m0 = ",m0]];
 roots2=findRootXi00m2[seeds2,\[Nu]p,\[Lambda]p,m0,roots1,roots3];
 If[M2!=#,Print["Unexpected number of roots on seedline 2. Expected: ",M2,", Found: ",#,"."];  Print[{\[Nu],\[Lambda]}];
  Print[seeds2]; Print[roots2]; Abort[];]&@ (Length@Union[roots2,SameTest->((sa[#]==sa[#2])&)]);
 If[debug,Print["Xi00Roots:2"];Print[seeds2];Print[roots2]];
 roots=Union[roots1,roots2,roots3,SameTest->((sa[#]==sa[#2])&)];
 Switch[m0,
 1, (* one extra seed *)
 If[M1+M2+Nz==Length[roots]+1,
  If[M1==0,s1=Null;r1=Null,s1=Last[seeds1];r1=Last[roots1]];
  If[M2==0,s2=Null;r2=Null,s2=Last[seeds2];r2=Last[roots2]];
  s3=First[seeds3];r3=First[roots3];
  dd=Abs@{r1-s1,r2-s2,r3-s3};
  Which[
   TrueQ[sa[r1]==sa[r2]],pp=Position[#,Max@#]&@ReplacePart[dd,3->0],
   TrueQ[sa[r1]==sa[r3]],pp=Position[#,Max@#]&@ReplacePart[dd,2->0],
   TrueQ[sa[r2]==sa[r3]],pp=Position[#,Max@#]&@ReplacePart[dd,1->0],
   True,Print["Unexpected sorting problem (m0 = ",m0,")."]; Print[{s1,s2,s3}]; Print[{r1,r2,r3}];Print[dd];Abort[];];
  If[debug,Print[{s1,s2,s3}];Print[{r1,r2,r3}];Print[dd];Print["pp = ",pp]];
  Switch[First@@pp,
   1,roots1=Most@roots1,
   2,roots2=Most@roots2 ,
   3,roots3=Rest@roots3
  ],
  Print["Unexpected number of roots. Expected: ",M1+M2+Nz-1,", Found: ",Length[roots]," (m0 = ",m0,")"];
  If[(M1+1)+(M2-1)+Nz==Length[roots],(* could be a roots = seeds case? *)
   Print["Could be a misclassified m0 = 2 case."]];sortFail=True;
 ],
 0, (* #roots == #seeds *)
 If[M1+M2+Nz!= Length[roots],
  Print["Unexpected number of roots. Expected: ",M1+M2+Nz,", Found: ",Length[roots]," (m0 = ",m0,")"];sortFail=True;],
 -1,  (* one extra root *)
 If[M1+M2+Nz==Length[roots],
  (* Manufacture seed for extra root by expanding Eq.(7) about O\[Alpha] *)
  If[debug,Print["Xi00Roots:4.0"]; Print["seeds1 = ",seeds1]; Print["roots1 = ",roots1]; Print["seeds2 = ",seeds2];
   Print["roots2 = ",roots2]; Print["seeds3 = ",seeds3]; Print["roots3 = ",roots3];];
  eq22bits[Z_,nu_,lambda_]:={I Exp[nu (rho[Z]+rho[lambda Z])],-Exp[nu (rho[Z]-rho[lambda Z])],Exp[nu (-rho[Z]+rho[lambda Z])]};
  O\[Alpha]=ptO\[Alpha][Arg[\[Nu]p],\[Lambda]p];If[debug,Print["O\[Alpha] = ",O\[Alpha]]];
  s0=(Z/.Solve[Normal[Quiet@Series[Plus@@eq22bits[Z,\[Nu]p,\[Lambda]p],{Z,O\[Alpha],1},Assumptions->Im[Z]>0]]==0,Z][[1]]);
  r0=(Z/.FindRoot[Xi00[Z,\[Nu]p,\[Lambda]p],{Z,s0},Evaluated->False]);
  If[debug,Print["Xi00Roots:4.1"];Print[{N@s0,r0,Abs[r0-s0]}];
   Print[{#,Abs@#,Position[Abs@#,Min@Abs@#]}&@(N@eq22bits[s0,\[Nu]p,\[Lambda]p])];
    Print[{#,Abs@#,Position[Abs@#,Min@Abs@#]}&@(N@eq22bits[r0,\[Nu]p,\[Lambda]p])];];
  If[Length[roots]==Length@Union[{r0},roots,SameTest->((sa[#]==sa[#2])&)],Print["Extra seed did not converge to central root."];
   Print[{\[Nu],\[Lambda]}]; Print[s0]; Print[r0];  Print[roots]; Abort[];,
   (* else sort central root into one of the seedlines *)
   If[Re[\[Nu]p]==0,s1=Null;,s1=If[Im[\[Nu]p]==0,Xi00Z0m1[M1+1,\[Nu]p,\[Lambda]p],Xi00Zm1[M1+1,\[Nu]p,\[Lambda]p]]];
   If[Im[\[Nu]p]==0,s2=Null;,s2=Xi00Zm2[M2+1,\[Nu]p,\[Lambda]p]];
   If[(M3-1)==0,s3=Null;,s3=If[Im[\[Nu]p]==0,Xi00Z0m3[M3-1,\[Nu]p,\[Lambda]p],Xi00Zm3[M3-1,\[Nu]p,\[Lambda]p,1]]];
   dd=Abs[r0-{s1,s2,s3}];
   If[debug,Print[dd]; Print[{s1,s2,s3}] ;Print[{TrueQ[(s1==Null)&&(s3==Null)],TrueQ[(s2==Null)&&(s3==Null)],TrueQ[s1==Null],TrueQ[s2==Null],TrueQ[s3==Null]}]];
   Which[
    TrueQ[(s1==Null)&&(s3==Null)],pp=2,
    TrueQ[(s2==Null)&&(s3==Null)],pp=1,
    TrueQ[s1==Null],pp=Position[#,Min@#]&@ReplacePart[dd,1->10^10],
    TrueQ[s2==Null],pp=Position[#,Min@#]&@ReplacePart[dd,2->10^10],
    TrueQ[s3==Null],pp=Position[#,Min@#]&@ReplacePart[dd,3->10^10],
    True,pp=Position[dd,Min@dd]
   ];
   If[debug,Print[{s1,s2,s3}];Print[dd];Print["pp = ",pp]];
   Switch[First@@pp,
    1,roots1=Append[roots1,r0],
    2,roots2=Append[roots2,r0],
    3,roots3=Prepend[roots3,r0]
   ]
  ],
  Print["Unexpected number of roots. Expected: ",M1+M2+Nz,", Found: ",Length[roots]," (m0 = ",m0,")"]; sortFail=True;
 ],
 _,Print["Unexpected value for m0: ",{m0,{\[Nu]p,\[Lambda]p,M1,M2,M3}}, ". Aborting."];Abort[];
 ];
 If[sortFail,Print["Check whether the seeds have converged to the expected root. (m0 = ",m0,")"]; Print[{\[Nu]p,\[Lambda]p}];
  Print["seeds1 = ",seeds1]; Print["roots1 = ",roots1]; Print["seeds2 = ",seeds2]; Print["roots2 = ",roots2];
   Print["seeds3 = ",seeds3]; Print["roots3 = ",roots3]; Print["filtered roots = ",roots]; Abort[];];
 If[debug,Print["Xi00Roots:5"]];
 roots={roots1,roots2,roots3};
 Which[Arg[\[Nu]]==0,roots=Re@roots,Arg[\[Nu]]>0,roots=Conjugate[roots],True,roots];
 If[0<\[Lambda]<1,roots=roots/\[Lambda],roots]
]/;ReQ[\[Lambda]]&&(\[Lambda]>1)

BesselX00[z_,\[Nu]_,\[Lambda]_]:=BesselJ[\[Nu],z]BesselY[\[Nu],\[Lambda] z]-BesselY[\[Nu],z]BesselJ[\[Nu],\[Lambda] z]

BesselX00Zeros[\[Nu]_,\[Lambda]_,Nz_:10]:=
Block[{\[Nu]p,\[Lambda]p,seeds,roots,z,i,debug=False},
 \[Lambda]p=If[0<\[Lambda]<1,1/\[Lambda],\[Lambda]];
 \[Nu]p=Which[-\[Pi]<Arg[\[Nu]]<-\[Pi]/2,Conjugate[-\[Nu]],\[Pi]/2<Arg[\[Nu]]<=\[Pi],-\[Nu],0<Arg[\[Nu]]<=\[Pi]/2,Conjugate[\[Nu]],
  -\[Pi]/2<=Arg[\[Nu]]<=0,\[Nu],True,Print["Unexpected value for Arg[\[Nu]]."];Abort[]];
 If[debug,Print["BesselX00Zeros:1"];Print[\[Nu]," \[Equivalent] ",\[Nu]p,", ",\[Lambda]," \[Equivalent] ",\[Lambda]p,", ",Nz]];
 seeds=\[Nu]p Xi00Roots[\[Nu]p,\[Lambda]p,Nz];
 If[debug,Print["BesselX00Zeros:2"];Print[seeds];Print[\[Nu]p seeds]];
 roots=Table[z/.FindRoot[BesselX00[z,\[Nu]p,\[Lambda]p],{z,#},Evaluated->False]&/@seeds[[i]],{i,1,3}];
 If[debug,Print["BesselX00Zeros:3"];Print[roots]];
 Which[Arg[\[Nu]]==0,roots=Re@roots,-\[Pi]<Arg[\[Nu]]<-\[Pi]/2,roots=-Conjugate[roots],\[Pi]/2<Arg[\[Nu]]<=\[Pi],roots=-roots,
  0<Arg[\[Nu]]<=\[Pi]/2,roots=Conjugate[roots],-\[Pi]/2<=Arg[\[Nu]]<=0,roots];
 If[0<\[Lambda]<1,roots=roots/\[Lambda],roots]
]/;(ReQ[\[Lambda]]&&\[Lambda]>0 &&\[Lambda]!=1)


(* ::Subsubsection:: *)
(*Definitions for X01*)

Xi01[Z_,\[Nu]_,\[Lambda]_]:=AiryAi[\[Nu]^(2/3)zi[Z]]AiryBiPrime[\[Nu]^(2/3)zi[\[Lambda] Z]]-
  AiryBi[\[Nu]^(2/3)zi[Z]]AiryAiPrime[\[Nu]^(2/3)zi[\[Lambda] Z]]

Xi01Z0m1[m1_(*Integer*),\[Nu]_(*Complex*),\[Lambda]_Real]:=
  zed[(3/2Abs[m1-3/4]\[Pi]/Abs[\[Nu]])^(2/3)Exp[(-\[Pi]+2Abs[Arg[\[Nu]]]/3)I]]/\[Lambda]/;
   ((*(1\[LessEqual]m1\[LessEqual]M01m1[\[Nu],\[Lambda]])&&*)\[Lambda]>1)

Xi01Z0m2[m2_(*Integer*),\[Nu]_(*Complex*),\[Lambda]_Real]:=
  zed[(3/2Abs[m2-1/4]\[Pi]/Abs[\[Nu]])^(2/3)Exp[(-\[Pi]/3+2Abs[Arg[\[Nu]]]/3)I]]/;
   ((*(1\[LessEqual]m2\[LessEqual]M01m2[\[Nu],\[Lambda]])&&*)\[Lambda] > 1)

Xi01Z0m3[m3_(*Integer*), \[Nu]_Complex,\[Lambda]_Real]:=
Block[{Z,\[Nu]p,Zp},
 If[Arg[\[Nu]]<=0,\[Nu]p=\[Nu],\[Nu]p=Conjugate[\[Nu]]];
 Z=Zp/.FindRoot[\[Nu]p(\[Rho][Zp]-\[Rho][\[Lambda] Zp])+(m3-1/2) \[Pi] I,{Zp,1+I},Evaluated->False];
 If[Arg[\[Nu]]<=0,Z,Conjugate[Z]]
]/;((*(M01m3[\[Nu],\[Lambda]]\[LessEqual]m3)&&*)\[Lambda]>1)

Xi01Z0m3[m3_(*Integer*), \[Nu]_Real,\[Lambda]_Real]:=
Block[{Z},
 Re@Z/.FindRoot[\[Nu](\[Rho][Z]-\[Rho][\[Lambda] Z])+(m3-1/2) \[Pi] I,{Z,1},Evaluated->False]
]/;((*(M01m3[\[Nu],\[Lambda]]\[LessEqual]m3)&&*)\[Lambda]>1)

M01m1[\[Nu]_,\[Lambda]_]:=
  (Floor@Re@(-I # \[Rho][\[Lambda] ptO\[Alpha][Arg[#],\[Lambda]]]/\[Pi]+3/4))&@If[Arg[\[Nu]]>0,Conjugate[\[Nu]],\[Nu]]

M01m2[\[Nu]_,\[Lambda]_]:=
  (Floor@Re@(I # \[Rho][ptO\[Alpha][Arg[#],\[Lambda]]]/\[Pi]+1/4))&@If[Arg[\[Nu]]>0,Conjugate[\[Nu]],\[Nu]]

M01m3[\[Nu]_,\[Lambda]_]:=
Block[{\[Nu]p=If[Arg[\[Nu]]>0,Conjugate[\[Nu]],\[Nu]]},
 ((Ceiling@Re@(I \[Nu]p (\[Rho][#]-\[Rho][\[Lambda] #])/\[Pi]+1/2))&@ptO\[Alpha][Arg[\[Nu]p],\[Lambda]])
]

Xi01Zm1[m1_,\[Nu]_,\[Lambda]_Real,N_Integer:1,showAllZn_Symbol:False]:=
Block[{z,zz,z0rule,znrules,p,prules,n,lhs,eqns,soln,\[Epsilon],N\[Epsilon]=If[N==0,0,Floor[(N-1)/2]],z0=Xi01Z0m1[m1,\[Nu],\[Lambda]]},
 If[N==0,Return[z0],z0rule=z[0]->z0];
 prules=Join[{p[z[0]]->\[Rho][z0],p[\[Lambda] z[0]]->\[Rho][\[Lambda] z0]},
  Table[Derivative[n][p][z[0]]->Derivative[n][rho][z0],{n,1,N}],
   Table[Derivative[n][p][\[Lambda] z[0]]->Derivative[n][rho][\[Lambda] z0],{n,1,N}]];
 zz=Sum[\[Epsilon]^n z[n],{n,0,N}]+O[\[Epsilon]]^(N+1);
 lhs=Series[2Exp[I \[Pi]/4+\[Nu] p[zz]]Cos[\[Pi]/4-I \[Nu] p[\[Lambda] zz]]+\[Epsilon] Exp[-\[Nu](p[zz]-p[\[Lambda] zz])],{\[Epsilon],0,N}];
 eqns=(Map[#==0&,CoefficientList[lhs,\[Epsilon]]][[2;;N+1]])/.prules//Simplify;
 soln=Solve[eqns,Table[z[n],{n,1,N}]];
 znrules=Flatten[Append[soln,{z0rule}]];
 If[N\[Epsilon]==0,
  If[showAllZn,
   {Normal[zz/.znrules]/.{\[Epsilon]->1},Table[z[n],{n,0,N}]/.znrules},
   Normal[zz/.znrules]/.{\[Epsilon]->1}],
  If[showAllZn,
   {PadeApproximant[(zz/.znrules),{\[Epsilon],0,N\[Epsilon]}]/.{\[Epsilon]->1},Table[z[n],{n,0,N}]/.znrules},
   PadeApproximant[(zz/.znrules),{\[Epsilon],0,N\[Epsilon]}]/.{\[Epsilon]->1}]
 ]
]/;((*(1\[LessEqual]m1\[LessEqual]M01m1[\[Nu],\[Lambda]])&&*) N>=0&&\[Lambda]>1)

Xi01Zm2[m2_,\[Nu]_,\[Lambda]_Real,N_Integer:1,showAllZn_Symbol:False]:=
Block[{z,zz,z0rule,znrules,p,prules,n,lhs,eqns,soln,\[Epsilon],N\[Epsilon]=If[N==0,0,Floor[(N-1)/2]],z0=Xi01Z0m2[m2,\[Nu],\[Lambda]]},
 If[N==0,Return[z0],z0rule=z[0]->z0];
 prules=Join[{p[z[0]]->\[Rho][z0],p[\[Lambda] z[0]]->\[Rho][\[Lambda] z0]},
  Table[Derivative[n][p][z[0]]->Derivative[n][rho][z0],{n,1,N}],
   Table[Derivative[n][p][\[Lambda] z[0]]->Derivative[n][rho][\[Lambda] z0],{n,1,N}]];
 zz=Sum[\[Epsilon]^n z[n],{n,0,N}]+O[\[Epsilon]]^(N+1);
 lhs=Series[2 Exp[I \[Pi]/4+\[Nu] p[\[Lambda] zz]]Cos[\[Pi]/4-I \[Nu] p[zz]]+\[Epsilon] Exp[\[Nu](p[zz]-p[\[Lambda] zz])],{\[Epsilon],0,N}];
 eqns=(Map[#==0&,CoefficientList[lhs,\[Epsilon]]][[2;;N+1]])/.prules//Simplify;
 soln=Solve[eqns,Table[z[n],{n,1,N}]];
 znrules=Flatten[Append[soln,{z0rule}]];
 If[N\[Epsilon]==0,
  If[showAllZn,
   {Normal[zz/.znrules]/.{\[Epsilon]->1},Table[z[n],{n,0,N}]/.znrules},
   Normal[zz/.znrules]/.{\[Epsilon]->1}],
  If[showAllZn,
   {PadeApproximant[(zz/.znrules),{\[Epsilon],0,N\[Epsilon]}]/.{\[Epsilon]->1},Table[z[n],{n,0,N}]/.znrules},
   PadeApproximant[(zz/.znrules),{\[Epsilon],0,N\[Epsilon]}]/.{\[Epsilon]->1}]
 ]
]/;((*(1\[LessEqual]m2\[LessEqual]M01m[\[Nu],\[Lambda]])&&*) N>=0&&\[Lambda]>1)

Xi01Zm3[m3_,\[Nu]_,\[Lambda]_Real,N_Integer:1,showAllZn_Symbol:False]:=
Block[{z,zz,z0rule,znrules,p,prules,n,lhs,eqns,soln,\[Epsilon],N\[Epsilon]=If[N==0,0,Floor[(N-1)/2]],z0=Xi01Z0m3[m3,\[Nu],\[Lambda]]},
 If[N==0,Return[z0],z0rule=z[0]->z0];
 prules=Join[{p[z[0]]->\[Rho][z0],p[\[Lambda] z[0]]->\[Rho][\[Lambda] z0]},
  Table[Derivative[n][p][z[0]]->Derivative[n][rho][z0],{n,1,N}],
   Table[Derivative[n][p][\[Lambda] z[0]]->Derivative[n][rho][\[Lambda] z0],{n,1,N}]];
 zz=Sum[\[Epsilon]^n z[n],{n,0,N}]+O[\[Epsilon]]^(N+1);
 lhs=Series[2Cos[-I \[Nu](p[zz]- p[\[Lambda] zz])]+\[Epsilon] I Exp[\[Nu](p[zz]+p[\[Lambda] zz])],{\[Epsilon],0,N}];
 eqns=(Map[#==0&,CoefficientList[lhs,\[Epsilon]]][[2;;N+1]])/.prules//Simplify;
 soln=Solve[eqns,Table[z[n],{n,1,N}]];
 znrules=Flatten[Append[soln,{z0rule}]];
 If[N\[Epsilon]==0,
  If[showAllZn,
   {Normal[zz/.znrules]/.{\[Epsilon]->1},Table[z[n],{n,0,N}]/.znrules},
   Normal[zz/.znrules]/.{\[Epsilon]->1}],
  If[showAllZn,
   {PadeApproximant[(zz/.znrules),{\[Epsilon],0,N\[Epsilon]}]/.{\[Epsilon]->1},Table[z[n],{n,0,N}]/.znrules},
   PadeApproximant[(zz/.znrules),{\[Epsilon],0,N\[Epsilon]}]/.{\[Epsilon]->1}]
 ]
]/;((*(M01m3[\[Nu],\[Lambda]]\[LessEqual]m3)&&*) N>=0&&\[Lambda]>1)

findRootXi01m1[seeds_,\[Nu]_,\[Lambda]_]:=
Block[{Z,Z0,Z1,Z2,R,root,i,j,Nseeds=Length[seeds],m1,mm,debug=False},
 Table[
  Z0=seeds[[i]];
  root=Z/.FindRoot[Xi01[Z,\[Nu],\[Lambda]],{Z,Z0},Evaluated->False];
  If[i==1,Z1=1/\[Lambda],Z1=seeds[[i-1]]];
  If[i==Nseeds,Z2=Xi01Z0m1[Nseeds+1,\[Nu],\[Lambda]],Z2=seeds[[i+1]]];
  R=Min@Abs[Z0-{Z1,Z2}];
  If[debug,Print[i,"    ",{Z1,Z0,Z2},"     ",R];Print[root]];
  j=0;m1=i;mm=m1;
  While[(Abs[Z0-root]>=0.95R)&&mm>If[m1==1,3/4,m1-1],
   j=j+1;mm=m1-j*.001;Z0=Xi01Zm1[mm,\[Nu],\[Lambda],1];
   root=Z/.FindRoot[Xi01[Z,\[Nu],\[Lambda]],{Z,Z0},Evaluated->False,MaxIterations->1000];
   If[debug,Print["   ",j,"   ",mm,"    ",Z0,"    ",root]]
  ];
 root,{i,1,Nseeds}]
]

findRootXi01m2[seeds_,\[Nu]_,\[Lambda]_,m0_,roots1_,roots3_]:=
Block[{Z,Z0,Z1,Z2,R,root,i,j,Nseeds=Length[seeds],M1=Length[roots1],m2,mm,sa,test,roots13,r1,r3,debug=False},
 sa[Z_,a_:10]:=SetAccuracy[Z,a];
 Which[
  (m0==-1)||(m0==0),
   roots13=Join[roots1,roots3];test[rts_,rt_]:=ContainsAny[rts,{rt},SameTest->((sa[#]==sa[#2])&)];,
  m0==1,
   If[M1==0,r1=Null,r1=Last[roots1]];r3=First[roots3];
   If[TrueQ[sa[r1]==sa[r3]],
   roots13=Join[roots1,Rest[roots3]];
   test[rts_,rt_]:=ContainsAny[rts,{rt},SameTest->((sa[#]==sa[#2])&)];,
   roots13=Join[Most[roots1],Rest[roots3]];
   test[rts_,rt_]:=(ContainsAny[rts,{rt},SameTest->((sa[#]==sa[#2])&)]|| !(TrueQ[sa[rt]==sa[r1]]||(sa[rt]==sa[r3])));
 ];];
 Table[
  Z0=seeds[[i]];
  root=Z/.FindRoot[Xi01[Z,\[Nu],\[Lambda]],{Z,Z0},Evaluated->False];
  If[i==1,Z1=1+0.I,Z1=seeds[[i-1]]];
  If[i==Nseeds,Z2=ptO\[Alpha][Arg[\[Nu]],\[Lambda]],Z2=seeds[[i+1]]];
  R=Min@Abs[Z0-{Z1,Z2}];
  If[debug,Print[i,"    ",{Z2,Z0,Z1},"     ",R];Print[root]];
  j=0;m2=i;mm=m2;
  While[(Abs[Z0-root]>=0.95R)&&mm>If[m2==1,1/4,m2-1]&& test[roots13,root],
   j=j+1;mm=m2-j*.001;Z0=Xi01Zm2[mm,\[Nu],\[Lambda],1];
   root=Z/.FindRoot[Xi01[Z,\[Nu],\[Lambda]],{Z,Z0},Evaluated->False,MaxIterations->1000];
   If[debug,Print["   ",j,"   ",mm,"    ",Z0,"    ",root]];
  ];
  root,{i,1,Nseeds}]
]

Xi01Roots[\[Nu]_(*Complex*),\[Lambda]_(*Real*),Nz_Integer:10]:=
Block[{\[Nu]p,\[Lambda]p,seeds1,seeds2,seeds3,roots1,roots2,roots3,Z,m1,m2,m3,m0,M1,M2,M3,
       s1,s2,s3,r1,r2,r3,dd,pp,roots,sa,sortFail=False,O\[Alpha],s0,r0,eq22bits,debug=False},
 sa[Z_,a_:10]:=SetAccuracy[Z,a];
 \[Nu]p=If[Arg[\[Nu]]>0,Conjugate[\[Nu]],\[Nu]];
 \[Lambda]p=If[0<\[Lambda]<1,1/\[Lambda],\[Lambda]];
 If[debug,Print["Xi01Roots:0"];Print["\[Nu]: ",\[Nu]," \[Equivalent] ",\[Nu]p,", \[Lambda]: ",\[Lambda]," \[Equivalent] ",\[Lambda]p]];
 M1=M01m1[\[Nu]p,\[Lambda]p];
 If[Abs@Im[\[Nu]p]<0.1,
  seeds1=Table[Xi01Z0m1[m1,\[Nu]p,\[Lambda]p],{m1,1,M1}],
  seeds1=Table[Xi01Zm1[m1,\[Nu]p,\[Lambda]p,1],{m1,1,M1}]];
 roots1=findRootXi01m1[Chop@seeds1,\[Nu]p,\[Lambda]p];
 If[M1!=#,Print["Unexpected number of roots on seedline 1. Expected: ",M1,", Found: ",#,"."]; Print[{\[Nu],\[Lambda]}];
  Print[seeds1]; Print[roots1]; Abort[];]&@(Length@Union[roots1,SameTest->((sa[#]==sa[#2])&)]);
 If[debug,Print["Xi01Roots:1"];Print[seeds1];Print[roots1]];
 M2=M01m2[\[Nu]p,\[Lambda]p];
 seeds2=Table[Xi01Zm2[m2,\[Nu]p,\[Lambda]p,5],{m2,1,M2}];
 M3=M01m3[\[Nu]p,\[Lambda]p];
 If[debug,Print["Xi01Roots:3.0"];Print[{M3,M3+(Nz-1),Nz}]];
 If[Abs@Im[\[Nu]p]<0.1,
  seeds3=Table[Xi01Z0m3[m3,\[Nu]p,\[Lambda]p],{m3,M3,M3+(Nz-1)}],
  seeds3=Table[Xi01Zm3[m3,\[Nu]p,\[Lambda]p,1],{m3,M3,M3+(Nz-1)}]];
 roots3=Z/.FindRoot[Xi01[Z,\[Nu]p,\[Lambda]p],{Z,#},Evaluated->False]&/@seeds3;
 If[Nz!=#,Print["Unexpected number of roots on seedline 3. Expected: ",Nz,", Found: ",#,"."]; Print[{\[Nu],\[Lambda]}];
  Print[seeds3]; Print[roots3]; Abort[];]&@(Length@Union[roots3,SameTest->((sa[#]==sa[#2])&)]);
 If[debug,Print["Xi01Roots:3.1"];Print[seeds3];Print[roots3]];
 m0=M1+M2-M3+1;If[debug,Print[{M1,M2,M3}," m0 = ",m0]];
 roots2=findRootXi01m2[seeds2,\[Nu]p,\[Lambda]p,m0,roots1,roots3];
 If[M2!=#,Print["Unexpected number of roots on seedline 2. Expected: ",M2,", Found: ",#,"."];  Print[{\[Nu],\[Lambda]}];
  Print[seeds2]; Print[roots2]; Abort[];]&@ (Length@Union[roots2,SameTest->((sa[#]==sa[#2])&)]);
 If[debug,Print["Xi01Roots:2"];Print[seeds2];Print[roots2]];
 roots=Union[roots1,roots2,roots3,SameTest->((sa[#]==sa[#2])&)];
 Switch[m0,
  1, (* one extra seed *)
  If[M1+M2+Nz==Length[roots]+1,
   If[M1==0,s1=Null;r1=Null,s1=Last[seeds1];r1=Last[roots1]];
   If[M2==0,s2=Null;r2=Null,s2=Last[seeds2];r2=Last[roots2]];
   s3=First[seeds3];r3=First[roots3];
   dd=Abs@{r1-s1,r2-s2,r3-s3};
   Which[
    TrueQ[sa[r1]==sa[r2]],pp=Position[#,Max@#]&@ReplacePart[dd,3->0],
    TrueQ[sa[r1]==sa[r3]],pp=Position[#,Max@#]&@ReplacePart[dd,2->0],
    TrueQ[sa[r2]==sa[r3]],pp=Position[#,Max@#]&@ReplacePart[dd,1->0],
    True,Print["Unexpected sorting problem (m0 = ",m0,")."]; Print[{s1,s2,s3}]; Print[{r1,r2,r3}];Print[dd];Abort[];];
   If[debug,Print[{s1,s2,s3}];Print[{r1,r2,r3}];Print[dd];Print["pp = ",pp]];
   Switch[First@@pp,
    1,roots1=Most@roots1,
    2,roots2=Most@roots2 ,
    3,roots3=Rest@roots3
   ],
   Print["Unexpected number of roots. Expected: ",M1+M2+Nz-1,", Found: ",Length[roots]," (m0 = ",m0,")"];
   If[(M1+1)+(M2-1)+Nz==Length[roots],(* could be a roots = seeds case? *)
    Print["Could be a misclassified m0 = 0 case."]];sortFail=True;
  ],
  0, (* #roots == #seeds *)
  If[M1+M2+Nz!= Length[roots],
   Print["Unexpected number of roots. Expected: ",M1+M2+Nz,", Found: ",Length[roots]," (m0 = ",m0,")"];sortFail=True;],
  -1,  (* one extra root *)
  If[M1+M2+Nz==Length[roots],
  (* Manufacture seed for extra root by expanding Eq.(7) about O\[Alpha] *)
  If[debug,Print["Xi01Roots:4.0"]; Print["seeds1 = ",seeds1]; Print["roots1 = ",roots1]; Print["seeds2 = ",seeds2];
   Print["roots2 = ",roots2]; Print["seeds3 = ",seeds3]; Print["roots3 = ",roots3];];
  eq22bits[Z_,nu_,lambda_]:={I Exp[nu (rho[Z]+rho[lambda Z])],Exp[nu (rho[Z]-rho[lambda Z])],Exp[nu (-rho[Z]+rho[lambda Z])]};
  O\[Alpha]=ptO\[Alpha][Arg[\[Nu]p],\[Lambda]p];If[debug,Print["O\[Alpha] = ",O\[Alpha]]];
  s0=(Z/.Solve[Normal[Quiet@Series[Plus@@eq22bits[Z,\[Nu]p,\[Lambda]p],{Z,O\[Alpha],1},Assumptions->Im[Z]>0]]==0,Z][[1]]);
  r0=(Z/.FindRoot[Xi01[Z,\[Nu]p,\[Lambda]p],{Z,s0},Evaluated->False]);
  If[debug,Print["Xi01Roots:4.1"];Print[{N@s0,r0,Abs[r0-s0]}];
   Print[{#,Abs@#,Position[Abs@#,Min@Abs@#]}&@(N@eq22bits[s0,\[Nu]p,\[Lambda]p])];
   Print[{#,Abs@#,Position[Abs@#,Min@Abs@#]}&@(N@eq22bits[r0,\[Nu]p,\[Lambda]p])];];
  If[Length[roots]==Length@Union[{r0},roots,SameTest->((sa[#]==sa[#2])&)],Print["Extra seed did not converge to central root."];
   Print[{\[Nu],\[Lambda]}]; Print[s0]; Print[r0];  Print[roots]; Abort[];,
   (* else sort central root into one of the seedlines *)
   If[Re[\[Nu]p]==0,s1=Null;,s1=If[Im[\[Nu]p]==0,Xi01Z0m1[M1+1,\[Nu]p,\[Lambda]p],Xi01Zm1[M1+1,\[Nu]p,\[Lambda]p]]];
   If[Im[\[Nu]p]==0,s2=Null;,s2=Xi01Zm2[M2+1,\[Nu]p,\[Lambda]p]];
   If[(M3-1)==0,s3=Null;,s3=If[Im[\[Nu]p]==0,Xi01Z0m3[M3-1,\[Nu]p,\[Lambda]p],Xi01Zm3[M3-1,\[Nu]p,\[Lambda]p,1]]];
   dd=Abs[r0-{s1,s2,s3}];
   If[debug,Print[dd]; Print[{s1,s2,s3}] ;Print[{TrueQ[(s1==Null)&&(s3==Null)],TrueQ[(s2==Null)&&(s3==Null)],TrueQ[s1==Null],TrueQ[s2==Null],TrueQ[s3==Null]}]];
   Which[
    TrueQ[(s1==Null)&&(s3==Null)],pp=2,
    TrueQ[(s2==Null)&&(s3==Null)],pp=1,
    TrueQ[s1==Null],pp=Position[#,Min@#]&@ReplacePart[dd,1->10^10],
    TrueQ[s2==Null],pp=Position[#,Min@#]&@ReplacePart[dd,2->10^10],
    TrueQ[s3==Null],pp=Position[#,Min@#]&@ReplacePart[dd,3->10^10],
    True,pp=Position[dd,Min@dd]
   ];
   If[debug,Print[{s1,s2,s3}];Print[dd];Print["pp = ",pp]];
   Switch[First@@pp,
    1,roots1=Append[roots1,r0],
    2,roots2=Append[roots2,r0],
    3,roots3=Prepend[roots3,r0]
   ]
  ],
  Print["Unexpected number of roots. Expected: ",M1+M2+Nz,", Found: ",Length[roots]," (m0 = ",m0,")"]; sortFail=True;
  ],
  _,Print["Unexpected value for m0: ",{m0,{\[Nu]p,\[Lambda]p,M1,M2,M3}}, ". Aborting."];Abort[];
 ];
 If[sortFail,Print["Check whether the seeds have converged to the expected root. (m0 = ",m0,")"]; Print[{\[Nu]p,\[Lambda]p}];
  Print["seeds1 = ",seeds1]; Print["roots1 = ",roots1]; Print["seeds2 = ",seeds2]; Print["roots2 = ",roots2];
   Print["seeds3 = ",seeds3]; Print["roots3 = ",roots3]; Print["filtered roots = ",roots]; Abort[];];
 If[debug,Print["Xi01Roots:5"]];
 roots={roots1,roots2,roots3};
 Which[Arg[\[Nu]]==0,roots=Re@roots,Arg[\[Nu]]>0,roots=Conjugate[roots],True,roots];
 If[0<\[Lambda]<1,roots=roots/\[Lambda],roots]
]/;ReQ[\[Lambda]]&&(\[Lambda]>1)

BesselX01[z_,\[Nu]_,\[Lambda]_]:=
  (BesselJ[\[Nu],z](BesselY[-1+\[Nu],\[Lambda] z]-BesselY[1+\[Nu],\[Lambda] z])-
    BesselY[\[Nu],z](BesselJ[-1+\[Nu],\[Lambda] z]-BesselJ[1+\[Nu],\[Lambda] z]))/2

BesselX01Zeros[\[Nu]_,\[Lambda]_,Nz_:10]:=
Block[{\[Nu]p,\[Lambda]p,seeds,roots,z,i,debug=False},
 \[Lambda]p=If[0<\[Lambda]<1,1/\[Lambda],\[Lambda]];
 \[Nu]p=Which[-\[Pi]<Arg[\[Nu]]<-\[Pi]/2,Conjugate[-\[Nu]],\[Pi]/2<Arg[\[Nu]]<=\[Pi],-\[Nu],
  0<Arg[\[Nu]]<=\[Pi]/2,Conjugate[\[Nu]],-\[Pi]/2<=Arg[\[Nu]]<=0,\[Nu],True,Print["Unexpected value for Arg[\[Nu]]."];Abort[]];
 If[debug,Print["BesselX01Zeros:1"];Print[\[Nu]," \[Equivalent] ",\[Nu]p,", ",\[Lambda]," \[Equivalent] ",\[Lambda]p,", ",Nz]];
 seeds=\[Nu]p Xi01Roots[\[Nu]p,\[Lambda]p,Nz];
 If[debug,Print["BesselX01Zeros:2"];Print[seeds];Print[\[Nu]p seeds]];
 roots=Table[z/.FindRoot[BesselX01[z,\[Nu]p,\[Lambda]p],{z,#},Evaluated->False]&/@seeds[[i]],{i,1,3}];
 If[debug,Print["BesselX01Zeros:3"];Print[roots]];
 Which[Arg[\[Nu]]==0,roots=Re@roots,-\[Pi]<Arg[\[Nu]]<-\[Pi]/2,roots=-Conjugate[roots],
  \[Pi]/2<Arg[\[Nu]]<=\[Pi],roots=-roots,0<Arg[\[Nu]]<=\[Pi]/2,roots=Conjugate[roots],-\[Pi]/2<=Arg[\[Nu]]<=0,roots];
 If[0<\[Lambda]<1,roots=roots/\[Lambda],roots]
]/;(ReQ[\[Lambda]]&&\[Lambda]>0 &&\[Lambda]!=1)


(* ::Subsubsection:: *)
(*Definitions for X10*)

Xi10[Z_,\[Nu]_,\[Lambda]_]:=AiryAiPrime[\[Nu]^(2/3)zi[Z]]AiryBi[\[Nu]^(2/3)zi[\[Lambda] Z]]-
  AiryBiPrime[\[Nu]^(2/3)zi[Z]]AiryAi[\[Nu]^(2/3)zi[\[Lambda] Z]]

Xi10Z0m1[m1_(*Integer*),\[Nu]_(*Complex*),\[Lambda]_Real]:=
  zed[(3/2Abs[m1-1/4]\[Pi]/Abs[\[Nu]])^(2/3)Exp[(-\[Pi]+2Abs[Arg[\[Nu]]]/3)I]]/\[Lambda]/;
   ((*(1\[LessEqual]m1\[LessEqual]M10m1[\[Nu],\[Lambda]]+1)&&*)\[Lambda]>1)

Xi10Z0m2[m2_(*Integer*),\[Nu]_(*Complex*),\[Lambda]_Real]:=
  zed[(3/2Abs[m2-3/4]\[Pi]/Abs[\[Nu]])^(2/3)Exp[(-\[Pi]/3+2Abs[Arg[\[Nu]]]/3)I]]/;
   ((*(1\[LessEqual]m2\[LessEqual]M10m2[\[Nu],\[Lambda]])&&*)\[Lambda] > 1)

Xi10Z0m3[m3_(*Integer*), \[Nu]_Complex,\[Lambda]_Real]:=
Block[{Z,\[Nu]p,Zp},
 If[Arg[\[Nu]]<=0,\[Nu]p=\[Nu],\[Nu]p=Conjugate[\[Nu]]];
 Z=Zp/.FindRoot[\[Nu]p(\[Rho][Zp]-\[Rho][\[Lambda] Zp])+(m3-1/2) \[Pi] I,{Zp,1+I},Evaluated->False];
 If[Arg[\[Nu]]<=0,Z,Conjugate[Z]]
]/;((*(M10m3[\[Nu],\[Lambda]]\[LessEqual]m3)&&*)\[Lambda]>1)

Xi10Z0m3[m3_(*Integer*), \[Nu]_Real,\[Lambda]_Real]:=
Block[{Z},
 Re@Z/.FindRoot[\[Nu](\[Rho][Z]-\[Rho][\[Lambda] Z])+(m3-1/2) \[Pi] I,{Z,1},Evaluated->False]
]/;((*(M10m3[\[Nu],\[Lambda]]\[LessEqual]m3)&&*)\[Lambda]>1)

M10m1[\[Nu]_,\[Lambda]_]:=
  (Floor@Re@(-I # \[Rho][\[Lambda] ptO\[Alpha][Arg[#],\[Lambda]]]/\[Pi]+1/4))&@If[Arg[\[Nu]]>0,Conjugate[\[Nu]],\[Nu]]

M10m2[\[Nu]_,\[Lambda]_]:=
  (Floor@Re@(I # \[Rho][ptO\[Alpha][Arg[#],\[Lambda]]]/\[Pi]+3/4))&@If[Arg[\[Nu]]>0,Conjugate[\[Nu]],\[Nu]]

M10m3[\[Nu]_,\[Lambda]_]:=
Block[{\[Nu]p=If[Arg[\[Nu]]>0,Conjugate[\[Nu]],\[Nu]]},
 ((Ceiling@Re@(I \[Nu]p (\[Rho][#]-\[Rho][\[Lambda] #])/\[Pi]+1/2))&@ptO\[Alpha][Arg[\[Nu]p],\[Lambda]])
]

Xi10Zm1[m1_,\[Nu]_,\[Lambda]_Real,N_Integer:1,showAllZn_Symbol:False]:=
Block[{z,zz,z0rule,znrules,p,prules,n,lhs,eqns,soln,\[Epsilon],N\[Epsilon]=If[N==0,0,Floor[(N-1)/2]],z0=Xi10Z0m1[m1,\[Nu],\[Lambda]]},
 If[N==0,Return[z0],z0rule=z[0]->z0];
 prules=Join[{p[z[0]]->\[Rho][z0],p[\[Lambda] z[0]]->\[Rho][\[Lambda] z0]},
  Table[Derivative[n][p][z[0]]->Derivative[n][rho][z0],{n,1,N}],
   Table[Derivative[n][p][\[Lambda] z[0]]->Derivative[n][rho][\[Lambda] z0],{n,1,N}]];
 zz=Sum[\[Epsilon]^n z[n],{n,0,N}]+O[\[Epsilon]]^(N+1);
 lhs=Series[2I Exp[I \[Pi]/4+\[Nu] p[zz]]Sin[\[Pi]/4-I \[Nu] p[\[Lambda] zz]]-\[Epsilon] Exp[-\[Nu](p[zz]-p[\[Lambda] zz])],{\[Epsilon],0,N}];
 eqns=(Map[#==0&,CoefficientList[lhs,\[Epsilon]]][[2;;N+1]])/.prules//Simplify;
 soln=Solve[eqns,Table[z[n],{n,1,N}]];
 znrules=Flatten[Append[soln,{z0rule}]];
 If[N\[Epsilon]==0,
  If[showAllZn,
   {Normal[zz/.znrules]/.{\[Epsilon]->1},Table[z[n],{n,0,N}]/.znrules},
   Normal[zz/.znrules]/.{\[Epsilon]->1}],
  If[showAllZn,
   {PadeApproximant[(zz/.znrules),{\[Epsilon],0,N\[Epsilon]}]/.{\[Epsilon]->1},Table[z[n],{n,0,N}]/.znrules},
   PadeApproximant[(zz/.znrules),{\[Epsilon],0,N\[Epsilon]}]/.{\[Epsilon]->1}]
 ]
]/;((*(1\[LessEqual]m1\[LessEqual]M10m1[\[Nu],\[Lambda]])&&*)N>=0&&\[Lambda]>1)

Xi10Zm2[m2_,\[Nu]_,\[Lambda]_Real,N_Integer:1,showAllZn_Symbol:False]:=
Block[{z,zz,z0rule,znrules,p,prules,n,lhs,eqns,soln,\[Epsilon],N\[Epsilon]=If[N==0,0,Floor[(N-1)/2]],z0=Xi10Z0m2[m2,\[Nu],\[Lambda]]},
 If[N==0,Return[z0],z0rule=z[0]->z0];
 prules=Join[{p[z[0]]->\[Rho][z0],p[\[Lambda] z[0]]->\[Rho][\[Lambda] z0]},
  Table[Derivative[n][p][z[0]]->Derivative[n][rho][z0],{n,1,N}],
   Table[Derivative[n][p][\[Lambda] z[0]]->Derivative[n][rho][\[Lambda] z0],{n,1,N}]];
 zz=Sum[\[Epsilon]^n z[n],{n,0,N}]+O[\[Epsilon]]^(N+1);
 lhs=Series[2 I Exp[I \[Pi]/4+\[Nu] p[\[Lambda] zz]]Sin[\[Pi]/4-I \[Nu] p[zz]]-\[Epsilon] Exp[\[Nu](p[zz]-p[\[Lambda] zz])],{\[Epsilon],0,N}];
 eqns=(Map[#==0&,CoefficientList[lhs,\[Epsilon]]][[2;;N+1]])/.prules//Simplify;
 soln=Solve[eqns,Table[z[n],{n,1,N}]];
 znrules=Flatten[Append[soln,{z0rule}]];
 If[N\[Epsilon]==0,
  If[showAllZn,
   {Normal[zz/.znrules]/.{\[Epsilon]->1},Table[z[n],{n,0,N}]/.znrules},
   Normal[zz/.znrules]/.{\[Epsilon]->1}],
  If[showAllZn,
   {PadeApproximant[(zz/.znrules),{\[Epsilon],0,N\[Epsilon]}]/.{\[Epsilon]->1},Table[z[n],{n,0,N}]/.znrules},
   PadeApproximant[(zz/.znrules),{\[Epsilon],0,N\[Epsilon]}]/.{\[Epsilon]->1}]
 ]
]/;((*(1\[LessEqual]m2\[LessEqual]M10m2[\[Nu],\[Lambda]])&&*)N>=0&&\[Lambda]>1)

Xi10Zm3[m3_,\[Nu]_,\[Lambda]_Real,N_Integer:1,showAllZn_Symbol:False]:=
Block[{z,zz,z0rule,znrules,p,prules,n,lhs,eqns,soln,\[Epsilon],N\[Epsilon]=If[N==0,0,Floor[(N-1)/2]],z0=Xi10Z0m3[m3,\[Nu],\[Lambda]]},
 If[N==0,Return[z0],z0rule=z[0]->z0];
 prules=Join[{p[z[0]]->\[Rho][z0],p[\[Lambda] z[0]]->\[Rho][\[Lambda] z0]},
  Table[Derivative[n][p][z[0]]->Derivative[n][rho][z0],{n,1,N}],
   Table[Derivative[n][p][\[Lambda] z[0]]->Derivative[n][rho][\[Lambda] z0],{n,1,N}]];
 zz=Sum[\[Epsilon]^n z[n],{n,0,N}]+O[\[Epsilon]]^(N+1);
 lhs=Series[-2Cos[-I \[Nu](p[zz]- p[\[Lambda] zz])]+\[Epsilon] I Exp[\[Nu](p[zz]+p[\[Lambda] zz])],{\[Epsilon],0,N}];
 eqns=(Map[#==0&,CoefficientList[lhs,\[Epsilon]]][[2;;N+1]])/.prules//Simplify;
 soln=Solve[eqns,Table[z[n],{n,1,N}]];
 znrules=Flatten[Append[soln,{z0rule}]];
 If[N\[Epsilon]==0,
  If[showAllZn,
   {Normal[zz/.znrules]/.{\[Epsilon]->1},Table[z[n],{n,0,N}]/.znrules},
   Normal[zz/.znrules]/.{\[Epsilon]->1}],
  If[showAllZn,
   {PadeApproximant[(zz/.znrules),{\[Epsilon],0,N\[Epsilon]}]/.{\[Epsilon]->1},Table[z[n],{n,0,N}]/.znrules},
   PadeApproximant[(zz/.znrules),{\[Epsilon],0,N\[Epsilon]}]/.{\[Epsilon]->1}]
 ]
]/;((*(M10m3[\[Nu],\[Lambda]]\[LessEqual]m3)&&*)N>=0&&\[Lambda]>1)

findRootXi10m1[seeds_,\[Nu]_,\[Lambda]_]:=
Block[{Z,Z0,Z1,Z2,R,root,i,j,Nseeds=Length[seeds],m1,mm,debug=False},
 Table[
  Z0=seeds[[i]];
  root=Z/.FindRoot[Xi10[Z,\[Nu],\[Lambda]],{Z,Z0},Evaluated->False];
  If[i==1,Z1=1/\[Lambda],Z1=seeds[[i-1]]];
  If[i==Nseeds,Z2=Xi10Z0m1[Nseeds+1,\[Nu],\[Lambda]],Z2=seeds[[i+1]]];
  R=Min@Abs[Z0-{Z1,Z2}];
  If[debug,Print[i,"    ",{Z1,Z0,Z2},"     ",R];Print[root]];
  j=0;m1=i;mm=m1;
  While[(Abs[Z0-root]>=0.95R)&&mm>If[m1==1,1/4,m1-1],
   j=j+1;mm=m1-j*.001;Z0=Xi10Zm1[mm,\[Nu],\[Lambda],1];
   root=Z/.FindRoot[Xi10[Z,\[Nu],\[Lambda]],{Z,Z0},Evaluated->False,MaxIterations->1000];
   If[debug,Print["   ",j,"   ",mm,"    ",Z0,"    ",root]]
  ];
  root,{i,1,Nseeds}]
]

findRootXi10m2[seeds_,\[Nu]_,\[Lambda]_,m0_,roots1_,roots3_]:=
Block[{Z,Z0,Z1,Z2,R,root,i,j,Nseeds=Length[seeds],M1=Length[roots1],m2,mm,sa,test,roots13,r1,r3,debug=False},
 sa[Z_,a_:10]:=SetAccuracy[Z,a];
 Which[
  (m0==-1)||(m0==0),
   roots13=Join[roots1,roots3];test[rts_,rt_]:=ContainsAny[rts,{rt},SameTest->((sa[#]==sa[#2])&)];,
  m0==1,
   If[M1==0,r1=Null,r1=Last[roots1]];r3=First[roots3];
   If[TrueQ[sa[r1]==sa[r3]],
   roots13=Join[roots1,Rest[roots3]];
   test[rts_,rt_]:=ContainsAny[rts,{rt},SameTest->((sa[#]==sa[#2])&)];,
   roots13=Join[Most[roots1],Rest[roots3]];
   test[rts_,rt_]:=(ContainsAny[rts,{rt},SameTest->((sa[#]==sa[#2])&)]|| !(TrueQ[sa[rt]==sa[r1]]||(sa[rt]==sa[r3])));
 ];];
 Table[
  Z0=seeds[[i]];
  root=Z/.FindRoot[Xi10[Z,\[Nu],\[Lambda]],{Z,Z0},Evaluated->False];
  If[i==1,Z1=1+0.I,Z1=seeds[[i-1]]];
  If[i==Nseeds,Z2=ptO\[Alpha][Arg[\[Nu]],\[Lambda]],Z2=seeds[[i+1]]];
  R=Min@Abs[Z0-{Z1,Z2}];
  If[debug,Print[i,"    ",{Z2,Z0,Z1},"     ",R];Print[root]];
  j=0;m2=i;mm=m2;
  While[(Abs[Z0-root]>=0.95R)&&mm>If[m2==1,3/4,m2-1]&& test[roots13,root],
   j=j+1;mm=m2-j*.001;Z0=Xi10Zm2[mm,\[Nu],\[Lambda],1];
   root=Z/.FindRoot[Xi10[Z,\[Nu],\[Lambda]],{Z,Z0},Evaluated->False,MaxIterations->1000];
   If[debug,Print["   ",j,"   ",mm,"    ",Z0,"    ",root]];
  ];
  root,{i,1,Nseeds}]
]

Xi10Roots[\[Nu]_(*Complex*),\[Lambda]_(*Real*),Nz_Integer:10]:=
Block[{\[Nu]p,\[Lambda]p,seeds1,seeds2,seeds3,roots1,roots2,roots3,Z,m1,m2,m3,m0,M1,M2,M3,
      s1,s2,s3,r1,r2,r3,dd,pp,roots,sa,sortFail=False,O\[Alpha],s0,r0,eq22bits,debug=False},
 sa[Z_,a_:10]:=SetAccuracy[Z,a];
 \[Nu]p=If[Arg[\[Nu]]>0,Conjugate[\[Nu]],\[Nu]];
 \[Lambda]p=If[0<\[Lambda]<1,1/\[Lambda],\[Lambda]];
 If[debug,Print["Xi10Roots:0"];Print["\[Nu]: ",\[Nu]," \[Equivalent] ",\[Nu]p,", \[Lambda]: ",\[Lambda]," \[Equivalent] ",\[Lambda]p]];
 M1=M10m1[\[Nu]p,\[Lambda]p];
 seeds1=Table[Xi10Zm1[m1,\[Nu]p,\[Lambda]p,1],{m1,1,M1}];
 roots1=findRootXi10m1[Chop@seeds1,\[Nu]p,\[Lambda]p];
 If[M1!=#,Print["Unexpected number of roots on seedline 1. Expected: ",M1,", Found: ",#,"."]; Print[{\[Nu],\[Lambda]}];
  Print[seeds1]; Print[roots1]; Abort[];]&@(Length@Union[roots1,SameTest->((sa[#]==sa[#2])&)]);
 If[debug,Print["Xi10Roots:1"];Print[seeds1];Print[roots1]];
 M2=M10m2[\[Nu]p,\[Lambda]p];
 seeds2=Table[Xi10Zm2[m2,\[Nu]p,\[Lambda]p,5],{m2,1,M2}];
 M3=M10m3[\[Nu]p,\[Lambda]p];
 If[debug,Print["Xi10Roots:3.0"];Print[{M3,M3+(Nz-1),Nz}]];
 If[Abs@Im[\[Nu]p]<0.1,
  seeds3=Table[Xi10Z0m3[m3,\[Nu]p,\[Lambda]p],{m3,M3,M3+(Nz-1)}],
  seeds3=Table[Xi10Zm3[m3,\[Nu]p,\[Lambda]p,1],{m3,M3,M3+(Nz-1)}]];
 roots3=Z/.FindRoot[Xi10[Z,\[Nu]p,\[Lambda]p],{Z,#},Evaluated->False]&/@seeds3;
 If[Nz!=#,Print["Unexpected number of roots on seedline 3. Expected: ",Nz,", Found: ",#,"."]; Print[{\[Nu],\[Lambda]}];
  Print[seeds3]; Print[roots3]; Abort[];]&@(Length@Union[roots3,SameTest->((sa[#]==sa[#2])&)]);
 If[debug,Print["Xi10Roots:3.1"];Print[seeds3];Print[roots3]];
 m0=M1+M2-M3+1;If[debug,Print[{M1,M2,M3}," m0 = ",m0]];
 roots2=findRootXi10m2[seeds2,\[Nu]p,\[Lambda]p,m0,roots1,roots3];
 If[M2!=#,Print["Unexpected number of roots on seedline 2. Expected: ",M2,", Found: ",#,"."];  Print[{\[Nu],\[Lambda]}];
  Print[seeds2]; Print[roots2]; Abort[];]&@ (Length@Union[roots2,SameTest->((sa[#]==sa[#2])&)]);
 If[debug,Print["Xi10Roots:2"];Print[seeds2];Print[roots2]];
 roots=Union[roots1,roots2,roots3,SameTest->((sa[#]==sa[#2])&)];
 Switch[m0,
  1, (* one extra seed *)
  If[M1+M2+Nz==Length[roots]+1,
   If[M1==0,s1=Null;r1=Null,s1=Last[seeds1];r1=Last[roots1]];
   If[M2==0,s2=Null;r2=Null,s2=Last[seeds2];r2=Last[roots2]];
   s3=First[seeds3];r3=First[roots3];
   dd=Abs@{r1-s1,r2-s2,r3-s3};
   Which[
    TrueQ[sa[r1]==sa[r2]],pp=Position[#,Max@#]&@ReplacePart[dd,3->0],
    TrueQ[sa[r1]==sa[r3]],pp=Position[#,Max@#]&@ReplacePart[dd,2->0],
    TrueQ[sa[r2]==sa[r3]],pp=Position[#,Max@#]&@ReplacePart[dd,1->0],
    True,Print["Unexpected sorting problem (m0 = ",m0,")."]; Print[{s1,s2,s3}]; Print[{r1,r2,r3}];Print[dd];Abort[];];
   If[debug,Print[{s1,s2,s3}];Print[{r1,r2,r3}];Print[dd];Print["pp = ",pp]];
   Switch[First@@pp,
    1,roots1=Most@roots1,
    2,roots2=Most@roots2 ,
    3,roots3=Rest@roots3
   ],
   Print["Unexpected number of roots. Expected: ",M1+M2+Nz-1,", Found: ",Length[roots]," (m0 = ",m0,")"];
   If[(M1+1)+(M2-1)+Nz==Length[roots],(* could be a roots = seeds case? *)
    Print["Could be a misclassified m0 = 2 case."]];sortFail=True;
  ],
  0, (* #roots == #seeds *)
  If[M1+M2+Nz!= Length[roots],
   Print["Unexpected number of roots. Expected: ",M1+M2+Nz,", Found: ",Length[roots]," (m0 = ",m0,")"];sortFail=True;],
  -1,  (* one extra root *)
  If[M1+M2+Nz==Length[roots],
   (* Manufacture seed for extra root by expanding Eq.(7) about O\[Alpha] *)
   If[debug,Print["Xi10Roots:4.0"]; Print["seeds1 = ",seeds1]; Print["roots1 = ",roots1]; Print["seeds2 = ",seeds2];
    Print["roots2 = ",roots2]; Print["seeds3 = ",seeds3]; Print["roots3 = ",roots3];];
   eq22bits[Z_,nu_,lambda_]:={I Exp[nu (rho[Z]+rho[lambda Z])],-Exp[nu (rho[Z]-rho[lambda Z])],-Exp[nu (-rho[Z]+rho[lambda Z])]};
   O\[Alpha]=ptO\[Alpha][Arg[\[Nu]p],\[Lambda]p];If[debug,Print["O\[Alpha] = ",O\[Alpha]]];
   s0=(Z/.Solve[Normal[Quiet@Series[Plus@@eq22bits[Z,\[Nu]p,\[Lambda]p],{Z,O\[Alpha],1},Assumptions->Im[Z]>0]]==0,Z][[1]]);
   r0=(Z/.FindRoot[Xi10[Z,\[Nu]p,\[Lambda]p],{Z,s0},Evaluated->False]);
   If[debug,Print["Xi10Roots:4.1"];Print[{N@s0,r0,Abs[r0-s0]}];
    Print[{#,Abs@#,Position[Abs@#,Min@Abs@#]}&@(N@eq22bits[s0,\[Nu]p,\[Lambda]p])];
     Print[{#,Abs@#,Position[Abs@#,Min@Abs@#]}&@(N@eq22bits[r0,\[Nu]p,\[Lambda]p])];];
   If[Length[roots]==Length@Union[{r0},roots,SameTest->((sa[#]==sa[#2])&)],Print["Extra seed did not converge to central root."];
    Print[{\[Nu],\[Lambda]}]; Print[s0]; Print[r0];  Print[roots]; Abort[];,
   (* else sort central root into one of the seedlines *)
   If[Re[\[Nu]p]==0,s1=Null;,s1=If[Im[\[Nu]p]==0,Xi10Z0m1[M1+1,\[Nu]p,\[Lambda]p],Xi10Zm1[M1+1,\[Nu]p,\[Lambda]p]]];
   If[Im[\[Nu]p]==0,s2=Null;,s2=Xi01Zm2[M2+1,\[Nu]p,\[Lambda]p]];
   If[(M3-1)==0,s3=Null;,s3=If[Im[\[Nu]p]==0,Xi10Z0m3[M3-1,\[Nu]p,\[Lambda]p],Xi10Zm3[M3-1,\[Nu]p,\[Lambda]p,1]]];
   dd=Abs[r0-{s1,s2,s3}];
   If[debug,Print[dd]; Print[{s1,s2,s3}] ;Print[{TrueQ[(s1==Null)&&(s3==Null)],TrueQ[(s2==Null)&&(s3==Null)],TrueQ[s1==Null],TrueQ[s2==Null],TrueQ[s3==Null]}]];
   Which[
    TrueQ[(s1==Null)&&(s3==Null)],pp=2,
    TrueQ[(s2==Null)&&(s3==Null)],pp=1,
    TrueQ[s1==Null],pp=Position[#,Min@#]&@ReplacePart[dd,1->10^10],
    TrueQ[s2==Null],pp=Position[#,Min@#]&@ReplacePart[dd,2->10^10],
    TrueQ[s3==Null],pp=Position[#,Min@#]&@ReplacePart[dd,3->10^10],
    True,pp=Position[dd,Min@dd]
   ];
   If[debug,Print[{s1,s2,s3}];Print[dd];Print["pp = ",pp]];
   Switch[First@@pp,
    1,roots1=Append[roots1,r0],
    2,roots2=Append[roots2,r0],
    3,roots3=Prepend[roots3,r0]
   ]
  ],
  Print["Unexpected number of roots. Expected: ",M1+M2+Nz,", Found: ",Length[roots]," (m0 = ",m0,")"]; sortFail=True;
  ],
  _,Print["Unexpected value for m0: ",{m0,{\[Nu]p,\[Lambda]p,M1,M2,M3}}, ". Aborting."];Abort[];
 ];
 If[sortFail,Print["Check whether the seeds have converged to the expected root. (m0 = ",m0,")"]; Print[{\[Nu]p,\[Lambda]p}];
  Print["seeds1 = ",seeds1]; Print["roots1 = ",roots1]; Print["seeds2 = ",seeds2]; Print["roots2 = ",roots2];
   Print["seeds3 = ",seeds3]; Print["roots3 = ",roots3]; Print["filtered roots = ",roots]; Abort[];];
 If[debug,Print["Xi10Roots:5"]];
 roots={roots1,roots2,roots3};
 Which[Arg[\[Nu]]==0,roots=Re@roots,Arg[\[Nu]]>0,roots=Conjugate[roots],True,roots];
 If[0<\[Lambda]<1,roots=roots/\[Lambda],roots]
]/;ReQ[\[Lambda]]&&(\[Lambda]>1)

BesselX10[z_,\[Nu]_,\[Lambda]_]:=
  ((BesselJ[-1+\[Nu],z]-BesselJ[1+\[Nu],z])BesselY[\[Nu],\[Lambda] z]-(BesselY[-1+\[Nu],z]-BesselY[1+\[Nu],z])BesselJ[\[Nu],\[Lambda] z])/2

BesselX10Zeros[\[Nu]_,\[Lambda]_,Nz_:10]:=
Block[{\[Nu]p,\[Lambda]p,seeds,roots,z,i,debug=False},
 \[Lambda]p=If[0<\[Lambda]<1,1/\[Lambda],\[Lambda]];
 \[Nu]p=Which[-\[Pi]<Arg[\[Nu]]<-\[Pi]/2,Conjugate[-\[Nu]],\[Pi]/2<Arg[\[Nu]]<=\[Pi],-\[Nu],
  0<Arg[\[Nu]]<=\[Pi]/2,Conjugate[\[Nu]],-\[Pi]/2<=Arg[\[Nu]]<=0,\[Nu],True,Print["Unexpected value for Arg[\[Nu]]."];Abort[]];
 If[debug,Print["BesselX10Zeros:1"];Print[\[Nu]," \[Equivalent] ",\[Nu]p,", ",\[Lambda]," \[Equivalent] ",\[Lambda]p,", ",Nz]];
 seeds=\[Nu]p Xi10Roots[\[Nu]p,\[Lambda]p,Nz];
 If[debug,Print["BesselX10Zeros:2"];Print[seeds];Print[\[Nu]p seeds]];
 roots=Table[z/.FindRoot[BesselX10[z,\[Nu]p,\[Lambda]p],{z,#},Evaluated->False]&/@seeds[[i]],{i,1,3}];
 If[debug,Print["BesselX10Zeros:3"];Print[roots]];
 Which[Arg[\[Nu]]==0,roots=Re@roots,-\[Pi]<Arg[\[Nu]]<-\[Pi]/2,roots=-Conjugate[roots],\[Pi]/2<Arg[\[Nu]]<=\[Pi],roots=-roots,
  0<Arg[\[Nu]]<=\[Pi]/2,roots=Conjugate[roots],-\[Pi]/2<=Arg[\[Nu]]<=0,roots];
 If[0<\[Lambda]<1,roots=roots/\[Lambda],roots]
]/;(ReQ[\[Lambda]]&&\[Lambda]>0 &&\[Lambda]!=1)


End[]

Protect[BesselX00Zeros];
Protect[BesselX01Zeros];
Protect[BesselX10Zeros];
Protect[BesselX11Zeros];
Protect[BesselX00];
Protect[BesselX01];
Protect[BesselX10];
Protect[BesselX11];

SetAttributes[BesselX00Zeros,{ReadProtected,NumericFunction}];
SetAttributes[BesselX01Zeros,{ReadProtected,NumericFunction}];
SetAttributes[BesselX10Zeros,{ReadProtected,NumericFunction}];
SetAttributes[BesselX11Zeros,{ReadProtected,NumericFunction}];
SetAttributes[BesselX00,{ReadProtected,NumericFunction}];
SetAttributes[BesselX01,{ReadProtected,NumericFunction}];
SetAttributes[BesselX10,{ReadProtected,NumericFunction}];
SetAttributes[BesselX11,{ReadProtected,NumericFunction}];

EndPackage[]
