(* ::Package:: *)

BeginPackage["QHD`"]

QHDFormat::usage = "QHDFormat[stringlist_,MaxStringLength_] regularizes the output of strings, where \"MaxStringLength\" is used to restrict the space for the words."

GaussianIntegrate::usage = "GaussianIntegrate[variableList_, A_, B_, function_, end_] calculates the integral concerned with a gaussian function. \"A\" refers to the coefficient matrix for the variables in the second order, whilc \"B\" stands for the coefficient vector for the variables in the first order. Constants in the power cannot be calculated. 
Theoretically any arbitrary function before the gaussian function can be put into calculation, but a polynomial function is favored if a analytical result is needed. \"end\" defines the expanded series of exponential operator. Theoretically a higher value of it helps the program to return a better result."

ParallelGaussianIntegrate::usage = "Paralled version of GaussianIntegrate, which allows paralleled calculation of gaussian integrals."

SelfConvergingGaussianIntegrate::usage="SelfConvergingGaussianIntegrate[variableList_, A_, B_, function_, precisiongoal_] is designed for a converging gaussian integral of a function that cannot be derived to vanishment."

DerivativeExpand::usage= "DerivativeExpand[function1_, function2_, variablelist_, power_] helps calculate nested Poisson brackets, namely (\!\(\*FractionBox[\(\[PartialD]\\\ #1\), \(\[PartialD]x\)]\)\!\(\*FractionBox[\(\[PartialD]\\\ #2\), \(\[PartialD]p\)]\)-\!\(\*FractionBox[\(\[PartialD]\\\ #2\), \(\[PartialD]x\)]\)\!\(\*FractionBox[\(\[PartialD]\\\ #1\), \(\[PartialD]p\)]\)\!\(\*SuperscriptBox[\()\), \(n\)]\)[function1,function2],
which is intended to help calculate Moyal bracket. The variablelist must be a second-ordered list with all pairs of conjugated variables put in a list."

MoyalBracket::usage= "MoyalBracket[A_, H_, variableList_, end_] calculates the moyal bracket of two functions, A being an arbotrary function of observables, while H stands for the Hamiltonian function for the system.
end refers to the cutting point for the series of expansion. The variablelist must be a second-ordered list with all pairs of conjugated variables put in a list. h, namely Planck's constant, needs to be stated in numerical calculations."

ParallelMoyalBracket::usage= "Paralleled version of MoyalBracket, which allows paralleled calculation of moyal brackets."

UnitaryDiagonalization::usage= "UnitaryDiagonalization[covariancematrix_,variablelist_] is designed for the diagonization of variables that has corelation with each other. It returns the substitution rule for the variables."

UnitaryDiagonalizationMatrix::usage="UnitaryDiagonalizationMatrix[covariancematrix_] simply returns the unitary trasnformation for a Covariance Matrix that makes it diagonalized, meaning there's no correlation between the transformed variables."

QHDInput::usage="QHDInput[inputtype_,basicsets_,basicsetscut_,Hamiltonian_, InitialA_, InitialB_, variable_, timeinterval_, timegoal_, grade_, Integrateend_, filename_] is intended to help create an input file for QHD programs, which requires to input the type of the data (now being either the means of several functions or the coefficients of the basic sets), the type of sets(now either the polynomial sets or exponential sets.
	
	inputtype: \"Mean\",\"Coefficient\"

	basicsets:\"Polynomial\",\"Exponential\"

	A few others need to be set too, including the form of Hamiltonian function which can be written in traditional form as normally done in handwriting. It shall be mentioned that the variables must be written in a matrix form, with all the pairs of conjugates listed."

QHDInputRead::usage="QHDInputRead[filename_] is aimed to simply read a file created by QHDInput, and returns a list composed of the titles and their values."

PolynomialQHD::usage="PolynomialQHD[Hamiltonian_,InitialA_,InitialB_,meanlist_,variable_,timeinterval_,timegoal_,grade_,Integradeend_,maxstringlength_,outputfile_] contains the whole algorithm for simulating QHD, where the grades of the basic sets, the cutting of the integrals as well as the maximum length of words in the output need to be set manually."

IntegratedQHD::usage="IntegratedQHD[inputfilename_,outputfilename_] simply integrates all methods of QHD computation. The method needs to be instructed in the input file."


Begin["`Private`"]
h=1;

QHDFormat[stringlist_,MaxStringLength_] := (StringJoin @@ (# <> (StringJoin @@ Table[" ", MaxStringLength - StringLength[#]]) & /@ stringlist)) <> "\n";

GaussianIntegrate[variableList_,A_,B_,function_,end_]:=Module[{inverse,b,modifiedfunction},
inverse=Inverse[A];
b=inverse.B;
Evaluate[Sqrt[(2Pi)^Length[variableList]/Det[A]]*Exp[0.5*B.inverse.B]*

Sum[1/n!*Nest[Function[a,0.5*Flatten[inverse].(D[a,##]&@@@Tuples[variableList,2])],(function/.Rule@@@Transpose[{variableList,variableList+b}]),n],{n,0,end}]

]/.Rule@@@Transpose[{variableList,Table[0,Length[variableList]]}]
];

ParallelGaussianIntegrate[variableList_,A_,B_,function_,end_]:=Module[{inverse,b,modifiedfunction},
inverse=Inverse[A];
b=inverse.B;
Evaluate[Sqrt[(2Pi)^Length[variableList]/Det[A]]*Exp[0.5*B.inverse.B]*

ParallelSum[1/n!*Nest[Function[a,0.5*Flatten[inverse].(D[a,##]&@@@Tuples[variableList,2])],(function/.Rule@@@Transpose[{variableList,variableList+b}]),n],{n,0,end}]

]/.Rule@@@Transpose[{variableList,Table[0,Length[variableList]]}
]
];

SelfConvergingGaussianIntegrate[variableList_, A_, B_, function_, precisiongoal_] := 
Module[{inverse, b, modifiedfunction, temp, n = 0,result=0}, 
  inverse = Inverse[A]; 
  b = inverse.B; 
  temp=Evaluate[
      Sqrt[(2 Pi)^Length[variableList]/Det[A]]
        *Exp[0.5*B.inverse.B]*1/n!*
        Nest[Function[a, 0.5*Flatten[inverse].(D[a, ##] & @@@ Tuples[variableList, 2])], (function /. Rule @@@ Transpose[{variableList, variableList + b}]), n]
        ] /. Rule @@@ Transpose[{variableList, Table[0, Length[variableList]]}];
        Echo@temp;
  While[
    Abs[temp]>10^(-precisiongoal),
    result+=temp;
    n++;
    temp=Evaluate[
      Sqrt[(2 Pi)^Length[variableList]/Det[A]]
        *Exp[0.5*B.inverse.B]*1/n!*
        Nest[Function[a, 0.5*Flatten[inverse].(D[a, ##] & @@@ Tuples[variableList, 2])], (function /. Rule @@@ Transpose[{variableList, variableList + b}]), n]
        ] /. Rule @@@ Transpose[{variableList, Table[0, Length[variableList]]}];
    ];

    result
];

DerivativeExpand[function1_,function2_,variablelist_,power_]:=Module[{A,B,a,b,x,Regularizer,Coefficientlist,listtemp},
Regularizer[SubList_]:=Module[{temp},
temp=If[Nor[MatchQ[#,_Integer],MatchQ[#,_List]],{#,1},#]&/@SubList;
If[Not[MatchQ[First[temp],_Integer]],PrependTo[temp,1]];
temp
];
listtemp=Regularizer/@ReleaseHold[Hold[Evaluate@Expand[Total[A[#1]*B[#2]-A[#2]B[#1]&@@@variablelist]^power]]/.{Power->List,Times->List,Plus->List}]/.{A[x_]->{a,x},B[x_]->{b,x}};
Coefficientlist=First/@listtemp;
listtemp=Map[Drop[#,1]&,Map[GatherBy[#,First]&,Map[Flatten,Drop[#,1]&/@listtemp,{2}]],{3}];
(*Now the first element in each list refers to the derivatives of the first function*)
Coefficientlist.(Fold[D,function1,#[[1]]]*Fold[D,function2,#[[2]]]&/@listtemp)
];
MoyalBracket[A_,H_,variableList_,end_]:=Sum[(-1)^j/(2j+1)!*(h/4Pi)^(2j)*DerivativeExpand[A,H,variableList,2j+1],{j,0,end}];

ParallelMoyalBracket[A_,H_,variableList_,end_]:=ParallelSum[(-1)^j/(2j+1)!*(h/4Pi)^(2j)*DerivativeExpand[A,H,variableList,2j+1],{j,0,end}];

UnitaryDiagonalization[covariancematrix_,variablelist_]:=Rule@@@Transpose[{variablelist,Orthogonalize[Eigenvectors[covariancematrix]].variablelist}];

UnitaryDiagonalizationMatrix[covariancematrix_]:=Orthogonalize[Eigenvectors[covariancematrix]];


QHDInput[inputtype_,basicsets_,basicsetscut_,Hamiltonian_, InitialA_, InitialB_, variable_, timeinterval_, timegoal_, grade_, Integrateend_, filename_] :=
  Module[{string1, string2, Format, MaxStringLength = 20},
   Format[stringlist_] := (StringJoin @@ (# <> (StringJoin @@ Table[" ", MaxStringLength - StringLength[#]]) & /@ stringlist)) <> "\n";
   string1 = StringJoin @@ (Format /@ {
            {"Input type:",inputtype},
            {"Basic sets type:",basicsets},
            {"Basic sets cut:",ToString@basicsetscut},
            {"Hamiltonian:", ToString@FullForm[Hamiltonian]},
            {"Initial A:", ToString@FullForm[InitialA]},
            {"Initial B:", ToString@FullForm[InitialB]},
            {"variable:", ToString@variable},
            {"Time interval:", ToString@timeinterval},
            {"Time Goal:", ToString@timegoal},
            {"QHD grade:",ToString@grade},
            {"Integration Cut:", ToString@Integrateend}
            });
   
   string2 = Which[
     inputtype== "Mean",StringJoin @@ (Format /@ Transpose[{(ToString[#, InputForm] <> ":") & /@ Sort@DeleteDuplicates[Times @@@ Subsets[Flatten[Table[variable, basicsetscut]]]]}]),

     inputtype=="Coefficient", 
        Which[
          basicsets=="Polynomial",StringJoin @@ (Format /@ Transpose[{(ToString[#, InputForm] <> ":") & /@ Sort@DeleteDuplicates[Times @@@ Subsets[Flatten[Table[variable, basicsetscut]]]]}]),

          basicsets=="Exponential",StringJoin @@ (Format /@ Transpose[{(ToString[#, InputForm] <> ":") & /@ Exp /@Sort@DeleteDuplicates[Plus @@@ Subsets[Flatten[Table[variable, basicsetscut]]]]}])

          ]

   ];

 Export[NotebookDirectory[] <> filename, string1 <> "\n" <> string2]
];


QHDInputRead[filename_] := Module[{}, 
	Transpose[{
		First/@Select[StringSplit[#, ":"] & /@ StringSplit[Import[filename], "\n"], Length[#] > 1 &],
			MapAt[ToString,Transpose[
				Map[ToExpression, 
					Select[StringSplit[#, ":"] & /@ StringSplit[Import[filename], "\n"], Length[#] > 1 &], 
					{2}]
				][[2]],
			{{1},{2}}]
	}]
];


PolynomialQHD[Hamiltonian_,InitialA_,InitialB_,meanlist_,variable_,timeinterval_,timegoal_,grade_,Integradeend_,maxstringlength_,outputfile_]:=
Module[{fileoperator,MoyalBracketList,PolynomialList,Normalizer,Tuner,TunedMeanList,VarianceList,TemporaryCoefficientList,TemporaryA,TemporaryB,EquationList,TunerList,UnitaryTransformation,Polynomialfunction,coefficient,equationlist,CoefficientResult,timecounter=0,MeanChangeList,coefficientlist,polynomialfunction,tempseries,coefficients,result},
Tuner[polynomialfunction_]:=If[MatchQ[polynomialfunction/.{Times->List,Power->List},_List],Times@@(ParallelMap[1/#!&,Cases[Flatten[polynomialfunction/.{Times->List,Power->List}],_Integer]]),1];

fileoperator = OpenWrite[outputfile];

PolynomialList=Sort@DeleteDuplicates[Times@@@Subsets[Flatten[Table[variable,grade]]]];
TunerList=Tuner/@PolynomialList;
TunedMeanList=TunerList*meanlist;
PolynomialList=TunerList*PolynomialList;

MoyalBracketList=ParallelMap[MoyalBracket[#,Hamiltonian,variable,grade]&,PolynomialList];

TemporaryA=InitialA;
TemporaryB=InitialB;



coefficientlist=Table[tempseries[q],{q,Length[PolynomialList]}];

polynomialfunction=coefficientlist.PolynomialList;


result=Solve[ParallelMap[GaussianIntegrate[Flatten[variable],TemporaryA,TemporaryB,polynomialfunction*#,Integradeend]&,PolynomialList]==TunedMeanList,coefficientlist];
   
TemporaryCoefficientList=Flatten[Table[tempseries[q],{q,Length[PolynomialList]}]/.result];

Normalizer=1/GaussianIntegrate[Flatten[variable],TemporaryA,TemporaryB,TemporaryCoefficientList.PolynomialList,Integradeend];


WriteString[fileoperator,
	"Initial Status:"<>"\n"<>
	StringJoin@@(QHDFormat[#,maxstringlength]&/@(Transpose@{ToString[#,InputForm]&/@PolynomialList,ToString[#,InputForm]&/@TemporaryCoefficientList}))<>"\n"<>"Mean:"<>"\n"<>StringJoin@@(QHDFormat[#,maxstringlength]&/@(Transpose@{ToString[#,InputForm]&/@PolynomialList,ToString[#,InputForm]&/@TunedMeanList}))<>"\n"
];


While[timecounter<=timegoal,
timecounter+=timeinterval;

MeanChangeList=timeinterval*ParallelMap[GaussianIntegrate[Flatten[variable],TemporaryA,TemporaryB,Normalizer*(TemporaryCoefficientList.PolynomialList)*#,Integradeend]&, MoyalBracketList];
TunedMeanList+=MeanChangeList;

coefficientlist=Table[tempseries[q],{q,Length[PolynomialList]}];
polynomialfunction=coefficientlist.PolynomialList;

equations=ParallelMap[GaussianIntegrate[Flatten[variable],TemporaryA,TemporaryB,Normalizer*polynomialfunction*#,Integradeend]&,PolynomialList];


TemporaryCoefficientList=Flatten[Table[tempseries[q],{q,Length[PolynomialList]}]/.result];


Normalizer=1/GaussianIntegrate[Flatten[variable],TemporaryA,TemporaryB,TemporaryCoefficientList.PolynomialList,Integradeend];

WriteString[fileoperator,
	QHDFormat[{"t:",ToString@timecounter},maxstringlength]<>
	StringJoin@@(QHDFormat[#,maxstringlength]&/@(Transpose@{ToString[#,InputForm]&/@PolynomialList,ToString[#,InputForm]&/@TemporaryCoefficientList}))<>"\n"<>"Mean:"<>"\n"<>StringJoin@@(QHDFormat[#,maxstringlength]&/@(Transpose@{ToString[#,InputForm]&/@PolynomialList,ToString[#,InputForm]&/@TunedMeanList}))<>"\n"
];
  ];


Close[fileoperator];

];

CoefficientPolynomialQHD[Hamiltonian_,InitialA_,InitialB_,Coefficient_,variable_,timeinterval_,timegoal_,grade_,Integradeend_,maxstringlength_,outputfile_]:=
Module[{fileoperator,MoyalBracketList,PolynomialList,Normalizer,Tuner,TunedMeanList,VarianceList,TemporaryCoefficientList,TemporaryA,TemporaryB,EquationList,TunerList,UnitaryTransformation,Polynomialfunction,coefficient,equationlist,CoefficientResult,timecounter=0,MeanChangeList,coefficientlist,polynomialfunction,tempseries,coefficients,result},
Tuner[polynomialfunction_]:=If[MatchQ[polynomialfunction/.{Times->List,Power->List},_List],Times@@(ParallelMap[1/#!&,Cases[Flatten[polynomialfunction/.{Times->List,Power->List}],_Integer]]),1];

fileoperator = OpenWrite[outputfile];

PolynomialList=Sort@DeleteDuplicates[Times@@@Subsets[Flatten[Table[variable,grade]]]];
TunerList=Tuner/@PolynomialList;
PolynomialList=TunerList*PolynomialList;

MoyalBracketList=ParallelMap[MoyalBracket[#,Hamiltonian,variable,grade]&,PolynomialList];

TemporaryA=InitialA;
TemporaryB=InitialB;



coefficientlist=Table[tempseries[q],{q,Length[PolynomialList]}];

polynomialfunction=coefficientlist.PolynomialList;


result=Solve[ParallelMap[GaussianIntegrate[Flatten[variable],TemporaryA,TemporaryB,polynomialfunction*#,Integradeend]&,PolynomialList]==TunedMeanList,coefficientlist];
   
TemporaryCoefficientList=Flatten[Table[tempseries[q],{q,Length[PolynomialList]}]/.result];

Normalizer=1/GaussianIntegrate[Flatten[variable],TemporaryA,TemporaryB,TemporaryCoefficientList.PolynomialList,Integradeend];


WriteString[fileoperator,
  "Initial Status:"<>"\n"<>
  StringJoin@@(QHDFormat[#,maxstringlength]&/@(Transpose@{ToString[#,InputForm]&/@PolynomialList,ToString[#,InputForm]&/@TemporaryCoefficientList}))<>"\n"<>"Mean:"<>"\n"<>StringJoin@@(QHDFormat[#,maxstringlength]&/@(Transpose@{ToString[#,InputForm]&/@PolynomialList,ToString[#,InputForm]&/@TunedMeanList}))<>"\n"
];


While[timecounter<=timegoal,
timecounter+=timeinterval;

MeanChangeList=timeinterval*ParallelMap[GaussianIntegrate[Flatten[variable],TemporaryA,TemporaryB,Normalizer*(TemporaryCoefficientList.PolynomialList)*#,Integradeend]&, MoyalBracketList];
TunedMeanList+=MeanChangeList;

coefficientlist=Table[tempseries[q],{q,Length[PolynomialList]}];
polynomialfunction=coefficientlist.PolynomialList;

equations=ParallelMap[GaussianIntegrate[Flatten[variable],TemporaryA,TemporaryB,Normalizer*polynomialfunction*#,Integradeend]&,PolynomialList];


TemporaryCoefficientList=Flatten[Table[tempseries[q],{q,Length[PolynomialList]}]/.result];


Normalizer=1/GaussianIntegrate[Flatten[variable],TemporaryA,TemporaryB,TemporaryCoefficientList.PolynomialList,Integradeend];

WriteString[fileoperator,
  QHDFormat[{"t:",ToString@timecounter},maxstringlength]<>
  StringJoin@@(QHDFormat[#,maxstringlength]&/@(Transpose@{ToString[#,InputForm]&/@PolynomialList,ToString[#,InputForm]&/@TemporaryCoefficientList}))<>"\n"<>"Mean:"<>"\n"<>StringJoin@@(QHDFormat[#,maxstringlength]&/@(Transpose@{ToString[#,InputForm]&/@PolynomialList,ToString[#,InputForm]&/@TunedMeanList}))<>"\n"
];
  ];


Close[fileoperator];

];




IntegratedQHD[inputfilename_,outputfilename_] := Module[{Parameters,Coefficients,inputtype,basicssets,basicsetscut,hamiltonian,initalA,initialB,meanlist,variable,timeinterval,timegoal,grade,integradeend},
	
	Parameters = Rule@@@QHDInputRead[inputfilename];
	

	{inputtype,basicsets,basicsetscut,hamiltonian,initialA,initialB,variable,timeinterval,timegoal,grade,integradeend}=({"Input type","Basic sets type","Basic sets cut","Hamiltonian","Initial A","Initial B","variable","Time interval","Time Goal","QHD grade","Integration Cut"}/.Parameters);



	meanlist=Drop[Transpose[QHDInputRead[inputfilename]][[2]],11];

	(*string1 = StringJoin @@ (Format /@ {
            {"Input type:",inputtype},
            {"Basic sets type:",basicsets},
            {"Basic sets cut:",ToString@basicsetscut},
            {"Hamiltonian:", ToString@FullForm[Hamiltonian]},
            {"Initial A:", ToString@FullForm[InitialA]},
            {"Initial B:", ToString@FullForm[InitialB]},
            {"variable:", ToString@variable},
            {"Time interval:", ToString@timeinterval},
            {"Time Goal:", ToString@timegoal},
            {"QHD grade:",ToString@grade},
            {"Integration Cut:", ToString@Integrateend}
            });*)

	If[inputtype=="Mean"&&basicsets=="Polynomial",
		PolynomialQHD[hamiltonian,initialA,initialB,meanlist,variable,timeinterval,timegoal,grade,integradeend,30,outputfilename]
	]

];




End[]


EndPackage[]
