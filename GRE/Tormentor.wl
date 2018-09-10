(* ::Package:: *)

BeginPackage["GRETormentor`"]

DatabaseFormatize::usage="DatabaseFormatize[excel_] Formatizes the word list to be analyzed by the Initializer Function."

Initializer::usage="Reads a word list in xls or xlsx format, and transforms the list into a database to be manipulated by the main program."

Tormentor::usage="Torments your brain with words."

Begin["`Private`"]

DatabaseFormatize[excel_]:=
	Module[{meta,words,sticker,expression,permutedlist,list},
	meta=Drop[Import["GRE3000Template.xlsx"][[1]], 3];
	sticker=If[Not@MatchQ[#, _Real], StringRiffle[StringSplit[StringReplace[#,"||"->"\n"], "@"], "\n"], #] &;

	permutedlist=Map[Permute[Permute[Drop[Drop[#, {4}], -3], Cycles[{{2, 3, 4, 5}}]], Cycles[{{3, 4, 5}}]] &,meta];
	list=Map[sticker,permutedlist,{2}];
	
	words=list[[All,2]];

	expression="Expression: "<>#3<>"\n\n"<>"Antonyms: "<>#4<>"\n\n"<>"Synonyms: "<>#5<>"\n\n"<>"Examples: "<>#6&@@@list;

  	
	Export["formatized"<>excel,Thread[List[words,expression]]]
  ];


Initializer[excel_]:=Export["database.xlsx",Transpose[(Join@@{{Range[Length[#]]} , Transpose@PadRight[#,{Length[#],4}]})&[Map[ToString,Join@@Import[excel],{2}]]]];

Tormentor[weightfunction_, fontsize_] := 
  Module[{flag = 1, expressionflag = 0, database, algorithm, word, 
    wordcount = 0},
   database = First@Import["database.xlsx"];
   While[flag == 1,
    word = 
     RandomChoice[
      weightfunction /@ (#[[4]] & /@ database) -> database];
    wordcount++;
    DialogInput[
     {TextCell[Style[word[[2]], fontsize]],
      	Button[Style[#, fontsize] &@"WTF?", 
       database = MapAt[# + 1 &, database, {word[[1]], 4}]; 
       database = MapAt[# + 1 &, database, {word[[1]], 5}]; 
       DialogReturn[]],
      Button[Style[#, fontsize] &@"Rubbish!", 
       database = MapAt[# - 1 &, database, {word[[1]], 4}];
       	database = MapAt[# + 1 &, database, {word[[1]], 5}];
       	DialogReturn[]],
      Button[Style[#, fontsize] &@"WHAT'S THIS?",
       expressionflag=1;DialogReturn[];
       ],
      Button[Style[#, fontsize] &@"GG", flag = 0; DialogReturn[];]}
     ];
    If[expressionflag == 1, 
     DialogInput[{TextCell[Style[#, fontsize] &@(word[[2]] <> ":")], TextCell[Style[#, fontsize] &@word[[3]]],
       Row[{Button[Style[#, fontsize] &@"WTF?", 
          database = MapAt[# + 1 &, database, {word[[1]], 4}]; 
          database = MapAt[# + 1 &, database, {word[[1]], 5}]; 
          DialogReturn[]],
         Button[Style[#, fontsize] &@"I fucking know!", 
          database = MapAt[# - 1 &, database, {word[[1]], 4}];
          	database = MapAt[# + 1 &, database, {word[[1]], 5}];
          	DialogReturn[]]}, Spacer[2]]}]; expressionflag = 0];
    ];
   Print["Today I used " <> ToString@wordcount <> 
     " words to torment myself."];
   Export["database.xlsx", database]

   ];






End[]


EndPackage[]
