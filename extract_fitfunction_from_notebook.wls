#!/usr/bin/env wolframscript
(* Extract Run2 and Run3 from StanCodeGeneratorBasis.nb into a standalone .wls. *)

ClearAll[heldText, heldBodyText, shiftedInputText, massInputText];

heldText[expr_] := StringReplace[
  ToString[expr, InputForm, PageWidth -> Infinity],
  "Stan_intfunctions/" -> "Stan_inputs/"
];

heldBodyText[expr_] := Module[{held, text},
  held = expr /. HoldComplete[x___] :> HoldComplete[CompoundExpression[x]];
  text = heldText[held];
  If[StringStartsQ[text, "HoldComplete["] && StringEndsQ[text, "]"],
    StringTake[text, {14, -2}],
    text
  ]
];

shiftedInputText[inputs_, indices_List] :=
  StringRiffle[heldBodyText /@ inputs[[indices]], "\n\n"];

massInputText[] := StringRiffle[
  {
    "GF = 1.166371/10^5;",
    "mb = 4.18;",
    "mc = 0.915;",
    "mtau = 1.77686;",
    "mB = 5.27966;",
    "mBst = 5.32470;",
    "mD = 1.86966;",
    "mDst = 2.01026;",
    "qDstMAX = (mB - mDst)^2;",
    "rD = mD/mB;",
    "rDst = mDst/mB;",
    "QpQmD[qSQ_] := ((mB + mD)^2 - qSQ)*((mB - mD)^2 - qSQ);",
    "QpQmDst[qSQ_] := ((mB + mDst)^2 - qSQ)*((mB - mDst)^2 - qSQ);",
    "qSQfuncDst[w_] := mB^2 + mDst^2 - 2*mB*mDst*w;",
    "qSQfuncD[w_] := mB^2 + mD^2 - 2*mB*mD*w;",
    "wFuncDst[qSQ_] := (mB^2 + mDst^2 - qSQ)/(2*mB*mDst);",
    "wFuncD[qSQ_] := (mB^2 + mD^2 - qSQ)/(2*mB*mD);",
    "wMax[\"D\"] = w /. NSolve[qSQfuncD[w] == 0, w][[1]];",
    "wMax[\"Dst\"] = w /. NSolve[qSQfuncDst[w] == 0, w][[1]];",
    "wMaxTau[\"D\"] = 1.429071189096182;",
    "wMaxTau[\"Dst\"] = 1.3548194509494718;",
    "BrDstDpi = 0.677;",
    "SecondToGeV = 1/(6.582119569*10^-25);",
    "LifeTime[\"Bp\"] = 1.638*10^-12*SecondToGeV;"
  },
  "\n"
];

Module[
  {
    nbFile, outFile, inputs, setupIndices, ubIndices, bszIndices, bglIndices, hqetIndices,
    header, body
  },
  nbFile = If[Length[$ScriptCommandLine] >= 2,
    $ScriptCommandLine[[2]],
    "StanCodeGeneratorBasis.nb"
  ];
  outFile = If[Length[$ScriptCommandLine] >= 3,
    $ScriptCommandLine[[3]],
    "Stan_FitFunction.wls"
  ];

  inputs = NotebookImport[nbFile, "Input"];

  (* Input indices corresponding to StanCodeGeneratorBasis.nb:
     Run2/Run3 setup: 30-57
     Integration/save: BSZ 58-59, BGL 60-61, HQET 62-63
     HQET UB setup: 64-68
  *)
  setupIndices = Range[30, 57];
  ubIndices = Range[64, 68];
  bszIndices = {58, 59};
  bglIndices = {60, 61};
  hqetIndices = {62, 63};

  header = StringRiffle[
    {
      "#!/usr/bin/env wolframscript",
      "(* ::Package:: *)",
      "",
      "(* Standalone copy of StanCodeGeneratorBasis.nb Run2 and Run3.",
      "   Runtime does not load or reference the notebook.",
      "   Integration outputs are read/written under Stan_inputs/. *)",
      "",
      "ClearAll[StanFitFunctionSetup, StanFitFunctionUBSetup, StanFitFunctionIntegrate, StanGenerateFitFunctions];",
      "",
      ""
    },
    "\n"
  ];

  body = StringRiffle[
    {
      "StanFitFunctionSetup[] := (",
      massInputText[],
      shiftedInputText[inputs, Rest[setupIndices]],
      ");",
      "",
      "StanFitFunctionUBSetup[] := (",
      shiftedInputText[inputs, ubIndices],
      ");",
      "",
      "StanFitFunctionIntegrate[\"BSZ\"] := (",
      shiftedInputText[inputs, bszIndices],
      ");",
      "",
      "StanFitFunctionIntegrate[\"BGL\"] := (",
      shiftedInputText[inputs, bglIndices],
      ");",
      "",
      "StanFitFunctionIntegrate[\"HQET\"] := (",
      shiftedInputText[inputs, hqetIndices],
      ");",
      "",
      "StanFitFunctionIntegrate[ffModel_String] := (",
      "  Message[StanFitFunctionIntegrate::unkff, ffModel];",
      "  $Failed",
      ");",
      "",
      "StanFitFunctionIntegrate::unkff = \"Unknown FF model `1`. Expected BSZ, BGL, or HQET.\";",
      "",
      "Options[StanGenerateFitFunctions] = {",
      "  \"FFModel\" -> \"BGL\",",
      "  \"FFnf\" -> 2,",
      "  \"NPModel\" -> \"SM\",",
      "  \"HQETModel\" -> \"3/2/1\",",
      "  \"CutControl\" -> 10^-4,",
      "  \"InputDirectory\" -> \"Stan_inputs\"",
      "};",
      "",
      "StanGenerateFitFunctions[OptionsPattern[]] := Module[{ffModel},",
      "  ffModel = OptionValue[\"FFModel\"];",
      "  If[! DirectoryQ[OptionValue[\"InputDirectory\"]],",
      "    CreateDirectory[OptionValue[\"InputDirectory\"], CreateIntermediateDirectories -> True]",
      "  ];",
      "  StanFitFunctionSetup[];",
      "  If[ffModel === \"HQET\", StanFitFunctionUBSetup[]];",
      "  StanFitFunctionIntegrate[ffModel];",
      "];",
      "",
      "If[Length[Rest[$ScriptCommandLine]] > 0 &&",
      "    StringMatchQ[$InputFileName, ___ ~~ \"Stan_FitFunction.wls\"] &&",
      "    MemberQ[{\"BSZ\", \"BGL\", \"HQET\"}, First[Rest[$ScriptCommandLine]]],",
      "  StanGenerateFitFunctions[\"FFModel\" -> First[Rest[$ScriptCommandLine]]];",
      "];"
    },
    "\n"
  ];

  Export[outFile, header <> body <> "\n", "Text"];
  Print["Exported ", outFile];
];
