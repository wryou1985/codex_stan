#!/usr/bin/env wolframscript
(* Extract data-point setup from StanCodeGeneratorBasis.nb into a standalone .wls. *)

ClearAll[heldText, heldBodyText, shiftedInputText, kinematicsText, hpqcdFitText, belleFitText, associationText, apiText];

heldText[expr_] := ToString[expr, InputForm, PageWidth -> Infinity];

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

kinematicsText[] := StringRiffle[
  {
    "mB = 5.27966;",
    "mD = 1.86966;",
    "mDst = 2.01026;",
    "qSQfuncDst[w_] := mB^2 + mDst^2 - 2*mB*mDst*w;",
    "qSQfuncD[w_] := mB^2 + mD^2 - 2*mB*mD*w;",
    "wMax[\"D\"] = w /. NSolve[qSQfuncD[w] == 0, w][[1]];",
    "wMax[\"Dst\"] = w /. NSolve[qSQfuncDst[w] == 0, w][[1]];"
  },
  "\n"
];

hpqcdFitText[] := StringRiffle[
  {
    "StanLoadHPQCD23BGLFit[fitFile_] := Module[{file, importResult, modelParameter, fitResult, nList, fitResultList, bglFFList},",
    "  file = fitFile;",
    "  If[file === Automatic,",
    "    file = FileNameJoin[{ParentDirectory[Directory[]], \"mycode_FitAnalysisCode_v1.2\", \"Stan_output\", \"StanFit_HPQCD23_full_1001_1917.csv\"}]",
    "  ];",
    "  If[! FileExistsQ[file],",
    "    Print[\"Warning: HPQCD23 BGL fit CSV was not found: \", file];",
    "    Return[$Failed]",
    "  ];",
    "  importResult = Import[file];",
    "  modelParameter = importResult[[1]];",
    "  fitResult = Delete[importResult, 1];",
    "  nList = Length[fitResult];",
    "  Do[fitResultList[modelParameter[[i]]] = Table[fitResult[[n, i]], {n, nList}], {i, Length[modelParameter]}];",
    "  bglFFList = modelParameter[[36 ;; 70]];",
    "  bglHPQCD23cent = Table[Mean[fitResultList[para]], {para, bglFFList}];",
    "  bglHPQCD23centNoTensor = Table[Mean[fitResultList[para]], {para, bglFFList[[1 ;; 20]]}];",
    "  bglHPQCD23cov = Table[Covariance[fitResultList[para1], fitResultList[para2]], {para1, bglFFList}, {para2, bglFFList}];",
    "  bglHPQCD23covNoTensor = Table[Covariance[fitResultList[para1], fitResultList[para2]], {para1, bglFFList[[1 ;; 20]]}, {para2, bglFFList[[1 ;; 20]]}];",
    "  True",
    "];"
  },
  "\n"
];

belleFitText[] := StringRiffle[
  {
    "StanBelleFitOptions[] := {\"sc\", \"re_small\", \"re_large\", \"de\", \"pdg\", \"rms\"};",
    "",
    "StanBlockDiagonalMatrix[mats_List] := Module[{n},",
    "  n = Length[mats];",
    "  ArrayFlatten[",
    "    Table[",
    "      If[i == j,",
    "        mats[[i]],",
    "        ConstantArray[0, {Length[mats[[i]]], Length[mats[[j]]]}]",
    "      ],",
    "      {i, n}, {j, n}",
    "    ]",
    "  ]",
    "];",
    "",
    "StanLoadBelleFitInputFile[inputDirectory_, dist_String, option_String] := Module[{file, data, key, covKey},",
    "  file = FileNameJoin[{inputDirectory, \"BelleFit_\" <> dist <> \"_inputs.json\"}];",
    "  If[! FileExistsQ[file],",
    "    Print[\"Warning: BelleFit input JSON was not found: \", file];",
    "    Return[$Failed]",
    "  ];",
    "  data = Import[file, \"RawJSON\"];",
    "  key = \"BelleFit_\" <> dist <> \"_\" <> option;",
    "  covKey = key <> \"_cov\";",
    "  If[! KeyExistsQ[data, key] || ! KeyExistsQ[data, covKey],",
    "    Print[\"Warning: BelleFit option keys were not found in \", file, \". Required keys: \", key, \", \", covKey];",
    "    Return[$Failed]",
    "  ];",
    "  <|\"cent\" -> data[key], \"cov\" -> data[covKey]|>",
    "];",
    "",
    "StanLoadBelleFITData[inputDirectory_, belleFitOption_] := Module[{option, dstDists, dstData, dData, zeroDst, zeroD},",
    "  option = belleFitOption;",
    "  If[! MemberQ[StanBelleFitOptions[], option],",
    "    Print[\"Warning: BelleFitOption is not valid: \", option, \". Available options are: \", StringRiffle[StanBelleFitOptions[], \", \"]];",
    "    Return[$Failed]",
    "  ];",
    "  dstDists = {\"w\", \"cosL\", \"cosV\", \"chi\"};",
    "  dstData = Association@Table[dist -> StanLoadBelleFitInputFile[inputDirectory, dist, option], {dist, dstDists}];",
    "  dData = StanLoadBelleFitInputFile[inputDirectory, \"Dmode\", option];",
    "  If[MemberQ[Values[dstData], $Failed] || dData === $Failed, Return[$Failed]];",
    "  BelleFITcent = Join @@ (dstData[#][\"cent\"] & /@ dstDists);",
    "  BelleFITcov = StanBlockDiagonalMatrix[dstData[#][\"cov\"] & /@ dstDists];",
    "  BelleFITDcent = dData[\"cent\"];",
    "  BelleFITDcov = dData[\"cov\"];",
    "  zeroDst = ConstantArray[0, {Length[BelleFITcent], Length[BelleFITcent]}];",
    "  zeroD = ConstantArray[0, {Length[BelleFITDcent], Length[BelleFITDcent]}];",
    "  BelleSCATTERcov = zeroDst;",
    "  BelleSCATTERDcov = zeroD;",
    "  BelleRawFITcent = BelleFITcent;",
    "  BelleRawFITcov = BelleFITcov;",
    "  BelleRawFITDcent = BelleFITDcent;",
    "  BelleRawFITDcov = BelleFITDcov;",
    "  BellFITdata = <|",
    "    \"BelleFITcent_Dst\" -> BelleFITcent, \"BelleFITcov_Dst\" -> BelleFITcov + BelleSCATTERcov,",
    "    \"BelleFITcent_D\" -> BelleFITDcent, \"BelleFITcov_D\" -> BelleFITDcov + BelleSCATTERDcov,",
    "    \"BelleRawFITcent_Dst\" -> BelleRawFITcent, \"BelleRawFITcov_Dst\" -> BelleRawFITcov,",
    "    \"BelleRawFITcent_D\" -> BelleRawFITDcent, \"BelleRawFITcov_D\" -> BelleRawFITDcov",
    "  |>;",
    "  True",
    "];"
  },
  "\n"
];

associationText[] := StringRiffle[
  {
    "StanBuildDataPointAssociations[] := Module[{},",
    "  MILC15data = <|\"MILC15cent\" -> MILC15cent, \"MILC15cov\" -> MILC15cov|>;",
    "  HPQCD15data = <|\"HPQCD15cent\" -> HPQCD15cent, \"HPQCD15cov\" -> HPQCD15cov|>;",
    "  bglMILC21data = <|\"bglMILC21cent\" -> bglMILC21cent, \"bglMILC21cov\" -> bglMILC21cov|>;",
    "  bglJLQCD23data = <|\"bglJLQCD23cent\" -> bglJLQCD23cent, \"bglJLQCD23cov\" -> bglJLQCD23cov|>;",
    "  If[ValueQ[bglHPQCD23cent] && ValueQ[bglHPQCD23cov],",
    "    bglHPQCD23data = <|\"bglHPQCD23cent\" -> bglHPQCD23cent, \"bglHPQCD23cov\" -> bglHPQCD23cov*maskHPQCD, \"bglHPQCD23centNoTensor\" -> bglHPQCD23centNoTensor, \"bglHPQCD23covNoTensor\" -> bglHPQCD23covNoTensor*maskHPQCD20|>",
    "  ];",
    "  hqetHPQCD23data = <|\"hqetHPQCD23cent\" -> hqetHPQCD23cent, \"hqetHPQCD23cov\" -> hqetHPQCD23cov*maskHPQCD + IdentityMatrix[35]/10^5|>;",
    "  bszLCSR18data = <|\"bszLCSR18cent\" -> bszLCSR18cent, \"bszLCSR18cov\" -> bszLCSR18cov*maskLCSR18, \"bszLCSR18centNoTensor\" -> bszLCSR18centNoTensor, \"bszLCSR18covNoTensor\" -> bszLCSR18covNoTensor*maskLCSR18NoTensor|>;",
    "  bszLCSR23data = <|\"bszLCSR23cent\" -> bszLCSR23cent, \"bszLCSR23cov\" -> bszLCSR23cov*maskLCSR23|>;",
    "  BrRatio2025data = <|\"BrRatio2025cent\" -> BrRatio2025cent, \"BrRatio2025cov\" -> BrRatio2025cov|>;",
    "];"
  },
  "\n"
];

apiText[] := StringRiffle[
  {
    "StanDataPointSymbols[] := {",
    "  \"MILC15data\", \"HPQCD15data\", \"bglMILC21data\", \"bglJLQCD23data\",",
    "  \"bglHPQCD23data\", \"hqetHPQCD23data\", \"bszLCSR18data\",",
    "  \"bszLCSR23data\", \"BrRatio2025data\", \"BellFITdata\"",
    "};",
    "",
    "StanAvailableDataPointSymbols[] := Select[StanDataPointSymbols[], ValueQ[ToExpression[#]] &];",
    "",
    "Options[StanDataPointSetup] = {",
    "  \"HPQCD23FitFile\" -> Automatic,",
    "  \"InputDirectory\" -> \"Stan_inputs\",",
    "  \"BelleFitOption\" -> \"sc\"",
    "};",
    "StanDataPointSetup[OptionsPattern[]] := Module[{},",
    "  StanDataPointKinematics[];",
    "  StanDataPointRawSetup[];",
    "  StanLoadHPQCD23BGLFit[OptionValue[\"HPQCD23FitFile\"]];",
    "  If[StanLoadBelleFITData[OptionValue[\"InputDirectory\"], OptionValue[\"BelleFitOption\"]] === $Failed,",
    "    Return[$Failed]",
    "  ];",
    "  StanBuildDataPointAssociations[];",
    "];",
    "",
    "Options[StanExportDataPointInputs] = {\"OutputDirectory\" -> \"Stan_inputs\", \"FilePrefix\" -> \"standata\"};",
    "StanExportDataPointInputs[OptionsPattern[]] := Module[{outputDir, prefix, symbols, exported},",
    "  outputDir = OptionValue[\"OutputDirectory\"];",
    "  prefix = OptionValue[\"FilePrefix\"];",
    "  If[! DirectoryQ[outputDir], CreateDirectory[outputDir, CreateIntermediateDirectories -> True]];",
    "  symbols = StanAvailableDataPointSymbols[];",
    "  exported = Association@Table[",
    "    With[{name = symbol, file = FileNameJoin[{outputDir, prefix <> \"_\" <> symbol <> \".json\"}]},",
    "      Export[file, ToExpression[symbol], \"JSON\"];",
    "      name -> file",
    "    ],",
    "    {symbol, symbols}",
    "  ];",
    "  Export[FileNameJoin[{outputDir, prefix <> \"_available_data.json\"}], Keys[exported], \"JSON\"];",
    "  exported",
    "];",
    "",
    "Options[StanPrepareDataPointInputs] = Join[Options[StanDataPointSetup], Options[StanExportDataPointInputs]];",
    "StanPrepareDataPointInputs[OptionsPattern[]] := Module[{},",
    "  If[",
    "  StanDataPointSetup[",
    "    \"HPQCD23FitFile\" -> OptionValue[\"HPQCD23FitFile\"],",
    "    \"InputDirectory\" -> OptionValue[\"InputDirectory\"],",
    "    \"BelleFitOption\" -> OptionValue[\"BelleFitOption\"]",
    "  ] === $Failed,",
    "    Return[$Failed]",
    "  ];",
    "  StanExportDataPointInputs[\"OutputDirectory\" -> OptionValue[\"OutputDirectory\"], \"FilePrefix\" -> OptionValue[\"FilePrefix\"]]",
    "];",
    "",
    "If[Length[$ScriptCommandLine] > 0 &&",
    "    StringMatchQ[FileNameTake[First[$ScriptCommandLine]], \"Stan_DataPoint.wls\"],",
    "  StanPrepareDataPointInputs[];",
    "];"
  },
  "\n"
];

Module[{nbFile, outFile, inputs, header, body},
  nbFile = If[Length[$ScriptCommandLine] >= 2, $ScriptCommandLine[[2]], "StanCodeGeneratorBasis.nb"];
  outFile = If[Length[$ScriptCommandLine] >= 3, $ScriptCommandLine[[3]], "Stan_DataPoint.wls"];
  inputs = NotebookImport[nbFile, "Input"];
  header = StringRiffle[
    {
      "#!/usr/bin/env wolframscript",
      "(* ::Package:: *)",
      "",
      "(* Standalone copy of StanCodeGeneratorBasis.nb data-point setup.",
      "   Runtime does not load or reference the notebook.",
      "   Data inputs are exported under Stan_inputs/. *)",
      "",
      "ClearAll[StanDataPointKinematics, StanDataPointRawSetup, StanLoadHPQCD23BGLFit,",
      "  StanBelleFitOptions, StanBlockDiagonalMatrix, StanLoadBelleFitInputFile, StanLoadBelleFITData,",
      "  StanBuildDataPointAssociations, StanDataPointSymbols, StanAvailableDataPointSymbols,",
      "  StanDataPointSetup, StanExportDataPointInputs, StanPrepareDataPointInputs];",
      ""
    },
    "\n"
  ];
  body = StringRiffle[
    {
      "StanDataPointKinematics[] := (",
      kinematicsText[],
      ");",
      "",
      "StanDataPointRawSetup[] := (",
      shiftedInputText[inputs, Join[Range[2, 10], Range[12, 15], Range[17, 27], {16}]],
      ");",
      "",
      hpqcdFitText[],
      "",
      belleFitText[],
      "",
      associationText[],
      "",
      apiText[]
    },
    "\n"
  ];
  Export[outFile, header <> body <> "\n", "Text"];
  Print["Exported ", outFile];
];
