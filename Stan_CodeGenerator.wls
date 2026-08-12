#!/usr/bin/env wolframscript
(* ::Package:: *)

(* Modular Stan code builder for fit procedure (B-PDG).

   Intended use:
   1. Stan_FitAnalysis.py writes a JSON configuration file.
   2. Stan_FitAnalysis.py runs this script with wolframscript.
   3. This script loads Stan_FitFunction.wls, prepares the selected fit
      functions, and exports the .json/.stan/_param.json files for one fit
      model directory.

   This file does not edit or load notebooks at runtime.
*)

ClearAll[
  CodexNL,
  CodexJoinLines,
  CodexStringList,
  CodexRequireSymbols,
  CodexSymbolValue,
  CodexStanVectorAssignments,
  CodexNPParameterBlock,
  CodexNPPriorBlock,
  CodexParamListForPython,
  CodexDefaultSelectedData,
  CodexChiSqDataKeys,
  CodexRequiredDataKeys,
  CodexCasePDGDataModules,
  CodexCasePDGOutputModules,
  CodexResolveOutputKeys,
  CodexBuildStandata,
  CodexModuleDataAssoc,
  CodexBuildStanDataBlock,
  CodexBuildParameterBlock,
  CodexValidUBModels,
  CodexEffectiveUBModel,
  CodexSetPaperPriors,
  CodexUBDataDecl,
  CodexUBFunctionDecl,
  CodexUBFunctionBody,
  CodexUBModelTerm,
  CodexBuildTransformedParameterBlock,
  CodexBuildModelBlock,
  CodexFullChiSqLines,
  CodexBuildGeneratedQuantitiesBlock,
  CodexBuildCasePDGStanCode,
  CodexExportCasePDGStanFit,
  CodexGenerateCasePDGStanFit,
  CodexApplyGeneratorConfig,
  CodexGeneratorCLI
];

CodexNL = "\n";

CodexJoinLines[parts_List] := StringRiffle[DeleteCases[Flatten[parts], ""], CodexNL];

CodexStringList[x_String] := x;
CodexStringList[x_List] := StringJoin[x];
CodexStringList[x_] := ToString[x];

CodexRequireSymbols[symbolNames_List] := Module[{missing},
  missing = Select[symbolNames, ! ValueQ[ToExpression[#]] &];
  If[missing =!= {},
    Message[CodexRequireSymbols::missing, StringRiffle[missing, ", "]];
    Return[$Failed];
  ];
  True
];

CodexRequireSymbols::missing =
  "Required symbols are not defined in the current kernel: `1`. Run the original generator setup first.";

CodexSymbolValue[name_String] := ToExpression[name];

CodexStanVectorAssignments[name_String, values_List] :=
  StringJoin[
    Table[
      name <> "[" <> ToString[i] <> "]=" <> CodexStringList[values[[i]]] <> ";" <> CodexNL,
      {i, Length[values]}
    ]
  ];

CodexNPParameterBlock[] := Module[{},
  If[StringContainsQ[NPmodel, "SM"],
    "",
    If[StringContainsQ[NPmodel, "T"],
      "real<lower=-0.5, upper=0.5> C" <> NPmodel <> ";" <> CodexNL <>
        CodexStringList[ParaDefNP[NPmodel, FFmodel, FFnf]],
      "real<lower=-0.4, upper=0.3> C" <> NPmodel <> ";"
    ]
  ]
];

CodexNPPriorBlock[] := Module[{},
  If[StringContainsQ[NPmodel, "SM"],
    "",
    If[StringContainsQ[NPmodel, "T"],
      "C" <> NPmodel <> " ~ normal(0, 0.5);" <> CodexNL <>
        CodexStringList[ParaPriorNP[NPmodel, FFprior, FFmodel, FFnf]],
      "C" <> NPmodel <> " ~ normal(0, 0.5);" <> CodexNL
    ]
  ]
];

CodexParamListForPython[] := Module[{paramList},
  If[! ValueQ[ParamList],
    If[CodexRequireSymbols[{"FFlist", "NPpara", "NPmodel"}] === $Failed, Return[$Failed]];
    paramList = DeleteDuplicates[Join[FFlist, NPpara[NPmodel]]],
    paramList = ParamList
  ];
  ToString /@ DeleteCases[paramList, Vcb]
];

CodexDefaultSelectedData[] := {
  "MILC15",
  "MILC21",
  "JLQCD23",
  "HPQCD23",
  "LCSR18",
  "LCSR23",
  "BelleFITD",
  "BelleFITDstW",
  "BelleFITDstCosL",
  "BelleFITDstCosV",
  "BelleFITDstChi"
};

CodexChiSqDataKeys[] := {
  "MILC15",
  "MILC21",
  "JLQCD23",
  "HPQCD23",
  "BelleFITD",
  "BelleFITDstW",
  "BelleFITDstCosL",
  "BelleFITDstCosV",
  "BelleFITDstChi",
  "LCSR18",
  "LCSR23"
};

CodexRequiredDataKeys[selectedData_List] :=
  DeleteDuplicates[Join[selectedData, CodexChiSqDataKeys[], {"BrRatio2025"}]];

CodexCasePDGDataModules[] := Module[
  {hpqcdDecl, hpqcdFuncDecl, hpqcdFuncBody, hpqcdModel, hpqcdChi,
   lcsr18Decl, lcsr18FuncDecl, lcsr18FuncBody, lcsr18Model, lcsr18Chi},

  hpqcdDecl = If[ValueQ[HPQCD23data],
    CodexStringList[HPQCD23data],
    If[ValueQ[NPmodel] && NPmodel === "T",
      "vector[35] bglHPQCD23cent;" <> CodexNL <> "cov_matrix[35] bglHPQCD23cov;",
      "vector[20] bglHPQCD23centNoTensor;" <> CodexNL <> "cov_matrix[20] bglHPQCD23covNoTensor;"
    ]
  ];
  hpqcdFuncDecl = "vector[" <> ToString[If[ValueQ[HPQCD23number], HPQCD23number, If[ValueQ[NPmodel] && NPmodel === "T", 35, 20]]] <> "] HPQCD23func;";
  hpqcdFuncBody = If[ValueQ[HPQCD23func],
    CodexStringList[HPQCD23func],
    ""
  ];
  hpqcdModel = If[ValueQ[HPQCD23fit],
    CodexStringList[HPQCD23fit],
    If[ValueQ[NPmodel] && NPmodel === "T",
      "bglHPQCD23cent ~ multi_normal(HPQCD23func, bglHPQCD23cov);",
      "bglHPQCD23centNoTensor ~ multi_normal(HPQCD23func, bglHPQCD23covNoTensor);"
    ]
  ];
  hpqcdChi = If[ValueQ[NPmodel] && NPmodel === "T",
    "dot_product(bglHPQCD23cent - HPQCD23func, inverse(bglHPQCD23cov) * (bglHPQCD23cent - HPQCD23func))",
    "dot_product(bglHPQCD23centNoTensor - HPQCD23func, inverse(bglHPQCD23covNoTensor) * (bglHPQCD23centNoTensor - HPQCD23func))"
  ];

  lcsr18Decl = If[ValueQ[LCSR18data],
    CodexStringList[LCSR18data],
    If[ValueQ[NPmodel] && NPmodel === "T",
      "vector[37] bszLCSR18cent;" <> CodexNL <> "cov_matrix[37] bszLCSR18cov;",
      "vector[22] bszLCSR18centNoTensor;" <> CodexNL <> "cov_matrix[22] bszLCSR18covNoTensor;"
    ]
  ];
  lcsr18FuncDecl = "vector[" <> ToString[If[ValueQ[LCSR18number], LCSR18number, If[ValueQ[NPmodel] && NPmodel === "T", 37, 22]]] <> "] LCSR18func;";
  lcsr18FuncBody = If[ValueQ[LCSR18func],
    CodexStringList[LCSR18func],
    ""
  ];
  lcsr18Model = If[ValueQ[LCSR18fit],
    CodexStringList[LCSR18fit],
    If[ValueQ[NPmodel] && NPmodel === "T",
      "bszLCSR18cent ~ multi_normal(LCSR18func, bszLCSR18cov);",
      "bszLCSR18centNoTensor ~ multi_normal(LCSR18func, bszLCSR18covNoTensor);"
    ]
  ];
  lcsr18Chi = If[ValueQ[NPmodel] && NPmodel === "T",
    "dot_product(bszLCSR18cent - LCSR18func, inverse(bszLCSR18cov) * (bszLCSR18cent - LCSR18func))",
    "dot_product(bszLCSR18centNoTensor - LCSR18func, inverse(bszLCSR18covNoTensor) * (bszLCSR18centNoTensor - LCSR18func))"
  ];

  <|
    "BrRatio2025" -> <|
      "DataSymbols" -> {"BrRatio2025data"},
      "DataDecl" -> "vector[2] BrRatio2025cent;" <> CodexNL <> "cov_matrix[2] BrRatio2025cov;",
      "FunctionDecl" -> "",
      "FunctionBody" -> "",
      "ModelTerm" -> "",
      "ChiSqExpr" -> Missing["NotFit"],
      "RawOutputKey" -> Nothing
    |>,
    "MILC15" -> <|
      "DataSymbols" -> {"MILC15data"},
      "DataDecl" -> "vector[6] MILC15cent;" <> CodexNL <> "cov_matrix[6] MILC15cov;",
      "FunctionDecl" -> "vector[6] MILC15func;",
      "FunctionBody" -> CodexStanVectorAssignments["MILC15func", MILC15funclist],
      "ModelTerm" -> "MILC15cent ~ multi_normal(MILC15func, MILC15cov);",
      "ChiSqExpr" -> "dot_product(MILC15cent - MILC15func, inverse(MILC15cov) * (MILC15cent - MILC15func))"
    |>,
    "MILC21" -> <|
      "DataSymbols" -> {"bglMILC21data"},
      "DataDecl" -> "vector[12] bglMILC21cent;" <> CodexNL <> "cov_matrix[12] bglMILC21cov;",
      "FunctionDecl" -> "vector[12] MILC21func;",
      "FunctionBody" -> CodexStanVectorAssignments["MILC21func", MILC21funclist],
      "ModelTerm" -> "bglMILC21cent ~ multi_normal(MILC21func, bglMILC21cov);",
      "ChiSqExpr" -> "dot_product(bglMILC21cent - MILC21func, inverse(bglMILC21cov) * (bglMILC21cent - MILC21func))"
    |>,
    "JLQCD23" -> <|
      "DataSymbols" -> {"bglJLQCD23data"},
      "DataDecl" -> "vector[12] bglJLQCD23cent;" <> CodexNL <> "cov_matrix[12] bglJLQCD23cov;",
      "FunctionDecl" -> "vector[12] JLQCD23func;",
      "FunctionBody" -> CodexStanVectorAssignments["JLQCD23func", JLQCD23funclist],
      "ModelTerm" -> "bglJLQCD23cent ~ multi_normal(JLQCD23func, bglJLQCD23cov);",
      "ChiSqExpr" -> "dot_product(bglJLQCD23cent - JLQCD23func, inverse(bglJLQCD23cov) * (bglJLQCD23cent - JLQCD23func))"
    |>,
    "HPQCD23" -> <|
      "DataSymbols" -> {"bglHPQCD23cent", "bglHPQCD23cov", "bglHPQCD23centNoTensor", "bglHPQCD23covNoTensor"},
      "DataAssoc" -> If[ValueQ[NPmodel] && NPmodel === "T",
        <|"bglHPQCD23cent" -> bglHPQCD23cent, "bglHPQCD23cov" -> bglHPQCD23cov|>,
        <|"bglHPQCD23centNoTensor" -> bglHPQCD23centNoTensor, "bglHPQCD23covNoTensor" -> bglHPQCD23covNoTensor|>
      ],
      "DataDecl" -> hpqcdDecl,
      "FunctionDecl" -> hpqcdFuncDecl,
      "FunctionBody" -> hpqcdFuncBody,
      "ModelTerm" -> hpqcdModel,
      "ChiSqExpr" -> hpqcdChi
    |>,
    "LCSR18" -> <|
      "DataSymbols" -> {"bszLCSR18cent", "bszLCSR18cov", "bszLCSR18centNoTensor", "bszLCSR18covNoTensor"},
      "DataAssoc" -> If[ValueQ[NPmodel] && NPmodel === "T",
        <|"bszLCSR18cent" -> bszLCSR18cent, "bszLCSR18cov" -> bszLCSR18cov|>,
        <|"bszLCSR18centNoTensor" -> bszLCSR18centNoTensor, "bszLCSR18covNoTensor" -> bszLCSR18covNoTensor|>
      ],
      "DataDecl" -> lcsr18Decl,
      "FunctionDecl" -> lcsr18FuncDecl,
      "FunctionBody" -> lcsr18FuncBody,
      "ModelTerm" -> lcsr18Model,
      "ChiSqExpr" -> lcsr18Chi
    |>,
    "LCSR23" -> <|
      "DataSymbols" -> {"bszLCSR23data"},
      "DataDecl" -> "vector[34] bszLCSR23cent;" <> CodexNL <> "cov_matrix[34] bszLCSR23cov;",
      "FunctionDecl" -> "vector[34] LCSR23func;",
      "FunctionBody" -> CodexStanVectorAssignments["LCSR23func", LCSR23funclist],
      "ModelTerm" -> "bszLCSR23cent ~ multi_normal(LCSR23func, bszLCSR23cov);",
      "ChiSqExpr" -> "dot_product(bszLCSR23cent - LCSR23func, inverse(bszLCSR23cov) * (bszLCSR23cent - LCSR23func))"
    |>,
    "BelleFITD" -> <|
      "DataSymbols" -> {"BelleFITDcent", "BelleFITDcov", "BelleSCATTERDcov"},
      "DataAssoc" -> <|"BelleFITcent_D" -> BelleFITDcent, "BelleFITcov_D" -> BelleFITDcov + BelleSCATTERDcov|>,
      "DataDecl" -> "vector[10] BelleFITcent_D;" <> CodexNL <> "cov_matrix[10] BelleFITcov_D;",
      "FunctionDecl" -> "vector[10] BelleD_func;",
      "FunctionBody" -> StringJoin[Table["BelleD_func[" <> ToString[i] <> "]=(" <> BelleFITfunclist["D"][[i]] <> ")/(" <> BelleFITfunctotal["D"] <> ");" <> CodexNL, {i, 10}]],
      "ModelTerm" -> "BelleFITcent_D ~ multi_normal(BelleD_func, BelleFITcov_D);",
      "ChiSqExpr" -> "dot_product(BelleFITcent_D - BelleD_func, inverse(BelleFITcov_D) * (BelleFITcent_D - BelleD_func))",
      "RawOutputKey" -> "RawBelleD"
    |>,
    "BelleFITDstW" -> <|
      "DataSymbols" -> {"BelleFITcent", "BelleFITcov", "BelleSCATTERcov"},
      "DataAssoc" -> <|
        "BelleFITcent_w" -> BelleFITcent[[1 ;; 10]],
        "BelleFITcov_w" -> (BelleFITcov + BelleSCATTERcov)[[1 ;; 10, 1 ;; 10]]
      |>,
      "DataDecl" -> "vector[10] BelleFITcent_w;" <> CodexNL <> "cov_matrix[10] BelleFITcov_w;",
      "FunctionDecl" -> "vector[10] BelleDst_w_func;",
      "FunctionBody" -> StringJoin[Table["BelleDst_w_func[" <> ToString[i] <> "]=(" <> BelleFITfunclist["w"][[i]] <> ")/(" <> BelleFITfunctotal["Dst"] <> ");" <> CodexNL, {i, 10}]],
      "ModelTerm" -> "BelleFITcent_w ~ multi_normal(BelleDst_w_func, BelleFITcov_w);",
      "ChiSqExpr" -> "dot_product(BelleFITcent_w - BelleDst_w_func, inverse(BelleFITcov_w) * (BelleFITcent_w - BelleDst_w_func))",
      "RawOutputKey" -> "RawBelleDstW"
    |>,
    "BelleFITDstCosL" -> <|
      "DataSymbols" -> {"BelleFITcent", "BelleFITcov", "BelleSCATTERcov"},
      "DataAssoc" -> <|
        "BelleFITcent_cosL" -> BelleFITcent[[11 ;; 20]],
        "BelleFITcov_cosL" -> (BelleFITcov + BelleSCATTERcov)[[11 ;; 20, 11 ;; 20]]
      |>,
      "DataDecl" -> "vector[10] BelleFITcent_cosL;" <> CodexNL <> "cov_matrix[10] BelleFITcov_cosL;",
      "FunctionDecl" -> "vector[10] BelleDst_cosL_func;",
      "FunctionBody" -> StringJoin[Table["BelleDst_cosL_func[" <> ToString[i] <> "]=(" <> BelleFITfunclist["cosL"][[i]] <> ")/(" <> BelleFITfunctotal["Dst"] <> ");" <> CodexNL, {i, 10}]],
      "ModelTerm" -> "BelleFITcent_cosL ~ multi_normal(BelleDst_cosL_func, BelleFITcov_cosL);",
      "ChiSqExpr" -> "dot_product(BelleFITcent_cosL - BelleDst_cosL_func, inverse(BelleFITcov_cosL) * (BelleFITcent_cosL - BelleDst_cosL_func))",
      "RawOutputKey" -> "RawBelleDstCosL"
    |>,
    "BelleFITDstCosV" -> <|
      "DataSymbols" -> {"BelleFITcent", "BelleFITcov", "BelleSCATTERcov"},
      "DataAssoc" -> <|
        "BelleFITcent_cosV" -> BelleFITcent[[21 ;; 30]],
        "BelleFITcov_cosV" -> (BelleFITcov + BelleSCATTERcov)[[21 ;; 30, 21 ;; 30]]
      |>,
      "DataDecl" -> "vector[10] BelleFITcent_cosV;" <> CodexNL <> "cov_matrix[10] BelleFITcov_cosV;",
      "FunctionDecl" -> "vector[10] BelleDst_cosV_func;",
      "FunctionBody" -> StringJoin[Table["BelleDst_cosV_func[" <> ToString[i] <> "]=(" <> BelleFITfunclist["cosV"][[i]] <> ")/(" <> BelleFITfunctotal["Dst"] <> ");" <> CodexNL, {i, 10}]],
      "ModelTerm" -> "BelleFITcent_cosV ~ multi_normal(BelleDst_cosV_func, BelleFITcov_cosV);",
      "ChiSqExpr" -> "dot_product(BelleFITcent_cosV - BelleDst_cosV_func, inverse(BelleFITcov_cosV) * (BelleFITcent_cosV - BelleDst_cosV_func))",
      "RawOutputKey" -> "RawBelleDstCosV"
    |>,
    "BelleFITDstChi" -> <|
      "DataSymbols" -> {"BelleFITcent", "BelleFITcov", "BelleSCATTERcov"},
      "DataAssoc" -> <|
        "BelleFITcent_chi" -> BelleFITcent[[31 ;; 40]],
        "BelleFITcov_chi" -> (BelleFITcov + BelleSCATTERcov)[[31 ;; 40, 31 ;; 40]]
      |>,
      "DataDecl" -> "vector[10] BelleFITcent_chi;" <> CodexNL <> "cov_matrix[10] BelleFITcov_chi;",
      "FunctionDecl" -> "vector[10] BelleDst_chi_func;",
      "FunctionBody" -> StringJoin[Table["BelleDst_chi_func[" <> ToString[i] <> "]=(" <> BelleFITfunclist["chi"][[i]] <> ")/(" <> BelleFITfunctotal["Dst"] <> ");" <> CodexNL, {i, 10}]],
      "ModelTerm" -> "BelleFITcent_chi ~ multi_normal(BelleDst_chi_func, BelleFITcov_chi);",
      "ChiSqExpr" -> "dot_product(BelleFITcent_chi - BelleDst_chi_func, inverse(BelleFITcov_chi) * (BelleFITcent_chi - BelleDst_chi_func))",
      "RawOutputKey" -> "RawBelleDstChi"
    |>
  |>
];

CodexCasePDGOutputModules[] := <|
  "VcbFromBr" -> <|
    "Decl" -> CodexJoinLines[{
      "vector[2] rawBrRatio_func;",
      "vector[2] BrRatio_func;",
      "vector[2] BrRatio_data;",
      "vector[2] Vcb_Br;",
      "real Vcb_mean;"
    }],
    "Body" -> StringJoin[
      Table["rawBrRatio_func[" <> ToString[i] <> "]=" <> BelleFITfunclist["Br"][[i]] <> ";" <> CodexNL, {i, 2}],
      CodexNL,
      "BrRatio_data = multi_normal_rng(BrRatio2025cent , BrRatio2025cov);" <> CodexNL,
      "Vcb_Br = sqrt(BrRatio_data ./ (rawBrRatio_func));" <> CodexNL,
      "Vcb_mean = mean(sqrt(BrRatio2025cent ./ (rawBrRatio_func)));" <> CodexNL,
      "BrRatio_func = rawBrRatio_func * square(Vcb_mean);" <> CodexNL
    ]
  |>,
  "RD" -> <|
    "Decl" -> "vector[2] RD_func;",
    "Body" -> StringJoin[
      Table[
        "RD_func[" <> ToString[i] <> "]= (" <> BelleFITfunclist["BrTau"][[i]] <> ") / (" <> BelleFITfunclist["Br"][[i]] <> ");" <> CodexNL,
        {i, 2}
      ]
    ]
  |>,
  "RawBelleD" -> <|
    "Decl" -> "vector[10] rawBelleD_func;",
    "Body" -> StringJoin[Table["rawBelleD_func[" <> ToString[i] <> "]= " <> BelleFITfunclist["D"][[i]] <> ";" <> CodexNL, {i, 10}]]
  |>,
  "RawBelleDstW" -> <|
    "Decl" -> "vector[10] rawBelleDst_w_func;",
    "Body" -> StringJoin[Table["rawBelleDst_w_func[" <> ToString[i] <> "]= " <> BelleFITfunclist["w"][[i]] <> ";" <> CodexNL, {i, 10}]]
  |>,
  "RawBelleDstCosL" -> <|
    "Decl" -> "vector[10] rawBelleDst_cosL_func;",
    "Body" -> StringJoin[Table["rawBelleDst_cosL_func[" <> ToString[i] <> "]= " <> BelleFITfunclist["cosL"][[i]] <> ";" <> CodexNL, {i, 10}]]
  |>,
  "RawBelleDstCosV" -> <|
    "Decl" -> "vector[10] rawBelleDst_cosV_func;",
    "Body" -> StringJoin[Table["rawBelleDst_cosV_func[" <> ToString[i] <> "]= " <> BelleFITfunclist["cosV"][[i]] <> ";" <> CodexNL, {i, 10}]]
  |>,
  "RawBelleDstChi" -> <|
    "Decl" -> "vector[10] rawBelleDst_chi_func;",
    "Body" -> StringJoin[Table["rawBelleDst_chi_func[" <> ToString[i] <> "]= " <> BelleFITfunclist["chi"][[i]] <> ";" <> CodexNL, {i, 10}]]
  |>
|>;

CodexResolveOutputKeys[selectedData_List, outputSpec_] := Module[
  {modules, rawKeys, base},
  modules = CodexCasePDGDataModules[];
  rawKeys = DeleteMissing[Lookup[modules /@ selectedData, "RawOutputKey", Missing[]]];
  base = If[outputSpec === Automatic,
    Join[{"VcbFromBr", "RD"}, rawKeys],
    outputSpec
  ];
  DeleteDuplicates[Join[{"VcbFromBr"}, base]]
];

CodexModuleDataAssoc[module_Association] := Module[{dataAssoc},
  dataAssoc = Lookup[module, "DataAssoc", Missing["UseDataSymbols"]];
  If[dataAssoc === Missing["UseDataSymbols"],
    Join @@ (CodexSymbolValue /@ Lookup[module, "DataSymbols"]),
    dataAssoc
  ]
];

CodexBuildStandata[selectedData_List] := Module[
  {modules, selectedModules, allData, dataSymbols},
  modules = CodexCasePDGDataModules[];
  allData = CodexRequiredDataKeys[selectedData];
  selectedModules = modules /@ allData;
  dataSymbols = DeleteDuplicates[Flatten[Lookup[selectedModules, "DataSymbols"]]];
  If[CodexRequireSymbols[dataSymbols] === $Failed, Return[$Failed]];
  Join @@ (CodexModuleDataAssoc /@ selectedModules)
];

CodexBuildStanDataBlock[selectedData_List] := Module[
  {modules, allData, ubDecl},
  modules = CodexCasePDGDataModules[];
  allData = CodexRequiredDataKeys[selectedData];
  ubDecl = CodexUBDataDecl[CodexEffectiveUBModel[If[ValueQ[UBModel], UBModel, "no-UB"]]];
  "data {" <> CodexNL <>
    CodexJoinLines[Join[Lookup[modules /@ allData, "DataDecl"], {ubDecl}]] <>
    CodexNL <> "}"
];

CodexValidUBModels[] := {"no-UB", "soft-UB", "hard-UB"};

CodexEffectiveUBModel[ubModel_String] := Module[{},
  If[ValueQ[FFmodel] && FFmodel =!= "HQET", Return["no-UB"]];
  ubModel
];

CodexSetPaperPriors[] := Module[{},
  If[! ValueQ[FFmodel], Return[$Failed]];
  FFprior = Switch[
    FFmodel,
    "BSZ", "wide2",
    "BGL", "simple",
    "HQET", "QCDSR",
    _, Message[CodexSetPaperPriors::unkff, FFmodel]; Return[$Failed]
  ];
  True
];

CodexSetPaperPriors::unkff =
  "Unknown FFmodel `1`. Expected BSZ, BGL, or HQET.";

CodexUBDataDecl[ubModel_String] := Module[{},
  Switch[
    ubModel,
    "no-UB", "",
    "soft-UB", "vector[4] UBmax;" <> CodexNL <> "vector<lower=0>[4] UB_sigma;",
    "hard-UB", "vector[4] UBmax;",
    _, Message[CodexUBModel::unk, ubModel, StringRiffle[CodexValidUBModels[], ", "]]; $Failed
  ]
];

CodexUBModel::unk =
  "Unknown UBModel `1`. Available UB models are: `2`.";

CodexUBFunctionDecl[] := Module[{},
  If[ValueQ[UBdef] && StringLength[StringTrim[CodexStringList[UBdef]]] > 0,
    CodexStringList[UBdef],
    If[ValueQ[UBfunc] && StringLength[StringTrim[CodexStringList[UBfunc]]] > 0,
      "vector[4] UBfunc;",
      ""
    ]
  ]
];

CodexUBFunctionBody[] := Module[{},
  If[ValueQ[UBfunc], CodexStringList[UBfunc], ""]
];

CodexUBModelTerm[ubModel_String] := Module[{},
  If[! ValueQ[UBfunc] || StringLength[StringTrim[CodexStringList[UBfunc]]] == 0,
    Return[""]
  ];
  Switch[
    ubModel,
    "no-UB", "",
    "soft-UB",
      CodexJoinLines[{
        "for (k in 1:4) {",
        "  target += normal_lccdf(UBfunc[k] | UBmax[k], UB_sigma[k]);",
        "}"
      }],
    "hard-UB",
      CodexJoinLines[{
        "for (k in 1:4) {",
        "  if (UBfunc[k] > UBmax[k]) {",
        "    target += negative_infinity();",
        "  }",
        "}"
      }],
    _, Message[CodexUBModel::unk, ubModel, StringRiffle[CodexValidUBModels[], ", "]]; $Failed
  ]
];

CodexBuildParameterBlock[] := Module[{},
  If[! ValueQ[FFrange], FFrange = "full"];
  If[CodexRequireSymbols[{"NPmodel", "FFmodel", "FFnf", "FFrange", "ParaDef"}] === $Failed, Return[$Failed]];
  "parameters {" <> CodexNL <>
    CodexNPParameterBlock[] <>
    CodexStringList[ParaDef[FFrange, FFmodel, FFnf]] <>
    CodexNL <> "}"
];

CodexBuildTransformedParameterBlock[selectedData_List] := Module[
  {modules, functionData, decls, bodies, ubDecl, ubBody},
  modules = CodexCasePDGDataModules[];
  functionData = DeleteCases[CodexRequiredDataKeys[selectedData], "BrRatio2025"];
  decls = Lookup[modules /@ functionData, "FunctionDecl"];
  bodies = Lookup[modules /@ functionData, "FunctionBody"];
  ubDecl = CodexUBFunctionDecl[];
  ubBody = CodexUBFunctionBody[];
  "transformed parameters {" <> CodexNL <>
    CodexJoinLines[Join[decls, {ubDecl}]] <> CodexNL <>
    CodexStringList[StringJoin[Join[bodies, {ubBody}]]] <>
    "}"
];

CodexBuildModelBlock[selectedData_List] := Module[
  {modules, modelTerms, ubTerm},
  If[CodexSetPaperPriors[] === $Failed, Return[$Failed]];
  If[CodexRequireSymbols[{"FFprior", "FFmodel", "FFnf", "ParaPrior"}] === $Failed, Return[$Failed]];
  modules = CodexCasePDGDataModules[];
  modelTerms = DeleteCases[Lookup[modules /@ selectedData, "ModelTerm"], ""];
  ubTerm = CodexUBModelTerm[CodexEffectiveUBModel[If[ValueQ[UBModel], UBModel, "no-UB"]]];
  If[ubTerm === $Failed, Return[$Failed]];
  "model {" <> CodexNL <>
    CodexNPPriorBlock[] <>
    CodexStringList[ParaPrior[FFprior, FFmodel, FFnf]] <> CodexNL <>
    CodexJoinLines[Join[modelTerms, {ubTerm}]] <>
    CodexNL <> "}"
];

CodexFullChiSqLines[] := Module[
  {modules, lcsr18Diag, lcsr23Diag},
  modules = CodexCasePDGDataModules[];
  lcsr18Diag = If[ValueQ[NPmodel] && NPmodel === "T",
    "dot_self((bszLCSR18cent - LCSR18func) ./ sqrt(diagonal(bszLCSR18cov)))",
    "dot_self((bszLCSR18centNoTensor - LCSR18func) ./ sqrt(diagonal(bszLCSR18covNoTensor)))"
  ];
  lcsr23Diag = "dot_self((bszLCSR23cent - LCSR23func) ./ sqrt(diagonal(bszLCSR23cov)))";
  {
    "chi_sq[1] = " <> modules["MILC15"]["ChiSqExpr"] <> ";",
    "chi_sq[2] = " <> modules["MILC21"]["ChiSqExpr"] <> ";",
    "chi_sq[3] = " <> modules["JLQCD23"]["ChiSqExpr"] <> ";",
    "chi_sq[4] = " <> modules["HPQCD23"]["ChiSqExpr"] <> ";",
    "chi_sq[5] = " <> modules["BelleFITD"]["ChiSqExpr"] <> ";",
    "chi_sq[6] = " <> modules["BelleFITDstW"]["ChiSqExpr"] <> ";",
    "chi_sq[7] = " <> modules["BelleFITDstCosL"]["ChiSqExpr"] <> ";",
    "chi_sq[8] = " <> modules["BelleFITDstCosV"]["ChiSqExpr"] <> ";",
    "chi_sq[9] = " <> modules["BelleFITDstChi"]["ChiSqExpr"] <> ";",
    "chi_sq[10] = dot_product(BrRatio2025cent - BrRatio_func, inverse(BrRatio2025cov) * (BrRatio2025cent - BrRatio_func));",
    "chi_sq[11] = " <> modules["LCSR18"]["ChiSqExpr"] <> ";",
    "chi_sq[12] = " <> modules["LCSR23"]["ChiSqExpr"] <> ";",
    "chi_sq[13] = " <> lcsr18Diag <> ";",
    "chi_sq[14] = " <> lcsr23Diag <> ";",
    "chi_sq_total = sum(chi_sq[1:12]);"
  }
];

CodexBuildGeneratedQuantitiesBlock[selectedData_List, outputSpec_: Automatic] := Module[
  {outputs, outputKeys, outputAssoc, decls, bodies},
  outputs = CodexCasePDGOutputModules[];
  outputKeys = CodexResolveOutputKeys[selectedData, outputSpec];
  outputAssoc = outputs /@ outputKeys;
  decls = Join[
    Lookup[outputAssoc, "Decl"],
    {"vector[14] chi_sq;", "real chi_sq_total;"}
  ];
  bodies = Join[
    Lookup[outputAssoc, "Body"],
    CodexFullChiSqLines[]
  ];
  "generated quantities {" <> CodexNL <>
    CodexJoinLines[decls] <> CodexNL <>
    CodexJoinLines[bodies] <>
    CodexNL <> "}"
];

Options[CodexBuildCasePDGStanCode] = {
  "SelectedData" -> Automatic,
  "GeneratedOutputs" -> Automatic,
  "UBModel" -> Automatic
};

CodexBuildCasePDGStanCode[OptionsPattern[]] := Module[
  {selectedData, generatedOutputs, ubModel, knownKeys, unknown},
  selectedData = Replace[OptionValue["SelectedData"], Automatic :> CodexDefaultSelectedData[]];
  generatedOutputs = OptionValue["GeneratedOutputs"];
  ubModel = CodexEffectiveUBModel[Replace[OptionValue["UBModel"], Automatic :> If[ValueQ[UBModel], UBModel, "no-UB"]]];
  If[! MemberQ[CodexValidUBModels[], ubModel],
    Message[CodexUBModel::unk, ubModel, StringRiffle[CodexValidUBModels[], ", "]];
    Return[$Failed];
  ];
  UBModel = ubModel;
  knownKeys = Keys[CodexCasePDGDataModules[]];
  unknown = Complement[selectedData, knownKeys];
  If[unknown =!= {},
    Message[CodexBuildCasePDGStanCode::unkdata, StringRiffle[unknown, ", "], StringRiffle[knownKeys, ", "]];
    Return[$Failed];
  ];
  StringRiffle[
    {
      CodexBuildStanDataBlock[selectedData],
      CodexBuildParameterBlock[],
      CodexBuildTransformedParameterBlock[selectedData],
      CodexBuildModelBlock[selectedData],
      CodexBuildGeneratedQuantitiesBlock[selectedData, generatedOutputs]
    },
    CodexNL <> CodexNL
  ]
];

CodexBuildCasePDGStanCode::unkdata =
  "Unknown SelectedData key(s): `1`. Available keys are: `2`.";

Options[CodexExportCasePDGStanFit] = {
  "SelectedData" -> Automatic,
  "GeneratedOutputs" -> Automatic,
  "OutputRoot" -> "Stan_inputs",
  "UBModel" -> Automatic
};

CodexExportCasePDGStanFit[analysisName_String, OptionsPattern[]] := Module[
  {selectedData, generatedOutputs, outputRoot, ubModel, fitDir, stanCode, standata, params,
   dataFile, stanFile, paramFile},
  selectedData = Replace[OptionValue["SelectedData"], Automatic :> CodexDefaultSelectedData[]];
  generatedOutputs = OptionValue["GeneratedOutputs"];
  outputRoot = OptionValue["OutputRoot"];
  ubModel = CodexEffectiveUBModel[Replace[OptionValue["UBModel"], Automatic :> If[ValueQ[UBModel], UBModel, "no-UB"]]];
  stanCode = CodexBuildCasePDGStanCode[
    "SelectedData" -> selectedData,
    "GeneratedOutputs" -> generatedOutputs,
    "UBModel" -> ubModel
  ];
  If[stanCode === $Failed, Return[$Failed]];
  standata = CodexBuildStandata[selectedData];
  If[standata === $Failed, Return[$Failed]];
  params = CodexParamListForPython[];
  If[params === $Failed, Return[$Failed]];

  fitDir = FileNameJoin[{outputRoot, analysisName}];
  If[! DirectoryQ[fitDir], CreateDirectory[fitDir, CreateIntermediateDirectories -> True]];
  dataFile = FileNameJoin[{fitDir, "standata_" <> analysisName <> ".json"}];
  stanFile = FileNameJoin[{fitDir, "stancode_" <> analysisName <> ".stan"}];
  paramFile = FileNameJoin[{fitDir, "stancode_" <> analysisName <> "_param.json"}];

  Export[dataFile, standata, "JSON"];
  Export[stanFile, stanCode, "Text"];
  Export[paramFile, params, "JSON"];

  <|
    "AnalysisName" -> analysisName,
    "FitDirectory" -> fitDir,
    "DataFile" -> dataFile,
    "StanFile" -> stanFile,
    "ParamFile" -> paramFile,
    "SelectedData" -> selectedData,
    "GeneratedOutputs" -> CodexResolveOutputKeys[selectedData, generatedOutputs],
    "UBModel" -> ubModel
  |>
];

Options[CodexGenerateCasePDGStanFit] = Options[CodexExportCasePDGStanFit];

CodexGenerateCasePDGStanFit[analysisName_String, opts : OptionsPattern[]] :=
  CodexExportCasePDGStanFit[analysisName, opts];

CodexApplyGeneratorConfig[config_Association] := Module[
  {setupFiles, dataPointFile, ffFunctionFile, selectedData, generatedOutputs},
  SetDirectory[Lookup[config, "WorkingDirectory", Directory[]]];
  setupFiles = Lookup[config, "SetupFiles", {}];
  Scan[Get, setupFiles];

  dataPointFile = Lookup[config, "DataPointFile", "Stan_DataPoint.wls"];
  If[FileExistsQ[dataPointFile], Get[dataPointFile]];

  ffFunctionFile = Lookup[config, "FitFunctionFile", "Stan_FitFunction.wls"];
  If[FileExistsQ[ffFunctionFile], Get[ffFunctionFile]];

  NPmodel = Lookup[config, "NPModel", If[ValueQ[NPmodel], NPmodel, "SM"]];
  FFmodel = Lookup[config, "FFModel", If[ValueQ[FFmodel], FFmodel, "BGL"]];
  FFnf = Lookup[config, "FFnf", If[ValueQ[FFnf], FFnf, 2]];
  FFrange = Lookup[config, "FFrange", If[ValueQ[FFrange], FFrange, "full"]];
  UBModel = Lookup[config, "UBModel", If[ValueQ[UBModel], UBModel, "no-UB"]];
  CutControl = Lookup[config, "CutControl", If[ValueQ[CutControl], CutControl, 10^-4]];

  If[Length[DownValues[StanPrepareDataPointInputs]] > 0,
    StanPrepareDataPointInputs[
      "OutputDirectory" -> Lookup[config, "InputDirectory", "Stan_inputs"],
      "FilePrefix" -> Lookup[config, "DataPointFilePrefix", "standata"],
      "InputDirectory" -> Lookup[config, "InputDirectory", "Stan_inputs"],
      "HPQCD23FitFile" -> Lookup[config, "HPQCD23FitFile", Automatic],
      "BelleFitOption" -> Lookup[config, "BelleFitOption", "sc"]
    ];
  ];

  If[Length[DownValues[StanGenerateFitFunctions]] > 0,
    StanGenerateFitFunctions[
      "FFModel" -> FFmodel,
      "FFnf" -> FFnf,
      "NPModel" -> NPmodel,
      "HQETModel" -> Lookup[config, "HQETModel", Missing["NotUsed"]],
      "CutControl" -> CutControl,
      "InputDirectory" -> Lookup[config, "InputDirectory", "Stan_inputs"]
    ];
  ];

  selectedData = Replace[Lookup[config, "SelectedData", Automatic], "Automatic" -> Automatic];
  generatedOutputs = Replace[Lookup[config, "GeneratedOutputs", Automatic], "Automatic" -> Automatic];

  CodexGenerateCasePDGStanFit[
    Lookup[config, "AnalysisName"],
    "SelectedData" -> selectedData,
    "GeneratedOutputs" -> generatedOutputs,
    "OutputRoot" -> Lookup[config, "OutputRoot", "Stan_inputs"],
    "UBModel" -> UBModel
  ]
];

CodexGeneratorCLI[] := Module[{args, config, result},
  args = Rest[$ScriptCommandLine];
  If[Length[args] == 0, Return[Null]];
  config = Import[First[args], "RawJSON"];
  If[! AssociationQ[config],
    Print["Generator config must be a JSON object: ", First[args]];
    Exit[2];
  ];
  result = CodexApplyGeneratorConfig[config];
  If[result === $Failed,
    Print["Stan generation failed."];
    Exit[1];
  ];
  Print[ExportString[result, "JSON"]];
];

If[
  StringMatchQ[$InputFileName, ___ ~~ "Stan_CodeGenerator.wls"] &&
    Length[Rest[$ScriptCommandLine]] > 0 &&
    FileExistsQ[First[Rest[$ScriptCommandLine]]] &&
    StringEndsQ[First[Rest[$ScriptCommandLine]], ".json"],
  CodexGeneratorCLI[];
];
