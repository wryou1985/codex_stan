#!/usr/bin/env wolframscript
(* ::Package:: *)

(* Standalone copy of StanCodeGeneratorBasis.nb Run2 and Run3.
   Runtime does not load or reference the notebook.
   Integration outputs are read/written under Stan_inputs/. *)

ClearAll[
  StanFitFunctionSetup,
  StanFitFunctionUBSetup,
  StanFitFunctionIntegrate,
  StanFitFunctionCommonSetup,
  StanPrepareFitFunctionLists,
  StanGenerateFitFunctions
];

StanFitFunctionSetup[] := (
GF = 1.166371/10^5;
mb = 4.18;
mc = 0.915;
mtau = 1.77686;
mB = 5.27966;
mBst = 5.32470;
mD = 1.86966;
mDst = 2.01026;
qDstMAX = (mB - mDst)^2;
rD = mD/mB;
rDst = mDst/mB;
QpQmD[qSQ_] := ((mB + mD)^2 - qSQ)*((mB - mD)^2 - qSQ);
QpQmDst[qSQ_] := ((mB + mDst)^2 - qSQ)*((mB - mDst)^2 - qSQ);
qSQfuncDst[w_] := mB^2 + mDst^2 - 2*mB*mDst*w;
qSQfuncD[w_] := mB^2 + mD^2 - 2*mB*mD*w;
wFuncDst[qSQ_] := (mB^2 + mDst^2 - qSQ)/(2*mB*mDst);
wFuncD[qSQ_] := (mB^2 + mD^2 - qSQ)/(2*mB*mD);
wMax["D"] = w /. NSolve[qSQfuncD[w] == 0, w][[1]];
wMax["Dst"] = w /. NSolve[qSQfuncDst[w] == 0, w][[1]];
wMaxTau["D"] = 1.429071189096182;
wMaxTau["Dst"] = 1.3548194509494718;
BrDstDpi = 0.677;
SecondToGeV = 1/(6.582119569*10^-25);
LifeTime["Bp"] = 1.638*10^-12*SecondToGeV;
(FFc["BSZ_fp"] = {afp0, afp1, afp2, afp3}; ); Null; (FFc["BSZ_f0"] = {af00, af01, af02, af03}; ); Null; (FFc["BSZ_fT"] = {afT0, afT1, afT2, afT3}; ); Null; (FFc["BSZ_A1"] = {aA10, aA11, aA12, aA13}; ); Null; (FFc["BSZ_A12"] = {aA120, aA121, aA122, aA123}; ); Null; (FFc["BSZ_A0"] = {aA00, aA01, aA02, aA03}; ); Null; (FFc["BSZ_V"] = {aV0, aV1, aV2, aV3}; ); Null; (FFc["BSZ_T1"] = {aT10, aT11, aT12, aT13}; ); Null; (FFc["BSZ_T2"] = {aT20, aT21, aT22, aT23}; ); Null; (FFc["BSZ_T23"] = {aT230, aT231, aT232, aT233}; ); Null; (MR["BSZ_A0"] = 6.275; ); Null; (MR["BSZ_fp"] = MR["BSZ_fT"] = MR["BSZ_V"] = MR["BSZ_T1"] = 6.33; ); Null; (MR["BSZ_f0"] = 6.42; ); Null; (MR["BSZ_A1"] = MR["BSZ_A12"] = MR["BSZ_T2"] = MR["BSZ_T23"] = 6.767; ); Null; (tpD = (mB + mD)^2; ); Null; (t0DBSZ = (mB + mD)*(Sqrt[mB] - Sqrt[mD])^2; ); Null; (zD["BSZ", qSQ_] = (Sqrt[tpD - qSQ] - Sqrt[tpD - t0DBSZ])/(Sqrt[tpD - qSQ] + Sqrt[tpD - t0DBSZ]); ); Null; (tpDst = (mB + mDst)^2; ); Null; (t0DstBSZ = (mB + mDst)*(Sqrt[mB] - Sqrt[mDst])^2; ); Null; (zDst["BSZ", qSQ_] = (Sqrt[tpDst - qSQ] - Sqrt[tpDst - t0DstBSZ])/(Sqrt[tpDst - qSQ] + Sqrt[tpDst - t0DstBSZ]); ); Null; zBSZ[FF_, qSQ_] := If[StringContainsQ[FF, "f"], zD["BSZ", qSQ], zDst["BSZ", qSQ]]; Null; BSZff[FFmode_, NF_, qSQ_] := (1/(1 - qSQ/MR[FFmode]^2))*Sum[FFc[FFmode][[n + 1]]*(zBSZ[FFmode, qSQ] - zBSZ[FFmode, 0])^n, {n, 0, NF}]; Null; BSZff["BSZ_A2", NF_, qSQ_] := (1/QpQmDst[qSQ])*((mB + mDst)^2*(mB^2 - mDst^2 - qSQ)*BGLff["BSZ_A1", NF, qSQ] - 16*mB*mDst^2*(mB + mDst)*BGLff["BSZ_A12", NF, qSQ]); Null; BSZff["BSZ_fS", NF_, qSQ_] := ((mB - mD)/(mb - mc))*BSZff["BSZ_f0", NF, qSQ]; Null; BSZff["BSZ_P", NF_, qSQ_] := (-((2*mDst)/(mb + mc)))*BSZff["BSZ_A0", NF, qSQ]; Null; (KineticBSZ = Join[Join[{af00 -> afp0, aA00 -> 3.562520345*aA120, aT20 -> aT10}, Flatten[NSolve[BSZff["BSZ_A1", 2, qDstMAX] == ((8*mB*mDst)/(mB^2 - mDst^2))*BSZff["BSZ_A12", 2, qDstMAX], aA120]]], Flatten[NSolve[BSZff["BSZ_T2", 2, qDstMAX] == ((4*mB*mDst)/(mB + mDst)^2)*BSZff["BSZ_T23", 2, qDstMAX], aT230]]]; )

; 

(BSZff["BSZ_f", NF_, qSQ_] := (mB + mDst)*BSZff["BSZ_A1", NF, qSQ]; ); Null; (BSZff["BSZ_F1", NF_, qSQ_] := 8*mB*mDst*BSZff["BSZ_A12", NF, qSQ]; ); Null; (BSZff["BSZ_F2", NF_, qSQ_] := 2*BSZff["BSZ_A0", NF, qSQ]; ); Null; (BSZff["BSZ_g", NF_, qSQ_] := (2/(mB + mDst))*BSZff["BSZ_V", NF, qSQ]; )

; 

(FFc["BGL_fp"] = {bfp0, bfp1, bfp2, bfp3}; ); Null; (FFc["BGL_f0"] = {bf00, bf01, bf02, bf03}; ); Null; (FFc["BGL_fT"] = {bfT0, bfT1, bfT2, bfT3}; ); Null; (FFc["BGL_f"] = {bf0, bf1, bf2, bf3}; ); Null; (FFc["BGL_F1"] = {bF10, bF11, bF12, bF13}; ); Null; (FFc["BGL_F2"] = {bF20, bF21, bF22, bF23}; ); Null; (FFc["BGL_g"] = {bg0, bg1, bg2, bg3}; ); Null; (FFc["BGL_T1"] = {bT10, bT11, bT12, bT13}; ); Null; (FFc["BGL_T2"] = {bT20, bT21, bT22, bT23}; ); Null; (FFc["BGL_T23"] = {bT230, bT231, bT232, bT233}; ); Null; (MR["BGL_fT"] = MR["BGL_T1"] = 6.33; ); Null; (MR["BGL_T2"] = MR["BGL_T23"] = 6.767; ); Null; (zF0 = {-0.433, -0.819}; ); Null; (zFp = {-0.308, -0.555, -0.646}; ); Null; (zF = {-0.401861, -0.405917, -0.636689, -0.641829}; ); Null; (zF1 = {-0.401861, -0.405917, -0.636689, -0.641829}; ); Null; (zF2 = {-0.274321, -0.442859, -0.790626}; ); Null; (zG = {-0.286001, -0.479425, -0.536533, -0.88974}; ); Null; (PF["BGL_fp", z_] := Product[(z - zn)/(1 - z*zn), {zn, zFp}]; ); Null; (PF["BGL_f0", z_] := Product[(z - zn)/(1 - z*zn), {zn, zF0}]; ); Null; (PF["BGL_f", z_] := Product[(z - zn)/(1 - z*zn), {zn, zF}]; ); Null; (PF["BGL_F1", z_] := Product[(z - zn)/(1 - z*zn), {zn, zF1}]; ); Null; (PF["BGL_F2", z_] := Product[(z - zn)/(1 - z*zn), {zn, zF2}]; ); Null; (PF["BGL_g", z_] := Product[(z - zn)/(1 - z*zn), {zn, zG}]; ); Null; (nI = 2.6; ); Null; (χT1m = 5.131/10^4; ); Null; (χT1p = 3.894/10^4; ); Null; (χL1p = 19.421/10^3; ); Null; (ϕF["BGL_fp", z_] = (12.43*(1 + z)^2*(1 - z)^(1/2))/((1 + rD)*(1 - z) + 2*Sqrt[rD]*(1 + z))^5; ); Null; (ϕF["BGL_f0", z_] = (10.11*(1 + z)*(1 - z)^(3/2))/((1 + rD)*(1 - z) + 2*Sqrt[rD]*(1 + z))^4; ); Null; (ϕF["BGL_f", z_] = (((4*rDst)/mB^2)*Sqrt[nI/(3*Pi*χT1p)]*(1 + z)*(1 - z)^(3/2))/((1 + rDst)*(1 - z) + 2*Sqrt[rDst]*(1 + z))^4; ); Null; (ϕF["BGL_F1", z_] = (((4*rDst)/mB^3)*Sqrt[nI/(6*Pi*χT1p)]*(1 + z)*(1 - z)^(5/2))/((1 + rDst)*(1 - z) + 2*Sqrt[rDst]*(1 + z))^5; ); Null; (ϕF["BGL_F2", z_] = (Sqrt[nI/(Pi*χL1p)]*8*Sqrt[2]*rDst^2*(1 + z)^2)/((1 - z)^2^(-1)*((1 + rDst)*(1 - z) + 2*Sqrt[rDst]*(1 + z))^4); ); Null; (ϕF["BGL_g", z_] = (2^4*Sqrt[nI/(3*Pi*χT1m)]*rDst^2*(1 + z)^2)/((1 - z)^2^(-1)*((1 + rDst)*(1 - z) + 2*Sqrt[rDst]*(1 + z))^4); ); Null; (tpD = (mB + mD)^2; ); Null; (t0DBGL = (mB - mD)^2; ); Null; (zD["BGL", qSQ_] = (Sqrt[tpD - qSQ] - Sqrt[tpD - t0DBGL])/(Sqrt[tpD - qSQ] + Sqrt[tpD - t0DBGL]); ); Null; (tpDst = (mB + mDst)^2; ); Null; (t0DstBGL = (mB - mDst)^2; ); Null; (zDst["BGL", qSQ_] = (Sqrt[tpDst - qSQ] - Sqrt[tpDst - t0DstBGL])/(Sqrt[tpDst - qSQ] + Sqrt[tpDst - t0DstBGL]); ); Null; zBGL[FF_, qSQ_] := If[StringContainsQ[FF, {"fp", "f0", "fT"}], zD["BGL", qSQ], zDst["BGL", qSQ]]; Null; (BGLlist = {"BGL_fp", "BGL_f0", "BGL_f", "BGL_F1", "BGL_F2", "BGL_g"}; ); Null; BGLff[FFmode_, NF_, qSQ_] := If[StringMatchQ[FFmode, BGLlist], (1/(PF[FFmode, zBGL[FFmode, qSQ]]*ϕF[FFmode, zBGL[FFmode, qSQ]]))*Sum[FFc[FFmode][[n + 1]]*zBGL[FFmode, qSQ]^n, {n, 0, NF}], (1/(1 - qSQ/MR[FFmode]^2))*Sum[FFc[FFmode][[n + 1]]*(zBSZ[FFmode, qSQ] - zBSZ[FFmode, 0])^n, {n, 0, NF}]]; Null; BGLff["BGL_A0", NF_, qSQ_] := (1/2)*BGLff["BGL_F2", NF, qSQ]; Null; BGLff["BGL_A1", NF_, qSQ_] := (1/(mB + mDst))*BGLff["BGL_f", NF, qSQ]; Null; BGLff["BGL_A2", NF_, qSQ_] := ((mB + mDst)/QpQmDst[qSQ])*((mB^2 - mDst^2 - qSQ)*BGLff["BGL_f", NF, qSQ] - 2*mDst*BGLff["BGL_F1", NF, qSQ]); Null; BGLff["BGL_A12", NF_, qSQ_] := BGLff["BGL_F1", NF, qSQ]/(8*mB*mDst); Null; BGLff["BGL_V", NF_, qSQ_] := ((mB + mDst)/2)*BGLff["BGL_g", NF, qSQ]; Null; BGLff["BGL_fS", NF_, qSQ_] := ((mB - mD)/(mb - mc))*BGLff["BGL_f0", NF, qSQ]; Null; BGLff["BGL_P", NF_, qSQ_] := (-((2*mDst)/(mb + mc)))*BGLff["BGL_A0", NF, qSQ]

CompoundExpression[KineticBGL = Join[Join[Join[{bT20 -> bT10, bF10 -> 0.16745498661163563*bf0}, Flatten[NSolve[BGLff["BGL_T2", 2, qDstMAX] == ((4*mB*mDst)/(mB + mDst)^2)*BGLff["BGL_T23", 2, qDstMAX], bT230]]], Flatten[NSolve[Limit[BGLff["BGL_fp", 2, qSQ], qSQ -> 0] == Limit[BGLff["BGL_f0", 2, qSQ], qSQ -> 0], bf00]]], Flatten[NSolve[Limit[BGLff["BGL_A0", 2, qSQ], qSQ -> 0] == ((8*mB*mDst)/(mB^2 - mDst^2))*Limit[BGLff["BGL_A12", 2, qSQ], qSQ -> 0], bF20]]]; ]

; 

; 

; 

(FFc["HQET_xi"] = {cxi0, cxi1, cxi2, cxi3}; ); Null; (FFc["HQET_et"] = {cet0, cet1, cet2}; ); Null; (FFc["HQET_ci2"] = {cci20, cci21, cci22}; ); Null; (FFc["HQET_ci3"] = {cci30, cci31, cci32}; ); Null; (FFc["HQET_el1"] = {cel10, cel11}; ); Null; (FFc["HQET_el2"] = {cel20, cel21}; ); Null; (FFc["HQET_el3"] = {cel30, cel31}; ); Null; (FFc["HQET_el4"] = {cel40, cel41}; ); Null; (FFc["HQET_el5"] = {cel50, cel51}; ); Null; (FFc["HQET_el6"] = {cel60, cel61}; ); Null; (cxi0 = 1; ); Null; (cci30 = 0; ); Null; (zHQET[w_] := (Sqrt[w + 1] - Sqrt[2])/(Sqrt[w + 1] + Sqrt[2]); ); Null; (HQETlist = {"HQET_xi", "HQET_et", "HQET_ci2", "HQET_ci3", "HQET_el1", "HQET_el2", "HQET_el3", "HQET_el4", "HQET_el5", "HQET_el6"}; ); Null; (Protect[z]; ); Null; (HQETff["HQET_xi", w_] = cxi0 + 8*cxi1*z + 16*(cxi1 + 2*cxi2)*z^2 + (8/3)*(9*cxi1 + 48*cxi2 + 32*cxi3)*z^3; ); Null; (HQETff["HQET_et", w_] = cet0 + 8*cet1*z + 16*(cet1 + 2*cet2)*z^2; ); Null; (HQETff["HQET_ci2", w_] = cci20 + 8*cci21*z + 16*(cci21 + 2*cci22)*z^2; ); Null; (HQETff["HQET_ci3", w_] = cci30 + 8*cci31*z + 16*(cci31 + 2*cci32)*z^2; ); Null; (HQETff["HQET_el1", w_] = cel10 + 8*cel11*z; ); Null; (HQETff["HQET_el2", w_] = cel20 + 8*cel21*z; ); Null; (HQETff["HQET_el3", w_] = cel30 + 8*cel31*z; ); Null; (HQETff["HQET_el4", w_] = cel40 + 8*cel41*z; ); Null; (HQETff["HQET_el5", w_] = cel50 + 8*cel51*z; ); Null; (HQETff["HQET_el6", w_] = cel60 + 8*cel61*z; ); Null; (wP[w_] := w + Sqrt[-1 + w^2]; ); Null; (wM[w_] := w - Sqrt[-1 + w^2]; ); Null; (rw[w_] := Log[wP[w]]/Sqrt[-1 + w^2]; ); Null; (Li2[x_] := PolyLog[2, x]; ); Null; (Ow[w_] := 1 + (w*(2*Li2[1 - zcb*wM[w]] - Li2[1 - wM[w]^2] - 2*Li2[1 - zcb*wP[w]] + Li2[1 - wP[w]^2]))/(2*Sqrt[-1 + w^2]) - w*Log[zcb]*rw[w]; ); Null; (alCV1[w_] := (-12*(w - wcb)*zcb - (-1 + zcb^2)*Log[zcb] + 4*(w - wcb)*zcb*Ow[w] + 2*(1 + w)*(-1 + (-1 + 3*w)*zcb - zcb^2)*rw[w])/(6*(w - wcb)*zcb) + ScaleFactor["CV1", w]; ); Null; (alCV2[w_] := -(zcb*(-2*(w - wcb)*(-1 + zcb) + (3 - 2*w - (-2 + 4*w)*zcb + zcb^2)*Log[zcb]) + (2 - (-1 + 5*w + 2*w^2)*zcb + (2*w + 4*w^2)*zcb^2 - (1 + w)*zcb^3)*rw[w])/(6*(w - wcb)^2*zcb^2); ); Null; (alCV3[w_] := (-2*(w - wcb)*(-1 + zcb)*zcb + (1 - (-2 + 4*w)*zcb + (3 - 2*w)*zcb^2)*Log[zcb] + (1 + w - (2*w + 4*w^2)*zcb + (-1 + 5*w + 2*w^2)*zcb^2 - 2*zcb^3)*rw[w])/(6*(w - wcb)^2*zcb); ); Null; (alCA1[w_] := (-12*(w - wcb)*zcb - (-1 + zcb^2)*Log[zcb] + 4*(w - wcb)*zcb*Ow[w] + 2*(-1 + w)*(-1 + (1 + 3*w)*zcb - zcb^2)*rw[w])/(6*(w - wcb)*zcb) + ScaleFactor["CA1", w]; ); Null; (alCA2[w_] := -(zcb*(-2*(w - wcb)*(1 + zcb) + (3 + 2*w - (2 + 4*w)*zcb + zcb^2)*Log[zcb]) + (2 + (-1 - 5*w + 2*w^2)*zcb + (-2*w + 4*w^2)*zcb^2 + (1 - w)*zcb^3)*rw[w])/(6*(w - wcb)^2*zcb^2); ); Null; (alCA3[w_] := (-2*(w - wcb)*zcb*(1 + zcb) - (1 - (2 + 4*w)*zcb + (3 + 2*w)*zcb^2)*Log[zcb] + (1 - w + (-2*w + 4*w^2)*zcb + (-1 - 5*w + 2*w^2)*zcb^2 + 2*zcb^3)*rw[w])/(6*(w - wcb)^2*zcb); ); Null; (alCS[w_] := ((-1 + zcb^2)*Log[zcb] + 2*(w - wcb)*zcb*Ow[w] - (-1 + w)*(1 + zcb)^2*rw[w])/(3*(w - wcb)*zcb) + ScaleFactor["CS", w]; ); Null; (alCP[w_] := ((-1 + zcb^2)*Log[zcb] + 2*(w - wcb)*zcb*Ow[w] - (1 + w)*(-1 + zcb)^2*rw[w])/(3*(w - wcb)*zcb) + ScaleFactor["CP", w]; ); Null; (alCT1[w_] := (-6*(w - wcb)*zcb - (-1 + zcb^2)*Log[zcb] + 2*(w - wcb)*zcb*Ow[w] + (-1 + w)*(-1 + (2 + 4*w)*zcb - zcb^2)*rw[w])/(3*(w - wcb)*zcb) + ScaleFactor["CT1", w]; ); Null; (alCT2[w_] := (2*(zcb*Log[zcb] + (1 - w*zcb)*rw[w]))/(3*(w - wcb)*zcb); ); Null; (alCT3[w_] := (2*(Log[zcb] + (w - zcb)*rw[w]))/(3*(w - wcb)); ); Null; (ScaleChoice[scale_] := Which[scale == "mbmc", Module[{}, ScaleFactor["CV1", w_] = 0; ScaleFactor["CA1", w_] = 0; ScaleFactor["CS", w_] = 0; ScaleFactor["CP", w_] = 0; ScaleFactor["CT1", w_] = 0; LigetiFit = {als -> 0.08276057040778559, epc -> 0.1807, epb -> 0.0522}; zcb = 1.3818/4.7818; wcb = (zcb^(-1) + zcb)/2; ], scale == "mbmb", Module[{}, ScaleFactor["CV1", w_] = (-2*Log[zcb]*(-1 + w*rw[w]))/3; ScaleFactor["CA1", w_] = (-2*Log[zcb]*(-1 + w*rw[w]))/3; ScaleFactor["CS", w_] = -(Log[zcb]*(1 + 2*w*rw[w]))/3; ScaleFactor["CP", w_] = -(Log[zcb]*(1 + 2*w*rw[w]))/3; ScaleFactor["CT1", w_] = -(Log[zcb]*(-3 + 2*w*rw[w]))/3; LigetiFit = {als -> 0.0716197243913529, epc -> 0.1807, epb -> 0.0522}; zcb = 1.3818/4.7818; wcb = (zcb^(-1) + zcb)/2; ]]; ); Null; (ScaleChoice["mbmb"]; ); Null; (HQETff["HQET_Δ1", w_] := -4*(w - 1)*HQETff["HQET_ci2", w] + 12*HQETff["HQET_ci3", w]; ); Null; (HQETff["HQET_Δ2", w_] := -4*HQETff["HQET_ci3", w]; ); Null; (Cgc = (-(3/2))*(Log[zcb] - 13/9); ); Null; (DalCV1[w_] = D[alCV1[w], w]; ); Null; (DalCV2[w_] = D[alCV2[w], w]; ); Null; (DalCV3[w_] = D[alCV3[w], w]; ); Null; (DalCA1[w_] = D[alCA1[w], w]; ); Null; (DalCA2[w_] = D[alCA2[w], w]; ); Null; (DalCA3[w_] = D[alCA3[w], w]; ); Null; (DalCT1[w_] = D[alCT1[w], w]; ); Null; (DalCT2[w_] = D[alCT2[w], w]; ); Null; (DalCT3[w_] = D[alCT3[w], w]; ); Null; (DalCS[w_] = D[alCS[w], w]; ); Null; (DalCP[w_] = D[alCP[w], w]; ); Null; (alCV1ap[w_] = Evaluate[Chop[Normal[Series[alCV1[w], {w, 1, 5}]]]]; ); Null; (alCV2ap[w_] = Evaluate[Chop[Normal[Series[alCV2[w], {w, 1, 5}]]]]; ); Null; (alCV3ap[w_] = Evaluate[Chop[Normal[Series[alCV3[w], {w, 1, 5}]]]]; ); Null; (alCA1ap[w_] = Evaluate[Chop[Normal[Series[alCA1[w], {w, 1, 5}]]]]; ); Null; (alCA2ap[w_] = Evaluate[Chop[Normal[Series[alCA2[w], {w, 1, 5}]]]]; ); Null; (alCA3ap[w_] = Evaluate[Chop[Normal[Series[alCA3[w], {w, 1, 5}]]]]; ); Null; (alCSap[w_] = Evaluate[Chop[Normal[Series[alCS[w], {w, 1, 5}]]]]; ); Null; (alCPap[w_] = Evaluate[Chop[Normal[Series[alCP[w], {w, 1, 5}]]]]; ); Null; (alCT1ap[w_] = Evaluate[Chop[Normal[Series[alCT1[w], {w, 1, 5}]]]]; ); Null; (alCT2ap[w_] = Evaluate[Chop[Normal[Series[alCT2[w], {w, 1, 5}]]]]; ); Null; (alCT3ap[w_] = Evaluate[Chop[Normal[Series[alCT3[w], {w, 1, 5}]]]]; ); Null; (DalCV1ap[w_] = Evaluate[Chop[Normal[Series[D[alCV1[w], w], {w, 1, 6}]]]]; ); Null; (DalCV2ap[w_] = Evaluate[Chop[Normal[Series[D[alCV2[w], w], {w, 1, 6}]]]]; ); Null; (DalCV3ap[w_] = Evaluate[Chop[Normal[Series[D[alCV3[w], w], {w, 1, 6}]]]]; ); Null; (DalCA1ap[w_] = Evaluate[Chop[Normal[Series[D[alCA1[w], w], {w, 1, 6}]]]]; ); Null; (DalCA2ap[w_] = Evaluate[Chop[Normal[Series[D[alCA2[w], w], {w, 1, 6}]]]]; ); Null; (DalCA3ap[w_] = Evaluate[Chop[Normal[Series[D[alCA3[w], w], {w, 1, 6}]]]]; ); Null; (DalCSap[w_] = Evaluate[Chop[Normal[Series[D[alCS[w], w], {w, 1, 6}]]]]; ); Null; (DalCPap[w_] = Evaluate[Chop[Normal[Series[D[alCP[w], w], {w, 1, 6}]]]]; ); Null; (DalCT1ap[w_] = Evaluate[Chop[Normal[Series[D[alCT1[w], w], {w, 1, 6}]]]]; ); Null; (DalCT2ap[w_] = Evaluate[Chop[Normal[Series[D[alCT2[w], w], {w, 1, 6}]]]]; ); Null; (DalCT3ap[w_] = Evaluate[Chop[Normal[Series[D[alCT3[w], w], {w, 1, 6}]]]]; ); Null; (HQETreplace["IW_basis"] = {funcL1[w] -> -4*(w - 1)*HQETff["HQET_ci2", w] + 12*HQETff["HQET_ci3", w], funcL2[w] -> -4*HQETff["HQET_ci3", w], funcL3[w] -> 4*HQETff["HQET_ci2", w], funcL4[w] -> 2*HQETff["HQET_et", w] - 1, funcL5[w] -> -1, funcL6[w] -> (-(2/(1 + w)))*(1 + HQETff["HQET_et", w]), funcalCV1[w] -> alCV1[w], funcalCV2[w] -> alCV2[w], funcalCV3[w] -> alCV3[w], funcalCA1[w] -> alCA1[w], funcalCA2[w] -> alCA2[w], funcalCA3[w] -> alCA3[w], funcalCP[w] -> alCP[w], funcalCS[w] -> alCS[w], funcalCT1[w] -> alCT1[w], funcalCT2[w] -> alCT2[w], funcalCT3[w] -> alCT3[w], funcEl1[w] -> HQETff["HQET_el1", w], funcEl2[w] -> HQETff["HQET_el2", w], funcEl3[w] -> HQETff["HQET_el3", w], funcEl4[w] -> HQETff["HQET_el4", w], funcEl5[w] -> HQETff["HQET_el5", w], funcEl6[w] -> HQETff["HQET_el6", w]}; ); Null; (HQETreplace["als/mc_full"] = {funcalhp[w] -> alCV1[w]*funcL1[w] + ((w + 1)/2)*(alCV2[w] + alCV3[w])*(funcL1[w] - ((w - 1)/(w + 1))*funcL4[w]) + 2*(w - 1)*(DalCV1[w] + ((w + 1)/2)*(DalCV2[w] + DalCV3[w])) - (w - 1)*alCV3[w]*funcL5[w] + Cgc*HQETff["HQET_Δ1", w], funcalhm[w] -> ((w + 1)/2)*((alCV2[w] - alCV3[w])*(funcL1[w] - ((w - 1)/(w + 1))*funcL4[w]) + 2*(w - 1)*(DalCV2[w] - DalCV3[w])) + alCV1[w]*funcL4[w] + (w + 1)*alCV3[w]*funcL5[w], funcalhV[w] -> alCV1[w]*(funcL2[w] - funcL5[w]) - alCV3[w]*(funcL4[w] - funcL5[w]) + 2*(w - 1)*DalCV1[w] + Cgc*HQETff["HQET_Δ2", w], funcalhA1[w] -> alCA1[w]*(funcL2[w] - ((w - 1)/(w + 1))*funcL5[w]) + 2*(w - 1)*DalCA1[w] + ((w - 1)/(w + 1))*alCA3[w]*(funcL4[w] - funcL5[w]) + Cgc*HQETff["HQET_Δ2", w], funcalhA2[w] -> alCA1[w]*(funcL3[w] + funcL6[w]) + alCA2[w]*(funcL2[w] + (w - 1)*funcL3[w] + funcL5[w] - (w + 1)*funcL6[w]) - (1/(w + 1))*alCA3[w]*(funcL4[w] - 3*funcL5[w]) + 2*(w - 1)*DalCA2[w] + Cgc*funcL3[w], funcalhA3[w] -> alCA1[w]*(funcL2[w] - funcL3[w] - funcL5[w] + funcL6[w]) + (w/(w + 1))*alCA3[w]*(funcL4[w] - 3*funcL5[w]) + alCA3[w]*(funcL2[w] + (w - 1)*funcL3[w] + funcL5[w] - (w + 1)*funcL6[w]) + 2*(w - 1)*(DalCA1[w] + DalCA3[w]) + Cgc*(HQETff["HQET_Δ2", w] - funcL3[w]), funcalhS[w] -> alCS[w]*(funcL1[w] - ((w - 1)/(w + 1))*funcL4[w]) + 2*(w - 1)*DalCS[w] + Cgc*HQETff["HQET_Δ1", w], funcalhP[w] -> alCP[w]*(funcL2[w] + (w - 1)*funcL3[w] + funcL5[w] - (w + 1)*funcL6[w]) + 2*(w - 1)*DalCP[w] + Cgc*(HQETff["HQET_Δ2", w] + (w - 1)*funcL3[w]), funcalhT[w] -> (alCT1[w] - alCT2[w] + alCT3[w])*funcL1[w] - alCT1[w]*funcL4[w] + 2*(w - 1)*(DalCT1[w] - DalCT2[w] + DalCT3[w]) + (alCT2[w] + alCT3[w])*funcL4[w] - 2*alCT3[w]*funcL5[w] + Cgc*HQETff["HQET_Δ1", w], funcalhT1[w] -> alCT1[w]*funcL2[w] + ((w - 1)/2)*(alCT2[w] - alCT3[w])*(funcL2[w] - funcL5[w]) + 2*(w - 1)*(DalCT1[w] + ((w - 1)/2)*(DalCT2[w] - DalCT3[w])) + (w - 1)*alCT3[w]*funcL5[w] + Cgc*HQETff["HQET_Δ2", w], funcalhT2[w] -> alCT1[w]*funcL5[w] + ((w + 1)/2)*(alCT2[w] + alCT3[w])*(funcL2[w] - funcL5[w]) + alCT3[w]*(funcL4[w] - w*funcL5[w]) + (w^2 - 1)*(DalCT2[w] + DalCT3[w]), funcalhT3[w] -> (-alCT1[w])*(funcL3[w] - funcL6[w]) + alCT2[w]*(funcL2[w] - funcL5[w]) - (1/(w + 1))*alCT3[w]*(funcL4[w] - 3*funcL5[w]) + 2*(w - 1)*DalCT2[w] - Cgc*funcL3[w]}; ); Null; (HQETreplace["als/mc_app"] = {funcalhp[w] -> alCV1ap[w]*funcL1[w] + ((w + 1)/2)*(alCV2ap[w] + alCV3ap[w])*(funcL1[w] - ((w - 1)/(w + 1))*funcL4[w]) + 2*(w - 1)*(DalCV1ap[w] + ((w + 1)/2)*(DalCV2ap[w] + DalCV3ap[w])) - (w - 1)*alCV3ap[w]*funcL5[w] + Cgc*HQETff["HQET_Δ1", w], funcalhm[w] -> ((w + 1)/2)*((alCV2ap[w] - alCV3ap[w])*(funcL1[w] - ((w - 1)/(w + 1))*funcL4[w]) + 2*(w - 1)*(DalCV2ap[w] - DalCV3ap[w])) + alCV1ap[w]*funcL4[w] + (w + 1)*alCV3ap[w]*funcL5[w], funcalhV[w] -> alCV1ap[w]*(funcL2[w] - funcL5[w]) - alCV3ap[w]*(funcL4[w] - funcL5[w]) + 2*(w - 1)*DalCV1ap[w] + Cgc*HQETff["HQET_Δ2", w], funcalhA1[w] -> alCA1ap[w]*(funcL2[w] - ((w - 1)/(w + 1))*funcL5[w]) + 2*(w - 1)*DalCA1ap[w] + ((w - 1)/(w + 1))*alCA3ap[w]*(funcL4[w] - funcL5[w]) + Cgc*HQETff["HQET_Δ2", w], funcalhA2[w] -> alCA1ap[w]*(funcL3[w] + funcL6[w]) + alCA2ap[w]*(funcL2[w] + (w - 1)*funcL3[w] + funcL5[w] - (w + 1)*funcL6[w]) - (1/(w + 1))*alCA3ap[w]*(funcL4[w] - 3*funcL5[w]) + 2*(w - 1)*DalCA2ap[w] + Cgc*funcL3[w], funcalhA3[w] -> alCA1ap[w]*(funcL2[w] - funcL3[w] - funcL5[w] + funcL6[w]) + (w/(w + 1))*alCA3ap[w]*(funcL4[w] - 3*funcL5[w]) + alCA3ap[w]*(funcL2[w] + (w - 1)*funcL3[w] + funcL5[w] - (w + 1)*funcL6[w]) + 2*(w - 1)*(DalCA1ap[w] + DalCA3ap[w]) + Cgc*(HQETff["HQET_Δ2", w] - funcL3[w]), funcalhS[w] -> alCSap[w]*(funcL1[w] - ((w - 1)/(w + 1))*funcL4[w]) + 2*(w - 1)*DalCSap[w] + Cgc*HQETff["HQET_Δ1", w], funcalhP[w] -> alCPap[w]*(funcL2[w] + (w - 1)*funcL3[w] + funcL5[w] - (w + 1)*funcL6[w]) + 2*(w - 1)*DalCPap[w] + Cgc*(HQETff["HQET_Δ2", w] + (w - 1)*funcL3[w]), funcalhT[w] -> (alCT1ap[w] - alCT2ap[w] + alCT3ap[w])*funcL1[w] - alCT1ap[w]*funcL4[w] + 2*(w - 1)*(DalCT1ap[w] - DalCT2ap[w] + DalCT3ap[w]) + (alCT2ap[w] + alCT3ap[w])*funcL4[w] - 2*alCT3ap[w]*funcL5[w] + Cgc*HQETff["HQET_Δ1", w], funcalhT1[w] -> alCT1ap[w]*funcL2[w] + ((w - 1)/2)*(alCT2ap[w] - alCT3ap[w])*(funcL2[w] - funcL5[w]) + 2*(w - 1)*(DalCT1ap[w] + ((w - 1)/2)*(DalCT2ap[w] - DalCT3ap[w])) + (w - 1)*alCT3[w]*funcL5[w] + Cgc*HQETff["HQET_Δ2", w], funcalhT2[w] -> alCT1ap[w]*funcL5[w] + ((w + 1)/2)*(alCT2ap[w] + alCT3ap[w])*(funcL2[w] - funcL5[w]) + alCT3ap[w]*(funcL4[w] - w*funcL5[w]) + (w^2 - 1)*(DalCT2ap[w] + DalCT3ap[w]), funcalhT3[w] -> (-alCT1ap[w])*(funcL3[w] - funcL6[w]) + alCT2ap[w]*(funcL2[w] - funcL5[w]) - (1/(w + 1))*alCT3ap[w]*(funcL4[w] - 3*funcL5[w]) + 2*(w - 1)*DalCT2ap[w] - Cgc*funcL3[w]}; ); Null; (HQETff["HQET_hp", NF_, w_] := HQETff["HQET_xi", w]*(1 + als*(funcalCV1[w] + ((1 + w)/2)*(funcalCV2[w] + funcalCV3[w])) + (epb + epc)*funcL1[w] + epc^2*funcEl1[w] + als*epc*funcalhp[w]); ); Null; (HQETff["HQET_hm", NF_, w_] := HQETff["HQET_xi", w]*(0 + als*((1 + w)/2)*(funcalCV2[w] - funcalCV3[w]) - (epb - epc)*funcL4[w] + epc^2*funcEl4[w] + als*epc*funcalhm[w]); ); Null; (HQETff["HQET_hV", NF_, w_] := HQETff["HQET_xi", w]*(1 + als*funcalCV1[w] + epb*(funcL1[w] - funcL4[w]) + epc*(funcL2[w] - funcL5[w]) + epc^2*(funcEl2[w] - funcEl5[w]) + als*epc*funcalhV[w]); ); Null; (HQETff["HQET_hA1", NF_, w_] := HQETff["HQET_xi", w]*(1 + als*funcalCA1[w] + epb*(funcL1[w] - ((w - 1)/(w + 1))*funcL4[w]) + epc*(funcL2[w] - ((w - 1)/(w + 1))*funcL5[w]) + epc^2*(funcEl2[w] - ((w - 1)/(w + 1))*funcEl5[w]) + als*epc*funcalhA1[w]); ); Null; (HQETff["HQET_hA2", NF_, w_] := HQETff["HQET_xi", w]*(0 + als*funcalCA2[w] + epc*(funcL3[w] + funcL6[w]) + epc^2*(funcEl3[w] + funcEl6[w]) + als*epc*funcalhA2[w]); ); Null; (HQETff["HQET_hA3", NF_, w_] := HQETff["HQET_xi", w]*(1 + als*(funcalCA1[w] + funcalCA3[w]) + epb*(funcL1[w] - funcL4[w]) + epc*(funcL2[w] - funcL3[w] - funcL5[w] + funcL6[w]) + epc^2*(funcEl2[w] - funcEl3[w] - funcEl5[w] + funcEl6[w]) + als*epc*funcalhA3[w]); ); Null; (HQETff["HQET_hS", NF_, w_] := HQETff["HQET_xi", w]*(1 + als*funcalCS[w] + (epc + epb)*(funcL1[w] - ((w - 1)/(w + 1))*funcL4[w]) + epc^2*(funcEl1[w] - ((w - 1)/(w + 1))*funcEl4[w]) + als*epc*funcalhS[w]); ); Null; (HQETff["HQET_hP", NF_, w_] := HQETff["HQET_xi", w]*(1 + als*funcalCP[w] + epc*(funcL2[w] + (w - 1)*funcL3[w] + funcL5[w] - (w + 1)*funcL6[w]) + epb*(funcL1[w] - funcL4[w]) + epc^2*(funcEl2[w] + (w - 1)*funcEl3[w] + funcEl5[w] - (w + 1)*funcEl6[w]) + als*epc*funcalhP[w]); ); Null; (HQETff["HQET_hT", NF_, w_] := HQETff["HQET_xi", w]*(1 + als*(funcalCT1[w] - funcalCT2[w] + funcalCT3[w]) + (epb + epc)*(funcL1[w] - funcL4[w]) + epc^2*(funcEl1[w] - funcEl4[w]) + als*epc*funcalhT[w]); ); Null; (HQETff["HQET_hT1", NF_, w_] := HQETff["HQET_xi", w]*(1 + als*(funcalCT1[w] + ((w - 1)/2)*(funcalCT2[w] - funcalCT3[w])) + epb*funcL1[w] + epc*funcL2[w] + epc^2*funcEl2[w] + als*epc*funcalhT1[w]); ); Null; (HQETff["HQET_hT2", NF_, w_] := HQETff["HQET_xi", w]*(0 + als*((w + 1)/2)*(funcalCT2[w] + funcalCT3[w]) - epb*funcL4[w] + epc*funcL5[w] + epc^2*funcEl5[w] + als*epc*funcalhT2[w]); ); Null; (HQETff["HQET_hT3", NF_, w_] := HQETff["HQET_xi", w]*(0 + als*funcalCT2[w] + epc*(-funcL3[w] + funcL6[w]) + epc^2*(-funcEl3[w] + funcEl6[w]) + als*epc*funcalhT3[w]); ); Null; (HQETff["HQET_fp", NF_, w_] := (1/(2*Sqrt[rD]))*((1 + rD)*HQETff["HQET_hp", NF, w] - (1 - rD)*HQETff["HQET_hm", NF, w]); ); Null; (HQETff["HQET_f0", NF_, w_] := Sqrt[rD]*(((w + 1)/(1 + rD))*HQETff["HQET_hp", NF, w] - ((w - 1)/(1 - rD))*HQETff["HQET_hm", NF, w]); ); Null; (HQETff["HQET_fS", NF_, w_] := (Sqrt[rD]/(1 + rD))*(w + 1)*HQETff["HQET_hS", NF, w]; ); Null; (HQETff["HQET_fT", NF_, w_] := ((1 + rD)/(2*Sqrt[rD]))*HQETff["HQET_hT", NF, w]; ); Null; (HQETff["HQET_A1", NF_, w_] := (Sqrt[rDst]/(1 + rDst))*(w + 1)*HQETff["HQET_hA1", NF, w]; ); Null; (HQETff["HQET_A2", NF_, w_] := ((1 + rDst)/(2*Sqrt[rDst]))*(rDst*HQETff["HQET_hA2", NF, w] + HQETff["HQET_hA3", NF, w]); ); Null; (HQETff["HQET_A0", NF_, w_] := (1/(2*Sqrt[rDst]))*((w + 1)*HQETff["HQET_hA1", NF, w] + (w*rDst - 1)*HQETff["HQET_hA2", NF, w] + (rDst - w)*HQETff["HQET_hA3", NF, w]); ); Null; HQETff["HQET_A12", NF_, w_] := (-((w + 1)/(8*Sqrt[rDst])))*((rDst - w)*HQETff["HQET_hA1", NF, w] + (w - 1)*(rDst*HQETff["HQET_hA2", NF, w] + HQETff["HQET_hA3", NF, w])); Null; (HQETff["HQET_V", NF_, w_] := ((1 + rDst)/(2*Sqrt[rDst]))*HQETff["HQET_hV", NF, w]; ); Null; (HQETff["HQET_P", NF_, w_] := (-Sqrt[rDst])*HQETff["HQET_hP", NF, w]; ); Null; (HQETff["HQET_T1", NF_, w_] := (1/(2*Sqrt[rDst]))*((1 + rDst)*HQETff["HQET_hT1", NF, w] - (1 - rDst)*HQETff["HQET_hT2", NF, w]); ); Null; (HQETff["HQET_T2", NF_, w_] := Sqrt[rDst]*(((w + 1)/(1 + rDst))*HQETff["HQET_hT1", NF, w] - ((w - 1)/(1 - rDst))*HQETff["HQET_hT2", NF, w]); ); Null; (HQETff["HQET_T23", NF_, w_] := ((1 + rDst)/(4*Sqrt[rDst]))*((w + 1)*HQETff["HQET_hT1", NF, w] + (w - 1)*HQETff["HQET_hT2", NF, w] - (w^2 - 1)*HQETff["HQET_hT3", NF, w]); ); Null; (HQETff["HQET_f", NF_, w_] := (mB + mDst)*HQETff["HQET_A1", NF, w]; ); Null; (HQETff["HQET_F1", NF_, w_] := 8*mB*mDst*HQETff["HQET_A12", NF, w]; ); Null; (HQETff["HQET_F2", NF_, w_] := 2*HQETff["HQET_A0", NF, w]; ); Null; (HQETff["HQET_g", NF_, w_] := (2/(mB + mDst))*HQETff["HQET_V", NF, w]; ); Null; (HQETallreplace[w_] = Join[Join[Join[LigetiFit, HQETreplace["IW_basis"]], HQETreplace["als/mc_full"] /. HQETreplace["IW_basis"]], {z -> zHQET[w]}]; ); Null; (DeltaKinetic["f"] = Expand[Limit[HQETff["HQET_f0", 2, w] //. HQETallreplace[w], w -> 1.503556184313367]] - Expand[Limit[HQETff["HQET_fp", 2, w] //. HQETallreplace[w], w -> 1.503556184313367]]; ); Null; (DeltaKinetic["A"] = Expand[Limit[HQETff["HQET_A0", 2, w] //. HQETallreplace[w], w -> 1.503556184313367]] - Expand[((8*mB*mDst)/(mB^2 - mDst^2))*Limit[HQETff["HQET_A12", 2, w] //. HQETallreplace[w], w -> 1.503556184313367]]; ); Null; (DeltaKinetic["T"] = Expand[Limit[HQETff["HQET_T1", 2, w] //. HQETallreplace[w], w -> 1.503556184313367]] - Expand[Limit[HQETff["HQET_T2", 2, w] //. HQETallreplace[w], w -> 1.503556184313367]]; )

; 

(SpecifyModel[FFmodel_, NF_] := Module[{FF, FFlist}, FFlist = {fp, f0, fS, fT, A0, A1, A12, V, P, T1, T2, T23}; Clear[FFreplace]; FFreplace = Table[n -> ToExpression[StringJoin[FFmodel, "ff"]][StringJoin[FFmodel, "_", ToString[n]], NF, If[FFmodel == "HQET", w, qSQ]], {n, FFlist}]; Print[StringJoin["Notice: (", FFmodel, ", NF=", ToString[NF], ") is chosen. FF is stored in 'FFreplace'."]]; ]; ); Null; (NF2off := If[FFnf == 1, {bf02 -> 0, bfp2 -> 0, bfT2 -> 0, bf2 -> 0, bF12 -> 0, bF22 -> 0, bg2 -> 0, bT12 -> 0, bT22 -> 0, bT232 -> 0, af02 -> 0, afp2 -> 0, afT2 -> 0, aA02 -> 0, aA12 -> 0, aA122 -> 0, aV2 -> 0, aT12 -> 0, aT22 -> 0, aT232 -> 0, cxi3 -> 0, cet2 -> 0, cci22 -> 0, cci32 -> 0, cel11 -> 0, cel21 -> 0, cel31 -> 0, cel41 -> 0, cel51 -> 0, cel61 -> 0}, {}]; ); Null; (SpecifyKinetic[choice_] := Switch[choice, "On", KineticChoice = {af00 -> Nothing, aA00 -> Nothing, aA120 -> Nothing, aT20 -> Nothing, aT230 -> Nothing, bf00 -> Nothing, bF20 -> Nothing, bF10 -> Nothing, bT20 -> Nothing, bT230 -> Nothing}; KineticReplace = Join[Join[Join[NF2off, KineticBSZ /. NF2off], KineticBGL /. NF2off], HQETallreplace[w]]; , "Off", KineticChoice = {}; KineticReplace = Join[NF2off, HQETallreplace[w]]; , _, "Invalid choice"]; )

CompoundExpression[Protect[fp, f0, fS, fT, A0, A1, A2, A12, V, P, T1, T2, T23]; ]

CompoundExpression[(HV = (Sqrt[QpQmD[qSQ]]/Sqrt[qSQ])*fp; )*(HS = (mB + mD)*fS; )*(HT = (-(Sqrt[QpQmD[qSQ]]/(mB + mD)))*fT; )*(HVp = (mB + mDst)*A1 - (Sqrt[QpQmDst[qSQ]]/(mB + mDst))*V; )*(HVm = (mB + mDst)*A1 + (Sqrt[QpQmDst[qSQ]]/(mB + mDst))*V; )*(HV0 = (-((8*mB*mDst)/Sqrt[qSQ]))*A12; )*(HP = (Sqrt[QpQmDst[qSQ]]/(2*mDst))*P; )*(HTp = (1/Sqrt[qSQ])*(Sqrt[QpQmDst[qSQ]]*T1 + (mB^2 - mDst^2)*T2); )*(HTm = (1/Sqrt[qSQ])*(Sqrt[QpQmDst[qSQ]]*T1 - (mB^2 - mDst^2)*T2); )*(HT0 = -4*((mB*mDst)/(mB + mDst))*T23; )*(HDVt = ((mB^2 - mD^2)/Sqrt[qSQ])*f0; )*(HDstVt = (-(Sqrt[QpQmDst[qSQ]]/Sqrt[qSQ]))*A0; )]

CompoundExpression[HHreplace = {HHVp -> HVp, HHVm -> HVm, HHV0 -> HV0, HHP -> HP, HHTp -> HTp, HHTm -> HTm, HHT0 -> HT0}; ]

(Hcal[1] = (HHVp - CVR*HHVm)^2; ); Null; (Hcal[2] = (HHVm - CVR*HHVp)^2; ); Null; (Hcal[3] = (1 - CVR)^2*HHV0^2; ); Null; (Hcal[4] = (CSR - CSL)^2*HHP^2; ); Null; (Hcal[5] = CT^2*(HHTp^2 + HHTm^2); ); Null; (Hcal[6] = CT^2*HHT0^2; ); Null; (Hcal[7] = (HHVp - CVR*HHVm)*(HHVm - CVR*HHVp) + 16*CT^2*HHTp*HHTm; ); Null; (Hcal[8] = (1 - CVR)*HHV0*(HHVp - CVR*HHVm) - 8*CT^2*HHT0*(HHTp - HHTm) + 2*CT*(CSR - CSL)*HHP*(HHTp - HHTm); ); Null; (Hcal[9] = (1 - CVR)*HHV0*(HHVm - CVR*HHVp) - 8*CT^2*HHT0*(HHTp - HHTm) - 2*CT*(CSR - CSL)*HHP*(HHTp - HHTm); ); Null; (Hcal[10] = (CSR - CSL)*CT*HHP*HHT0; ); Null; (SinReplace = {sinL -> Sqrt[1 - cosL^2], sinV -> Sqrt[1 - cosV^2]}; ); Null; (Jcal[1] = (1 - cosL)^2*sinV^2 /. SinReplace; ); Null; (Jcal[2] = (1 + cosL)^2*sinV^2 /. SinReplace; ); Null; (Jcal[3] = (Plus[4])*sinL^2*cosV^2 /. SinReplace; ); Null; (Jcal[4] = (Plus[4])*cosV^2 /. SinReplace; ); Null; (Jcal[5] = (Plus[16])*sinL^2*sinV^2 /. SinReplace; ); Null; (Jcal[6] = (Plus[64])*cosL^2*cosV^2 /. SinReplace; ); Null; (Jcal[7] = -2*sinL^2*sinV^2*Cos[2*χ] /. SinReplace; ); Null; (Jcal[8] = (Plus[4])*(1 - cosL)*sinL*cosV*sinV*Cos[χ] /. SinReplace; ); Null; (Jcal[9] = -4*(1 + cosL)*sinL*cosV*sinV*Cos[χ] /. SinReplace; ); Null; (Jcal[10] = (Plus[32])*cosL*cosV^2 /. SinReplace; ); Null; (AngularAmp[cosL_, cosV_, χ_] = Sum[Hcal[i]*Jcal[i], {i, 10}]; )

(SingleAngularAmp["w"] = Integrate[AngularAmp[cosL, cosV, χ], {cosL, -1, 1}, {cosV, -1, 1}, {χ, 0, 2*Pi}]; ); Null; (SingleAngularAmp["w_cosL"] = Integrate[AngularAmp[cosL, cosV, χ], {cosV, -1, 1}, {χ, 0, 2*Pi}]; ); Null; (SingleAngularAmp["w_cosV"] = Integrate[AngularAmp[cosL, cosV, χ], {cosL, -1, 1}, {χ, 0, 2*Pi}]; ); Null; (SingleAngularAmp["w_χ"] = Integrate[AngularAmp[cosL, cosV, χ], {cosL, -1, 1}, {cosV, -1, 1}]; )

(etaEM = 1.00662; ); Null; (DecayRate["D", "w"] := 2*mB*mD*((etaEM^2*GF^2)/(192*Pi^3*mB^3))*qSQ*Sqrt[QpQmD[qSQ]]*((1 + CVR)^2*HV^2 + (3/2)*(CSL + CSR)^2*HS^2 + 8*CT^2*HT^2) /. FFreplace /. qSQ -> qSQfuncD[w]; ); Null; (DecayRate["Dst", var_] := 2*BrDstDpi*mB*mDst*((3*etaEM^2*GF^2)/(4096*Pi^4*mB^3))*qSQ*Sqrt[QpQmDst[qSQ]]*SingleAngularAmp[var] /. HHreplace /. FFreplace /. qSQ -> qSQfuncDst[w]; )

(DecayRateTauSM["D", "w"] := 2*mB*mD*((etaEM^2*GF^2)/(192*Pi^3*mB^3))*(1 - mtau^2/qSQ)^2*qSQ*Sqrt[QpQmD[qSQ]]*((1 + mtau^2/(2*qSQ))*HV^2 + (3/2)*(mtau^2/qSQ)*HDVt^2) /. FFreplace /. qSQ -> qSQfuncD[w]; ); Null; (DecayRateTauSM["Dst", "w"] := 2*mB*mDst*((etaEM^2*GF^2)/(192*Pi^3*mB^3))*(1 - mtau^2/qSQ)^2*qSQ*Sqrt[QpQmDst[qSQ]]*((1 + mtau^2/(2*qSQ))*(HVp^2 + HVm^2 + HV0^2) + (3/2)*(mtau^2/qSQ)*HDstVt^2) /. HHreplace /. FFreplace /. qSQ -> qSQfuncDst[w]; )

(ObsFunc["Belle15", "w"] := ObsScale["Belle15"]*DecayRate["D", "w"]; ); Null; (ObsFunc["BelleII25", "w"] := ObsScale["BelleII25"]*DecayRate["D", "w"]; ); Null; (ObsFunc["Belle17", "w"] := ObsScale["Belle17"]*DecayRate["Dst", "w"]; ); Null; (ObsFunc["Belle17", "cosL"] := ObsScale["Belle17"]*DecayRate["Dst", "cosL"]; ); Null; (ObsFunc["Belle17", "cosV"] := ObsScale["Belle17"]*DecayRate["Dst", "cosV"]; ); Null; (ObsFunc["Belle17", "χ"] := ObsScale["Belle17"]*DecayRate["Dst", "χ"]; ); Null; (ObsFunc["Belle23", "w"] := ObsScale["Belle23"]*DecayRate["Dst", "w"]; ); Null; (ObsFunc["Belle23", "cosL"] := ObsScale["Belle23"]*DecayRate["Dst", "cosL"]; ); Null; (ObsFunc["Belle23", "cosV"] := ObsScale["Belle23"]*DecayRate["Dst", "cosV"]; ); Null; (ObsFunc["Belle23", "χ"] := ObsScale["Belle23"]*DecayRate["Dst", "χ"]; ); Null; (ObsFunc["BelleII23", "w"] := ObsScale["BelleII23"]*DecayRate["Dst", "w"]; ); Null; (ObsFunc["BelleII23", "cosL"] := ObsScale["BelleII23"]*DecayRate["Dst", "cosL"]; ); Null; (ObsFunc["BelleII23", "cosV"] := ObsScale["BelleII23"]*DecayRate["Dst", "cosV"]; ); Null; (ObsFunc["BelleII23", "χ"] := ObsScale["BelleII23"]*DecayRate["Dst", "χ"]; ); Null; (ObsFunc["BrRatio2025", "D"] := LifeTime["Bp"]*DecayRate["D", "w"]; ); Null; (ObsFunc["BrRatio2025", "Dst"] := LifeTime["Bp"]*(DecayRate["Dst", "w"]/BrDstDpi); ); Null; (ObsFunc["BrRatioTau", "D"] := LifeTime["Bp"]*DecayRateTauSM["D", "w"]; ); Null; (ObsFunc["BrRatioTau", "Dst"] := LifeTime["Bp"]*DecayRateTauSM["Dst", "w"]; )

CompoundExpression[If[FileExistsQ["Stan_inputs/IntFunctions_HQET.ma"], Print["Files already exist. Skip this step."], FFmodel = "HQET"; SpecifyModel[FFmodel, "irelevant"]; ExpandDecayRate["D", "w"] = Expand[DecayRate["D", "w"]]; ExpandDecayRate["Dst", "w"] = Expand[DecayRate["Dst", "w"]]; ExpandDecayRate["Dst", "w_cosL"] = Expand[DecayRate["Dst", "w_cosL"]]; ExpandDecayRate["Dst", "w_cosV"] = Expand[DecayRate["Dst", "w_cosV"]]; ExpandDecayRate["Dst", "w_χ"] = Expand[DecayRate["Dst", "w_χ"]]; ExpandObsFunc["Belle15", "w"] := ObsScale["Belle15"]*ExpandDecayRate["D", "w"]; ExpandObsFunc["BelleII25", "w"] := ObsScale["BelleII25"]*ExpandDecayRate["D", "w"]; ExpandObsFunc["Belle17", "w"] := ObsScale["Belle17"]*ExpandDecayRate["Dst", "w"]; ExpandObsFunc["Belle23", "w"] := ObsScale["Belle23"]*ExpandDecayRate["Dst", "w"]; ExpandObsFunc["BelleII23", "w"] := ObsScale["BelleII23"]*ExpandDecayRate["Dst", "w"]; ExpandObsFunc["BrRatio2025", "D"] := LifeTime["Bp"]*ExpandDecayRate["D", "w"]; ExpandObsFunc["BrRatio2025", "Dst"] := LifeTime["Bp"]*(ExpandDecayRate["Dst", "w"]/BrDstDpi); ExpandDecayRateTauSM["D", "w"] = Expand[DecayRateTauSM["D", "w"]]; ExpandDecayRateTauSM["Dst", "w"] = Expand[DecayRateTauSM["Dst", "w"]]; ExpandObsFunc["BrRatioTau", "D"] := LifeTime["Bp"]*ExpandDecayRateTauSM["D", "w"]; ExpandObsFunc["BrRatioTau", "Dst"] := LifeTime["Bp"]*ExpandDecayRateTauSM["Dst", "w"]; Clear[FFmodel]; ]; ]

; 

CompoundExpression[ParallelEvaluate[Off[NIntegrate::ncvb]; Off[NIntegrate::slwcon]; ]]

CompoundExpression[IntPerform[f_, variable_, x_, {x1_, x2_}] := NIntegrate[f/Times @@ (variable^Exponent[f, variable]), {x, x1, x2}, PrecisionGoal -> 8, AccuracyGoal -> 8, MaxRecursion -> 6, Method -> {"LocalAdaptive", "SymbolicProcessing" -> 0}]*Times @@ (variable^Exponent[f, variable])]

(wDstMax = w /. NSolve[qSQfuncDst[w] == 0, w][[1]]; ); Null; PreIntPerform := Module[{IntVariable, IntCoefficient}, Do[Clear[IntVariable, IntCoefficient]; IntVariable = Variables[DecayRate["Dst", StringJoin["w_", dist]]] /. {Cos[2*χ] -> Cos2Chi, Cos[χ] -> CosChi, w -> Nothing}; IntCoefficient = MonomialList[DecayRate["Dst", StringJoin["w_", dist]] /. {Cos[2*χ] -> Cos2Chi, Cos[χ] -> CosChi}, IntVariable]; DecayRate["Dst", dist] = Total[ParallelMap[IntPerform[#1, IntVariable, w, {1.0001, wDstMax - 0.0001}] & , IntCoefficient]] /. {Cos2Chi -> Cos[2*χ], CosChi -> Cos[χ]}; , {dist, {"cosL", "cosV", "χ"}}]; ]

CompoundExpression[ObsIntPerform[ObsName_, Dist_] := Module[{NameBin, IntVariable, IntBin, IntCoefficient}, If[StringContainsQ[ObsName, "BrRatio"], IntVariable = Variables[ObsFunc[ObsName, Dist]] /. w -> Nothing; IntCoefficient = MonomialList[ObsFunc[ObsName, Dist], IntVariable]; IntBin = {1, If[ObsName == "BrRatio2025", wMax[Dist] - 0.0001, wMaxTau[Dist] - 0.0001]}; Obs[ObsName, Dist] = Total[ParallelMap[IntPerform[#1, IntVariable, w, IntBin] & , IntCoefficient]]; , NameBin = Bin[ObsName, Dist]; IntVariable = Variables[ObsFunc[ObsName, Dist]] /. ToExpression[Dist] -> Nothing; IntCoefficient = MonomialList[ObsFunc[ObsName, Dist], IntVariable]; Do[IntBin = NameBin[[n ;; n + 1]]; Obs[ObsName, Dist, n] = Total[ParallelMap[IntPerform[#1, IntVariable, ToExpression[Dist], IntBin] & , IntCoefficient]]; , {n, Length[NameBin] - 1}]; ]]]

(IWcutNLO = {cci20 -> 0, cci21 -> 0, cci22 -> 0, cci31 -> 0, cci32 -> 0, cet0 -> 0, cet1 -> 0, cet2 -> 0}; ); Null; (IWcutNNLO = {cel10 -> 0, cel11 -> 0, cel20 -> 0, cel21 -> 0, cel30 -> 0, cel31 -> 0, cel40 -> 0, cel41 -> 0, cel50 -> 0, cel51 -> 0, cel60 -> 0, cel61 -> 0}; ); Null; (OrderReplaceAll[order_] := Which[order == "LO", {Nals -> 0, Nepc -> 0, Nepb -> 0}, order == "als-1", {Nals -> 1, Nepc -> 0, Nepb -> 0}, order == "als-2", {Nals -> 2, Nepc -> 0, Nepb -> 0}, order == "epc-1", {Nals -> 0, Nepc -> 1, Nepb -> 0}, order == "epc-2", {Nals -> 0, Nepc -> 2, Nepb -> 0}, order == "epc-3", {Nals -> 0, Nepc -> 3, Nepb -> 0}, order == "epc-4", {Nals -> 0, Nepc -> 4, Nepb -> 0}, order == "epb-1", {Nals -> 0, Nepc -> 0, Nepb -> 1}, order == "epb-2", {Nals -> 0, Nepc -> 0, Nepb -> 2}, order == "epb-epc-1", {Nals -> 0, Nepc -> 1, Nepb -> 1}, order == "als-epc-1", {Nals -> 1, Nepc -> 1, Nepb -> 0}, order == "als-epc-2", {Nals -> 2, Nepc -> 1, Nepb -> 0}, order == "als-epc-3", {Nals -> 1, Nepc -> 2, Nepb -> 0}, order == "als-epc-4", {Nals -> 1, Nepc -> 3, Nepb -> 0}, order == "als-epb-1", {Nals -> 1, Nepc -> 0, Nepb -> 1}]; ); Null; (OrderListAll = {"LO", "als-1", "als-2", "epc-1", "epc-2", "epc-3", "epc-4", "epb-1", "epb-2", "epb-epc-1", "als-epc-1", "als-epc-2", "als-epc-3", "als-epc-4", "als-epb-1"}; ); Null; (NzOrder[order_] := If[StringContainsQ[order, "als"], 4, 6]; ); Null; (ExpandObs[obs_, order_] := Coefficient[Coefficient[Coefficient[obs, als, Nals /. OrderReplaceAll[order]], epc, Nepc /. OrderReplaceAll[order]], epb, Nepb /. OrderReplaceAll[order]]*als^Nals*epc^Nepc*epb^Nepb /. OrderReplaceAll[order]; ); Null; (HQETset[model_] := Which[model == "3/2/1", FFnf = 2; HQETname = "321"; , model == "2/1/0", FFnf = 1; HQETname = "210"; , model =!= "3/2/1" || "2/1/0", Print["HQET model is not valid"]; ]; ); Null; (time := SessionTime[]; ); Null; PreIntPerformHQET := Module[{OrderDecayRateFunc, IntVariable, IntCoefficient}, Do[Do[Clear[IntVariable, IntCoefficient, TimeStart]; TimeStart = time; Print[StringJoin["Expanding (w_", dist, ",", order, ") term"]]; OrderDecayRateFunc["Dst", StringJoin["w_", dist], order] = Sum[Coefficient[ExpandObs[ExpandDecayRate["Dst", StringJoin["w_", dist]], order] /. HQETreplace["als/mc_app"] /. HQETreplace["IW_basis"], z, i]*z^i, {i, 0, NzOrder[order]}] /. {z -> zHQET[w]}; IntVariable = Variables[OrderDecayRateFunc["Dst", StringJoin["w_", dist], order]] /. {Cos[2*χ] -> Cos2Chi, Cos[χ] -> CosChi, w -> Nothing}; IntCoefficient = MonomialList[OrderDecayRateFunc["Dst", StringJoin["w_", dist], order] /. {Cos[2*χ] -> Cos2Chi, Cos[χ] -> CosChi}, IntVariable]; Print[StringJoin["Integrating (w_", dist, ",", order, ") term"]]; OrderDecayRate["Dst", dist, order] = Total[ParallelMap[IntPerform[#1, IntVariable, w, {1, wDstMax - 0.001}] & , IntCoefficient]] /. {Cos2Chi -> Cos[2*χ], CosChi -> Cos[χ]}; Print[StringJoin["Calculation time: ", ToString[time - TimeStart], "sec."]]; , {order, OrderListAll}]; DecayRate["Dst", dist] = Sum[OrderDecayRate["Dst", dist, OrderName], {OrderName, OrderListAll}]; , {dist, {"cosL", "cosV", "χ"}}]; ]; Null; (ObsIntPerformHQET[ObsName_, Dist_] := Module[{OrderObsFunc, NameBin, IntVariable, IntBin, IntCoefficient}, Print[StringJoin["[Integration for ", ObsName, "]"]]; TimeStart = time; If[Dist == "cosL" || Dist == "cosV" || Dist == "χ", Print["Please use 'ObsIntPerform' for (cosL, cosV, χ). Use 'ObsIntPerformHQET' for w"]; , If[StringContainsQ[ObsName, "BrRatio"], Do[OrderObsFunc[ObsName, Dist, order] = Sum[Coefficient[ExpandObs[ExpandObsFunc[ObsName, Dist], order] /. HQETreplace["als/mc_app"] /. HQETreplace["IW_basis"], z, i]*z^i, {i, 0, NzOrder[order]}] /. {z -> zHQET[w]}; IntVariable = Variables[OrderObsFunc[ObsName, Dist, order]] /. w -> Nothing; IntCoefficient = MonomialList[OrderObsFunc[ObsName, Dist, order], IntVariable]; IntBin = {1, If[ObsName == "BrRatio2025", wMax[Dist] - 0.001, wMaxTau[Dist] - 0.001]}; OrderObs[ObsName, Dist, order] = Total[ParallelMap[IntPerform[#1, IntVariable, w, IntBin] & , IntCoefficient]]; , {order, OrderListAll}]; Obs[ObsName, Dist] = Sum[OrderObs[ObsName, Dist, order], {order, OrderListAll}]; , Do[Print[StringJoin["Expanding ", order, " term"]]; OrderObsFunc[ObsName, Dist, order] = Sum[Coefficient[ExpandObs[ExpandObsFunc[ObsName, Dist], order] /. HQETreplace["als/mc_app"] /. HQETreplace["IW_basis"], z, i]*z^i, {i, 0, NzOrder[order]}] /. {z -> zHQET[w]}; NameBin = Bin[ObsName, Dist]; IntVariable = Variables[OrderObsFunc[ObsName, Dist, order]] /. ToExpression[Dist] -> Nothing; IntCoefficient = MonomialList[OrderObsFunc[ObsName, Dist, order], IntVariable]; Print[StringJoin["Integrating ", order, " term"]]; Do[IntBin = NameBin[[n ;; n + 1]]; OrderObs[ObsName, Dist, order, n] = Total[ParallelMap[IntPerform[#1, IntVariable, ToExpression[Dist], IntBin] & , IntCoefficient]]; , {n, Length[NameBin] - 1}]; , {order, OrderListAll}]; Do[Obs[ObsName, Dist, n] = Sum[OrderObs[ObsName, Dist, order, n], {order, OrderListAll}]; , {n, Length[NameBin] - 1}]]]; Print[StringJoin["Calculation time: ", ToString[time - TimeStart], "sec."]]; ]; )
);

StanFitFunctionUBSetup[] := (
CompoundExpression[Clear[FFmodel]; ]

(BcN[state_] := Which[state == "0+", 2, state == "1-", 4, state == "0-", 3, state == "1+", 4]; ); Null; (BcMn["0+", 1] = 6.714; ); Null; (BcMn["0+", 2] = 7.107; ); Null; (BcMn["1-", 1] = 6.326; ); Null; (BcMn["1-", 2] = 6.89; ); Null; (BcMn["1-", 3] = 7.02; ); Null; (BcMn["1-", 4] = 7.252; ); Null; (BcMn["0-", 1] = 6.271; ); Null; (BcMn["0-", 2] = 6.871; ); Null; (BcMn["0-", 3] = 7.239; ); Null; (BcMn["1+", 1] = 6.757; ); Null; (BcMn["1+", 2] = 6.776; ); Null; (BcMn["1+", 3] = 7.134; ); Null; (BcMn["1+", 4] = 7.142; ); Null; (nf = 2.5; ); Null; (mBINPUT = 5.2793; ); Null; (mDINPUT = 1.8672; ); Null; (mBstarINPUT = 5.3247; ); Null; (mDstarINPUT = 2.0085; ); Null; (rDBINPUT = mDINPUT/mBINPUT; ); Null; (rDstarBINPUT = mDstarINPUT/mBINPUT; ); Null; (rDBstarINPUT = mDINPUT/mBstarINPUT; ); Null; (rDstarBstarINPUT = mDstarINPUT/mBstarINPUT; ); Null; (Betaj = (m1 + m2)/(2*Sqrt[m1*m2]); ); Null; (poleZj = (Sqrt[(m1 + m2)^2 - m3^2] - 2*Sqrt[m1*m2])/(Sqrt[(m1 + m2)^2 - m3^2] + 2*Sqrt[m1*m2]); ); Null; (ReplaceBD[state_, n_] := {m1 -> mBINPUT, m2 -> mDINPUT, m3 -> BcMn[state, n]}; ); Null; (ReplaceBstarD[state_, n_] := {m1 -> mBstarINPUT, m2 -> mDINPUT, m3 -> BcMn[state, n]}; ); Null; (ReplaceBDstar[state_, n_] := {m1 -> mBINPUT, m2 -> mDstarINPUT, m3 -> BcMn[state, n]}; ); Null; (ReplaceBstarDstar[state_, n_] := {m1 -> mBstarINPUT, m2 -> mDstarINPUT, m3 -> BcMn[state, n]}; ); Null; (βS[j_] := Which[j == 1, Betaj /. ReplaceBD["0+", 1], j == 2 || j == 3, Betaj /. ReplaceBstarDstar["0+", 1]]; ); Null; (βP[j_] := Which[j == 1, Betaj /. ReplaceBDstar["0-", 1], j == 2, Betaj /. ReplaceBstarD["0-", 1], j == 3, Betaj /. ReplaceBstarDstar["0-", 1]]; ); Null; (βV[j_] := Which[j == 1, Betaj /. ReplaceBD["1-", 1], j == 2 || j == 3 || j == 6 || j == 7, Betaj /. ReplaceBstarDstar["1-", 1], j == 4, Betaj /. ReplaceBDstar["1-", 1], j == 5, Betaj /. ReplaceBstarD["1-", 1]]; ); Null; (βA[j_] := Which[j == 1 || j == 5, Betaj /. ReplaceBDstar["1+", 1], j == 2 || j == 6, Betaj /. ReplaceBstarD["1+", 1], j == 3 || j == 4 || j == 7, Betaj /. ReplaceBstarDstar["1+", 1]]; ); Null; (Zjn["0+", j_, n_] := Which[j == 1, poleZj /. ReplaceBD["0+", n], j == 2 || j == 3, poleZj /. ReplaceBstarDstar["0+", n]]; ); Null; (Zjn["0-", j_, n_] := Which[j == 1, poleZj /. ReplaceBDstar["0-", n], j == 2, poleZj /. ReplaceBstarD["0-", n], j == 3, poleZj /. ReplaceBstarDstar["0-", n]]; ); Null; (Zjn["1-", j_, n_] := Which[j == 1, poleZj /. ReplaceBD["1-", n], j == 2 || j == 3 || j == 6 || j == 7, poleZj /. ReplaceBstarDstar["1-", n], j == 4, poleZj /. ReplaceBDstar["1-", n], j == 5, poleZj /. ReplaceBstarD["1-", n]]; ); Null; (Zjn["1+", j_, n_] := Which[j == 1 || j == 5, poleZj /. ReplaceBDstar["1+", n], j == 2 || j == 6, poleZj /. ReplaceBstarD["1+", n], j == 3 || j == 4 || j == 7, poleZj /. ReplaceBstarDstar["1+", n]]; ); Null; (wTOz[z_] := 2*((1 + z)/(1 - z))^2 - 1; ); Null; (HQETff["HQET_hVbar", NF_, w_] := HQETff["HQET_xi", w]*(1 + als*funcalCV1[w] + epc*(funcL1[w] - funcL4[w]) + epb*(funcL2[w] - funcL5[w]) + epc^2*(funcEl1[w] - funcEl4[w])); ); Null; (HQETff["HQET_hA1bar", NF_, w_] := HQETff["HQET_xi", w]*(1 + als*funcalCA1[w] + epc*(funcL1[w] - funcL4[w]*((w - 1)/(w + 1))) + epb*(funcL2[w] - funcL5[w]*((w - 1)/(w + 1))) + epc^2*(funcEl1[w] - funcEl4[w]*((w - 1)/(w + 1)))); ); Null; (HQETff["HQET_hA2bar", NF_, w_] := HQETff["HQET_xi", w]*((-als)*funcalCA3[w] + epb*(funcL3[w] + funcL6[w])); ); Null; (HQETff["HQET_hA3bar", NF_, w_] := HQETff["HQET_xi", w]*(1 + als*(funcalCA1[w] - funcalCA2[w]) + epc*(funcL1[w] - funcL4[w]) + epb*(funcL2[w] - funcL3[w] - funcL5[w] + funcL6[w]) + epc^2*(funcEl1[w] - funcEl4[w])); ); Null; (HQETff["HQET_h1", NF_, w_] := HQETff["HQET_xi", w]*(1 + als*(funcalCV1[w] + ((1 + w)*(funcalCV2[w] + funcalCV3[w]))/2) + (epb + epc)*funcL2[w] + epc^2*funcEl2[w]); ); Null; (HQETff["HQET_h2", NF_, w_] := HQETff["HQET_xi", w]*((als*(1 + w)*(funcalCV2[w] - funcalCV3[w]))/2 + (epc - epb)*funcL5[w] + epc^2*funcEl5[w]); ); Null; (HQETff["HQET_h3", NF_, w_] := HQETff["HQET_xi", w]*(1 + als*funcalCV1[w] + epc*(funcL2[w] + (w - 1)*funcL3[w] + funcL5[w] - (w + 1)*funcL6[w]) + epb*(funcL2[w] - funcL5[w]) + epc^2*(funcEl2[w] + (w - 1)*funcEl3[w] + funcEl5[w] - (w + 1)*funcEl6[w])); ); Null; (HQETff["HQET_h4", NF_, w_] := HQETff["HQET_xi", w]*(1 + als*funcalCV1[w] + epc*(funcL2[w] - funcL5[w]) + epb*(funcL2[w] + (w - 1)*funcL3[w] + funcL5[w] - (w + 1)*funcL6[w]) + epc^2*(funcEl2[w] - funcEl5[w])); ); Null; (HQETff["HQET_h5", NF_, w_] := HQETff["HQET_xi", w]*((-als)*funcalCV2[w] + epc*(funcL3[w] - funcL6[w]) + epc^2*(funcEl3[w] - funcEl6[w])); ); Null; (HQETff["HQET_h6", NF_, w_] := HQETff["HQET_xi", w]*((-als)*funcalCV3[w] + epb*(funcL3[w] - funcL6[w])); ); Null; (HQETff["HQET_h7", NF_, w_] := HQETff["HQET_xi", w]*(1 + als*(funcalCA1[w] + (w - 1)*((funcalCA2[w] - funcalCA3[w])/2)) + (epb + epc)*funcL2[w] + epc^2*funcEl2[w]); ); Null; (HQETff["HQET_h8", NF_, w_] := HQETff["HQET_xi", w]*(als*(w + 1)*((funcalCA2[w] + funcalCA3[w])/2) + (epc - epb)*funcL5[w] + epc^2*funcEl5[w]); ); Null; (HQETff["HQET_h9", NF_, w_] := HQETff["HQET_xi", w]*((-als)*funcalCA2[w] + epc*(funcL3[w] - funcL6[w]) + epc^2*(funcEl3[w] - funcEl6[w])); ); Null; (HQETff["HQET_h10", NF_, w_] := HQETff["HQET_xi", w]*(als*funcalCA3[w] + epb*(funcL3[w] - funcL6[w])); ); Null; (Sj[1, w_] := HQETff["HQET_hp", NF, w] - ((1 + rDBINPUT)/(1 - rDBINPUT))*((w - 1)/(w + 1))*HQETff["HQET_hm", NF, w]; ); Null; (Sj[2, w_] := HQETff["HQET_h1", NF, w] - ((1 + rDstarBstarINPUT)/(1 - rDstarBstarINPUT))*((w - 1)/(w + 1))*HQETff["HQET_h2", NF, w]; ); Null; (Sj[3, w_] := w*(HQETff["HQET_h1", NF, w] - ((1 + rDstarBstarINPUT)/(1 - rDstarBstarINPUT))*((w - 1)/(w + 1))*HQETff["HQET_h2", NF, w]) + ((w - 1)/(1 - rDstarBstarINPUT))*(rDstarBstarINPUT*HQETff["HQET_h3", NF, w] - HQETff["HQET_h4", NF, w] + (1 - rDstarBstarINPUT*w)*HQETff["HQET_h5", NF, w] + (w - rDstarBstarINPUT)*HQETff["HQET_h6", NF, w]); ); Null; (Sj[4, w_] := 0; ); Null; (Sj[5, w_] := 0; ); Null; (Sj[6, w_] := 0; ); Null; (Sj[7, w_] := 0; ); Null; (Sjz[n_, z_] := Sj[n, w] /. {als -> 0} //. HQETreplace["IW_basis"] /. w -> wTOz[z]; ); Null; (Pj[1, w_] := (1/(1 + rDstarBINPUT))*((w + 1)*HQETff["HQET_hA1", NF, w] - (1 - rDstarBINPUT*w)*HQETff["HQET_hA2", NF, w] - (w - rDstarBINPUT)*HQETff["HQET_hA3", NF, w]); ); Null; (Pj[2, w_] := (1/(1 + rDBstarINPUT))*(rDBstarINPUT*(w + 1)*HQETff["HQET_hA1bar", NF, w] - (rDBstarINPUT - w)*HQETff["HQET_hA2bar", NF, w] - (rDBstarINPUT*w - 1)*HQETff["HQET_hA3bar", NF, w]); ); Null; (Pj[3, w_] := HQETff["HQET_h7", NF, w] - ((1 - rDstarBstarINPUT)/(1 + rDstarBstarINPUT))*HQETff["HQET_h8", NF, w]; ); Null; (Pj[4, w_] := 0; ); Null; (Pj[5, w_] := 0; ); Null; (Pj[6, w_] := 0; ); Null; (Pj[7, w_] := 0; ); Null; (Pjz[n_, z_] := Pj[n, w] /. {als -> 0} //. HQETreplace["IW_basis"] /. w -> wTOz[z]; ); Null; (Vj[1, w_] := HQETff["HQET_hp", NF, w] - ((1 - rDBINPUT)/(1 + rDBINPUT))*HQETff["HQET_hm", NF, w]; ); Null; (Vj[2, w_] := HQETff["HQET_h1", NF, w] - ((1 - rDstarBstarINPUT)/(1 + rDstarBstarINPUT))*HQETff["HQET_h2", NF, w]; ); Null; (Vj[3, w_] := w*(HQETff["HQET_h1", NF, w] - ((1 - rDstarBstarINPUT)/(1 + rDstarBstarINPUT))*HQETff["HQET_h2", NF, w]) + (1/(1 + rDstarBstarINPUT))*((1 - rDstarBstarINPUT*w)*HQETff["HQET_h3", NF, w] + (rDstarBstarINPUT - w)*HQETff["HQET_h4", NF, w] + (w^2 - 1)*(rDstarBstarINPUT*HQETff["HQET_h5", NF, w] + HQETff["HQET_h6", NF, w])); ); Null; (Vj[4, w_] := HQETff["HQET_hV", NF, w]; ); Null; (Vj[5, w_] := HQETff["HQET_hVbar", NF, w]; ); Null; (Vj[6, w_] := HQETff["HQET_h3", NF, w]; ); Null; (Vj[7, w_] := HQETff["HQET_h4", NF, w]; ); Null; (Vjz[n_, z_] := Vj[n, w] /. {als -> 0} //. HQETreplace["IW_basis"] /. w -> wTOz[z]; ); Null; (Aj[1, w_] := HQETff["HQET_hA1", NF, w]; ); Null; (Aj[2, w_] := HQETff["HQET_hA1bar", NF, w]; ); Null; (Aj[3, w_] := HQETff["HQET_h7", NF, w] - ((w - 1)/(w + 1))*HQETff["HQET_h8", NF, w] + (w - 1)*HQETff["HQET_h10", NF, w]; ); Null; (Aj[4, w_] := HQETff["HQET_h7", NF, w] + ((w - 1)/(w + 1))*HQETff["HQET_h8", NF, w] + (w - 1)*HQETff["HQET_h9", NF, w]; ); Null; (Aj[5, w_] := (1/(1 - rDstarBINPUT))*((w - rDstarBINPUT)*HQETff["HQET_hA1", NF, w] - (w - 1)*(rDstarBINPUT*HQETff["HQET_hA2", NF, w] + HQETff["HQET_hA3", NF, w])); ); Null; (Aj[6, w_] := (1/(1 - rDBstarINPUT))*((1 - rDBstarINPUT*w)*HQETff["HQET_hA1bar", NF, w] + (w - 1)*(HQETff["HQET_hA2bar", NF, w] + rDBstarINPUT*HQETff["HQET_hA3bar", NF, w])); ); Null; (Aj[7, w_] := HQETff["HQET_h7", NF, w] - ((1 + rDstarBstarINPUT)/(1 - rDstarBstarINPUT))*((w - 1)/(w + 1))*HQETff["HQET_h8", NF, w]; ); Null; (Ajz[n_, z_] := Aj[n, w] /. {als -> 0} //. HQETreplace["IW_basis"] /. w -> wTOz[z]; ); Null; (PreFactor["outer", OverAll_, β_, n_, m_] := (Sqrt[2*Pi*OverAll]*(4^(2 - n)/(β + 1)^m)*(1 + z)^(n + 2)*(1 - z)^(m - 9/2))/(1 - ((β - 1)/(β + 1))*z)^m; ); Null; (PreFactor["inner", state_, j_] := Product[If[Element[Zjn[state, j, n], Reals], (z - Zjn[state, j, n])/(1 - z*Zjn[state, j, n]), 1], {n, BcN[state]}]; ); Null; (ChiFuncZplane["0+", j_, z_] := Sjz[j, z]*PreFactor["inner", "0+", j]*PreFactor["outer", (nf/(64*Pi^2))*(1 + If[j == 2, 1, 0])*(βS[j]^2 - 1), βS[j], 1, 4]; ); Null; (ChiFuncZplane["1-", j_, z_] := Vjz[j, z]*PreFactor["inner", "1-", j]*Which[j == 1 || j == 3, PreFactor["outer", (nf/(3072*Pi^2))*βV[j]^2, βV[j], 0, 5], j == 2, PreFactor["outer", (nf/(3072*Pi^2))*2*βV[j]^2, βV[j], 0, 5], j == 4 || j == 5 || j == 6 || j == 7, PreFactor["outer", (nf/(3072*Pi^2))*2, βV[j], 0, 4]]; ); Null; (ChiFuncZplane["1-CLN", j_, z_] := Vjz[j, z]*PreFactor["inner", "1-", j]*Which[j == 1 || j == 3, PreFactor["outer", (nf/(3072*Pi^2))*βV[j]^2, βV[j], 0, 5], j == 2, PreFactor["outer", (nf/(3072*Pi^2))*2*βV[j]^2, βV[j], 0, 5], j == 4 || j == 5 || j == 6 || j == 7, PreFactor["outer", (nf/(3072*Pi^2))*2, βV[j], 0, 5]]; ); Null; (ChiFuncZplane["0-", j_, z_] := Pjz[j, z]*PreFactor["inner", "0-", j]*PreFactor["outer", (nf/(256*Pi^2))*(1 + If[j == 3, 1, 0])*βP[j]^2, βP[j], 0, 4]; ); Null; (ChiFuncZplane["1+", j_, z_] := Ajz[j, z]*PreFactor["inner", "1+", j]*Which[j == 1 || j == 2 || j == 3 || j == 4, PreFactor["outer", (nf/(768*Pi^2))*2, βA[j], 1, 4], j == 5 || j == 6, PreFactor["outer", (nf/(768*Pi^2))*(βA[j]^2 - 1), βA[j], 1, 5], j == 7, PreFactor["outer", (nf/(768*Pi^2))*2*(βA[j]^2 - 1), βA[j], 1, 5]]; )

(FFmodel = "HQET"; ); Null; (HQETmodel = "3/2/1"; ); Null; (HQETset[HQETmodel]; ); Null; (SubIWcut = {cci20 -> 0, cci21 -> 0, cci22 -> 0, cci31 -> 0, cci32 -> 0, cet0 -> 0, cet1 -> 0, cet2 -> 0}; ); Null; (SubSubIWcut = {cel10 -> 0, cel11 -> 0, cel20 -> 0, cel21 -> 0, cel30 -> 0, cel31 -> 0, cel40 -> 0, cel41 -> 0, cel50 -> 0, cel51 -> 0, cel60 -> 0, cel61 -> 0}; ); Null; (AllIWcut = Join[Join[{cxi1 -> 0, cxi2 -> 0, cxi3 -> 0}, SubIWcut], SubSubIWcut]; ); Null; (Do[ChiFuncDerivativeFunc[state, n, j, z_] = (1/n!)*D[ChiFuncZplane[state, j, z] /. HQETreplace["IW_basis"] /. {als -> 0}, {z, n}]; ChiFuncDerivative[state, n, j] = Limit[ChiFuncDerivativeFunc[state, n, j, z], z -> 0]; , {state, {"0+", "0-", "1+", "1-"}}, {n, 0, 4}, {j, 7}]; ); Null; (EpsilonExpandFunc[state_, n_, j_, Nc_, Nb_] := Coefficient[Coefficient[Expand[ChiFuncDerivative[state, n, j]], epc, Nc], epb, Nb]*epc^Nc*epb^Nb; ); Null; (EpsilonExpand[state_, Nc_, Nb_, Nexpand_] := Sum[EpsilonExpandFunc[state, n, j, Nc, Nb]^2, {n, 0, Nexpand}, {j, 7}]; ); Null; UBfunction["HQET321", state_] := EpsilonExpand[state, 0, 0, 3] + EpsilonExpand[state, 1, 0, 2] + EpsilonExpand[state, 0, 1, 2] + (EpsilonExpand[state, 2, 0, 1] /. SubIWcut); Null; UBfunction["HQET210", state_] := EpsilonExpand[state, 0, 0, 2] + EpsilonExpand[state, 1, 0, 1] + EpsilonExpand[state, 0, 1, 1] + (EpsilonExpand[state, 2, 0, 0] /. SubIWcut); Null; (Do[UBformula["HQET321", state] = UBfunction["HQET321", state] /. LigetiFit; UBformula["HQET210", state] = UBfunction["HQET210", state] /. LigetiFit; , {state, {"0+", "0-", "1+", "1-"}}]; )

; 

; 
);

StanFitFunctionIntegrate["BSZ"] := (
CompoundExpression[Clear[FFmodel, Obs]; ]

CompoundExpression[If[FileExistsQ["Stan_inputs/IntFunctions_BSZ.ma"], Print["Already saved. Skip the integration step."], FFmodel = "BSZ"; FFnf = 2; SpecifyModel[FFmodel, FFnf]; PreIntPerform; ObsIntPerform["Belle15", "w"]; ObsIntPerform["BelleII25", "w"]; ObsIntPerform["Belle17", "w"]; ObsIntPerform["Belle17", "cosL"]; ObsIntPerform["Belle17", "cosV"]; ObsIntPerform["Belle17", "χ"]; Do[Obs["Belle18", dist, n] = Expand[(ObsScale["Belle18"]/ObsScale["Belle17"])*Obs["Belle17", dist, n]]; , {dist, {"w", "cosL", "cosV", "χ"}}, {n, 10}]; ObsIntPerform["Belle23", "w"]; ObsIntPerform["Belle23", "cosL"]; ObsIntPerform["Belle23", "cosV"]; ObsIntPerform["Belle23", "χ"]; ObsIntPerform["BelleII23", "w"]; ObsIntPerform["BelleII23", "cosL"]; ObsIntPerform["BelleII23", "cosV"]; ObsIntPerform["BelleII23", "χ"]; ObsIntPerform["BrRatio2025", "D"]; ObsIntPerform["BrRatio2025", "Dst"]; ObsIntPerform["BrRatioTau", "D"]; ObsIntPerform["BrRatioTau", "Dst"]; Save["Stan_inputs/IntFunctions_BSZ.ma", Obs]; ]]
);

StanFitFunctionIntegrate["BGL"] := (
CompoundExpression[Clear[FFmodel, Obs]; ]

CompoundExpression[If[FileExistsQ["Stan_inputs/IntFunctions_BGL.ma"], Print["Already saved. Skip the integration step."], FFmodel = "BGL"; FFnf = 2; SpecifyModel[FFmodel, FFnf]; PreIntPerform; ObsIntPerform["Belle15", "w"]; ObsIntPerform["BelleII25", "w"]; ObsIntPerform["Belle17", "w"]; ObsIntPerform["Belle17", "cosL"]; ObsIntPerform["Belle17", "cosV"]; ObsIntPerform["Belle17", "χ"]; Do[Obs["Belle18", dist, n] = Expand[(ObsScale["Belle18"]/ObsScale["Belle17"])*Obs["Belle17", dist, n]]; , {dist, {"w", "cosL", "cosV", "χ"}}, {n, 10}]; ObsIntPerform["Belle23", "w"]; ObsIntPerform["Belle23", "cosL"]; ObsIntPerform["Belle23", "cosV"]; ObsIntPerform["Belle23", "χ"]; ObsIntPerform["BelleII23", "w"]; ObsIntPerform["BelleII23", "cosL"]; ObsIntPerform["BelleII23", "cosV"]; ObsIntPerform["BelleII23", "χ"]; ObsIntPerform["BrRatio2025", "D"]; ObsIntPerform["BrRatio2025", "Dst"]; ObsIntPerform["BrRatioTau", "D"]; ObsIntPerform["BrRatioTau", "Dst"]; Save["Stan_inputs/IntFunctions_BGL.ma", Obs]; ]]
);

StanFitFunctionIntegrate["HQET"] := (
CompoundExpression[Clear[FFmodel, HQETmodel, HQETname]; ]

CompoundExpression[If[FileExistsQ["Stan_inputs/IntFunctions_HQET.ma"], Print["Already saved. Skip the integration step."], FFmodel = "HQET"; HQETmodel = "3/2/1"; HQETset[HQETmodel]; SpecifyModel[FFmodel, HQETmodel]; PreIntPerformHQET; ObsIntPerformHQET["Belle15", "w"]; ObsIntPerformHQET["BelleII25", "w"]; ObsIntPerformHQET["Belle17", "w"]; ObsIntPerform["Belle17", "cosL"]; ObsIntPerform["Belle17", "cosV"]; ObsIntPerform["Belle17", "χ"]; Do[Obs["Belle18", dist, n] = Expand[(ObsScale["Belle18"]/ObsScale["Belle17"])*Obs["Belle17", dist, n]]; , {dist, {"w", "cosL", "cosV", "χ"}}, {n, 10}]; ObsIntPerformHQET["Belle23", "w"]; ObsIntPerform["Belle23", "cosL"]; ObsIntPerform["Belle23", "cosV"]; ObsIntPerform["Belle23", "χ"]; ObsIntPerformHQET["BelleII23", "w"]; ObsIntPerform["BelleII23", "cosL"]; ObsIntPerform["BelleII23", "cosV"]; ObsIntPerform["BelleII23", "χ"]; ObsIntPerformHQET["BrRatio2025", "D"]; ObsIntPerformHQET["BrRatio2025", "Dst"]; ObsIntPerformHQET["BrRatioTau", "D"]; ObsIntPerformHQET["BrRatioTau", "Dst"]; Save["Stan_inputs/IntFunctions_HQET.ma", Obs]; ]]
);

StanFitFunctionIntegrate[ffModel_String] := (
  Message[StanFitFunctionIntegrate::unkff, ffModel];
  $Failed
);

StanFitFunctionIntegrate::unkff = "Unknown FF model `1`. Expected BSZ, BGL, or HQET.";

StanFitFunctionCommonSetup[] := (
  FFpriorSet["simple", nf_] := Table[" ~ normal(0, 1);", {i, 30}];
  FFpriorSet["narrow", nf_] := Table[" ~ normal(0, 0.1);", {i, 30}];
  FFpriorSet["wide", nf_] := Table[" ~ normal(0, 10);", {i, 30}];
  FFpriorSet["wide2", nf_] := Table[" ~ normal(0, 100);", {i, 30}];

  FFparaList["BSZ", 2] = {"afp0", "afp1", "afp2", "af01", "af02", "aV0", "aV1", "aV2", "aA10", "aA11", "aA12", "aA121", "aA122", "aA01", "aA02"};
  FFparaListNP["T", "BSZ", 2] = {"afT0", "afT1", "afT2", "aT10", "aT11", "aT12", "aT21", "aT22", "aT231", "aT232"};
  FFparaList["BSZ", 1] = {"afp0", "afp1", "af01", "aV0", "aV1", "aA10", "aA11", "aA121", "aA01"};
  FFparaListNP["T", "BSZ", 1] = {"afT0", "afT1", "aT10", "aT11", "aT21", "aT231"};

  FFparaList["BGL", 2] = {"bfp0", "bfp1", "bfp2", "bf01", "bf02", "bg0", "bg1", "bg2", "bf0", "bf1", "bf2", "bF11", "bF12", "bF21", "bF22"};
  FFparaListNP["T", "BGL", 2] = {"bfT0", "bfT1", "bfT2", "bT10", "bT11", "bT12", "bT21", "bT22", "bT231", "bT232"};
  FFparaList["BGL", 1] = {"bfp0", "bfp1", "bf01", "bg0", "bg1", "bf0", "bf1", "bF11", "bF21"};
  FFparaListNP["T", "BGL", 1] = {"bfT0", "bfT1", "bT10", "bT11", "bT21", "bT231"};

  FFparaList["HQET", 2] = {"cxi1", "cxi2", "cxi3", "cet0", "cet1", "cet2", "cci20", "cci21", "cci22", "cci31", "cci32", "cel10", "cel11", "cel20", "cel21", "cel30", "cel31", "cel40", "cel41", "cel50", "cel51", "cel60", "cel61"};
  FFparaList["HQET", 1] = {"cxi1", "cxi2", "cet0", "cet1", "cci20", "cci21", "cci31", "cel10", "cel20", "cel30", "cel40", "cel50", "cel60"};

  ParaDef["full", ff_, nf_] := Table["real " <> FFparaList[ff, nf][[ii]] <> ";\n", {ii, Length[FFparaList[ff, nf]]}];
  ParaDefNP["T", ff_, nf_] := Table["real " <> FFparaListNP["T", ff, nf][[ii]] <> ";\n", {ii, Length[FFparaListNP["T", ff, nf]]}];
  ParaDefNP["T", "HQET", nf_] := {};
  ParaPrior[ffprior_, ff_, nf_] := Table[FFparaList[ff, nf][[ii]] <> FFpriorSet[ffprior, nf][[ii]] <> "\n", {ii, Length[FFparaList[ff, nf]]}];
  ParaPriorNP["T", ffprior_, ff_, nf_] := Table[FFparaListNP["T", ff, nf][[ii]] <> FFpriorSet[ffprior, nf][[ii]] <> "\n", {ii, Length[FFparaListNP["T", ff, nf]]}];
  ParaPriorNP["T", ffprior_, "HQET", nf_] := {};

  ParaPrior["QCDSR", "HQET", 1] = {
    "cxi1 ~ normal(0, 1);\n", "cxi2 ~ normal(0, 1);\n",
    "cet0 ~ normal(0.62, 0.12);\n", "cet1 ~ normal(0.04, 0.03);\n",
    "cci20 ~ normal(-0.06, 0.02);\n", "cci21 ~ normal(0.00, 0.02);\n",
    "cci31 ~ normal(0.04, 0.03);\n",
    "cel10 ~ normal(0, 1);\n", "cel20 ~ normal(0, 1);\n", "cel30 ~ normal(0, 1);\n",
    "cel40 ~ normal(0, 1);\n", "cel50 ~ normal(0, 1);\n", "cel60 ~ normal(0, 1);\n"
  };
  ParaPrior["QCDSR", "HQET", 2] = {
    "cxi1 ~ normal(0, 1);\n", "cxi2 ~ normal(0, 1);\n", "cxi3 ~ normal(0, 1);\n",
    "cet0 ~ normal(0.62, 0.12);\n", "cet1 ~ normal(0.04, 0.03);\n", "cet2 ~ normal(0.05, 0.07);\n",
    "cci20 ~ normal(-0.06, 0.02);\n", "cci21 ~ normal(0.00, 0.02);\n", "cci22 ~ normal(-0.01, 0.02);\n",
    "cci31 ~ normal(0.04, 0.03);\n", "cci32 ~ normal(-0.03, 0.05);\n",
    "cel10 ~ normal(0, 1);\n", "cel11 ~ normal(0, 1);\n", "cel20 ~ normal(0, 1);\n",
    "cel21 ~ normal(0, 1);\n", "cel30 ~ normal(0, 1);\n", "cel31 ~ normal(0, 1);\n",
    "cel40 ~ normal(0, 1);\n", "cel41 ~ normal(0, 1);\n", "cel50 ~ normal(0, 1);\n",
    "cel51 ~ normal(0, 1);\n", "cel60 ~ normal(0, 1);\n", "cel61 ~ normal(0, 1);\n"
  };

  DstFitData[dstData_] := Switch[dstData, "Wdata", "w", "CosLdata", "cosL", "CosVdata", "cosV", "Chidata", "\[Chi]"];
  DstFitList := {"w", "cosL", "cosV", "\[Chi]"};
  NPsetup[np_] := Switch[np, "SM", {CT -> 0, CVR -> 0, CSR -> 0, CSL -> 0}, "VR", {CT -> 0, CSR -> 0, CSL -> 0}, "SL", {CT -> 0, CVR -> 0, CSR -> 0}, "SR", {CT -> 0, CVR -> 0, CSL -> 0}, "T", {CVR -> 0, CSR -> 0, CSL -> 0}];
  NPpara[np_] := Switch[np, "SM", {Vcb}, "VR", {CVR, Vcb}, "SL", {CSL, Vcb}, "SR", {CSR, Vcb}, "T", {CT, Vcb}];
  NPparaSplit[np_] := Switch[np, "SM", {VcbD, VcbDst}, "VR", {CVR, VcbD, VcbDst}, "SL", {CSL, VcbD, VcbDst}, "SR", {CSR, VcbD, VcbDst}, "T", {CT, VcbD, VcbDst}];

  FFsetup[ffmodel_, np_] := Switch[
    ffmodel,
    "BSZ",
      If[np == "T", BSZfunc = {"fp", "f0", "fT", "V", "A0", "A1", "A12", "T1", "T2", "T23"}, BSZfunc = {"fp", "f0", "V", "A0", "A1", "A12"}];
      FFlist = Flatten[Table[Variables[BSZff[ffmodel <> "_" <> ff, FFnf, qSQ]], {ff, BSZfunc}]] //. Join[{qSQ -> Nothing}, KineticChoice];
      UBdata = ""; UBdef = ""; UBfunc = ""; UBfit = "";,
    "BGL",
      If[np == "T", BGLfunc = {"fp", "f0", "fT", "f", "F1", "F2", "g", "T1", "T2", "T23"}, BGLfunc = {"fp", "f0", "f", "F1", "F2", "g"}];
      FFlist = Flatten[Table[Variables[BGLff[ffmodel <> "_" <> ff, FFnf, qSQ]], {ff, BGLfunc}]] //. Join[{qSQ -> Nothing}, KineticChoice];
      UBdata = ""; UBdef = ""; UBfunc = ""; UBfit = "";,
    "HQET",
      FFlist = Select[Flatten[Table[FFc[list], {list, HQETlist}]] /. NF2off, ! NumberQ[#] &];
      UBdata = ""; UBdef = ""; UBfunc = ""; UBfit = "";
  ];

  TensorForLCSR18[np_] := Module[{},
    If[np == "T",
      LCSR18data = "vector[37] bszLCSR18cent;\n cov_matrix[37] bszLCSR18cov;";
      LCSR18number = 37;
      LCSR18func = Table["LCSR18func[" <> ToString[i] <> "]=" <> LCSR18funclist[[i]] <> ";\n", {i, 37}];
      LCSR18fit = "bszLCSR18cent ~ multi_normal(LCSR18func, bszLCSR18cov);";,
      LCSR18data = "vector[22] bszLCSR18centNoTensor;\n cov_matrix[22] bszLCSR18covNoTensor;";
      LCSR18number = 22;
      LCSR18func = Table["LCSR18func[" <> ToString[i] <> "]=" <> LCSR18funclistNoTensor[[i]] <> ";\n", {i, 22}];
      LCSR18fit = "bszLCSR18centNoTensor ~ multi_normal(LCSR18func, bszLCSR18covNoTensor);";
    ];
  ];

  TensorForHPQCD23[np_] := Module[{},
    If[np == "T",
      HPQCD23data = "vector[35] bglHPQCD23cent;\n cov_matrix[35] bglHPQCD23cov;";
      HPQCD23number = 35;
      HPQCD23func = Table["HPQCD23func[" <> ToString[i] <> "]=" <> HPQCD23funclist[[i]] <> ";\n", {i, 35}];
      HPQCD23fit = "bglHPQCD23cent ~ multi_normal(HPQCD23func, bglHPQCD23cov);";,
      HPQCD23data = "vector[20] bglHPQCD23centNoTensor;\n cov_matrix[20] bglHPQCD23covNoTensor;";
      HPQCD23number = 20;
      HPQCD23func = Table["HPQCD23func[" <> ToString[i] <> "]=" <> HPQCD23funclistNoTensor[[i]] <> ";\n", {i, 20}];
      HPQCD23fit = "bglHPQCD23centNoTensor ~ multi_normal(HPQCD23func, bglHPQCD23covNoTensor);";
    ];
  ];

  ChopSmallTerm[func_, cut_] := Module[{terms, vars, coeff, coeffs, maxcoeff, threshold},
    If[func === 0 || func === {0}, Return[0]];
    terms = If[Head[func] === Plus, List @@ func, {func}];
    vars = Variables[func];
    coeff[t_] := t/Times @@ (vars^Exponent[t, vars]);
    coeffs = N[coeff /@ terms];
    If[coeffs === {} || Max[Abs[coeffs]] == 0, Return[0]];
    maxcoeff = Max[Abs[coeffs]];
    threshold = maxcoeff*cut;
    Total[Select[terms, Abs[N[coeff[#]]] >= threshold &]]
  ];
  NLOpara = {cet1, cet2, cci20, cci21, cci22, cci31, cci32};
  NLOsq = (Times @@ (NLOpara^#) &) /@ Select[Tuples[Range[0, 2], Length[NLOpara]], Total[#] == 2 &];
  NLOcut = Table[i -> 0, {i, NLOsq}];
);

Options[StanPrepareFitFunctionLists] = {
  "FFModel" -> "BGL",
  "FFnf" -> 2,
  "NPModel" -> "SM",
  "HQETModel" -> "3/2/1",
  "CutControl" -> 10^-4,
  "InputDirectory" -> "Stan_inputs"
};

StanPrepareFitFunctionLists[OptionsPattern[]] := Module[
  {inputDirectory, intFile, ffExpr, qForD, qForDst},
  FFmodel = OptionValue["FFModel"];
  FFnf = OptionValue["FFnf"];
  NPmodel = OptionValue["NPModel"];
  HQETmodel = OptionValue["HQETModel"];
  CutControl = OptionValue["CutControl"];
  inputDirectory = OptionValue["InputDirectory"];

  If[FFmodel == "HQET",
    HQETset[HQETmodel];
    SpecifyModel[FFmodel, HQETmodel];,
    SpecifyModel[FFmodel, FFnf];
  ];
  SpecifyKinetic["On"];

  intFile = FileNameJoin[{inputDirectory, "IntFunctions_" <> FFmodel <> ".ma"}];
  If[! FileExistsQ[intFile],
    Message[StanPrepareFitFunctionLists::noint, intFile];
    Return[$Failed];
  ];
  Get[intFile];

  FFsetup[FFmodel, NPmodel];
  ParamList = DeleteDuplicates[Join[FFlist, NPpara[NPmodel]]];
  StanForm = Table[
    ToString[ParamList[[i]]] <> "**2" -> ToString[ParamList[[i]]] <> "*" <> ToString[ParamList[[i]]],
    {i, Length[ParamList]}
  ];
  MathToStan[func_] := StringReplace[
    StringReplace[
      ToString[FortranForm[N[func]], PageWidth -> Infinity],
      WhitespaceCharacter.. -> ""
    ],
    StanForm
  ];
  ffExpr[ff_, arg_] := ToExpression[FFmodel <> "ff"][FFmodel <> "_" <> ff, FFnf, arg];
  qForD[wval_] := If[FFmodel == "HQET", wval, qSQfuncD[wval]];
  qForDst[wval_] := If[FFmodel == "HQET", wval, qSQfuncDst[wval]];

  MILC15funclist = Flatten[Table[
    MathToStan[ChopSmallTerm[Expand[ffExpr[ff, qForD[wval]] /. NLOcut //. KineticReplace], CutControl]],
    {ff, {"fp", "f0"}}, {wval, Bin["MILC15", "w"]}
  ]];
  MILC21funclist = Flatten[Table[
    MathToStan[ChopSmallTerm[Expand[ffExpr[ff, qForDst[wval]] /. NLOcut //. KineticReplace], CutControl]],
    {ff, {"f", "F1", "F2", "g"}}, {wval, Bin["MILC21", "w"]}
  ]];
  JLQCD23funclist = Flatten[Table[
    MathToStan[ChopSmallTerm[Expand[ffExpr[ff, qForDst[wval]] /. NLOcut //. KineticReplace], CutControl]],
    {ff, {"f", "F1", "F2", "g"}}, {wval, Bin["bglJLQCD23", "w"]}
  ]];
  HPQCD23funclist = Flatten[Table[
    MathToStan[ChopSmallTerm[Expand[Limit[ffExpr[ff, qForDst[w]], w -> wlim] /. NLOcut //. KineticReplace], CutControl]],
    {ff, {"f", "F1", "F2", "g", "T1", "T2", "T23"}}, {wlim, Bin["bglHPQCD23", "w"]}
  ]];
  HPQCD23funclistNoTensor = Delete[HPQCD23funclist, Table[{i}, {i, 21, 35}]];

  FForder["LCSR18"] = {"fp", "f0", "fT", "A0", "A1", "A2", "V", "T1", "T2", "T23"};
  LCSR18funclist = Delete[
    Flatten[Table[
      MathToStan[ChopSmallTerm[
        Expand[Limit[
          ffExpr[ff, If[FFmodel == "HQET", w, qSQ]] /. NLOcut //. KineticReplace,
          If[FFmodel == "HQET", w -> If[StringContainsQ[ff, "f"], wFuncD[qSQlim], wFuncDst[qSQlim]], qSQ -> qSQlim]
        ]],
        CutControl
      ]],
      {ff, FForder["LCSR18"]}, {qSQlim, Bin["LCSR18", "qSQ"]}
    ]],
    {{8}, {16}, {36}}
  ];
  LCSR18funclistNoTensor = Delete[LCSR18funclist, Table[{i}, {i, TensorList}]];
  TensorForLCSR18[NPmodel];
  TensorForHPQCD23[NPmodel];

  EpDst[qSQ_] = (mB^2 + mDst^2 - qSQ)/(2*mB);
  EpD[qSQ_] = (mB^2 + mD^2 - qSQ)/(2*mB);
  LCSR23factor["fp", qSQ_] = 1;
  LCSR23factor["f0", qSQ_] = mB/(2*EpD[qSQ]);
  LCSR23factor["V", qSQ_] = mB/(mB + mDst);
  LCSR23factor["A0", qSQ_] = mDst/EpDst[qSQ];
  LCSR23factor["A1", qSQ_] = (mB + mDst)/(2*EpDst[qSQ]);
  LCSR23factor["A2", qSQ_] = (mB - mDst)/mB;
  FForder["LCSR23"] = {"fp", "f0", "V", "A0", "A1"};
  LCSR23funclist = Delete[
    If[FFmodel == "HQET",
      Join[
        Flatten[Table[
          MathToStan[ChopSmallTerm[Expand[LCSR23factor[ff, qSQlim]*Limit[HQETff["HQET_" <> ff, FFnf, w] /. NLOcut //. KineticReplace, w -> If[StringContainsQ[ff, "f"], wFuncD[qSQlim], wFuncDst[qSQlim]]]], CutControl]],
          {ff, FForder["LCSR23"]}, {qSQlim, Bin["LCSR23", "qSQ"]}
        ]],
        Table[
          MathToStan[ChopSmallTerm[Expand[Limit[LCSR23factor["A1", qSQ]*HQETff["HQET_A1", FFnf, w] - LCSR23factor["A2", qSQ]*HQETff["HQET_A2", FFnf, w] /. NLOcut //. KineticReplace, w -> wFuncDst[qSQ]]], CutControl]],
          {qSQ, Bin["LCSR23", "qSQ"]}
        ]
      ],
      Join[
        Flatten[Table[
          MathToStan[ChopSmallTerm[Expand[LCSR23factor[ff, qSQlim]*Limit[ffExpr[ff, qSQ] /. NLOcut //. KineticReplace, qSQ -> qSQlim]], CutControl]],
          {ff, FForder["LCSR23"]}, {qSQlim, Bin["LCSR23", "qSQ"]}
        ]],
        Table[
          MathToStan[ChopSmallTerm[Expand[LCSR23factor["A1", qSQlim]*Limit[ffExpr["A1", qSQ] //. KineticReplace, qSQ -> qSQlim] - LCSR23factor["A2", qSQlim]*Limit[ffExpr["A2", qSQ] //. KineticReplace, qSQ -> qSQlim] /. NLOcut //. KineticReplace], CutControl]],
          {qSQlim, Bin["LCSR23", "qSQ"]}
        ]
      ]
    ],
    {{10}, {22}}
  ];

  If[FFmodel == "HQET",
    UBfunclist = Flatten[Table[MathToStan[ChopSmallTerm[Expand[UBformula[FFmodel <> HQETname, state]], CutControl]], {state, {"0+", "0-", "1+", "1-"}}]];
    UBfunc = Table["UBfunc[" <> ToString[i] <> "]=" <> UBfunclist[[i]] <> ";\n", {i, 4}];,
    UBfunc = "";
  ];

  BelleFITfunclist["Br"] = Flatten[Table[
    MathToStan[ChopSmallTerm[Expand[Obs["BrRatio2025", bin] /. NLOcut //. KineticReplace /. NPsetup[NPmodel]], CutControl]],
    {bin, {"D", "Dst"}}
  ]];
  BelleFITfunclist["BrTau"] = Flatten[Table[
    MathToStan[ChopSmallTerm[Expand[Obs["BrRatioTau", bin] /. NLOcut //. KineticReplace /. NPsetup[NPmodel]], CutControl]],
    {bin, {"D", "Dst"}}
  ]];

  ReRescale[bin_] := If[bin == 10, wMax["D"] - 1.54, 0.06];
  BDlnuRescale["Belle15", bin_] := If[bin == 10, 0.06/(wMax["D"] - 1.54), 1];
  BelleFITfunclist["D"] = Flatten[Table[
    MathToStan[ChopSmallTerm[Expand[ReRescale[bin]*BDlnuRescale["Belle15", bin]*Obs["Belle15", "w", bin] /. NLOcut //. KineticReplace /. NPsetup[NPmodel]], CutControl]],
    {bin, 10}
  ]];
  BelleFITfunctotal["D"] = MathToStan[ChopSmallTerm[Expand[Sum[ReRescale[bin]*BDlnuRescale["Belle15", bin]*Obs["Belle15", "w", bin] /. NLOcut //. KineticReplace /. NPsetup[NPmodel], {bin, 10}]], CutControl]];
  BelleFITfunclist["w"] = Flatten[Table[MathToStan[ChopSmallTerm[Expand[Obs["Belle17", "w", bin] /. NLOcut //. KineticReplace /. NPsetup[NPmodel]], CutControl]], {bin, 10}]];
  BelleFITfunclist["cosL"] = Flatten[Table[MathToStan[ChopSmallTerm[Expand[Obs["Belle17", "cosL", bin] /. NLOcut //. KineticReplace /. NPsetup[NPmodel]], CutControl]], {bin, 10}]];
  BelleFITfunclist["cosV"] = Flatten[Table[MathToStan[ChopSmallTerm[Expand[Obs["Belle17", "cosV", bin] /. NLOcut //. KineticReplace /. NPsetup[NPmodel]], CutControl]], {bin, 10}]];
  BelleFITfunclist["chi"] = Flatten[Table[MathToStan[ChopSmallTerm[Expand[Obs["Belle17", "\[Chi]", bin] /. NLOcut //. KineticReplace /. NPsetup[NPmodel]], CutControl]], {bin, 10}]];
  BelleFITfunctotal["Dst"] = MathToStan[ChopSmallTerm[Expand[Sum[Obs["Belle17", "w", bin] /. NLOcut //. KineticReplace /. NPsetup[NPmodel], {bin, 10}]], CutControl]];
  BelleFITfunclist["Dst"] = Join[BelleFITfunclist["w"], BelleFITfunclist["cosL"], BelleFITfunclist["cosV"], BelleFITfunclist["chi"]];
  True
];

StanPrepareFitFunctionLists::noint = "Integrated function file was not found: `1`.";

Options[StanGenerateFitFunctions] = {
  "FFModel" -> "BGL",
  "FFnf" -> 2,
  "NPModel" -> "SM",
  "HQETModel" -> "3/2/1",
  "CutControl" -> 10^-4,
  "InputDirectory" -> "Stan_inputs"
};

StanGenerateFitFunctions[OptionsPattern[]] := Module[{ffModel},
  ffModel = OptionValue["FFModel"];
  If[! DirectoryQ[OptionValue["InputDirectory"]],
    CreateDirectory[OptionValue["InputDirectory"], CreateIntermediateDirectories -> True]
  ];
  StanFitFunctionSetup[];
  StanFitFunctionCommonSetup[];
  If[ffModel === "HQET", StanFitFunctionUBSetup[]];
  StanFitFunctionIntegrate[ffModel];
  StanPrepareFitFunctionLists[
    "FFModel" -> ffModel,
    "FFnf" -> OptionValue["FFnf"],
    "NPModel" -> OptionValue["NPModel"],
    "HQETModel" -> OptionValue["HQETModel"],
    "CutControl" -> OptionValue["CutControl"],
    "InputDirectory" -> OptionValue["InputDirectory"]
  ];
];

If[Length[Rest[$ScriptCommandLine]] > 0 &&
    StringMatchQ[$InputFileName, ___ ~~ "Stan_FitFunction.wls"] &&
    MemberQ[{"BSZ", "BGL", "HQET"}, First[Rest[$ScriptCommandLine]]],
  StanGenerateFitFunctions["FFModel" -> First[Rest[$ScriptCommandLine]]];
];
