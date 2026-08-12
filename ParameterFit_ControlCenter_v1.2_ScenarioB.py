#source stan_env/bin/activate (mac)
#source stan_env/Scripts/activate (win)
from cmdstanpy import CmdStanModel
from pathlib import Path
import json
import os
import arviz as az
import numpy as np
import pandas as pd
from datetime import datetime
#from Rhat import RhatSum

# Stan model
output = "full"
name_analysis = "CasePDG-UB_SM_BGLnf2_simple"
#name_analysis = "CasePDG-UB_SM_BGLnf2_simple"
#name_analysis = "CasePDG-UB_VR_HQET321_QCDSR"

name_path = Path.cwd()
name_model = "stancode_" + name_analysis + ".stan"
name_param = "stancode_" + name_analysis + "_param.json"
name_data = "standata_" + name_analysis + ".json"
name_output = "StanFit_" + name_analysis

# Set your directory for the above three files (设置上述三个文件的目录)
name_fullmodel = name_path / "Stan_model" / name_model
name_fullparam = name_path / "Stan_model" / name_param
name_fulldata = name_path / "Stan_data" / name_data

##############################
# Input data
with open(name_fulldata, "r") as datataken:
    data = json.load(datataken)

MILC15cent = np.array(data["MILC15cent"])
MILC15cov = np.array(data["MILC15cov"])

bglMILC21cent = np.array(data["bglMILC21cent"])
bglMILC21cov = np.array(data["bglMILC21cov"])

bglJLQCD23cent = np.array(data["bglJLQCD23cent"])
bglJLQCD23cov = np.array(data["bglJLQCD23cov"])

bglHPQCD23cent = np.array(data["bglHPQCD23cent"])
bglHPQCD23cov = np.array(data["bglHPQCD23cov"])

bglHPQCD23centNoTensor = np.array(data["bglHPQCD23centNoTensor"])
bglHPQCD23covNoTensor = np.array(data["bglHPQCD23covNoTensor"])

bszLCSR18cent = np.array(data["bszLCSR18cent"])
bszLCSR18cov = np.array(data["bszLCSR18cov"])

bszLCSR18centNoTensor = np.array(data["bszLCSR18centNoTensor"])
bszLCSR18covNoTensor = np.array(data["bszLCSR18covNoTensor"])

bszLCSR23cent = np.array(data["bszLCSR23cent"])
bszLCSR23cov = np.array(data["bszLCSR23cov"])

BrRatio2025cent = np.array(data["BrRatio2025cent"])
BrRatio2025cov = np.array(data["BrRatio2025cov"])

BelleFITcent = np.array(data["BelleFITcent_Dst"])
BelleFITcov = np.array(data["BelleFITcov_Dst"])
BelleRawFITcent = np.array(data["BelleRawFITcent_Dst"])
BelleRawFITcov = np.array(data["BelleRawFITcov_Dst"])

BelleFITDcent = np.array(data["BelleFITcent_D"])
BelleFITDcov = np.array(data["BelleFITcov_D"])
BelleRawFITDcent = np.array(data["BelleRawFITcent_D"])
BelleRawFITDcov = np.array(data["BelleRawFITcov_D"])


stan_data = {
    'MILC15cent': MILC15cent,
    'MILC15cov': MILC15cov,

    'bglMILC21cent': bglMILC21cent,
    'bglMILC21cov': bglMILC21cov,

    'bglJLQCD23cent': bglJLQCD23cent,
    'bglJLQCD23cov': bglJLQCD23cov,

    'bszLCSR18cent': bszLCSR18cent,
    'bszLCSR18cov': bszLCSR18cov,
    'bszLCSR18centNoTensor': bszLCSR18centNoTensor,
    'bszLCSR18covNoTensor': bszLCSR18covNoTensor,

    'bglHPQCD23cent': bglHPQCD23cent,
    'bglHPQCD23cov': bglHPQCD23cov,
    'bglHPQCD23centNoTensor': bglHPQCD23centNoTensor,
    'bglHPQCD23covNoTensor': bglHPQCD23covNoTensor,

    'bszLCSR23cent': bszLCSR23cent,
    'bszLCSR23cov': bszLCSR23cov,

    'BrRatio2025cent': BrRatio2025cent,
    'BrRatio2025cov': BrRatio2025cov,

    'BelleFITcent_w': BelleFITcent[0:10],
    'BelleFITcent_cosL': BelleFITcent[10:20],
    'BelleFITcent_cosV': BelleFITcent[20:30],
    'BelleFITcent_chi': BelleFITcent[30:40],
    'BelleFITcov_w': BelleFITcov[0:10, 0:10],
    'BelleFITcov_cosL': BelleFITcov[10:20, 10:20],
    'BelleFITcov_cosV': BelleFITcov[20:30, 20:30],
    'BelleFITcov_chi': BelleFITcov[30:40, 30:40],
    'BelleRawFITcent_Dst': BelleRawFITcent,
    'BelleRawFITcov_Dst': BelleRawFITcov,

    'BelleFITcent_D': BelleFITDcent,
    'BelleFITcov_D': BelleFITDcov,
    'BelleRawFITcent_D': BelleRawFITDcent,
    'BelleRawFITDcov_D': BelleRawFITDcov,

    'UBcent': np.array([0.0, 0.0, 0.0, 0.0]),
    'UBcov': np.diag([4.1e-5, 1.9e-4, 1.65e-5, 1.41e-5]),
    'KNcent': np.array([0.0, 0.0]),
    'KNcov': np.diag([2.5e-3, 2.5e-3])
}

print("-------------------------------------------------------")
print("-------------------------------------------------------")
print(" [Fit Analysis Code for fit scenario(B):  version 1.2] ")
print("-------------------------------------------------------")
print("-------------------------------------------------------")
print(" ")
print(" ")

##############################
# Define parameters
with open(name_fullparam, "r") as paramtaken:
    FitPara = json.load(paramtaken)

##############################
# Compile stan model
model = CmdStanModel(stan_file=name_fullmodel)

##############################
# Pre-fit run
print("-------------")
print(" Pre fit run ")
print("-------------")

pre_fit = model.sample(
    data=stan_data,
    chains = 10,
    iter_warmup=500,
    iter_sampling=200,
    adapt_delta=0.95,
    max_treedepth=12
)
pre_df = pre_fit.draws_pd()

# Good chain check
for c in sorted(pre_df["chain__"].unique()):
    sub = pre_df[pre_df["chain__"] == c]["lp__"]
    print(f"Chain {c}:  sample_length = {len(sub)}, lp_max   = {sub.max():.2f}")

threshold = pre_df["lp__"].max() - 200
good_chains = []

for c in sorted(pre_df["chain__"].unique()):
    sub = pre_df[pre_df["chain__"] == c]["lp__"]
    if sub.mean() > threshold:
        good_chains.append(int(c))

print(f"Good chains: {[int(x) for x in good_chains]}")

# Set initial inputs (to avoid bad chain)
good_df = pre_df[pre_df["lp__"] > threshold]

good_init = good_df.iloc[0][FitPara].to_dict()
inits_set = [good_init] * 3

##############################
# Fit analysis run
print("-------------------")
print(" Ffit analysis run ")
print("-------------------")

fit = model.sample(
    data=stan_data,
    chains=3,
    iter_warmup=3000,
    iter_sampling=6000,
    adapt_delta=0.99,
    max_treedepth=15,
    seed=42,
    inits=inits_set
)

##############################
# Collect fit results
df = fit.draws_pd()

print("------------------")
print(" Bad chain check: ")
print("------------------")

for c in sorted(df["chain__"].unique()):
    sub = df[df["chain__"] == c]["lp__"]
    print(f"Chain {c}:  sample_length = {len(sub)}, lp_max   = {sub.max():.2f}")

# Reject bad chains
good_chains = []
for c in sorted(df["chain__"].unique()):
    sub = df[df["chain__"] == c]["lp__"]
    if sub.mean() > df["lp__"].max() - 200:
        good_chains.append(c)
df_good = df[df["chain__"].isin(good_chains)]

print(f"Good chains: {[int(x) for x in good_chains]}")
print(" ")
print("Bad chains are rejected in output. 2 good chains are needed at least.")
print(" ")

obs_option = {
    "full": ["MILC15func", "MILC21func", "HPQCD23func", "JLQCD23func", "LCSR18func", "LCSR23func", 
             "rawBrRatio_func", "rawBelleD_func", "rawBelleDst_w_func", "rawBelleDst_cosL_func", "rawBelleDst_cosV_func", "rawBelleDst_chi_func", 
             "BrRatio_func", "BelleD_func", "BelleDst_w_func", "BelleDst_cosL_func", "BelleDst_cosV_func", "BelleDst_chi_func", 
             "RD_func", "chi_sq", "Vcb_Br"],     
    "light": ["MILC15func", "MILC21func", "HPQCD23func", "JLQCD23func", "LCSR18func", "LCSR23func"]
}

fitted_obs = obs_option.get(output, [])
#FitObs = [col for col in df.columns if any(v in col for v in fitted_obs)]
FitObs = [
    col for col in df.columns
    if any(col == v or col.startswith(v + "[") for v in fitted_obs)
]
FitAll = FitPara + FitObs + ["lp__", "chi_sq_total", "Vcb_mean"]
params_df = df_good[FitAll]


# CSVに保存
timestamp = datetime.today().strftime('%m%d_%H%M')
filename = f"{name_output}_{timestamp}.csv"
name_fulloutput = name_path / "Stan_output" / filename

params_df.to_csv(name_fulloutput, index=False)


##############################
# Display the result
print("---------------------")
print(" Display Fit Result: ")
print("---------------------")
print(params_df.head())
print(" ")

# rhatを表示
print("-------------------------------")
print(" R-hat values (FF parameters): ")
print("-------------------------------")

idata = az.from_cmdstanpy(fit)
good_chains_idata = [int(c) - 1 for c in good_chains]
good_idata = idata.sel(chain=good_chains_idata)
summ = az.summary(good_idata, var_names=FitPara, round_to=4)
print(summ)
print(" ")

if output == "full":
    print("-------------------------------------")
    print(" R-hat values (RD, Br, Vcb, chi_sq): ")
    print("-------------------------------------")

    summ2 = az.summary(good_idata, var_names=["RD_func", "BrRatio_func", "Vcb_Br", "Vcb_mean", "chi_sq_total"], round_to=4)
    print(summ2)
print(" ")

print("----------------------")
print("  Fit analysis done!  ")
print("----------------------")



