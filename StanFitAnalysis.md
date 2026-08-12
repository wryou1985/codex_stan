# Stan Fit Analysis 手順メモ

このメモは、`Stan_FitAnalysis.py` を正規入口として Stan fit analysis を実行するための手順です。別スレッドで作業を再開するときは、inline Python runner ではなく、このファイルと `Stan_FitAnalysis.py` を参照してください。

## 実行場所と環境

作業ディレクトリは必ず `mycode_MetaFitAnalysis` にする。

```bash
cd /Users/ryoutaro/Documents/Codex/Project_VcbFit/mycode_MetaFitAnalysis
source /Users/ryoutaro/stan_env/bin/activate
python Stan_FitAnalysis.py
```

`Stan_FitAnalysis.py` は Python/cmdstanpy で実行する。通常は `wolframscript` を呼ばない。

## 解析設定

実行前に `Stan_FitAnalysis.py` 冒頭の control block を直接編集する。

```python
ANALYSIS_NAME = "CasePDG-noUB_SM_BGLnf1"
FF_MODEL = "BGL"
FF_NF = 1
NP_MODEL = "SM"
HQET_MODEL = "3/2/1"
UB_MODEL = "no-UB"
CUT_CONTROL = 1.0e-4
BELLE_FIT_OPTION = "sc"
```

主な意味は以下。

- `ANALYSIS_NAME`: 出力フォルダ名と最終CSV名に使う解析名。ケースごとに必ず変える。
- `FF_MODEL`: `"BGL"`, `"BSZ"`, `"HQET"` のいずれか。
- `FF_NF`: form factor 展開次数。例: `1`, `2`。
- `NP_MODEL`: 例: `"SM"`。
- `HQET_MODEL`: HQET の場合に使う模型指定。例: `"3/2/1"`, `"2/1/0"`。BGL/BSZ では実質使わない。
- `UB_MODEL`: `"no-UB"`, `"soft-UB"`, `"hard-UB"`。`effective_ub_model()` により、UB は HQET の場合だけ有効になり、BGL/BSZ では自動的に `"no-UB"` になる。
- `CUT_CONTROL`: fit function の小さい項を落とす threshold。
- `BELLE_FIT_OPTION`: BelleFIT distribution data の選択。`"sc"`, `"re_small"`, `"re_large"`, `"de"`, `"pdg"`, `"rms"`。

取り込むデータセットは `SELECTED_DATA` で指定する。fit procedure (B-PDG) の標準セットは現在のリストに対応している。

```python
SELECTED_DATA = [
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
    "BelleFITDstChi",
]
```

`BrRatio2025` は通常の selected data には入れない。generator 側で常に generated quantities の `Vcb_Br`, `Vcb_mean`, `BrRatio_func` などに使われる設計。

## MCMC 設定

pre-fit は本番fitの初期値を探すための短い run。

```python
PRE_FIT_SETTINGS = {
    "chains": 10,
    "iter_warmup": 500,
    "iter_sampling": 200,
    "adapt_delta": 0.95,
    "max_treedepth": 12,
}
```

遅いpre-fit chainを捨てるため、pre-fitだけ chain ごとに独立実行される。

```python
PRE_FIT_CHAIN_TIMEOUT = 600.0
PRE_FIT_PARALLEL_CHAINS = None
PRE_FIT_MIN_SUCCESSFUL_CHAINS = 1
PRE_FIT_LP_THRESHOLD_DROP = 200.0
```

- `PRE_FIT_CHAIN_TIMEOUT`: 1 chain あたりの制限時間。秒単位。`None` にすると無効。
- `PRE_FIT_PARALLEL_CHAINS`: 同時実行するpre-fit chain数。`None` なら `min(chains, os.cpu_count())`。
- `PRE_FIT_MIN_SUCCESSFUL_CHAINS`: timeoutされず成功する必要がある最小chain数。
- `PRE_FIT_LP_THRESHOLD_DROP`: pre-fitの最大 `lp__` からこの範囲内を「良い初期値候補」とみなす幅。現在は候補の中から最大 `lp__` の1点を本番fit初期値に選ぶ。

本番fitは以下の設定。

```python
FIT_SETTINGS = {
    "chains": 3,
    "iter_warmup": 3000,
    "iter_sampling": 6000,
    "adapt_delta": 0.99,
    "max_treedepth": 15,
    "seed": 42,
}
```

## 実行時の内部フロー

`python Stan_FitAnalysis.py` は `run_fit()` を実行する。

1. `build_generator_config()` で解析設定を作る。
2. `Stan_inputs/config_<ANALYSIS_NAME>.json` を保存する。
3. `generate_from_config()` で `.stan`, `standata_*.json`, `*_param.json` を生成する。
4. 生成物は `Stan_outputs/<ANALYSIS_NAME>/` に保存される。
5. `standata_*.json` を読み込む。
6. `CmdStanModel(stan_file=stan_file)` で Stan model を用意する。必要なら cmdstanpy がコンパイルする。
7. `run_pre_fit()` でpre-fitを chain ごとに走らせる。timeoutしたchainは捨てる。
8. pre-fit結果から最大 `lp__` の点を本番fitの初期値に使う。
9. `FIT_SETTINGS` で本番MCMCを走らせる。
10. `good_chain_ids()` で `lp__` が悪いchainを落とす。
11. FF parameters, selected observables, `chi_sq_total`, `Vcb_mean` などを抽出する。
12. `Stan_outputs/StanFit_<ANALYSIS_NAME>_<MMDD_HHMM>.csv` に最終結果CSVを保存する。

## 出力ファイル

代表的な出力は以下。

```text
Stan_inputs/config_<ANALYSIS_NAME>.json
Stan_outputs/<ANALYSIS_NAME>/standata_<ANALYSIS_NAME>.json
Stan_outputs/<ANALYSIS_NAME>/stancode_<ANALYSIS_NAME>.stan
Stan_outputs/<ANALYSIS_NAME>/stancode_<ANALYSIS_NAME>_param.json
Stan_outputs/<ANALYSIS_NAME>/cmdstan_prefit_<MMDD_HHMM>/chain_<N>/...
Stan_outputs/StanFit_<ANALYSIS_NAME>_<MMDD_HHMM>.csv
```

`StanFit_*.csv` が解析後に見る主な結果ファイル。pre-fit の raw CmdStan CSV は `Stan_outputs/<ANALYSIS_NAME>/cmdstan_prefit_<MMDD_HHMM>/chain_<N>/` 以下に残る。

## 正しい実行例

BGL nf=1 SM を実行する場合:

```python
ANALYSIS_NAME = "CasePDG-noUB_SM_BGLnf1"
FF_MODEL = "BGL"
FF_NF = 1
NP_MODEL = "SM"
HQET_MODEL = "3/2/1"
UB_MODEL = "no-UB"
BELLE_FIT_OPTION = "sc"
```

BSZ nf=1 SM を実行する場合:

```python
ANALYSIS_NAME = "CasePDG-noUB_SM_BSZnf1"
FF_MODEL = "BSZ"
FF_NF = 1
NP_MODEL = "SM"
HQET_MODEL = "3/2/1"
UB_MODEL = "no-UB"
BELLE_FIT_OPTION = "sc"
```

HQET2/1/0 SM を実行する場合:

```python
ANALYSIS_NAME = "CasePDG-noUB_SM_HQET210"
FF_MODEL = "HQET"
FF_NF = 1
NP_MODEL = "SM"
HQET_MODEL = "2/1/0"
UB_MODEL = "no-UB"
BELLE_FIT_OPTION = "sc"
```

UB付きHQETを実行する場合だけ `UB_MODEL` を `"soft-UB"` または `"hard-UB"` に変える。BGL/BSZでは `effective_ub_model()` により UB は無効化される。

## 注意

- `Stan_FitAnalysis.py` を正規に使う場合、inline Python で `CmdStanModel(..., exe_file=..., compile=False)` を直接呼ぶ必要はない。
- コンパイル済み executable を明示的に再利用するテストランは特殊用途。通常解析では `python Stan_FitAnalysis.py` を実行する。
- `ANALYSIS_NAME` を使い回すと `Stan_outputs/<ANALYSIS_NAME>/` の中身が上書き・混在しやすいので、ケースごとに名前を変える。
- pre-fit timeoutで捨てられたchainは初期値候補に使われない。これは悪い初期条件で極端に遅くなるchainを避けるため。
- full MCMCは時間がかかる。短い動作確認だけしたい場合は `PRE_FIT_SETTINGS` と `FIT_SETTINGS` を一時的に小さくする。ただし、その結果を物理的なfit結果として扱わない。
