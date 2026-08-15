# Extended Stan fit analysis

## 1. 目的

この文書は、既存のフィットコードを保存したまま追加した拡張版

- `Stan_CodeGenerator_extended.py`
- `Stan_FitAnalysis_extended.py`

の構造、設計方針、実行方法、出力をまとめたものである。別スレッドで解析を継続するときは、まずこの文書を参照する。

拡張版の主目的は、次の統計解析を同じ枠組みで実行することである。

1. 全データフィット
2. Leave-one-dataset-out（LODO）fit
3. whitened residual／eigenmode report
4. posterior predictive test（PPC）
5. node-splitting analysis
6. covariance-shrinkage scan

元の `Stan_CodeGenerator.py` と `Stan_FitAnalysis.py` は編集していない。また、既存の補助モジュール

- `Stan_DataPoint.py`
- `Stan_FitFunction.py`
- `Stan_IntFunctions.py`

はそのまま再利用する。

## 2. 設計方針

### 2.1 一つの Stan model を全フィットで共用する

Stan data に、各データセットを likelihood に含めるかどうかを指定する

```stan
array[11] int<lower=0, upper=1> use_dataset;
```

を追加した。Stan の `model` block では、概念的に

```stan
if (use_dataset[i] == 1)
  target += multi_normal_cholesky_lpdf(...);
```

のように各 likelihood を ON/OFF する。

このため、全データフィット、各 LODO fit、各データセットだけを用いる fit は、同一の Stan source と同一の executable を使う。フィットごとに変更するのは standata 中の `use_dataset` だけである。

したがって、ここでいう「再コンパイルしない」とは「再フィットしない」という意味ではない。異なる data mask を使う各 case について、pre-fit と MCMC sampling はそれぞれ必要である。

### 2.2 covariance matrix は Cholesky 分解して使う

各データセットの covariance matrix は `transformed data` block で Cholesky 分解し、likelihood には `multi_normal_cholesky_lpdf` を使う。PPC の replicated data の生成にも同じ Cholesky factor を使う。

したがって、入力 covariance matrix は正定値でなければならない。少なくとも現在の `BELLE_FIT_OPTION="sc"` には rank deficiency の問題があるため、production run の前に covariance の定義・処理方針を確定する必要がある。コードの既定値は元解析との整合性のため `sc` のままにしている。コンパイルおよび短い動作試験では `rms` を runtime override として用いた。

### 2.3 generated quantities は必要な量に限定する

元コードの `generated quantities` にあった `chi_sq`、`chi_sq_total`、`RD_func`、Belle の raw output などは拡張版では出力しない。

既存解析との接続に必要な $|V_{cb}|$ 関連量

- `rawBrRatio_func`
- `BrRatio_func`
- `BrRatio_data`
- `Vcb_Br`
- `Vcb_mean`

は残している。

さらに各データセット $g$ について、次を出力する。

- `T_obs_<dataset>`: 観測データと予測の二次形式
- `T_rep_<dataset>`: replicated data と予測の二次形式
- `log_lik_<dataset>`: 観測データの multivariate Gaussian log likelihood

予測ベクトル自体は `transformed parameters` にあり、runner 側で保存・診断に利用する。

## 3. データセットと mask の順序

`DATASET_ORDER` と `use_dataset` の位置は次のように対応する。

| `use_dataset` の位置 | データセット |
|---:|---|
| 1 | `MILC15` |
| 2 | `MILC21` |
| 3 | `JLQCD23` |
| 4 | `HPQCD23` |
| 5 | `LCSR18` |
| 6 | `LCSR23` |
| 7 | `BelleFITD` |
| 8 | `BelleFITDstW` |
| 9 | `BelleFITDstCosL` |
| 10 | `BelleFITDstCosV` |
| 11 | `BelleFITDstChi` |

全データ fit では全要素を 1 にし、`LODO_LCSR23` では 6 番目だけを 0 にする。`Only_LCSR23` では 6 番目だけを 1 にする。

BrRatio2025 の情報はこの dataset mask とは別に常に model に含まれる。

## 4. ファイルごとの役割

### 4.1 `Stan_CodeGenerator_extended.py`

このファイルは、拡張版の Stan source、standata、parameter information を生成する。

主な変更点は次の通りである。

- 全データセットを常に Stan data として宣言する。
- `use_dataset` を standata に追加する。
- 各データセットの likelihood を `use_dataset` で条件分岐する。
- covariance matrix の Cholesky factor を生成する。
- dataset ごとの `T_obs`、`T_rep`、`log_lik` を生成する。
- Stan source の内容が既存ファイルと同じ場合は書き直さない。

最後の処理により、別の invocation でも source の timestamp を不必要に更新せず、CmdStan の既存 executable を再利用しやすくしている。

### 4.2 `Stan_FitAnalysis_extended.py`

このファイルは解析 case を組み立て、Stan model を一度だけ読み込み、各 case の standata を変えて順番に fit する。また、fit 後の診断量を計算して CSV にまとめる。

通常の entry point は `run_extended_analysis()` である。互換用の `run_fit()` は全データ case だけを実行する。

## 5. 解析を制御する設定

主要な設定は `Stan_FitAnalysis_extended.py` の先頭付近にある。

```python
RUN_ALL_DATA = True
RUN_LODO = False
LODO_DATASETS = list(DATASET_ORDER)

RUN_NODE_SPLIT = False
NODE_SPLIT_DATASETS = []

RUN_COVARIANCE_SCAN = False
COVARIANCE_SCAN_DATASETS = []
COVARIANCE_ALPHAS = [0.0, 0.25, 0.5, 0.75, 1.0]

PHYSICS_OUTPUTS = ["Vcb_mean"]
```

既定状態では全データ fit だけを実行する。全データセットについて LODO を行うには、

```python
RUN_ALL_DATA = True
RUN_LODO = True
LODO_DATASETS = list(DATASET_ORDER)
```

とする。この場合、全データ fit 1 回と LODO fit 11 回の合計 12 回の MCMC fit を実行するが、Stan compilation は 1 回だけである。

一部のデータセットだけを調べる場合は、例えば

```python
RUN_LODO = True
LODO_DATASETS = ["LCSR23", "HPQCD23"]
```

とする。

## 6. case の構築と自動依存関係

runner は設定から必要な case を自動的に組み立てる。

### 6.1 全データ fit

`AllData` case を作り、全 dataset likelihood を有効にする。

LODO、node split、covariance scan のいずれかを有効にすると、比較の基準として `AllData` も自動的に含まれる。

### 6.2 LODO fit

対象 $g$ ごとに `LODO_<g>` case を作る。対象 $g$ の likelihood だけを除き、それ以外を使って posterior を求める。

全データ posterior と LODO posterior の物理量の変化、および除外したデータに対する予測診断を調べる。

### 6.3 node-splitting analysis

対象 $g$ ごとに次の二つを作る。

- `LODO_<g>`: $g$ 以外のデータから得た posterior
- `Only_<g>`: $g$ だけから得た posterior

`RUN_NODE_SPLIT=True` のとき、対応する LODO case は `RUN_LODO=False` でも自動的に作られる。

両 posterior から、`NODE_SPLIT_SEPARATOR_PREFIXES` で指定した共通の物理予測ベクトルを計算し、その差を比較する。既定の prefix は

```python
NODE_SPLIT_SEPARATOR_PREFIXES = [
    "MILC15func",
    "MILC21func",
    "JLQCD23func",
    "HPQCD23func",
    "LCSR18func",
    "LCSR23func",
]
```

である。

### 6.4 covariance-shrinkage scan

対象データセットの covariance を

\[
C(\alpha)=D\{(1-\alpha)I+\alpha R\}D
\]

で置き換える。ここで $D$ は各成分の標準偏差、$R$ は元の correlation matrix である。

- $\alpha=1$: 元の covariance
- $\alpha=0$: 同じ分散を保った対角 covariance

各 $\alpha$ について全データ fit をやり直す。$\alpha=1$ は baseline の `AllData` を再利用する。

また、対象データを除いた `LODO_<g>` posterior は $\alpha$ に依存しないため 1 回だけ fit し、その posterior predictive distribution を各 $C(\alpha)$ で再評価する。したがって LODO fit を $\alpha$ ごとに繰り返す必要はない。

`RUN_COVARIANCE_SCAN=True` のときも、必要な `LODO_<g>` は自動的に作られる。

## 7. fit 後に計算する診断量

### 7.1 whitened residual／eigenmode report

各 posterior sample と各データセットについて、通常の residual、標準偏差で割った bin pull、covariance の固有基底で whiten した eigenmode pull を計算する。

出力には、各 mode の固有値、最大 loadings を持つ bin、全 loadings、各 mode の $\chi^2$ 寄与などを含める。これにより、単に全体の不一致を示すだけでなく、どの相関付き線形結合が不一致を作っているかを調べる。

### 7.2 posterior predictive test

Stan が出力した `T_obs` と `T_rep` から、

\[
p_{\mathrm{PPC}}=\Pr(T^{\mathrm{rep}}\geq T^{\mathrm{obs}}\mid y)
\]

を推定する。

さらに Python 側で replicated residual を生成し、最大 absolute bin pull と最大 absolute eigenmode pull に基づく tail probability も計算する。個別データセットだけでなく、全データセットを横断した最大値について `GLOBAL_MAX` 行も出力する。

### 7.3 LODO influence

`PHYSICS_OUTPUTS` で指定した量について、全データ posterior と LODO posterior の

- median
- median の差
- 全データ posterior standard deviation で規格化した差
- credible interval width
- width ratio

を比較する。既定では `Vcb_mean` を比較する。

### 7.4 node split

`Only_<g>` と `LODO_<g>` の共通予測ベクトルの差 $\Delta$ を調べる。

各成分について差、標準偏差、Z score、Gaussian approximation による両側 p-value を出力する。また、多変量比較として

\[
Q=\Delta^{\mathsf T}V_\Delta^{+}\Delta
\]

を計算し、effective rank と $\chi^2$ approximation による p-value を出力する。$V_\Delta^{+}$ は数値的 rank を考慮した擬似逆行列である。

### 7.5 covariance scan

各 $\alpha$ で、物理量、PPC、log predictive density、whitened residual／eigenmode diagnostics がどのように変化するかを出力する。

変化が $\alpha$ に強く依存する場合は、不一致の評価が off-diagonal correlation の仮定に敏感であることを示す。対角化しても不一致が残る場合は、中心値または各点の分散だけでも tension が存在する可能性が高い。

## 8. 出力ディレクトリとファイル

generator の出力は解析名の下に置かれる。

```text
Stan_outputs/<ANALYSIS_NAME>/
  standata_*.json
  stancode_*.stan
  stancode_*_param.json
  stancode_*                 # compiled executable
```

各 runner invocation は timestamp 付きディレクトリを作る。

```text
Stan_outputs/<ANALYSIS_NAME>/extended_runs_<MMDD_HHMMSS>/
```

その中に case ごとのディレクトリを作る。

```text
AllData/
LODO_<dataset>/
Only_<dataset>/
CovShrink_<dataset>_a.../
```

各 case directory には、おおむね次が入る。

- `cmdstan_prefit/`: pre-fit の CmdStan raw output
- `cmdstan_fit/`: main fit の CmdStan raw output
- `StanFit_<ANALYSIS_NAME>_<case>.csv`: posterior samples

run directory 直下には解析結果を集約した CSV を出力する。

| ファイル | 内容 |
|---|---|
| `dataset_diagnostics.csv` | dataset 単位の $T$, PPC, log predictive density, covariance condition など |
| `bin_diagnostics.csv` | 通常の bin pull の posterior summary |
| `eigenmode_diagnostics.csv` | eigenvalue、loadings、eigenmode pull、$\chi^2$ 寄与 |
| `lodo_influence.csv` | 全データ fit と LODO fit の物理量比較 |
| `node_split_components.csv` | node split の成分ごとの差、Z score、p-value |
| `node_split_summary.csv` | node split の多変量 $Q$、rank、p-value |
| `covariance_scan_outputs.csv` | 各 $\alpha$ の fit 結果と診断量 |
| `covariance_scan_lodo.csv` | 一つの LODO posterior を各 $C(\alpha)$ で再評価した結果 |
| `run_metadata.json` | 実行設定、case、出力場所などの記録 |

解析モードを有効にしていない場合、そのモード専用 CSV は作られない。

Stan の CSV は `sig_figs=18` で保存する。これは、保存した予測値から residual の二次形式を再構成するときの丸め誤差を抑えるためである。

## 9. 代表的な実行設定

### 9.1 全データ fit のみ

```python
RUN_ALL_DATA = True
RUN_LODO = False
RUN_NODE_SPLIT = False
RUN_COVARIANCE_SCAN = False
```

### 9.2 全 11 データセットの LODO

```python
RUN_ALL_DATA = True
RUN_LODO = True
LODO_DATASETS = list(DATASET_ORDER)
RUN_NODE_SPLIT = False
RUN_COVARIANCE_SCAN = False
```

### 9.3 特定データセットの node split

```python
RUN_LODO = False
RUN_NODE_SPLIT = True
NODE_SPLIT_DATASETS = ["LCSR23"]
RUN_COVARIANCE_SCAN = False
```

この場合は `AllData`、`LODO_LCSR23`、`Only_LCSR23` が自動的に実行される。

### 9.4 特定データセットの covariance scan

```python
RUN_LODO = False
RUN_NODE_SPLIT = False
RUN_COVARIANCE_SCAN = True
COVARIANCE_SCAN_DATASETS = ["LCSR23"]
COVARIANCE_ALPHAS = [0.0, 0.25, 0.5, 0.75, 1.0]
```

`AllData`、`LODO_LCSR23`、および $\alpha<1$ の covariance scan fit が実行される。$\alpha=1$ は `AllData` を再利用する。

### 9.5 複数解析を同時に有効化する場合

LODO、node split、covariance scan が要求する case に重複があれば、自動的に一つにまとめられる。例えば node split と covariance scan の双方が `LODO_LCSR23` を必要としても、その fit は一度しか行わない。

## 10. 実行方法

プロジェクトの指定された Stan Python environment を有効にし、コードディレクトリから実行する。

```bash
cd /Users/ryoutaro/Documents/Codex/Project_VcbFit/mycode_MetaFitAnalysis
source /Users/ryoutaro/stan_env/bin/activate
python Stan_FitAnalysis_extended.py
```

`Stan_FitAnalysis_extended.py` は最初に Wolfram scripts を通じて Stan inputs/model を生成し、その後 `cmdstanpy` を実行する。Wolfram kernel が sandbox 内で license や shared-memory の問題を起こす場合は、物理コードを変更せず Wolfram step を sandbox 外で実行する。

## 11. 検証済みの範囲

現時点で次を確認済みである。

- 二つの extended Python files の syntax check
- 生成した Stan source の compilation
- 全データ mask と LCSR23 除外 mask で Stan source が同一であること
- 一つの `CmdStanModel`／executable を複数 case で共用できること
- `AllData + LODO_LCSR23` の短い end-to-end smoke test
- `AllData + LODO_LCSR23 + Only_LCSR23 + CovShrink_LCSR23` の短い smoke test
- 各 summary CSV の生成

ただし smoke test は 1 chain、warmup 10、sampling 10 という実行確認専用の設定であり、多数の divergence が出ている。これらの数値は物理解析には使用できない。

## 12. 既知の注意点

1. `BELLE_FIT_OPTION="sc"` の covariance は現在 rank deficient であり、正定値を仮定する Cholesky likelihood と両立しない。production run の前に扱いを決める。
2. LODO、node split、covariance scan の case 数を増やすと、コンパイルは一度でも MCMC fit の回数と計算時間は増える。
3. node split の Gaussian p-value は posterior difference の Gaussian approximation に基づく。非 Gaussian 性が強い場合は posterior sample 自体の分布も確認する。
4. PPC の tail probability は posterior predictive check であり、frequentist hypothesis test の一様な p-value と同じ解釈をしない。
5. eigenmode の符号は固有ベクトルの任意性を持つため、符号そのものではなく絶対 pull、寄与、loadings の組合せを見る。
6. 本格実行では R-hat、effective sample size、divergence、treedepth など通常の MCMC diagnostics を先に確認する。

## 13. 別スレッドで作業を再開するときの確認順序

1. この文書と `Stan_FitAnalysis_extended.py` の解析制御変数を読む。
2. 対象データセットと必要な解析モードだけを設定する。
3. `BELLE_FIT_OPTION="sc"` の rank deficiency をどう扱うか確認する。
4. 必要なら `rms` と短い sampler settings で smoke test を行う。
5. 生成された `run_metadata.json` で、実際に作られた case と設定を確認する。
6. production settings で実行し、MCMC diagnostics を確認してから統計診断 CSV を解釈する。
