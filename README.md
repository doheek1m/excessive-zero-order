# Excessive-Zero Order

R code and reproducibility notes for the manuscript
**Counting models with excessive zeros ensuring stochastic monotonicity**.

This project studies claim-count models for insurance data with many zero
outcomes. The actuarial question is whether a model preserves the
**credibility order**: after observing a larger claim history, the posterior
predictive distribution for future claims should not move downward.

The manuscript shows that standard excessive-zero models with *separate*
random effects for the zero and positive-count components can violate this
monotonicity property. It then proposes **comonotonic random-effect hurdle
models** that preserve stochastic monotonicity while retaining competitive
out-of-sample performance.

The repository contains two independent bodies of evidence:

1. **Simulation study** (`R_simulation/`) — closed-form Monte Carlo
   experiments that demonstrate the violation, and its dependence on the
   variance and correlation parameters, without any estimation step.
2. **Empirical study** (`model/`, `Benchmark_model/`) — Bayesian fits to the
   Wisconsin LGPIF collision data, with an out-of-sample validation and a
   credibility-order diagnostic.

---

## Repository structure

```text
.
├── README.md
├── R_simulation/
│   ├── table1_2.R
│   ├── Supplementary_Table_S1.R
│   ├── Supplementary_Table_S2.R
│   └── Supplementary_Table_S3.R
├── model/
│   ├── Model2_MVN_Poisson_hurdle_nimble_f.R
│   ├── Model3_Independent_Poisson_hurdle_nimble_condition_f.R
│   ├── Model4_Comonotonic_Poisson_hurdle_nimble_f.R
│   └── Model5_Comonotonic_NB_hurdle_nimble_f.R
├── Benchmark_model/
│   ├── BM1_2_Poisson_GLMM.R
│   └── benchmark_models_234.R
└── Data/
    ├── data.RData
    ├── dataout.RData
    └── Data_for_monotonicity_f.R
```

---

## Code files

### Simulation (`R_simulation/`)

These scripts are self-contained. They require no data file, draw directly
from the assumed random-effect distribution, and report Monte Carlo standard
errors alongside every reported mean.

| Path | Produces | What it shows |
|---|---|---|
| `R_simulation/table1_2.R` | Tables 1 and 2 | Model 2 (bivariate normal RE, Poisson-hurdle with `exp` link). Table 1 decreases `sigma1^2` over `{5, 2, 1, 0.1}`; Table 2 increases `sigma2^2` over `{0.01, 0.1, 1, 2}`. Compares `E[Y2 | Y1 = 0]` against `E[Y2 | Y1 = 1]` and locates the parameter regions where the credibility order fails. |
| `R_simulation/Supplementary_Table_S1.R` | Table S.1 | Same model, but comparing `E[Y2 | Y1 = 1]` against `E[Y2 | Y1 = 2]` as `rho` moves over `{-0.8, -0.5, 0, 0.5, 0.8}`, with `sigma1^2 = sigma2^2 = 0.5`, `mu = (0, -2)`. Shows the violation is not confined to the zero-versus-positive comparison. |
| `R_simulation/Supplementary_Table_S2.R` | Table S.2 | Bivariate-RE analogue of Model 4 — the softplus activation of Remark 1 under condition (14). Isolates the effect of the *link function* from the effect of the *comonotonic coupling*. |
| `R_simulation/Supplementary_Table_S3.R` | Table S.3 | Bivariate-RE analogue of Model S.1, the Poisson **zero-inflated** mixture under condition (S.2), with `c = 0`, `d = 2`. Confirms that the phenomenon is not specific to the hurdle construction. |

Simulation sizes and seeds are fixed inside each script
(`N = 2e6` for Tables 1, 2, S.2, S.3; `S = 3e5` for Table S.1;
seeds `123` and `456`).

### Empirical models (`model/`)

| Path | Model | Random effects | Credibility order |
|---|---|---|---|
| `model/Model2_MVN_Poisson_hurdle_nimble_f.R` | Model 2 — bivariate normal RE Poisson-hurdle | Separate hurdle and count effects, `exp` link | **Not** guaranteed |
| `model/Model4_Comonotonic_Poisson_hurdle_nimble_f.R` | Model 4 — comonotonic RE Poisson-hurdle | One shared latent effect, softplus link | Guaranteed |
| `model/Model5_Comonotonic_NB_hurdle_nimble_f.R` | Model 5 — comonotonic RE Negative-Binomial-hurdle | One shared latent effect | Guaranteed |


### Benchmarks (`Benchmark_model/`)

| Path | Model | Fitted by |
|---|---|---|
| `Benchmark_model/BM1_2_Poisson_GLMM.R` | BM 1 — Poisson GLMM (one random effect) | Bayesian MCMC (`nimble`) |
| `Benchmark_model/benchmark_models_234.R` | BM 2 — Poisson GLM, BM 3 — Poisson hurdle, BM 4 — Poisson zero-inflated | Maximum likelihood (`glm`, `glmmTMB`) |

### Data (`Data/`)

| Path | Contents |
|---|---|
| `Data/data.RData` | LGPIF training records, object `data` (2006–2010) |
| `Data/dataout.RData` | LGPIF validation records, object `dataout` (2011) |
| `Data/Data_for_monotonicity_f.R` | Builds every object the model scripts consume: `Y_mat`, `XX_train`, `wY_mat`, `ID_train`, the hurdle split `I_mat` / `N_mat` with their exposure masks, the test-side objects `Y_test`, `XX_test`, `ID_test`, and the counterfactual histories `*_c_0` / `*_c_1` used for the diagnostic. |

---

## Background

Insurance claim-frequency data often contain a high proportion of zeros.
Standard Poisson and Negative Binomial models may fit such data poorly
because they do not separately model the zero-generating mechanism.
Zero-inflated and hurdle models address this by adding a zero component, but
random-effect versions of these models can create undesirable posterior
credibility behaviour.

The paper formalises the following principle as stochastic monotonicity, or
the credibility order:

```text
A larger observed claim history should not imply a lower future claim-risk prediction.
```

Violating this order is a problem for experience rating, because an
additional claim could reduce a future premium-related predictive functional.

---

## Empirical data

The empirical study uses claim data from the Wisconsin Local Government
Property Insurance Fund (LGPIF), restricted to collision coverage for new and
old vehicles (`FreqCN + FreqCO`, `CoverageCN + CoverageCO`), keeping only
policies with positive coverage on both.

- **Training period:** 2006–2010 (five yearly observations per policy, unbalanced)
- **Validation period:** 2011
- **Response:** claim frequency `n`
- **Covariates:** entity type (city / county / school / town / village, one-hot)
  and a three-level coverage index formed from the empirical tertiles of
  `col.Cov`
- **Validation target:** one-step-ahead predictive mean for 2011
- **Panel size:** 409 policyholders in training; the test set is the subset of
  2011 records whose policy number appears in training

`Data/data.RData` and `Data/dataout.RData` are included in this repository.
Please cite the LGPIF source when reusing them.

---

## How to run

All commands assume the repository root as the working directory.

### 1. Simulation (no data needed, runs in minutes)

```bash
Rscript R_simulation/table1_2.R
Rscript R_simulation/Supplementary_Table_S1.R
Rscript R_simulation/Supplementary_Table_S2.R
Rscript R_simulation/Supplementary_Table_S3.R
```

### 2. Empirical models

Each model script begins by sourcing the data-preparation script, so run them
from the repository root:

```bash
Rscript Benchmark_model/BM1_2_Poisson_GLMM.R
Rscript Benchmark_model/benchmark_models_234.R

Rscript model/Model2_MVN_Poisson_hurdle_nimble_f.R
Rscript model/Model4_Comonotonic_Poisson_hurdle_nimble_f.R
Rscript model/Model5_Comonotonic_NB_hurdle_nimble_f.R
```

Each MCMC run uses `niter = 30000`, `nburnin = 10000`, `thin = 10`,
`nchains = 3` (the counterfactual refits inside `Model2_...R` use
`niter = 50000`). Expect a wall-clock time of roughly ten minutes per model
on a modern laptop; `Model2_...R` is longer because it fits three chains three
times over.

---

## Credibility-order diagnostic

For each policyholder the diagnostic holds `y_{i,1:4}` fixed, replaces the
most recent observation by two counterfactual values, refits, and compares the
resulting predictive functionals for 2011:

```text
Y_{i,5} = 0  versus  Y_{i,5} = 1
```

| Functional | Diagnostic inequality |
|---|---|
| Full coverage | `E[Y_{i,6} | history, Y_{i,5}=0] <= E[Y_{i,6} | history, Y_{i,5}=1]` |
| Deductible coverage | `E[(Y_{i,6}-d)_+ | history, Y_{i,5}=0] <= E[(Y_{i,6}-d)_+ | history, Y_{i,5}=1]`, for `d = 1, 2` |
| Limited coverage | `E[min(Y_{i,6}, d) | history, Y_{i,5}=0] <= E[min(Y_{i,6}, d) | history, Y_{i,5}=1]`, for `d = 1, 2` |


---

## Representative results

### Credibility-order violation rates (LGPIF)

| Model | Base `Y` | `(Y-d)+`, `d=1` | `(Y-d)+`, `d=2` | `min(Y,d)`, `d=1` | `min(Y,d)`, `d=2` |
|---|---:|---:|---:|---:|---:|
| Bivariate RE Poisson-hurdle (Model 2) | 9.29% | 10.76% | 12.72% | 0.00% | 0.98% |
| Comonotonic RE Poisson-hurdle (Model 4) | 0.00% | 0.00% | 0.00% | 0.00% | 0.00% |
| Comonotonic RE NB-hurdle (Model 5) | 0.00% | 0.00% | 0.00% | 0.00% | 0.00% |
| Poisson GLMM (BM 1) | 0.00% | 0.00% | 0.00% | 0.00% | 0.00% |


### Out-of-sample validation (2011)

Bayesian models are scored by posterior predictive log-density (LPD, summed
over test policies); the maximum-likelihood benchmarks are scored by plug-in
test log-likelihood. The two are **not** on the same scale and should be
compared within, not across, the two blocks.

| Model | Test log score | MSE | MAE | Guaranteed credibility order |
|---|---:|---:|---:|---|
| *Bayesian (posterior predictive LPD)* | | | | |
| BM 1 — Poisson GLMM | -208.6458 | 1.1171 | 0.6063 | Yes |
| Model 2 — bivariate RE Poisson-hurdle | -208.9472 | 1.0489 | 0.5990 | No |
| Model 4 — comonotonic RE Poisson-hurdle | -208.5055 | 0.9382 | 0.5994 | Yes |
| Model 5 — comonotonic RE NB-hurdle | -214.7405 | 0.9624 | 0.6060 | Yes |
| *Maximum likelihood (plug-in log-likelihood)* | | | | |
| BM 2 — Poisson GLM | -266.5104 | 2.0241 | 0.8575 | Trivial / no RE |
| BM 3 — Poisson hurdle | -301.3224 | 2.4239 | 0.9158 | Trivial / no RE |
| BM 4 — Poisson zero-inflated | -256.4331 | 1.9762 | 0.8396 | Trivial / no RE |


- Model 4 attains the best LPD and the lowest MSE among the random-effect
  models while guaranteeing the credibility order.
- Model 5 also guarantees the order, with MSE close to Model 4 and a weaker
  LPD.
- The differences in LPD among BM 1, Model 2 and Model 4 are small
  (under one nat over 253 test policies); the MSE gap is the more substantial
  comparison.
- The fixed-effect benchmarks are clearly worse on both MSE and MAE..
