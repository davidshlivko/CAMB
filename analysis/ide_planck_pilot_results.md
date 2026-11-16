# Planck 2018 IDE CMB pipeline pilot

This is a software validation record, not an IDE posterior constraint. The numerical
stability gate over the full prior range has **not** passed.

## Setup

- CAMB-IDE worktree, based on commit `5372c06180cf79fe085fe201f07a2431dcef18ac`;
  clean reference CAMB is a detached worktree at that commit.
- Python 3.11 and pinned Cobaya 3.6.2. The environment and likelihood data are under
  ignored `external/`. Cobaya downloaded the three native datasets from
  `CobayaSampler/planck_native_data` release `v1`.
- Likelihoods: `planck_2018_highl_plik.TTTEEE_lite_native`,
  `planck_2018_lowl.TT`, and `planck_2018_lowl.EE`.
- Fixed inputs: `ombh2=0.02237`, `mnu=0.06` eV, `nnu=3.046`,
  `As=2.10e-9`, `ns=0.9649`, `tau=0.0544`, `A_planck=1`.
  We verified the resolved CAMB parameters retain the first four physical density
  and reionization inputs at two distinct Cobaya sampled points.
- Linear lensed CMB only: `NonLinear_none`, `DoLateRadTruncation=False`,
  `lens_potential_accuracy=0`, explicit `max_eta_k=18000`. The three CAMB accuracy
  boosts in the table below are `AccuracyBoost`, `lSampleBoost`, and
  `lAccuracyBoost`, changed together unless stated otherwise.
- Likelihood evaluations and the pilot are reproducible with
  `analysis/ide_planck_pilot.py`; evaluation JSON and sampler files are under
  ignored `external/planck-pilot/`.

To recreate the isolated environment and likelihood data from this worktree:

```text
python3.11 -m venv external/planck-pilot-venv
external/planck-pilot-venv/bin/python -m pip install cobaya==3.6.2
external/planck-pilot-venv/bin/cobaya-install planck_2018_highl_plik.TTTEEE_lite_native planck_2018_lowl.TT planck_2018_lowl.EE --packages-path external/planck-packages --just-data --no-set-global
```

After building the CAMB extension, run a fixed-point likelihood or the pilot with:

```text
external/planck-pilot-venv/bin/python analysis/ide_planck_pilot.py evaluate --kind ide
external/planck-pilot-venv/bin/python analysis/ide_planck_pilot.py pilot --samples 12 --accuracy 3
```

## Same-likelihood baseline check

At `V0_over_rho100=0.31`, `omc_early=0.12`, and `beta=lam=0`, IDE resolves
`H0=67.30906608490525` km/s/Mpc and `omch2=0.12`. Clean CAMB used exactly those
physical inputs and the same three likelihood files. At default boosts, the total
log likelihoods are `-503.0054855883` (IDE) and `-503.0047646876` (clean),
giving `IDE-clean=-0.0007209007`. A repeated IDE evaluation returned identical
likelihood values.

| Point `(beta, lam)` | Boosts | Log likelihood change when tightened |
| --- | --- | ---: |
| IDE `(0, 0)` | 1 → 1.5 | -0.115754 |
| Clean ΛCDM `(0, 0)` | 1 → 1.5 | -0.118185 |
| IDE `(0, 0)` | 1.5 → 2 | -0.043690 |
| Clean ΛCDM `(0, 0)` | 1.5 → 2 | -0.041600 |
| IDE `(+0.01, 0.1)` | 1.5 → 2 | -0.094587 |
| IDE `(-0.01, 0.1)` | 1.5 → 2 | -0.026836 |
| IDE `(+0.01, 0.1)` | 3 → 4 | +0.056026 |
| IDE `(+0.05, 0.5)` | 2 → 3 | +1.121059 |
| IDE `(+0.05, 0.5)` | 3 → 4 | +0.046387 |
| IDE `(-0.05, 0.5)` | 2 → 3 | -0.258605 |
| IDE `(-0.05, 0.5)` | 3 → 4 | +0.099166 |
| IDE `(+0.1, 1)` | 3 → 4 | +2.108547 |
| IDE `(-0.1, 1)` | 3 → 4 | +0.623221 |

The default-accuracy shift at the ΛCDM limit is shared by clean CAMB. The
stronger-coupling points also have large derived changes in `H0` at fixed
`V0_over_rho100` and `omc_early`: about 59.7/73.2 km/s/Mpc for `beta=±0.05,
lam=0.5`, and 45.8/90.2 km/s/Mpc for `beta=±0.1, lam=1`. They remain finite
but are far from this fixed-input Planck likelihood's preferred region. At
`beta=+0.05, lam=0.5`, changing just `AccuracyBoost`, `lSampleBoost`, or
`lAccuracyBoost` from 2 to 3 moves log likelihood by `+0.3943`, `+0.5467`, or
`+0.3976`, respectively. No single boost accounts for the instability.

## Four-input pilot

Cobaya's MCMC sampler varied `V0_over_rho100 ∈ [0.20, 0.42]`,
`omc_early ∈ [0.10, 0.14]`, `beta ∈ [-0.1, 0.1]`, and `lam ∈ [0, 1]`.
It used the fixed inputs and three CAMB boosts of 3 above, seed `20261115`,
and started at `(0.31, 0.12, 0, 0.1)`.

The run completed **12 accepted steps** with acceptance rate **0.75** in
**131.7 seconds**. It raised **zero CAMB/likelihood computation errors**;
`stop_at_error=True` was active. The chain and full evaluation JSON remain under
ignored `external/planck-pilot/`. The chain is deliberately too short and the
full-range stability gate is open, so it cannot be used to infer an IDE posterior.
