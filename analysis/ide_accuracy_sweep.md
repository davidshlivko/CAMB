# IDE numerical convergence beyond the pilot

The reproducible spectrum sweep is `analysis/ide_accuracy_sweep.py`. It uses
the same fixed physical and primordial inputs as the Planck pilot, with
`V0_over_rho100=0.31`, `omc_early=0.12`, linear lensing, and
`max_eta_k=18000`. Its baseline has `AccuracyBoost=lSampleBoost=lAccuracyBoost=2`,
`a_init=1e-10`, `history_points=2049`, `integrate_tol=1e-10`, and a 0.06 eV
massive neutrino. Spectrum differences are measured over `30 <= ell <= 1000`;
TE is scaled by `sqrt(TT*EE)`. Run one point per process with, for example,
`OMP_NUM_THREADS=1 python analysis/ide_accuracy_sweep.py --point-index 1 --output external/ide_accuracy_sweep_1.json`.
Point indices 0 through 4 correspond to `(beta, lam)` values `(0,0)`,
`(+0.05,0.5)`, `(-0.05,0.5)`, `(+0.1,1)`, and `(-0.1,1)`. Output files under
`external/` are ignored by Git; the script writes partial JSON after each
variant to preserve progress if a later calculation fails.

The sweep changes one control at a time: `a_init` to `1e-12`, history samples
to 4097, background tolerance to `1e-11`, `neutrino_q_boost` to 2,
`MassiveNuMethod` to `Nu_int`, `TimeSwitchBoost` to 2, `IntTolBoost` to 2, and
the three main accuracy boosts together from 2 to 3. `TimeSwitchBoost`
controls perturbation switch timing, including massive-neutrino switches. It
does not directly control the photon tight-coupling switch.

CAMB hard-codes the tight-coupling switch coefficient `ep0=0.01` in
`fortran/equations.f90`. To test it independently, we copied the current
CAMB-IDE source to `/private/tmp/ide-tc-copy-20261115`, changed only `ep0` to
0.005 in that copy, and built it with an isolated `CAMB_BUILD_DIR`. For
`beta=+0.05`, `lam=0.5`, `ell_max=2500`, and the baseline boosts of 2,
the maximum absolute scaled differences over `30 <= ell <= 2500` were:

| Spectrum | Maximum | 95th percentile |
|---|---:|---:|
| Lensed TT | 1.0061e-4 | 7.43e-5 |
| Lensed EE | 1.1469e-4 | 8.20e-5 |
| Lensed TE / sqrt(TT EE) | 4.0101e-5 | 2.94e-5 |
| Lensing potential | 1.2577e-6 | 1.20e-6 |

This tests one representative coupled point; it does not replace a full
range photon-hierarchy convergence study. The main CAMB-IDE source was not
changed for this experiment.

## Spectrum sweep results

The serial `beta=+0.05`, `lam=0.5` point completed at `ell_max=1000` in about
105 seconds. Its baseline resolved `H0=59.72736` km/s/Mpc and
`omch2=0.11559814`. The table gives maximum absolute scaled differences
against that baseline over `30 <= ell <= 1000`:

| Changed control | Lensed TT | Lensed EE | TE / sqrt(TT EE) | Lensing potential |
|---|---:|---:|---:|---:|
| `a_init`: `1e-10` → `1e-12` | 1.86e-7 | 5.23e-7 | 2.42e-7 | 1.15e-5 |
| `history_points`: 2049 → 4097 | 7.63e-11 | 6.10e-13 | 7.39e-11 | 6.45e-13 |
| `integrate_tol`: `1e-10` → `1e-11` | 0 | 0 | 0 | 0 |
| `neutrino_q_boost`: 1 → 2 | 6.86e-6 | 8.02e-6 | 3.74e-6 | 1.25e-5 |
| `MassiveNuMethod`: `Nu_best` → `Nu_int` | 1.87e-12 | 3.42e-12 | 1.35e-12 | 0 |
| `TimeSwitchBoost`: 1 → 2 | 1.85e-8 | 3.39e-9 | 1.21e-8 | 8.57e-8 |
| `IntTolBoost`: 1 → 2 | 1.09e-5 | 7.45e-6 | 3.95e-6 | 3.65e-6 |
| All three main boosts: 2 → 3 | 7.25e-4 | 2.16e-4 | 2.78e-4 | 5.72e-4 |

The first 45-calculation concurrent batch was interrupted when high-boost
Planck runs exhausted available memory. This one-point run was completed
serially and saved as `external/ide_accuracy_sweep_1.json`. The main boost
change affects many CAMB controls, so it is a convergence comparison rather
than an isolated perturbation integrator test. The remaining four points
are not yet measured; no full-range claim follows from this table.

## Prior-edge Planck log likelihood

At `beta=+0.1`, `lam=1`, all three main boosts at 4 gave
`logL=-52471.5131449` in 132 s. The first boost-5 attempt hit CAMB's
`max_l_evolve=256` hard ceiling. The main source limit was raised to 512 and
built. A serial boost-5 rerun completed in 246.5 s with
`logL=-52470.5473846`, a change of `+0.9657603`. The full-prior absolute
`|Delta logL| < 0.1` gate therefore still fails at this point;
the 12-step sampler in `ide_planck_pilot_results.md` is a pipeline check only.

Holding `AccuracyBoost=lAccuracyBoost=4` and raising only `lSampleBoost`
from 4 to 6 gave `logL=-52471.5110227` in 136.7 s, just `+0.0021223`
relative to all boosts at 4. The remaining all-boost 4→5 shift is therefore
not explained by high-ell interpolation density alone.

Further one-control changes at the same point and on the same likelihood were:

| `AccuracyBoost` | `lSampleBoost` | `lAccuracyBoost` | Total logL | Change from `(4,4,4)` | Runtime |
|---:|---:|---:|---:|---:|---:|
| 4 | 4 | 4 | -52471.5131449 | 0 | 132.0 s |
| 4 | 6 | 4 | -52471.5110227 | +0.0021223 | 136.7 s |
| 5 | 4 | 4 | -52470.6120211 | +0.9011238 | 214.7 s |
| 4 | 4 | 5 | -52471.5380815 | -0.0249366 | 141.3 s |
| 5 | 5 | 5 | -52470.5473846 | +0.9657603 | 246.5 s |

The dominant remaining drift is associated with `AccuracyBoost`, which
changes several integration tolerances, time steps, source sampling, and
switches. It has not yet been reduced to one sub-control. These large
absolute log-likelihood changes occur at `logL` around `-5.25e4`, far from
the pilot's good-likelihood region, but the predefined absolute gate still
fails.

To isolate likely sub-controls, `analysis/ide_accuracy_likelihood.py` keeps
`AccuracyBoost=lSampleBoost=lAccuracyBoost=4` at the same `beta=+0.1`,
`lam=1` point and changes only one nested CAMB accuracy field per run. Its
calls reuse the pilot's fixed cosmology, Planck likelihood data, and
`max_eta_k=18000`. The isolated checks set `SourcekAccuracyBoost`,
`TimeStepBoost`, or `IntkAccuracyBoost` to 1.25, compared to their defaults
of 1. The last control multiplies the k-sample density for high-ell C_l
integration.

| Additional nested control | High-ell logL | LowE logL | LowT logL | Total logL | Change from default | Runtime |
|---|---:|---:|---:|---:|---:|---:|
| None (all main boosts 4) | -52262.3666198 | -198.4425013 | -10.7040239 | -52471.5131449 | 0 | 132.0 s |
| `SourcekAccuracyBoost=1.25` | -52262.4096673 | -198.4425013 | -10.7041662 | -52471.5563348 | -0.0431898 | 153.6 s |
| `TimeStepBoost=1.25` | -52262.3700322 | -198.4425013 | -10.7040595 | -52471.5165929 | -0.0034480 | 138.2 s |
| `IntkAccuracyBoost=1.25` | -52262.3834386 | -198.4425013 | -10.7040845 | -52471.5300244 | -0.0168795 | 132.7 s |

The three runs were serial with `OMP_NUM_THREADS=1`. Their full settings and
results are in `external/planck-pilot/edge_pos_sourcek_1p25.json`,
`external/planck-pilot/edge_pos_timestep_1p25.json`, and
`external/planck-pilot/edge_pos_intk_1p25.json`, which are ignored by Git.
None of these one-control changes accounts for the `+0.9011238` shift from
raising `AccuracyBoost` alone from 4 to 5. Other sub-controls or their
combined effects remain unresolved.

## Follow-up at beta = +0.07, lambda = 1

With the pilot's fixed `V0_over_rho100=0.31`, `omc_early=0.12`, and other
physical inputs, the resolved `H0` is 50.081557 km/s/Mpc. Using the same
Planck likelihood and changing all three main CAMB boosts together gave:

| Boosts | Total logL | Change from preceding setting | Runtime |
|---:|---:|---:|---:|
| 3 | -27568.4566088 | — | 40.1 s |
| 4 | -27567.6590456 | +0.7975632 | 130.3 s |
| 5 | -27567.0980048 | +0.5610408 | 252.9 s |

Both changes exceed the required absolute `|Delta logL| < 0.1`. The shifts
come almost entirely from the high-ell Plik component; lowE is unchanged
and lowT changes by less than 0.001. The runs were serial with
`OMP_NUM_THREADS=1`; full output is in the ignored
`external/planck-pilot/beta_007_lam1_boost{3,4,5}.json` files. This point
does **not** pass the pilot accuracy gate at boost 3 or 4.

### Retuned inputs at the same beta and lambda

A two-parameter Nelder-Mead search using the same Planck likelihood at
default CAMB accuracy found a better-fitting point at
`V0_over_rho100=0.9672031528954944` and
`omc_early=0.12113151435361558`. All other fixed pilot inputs, including
`beta=+0.07` and `lam=1`, were retained. This is a numerical test point,
not a posterior fit. It resolves to `H0=67.27603` km/s/Mpc and
`omch2=0.11130512`.

| Boosts | Total logL | Change from preceding setting | Runtime |
|---:|---:|---:|---:|
| 3 | -505.5138033 | — | 41.8 s |
| 4 | -505.4731714 | +0.0406318 | 141.9 s |
| 5 | -505.4808917 | -0.0077203 | 269.4 s |

Both changes satisfy `|Delta logL| < 0.1` at this retuned point. The runs
were serial with `OMP_NUM_THREADS=1`; full outputs are in the ignored
`external/planck-pilot/beta_007_lam1_fit_boost{3,4,5}.json` files. The
retuned potential amplitude is **outside** the pilot's fixed prior
`V0_over_rho100 in [0.20, 0.42]`. Thus this result demonstrates that the
accuracy failure at the original fixed-input point is not intrinsic to
`beta=0.07`, but it does not validate the original four-input prior or
establish a posterior constraint.
