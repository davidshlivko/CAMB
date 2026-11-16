# IDE scalar-mode equation residuals

Run `python analysis/ide_mode_residuals.py --max-residual 1e-4` after building CAMB-IDE. The script evaluates $k=0.005\,\mathrm{Mpc}^{-1}$ at conformal times 500, 2000, and 8000 Mpc, with $V_0/\rho_{100}=0.31$, $\omega_{c,\mathrm{early}}=0.12$, $\lambda=0.5$, and $\beta\in\{0,\pm0.05,\pm0.1\}$. It samples each integrated mode at $\tau\pm10^{-3}\tau$ and takes central differences. The background $\bar\varphi'$ is obtained independently by finite-differencing `get_ide_history` in $\ln a$.

All quantities below use CAMB's synchronous gauge. Let $C=8\pi Ga^2\bar\rho_c$, $U=8\pi Ga^2V$, $\pi=\delta\varphi'$, and $Q=\beta C(\bar\varphi'\delta_c+\pi)$. The script checks the CDM continuity and Euler equations, scalar Klein–Gordon equation, and the Einstein momentum equation $\eta_k'=\delta G_q/2$. For energy exchange it forms the scalar density and pressure perturbations $S=\bar\varphi'\pi-\lambda U\delta\varphi$ and $P=\bar\varphi'\pi+\lambda U\delta\varphi$, then checks

\[
(C\delta_c)' + \mathcal H C\delta_c + Ck(z+v_c)=-Q,
\]
\[
S'+\mathcal H S+3\mathcal H P+kz\bar\varphi'^2+k^2\bar\varphi'\delta\varphi=Q.
\]

Adding these gives a direct finite-difference check that the CDM and scalar exchange terms cancel. Each residual is divided by the sum of the absolute terms in its equation.

| $\beta$ | Largest residual, all equations | Einstein momentum | Scalar KG | Combined energy |
|---:|---:|---:|---:|---:|
| 0 | 8.93e-6 | 1.17e-6 | 7.98e-7 | 8.68e-6 |
| -0.05 | 5.74e-6 | 3.19e-6 | 4.99e-7 | 5.73e-6 |
| +0.05 | 5.74e-6 | 2.02e-6 | 4.33e-7 | 5.73e-6 |
| -0.10 | 6.55e-6 | 5.72e-6 | 5.11e-7 | 6.43e-6 |
| +0.10 | 5.66e-6 | 3.98e-6 | 4.70e-7 | 5.55e-6 |

Changing the finite-difference half-width from $10^{-3}\tau$ to $5\times10^{-4}\tau$ raises the largest residual to about $2.3\times10^{-5}$; increasing it to $2\times10^{-3}\tau$ gives about $2.3\times10^{-5}$. The smaller step amplifies quantization in CAMB's single-precision `delta_cdm` transfer output, while the larger step increases finite-difference truncation. The internal $z$ and shear identities are algebraic in CAMB and are therefore not counted as independent Einstein checks. These residuals establish consistency of the implemented evolution equations and do not establish CMB likelihood convergence or physical validity outside the tested modes.
