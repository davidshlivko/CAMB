"""Reproduce a Planck 2018 LCDM background and plot selected IDE histories.

Run from the CAMB-IDE worktree with:
    python -m analysis.planck_ide_background

Planck 2018 TT,TE,EE+lowE+lensing posterior means (Table 2 of Planck VI):
    H0 = 67.36 km/s/Mpc, omega_b = 0.02237, omega_c = 0.1200.
The base model fixes one 0.06 eV neutrino, N_eff = 3.046, and T_CMB = 2.7255 K.
https://doi.org/10.1051/0004-6361/201833910
"""

from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

import camb
from camb import constants

matplotlib.use("Agg")


PLANCK_H0 = 67.36
PLANCK_OMBH2 = 0.02237
PLANCK_OMCH2 = 0.1200
NEFF = 3.046
MNU = 0.06
TCMB = 2.7255
CASES = [
    (0.0, 0.0, "#333333", "Planck $\\Lambda$CDM limit"),
    (0.05, 0.30, "#0072B2", r"$\beta=+0.05,\ \lambda=0.30$"),
    (-0.05, 0.30, "#D55E00", r"$\beta=-0.05,\ \lambda=0.30$"),
    (0.10, 0.60, "#009E73", r"$\beta=+0.10,\ \lambda=0.60$"),
]


def planck_lcdm_background():
    pars = camb.CAMBparams()
    pars.set_cosmology(
        H0=PLANCK_H0,
        ombh2=PLANCK_OMBH2,
        omch2=PLANCK_OMCH2,
        mnu=MNU,
        nnu=NEFF,
        standard_neutrino_neff=NEFF,
        TCMB=TCMB,
    )
    return camb.get_background(pars, no_thermo=True)


def ide_background(v0, beta, lam):
    pars = camb.CAMBparams()
    pars.set_ide_cosmology(
        V0_over_rho100=v0,
        omc_early=PLANCK_OMCH2,
        beta=beta,
        lam=lam,
        ombh2=PLANCK_OMBH2,
        mnu=MNU,
        nnu=NEFF,
        standard_neutrino_neff=NEFF,
        TCMB=TCMB,
    )
    return camb.get_background(pars, no_thermo=True)


def main():
    outdir = Path(__file__).resolve().parent
    rho100_grho = 3 * (100 * 1000 / constants.c) ** 2
    standard = planck_lcdm_background()
    v0 = standard.grhov / rho100_grho
    a = np.geomspace(1e-8, 1.0, 650)

    fig, axes = plt.subplots(3, 1, figsize=(8.5, 10), sharex=True, constrained_layout=True)
    fig.suptitle("Exponential-mass IDE backgrounds at fixed Planck-calibrated inputs", fontsize=13)
    axes[0].set_title(
        rf"$V_0/\rho_{{\rm crit,100}}={v0:.9f}$, "
        rf"$\omega_{{c,\rm early}}={PLANCK_OMCH2:.4f}$; "
        r"$\omega_b=0.02237$, $m_\nu=0.06\,\mathrm{eV}$",
        fontsize=10,
        pad=9,
    )

    print(f"V0_over_rho100 = {v0:.12f}")
    print(f"omc_early = {PLANCK_OMCH2:.4f}")
    print("beta   lam       H0        omch2_0      phi_0")
    for beta, lam, color, label in CASES:
        bg = ide_background(v0, beta, lam)
        history = bg.get_ide_history(a)
        phi = history["phi"]
        rho_c = history["rho_c_over_rho100"]
        mass_ratio = rho_c * a**3 / PLANCK_OMCH2
        line_style = "--" if beta == lam == 0 else "-"
        line_width = 1.7 if beta == lam == 0 else 2.0
        axes[0].plot(a, phi, color=color, ls=line_style, lw=line_width, label=label)
        axes[1].plot(a, rho_c, color=color, ls=line_style, lw=line_width)
        axes[2].plot(a, mass_ratio, color=color, ls=line_style, lw=line_width)
        print(f"{beta:+.2f} {lam:.2f}  {bg.Params.H0:10.6f}  {bg.Params.omch2:10.7f}  {phi[-1]:+10.6f}")
        if beta == lam == 0:
            np.testing.assert_allclose(bg.Params.H0, PLANCK_H0, rtol=0, atol=1e-10)
            np.testing.assert_allclose(bg.Params.omch2, PLANCK_OMCH2, rtol=0, atol=1e-12)
            np.testing.assert_allclose(phi, 0, rtol=0, atol=1e-12)

    axes[0].axhline(0, color="0.75", lw=0.8)
    axes[0].set_ylabel(r"Dimensionless field $\varphi$")
    axes[0].legend(loc="upper left", frameon=False, fontsize=9)
    axes[1].set_yscale("log")
    axes[1].set_ylabel(r"$\rho_c/\rho_{\rm crit,100}$")
    axes[2].axhline(1, color="0.75", lw=0.8)
    axes[2].set_ylabel(r"$a^3\rho_c/(\omega_{c,\rm early}\rho_{\rm crit,100})$")
    axes[2].set_xlabel("Scale factor $a$")
    for ax in axes:
        ax.set_xscale("log")
        ax.set_xlim(a[0], a[-1])
        ax.grid(alpha=0.18, which="both")

    for extension in ("png", "pdf"):
        path = outdir / f"planck_ide_background.{extension}"
        fig.savefig(path, dpi=220, bbox_inches="tight")
        print(path)
    plt.close(fig)


if __name__ == "__main__":
    main()
