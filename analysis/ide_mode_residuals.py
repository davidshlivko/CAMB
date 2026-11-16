"""Finite-difference checks of the IDE scalar and CDM perturbation equations.

The variables returned by ``get_time_evolution`` are synchronous-gauge CAMB
variables. Each residual differentiates independently sampled mode solutions;
the checks therefore test the integrated solution, not only its RHS routine.
Run after building the CAMB-IDE Fortran library.
"""

import argparse
import importlib
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
camb = importlib.import_module("camb")
constants = importlib.import_module("camb.constants")

NAMES = [
    "a",
    "H",
    "etak",
    "delta_cdm",
    "delta_phi_ide",
    "delta_phi_prime_ide",
    "v_cdm_synch_ide",
    "z_synch_ide",
    "dgq_ide",
    "dgrho_ide",
    "sigma_synch_ide",
]
GRHO100 = 3 * (100 * 1000 / constants.c) ** 2


def normalized_residual(terms):
    """Residual divided by the sum of absolute terms, avoiding zero division."""
    scale = sum(abs(term) for term in terms)
    return abs(sum(terms)) / max(scale, 1e-100)


def run_point(beta, lam, times, k, step_fraction):
    pars = camb.set_params(
        dark_energy_model="IDEExponential",
        V0_over_rho100=0.31,
        omc_early=0.12,
        beta=beta,
        lam=lam,
        ombh2=0.02237,
        mnu=0.06,
        tau=0.0544,
        As=2.10e-9,
        ns=0.9649,
        lmax=80,
        lens_potential_accuracy=0,
    )
    result = camb.get_results(pars)
    halfsteps = step_fraction * times
    eta = np.column_stack((times - halfsteps, times, times + halfsteps)).ravel()
    mode_array = result.get_time_evolution(k, eta, vars=NAMES, lAccuracyBoost=1)
    mode = {name: mode_array[:, i] for i, name in enumerate(NAMES)}
    a = mode["a"]
    phi = result.get_ide_history(a)["phi"]
    log_a_step = 1e-4
    a_bracket = np.column_stack((a * np.exp(-log_a_step), a * np.exp(log_a_step)))
    phi_bracket = result.get_ide_history(a_bracket.ravel())["phi"].reshape(-1, 2)
    phi_dot = (phi_bracket[:, 1] - phi_bracket[:, 0]) / (2 * log_a_step) * mode["H"]
    cdm_density = result.get_background_densities(a, vars=["cdm"])["cdm"] / a**2
    potential = GRHO100 * pars.DarkEnergy.V0_over_rho100 * a**2 * np.exp(-lam * phi)

    def dtime(values, i):
        return (values[3 * i + 2] - values[3 * i]) / (2 * halfsteps[i])

    rows = []
    for i, t in enumerate(times):
        j = 3 * i + 1
        h, delta_c = mode["H"][j], mode["delta_cdm"][j]
        delta_phi, pi = mode["delta_phi_ide"][j], mode["delta_phi_prime_ide"][j]
        vc, z = mode["v_cdm_synch_ide"][j], mode["z_synch_ide"][j]
        c, u, phi_prime = cdm_density[j], potential[j], phi_dot[j]
        scalar_drho = phi_dot * mode["delta_phi_prime_ide"] - lam * potential * mode["delta_phi_ide"]
        scalar_dp = phi_prime * pi + lam * u * delta_phi
        cdm_drho = cdm_density * mode["delta_cdm"]
        exchange = beta * c * (phi_prime * delta_c + pi)

        residuals = {
            "cdm_continuity": normalized_residual([dtime(mode["delta_cdm"], i), k * (z + vc), beta * pi]),
            "cdm_euler": normalized_residual(
                [dtime(mode["v_cdm_synch_ide"], i), (h - beta * phi_prime) * vc, beta * k * delta_phi]
            ),
            "scalar_kg": normalized_residual(
                [
                    dtime(mode["delta_phi_prime_ide"], i),
                    2 * h * pi,
                    (k**2 + lam**2 * u) * delta_phi,
                    k * z * phi_prime,
                    -beta * c * delta_c,
                ]
            ),
            "einstein_momentum": normalized_residual([dtime(mode["etak"], i), -0.5 * mode["dgq_ide"][j]]),
            "cdm_energy_exchange": normalized_residual(
                [dtime(cdm_drho, i), h * c * delta_c, c * k * (z + vc), exchange]
            ),
            "scalar_energy_exchange": normalized_residual(
                [
                    dtime(scalar_drho, i),
                    h * scalar_drho[j],
                    3 * h * scalar_dp,
                    k * z * phi_prime**2,
                    k**2 * phi_prime * delta_phi,
                    -exchange,
                ]
            ),
            "total_energy_exchange": normalized_residual(
                [
                    dtime(cdm_drho, i) + dtime(scalar_drho, i),
                    h * (c * delta_c + scalar_drho[j]),
                    3 * h * scalar_dp,
                    c * k * (z + vc),
                    k * z * phi_prime**2,
                    k**2 * phi_prime * delta_phi,
                ]
            ),
        }
        rows.append(
            {
                "eta_mpc": float(t),
                "a": float(a[j]),
                "exchange_source": float(exchange),
                **{key: float(value) for key, value in residuals.items()},
            }
        )

    keys = list(residuals)
    return {
        "beta": beta,
        "lam": lam,
        "k_mpc_inv": k,
        "finite_difference_fraction": step_fraction,
        "rows": rows,
        "max_normalized_residual": {key: max(row[key] for row in rows) for key in keys},
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--beta", type=float, nargs="+", default=[0.0, -0.05, 0.05, -0.1, 0.1])
    parser.add_argument("--lam", type=float, default=0.5)
    parser.add_argument("--times", type=float, nargs="+", default=[500.0, 2000.0, 8000.0])
    parser.add_argument("--k", type=float, default=0.005)
    parser.add_argument("--step-fraction", type=float, default=1e-3)
    parser.add_argument("--max-residual", type=float, help="Fail if any normalized residual exceeds this value")
    args = parser.parse_args()

    times = np.asarray(args.times, dtype=np.float64)
    if np.any(times <= 0) or np.any(np.diff(times) <= 0) or args.k <= 0 or args.step_fraction <= 0:
        parser.error("times must be positive and increasing; k and step fraction must be positive")
    results = [run_point(beta, args.lam, times, args.k, args.step_fraction) for beta in args.beta]
    print(json.dumps(results, indent=2, sort_keys=True))
    if args.max_residual is not None:
        worst = max(max(result["max_normalized_residual"].values()) for result in results)
        if worst > args.max_residual:
            raise SystemExit(f"maximum normalized residual {worst:.5g} exceeds {args.max_residual:.5g}")


if __name__ == "__main__":
    main()
