"""Compare IDE CMB spectra under independent numerical controls.

Run from the CAMB-IDE worktree with ``OMP_NUM_THREADS=1``. By default the
script checks one coupled point; run separate processes for the other points
to limit memory use. The script imports the local CAMB package and prints
JSON. The tight-coupling switch coefficient is hard-coded in CAMB;
``AccuracyBoost`` affects its high-k branch but also changes other settings,
so it is not an isolated tight-coupling test.
"""

import argparse
import gc
import importlib
import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
camb = importlib.import_module("camb")

POINTS = (
    (0.0, 0.0),
    (0.05, 0.5),
    (-0.05, 0.5),
    (0.1, 1.0),
    (-0.1, 1.0),
)


def calculate(beta, lam, lmax, changes):
    settings = {
        "a_init": 1e-10,
        "history_points": 2049,
        "integrate_tol": 1e-10,
        "accuracy": 2.0,
        "neutrino_q_boost": 1.0,
        "time_switch_boost": 1.0,
        "int_tol_boost": 1.0,
        "massive_nu_method": "Nu_best",
    }
    settings.update(changes)
    start = time.monotonic()
    pars = camb.set_params(
        dark_energy_model="IDEExponential",
        V0_over_rho100=0.31,
        omc_early=0.12,
        beta=beta,
        lam=lam,
        a_init=settings["a_init"],
        history_points=settings["history_points"],
        integrate_tol=settings["integrate_tol"],
        ombh2=0.02237,
        mnu=0.06,
        nnu=3.046,
        tau=0.0544,
        As=2.10e-9,
        ns=0.9649,
        lmax=lmax,
        max_eta_k=18000,
        lens_potential_accuracy=0,
        AccuracyBoost=settings["accuracy"],
        lSampleBoost=settings["accuracy"],
        lAccuracyBoost=settings["accuracy"],
    )
    pars.Accuracy.neutrino_q_boost = settings["neutrino_q_boost"]
    pars.Accuracy.TimeSwitchBoost = settings["time_switch_boost"]
    pars.Accuracy.IntTolBoost = settings["int_tol_boost"]
    pars.MassiveNuMethod = settings["massive_nu_method"]
    result = camb.get_results(pars)
    return {
        "cl": result.get_lensed_scalar_cls(lmax),
        "pp": result.get_lens_potential_cls(lmax)[:, 0],
        "H0": result.Params.H0,
        "omch2": result.Params.omch2,
        "elapsed_seconds": time.monotonic() - start,
        "settings": settings,
    }


def metrics(other, baseline, lmin=30):
    base = baseline["cl"][lmin:]
    trial = other["cl"][lmin:]
    tt_scale = np.maximum(np.abs(base[:, 0]), 1e-5 * np.max(np.abs(base[:, 0])))
    ee_scale = np.maximum(np.abs(base[:, 1]), 1e-5 * np.max(np.abs(base[:, 1])))
    te_scale = np.sqrt(np.maximum(base[:, 0] * base[:, 1], 0.0))
    te_scale = np.maximum(te_scale, 1e-5 * np.max(te_scale))
    pp_base = baseline["pp"][lmin:]
    pp_scale = np.maximum(np.abs(pp_base), 1e-5 * np.max(np.abs(pp_base)))
    scaled = {
        "TT": np.abs(trial[:, 0] - base[:, 0]) / tt_scale,
        "EE": np.abs(trial[:, 1] - base[:, 1]) / ee_scale,
        "TE_over_sqrt_TT_EE": np.abs(trial[:, 3] - base[:, 3]) / te_scale,
        "PP": np.abs(other["pp"][lmin:] - pp_base) / pp_scale,
    }
    return {
        name: {"max": float(np.max(values)), "p95": float(np.percentile(values, 95))} for name, values in scaled.items()
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--lmax", type=int, default=1000)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--point-index", type=int, choices=range(len(POINTS)), default=1)
    parser.add_argument("--all-points", action="store_true", help="Run all five points in one process")
    parser.add_argument("--variants", nargs="+")
    args = parser.parse_args()
    if args.lmax < 30:
        parser.error("lmax must be at least 30")
    variants = {
        "earlier_ai": {"a_init": 1e-12},
        "denser_history": {"history_points": 4097},
        "tighter_background": {"integrate_tol": 1e-11},
        "neutrino_q_x2": {"neutrino_q_boost": 2.0},
        "neutrino_full_integral": {"massive_nu_method": "Nu_int"},
        "later_switches": {"time_switch_boost": 2.0},
        "ode_tol_x2": {"int_tol_boost": 2.0},
        "all_accuracy_x1p5": {"accuracy": 3.0},
    }
    if args.variants:
        unknown = set(args.variants) - set(variants)
        if unknown:
            parser.error(f"unknown variants: {', '.join(sorted(unknown))}")
        variants = {name: variants[name] for name in args.variants}
    output = {"lmax": args.lmax, "points": []}
    points = POINTS if args.all_points else [POINTS[args.point_index]]
    for beta, lam in points:
        print(f"Running beta={beta:+.3f}, lam={lam:.3f}", file=sys.stderr, flush=True)
        baseline = calculate(beta, lam, args.lmax, {})
        row = {
            "beta": beta,
            "lam": lam,
            "baseline_H0": baseline["H0"],
            "baseline_omch2": baseline["omch2"],
            "baseline_seconds": baseline["elapsed_seconds"],
            "variants": {},
        }
        for name, changes in variants.items():
            print(f"  {name}", file=sys.stderr, flush=True)
            trial = calculate(beta, lam, args.lmax, changes)
            row["variants"][name] = {
                "settings": changes,
                "seconds": trial["elapsed_seconds"],
                "H0_difference": trial["H0"] - baseline["H0"],
                "spectra": metrics(trial, baseline),
            }
            del trial
            gc.collect()
            if args.output:
                args.output.parent.mkdir(parents=True, exist_ok=True)
                args.output.write_text(json.dumps(output | {"points": output["points"] + [row]}, indent=2) + "\n")
        output["points"].append(row)
    rendered = json.dumps(output, indent=2, sort_keys=True)
    print(rendered)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered + "\n")


if __name__ == "__main__":
    main()
