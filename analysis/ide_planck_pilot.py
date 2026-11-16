"""Reproducible Planck 2018 TT,TE,EE+lowE pipeline check for CAMB-IDE.

Run with the isolated Python 3.11 environment documented in
``analysis/ide_planck_pilot_results.md``.
The four-dimensional chain is a software pilot, not a posterior constraint.
"""

import argparse
import json
import math
import time
from importlib.metadata import version
from pathlib import Path

from cobaya.model import get_model
from cobaya.run import run

ROOT = Path(__file__).resolve().parents[1]
PACKAGES = ROOT / "external" / "planck-packages"
OUTPUT = ROOT / "external" / "planck-pilot"
LIKELIHOODS = (
    "planck_2018_highl_plik.TTTEEE_lite_native",
    "planck_2018_lowl.TT",
    "planck_2018_lowl.EE",
)
FIXED = {
    "ombh2": 0.02237,
    "mnu": 0.06,
    "nnu": 3.046,
    "As": 2.10e-9,
    "ns": 0.9649,
    "tau": 0.0544,
    "A_planck": 1.0,
}
BASELINE = {"V0_over_rho100": 0.31, "omc_early": 0.12, "beta": 0.0, "lam": 0.0}
PRIORS = {
    "V0_over_rho100": {"min": 0.20, "max": 0.42},
    "omc_early": {"min": 0.10, "max": 0.14},
    "beta": {"min": -0.1, "max": 0.1},
    "lam": {"min": 0.0, "max": 1.0},
}
PROPOSAL = {"V0_over_rho100": 0.002, "omc_early": 0.0003, "beta": 0.003, "lam": 0.01}


def make_info(camb_root, kind, point, *, accuracy=1.0, boosts=None, sample=False):
    """Use exactly the same likelihood data and CMB numerical settings for both models."""
    extra_args = {
        "lens_potential_accuracy": 0,
        "max_eta_k": 18000,
        "DoLateRadTruncation": False,
    }
    if accuracy != 1.0:
        extra_args.update(AccuracyBoost=accuracy, lSampleBoost=accuracy, lAccuracyBoost=accuracy)
    if boosts:
        extra_args.update(boosts)
    params = dict(FIXED)
    if kind == "ide":
        extra_args["dark_energy_model"] = "IDEExponential"
        if sample:
            for name, bounds in PRIORS.items():
                params[name] = {"prior": bounds, "ref": point[name], "proposal": PROPOSAL[name]}
        else:
            params.update(point)
        params.update(H0=None, omch2=None)
    else:
        if sample:
            raise ValueError("The four-input pilot samples only IDE")
        if "H0" not in point:
            raise ValueError("Clean CAMB requires the H0 derived from the IDE baseline")
        params.update(H0=point["H0"], omch2=point["omc_early"])
    return {
        "packages_path": str(PACKAGES),
        "theory": {"camb": {"path": str(camb_root), "extra_args": extra_args, "stop_at_error": True}},
        "likelihood": {name: None for name in LIKELIHOODS},
        "params": params,
    }


def evaluate(camb_root, kind, point, *, accuracy=1.0, boosts=None):
    info = make_info(camb_root, kind, point, accuracy=accuracy, boosts=boosts)
    started = time.monotonic()
    with get_model(info) as model:
        posterior = model.logposterior({}, cached=False)
        names = list(model.likelihood)
        loglikes = {name: float(value) for name, value in zip(names, posterior.loglikes)}
        derived_names = list(model.parameterization.derived_params())
        derived = {name: float(value) for name, value in zip(derived_names, posterior.derived)}
    if not loglikes or not all(math.isfinite(value) for value in loglikes.values()):
        raise RuntimeError(f"Planck likelihood returned non-finite values: {loglikes}")
    return {
        "kind": kind,
        "point": point,
        "accuracy_boost": accuracy,
        "numerical_settings": info["theory"]["camb"]["extra_args"],
        "loglikes": loglikes,
        "total_loglike": sum(loglikes.values()),
        "derived": derived,
        "elapsed_seconds": time.monotonic() - started,
        "cobaya_version": version("cobaya"),
        "camb_root": str(camb_root),
    }


def pilot(camb_root, samples, accuracy):
    initial = {**BASELINE, "lam": 0.1}
    info = make_info(camb_root, "ide", initial, sample=True, accuracy=accuracy)
    info["sampler"] = {
        "mcmc": {
            "max_samples": samples,
            "burn_in": 0,
            "learn_proposal": False,
            "measure_speeds": False,
            "proposal_scale": 0.6,
            "seed": 20261115,
            "output_every": 1,
        }
    }
    OUTPUT.mkdir(parents=True, exist_ok=True)
    info["output"] = str(OUTPUT / "ide_four_input_pilot")
    started = time.monotonic()
    _, sampler = run(info, no_mpi=True, force=True)
    accepted = int(sampler.n())
    acceptance = float(sampler.get_acceptance_rate())
    result = {
        "status": "completed",
        "accepted_steps": accepted,
        "acceptance_rate": acceptance,
        "elapsed_seconds": time.monotonic() - started,
        "computation_failures": 0,
        "cobaya_version": version("cobaya"),
        "accuracy_boost": accuracy,
        "initial_point": initial,
        "output_prefix": info["output"],
        "interpretation": "Pipeline check only; these samples are not an IDE posterior constraint.",
    }
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("evaluate", "pilot", "compare"))
    parser.add_argument("--camb-root", type=Path, default=ROOT)
    parser.add_argument("--kind", choices=("ide", "lcdm"), default="ide")
    parser.add_argument("--v0", type=float, default=BASELINE["V0_over_rho100"])
    parser.add_argument("--omc-early", type=float, default=BASELINE["omc_early"])
    parser.add_argument("--beta", type=float, default=BASELINE["beta"])
    parser.add_argument("--lam", type=float, default=BASELINE["lam"])
    parser.add_argument("--h0", type=float)
    parser.add_argument("--tight", action="store_true", help="Use 1.5 for all three CAMB accuracy boosts")
    parser.add_argument("--accuracy", type=float, help="Override AccuracyBoost, lSampleBoost and lAccuracyBoost")
    parser.add_argument("--accuracy-boost", type=float, help="Override AccuracyBoost only")
    parser.add_argument("--lsample-boost", type=float, help="Override lSampleBoost only")
    parser.add_argument("--laccuracy-boost", type=float, help="Override lAccuracyBoost only")
    parser.add_argument("--samples", type=int, default=12)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--reference", type=Path, help="JSON from clean CAMB at the same physical LCDM point")
    parser.add_argument("--ide", type=Path, help="JSON from IDE at beta=lambda=0")
    args = parser.parse_args()

    if args.action == "compare":
        if not args.reference or not args.ide:
            parser.error("compare requires --reference and --ide")
        reference = json.loads(args.reference.read_text())
        ide = json.loads(args.ide.read_text())
        differences = {name: ide["loglikes"][name] - value for name, value in reference["loglikes"].items()}
        result = {"loglike_differences": differences, "total_difference": sum(differences.values())}
    elif args.action == "pilot":
        result = pilot(args.camb_root.resolve(), args.samples, args.accuracy or 3.0)
    else:
        point = {"V0_over_rho100": args.v0, "omc_early": args.omc_early, "beta": args.beta, "lam": args.lam}
        if args.h0 is not None:
            point["H0"] = args.h0
        boosts = {
            name: value
            for name, value in (
                ("AccuracyBoost", args.accuracy_boost),
                ("lSampleBoost", args.lsample_boost),
                ("lAccuracyBoost", args.laccuracy_boost),
            )
            if value is not None
        }
        result = evaluate(
            args.camb_root.resolve(),
            args.kind,
            point,
            accuracy=args.accuracy or (1.5 if args.tight else 1.0),
            boosts=boosts,
        )
    rendered = json.dumps(result, indent=2, sort_keys=True)
    print(rendered)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered + "\n")


if __name__ == "__main__":
    main()
