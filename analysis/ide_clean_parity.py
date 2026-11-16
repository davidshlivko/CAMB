"""Compare the IDE Lambda-CDM limit with a separate, clean CAMB worktree.

Run, for example::

    OMP_NUM_THREADS=1 python analysis/ide_clean_parity.py \
        --clean-tree /private/tmp/camb-ide-clean-20261115

The script launches separate Python processes from each worktree, so they load
their own Fortran extensions. Both use the same physical densities, primordial
spectrum, RECFAST, scalar modes, linear lensing, and late-radiation setting.
"""

import argparse
import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

CASE_CODE = r"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path.cwd()))
import camb

kind, output_path, lmax_text, h0_text = sys.argv[1:]
lmax = int(lmax_text)
pars = camb.CAMBparams()
if kind == "ide":
    pars.set_ide_cosmology(
        V0_over_rho100=0.31,
        omc_early=0.12,
        beta=0.0,
        lam=0.0,
        ombh2=0.02237,
        mnu=0.06,
        tau=0.0544,
    )
else:
    pars.set_cosmology(H0=float(h0_text), ombh2=0.02237, omch2=0.12, mnu=0.06, tau=0.0544)
pars.InitPower.set_params(As=2.10e-9, ns=0.9649)
pars.set_accuracy(AccuracyBoost=1.0, lAccuracyBoost=1.0, DoLateRadTruncation=False)
pars.set_for_lmax(lmax, lens_potential_accuracy=0, nonlinear=False)
result = camb.get_results(pars)
z_distance = np.array([0.1, 0.5, 2.0, 10.0, 100.0, 1100.0, 2500.0])
z_thermal = np.geomspace(1.0, 5000.0, 500) - 1.0
thermal = result.get_background_redshift_evolution(
    z_thermal, vars=["x_e", "visibility"], format="array"
)
np.savez(
    output_path,
    camb_module=str(Path(camb.__file__).resolve()),
    H0=result.Params.H0,
    omch2=result.Params.omch2,
    hubble=result.hubble_parameter(z_distance),
    distance=result.angular_diameter_distance(z_distance),
    xe=thermal[:, 0],
    visibility=thermal[:, 1],
    unlensed=result.get_unlensed_scalar_cls(lmax),
    lensed=result.get_lensed_scalar_cls(lmax),
    pp=result.get_lens_potential_cls(lmax)[:, 0],
)
"""


def run_case(interpreter, tree, kind, output, lmax, h0=None):
    env = os.environ.copy()
    env.pop("PYTHONPATH", None)
    env["OMP_NUM_THREADS"] = "1"
    command = [str(interpreter), "-c", CASE_CODE, kind, str(output), str(lmax), "" if h0 is None else str(h0)]
    process = subprocess.run(
        command, cwd=tree, env=env, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, check=False
    )
    if process.returncode:
        raise RuntimeError(f"{kind} CAMB calculation failed in {tree}:\n{process.stdout}")


def max_scaled(actual, reference, floor_fraction=1e-4):
    scale = np.maximum(np.abs(reference), floor_fraction * np.max(np.abs(reference)))
    return float(np.max(np.abs(actual - reference) / scale))


def compare(ide, clean):
    report = {
        "H0_km_s_Mpc": float(ide["H0"]),
        "omch2_ide": float(ide["omch2"]),
        "hubble_max_relative": max_scaled(ide["hubble"], clean["hubble"], 0.0),
        "distance_max_relative": max_scaled(ide["distance"], clean["distance"], 0.0),
        "xe_max_scaled": max_scaled(ide["xe"], clean["xe"], 1e-6),
        "visibility_peak_normalized": float(
            np.max(np.abs(ide["visibility"] - clean["visibility"])) / np.max(clean["visibility"])
        ),
    }
    for name in ("unlensed", "lensed"):
        actual = ide[name][2:]
        reference = clean[name][2:]
        report[f"{name}_TT_max_scaled"] = max_scaled(actual[:, 0], reference[:, 0])
        report[f"{name}_EE_max_scaled"] = max_scaled(actual[:, 1], reference[:, 1])
        te_scale = np.sqrt(np.maximum(reference[:, 0] * reference[:, 1], 0.0))
        te_scale = np.maximum(te_scale, 1e-5 * np.max(te_scale))
        report[f"{name}_TE_over_sqrt_TT_EE"] = float(np.max(np.abs(actual[:, 3] - reference[:, 3]) / te_scale))
    report["PP_max_scaled"] = max_scaled(ide["pp"][2:], clean["pp"][2:])
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--clean-tree", required=True, type=Path, help="separate unmodified CAMB worktree with camblib")
    parser.add_argument("--ide-tree", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--lmax", type=int, default=1000)
    args = parser.parse_args()
    clean_tree = args.clean_tree.resolve()
    ide_tree = args.ide_tree.resolve()
    for tree in (clean_tree, ide_tree):
        if not (tree / "camb" / "camblib.so").exists():
            parser.error(f"{tree} has no built camb/camblib.so")
    if args.lmax < 2:
        parser.error("--lmax must be at least 2")

    with tempfile.TemporaryDirectory(prefix="ide-clean-parity-") as temp:
        ide_file = Path(temp) / "ide.npz"
        clean_file = Path(temp) / "clean.npz"
        run_case(sys.executable, ide_tree, "ide", ide_file, args.lmax)
        with np.load(ide_file) as ide:
            run_case(sys.executable, clean_tree, "clean", clean_file, args.lmax, float(ide["H0"]))
            with np.load(clean_file) as clean:
                ide_module = Path(str(ide["camb_module"])).resolve()
                clean_module = Path(str(clean["camb_module"])).resolve()
                if not ide_module.is_relative_to(ide_tree) or not clean_module.is_relative_to(clean_tree):
                    raise RuntimeError("The two subprocesses did not import CAMB from their requested worktrees")
                report = compare(ide, clean)
                report.update(ide_module=str(ide_module), clean_module=str(clean_module), lmax=args.lmax)
                print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
