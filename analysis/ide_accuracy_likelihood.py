"""Evaluate one CAMB accuracy sub-control with the Planck pilot likelihood.

For example::

    OMP_NUM_THREADS=1 python analysis/ide_accuracy_likelihood.py \
        --control SourcekAccuracyBoost --output external/planck-pilot/edge_pos_sourcek_1p25.json
"""

import argparse
import json
from pathlib import Path

from ide_planck_pilot import BASELINE, ROOT, evaluate

CONTROLS = ("SourcekAccuracyBoost", "TimeStepBoost", "IntkAccuracyBoost")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--control", required=True, choices=CONTROLS)
    parser.add_argument("--value", type=float, default=1.25)
    parser.add_argument("--accuracy", type=float, default=4.0)
    parser.add_argument("--beta", type=float, default=0.1)
    parser.add_argument("--lam", type=float, default=1.0)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    if args.value <= 0 or args.accuracy <= 0:
        parser.error("accuracy controls must be positive")

    point = {**BASELINE, "beta": args.beta, "lam": args.lam}
    result = evaluate(
        ROOT,
        "ide",
        point,
        accuracy=args.accuracy,
        boosts={f"Accuracy.{args.control}": args.value},
    )
    rendered = json.dumps(result, indent=2, sort_keys=True)
    print(rendered)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered + "\n")


if __name__ == "__main__":
    main()
