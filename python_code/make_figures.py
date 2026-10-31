"""Redraw figures from saved results without rerunning the analyses (each of 01-05 already draws its own
figures when it finishes).

Example:
    python make_figures.py                  # every figure
    python make_figures.py --step 04        # Fig 8 and S3 Fig
    python make_figures.py --only fig8,fig9
"""

from __future__ import annotations

import argparse
import sys

from lrps.figures import FIGURES, make_figures


def main():
    sys.stdout.reconfigure(encoding="utf-8")
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--step", choices=sorted(FIGURES), help="figures of one analysis script (default: all)")
    ap.add_argument("--only", help="comma-separated figure names, e.g. fig1,fig8")
    args = ap.parse_args()
    make_figures(args.step, only=args.only.split(",") if args.only else None)


if __name__ == "__main__":
    main()
