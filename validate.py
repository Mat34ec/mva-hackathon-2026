#!/usr/bin/env python3
"""
One-command validation of this submission's central method.

    python validate.py

Runs with **no patient data and no API key**. Builds a synthetic VCF with four
contigs of known composition, runs the unmodified analysis code from
scripts/aneu.py against it, and checks that the discriminator does what the
reports claim: separate a real mosaic chromosome gain from paralogous
mismapping on contigs that are indistinguishable by read depth.

Requires only pysam and numpy:

    pip install pysam numpy

Other entry points:

    python validate.py --panel     rebuild the 12-gene panel live from Ensembl
                                   (needs network, no data and no key)
    python validate.py --clean     remove generated fixtures and outputs
"""

import argparse
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))


def run(cmd, **kw):
    print("$ " + " ".join(cmd))
    return subprocess.run(cmd, check=False, **kw).returncode


def do_validate():
    missing = []
    for m in ("pysam", "numpy"):
        try:
            __import__(m)
        except ImportError:
            missing.append(m)
    if missing:
        print(f"Missing: {', '.join(missing)}\n  pip install {' '.join(missing)}")
        return 1
    return run([sys.executable, os.path.join(HERE, "tests", "test_discriminator.py")])


def do_panel():
    return run([sys.executable, os.path.join(HERE, "scripts", "make_panel_bed.py"),
                "--build", "GRCh38", "--pad", "5000", "--out", "mva_panel"])


def do_clean():
    import shutil
    gone = []
    for p in ("tests/synthetic.vcf.gz", "tests/synthetic.vcf.gz.tbi",
              "aneu_parts.json", "loh_bins.json"):
        f = os.path.join(HERE, p)
        if os.path.exists(f):
            os.remove(f)
            gone.append(p)
    for root, dirs, _ in os.walk(HERE):
        for d in list(dirs):
            if d == "__pycache__":
                shutil.rmtree(os.path.join(root, d), ignore_errors=True)
                dirs.remove(d)
    print("removed: " + (", ".join(gone) if gone else "nothing to remove"))
    return 0


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--panel", action="store_true", help="rebuild the panel from Ensembl")
    ap.add_argument("--clean", action="store_true", help="remove generated files")
    a = ap.parse_args()
    if a.clean:
        sys.exit(do_clean())
    if a.panel:
        sys.exit(do_panel())
    sys.exit(do_validate())
