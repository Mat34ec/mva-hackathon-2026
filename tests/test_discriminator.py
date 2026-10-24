#!/usr/bin/env python3
"""
Validate the copy-number / mismapping discriminator against known ground truth.

Runs with **no patient data and no API key**. Builds a synthetic VCF with four
contigs of known composition (tests/make_synthetic_vcf.py), runs the unmodified
analysis code from scripts/aneu.py against it, and asserts that the method
behaves as the reports claim.

The claim under test, from Track 1 report §5 and Track 2 report §7:

    Depth alone cannot separate a real mosaic chromosome gain from paralogous
    mismapping. The shape of the heterozygous allele-balance distribution can,
    because a gain widens it SYMMETRICALLY (the extra copy duplicates either
    allele with equal probability) while mismapping only ever adds
    reference-matching reads and so widens the LOW side alone.

The fixture gives chrGAIN and chrMAP the *identical* depth lift on purpose. If
the discriminator is doing real work it must separate them anyway; if it is
secretly just reading depth, this test fails.

    python tests/test_discriminator.py        # standalone, prints a table
    pytest tests/                             # or under pytest
"""

import os
import subprocess
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts"))

HERE = os.path.dirname(os.path.abspath(__file__))
FIXTURE = os.path.join(HERE, "synthetic.vcf.gz")

# Ground truth built into make_synthetic_vcf.py
TRUTH = {
    "chrNORM": "clean diploid",
    "chrGAIN": "mosaic trisomy, f=0.30",
    "chrMAP": "paralogous mismapping",
    "chrHEMI": "hemizygous",
}


def ensure_fixture():
    if not os.path.exists(FIXTURE):
        r = subprocess.run(
            [sys.executable, os.path.join(HERE, "make_synthetic_vcf.py"),
             "--out", FIXTURE],
            capture_output=True, text=True)
        if r.returncode != 0 or not os.path.exists(FIXTURE):
            # Surface the real reason. Swallowing this makes a build failure
            # look like a test failure, which wastes a reviewer's time.
            raise SystemExit(
                "Could not build the test fixture.\n"
                f"  exit code: {r.returncode}\n"
                f"  stdout: {r.stdout.strip()[-600:]}\n"
                f"  stderr: {r.stderr.strip()[-600:]}\n"
                "  Most likely cause: pysam is missing (pip install pysam numpy),\n"
                "  or the tests/ directory is not writable."
            )
    return FIXTURE


def analyse():
    import aneu
    vcf = ensure_fixture()
    res = aneu.run(list(TRUTH), vcf)
    out = {}
    for c, d in res.items():
        w = d["whole"]
        out[c] = {
            "dp": w["dp_med"],
            "lo": w["ab_out_lo"],        # fraction of het sites below AB 0.40
            "hi": w["ab_out_hi"],        # fraction above AB 0.60
            "outside": w["ab_outside"],  # both tails together
            "het_frac": w["het_frac"],
            "n": w["n_used"],
        }
        out[c]["sym"] = out[c]["hi"] / out[c]["lo"] if out[c]["lo"] else float("nan")
    base = out["chrNORM"]["dp"]
    for c in out:
        out[c]["ratio"] = out[c]["dp"] / base
    return out


def report(st):
    print(f"\n{'contig':<9} {'truth':<24} {'dp':>5} {'ratio':>6} "
          f"{'<0.40':>7} {'>0.60':>7} {'both':>7} {'hi/lo':>6}")
    for c, t in TRUTH.items():
        s = st[c]
        print(f"{c:<9} {t:<24} {s['dp']:>5.0f} {s['ratio']:>6.3f} "
              f"{s['lo']:>7.3f} {s['hi']:>7.3f} {s['outside']:>7.3f} "
              f"{s['sym']:>6.2f}")


def run_checks(st):
    checks = []

    def chk(name, ok, detail):
        checks.append((name, ok, detail))

    g, m, n, h = st["chrGAIN"], st["chrMAP"], st["chrNORM"], st["chrHEMI"]

    # 1. The method resolves a real one-copy change. chrX in a male is the
    #    equivalent built-in positive control in the real analysis.
    chk("hemizygous contig detected at ~0.50 depth ratio",
        0.45 <= h["ratio"] <= 0.56,
        f"ratio={h['ratio']:.3f}")

    # 2. Depth cannot discriminate: the fixture gave both the same lift.
    chk("depth alone CANNOT separate gain from mismapping",
        abs(g["ratio"] - m["ratio"]) < 0.05,
        f"gain={g['ratio']:.3f} map={m['ratio']:.3f}")

    # 3. Neither can total spread.
    chk("total spread alone CANNOT separate them either",
        abs(g["outside"] - m["outside"]) < 0.08,
        f"gain={g['outside']:.3f} map={m['outside']:.3f}")

    # 4. The discriminator does. A true gain lifts BOTH tails.
    chk("true gain raises the HIGH tail above baseline",
        g["hi"] > n["hi"] * 1.3,
        f"gain hi={g['hi']:.3f} vs norm hi={n['hi']:.3f}")

    # 5. Mismapping does not.
    chk("mismapping leaves the HIGH tail at baseline",
        m["hi"] <= n["hi"] * 1.3,
        f"map hi={m['hi']:.3f} vs norm hi={n['hi']:.3f}")

    # 6. The symmetry statistic separates them cleanly.
    chk("symmetry statistic separates gain from mismapping",
        g["sym"] > 0.70 and m["sym"] < 0.50,
        f"gain hi/lo={g['sym']:.2f}  map hi/lo={m['sym']:.2f}")

    # 7. Sanity: the clean contig looks clean.
    chk("clean diploid contig is symmetric and at baseline depth",
        0.95 <= n["ratio"] <= 1.05 and n["sym"] > 0.70,
        f"ratio={n['ratio']:.3f} hi/lo={n['sym']:.2f}")

    return checks


def main():
    st = analyse()
    report(st)
    checks = run_checks(st)
    print()
    bad = 0
    for name, ok, detail in checks:
        print(f"  [{'PASS' if ok else 'FAIL'}] {name}\n         {detail}")
        bad += not ok
    print()
    if bad:
        print(f"{bad} of {len(checks)} checks FAILED")
        return 1
    print(f"all {len(checks)} checks passed — "
          "the discriminator separates copy-number change from mismapping\n"
          "on contigs that are indistinguishable by depth.")
    return 0


# pytest entry points
def test_discriminator():
    st = analyse()
    for name, ok, detail in run_checks(st):
        assert ok, f"{name}: {detail}"


if __name__ == "__main__":
    sys.exit(main())
