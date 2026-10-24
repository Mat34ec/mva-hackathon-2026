#!/usr/bin/env python3
"""
Genome-wide scan for subclonal loss of heterozygosity, per 10 Mb bin.

Corrects a statistic error in loh.py: splitting var(z) by the sign of z
gives the variance of a truncated half-normal, which is not comparable to
the null value of 1. The correct side-specific statistic is E[z^2], which
equals 1 on BOTH sides under clean germline heterozygosity. So:

  E[z^2] up on both sides  -> copy-number imbalance (LOH or gain)
  E[z^2] up only for z<0   -> paralogous mismapping, reference-matching
                              reads added on one side only

Framing the question genome wide rather than at LZTR1 alone. A soft tissue
sarcoma almost always carries arm-level or chromosome-level copy change, so
if tumour DNA reached this blood sample at a detectable fraction it should
appear as an overdispersed bin SOMEWHERE, not only at the gene we care
about. A clean genome-wide result is therefore a much stronger statement
than a clean result at one locus, and it is the honest way to ask.
"""
import argparse
import json
import os
import sys

import numpy as np
import pysam

def resolve_vcf(cli_value=None):
    """Locate the VCF without hardcoding anyone's directory layout.

    Order: --vcf argument, then $MVA_VCF, then the dataset default. Keeping
    this resolvable is what lets the script run against the synthetic test
    fixture in tests/, so the method can be validated without access to the
    controlled-access patient data.
    """
    import glob
    if cli_value:
        return cli_value
    env = os.environ.get("MVA_VCF")
    if env:
        return env
    for guess in (
        "WGS_EX2312012_HGWCNDSX7.vcf.gz",
        os.path.join(os.environ.get("HOME", ""),
                     "mnt/MVA Hackathon 2026 - Dataset",
                     "WGS_EX2312012_HGWCNDSX7.vcf.gz"),
    ):
        if os.path.exists(guess):
            return guess
    hits = glob.glob("*.vcf.gz") + glob.glob("data/*.vcf.gz")
    if hits:
        return hits[0]
    raise SystemExit(
        "No VCF found. Pass --vcf PATH, or set MVA_VCF, or see DATA.md "
        "for how to obtain the dataset under its access terms."
    )
BIN = 10_000_000


def run(chroms, vcf):
    v = pysam.VariantFile(vcf)
    s = list(v.header.samples)[0]
    out = {}
    for c in chroms:
        acc = {}
        try:
            it = v.fetch(c)
        except ValueError:
            continue
        for rec in it:
            if rec.filter.keys() and "PASS" not in rec.filter.keys():
                continue
            if not rec.alts or len(rec.alts) != 1:
                continue
            if len(rec.ref) != 1 or len(rec.alts[0]) != 1:
                continue
            g = rec.samples[s]
            gt = g.get("GT")
            if gt is None or None in gt or len(set(gt)) < 2:
                continue
            gq = g.get("GQ")
            if gq is None or gq < 30:
                continue
            ad = g.get("AD")
            if not ad or len(ad) < 2:
                continue
            r, a = int(ad[0]), int(ad[1])
            dp = r + a
            if dp < 20 or dp > 200:
                continue
            z = (a - dp / 2) / np.sqrt(dp / 4)
            b = rec.pos // BIN
            if b not in acc:
                acc[b] = {"n": 0, "sz": 0.0, "z2p": 0.0, "np_": 0,
                          "z2n": 0.0, "nn": 0}
            d = acc[b]
            d["n"] += 1
            d["sz"] += z
            if z > 0:
                d["z2p"] += z * z
                d["np_"] += 1
            elif z < 0:
                d["z2n"] += z * z
                d["nn"] += 1
        out[c] = {
            str(b): {
                "n": d["n"],
                "mean_z": d["sz"] / d["n"],
                "ez2_pos": d["z2p"] / d["np_"] if d["np_"] else None,
                "ez2_neg": d["z2n"] / d["nn"] if d["nn"] else None,
            }
            for b, d in sorted(acc.items()) if d["n"] >= 300
        }
        print(f"  {c}: {len(out[c])} bins with >=300 het sites", flush=True)
    return out


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("chroms", nargs="+", help="contig names to scan")
    ap.add_argument("--vcf", default=None, help="VCF path (else $MVA_VCF)")
    ap.add_argument("--out", default=None, help="JSON output path")
    args = ap.parse_args()
    res = run(args.chroms, resolve_vcf(args.vcf))
    p = args.out or "loh_bins.json"
    allr = json.load(open(p)) if os.path.exists(p) else {}
    allr.update(res)
    json.dump(allr, open(p, "w"))
    print(f"stored {len(res)}; file holds {len(allr)}")
