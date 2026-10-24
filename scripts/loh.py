#!/usr/bin/env python3
"""
Test for subclonal loss of heterozygosity, i.e. tumour-derived DNA in the
blood sample.

Why this matters here. The mechanistic model says chromosomal instability
from BUB1B supplies the somatic second hit that a monoallelic inactivating
LZTR1 allele needs. That hit is in the tumour, and we have only germline
blood. But if any tumour DNA reached the blood, the lost LZTR1 allele
would show as allelic imbalance at heterozygous sites across the locus.
So the hypothesis is testable with the data in hand, and a negative is
worth having with a stated floor.

The statistic. At a true germline heterozygous site, alt-read count is
Binomial(DP, 0.5), so z = (alt - DP/2)/sqrt(DP/4) is standard normal and
var(z) = 1. Subclonal LOH at fraction f splits sites to AB 1/(2-f) and
(1-f)/(2-f), inflating var(z) above 1 SYMMETRICALLY, because which allele
survives is independent of which is reference.

Paralogous mismapping also inflates var(z), but only on the negative side,
because mismapped reads add reference-matching reads. So var(z) is split
by sign. Symmetric inflation is copy-number; one-sided is mapping. This is
the same discriminator that settled the aneuploidy question, reused here
on purpose so the two results rest on one logic.
"""
import json
import os

import numpy as np
import pysam

D = os.environ["HOME"] + "/mnt/MVA Hackathon 2026 - Dataset"
VCF = D + "/WGS_EX2312012_HGWCNDSX7.vcf.gz"

REGIONS = {
    "BUB1B":  ("15", 40158985, 40223137),
    "LZTR1":  ("22", 20979980, 21001037),
    # Controls on the same chromosomes, matched for size, away from the
    # pericentromeric bins the aneuploidy scan flagged.
    "ctrl_15": ("15", 60000000, 60064000),
    "ctrl_22": ("22", 37000000, 37021000),
    "ctrl_2":  ("2",  60000000, 60064000),
    "ctrl_7":  ("7",  80000000, 80064000),
}


def zs(chrom, start, end, v, s):
    out = []
    for rec in v.fetch(chrom, start, end):
        if rec.filter.keys() and "PASS" not in rec.filter.keys():
            continue
        if not rec.alts or len(rec.alts) != 1:
            continue
        if len(rec.ref) != 1 or len(rec.alts[0]) != 1:
            continue
        g = rec.samples[s]
        gt = g.get("GT")
        if gt is None or None in gt or len(set(gt)) < 2:
            continue                      # heterozygous sites only
        gq = g.get("GQ")
        if gq is None or gq < 30:
            continue
        ad = g.get("AD")
        if not ad or len(ad) < 2:
            continue
        r, a = int(ad[0]), int(ad[1])
        dp = r + a
        if dp < 20:                       # z is unreliable at low depth
            continue
        out.append(((a - dp / 2) / np.sqrt(dp / 4), dp, rec.pos, r, a))
    return out


def report(name, rows):
    if len(rows) < 8:
        print(f"  {name:<9} only {len(rows)} usable het sites, underpowered")
        return None
    z = np.array([r[0] for r in rows])
    pos_v = float(np.var(z[z > 0], ddof=1)) if (z > 0).sum() > 2 else float("nan")
    neg_v = float(np.var(z[z < 0], ddof=1)) if (z < 0).sum() > 2 else float("nan")
    print(f"  {name:<9} n={len(z):>4}  var(z)={z.var(ddof=1):>5.2f}  "
          f"mean(z)={z.mean():>+5.2f}  var+={pos_v:>5.2f}  var-={neg_v:>5.2f}  "
          f"medDP={int(np.median([r[1] for r in rows]))}")
    return {"n": len(z), "var": float(z.var(ddof=1)), "mean": float(z.mean()),
            "var_pos": pos_v, "var_neg": neg_v}


if __name__ == "__main__":
    v = pysam.VariantFile(VCF)
    s = list(v.header.samples)[0]
    print("var(z) = 1.0 is the null for clean germline heterozygosity.")
    print("Symmetric inflation (var+ and var- both up) = copy number.")
    print("One-sided (only var- up) = mismapping.\n")
    res = {}
    for name, (c, a, b) in REGIONS.items():
        res[name] = report(name, zs(c, a, b, v, s))

    print("\nPer-site detail at the two genes that matter:")
    for name in ("BUB1B", "LZTR1"):
        c, a, b = REGIONS[name]
        rows = zs(c, a, b, v, s)
        print(f"\n  {name}: {len(rows)} het sites")
        for z_, dp, pos, r, al in rows:
            flag = "  <<" if abs(z_) > 2.5 else ""
            print(f"    {pos:>10}  {r:>3},{al:<3} DP{dp:>4}  "
                  f"AB {al/dp:.3f}  z {z_:>+6.2f}{flag}")

    with open(os.environ["HOME"] + "/an/loh.json", "w") as fh:
        json.dump(res, fh, indent=1)
