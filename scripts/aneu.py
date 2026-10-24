#!/usr/bin/env python3
"""
Per-chromosome and binned test for mosaic aneuploidy, from depth and
allele balance jointly.

Why this and not the earlier shape test. The earlier run asked whether the
genome-wide allele-balance distribution was bimodal and found it unimodal
with a one-sided low tail that tracked segmental duplications, i.e.
mismapping. That was a global test. Aneuploidy is per chromosome, so a
global test can dilute a real single-chromosome event below notice. MVA is
characteristically whole-chromosome, so chromosome is the right unit.

The statistic. For a chromosome at copy number CN in a fraction f of cells:

  depth ratio      = 1 + f(CN-2)/2
  het allele balance splits into two modes, not one shifted mode:
      trisomy  -> (1+f)/(2+f) and 1/(2+f)
      monosomy -> 1/(2-f)     and (1-f)/(2-f)

So the MEAN allele balance stays at 0.5 by symmetry for both gain and
loss, and the signal is in the SPREAD. Testing the mean is the wrong test
and would return a false negative. We therefore track the robust spread
(IQR and the fraction of het sites outside 0.40-0.60) alongside depth, and
require the two to move together in a consistent direction before calling
anything. Chromosomes differ in GC content and mappability, which moves
depth on its own, so depth alone is not evidence.

Writes per-chromosome and 10 Mb-bin stats as JSON for later assembly.
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
DP_HIST_MAX = 300
AB_BINS = 200


def blank():
    return {
        "dp_hist": np.zeros(DP_HIST_MAX + 1, dtype=np.int64),
        "ab_hist": np.zeros(AB_BINS + 1, dtype=np.int64),
        "n_het": 0,
        "n_homalt": 0,
        "n_used": 0,
        "n_seen": 0,
    }


def add(acc, dp, ab, zyg):
    acc["n_seen"] += 1
    acc["n_used"] += 1
    acc["dp_hist"][min(dp, DP_HIST_MAX)] += 1
    if zyg == "het":
        acc["n_het"] += 1
        acc["ab_hist"][int(round(ab * AB_BINS))] += 1
    else:
        acc["n_homalt"] += 1


def hist_quantiles(h, qs, scale=1.0):
    """Quantiles straight from a count histogram, so nothing is stored."""
    tot = h.sum()
    if tot == 0:
        return [float("nan")] * len(qs)
    c = np.cumsum(h)
    out = []
    for q in qs:
        idx = int(np.searchsorted(c, q * tot))
        out.append(float(min(idx, len(h) - 1)) * scale)
    return out


def summarize(acc):
    if acc["n_used"] == 0:
        return None
    dp_q = hist_quantiles(acc["dp_hist"], [0.25, 0.5, 0.75])
    ab_q = hist_quantiles(acc["ab_hist"], [0.25, 0.5, 0.75], scale=1.0 / AB_BINS)
    ab = acc["ab_hist"]
    n_ab = int(ab.sum())
    lo = int(ab[: int(0.40 * AB_BINS)].sum())
    hi = int(ab[int(0.60 * AB_BINS) + 1 :].sum())
    centres = np.arange(len(ab)) / AB_BINS
    mean_ab = float((centres * ab).sum() / n_ab) if n_ab else float("nan")
    return {
        "n_used": acc["n_used"],
        "n_het": acc["n_het"],
        "n_homalt": acc["n_homalt"],
        "het_frac": acc["n_het"] / acc["n_used"],
        "dp_q25": dp_q[0], "dp_med": dp_q[1], "dp_q75": dp_q[2],
        "ab_q25": ab_q[0], "ab_med": ab_q[1], "ab_q75": ab_q[2],
        "ab_iqr": ab_q[2] - ab_q[0],
        "ab_mean": mean_ab,
        "ab_out_lo": lo / n_ab if n_ab else float("nan"),
        "ab_out_hi": hi / n_ab if n_ab else float("nan"),
        "ab_outside": (lo + hi) / n_ab if n_ab else float("nan"),
    }


def run(chroms, vcf):
    v = pysam.VariantFile(vcf)
    s = list(v.header.samples)[0]
    out = {}
    for c in chroms:
        whole, bins = blank(), {}
        try:
            it = v.fetch(c)
        except ValueError:
            print(f"  {c}: not in index, skipped", flush=True)
            continue
        for rec in it:
            # Biallelic PASS SNVs only. Indels and multiallelics carry
            # depth and allele-balance artefacts that would swamp a few
            # percent of mosaicism.
            if rec.filter.keys() and "PASS" not in rec.filter.keys():
                continue
            alts = rec.alts
            if not alts or len(alts) != 1:
                continue
            if len(rec.ref) != 1 or len(alts[0]) != 1:
                continue
            g = rec.samples[s]
            gt = g.get("GT")
            if gt is None or None in gt or all(a == 0 for a in gt):
                continue
            gq = g.get("GQ")
            if gq is None or gq < 30:
                continue
            ad = g.get("AD")
            if not ad or len(ad) < 2:
                continue
            ref_n, alt_n = int(ad[0]), int(ad[1])
            dp = ref_n + alt_n
            if dp < 10 or dp > DP_HIST_MAX:
                continue
            zyg = "het" if len(set(gt)) > 1 else "homalt"
            ab = alt_n / dp
            add(whole, dp, ab, zyg)
            b = rec.pos // BIN
            if b not in bins:
                bins[b] = blank()
            add(bins[b], dp, ab, zyg)
        w = summarize(whole)
        out[c] = {
            "whole": w,
            "bins": {str(k): summarize(a) for k, a in sorted(bins.items())},
        }
        if w:
            print(f"  {c}: {w['n_used']:>9,} sites  dp_med {w['dp_med']:>5.1f}  "
                  f"het_frac {w['het_frac']:.3f}  ab_med {w['ab_med']:.3f}  "
                  f"ab_iqr {w['ab_iqr']:.3f}  outside40-60 {w['ab_outside']:.3f}",
                  flush=True)
    return out


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("chroms", nargs="+", help="contig names to scan")
    ap.add_argument("--vcf", default=None, help="VCF path (else $MVA_VCF)")
    ap.add_argument("--out", default=None, help="JSON output path")
    args = ap.parse_args()
    res = run(args.chroms, resolve_vcf(args.vcf))
    path = args.out or "aneu_parts.json"
    allres = {}
    if os.path.exists(path):
        with open(path) as fh:
            allres = json.load(fh)
    allres.update(res)
    with open(path, "w") as fh:
        json.dump(allres, fh)
    print(f"\nstored {len(res)} chromosome(s); file now holds {len(allres)}")
