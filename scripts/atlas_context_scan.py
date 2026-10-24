#!/usr/bin/env python3
"""
Score every SNV the patient carries in the 1 Mb context window around
each gene, using AlphaGenome Atlas lookup.

Why 1 Mb and not the whole chromosome: AlphaGenome's context window is
1 Mb, fixed. A variant further than 500 kb from BUB1B cannot affect
BUB1B in this model, so a chromosome-wide scan would answer a different
question. This window is the complete answer to the one that matters,
and it is where both published MVA regulatory second alleles sat (44 kb
upstream in one case, a few kb in the other).

Input:  research/context_1mb_variants.tsv  (already extracted)
Output: research/atlas_context_scores.tsv  + a ranked summary

Usage:
    python atlas_context_scan.py
    python atlas_context_scan.py --window BUB1B_1Mb --top 40
"""

import argparse, csv, os, sys
import numpy as np

from alphagenome.atlas import atlas
from alphagenome.data import genome

IN = "context_1mb_variants.tsv"
OUT = "atlas_context_scores.tsv"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", default=IN)
    ap.add_argument("--out", default=OUT)
    ap.add_argument("--window", default=None, help="BUB1B_1Mb or LZTR1_1Mb")
    ap.add_argument("--scorer", default="AVI_SCORE")
    ap.add_argument("--batch", type=int, default=500)
    ap.add_argument("--top", type=int, default=30)
    args = ap.parse_args()

    key = os.environ.get("ALPHAGENOME_API_KEY", "").strip()
    if not key:
        sys.exit("ALPHAGENOME_API_KEY is not set.")

    rows = []
    with open(args.input) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["kind"] != "snv":
                continue
            if args.window and r["window"] != args.window:
                continue
            rows.append(r)
    if not rows:
        sys.exit(f"no SNVs matched in {args.input}")

    byw = {}
    for r in rows:
        byw.setdefault(r["window"], []).append(r)
    for w, v in byw.items():
        print(f"{w}: {len(v)} SNVs")
    print(f"scorer: {args.scorer}\n")

    c = atlas.create(key)
    scored = []
    for w, group in byw.items():
        for i in range(0, len(group), args.batch):
            chunk = group[i:i + args.batch]
            variants = [genome.Variant(
                chromosome="chr" + r["chrom"], position=int(r["pos"]),
                reference_bases=r["ref"], alternate_bases=r["alt"]) for r in chunk]
            try:
                res = c.query_variants(variants=variants,
                                       requested_scorers=[args.scorer],
                                       progress_bar=False)
            except Exception as e:
                print(f"  {w} batch {i//args.batch + 1} FAILED "
                      f"{type(e).__name__}: {e}")
                continue
            ad = res.get(args.scorer) if isinstance(res, dict) else res
            if ad is None:
                print(f"  {w} batch {i//args.batch + 1}: no result")
                continue
            x = np.asarray(ad.X).reshape(len(ad.obs), -1)
            obs = ad.obs["variant"].tolist()
            gene = ad.obs["gene_name"].tolist() if "gene_name" in ad.obs else None
            for j, v in enumerate(obs):
                pos = int(getattr(v, "position", 0))
                scored.append({
                    "window": w,
                    "chrom": str(getattr(v, "chromosome", "")),
                    "pos": pos,
                    "ref": str(getattr(v, "reference_bases", "")),
                    "alt": str(getattr(v, "alternate_bases", "")),
                    "gene": (gene[j] if gene and j < len(gene) else ""),
                    "avi": float(x[j, 0]),
                })
            print(f"  {w} batch {i//args.batch + 1}: "
                  f"{len(obs)} scored, cumulative {len(scored)}")

    if not scored:
        sys.exit("nothing scored")

    zyg = {(r["chrom"], int(r["pos"]), r["ref"], r["alt"]): r["zyg"] for r in rows}
    for s in scored:
        s["zyg"] = zyg.get((s["chrom"].replace("chr", ""), s["pos"],
                            s["ref"], s["alt"]), "")

    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(scored[0]), delimiter="\t")
        w.writeheader(); w.writerows(scored)
    print(f"\nwrote {args.out}  ({len(scored)} rows)")

    known = {(40209701, "T", "G"): "BUB1B L737Ter",
             (40220612, "T", "G"): "BUB1B N1002K",
             (20996720, "C", "G"): "LZTR1 Tyr748Ter"}
    scored.sort(key=lambda s: -abs(s["avi"]))
    print(f"\n=== top {args.top} by |AVI| across both 1 Mb windows ===")
    print(f"{'rank':>4}  {'position':<18} {'change':<8} {'zyg':<4} "
          f"{'AVI':>9}  gene")
    for i, s in enumerate(scored[:args.top], 1):
        tag = known.get((s["pos"], s["ref"], s["alt"]), "")
        mark = f"   <<< {tag}" if tag else ""
        print(f"{i:>4}  {s['chrom']}:{s['pos']:<11} {s['ref']}>{s['alt']:<5} "
              f"{s['zyg']:<4} {s['avi']:+9.4f}  {s['gene']}{mark}")

    vals = np.array([abs(s["avi"]) for s in scored])
    print(f"\nacross {len(vals):,} carried SNVs: median {np.median(vals):.4f}  "
          f"p95 {np.percentile(vals, 95):.4f}  p99 {np.percentile(vals, 99):.4f}  "
          f"max {vals.max():.4f}")
    print("For reference, the BUB1B saturation scan gave p99 = 1.125 and "
          "p99.9 = 2.205 over all possible changes in that gene.")


if __name__ == "__main__":
    main()
