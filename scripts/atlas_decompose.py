#!/usr/bin/env python3
"""
Decompose the AVI score for specific variants.

Why: AVI blends coding impact (AlphaMissense) with regulatory
predictions. N1002K scores at the 99th percentile in BUB1B, but that is
uninformative on its own, because an AlphaMissense of 0.92 would produce
a high AVI whether or not anything regulatory is happening. The thesis
says the damage is protein-structural, not regulatory. The feature
importances say which term is carrying the score, so this is the run
that actually tests it.

Usage:
    python atlas_decompose.py
    python atlas_decompose.py --variants chr15:40209701:T>G=L737Ter
"""

import argparse, os, sys
import numpy as np

from alphagenome.atlas import atlas
from alphagenome.data import genome

DEFAULT = [
    ("chr15:40209701:T>G", "BUB1B L737Ter"),
    ("chr15:40220612:T>G", "BUB1B N1002K"),
    ("chr22:20996720:C>G", "LZTR1 Y748Ter"),
]
SCORERS = ["AVI_SCORE", "AVI_SCORE_FEATURE_IMPORTANCE", "AVI_SCORE_MODEL_FEATURES"]
GENE_OF = {"chr15": "BUB1B", "chr22": "LZTR1", "chr5": "SLC34A1"}
# One gene per chromosome only holds while we work three loci. --gene
# overrides it, which is needed the moment two candidates share a
# chromosome. Passing the wrong gene name silently changes what the
# gene-centric scorers report, so it is worth being explicit.


def parse(spec):
    c, p, ch = spec.split(":")
    ref, alt = ch.split(">")
    return genome.Variant(chromosome=c, position=int(p),
                          reference_bases=ref, alternate_bases=alt)


def show(ad, label, scorer):
    print(f"\n--- {label} / {scorer} ---")
    if ad is None:
        print("    nothing returned")
        return
    x = np.asarray(ad.X)
    print(f"    X shape {x.shape}   obs{ad.obs.shape}  var{ad.var.shape}")
    print(f"    obs cols: {list(ad.obs.columns)}")
    print(f"    var cols: {list(ad.var.columns)}")
    # var rows name the features when this is a decomposition
    try:
        names = [str(v) for v in ad.var.iloc[:, 0].tolist()]
    except Exception:
        names = [f"f{i}" for i in range(x.shape[1] if x.ndim > 1 else 1)]
    flat = x.reshape(x.shape[0], -1) if x.ndim > 1 else x.reshape(1, -1)
    if flat.shape[0] == 1 and len(names) == flat.shape[1]:
        pairs = sorted(zip(names, flat[0]), key=lambda t: -abs(t[1]))
        print("    features by |contribution|:")
        for n, v in pairs[:25]:
            print(f"      {n:<46} {v:+.6g}")
    else:
        print(f"    values (first 12): {flat.ravel()[:12]}")
        print(f"    feature names (first 12): {names[:12]}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--variants", nargs="*", default=None,
                    help="chr:pos:REF>ALT=Label")
    ap.add_argument("--scorers", nargs="*", default=SCORERS)
    ap.add_argument("--gene", default=None,
                    help="override the chromosome->gene lookup")
    args = ap.parse_args()

    key = os.environ.get("ALPHAGENOME_API_KEY", "").strip()
    if not key:
        sys.exit("ALPHAGENOME_API_KEY is not set.")

    items = DEFAULT
    if args.variants:
        items = []
        for m in args.variants:
            spec, _, lab = m.partition("=")
            items.append((spec.strip(), lab.strip() or spec.strip()))

    c = atlas.create(key)
    for spec, label in items:
        v = parse(spec)
        gene = args.gene or GENE_OF.get(v.chromosome, None)
        print(f"\n================ {label}  {spec}  gene={gene} ================")
        for s in args.scorers:
            try:
                res = c.query_variant(variant=v, requested_scorers=[s],
                                      gene_names=[gene] if gene else None)
            except Exception as e:
                print(f"\n--- {label} / {s} --- FAILED {type(e).__name__}: {e}")
                continue
            ad = res.get(s) if isinstance(res, dict) else res
            show(ad, label, s)


if __name__ == "__main__":
    main()
