#!/usr/bin/env python3
"""
Score the two BUB1B variants with AlphaGenome.

  chr15:40209701 T>G   c.2210T>G  p.Leu737Ter   (allele 1, nonsense)
  chr15:40220612 T>G   c.3006T>G  p.Asn1002Lys  (allele 2, missense)

Why run this at all. For N1002K the result you WANT is negative: the
mechanism argument says the variant acts on protein folding, not on
transcription or splicing. A clean "no splicing effect, no expression
effect" is what licenses the protein-level argument, and a negative is
only worth anything if you actually ran the test. For L737Ter the
interesting readout is the opposite: a predicted drop in expression
would corroborate nonsense-mediated decay at the sequence level.

Build is GRCh38. AlphaGenome uses 'chr'-prefixed contigs; the source VCF
uses bare contigs, so the prefix is added explicitly below. Getting this
wrong returns nothing useful rather than an error.

Usage:
    export ALPHAGENOME_API_KEY='your-key'
    python3 score_bub1b_alphagenome.py --list-scorers
    python3 score_bub1b_alphagenome.py
    python3 score_bub1b_alphagenome.py --scorers RNA_SEQ SPLICE_SITES
"""

import argparse
import os
import sys
import traceback

# The quick start uses alphagenome.models.v1; older builds expose
# alphagenome.models directly. Try both rather than assume.
_IMPORT_NOTE = ""
try:
    from alphagenome.models.v1 import dna_client, variant_scorers
    _IMPORT_NOTE = "alphagenome.models.v1"
except ImportError:
    from alphagenome.models import dna_client, variant_scorers
    _IMPORT_NOTE = "alphagenome.models"
from alphagenome.data import genome

VARIANTS = [
    {"label": "L737Ter", "chrom": "chr15", "pos": 40209701, "ref": "T", "alt": "G",
     "note": "c.2210T>G nonsense, allele 1, expect NMD"},
    {"label": "N1002K", "chrom": "chr15", "pos": 40220612, "ref": "T", "alt": "G",
     "note": "c.3006T>G missense, allele 2, expect NO splicing effect"},
]

# Preference order. Intersected with whatever this build actually ships,
# so a renamed key degrades to a warning instead of a crash.
PREFERRED = [
    "RNA_SEQ", "SPLICE_SITES", "SPLICE_SITE_USAGE", "SPLICE_JUNCTIONS",
    "CAGE", "PROCAP", "POLYADENYLATION", "ATAC", "DNASE",
    "CHIP_HISTONE", "CHIP_TF", "CONTACT_MAPS",
]


def available_scorers():
    d = getattr(variant_scorers, "RECOMMENDED_VARIANT_SCORERS", None)
    if d is None:
        sys.exit("No RECOMMENDED_VARIANT_SCORERS in this build. Run with "
                 "--list-scorers after checking dir(variant_scorers).")
    return d


def tidy(ad, label, scorer_name, out_rows):
    """Flatten one AnnData of variant scores into long rows."""
    import numpy as np
    arr = np.asarray(ad.X)
    if arr.ndim == 1:
        arr = arr.reshape(-1, 1)
    obs = ad.obs.reset_index() if hasattr(ad, "obs") else None
    var = ad.var.reset_index() if hasattr(ad, "var") else None
    for i in range(arr.shape[0]):
        for j in range(arr.shape[1]):
            row = {"variant": label, "scorer": scorer_name,
                   "value": float(arr[i, j])}
            if obs is not None and i < len(obs):
                for c in obs.columns:
                    row["obs_" + str(c)] = obs.iloc[i][c]
            if var is not None and j < len(var):
                for c in var.columns:
                    row["var_" + str(c)] = var.iloc[j][c]
            out_rows.append(row)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--list-scorers", action="store_true",
                    help="print the scorer keys this build ships and exit")
    ap.add_argument("--scorers", nargs="*", default=None,
                    help="scorer keys to run (default: the preferred set)")
    ap.add_argument("--out", default="alphagenome_bub1b_scores.tsv")
    args = ap.parse_args()

    print(f"import path: {_IMPORT_NOTE}")
    scorers = available_scorers()

    if args.list_scorers:
        print(f"\n{len(scorers)} scorer keys in this build:")
        for k in sorted(scorers):
            print(f"  {k}")
        return

    key = os.environ.get("ALPHAGENOME_API_KEY", "").strip()
    if not key:
        sys.exit("ALPHAGENOME_API_KEY is not set. export it first; do not "
                 "paste the key into a shared file or a chat.")

    wanted = args.scorers if args.scorers else PREFERRED
    chosen = [k for k in wanted if k in scorers]
    missing = [k for k in wanted if k not in scorers]
    if missing:
        print(f"not in this build, skipping: {', '.join(missing)}")
    if not chosen:
        sys.exit("None of the requested scorers exist here. Run "
                 "--list-scorers to see what does.")
    print(f"scoring with: {', '.join(chosen)}\n")

    model = dna_client.create(key)
    out_rows, summary = [], []

    for v in VARIANTS:
        print(f"=== {v['label']}  {v['chrom']}:{v['pos']} {v['ref']}>{v['alt']}")
        print(f"    {v['note']}")
        variant = genome.Variant(
            chromosome=v["chrom"], position=v["pos"],
            reference_bases=v["ref"], alternate_bases=v["alt"],
        )
        interval = variant.reference_interval.resize(
            dna_client.SEQUENCE_LENGTH_1MB)

        for name in chosen:
            try:
                res = model.score_variant(
                    interval=interval, variant=variant,
                    variant_scorers=[scorers[name]],
                )
            except Exception as e:
                print(f"    {name:<20} FAILED  {type(e).__name__}: {e}")
                summary.append((v["label"], name, None, f"failed: {e}"))
                continue

            items = res if isinstance(res, (list, tuple)) else [res]
            before = len(out_rows)
            for ad in items:
                try:
                    tidy(ad, v["label"], name, out_rows)
                except Exception:
                    print(f"    {name:<20} scored but could not flatten:")
                    traceback.print_exc(limit=1)
            vals = [abs(r["value"]) for r in out_rows[before:]]
            peak = max(vals) if vals else float("nan")
            print(f"    {name:<20} ok  {len(vals):>5} tracks  "
                  f"peak |score| {peak:.4g}")
            summary.append((v["label"], name, peak, f"{len(vals)} tracks"))

    if out_rows:
        import pandas as pd
        df = pd.DataFrame(out_rows)
        df.to_csv(args.out, sep="\t", index=False)
        print(f"\nwrote {args.out}  ({len(df)} rows, {len(df.columns)} cols)")

    print("\n=== peak absolute score per variant and modality ===")
    for label, name, peak, note in summary:
        shown = "n/a" if peak is None else f"{peak:.4g}"
        print(f"  {label:<9} {name:<20} {shown:>12}   {note}")

    print("\nHow to read this. For N1002K a LOW peak on the splicing and "
          "RNA_SEQ scorers is the result that supports the mechanism "
          "argument, because it rules out a transcriptional or splicing "
          "explanation and leaves the protein-level one. Compare the two "
          "variants against each other rather than against an absolute "
          "threshold, and run SpliceAI on the same two variants so the "
          "splicing claim rests on two independent methods.")


if __name__ == "__main__":
    main()
