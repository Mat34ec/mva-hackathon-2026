#!/usr/bin/env python3
"""
AlphaGenome Atlas lookups for the BUB1B locus.

The Atlas holds precomputed scores for every possible single-nucleotide
change in the genome, so this is a lookup rather than inference. Two
things that buys us:

  1. Real calibration. Instead of hand-built controls, both patient
     variants get a percentile against every possible change in BUB1B.
  2. The non-coding door. The patient's rare intronic and regulatory
     variants across the locus can be screened in one batch, including
     c.2679-1026A>G, which the research bundle lists as unassessed.

Remember what this cannot do: AlphaGenome does not model protein
folding or thermodynamic stability. These results can only rule
alternatives out. They cannot support the destabilisation claim itself.

Usage:
    python atlas_bub1b.py --probe              # what Atlas offers, small query
    python atlas_bub1b.py --interval           # full BUB1B, percentiles
    python atlas_bub1b.py --variants my.tsv    # look up a list (chrom pos ref alt)
"""

import argparse, csv, os, sys
from collections import defaultdict

from alphagenome.atlas import atlas
from alphagenome.data import genome

GENE = "BUB1B"
CHROM = "chr15"
GENE_START, GENE_END = 40160983, 40221137          # GRCh38 gene body
REAL = {"L737Ter": (40209701, "T", "G"),
        "N1002K":  (40220612, "T", "G")}


def client():
    key = os.environ.get("ALPHAGENOME_API_KEY", "").strip()
    if not key:
        sys.exit("ALPHAGENOME_API_KEY is not set.")
    return atlas.create(key)


def probe(c):
    md = c.scorer_metadata()
    print(f"{len(md)} scorers in the Atlas:\n")
    for name in sorted(md):
        print(f"  {name}")
    print("\nsmall interval probe: 2 kb around N1002K, RNA_SEQ only")
    iv = genome.Interval(chromosome=CHROM, start=40219612, end=40221612)
    df = c.query_interval(interval=iv, requested_scorers=["AVI_SCORE"],
                          gene_names=[GENE], progress_bar=False)
    report(df, "probe")


def report(obj, label, depth=0, maxdepth=3):
    """Print the shape of whatever the Atlas returned, without assuming it."""
    pad = "  " * depth
    t = type(obj).__name__
    if isinstance(obj, dict):
        print(f"{pad}[{label}] dict, {len(obj)} keys")
        for k in list(obj)[:25]:
            print(f"{pad}  key: {k!r} -> {type(obj[k]).__name__}")
        if depth < maxdepth:
            for k, v in list(obj.items())[:4]:
                report(v, f"{label}[{k!r}]", depth + 1, maxdepth)
        return
    if isinstance(obj, (list, tuple)):
        print(f"{pad}[{label}] {t}, len={len(obj)}")
        if obj and depth < maxdepth:
            report(obj[0], f"{label}[0]", depth + 1, maxdepth)
        return
    print(f"{pad}[{label}] {t}")
    for attr in ("shape", "columns", "dtype"):
        if hasattr(obj, attr):
            v = getattr(obj, attr)
            print(f"{pad}  {attr}: {list(v) if attr == 'columns' else v}")
    if hasattr(obj, "obs") and hasattr(obj, "var"):
        print(f"{pad}  AnnData obs{obj.obs.shape} var{obj.var.shape}")
        print(f"{pad}  obs cols: {list(obj.obs.columns)[:12]}")
        print(f"{pad}  var cols: {list(obj.var.columns)[:12]}")
        try:
            print(obj.obs.head(4).to_string()[:900])
        except Exception:
            pass
        return
    if hasattr(obj, "head"):
        try:
            print(obj.head(6).to_string()[:1800])
        except Exception:
            pass


def pct(series, value):
    """Percentile of value within series, as a percentage."""
    import numpy as np
    a = np.asarray(sorted(abs(x) for x in series if x == x))
    if len(a) == 0:
        return float("nan")
    return 100.0 * (a < abs(value)).sum() / len(a)


def interval_run(c, scorers, out):
    iv = genome.Interval(chromosome=CHROM, start=GENE_START, end=GENE_END)
    print(f"querying {CHROM}:{GENE_START}-{GENE_END} "
          f"({(GENE_END - GENE_START)/1000:.1f} kb) for {GENE}")
    print(f"scorers: {', '.join(scorers)}\n")
    df = c.query_interval(interval=iv, requested_scorers=scorers,
                          gene_names=[GENE], progress_bar=True)
    report(df, "interval")
    try:
        df.to_csv(out, sep="\t", index=False)
        print(f"\nwrote {out}")
    except Exception as e:
        print(f"could not write {out}: {e}")
    return df


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--probe", action="store_true")
    ap.add_argument("--interval", action="store_true")
    ap.add_argument("--variants")
    ap.add_argument("--scorers", nargs="*",
                    default=["AVI_SCORE", "RNA_SEQ", "SPLICE_SITES", "SPLICE_SITE_USAGE",
                             "SPLICE_JUNCTIONS", "ATAC", "DNASE",
                             "CHIP_HISTONE", "CHIP_TF", "CAGE", "PROCAP",
                             "POLYADENYLATION", "CONTACT_MAPS"])
    ap.add_argument("--out", default="atlas_bub1b_interval.tsv")
    args = ap.parse_args()

    c = client()
    if args.probe:
        probe(c)
        return
    if args.interval:
        interval_run(c, args.scorers, args.out)
        return
    if args.variants:
        rows = [r for r in csv.reader(open(args.variants), delimiter="\t")
                if r and not r[0].startswith("#")]
        vs = []
        for r in rows:
            ch = r[0] if r[0].startswith("chr") else "chr" + r[0]
            vs.append(genome.Variant(chromosome=ch, position=int(r[1]),
                                     reference_bases=r[2], alternate_bases=r[3]))
        print(f"looking up {len(vs)} variants")
        df = c.query_variants(variants=vs, requested_scorers=args.scorers,
                              gene_names=[GENE], progress_bar=True)
        report(df, "variants")
        df.to_csv("atlas_variant_lookup.tsv", sep="\t", index=False)
        print("\nwrote atlas_variant_lookup.tsv")
        return
    ap.print_help()


if __name__ == "__main__":
    main()
