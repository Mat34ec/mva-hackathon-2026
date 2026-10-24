#!/usr/bin/env python3
"""
Place the patient's variants in the distribution of every possible
single-nucleotide change across their own gene, using AlphaGenome Atlas
precomputed scores.

This replaces hand-built controls. Rather than asking "is 0.0086 small",
it answers "N1002K sits at the Nth percentile of all 180,000 possible
changes in BUB1B", which is a statement a reviewer can act on.

Gene coordinates come from Ensembl GRCh38, not hardcoded. The interval
is queried in chunks so a 60 kb gene does not go out as one request.

Usage:
  python atlas_percentile.py --gene BUB1B \
      --mark chr15:40209701:T>G=L737Ter --mark chr15:40220612:T>G=N1002K
  python atlas_percentile.py --gene LZTR1 --mark chr22:20996720:C>G=Y748Ter
  python atlas_percentile.py --gene BUB1B --scorers AVI_SCORE RNA_SEQ SPLICE_SITES
"""

import argparse, json, os, sys, urllib.request
import numpy as np

from alphagenome.atlas import atlas
from alphagenome.data import genome

ENSEMBL = "https://rest.ensembl.org"


def gene_span(symbol):
    url = f"{ENSEMBL}/lookup/symbol/homo_sapiens/{symbol}?expand=0"
    req = urllib.request.Request(url, headers={"Content-Type": "application/json"})
    with urllib.request.urlopen(req, timeout=30) as r:
        g = json.loads(r.read().decode())
    chrom = str(g["seq_region_name"])
    return ("chr" + chrom if not chrom.startswith("chr") else chrom,
            int(g["start"]), int(g["end"]), g["id"])


def to_frame(ad, scorer):
    """AnnData -> rows of (key, chrom, pos, ref, alt, score).

    obs['variant'] holds genome.Variant OBJECTS, not strings; they merely
    render as 'chr15:40220125:A>C'. Read the attributes, and fall back to
    parsing str() only if the object shape ever changes.
    """
    x = np.asarray(ad.X).reshape(len(ad.obs), -1)
    out = []
    for i, v in enumerate(ad.obs["variant"].tolist()):
        c = p_ = ref = alt = None
        if hasattr(v, "chromosome") and hasattr(v, "position"):
            c = str(v.chromosome)
            p_ = int(v.position)
            ref = str(getattr(v, "reference_bases", "") or "")
            alt = str(getattr(v, "alternate_bases", "") or "")
        else:
            try:
                c, ps, ch = str(v).split(":")
                ref, alt = ch.split(">")
                p_ = int(ps)
            except ValueError:
                c, p_, ref, alt = "", 0, "", ""
        out.append((f"{c}:{p_}:{ref}>{alt}", c, p_, ref, alt, float(x[i, 0])))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gene", required=True)
    ap.add_argument("--scorers", nargs="*", default=["AVI_SCORE"])
    ap.add_argument("--mark", action="append", default=[],
                    help="chr:pos:REF>ALT=Label, repeatable")
    ap.add_argument("--chunk", type=int, default=10000, help="bp per request")
    ap.add_argument("--pad", type=int, default=2000)
    ap.add_argument("--out")
    args = ap.parse_args()

    key = os.environ.get("ALPHAGENOME_API_KEY", "").strip()
    if not key:
        sys.exit("ALPHAGENOME_API_KEY is not set.")

    chrom, start, end, ensid = gene_span(args.gene)
    start, end = start - args.pad, end + args.pad
    print(f"{args.gene} ({ensid}) {chrom}:{start}-{end}  "
          f"{(end - start)/1000:.1f} kb\nscorers: {', '.join(args.scorers)}")

    marks = {}
    for m in args.mark:
        spec, _, label = m.partition("=")
        marks[spec.strip()] = label.strip() or spec.strip()

    c = atlas.create(key)
    rows = {s: [] for s in args.scorers}
    chunks = [(p, min(p + args.chunk, end)) for p in range(start, end, args.chunk)]
    for i, (a, b) in enumerate(chunks, 1):
        iv = genome.Interval(chromosome=chrom, start=a, end=b)
        try:
            res = c.query_interval(interval=iv, requested_scorers=args.scorers,
                                   gene_names=[args.gene], progress_bar=False)
        except Exception as e:
            print(f"  chunk {i}/{len(chunks)} {a}-{b} FAILED {type(e).__name__}: {e}")
            continue
        for s in args.scorers:
            ad = res.get(s)
            if ad is not None:
                rows[s] += to_frame(ad, s)
        got = sum(len(rows[s]) for s in args.scorers)
        print(f"  chunk {i}/{len(chunks)} {a}-{b}  cumulative rows {got}")

    out = args.out or f"atlas_{args.gene}_scores.tsv"
    with open(out, "w") as fh:
        fh.write("scorer\tvariant\tchrom\tpos\tref\talt\tscore\n")
        for s in args.scorers:
            for v, ch, p, r, al, sc in rows[s]:
                fh.write(f"{s}\t{v}\t{ch}\t{p}\t{r}\t{al}\t{sc}\n")
    print(f"\nwrote {out}")

    for s in args.scorers:
        data = rows[s]
        if not data:
            print(f"\n=== {s}: nothing returned ===")
            continue
        vals = np.array([abs(d[5]) for d in data])
        print(f"\n=== {s} over {args.gene}: {len(vals):,} possible variants ===")
        qs = np.percentile(vals, [50, 90, 95, 99, 99.9])
        print(f"  |score| median {qs[0]:.5g}  p90 {qs[1]:.5g}  p95 {qs[2]:.5g}  "
              f"p99 {qs[3]:.5g}  p99.9 {qs[4]:.5g}  max {vals.max():.5g}")

        index = {d[0]: d for d in data}
        for spec, label in marks.items():
            hit = index.get(spec)
            if hit is None:
                alt = next((d for d in data if f"{d[1]}:{d[2]}:" in spec
                            and d[4] == spec.split(">")[-1]), None)
                hit = alt
            if hit is None:
                print(f"  {label:<10} NOT FOUND as {spec}")
                continue
            v = abs(hit[5])
            pctile = 100.0 * (vals < v).sum() / len(vals)
            rank = int((vals > v).sum()) + 1
            print(f"  {label:<10} score {hit[5]:+.6g}  |{v:.6g}|  "
                  f"percentile {pctile:.2f}  rank {rank:,} of {len(vals):,}")

        top = sorted(data, key=lambda d: -abs(d[5]))[:8]
        print(f"  most impactful positions in {args.gene}:")
        for d in top:
            print(f"    {d[0]:<26} {d[5]:+.6g}")


if __name__ == "__main__":
    main()
