#!/usr/bin/env python3
"""
Calibrate the two BUB1B variants against synthetic neighbours.

The question this answers: is L737Ter's splice-junction score actually
elevated, or is 0.125 simply what any single-base change scores at this
locus? The earlier run had no answer, because the only internal
background was the other genes in the window, and for L737Ter the
junction scorer returned BUB1B alone.

Controls are generated, not taken from the patient. For each real
variant we score:
  - the two other possible substitutions at the same position
    (separates an allele-specific effect from a position-specific one)
  - single-base changes at nearby offsets in the same exon
    (gives a local background distribution)

Reference bases come from Ensembl REST on GRCh38, so the controls are
real substitutions rather than invented ref/alt pairs.

Usage:
    export ALPHAGENOME_API_KEY=...        # or  $env:ALPHAGENOME_API_KEY
    python score_controls_alphagenome.py
    python score_controls_alphagenome.py --scorers SPLICE_JUNCTIONS SPLICE_SITES
    python score_controls_alphagenome.py --offsets -20 -10 10 20
"""

import argparse, json, os, statistics as st, sys, urllib.request
from collections import defaultdict

try:
    from alphagenome.models.v1 import dna_client, variant_scorers
except ImportError:
    from alphagenome.models import dna_client, variant_scorers
from alphagenome.data import genome

GENE = "BUB1B"
REAL = [
    {"label": "L737Ter", "chrom": "chr15", "pos": 40209701, "ref": "T", "alt": "G"},
    {"label": "N1002K",  "chrom": "chr15", "pos": 40220612, "ref": "T", "alt": "G"},
]
DEFAULT_OFFSETS = [-30, -20, -10, -5, 5, 10, 20, 30]
ENSEMBL = "https://rest.ensembl.org"


def ref_window(chrom, centre, flank):
    """Reference bases centred on a position, GRCh38, 1-based inclusive."""
    c = chrom[3:] if chrom.startswith("chr") else chrom
    url = (f"{ENSEMBL}/sequence/region/human/"
           f"{c}:{centre - flank}..{centre + flank}?coord_system_version=GRCh38")
    req = urllib.request.Request(url, headers={"Content-Type": "application/json"})
    with urllib.request.urlopen(req, timeout=30) as r:
        seq = json.loads(r.read().decode())["seq"].upper()
    return {centre - flank + i: b for i, b in enumerate(seq)}


def build_controls(v, offsets, flank):
    win = ref_window(v["chrom"], v["pos"], flank)
    got = win.get(v["pos"])
    if got and got != v["ref"]:
        print(f"  !! reference at {v['pos']} is {got}, the VCF says {v['ref']}. "
              "Check the build before trusting anything below.")
    out = []
    for alt in "ACGT":
        if alt not in (v["ref"], v["alt"]):
            out.append({"label": f"{v['label']}~same{v['ref']}>{alt}",
                        "chrom": v["chrom"], "pos": v["pos"],
                        "ref": v["ref"], "alt": alt, "kind": "same position"})
    for off in offsets:
        p = v["pos"] + off
        r = win.get(p)
        if not r or r not in "ACGT":
            continue
        alt = {"A": "G", "G": "A", "C": "T", "T": "C"}[r]   # transition
        out.append({"label": f"{v['label']}~{off:+d}", "chrom": v["chrom"],
                    "pos": p, "ref": r, "alt": alt, "kind": "neighbour"})
    return out


def score(model, scorers, name, v):
    variant = genome.Variant(chromosome=v["chrom"], position=v["pos"],
                             reference_bases=v["ref"], alternate_bases=v["alt"])
    interval = variant.reference_interval.resize(dna_client.SEQUENCE_LENGTH_1MB)
    res = model.score_variant(interval=interval, variant=variant,
                              variant_scorers=[scorers[name]])
    items = res if isinstance(res, (list, tuple)) else [res]
    vals = []
    import numpy as np
    for ad in items:
        arr = np.asarray(ad.X)
        if arr.ndim == 1:
            arr = arr.reshape(-1, 1)
        names = list(ad.obs["gene_name"]) if "gene_name" in getattr(ad, "obs", {}) else None
        for i in range(arr.shape[0]):
            if names and i < len(names) and names[i] != GENE:
                continue          # keep BUB1B tracks only
            for j in range(arr.shape[1]):
                vals.append(abs(float(arr[i, j])))
    return vals


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--scorers", nargs="*", default=["SPLICE_JUNCTIONS"])
    ap.add_argument("--offsets", nargs="*", type=int, default=DEFAULT_OFFSETS)
    ap.add_argument("--flank", type=int, default=40)
    ap.add_argument("--out", default="alphagenome_controls.tsv")
    args = ap.parse_args()

    key = os.environ.get("ALPHAGENOME_API_KEY", "").strip()
    if not key:
        sys.exit("ALPHAGENOME_API_KEY is not set.")
    scorers = variant_scorers.RECOMMENDED_VARIANT_SCORERS
    chosen = [s for s in args.scorers if s in scorers]
    if not chosen:
        sys.exit(f"none of {args.scorers} exist in this build")

    jobs = []
    for v in REAL:
        print(f"building controls for {v['label']} at {v['chrom']}:{v['pos']}")
        jobs.append(dict(v, kind="REAL"))
        jobs += build_controls(v, args.offsets, args.flank)
    print(f"\n{len(jobs)} variants x {len(chosen)} scorer(s) = "
          f"{len(jobs) * len(chosen)} calls\n")

    model = dna_client.create(key)
    rows, by = [], defaultdict(list)
    for name in chosen:
        for j in jobs:
            try:
                vals = score(model, scorers, name, j)
            except Exception as e:
                print(f"  {j['label']:<22} {name:<18} FAILED {type(e).__name__}: {e}")
                continue
            if not vals:
                print(f"  {j['label']:<22} {name:<18} no {GENE} track returned")
                continue
            med, mx = st.median(vals), max(vals)
            rows.append({"variant": j["label"], "kind": j["kind"], "scorer": name,
                         "chrom": j["chrom"], "pos": j["pos"], "ref": j["ref"],
                         "alt": j["alt"], "n_tracks": len(vals),
                         "median_abs": med, "max_abs": mx})
            by[name].append((med, j["label"], j["kind"]))
            print(f"  {j['label']:<22} {name:<18} n={len(vals):>6}  "
                  f"median {med:.5f}  max {mx:.5f}")

    if rows:
        import csv
        with open(args.out, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
            w.writeheader(); w.writerows(rows)
        print(f"\nwrote {args.out}")

    for name in chosen:
        vals = sorted(by[name], reverse=True)
        print(f"\n=== {name}: ranked by median |score| over {GENE} tracks ===")
        for i, (med, label, kind) in enumerate(vals, 1):
            mark = "  <<< PATIENT VARIANT" if kind == "REAL" else ""
            print(f"  {i:>2}. {med:.5f}  {label:<22} {kind}{mark}")
        ctrl = [m for m, _, k in vals if k != "REAL"]
        if ctrl:
            print(f"\n  controls: median {st.median(ctrl):.5f}, "
                  f"max {max(ctrl):.5f}, n={len(ctrl)}")
            for med, label, kind in vals:
                if kind == "REAL":
                    over = sum(1 for c in ctrl if c >= med)
                    print(f"  {label}: {over} of {len(ctrl)} controls score as high "
                          f"or higher  ->  empirical p ~ {(over + 1) / (len(ctrl) + 1):.3f}")


if __name__ == "__main__":
    main()
