#!/usr/bin/env python3
"""
Build the MVA candidate-gene panel BED from Ensembl REST.

Coordinates are fetched live rather than hardcoded, so the panel is
reproducible and auditable for the hackathon write-up. Emits both GRCh37
and GRCh38 by default, and both 'chr'-prefixed and bare-contig versions,
so the output works regardless of the VCF's naming convention.

Usage:
    python make_panel_bed.py                      # both builds, 5 kb pad
    python make_panel_bed.py --build GRCh38
    python make_panel_bed.py --pad 10000
    python make_panel_bed.py --out mva_panel

Outputs per build:
    <out>.<build>.nochr.bed     e.g.  15  40159161  40212812  BUB1B|T1
    <out>.<build>.chr.bed       e.g.  chr15 ...
    <out>.<build>.provenance.tsv    gene, ensembl id, locus, build, date
"""

import argparse
import json
import sys
import time
import urllib.error
import urllib.request
from datetime import datetime, timezone

# Tier 1 = numbered MVA types. Tier 2 = reported candidates.
# Tier 3 = exploratory. MVK deliberately excluded (mevalonic aciduria,
# unrelated disorder that shares the 'MVA' abbreviation).
PANEL = [
    ("BUB1B",    "T1"),
    ("CEP57",    "T1"),
    ("TRIP13",   "T1"),
    ("CENATAC",  "T1"),
    ("SLF2",     "T1"),
    ("SMC5",     "T1"),
    ("MAD1L1",   "T1"),
    ("MAD2L1BP", "T2"),
    ("BUB1",     "T2"),
    ("CEP192",   "T2"),
    ("TRIM37",   "T3"),
    ("BUB3",     "T3"),
]

# Fallback symbols if the primary lookup misses on a given build.
ALIASES = {
    "CENATAC":  ["CCDC84"],
    "SLF2":     ["FAM178A"],
    "MAD2L1BP": ["CMT2", "C17orf46"],
    "TRIP13":   [],
}

SERVERS = {
    "GRCh38": "https://rest.ensembl.org",
    "GRCh37": "https://grch37.rest.ensembl.org",
}

PRIMARY_CONTIGS = {str(c) for c in range(1, 23)} | {"X", "Y", "MT"}


def lookup(server, symbol, retries=3):
    """Return the Ensembl gene record for a symbol, or None."""
    url = f"{server}/lookup/symbol/homo_sapiens/{symbol}?expand=0"
    req = urllib.request.Request(url, headers={"Content-Type": "application/json"})
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(req, timeout=30) as resp:
                return json.loads(resp.read().decode())
        except urllib.error.HTTPError as e:
            if e.code == 400:          # symbol not found on this build
                return None
            if e.code == 429:          # rate limited
                time.sleep(2 * (attempt + 1))
                continue
            if attempt == retries - 1:
                raise
            time.sleep(1.5 * (attempt + 1))
        except urllib.error.URLError:
            if attempt == retries - 1:
                raise
            time.sleep(1.5 * (attempt + 1))
    return None


def resolve(server, symbol):
    """Try the primary symbol, then any known aliases."""
    rec = lookup(server, symbol)
    if rec:
        return rec, symbol
    for alias in ALIASES.get(symbol, []):
        rec = lookup(server, alias)
        if rec:
            return rec, alias
    return None, None


def build_panel(build, pad, out_prefix):
    server = SERVERS[build]
    rows, problems = [], []
    stamp = datetime.now(timezone.utc).strftime("%Y-%m-%d")

    print(f"\n=== {build} ===")
    for symbol, tier in PANEL:
        rec, used = resolve(server, symbol)
        time.sleep(0.12)   # stay well under Ensembl's rate limit

        if rec is None:
            problems.append(f"{symbol}: no Ensembl record on {build}")
            print(f"  MISS  {symbol}")
            continue

        contig = str(rec["seq_region_name"])
        if contig not in PRIMARY_CONTIGS:
            problems.append(
                f"{symbol}: maps to non-primary contig '{contig}' on {build}, "
                "check manually before trusting"
            )

        asm = rec.get("assembly_name", "")
        if asm and not asm.startswith(build):
            problems.append(f"{symbol}: server returned assembly '{asm}', expected {build}")

        # Ensembl is 1-based inclusive; BED is 0-based half-open.
        start = max(0, int(rec["start"]) - 1 - pad)
        end = int(rec["end"]) + pad

        rows.append({
            "symbol": symbol,
            "used_symbol": used,
            "ensembl_id": rec["id"],
            "contig": contig,
            "start": start,
            "end": end,
            "strand": "+" if rec.get("strand", 1) == 1 else "-",
            "biotype": rec.get("biotype", ""),
            "tier": tier,
        })
        note = "" if used == symbol else f"  (via alias {used})"
        print(f"  ok    {symbol:<9} {contig}:{start}-{end}{note}")

    if not rows:
        print(f"  nothing resolved for {build}, skipping file write")
        return problems

    def sort_key(r):
        c = r["contig"]
        return (int(c) if c.isdigit() else 100 + ord(c[0]), r["start"])

    rows.sort(key=sort_key)

    for prefixed in (False, True):
        suffix = "chr" if prefixed else "nochr"
        path = f"{out_prefix}.{build}.{suffix}.bed"
        with open(path, "w") as fh:
            for r in rows:
                contig = ("chr" + r["contig"]) if prefixed else r["contig"]
                fh.write(f"{contig}\t{r['start']}\t{r['end']}\t"
                         f"{r['symbol']}|{r['tier']}\t0\t{r['strand']}\n")
        print(f"  wrote {path}")

    prov = f"{out_prefix}.{build}.provenance.tsv"
    with open(prov, "w") as fh:
        fh.write("gene\tsymbol_queried\tensembl_id\tcontig\tstart_bed\tend_bed\t"
                 "strand\tbiotype\ttier\tbuild\tpad_bp\tretrieved_utc\tsource\n")
        for r in rows:
            fh.write(f"{r['symbol']}\t{r['used_symbol']}\t{r['ensembl_id']}\t"
                     f"{r['contig']}\t{r['start']}\t{r['end']}\t{r['strand']}\t"
                     f"{r['biotype']}\t{r['tier']}\t{build}\t{pad}\t{stamp}\t"
                     f"{SERVERS[build]}\n")
    print(f"  wrote {prov}")

    total_bp = sum(r["end"] - r["start"] for r in rows)
    print(f"  {len(rows)}/{len(PANEL)} genes, {total_bp/1e6:.2f} Mb "
          f"({100*total_bp/3.1e9:.4f}% of the genome)")

    return problems


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--build", choices=["GRCh37", "GRCh38", "both"], default="both")
    ap.add_argument("--pad", type=int, default=5000,
                    help="bp added either side of the gene body (default 5000, "
                         "catches promoter and UTR-proximal variants)")
    ap.add_argument("--out", default="mva_panel")
    args = ap.parse_args()

    builds = ["GRCh37", "GRCh38"] if args.build == "both" else [args.build]

    all_problems = []
    for b in builds:
        all_problems += build_panel(b, args.pad, args.out)

    if all_problems:
        print("\n!! review these before using the panel:")
        for p in all_problems:
            print(f"   - {p}")
        sys.exit(1)

    print("\nAll genes resolved cleanly.")


if __name__ == "__main__":
    main()
