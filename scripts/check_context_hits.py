#!/usr/bin/env python3
"""
Identify the top-scoring variants from the 1 Mb context scan.

The scan ranked 2,626 carried SNVs by AVI. The two BUB1B variants came
first and second, but several unidentified variants scored in the same
range as N1002K. Before claiming the regulatory door is shut, each needs
a gene, a consequence and a population frequency, because a high AVI on
a common polymorphism is background rather than a finding.

The one that matters most is chr22:21009455, about 8 kb past LZTR1: a
rare high-impact regulatory variant there would be a germline second hit.

Usage:
    python check_context_hits.py            # top 15
    python check_context_hits.py --top 30
"""

import argparse, csv, json, urllib.request

REST = "https://rest.ensembl.org"
KNOWN = {(40209701, "T", "G"): "BUB1B L737Ter",
         (40220612, "T", "G"): "BUB1B N1002K",
         (20996720, "C", "G"): "LZTR1 Tyr748Ter"}


def vep(lines):
    url = f"{REST}/vep/human/region?hgvs=1&canonical=1&mane=1&numbers=1"
    req = urllib.request.Request(
        url, data=json.dumps({"variants": lines}).encode(),
        headers={"Content-Type": "application/json", "Accept": "application/json"})
    with urllib.request.urlopen(req, timeout=180) as r:
        return json.loads(r.read().decode())


def max_af(rec):
    best, pop = None, ""
    for cv in rec.get("colocated_variants") or []:
        for allele, pops in (cv.get("frequencies") or {}).items():
            for p, v in pops.items():
                try:
                    v = float(v)
                except (TypeError, ValueError):
                    continue
                if best is None or v > best:
                    best, pop = v, p
    return best, pop


def rsid(rec):
    for cv in rec.get("colocated_variants") or []:
        if str(cv.get("id", "")).startswith("rs"):
            return cv["id"]
    return ""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--scores", default="atlas_context_scores.tsv")
    ap.add_argument("--top", type=int, default=15)
    args = ap.parse_args()

    rows = list(csv.DictReader(open(args.scores), delimiter="\t"))
    rows.sort(key=lambda r: -abs(float(r["avi"])))
    top = rows[:args.top]

    lines, meta = [], {}
    for r in top:
        c = r["chrom"].replace("chr", "")
        spec = f"{c} {r['pos']} . {r['ref']} {r['alt']} . . ."
        lines.append(spec)
        meta[spec] = r

    res = vep(lines)
    by_input = {x.get("input", ""): x for x in res}

    print(f"{'AVI':>8}  {'position':<18} {'chg':<6} {'zyg':<4} {'gnomAD AF':>11}  "
          f"gene / consequence")
    print("-" * 112)
    flagged = []
    for spec in lines:
        r = meta[spec]
        rec = by_input.get(spec)
        avi = float(r["avi"])
        tag = KNOWN.get((int(r["pos"]), r["ref"], r["alt"]), "")
        if rec is None:
            print(f"{avi:+8.4f}  {r['chrom']}:{r['pos']:<11} "
                  f"{r['ref']}>{r['alt']:<3} {r['zyg']:<4} {'?':>11}  NO VEP RESULT")
            continue
        af, pop = max_af(rec)
        afs = "absent" if af is None else f"{af:.5g}"
        genes, cons = set(), set()
        for t in rec.get("transcript_consequences", []) or []:
            if t.get("gene_symbol"):
                genes.add(t["gene_symbol"])
            for x in t.get("consequence_terms", []):
                cons.add(x)
        if not genes:
            for rf in rec.get("regulatory_feature_consequences", []) or []:
                genes.add("[regulatory:" + str(rf.get("biotype", "?")) + "]")
        sev = rec.get("most_severe_consequence", "")
        rare = (af is None or af < 0.01)
        mark = "  <<< RARE" if (rare and not tag) else (f"  <<< {tag}" if tag else "")
        if rare and not tag:
            flagged.append((avi, r, afs, sorted(genes), sev, rsid(rec)))
        print(f"{avi:+8.4f}  {r['chrom']}:{r['pos']:<11} "
              f"{r['ref']}>{r['alt']:<3} {r['zyg']:<4} {afs:>11}  "
              f"{','.join(sorted(genes)) or '-'} / {sev}{mark}")

    print(f"\n=== {len(flagged)} rare (AF < 1%) high-AVI variants other than the "
          f"three known ones ===")
    for avi, r, afs, genes, sev, rs in flagged:
        print(f"  {r['chrom']}:{r['pos']} {r['ref']}>{r['alt']} {r['zyg']}  "
              f"AVI {avi:+.4f}  AF {afs}  {rs}")
        print(f"     genes: {', '.join(genes) or '-'}   consequence: {sev}")
    if not flagged:
        print("  none. Everything else scoring high is common population "
              "variation, so the regulatory door is shut.")


if __name__ == "__main__":
    main()
