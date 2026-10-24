#!/usr/bin/env python3
"""
Annotate the indels at both loci.

The AlphaGenome saturation scans cover single-nucleotide changes only, so
every insertion and deletion the patient carries across BUB1B and LZTR1
is currently unassessed. That matters most for two adjacent heterozygous
deletions in LZTR1 at 20989703 and 20989712: if either is exonic and
frameshifting, and lies on the allele opposite p.Tyr748Ter, the patient
has biallelic LZTR1 rather than haploinsufficiency, which is a different
diagnosis with a different inheritance pattern.

Most are probably intronic. This checks rather than assumes.
"""

import json, urllib.request

REST = "https://rest.ensembl.org"

VARIANTS = [
    ("15 40180642 . CT C . . .",   "BUB1B del, LowQual"),
    ("15 40181093 . AT A . . .",   "BUB1B del"),
    ("15 40198239 . G GT . . .",   "BUB1B ins, hom"),
    ("15 40216568 . TA T . . .",   "BUB1B del, hom"),
    ("15 40216570 . TA T . . .",   "BUB1B del, hom"),
    ("22 20984622 . G GGT . . .",  "LZTR1 ins, hom"),
    ("22 20985059 . A AT . . .",   "LZTR1 ins, hom"),
    ("22 20986376 . T TGATA . . .", "LZTR1 ins, hom"),
    ("22 20986456 . T TGATA . . .", "LZTR1 ins, het"),
    ("22 20987388 . C CA . . .",   "LZTR1 ins, het"),
    ("22 20989703 . GGGCGCA G . . .", "LZTR1 6bp del, het  <<<"),
    ("22 20989712 . TAGAGGAGGTGAGGGGCACGGGGAGCCAGGGC T . . .",
     "LZTR1 31bp del, het  <<<"),
    ("22 20998799 . G GGCCAGGTGC . . .", "LZTR1 9bp ins, het"),
]


def vep(lines):
    url = f"{REST}/vep/human/region?hgvs=1&canonical=1&mane=1&numbers=1"
    req = urllib.request.Request(
        url, data=json.dumps({"variants": lines}).encode(),
        headers={"Content-Type": "application/json", "Accept": "application/json"})
    with urllib.request.urlopen(req, timeout=120) as r:
        return json.loads(r.read().decode())


def main():
    labels = dict(VARIANTS)
    res = vep([v for v, _ in VARIANTS])
    by_input = {r.get("input", ""): r for r in res}

    print(f"{'variant':<46} {'severe consequence':<26} MANE transcript detail")
    print("-" * 118)
    for spec, note in VARIANTS:
        rec = by_input.get(spec)
        if rec is None:
            print(f"{note:<46} NO RESULT")
            continue
        sev = rec.get("most_severe_consequence", "-")
        mane = None
        for t in rec.get("transcript_consequences", []) or []:
            if t.get("gene_symbol") in ("BUB1B", "LZTR1") and t.get("mane_select"):
                mane = t
                break
        if mane is None:
            for t in rec.get("transcript_consequences", []) or []:
                if t.get("gene_symbol") in ("BUB1B", "LZTR1") and t.get("canonical"):
                    mane = t
                    break
        if mane is None:
            print(f"{note:<46} {sev:<26} (no BUB1B/LZTR1 transcript)")
            continue
        terms = ",".join(mane.get("consequence_terms", []))
        where = f"exon {mane.get('exon')}" if mane.get("exon") else \
                (f"intron {mane.get('intron')}" if mane.get("intron") else "-")
        print(f"{note:<46} {sev:<26} {terms} | {where} | "
              f"{mane.get('hgvsc','-')} {mane.get('hgvsp','') or ''}")

    print("\nWhat to look for: any 'frameshift_variant', 'stop_gained', "
          "'splice_acceptor_variant' or 'splice_donor_variant' on a MANE "
          "transcript. Anything reported as intron_variant is not a second hit.")


if __name__ == "__main__":
    main()
