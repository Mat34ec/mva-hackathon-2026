#!/usr/bin/env python3
"""
Resolve the LZTR1 transcript discrepancy.

The pipeline annotated chr22:20996720 C>G as ENST00000646124.2
c.2244C>G p.Tyr748Ter, a stop_gained on the MANE transcript. AlphaGenome
Atlas returns PROTEIN_TERMINATION = 0 for the same variant and attributes
45% of its impact to splicing instead. Both cannot be right, and the
LZTR1 arm of the two-gene model depends on which.

This asks Ensembl VEP directly and prints every LZTR1 transcript
consequence, so the disagreement is resolved against the primary source
rather than by preferring whichever tool we like better.
"""

import json, urllib.request

REST = "https://rest.ensembl.org"
VARIANTS = [
    ("22 20996720 . C G . . .", "LZTR1 reported as p.Tyr748Ter"),
    ("15 40209701 . T G . . .", "BUB1B L737Ter, positive control"),
    ("15 40220612 . T G . . .", "BUB1B N1002K, missense control"),
]


def vep(lines):
    url = f"{REST}/vep/human/region?hgvs=1&canonical=1&mane=1&numbers=1"
    req = urllib.request.Request(
        url, data=json.dumps({"variants": lines}).encode(),
        headers={"Content-Type": "application/json", "Accept": "application/json"})
    with urllib.request.urlopen(req, timeout=120) as r:
        return json.loads(r.read().decode())


def main():
    res = vep([v for v, _ in VARIANTS])
    labels = dict(VARIANTS)
    for rec in res:
        inp = rec.get("input", "")
        print(f"\n================ {labels.get(inp, inp)} ================")
        print(f"input: {inp}")
        print(f"most_severe_consequence: {rec.get('most_severe_consequence')}")
        tcs = rec.get("transcript_consequences", []) or []
        genes = sorted({t.get("gene_symbol", "") for t in tcs})
        print(f"{len(tcs)} transcript consequences across genes: {genes}")

        def key(t):
            return (0 if t.get("mane_select") else 1,
                    0 if t.get("canonical") else 1,
                    t.get("transcript_id", ""))

        shown = 0
        for t in sorted(tcs, key=key):
            terms = t.get("consequence_terms", [])
            interesting = t.get("mane_select") or t.get("canonical") or \
                any(x in terms for x in ("stop_gained", "splice_donor_variant",
                                         "splice_acceptor_variant", "missense_variant"))
            if not interesting or shown >= 10:
                continue
            shown += 1
            print(f"\n  {t.get('gene_symbol')}  {t.get('transcript_id')}"
                  f"  biotype={t.get('biotype')}")
            print(f"    MANE_Select : {t.get('mane_select') or '-'}")
            print(f"    canonical   : {bool(t.get('canonical'))}")
            print(f"    consequence : {','.join(terms)}")
            print(f"    hgvsc       : {t.get('hgvsc','-')}")
            print(f"    hgvsp       : {t.get('hgvsp','-')}")
            print(f"    exon/intron : {t.get('exon','-')} / {t.get('intron','-')}")
            if t.get("protein_start"):
                print(f"    protein pos : {t.get('protein_start')}"
                      f"  aa {t.get('amino_acids','-')}")


if __name__ == "__main__":
    main()
