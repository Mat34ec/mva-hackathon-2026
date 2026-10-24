#!/usr/bin/env python3
"""
Classify the one candidate for the nephrocalcinosis.

Nephrocalcinosis is the single phenotype feature neither BUB1B nor LZTR1
explains. Across 32 renal stone and tubular genes the patient carries one
rare coding variant: SLC34A1 chr5:177386432 C>T, p.Ala133Val, het,
gnomAD 0.0041. Biallelic SLC34A1 loss causes infantile hypercalcaemia
type 2, which presents with nephrocalcinosis.

A single heterozygote in a recessive condition is a carrier finding, so
this only changes the phenotype table if ClinVar classifies the specific
allele as pathogenic rather than uncertain. Otherwise the row stays with
prematurity and NICU diuretics as the likelier explanation.
"""

import json, urllib.request

REST = "https://rest.ensembl.org"
VARIANTS = [
    ("5 177386432 . C T . . .", "SLC34A1 p.Ala133Val"),
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
    for rec in vep([v for v, _ in VARIANTS]):
        print(f"=== {labels.get(rec.get('input',''), rec.get('input'))} ===")
        print(f"most_severe_consequence: {rec.get('most_severe_consequence')}")

        sigs, rsids, afs = set(), set(), {}
        for cv in rec.get("colocated_variants") or []:
            for s in cv.get("clin_sig") or []:
                sigs.add(s)
            if str(cv.get("id", "")).startswith("rs"):
                rsids.add(cv["id"])
            for allele, pops in (cv.get("frequencies") or {}).items():
                for p, v in pops.items():
                    try:
                        afs[p] = float(v)
                    except (TypeError, ValueError):
                        pass

        print(f"rsID            : {', '.join(sorted(rsids)) or 'none'}")
        print(f"ClinVar clin_sig: {', '.join(sorted(sigs)) or 'NONE REPORTED'}")
        if afs:
            top = sorted(afs.items(), key=lambda kv: -kv[1])[:6]
            print("population AF   : " + ", ".join(f"{k} {v:.5g}" for k, v in top))

        for t in rec.get("transcript_consequences", []) or []:
            if t.get("gene_symbol") != "SLC34A1":
                continue
            if not (t.get("mane_select") or t.get("canonical")):
                continue
            print(f"\n  {t.get('transcript_id')}  MANE={t.get('mane_select') or '-'}")
            print(f"    consequence : {','.join(t.get('consequence_terms', []))}")
            print(f"    {t.get('hgvsc','-')}  {t.get('hgvsp','-')}")
            print(f"    exon        : {t.get('exon','-')}")
            for k in ("sift_prediction", "polyphen_prediction"):
                if t.get(k):
                    print(f"    {k}: {t[k]} ({t.get(k.replace('prediction','score'),'')})")

        print("\nHow to read this. Pathogenic or Likely pathogenic would make")
        print("SLC34A1 a named candidate contributor to the nephrocalcinosis,")
        print("still as a carrier state in a recessive condition. Uncertain")
        print("significance, benign, or nothing reported leaves prematurity and")
        print("NICU diuretic exposure as the likelier explanation, and the")
        print("phenotype table stays as written.")


if __name__ == "__main__":
    main()
