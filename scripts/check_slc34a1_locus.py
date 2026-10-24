#!/usr/bin/env python3
"""
Close the SLC34A1 question: is there a SECOND damaging allele at this locus?

Context. The genome-wide rare-coding filter surfaced one SLC34A1 variant,
rs148976897 c.398C>T p.Ala133Val, heterozygous, with CONFLICTING ClinVar
classifications. SLC34A1 disease is recessive, so one allele is a carrier
state and nothing more. The only way it becomes an explanation for this
child's nephrocalcinosis is if a second damaging allele sits in trans and
the rare-coding filter missed it, which it would if that allele were
common, non-coding, or an indel.

So this script annotates EVERY carried allele across the locus, not just
the rare ones, and reports which are protein-altering or splice-relevant.

Phase note, decided before running: the VCF carries no GATK phase tag
anywhere in this locus, so even a second damaging allele could not be
shown to be in trans from this data alone. A positive result here is a
candidate, never a diagnosis.

Run from PowerShell in the research folder, where Ensembl is reachable:
    python check_slc34a1_locus.py
"""

import json
import os
import sys
import time
import urllib.error
import urllib.request

SERVER = "https://rest.ensembl.org"

# Consequences that could plausibly make an allele the second hit.
INTERESTING = {
    "transcript_ablation", "splice_acceptor_variant", "splice_donor_variant",
    "stop_gained", "frameshift_variant", "stop_lost", "start_lost",
    "transcript_amplification", "inframe_insertion", "inframe_deletion",
    "missense_variant", "protein_altering_variant",
    "splice_region_variant", "splice_donor_5th_base_variant",
    "splice_donor_region_variant", "splice_polypyrimidine_tract_variant",
}


def load_locus():
    here = os.path.dirname(os.path.abspath(__file__))
    path = os.path.join(here, "slc34a1_locus.json")
    if not os.path.exists(path):
        sys.exit(f"missing {path}")
    with open(path) as fh:
        return json.load(fh)


def to_region(r):
    """Ensembl region format. SNV is pos pos ref/alt; indels need the
    VCF-style anchor base stripped and the span widened accordingly."""
    pos, ref, alt = r["pos"], r["ref"], r["alt"]
    if len(ref) == 1 and len(alt) == 1:
        return f"5:{pos}-{pos}:1/{alt}", ref
    if len(ref) > len(alt) and ref.startswith(alt):          # deletion
        start = pos + len(alt)
        end = pos + len(ref) - 1
        return f"5:{start}-{end}:1/-", ref
    if len(alt) > len(ref) and alt.startswith(ref):          # insertion
        start = pos + len(ref)
        return f"5:{start}-{start - 1}:1/{alt[len(ref):]}", ref
    return f"5:{pos}-{pos + len(ref) - 1}:1/{alt}", ref      # complex


def vep(regions, retries=3):
    body = json.dumps({"variants": regions}).encode()
    url = (f"{SERVER}/vep/human/region"
           "?canonical=1&mane=1&hgvs=1&af_gnomadg=1&variant_class=1")
    req = urllib.request.Request(
        url, data=body,
        headers={"Content-Type": "application/json",
                 "Accept": "application/json"})
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(req, timeout=120) as resp:
                return json.loads(resp.read().decode())
        except urllib.error.HTTPError as e:
            if e.code in (429, 503) and attempt < retries - 1:
                time.sleep(3 * (attempt + 1))
                continue
            sys.exit(f"VEP HTTP {e.code}: {e.read().decode()[:300]}")
        except urllib.error.URLError as e:
            if attempt < retries - 1:
                time.sleep(3 * (attempt + 1))
                continue
            sys.exit(f"VEP unreachable: {e}")


def main():
    rows = load_locus()
    regions, meta = [], {}
    for r in rows:
        reg, ref = to_region(r)
        regions.append(reg)
        meta[reg] = r

    print(f"annotating {len(regions)} carried alleles across the SLC34A1 locus\n")
    out = vep(regions)

    by_input = {o.get("input"): o for o in out}
    flagged, clean = [], 0

    for reg in regions:
        o = by_input.get(reg)
        r = meta[reg]
        if o is None:
            print(f"  !! no VEP result for {reg}")
            continue

        # Prefer the MANE/canonical SLC34A1 consequence; fall back to worst.
        picks = [c for c in o.get("transcript_consequences", [])
                 if c.get("gene_symbol") == "SLC34A1"
                 and (c.get("mane_select") or c.get("canonical"))]
        if not picks:
            picks = [c for c in o.get("transcript_consequences", [])
                     if c.get("gene_symbol") == "SLC34A1"]
        terms = set()
        for c in picks:
            terms.update(c.get("consequence_terms", []))
        if not terms:
            terms = set(o.get("most_severe_consequence", "").split(","))

        af = o.get("colocated_variants", [{}])
        maxaf, rsid, clin = 0.0, "", ""
        for cv in af:
            if str(cv.get("id", "")).startswith("rs"):
                rsid = rsid or cv["id"]
            for k, val in (cv.get("frequencies") or {}).items():
                for pop, f in val.items():
                    try:
                        maxaf = max(maxaf, float(f))
                    except (TypeError, ValueError):
                        pass
            if cv.get("clin_sig"):
                clin = ",".join(cv["clin_sig"])

        hit = terms & INTERESTING
        tag = "  <<<" if hit else ""
        line = (f"  {r['pos']:>10} {r['ref'][:6]}>{r['alt'][:6]:<6} {r['zyg']:<4} "
                f"AF={maxaf:<8.4g} {rsid:<13} {','.join(sorted(terms))[:46]}")
        if clin:
            line += f" | ClinVar: {clin}"
        print(line + tag)

        if hit:
            p = picks[0] if picks else {}
            flagged.append((r, sorted(hit), maxaf, rsid, clin, p))
        else:
            clean += 1

    print(f"\n{clean} alleles are non-coding or synonymous; "
          f"{len(flagged)} are protein-altering or splice-relevant")

    others = [f for f in flagged if f[0]["pos"] != 177386432]
    print("\n" + "=" * 70)
    if not others:
        print("VERDICT  No second protein-altering or splice-relevant allele.")
        print("         p.Ala133Val stands alone: a heterozygous carrier state")
        print("         at a conflicting-classification allele, in a recessive")
        print("         gene. It does NOT explain the nephrocalcinosis.")
    else:
        print("VERDICT  A second candidate allele is present. Inspect below.")
        print("         Remember: no phase tags exist in this locus, so trans")
        print("         configuration cannot be shown without parental samples.")
        for r, hit, maxaf, rsid, clin, p in others:
            print(f"\n  {r['pos']} {r['ref']}>{r['alt']} {r['zyg']}  AD={r['ad']}")
            print(f"    {', '.join(hit)}   AF={maxaf:.4g}  {rsid}")
            if p.get("hgvsc"):
                print(f"    {p['hgvsc']}  {p.get('hgvsp', '')}")
            if p.get("sift_prediction"):
                print(f"    SIFT {p['sift_prediction']} ({p.get('sift_score')})")
            if p.get("polyphen_prediction"):
                print(f"    PolyPhen {p['polyphen_prediction']} ({p.get('polyphen_score')})")
            if clin:
                print(f"    ClinVar: {clin}")
            if maxaf > 0.01:
                print("    NOTE: too common to be a pathogenic recessive allele")
    print("=" * 70)


if __name__ == "__main__":
    main()
