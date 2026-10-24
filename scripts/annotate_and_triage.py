#!/usr/bin/env python3
"""
Annotate panel variants via Ensembl VEP REST, then triage for a recessive
two-hit pattern.

Input:  panel_variants.tsv  (from bcftools query, see README block below)
Output: panel_annotated.tsv  every variant with consequence + population AF
        candidates.tsv       genes showing a plausible recessive pattern

Design notes that matter for correctness:
  * INFO/AF in the source VCF is per-sample allele fraction, NOT population
    frequency. Population AF here comes only from VEP's colocated_variants.
  * A variant ABSENT from gnomAD is treated as AF 0.0, not as unknown.
    Dropping unseen variants would discard exactly the class of ultra-rare
    allele that causes a severe recessive disorder.
  * Both compound-heterozygous and homozygous models are scored. Homozygous
    is common in consanguineous MVA pedigrees, so assuming compound het
    a priori would bias the search.

Usage:
    python3 annotate_and_triage.py
    python3 annotate_and_triage.py --af-max 0.001 --keep-nonpass
    python3 annotate_and_triage.py --minimal      # if rich plugins 400
"""

import argparse
import json
import sys
import time
import urllib.error
import urllib.request
from collections import defaultdict

REST = "https://rest.ensembl.org"
BATCH = 200

PANEL_TIER = {
    "BUB1B": "T1", "CEP57": "T1", "TRIP13": "T1", "CENATAC": "T1",
    "CCDC84": "T1", "SLF2": "T1", "FAM178A": "T1", "SMC5": "T1",
    "MAD1L1": "T1",
    "MAD2L1BP": "T2", "BUB1": "T2", "CEP192": "T2",
    "TRIM37": "T3", "BUB3": "T3",
}

TRUNCATING = {
    "transcript_ablation", "stop_gained", "frameshift_variant",
    "splice_acceptor_variant", "splice_donor_variant", "start_lost",
}
MISSENSE_LIKE = {
    "missense_variant", "inframe_deletion", "inframe_insertion",
    "protein_altering_variant", "stop_lost",
}
SPLICE_REGION = {
    "splice_region_variant", "splice_donor_5th_base_variant",
    "splice_polypyrimidine_tract_variant", "splice_donor_region_variant",
}

# rough severity ordering, most severe first
SEVERITY = (list(TRUNCATING) + list(MISSENSE_LIKE) + list(SPLICE_REGION) +
            ["synonymous_variant", "5_prime_UTR_variant", "3_prime_UTR_variant",
             "intron_variant", "upstream_gene_variant", "downstream_gene_variant",
             "intergenic_variant"])
SEV_RANK = {c: i for i, c in enumerate(SEVERITY)}


def severity_rank(consequences):
    return min((SEV_RANK.get(c, 999) for c in consequences), default=999)


def vep_post(variants, rich=True, retries=4):
    """POST a batch of 'CHROM POS . REF ALT . . .' strings to VEP REST."""
    params = {
        "hgvs": 1, "canonical": 1, "mane": 1, "numbers": 1,
        "variant_class": 1, "domains": 1,
    }
    if rich:
        params.update({"AlphaMissense": 1, "CADD": 1, "LOEUF": 1})

    qs = "&".join(f"{k}={v}" for k, v in params.items())
    url = f"{REST}/vep/human/region?{qs}"
    body = json.dumps({"variants": variants}).encode()
    req = urllib.request.Request(
        url, data=body,
        headers={"Content-Type": "application/json",
                 "Accept": "application/json"},
    )

    for attempt in range(retries):
        try:
            with urllib.request.urlopen(req, timeout=120) as resp:
                return json.loads(resp.read().decode())
        except urllib.error.HTTPError as e:
            if e.code == 400 and rich:
                print("  rich plugin options rejected, falling back to minimal",
                      file=sys.stderr)
                return vep_post(variants, rich=False, retries=retries)
            if e.code in (429, 503, 504):
                wait = 3 * (attempt + 1)
                print(f"  server busy ({e.code}), retrying in {wait}s",
                      file=sys.stderr)
                time.sleep(wait)
                continue
            if attempt == retries - 1:
                raise
            time.sleep(2 * (attempt + 1))
        except urllib.error.URLError:
            if attempt == retries - 1:
                raise
            time.sleep(3 * (attempt + 1))
    return []


def max_gnomad_af(colocated, alt):
    """Highest gnomAD population AF across any subpopulation.

    Returns (af, source_label). Absence from gnomAD returns (0.0, 'absent'),
    which is deliberate: unseen means ultra-rare, not unknown.
    """
    best, label, seen_any = 0.0, "absent", False
    for cv in colocated or []:
        freqs = cv.get("frequencies") or {}
        for allele, pops in freqs.items():
            if allele != alt:
                continue
            for pop, val in pops.items():
                if not pop.lower().startswith("gnomad"):
                    continue
                try:
                    v = float(val)
                except (TypeError, ValueError):
                    continue
                seen_any = True
                if v > best:
                    best, label = v, pop
    if not seen_any:
        return 0.0, "absent"
    return best, label


def clin_sig(colocated):
    sigs = set()
    for cv in colocated or []:
        for s in cv.get("clin_sig") or []:
            sigs.add(s)
    return ";".join(sorted(sigs)) if sigs else ""


def rsid(colocated):
    for cv in colocated or []:
        i = cv.get("id", "")
        if i.startswith("rs"):
            return i
    return ""


def pick_consequence(rec):
    """Best transcript consequence restricted to panel genes.

    Prefers MANE Select, then Ensembl canonical, then most severe.
    """
    tcs = [t for t in rec.get("transcript_consequences", [])
           if t.get("gene_symbol") in PANEL_TIER]
    if not tcs:
        return None

    def key(t):
        return (
            0 if t.get("mane_select") else 1,
            0 if t.get("canonical") else 1,
            severity_rank(t.get("consequence_terms", [])),
        )

    return sorted(tcs, key=key)[0]


def classify(terms):
    s = set(terms)
    if s & TRUNCATING:
        return "truncating"
    if s & MISSENSE_LIKE:
        return "missense"
    if s & SPLICE_REGION:
        return "splice_region"
    if "synonymous_variant" in s:
        return "synonymous"
    if s & {"5_prime_UTR_variant", "3_prime_UTR_variant"}:
        return "utr"
    if "intron_variant" in s:
        return "intronic"
    return "other"


def load_variants(path, keep_nonpass):
    rows, skipped = [], defaultdict(int)
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 12:
                skipped["malformed"] += 1
                continue
            chrom, pos, ref, alt, filt, qual, gt, ad, dp, gq, pgt, pid = f[:12]

            if alt in (".", "*"):
                skipped["no_alt"] += 1
                continue
            if gt in ("0/0", "0|0", "./.", ".|.", "."):
                skipped["not_carried"] += 1
                continue
            if not keep_nonpass and filt not in ("PASS", "."):
                skipped[f"filtered:{filt}"] += 1
                continue

            rows.append({
                "chrom": chrom, "pos": int(pos), "ref": ref, "alt": alt,
                "filter": filt, "qual": qual, "gt": gt, "ad": ad,
                "dp": dp, "gq": gq, "pgt": pgt, "pid": pid,
            })
    return rows, skipped


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", default="panel_variants.tsv")
    ap.add_argument("--af-max", type=float, default=0.005,
                    help="max gnomAD popmax AF to call a variant rare "
                         "(default 0.005; severe recessive disease)")
    ap.add_argument("--keep-nonpass", action="store_true",
                    help="keep variants failing the GATK FILTER column")
    ap.add_argument("--minimal", action="store_true",
                    help="skip AlphaMissense/CADD/LOEUF plugin options")
    args = ap.parse_args()

    rows, skipped = load_variants(args.input, args.keep_nonpass)
    print(f"loaded {len(rows)} carried variants from {args.input}")
    for k, v in sorted(skipped.items(), key=lambda x: -x[1]):
        print(f"  skipped {v:>7}  {k}")
    if not rows:
        sys.exit("nothing to annotate, check the extraction step")

    # VEP join key must match the 'input' echo exactly
    for r in rows:
        r["vep_in"] = f"{r['chrom']} {r['pos']} . {r['ref']} {r['alt']} . . ."

    annotations = {}
    batches = [rows[i:i + BATCH] for i in range(0, len(rows), BATCH)]
    print(f"annotating in {len(batches)} VEP batches")
    for i, batch in enumerate(batches, 1):
        results = vep_post([r["vep_in"] for r in batch], rich=not args.minimal)
        for rec in results:
            annotations[rec.get("input", "")] = rec
        print(f"  batch {i}/{len(batches)}  ({len(results)} returned)")
        time.sleep(0.4)   # stay under the published rate limit

    out_rows = []
    for r in rows:
        rec = annotations.get(r["vep_in"])
        if rec is None:
            r.update(gene="", tier="", consequence="unannotated", klass="unknown",
                     hgvsc="", hgvsp="", af=0.0, af_pop="no_vep_result",
                     clinvar="", rs="", am="", cadd="", loeuf="", mane="")
            out_rows.append(r)
            continue

        tc = pick_consequence(rec)
        col = rec.get("colocated_variants")
        af, af_pop = max_gnomad_af(col, r["alt"])

        if tc is None:
            r.update(gene="", tier="", consequence=rec.get("most_severe_consequence", ""),
                     klass="off_panel_transcript", hgvsc="", hgvsp="", mane="")
        else:
            terms = tc.get("consequence_terms", [])
            gene = tc.get("gene_symbol", "")
            r.update(
                gene=gene,
                tier=PANEL_TIER.get(gene, ""),
                consequence=",".join(terms),
                klass=classify(terms),
                hgvsc=tc.get("hgvsc", ""),
                hgvsp=tc.get("hgvsp", ""),
                mane=tc.get("mane_select", "") or "",
                am=tc.get("alphamissense", {}).get("am_pathogenicity", "")
                   if isinstance(tc.get("alphamissense"), dict) else "",
                cadd=tc.get("cadd_phred", ""),
                loeuf=tc.get("loeuf", ""),
            )
        r.update(af=af, af_pop=af_pop, clinvar=clin_sig(col), rs=rsid(col))
        out_rows.append(r)

    cols = ["chrom", "pos", "ref", "alt", "gene", "tier", "klass", "consequence",
            "hgvsc", "hgvsp", "gt", "af", "af_pop", "clinvar", "rs",
            "cadd", "am", "loeuf", "filter", "qual", "ad", "dp", "gq",
            "pgt", "pid", "mane"]
    with open("panel_annotated.tsv", "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in out_rows:
            fh.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")
    print(f"\nwrote panel_annotated.tsv ({len(out_rows)} variants)")

    # ---- triage -------------------------------------------------------
    DAMAGING = {"truncating", "missense", "splice_region"}
    by_gene = defaultdict(list)
    for r in out_rows:
        if not r["gene"]:
            continue
        if r["klass"] not in DAMAGING:
            continue
        if r["af"] > args.af_max:
            continue
        by_gene[r["gene"]].append(r)

    def is_hom(gt):
        a = gt.replace("|", "/").split("/")
        return len(a) == 2 and a[0] == a[1] and a[0] not in ("0", ".")

    candidates = []
    for gene, vs in by_gene.items():
        homs = [v for v in vs if is_hom(v["gt"])]
        hets = [v for v in vs if not is_hom(v["gt"])]
        trunc = [v for v in vs if v["klass"] == "truncating"]

        model, score = None, 0
        if homs:
            model = "homozygous"
            score = 100 + (20 if any(v["klass"] == "truncating" for v in homs) else 0)
        elif len(hets) >= 2:
            model = "possible compound het"
            score = 80 + (20 if trunc else 0) + min(len(hets), 4)
        elif len(hets) == 1:
            model = "single het only"
            score = 20 + (10 if trunc else 0)

        if model:
            tier_bonus = {"T1": 15, "T2": 8, "T3": 3}.get(vs[0]["tier"], 0)
            candidates.append((score + tier_bonus, gene, model, vs))

    candidates.sort(reverse=True, key=lambda x: x[0])

    with open("candidates.tsv", "w") as fh:
        fh.write("score\tgene\ttier\tmodel\tn_rare_damaging\tdetail\n")
        for score, gene, model, vs in candidates:
            detail = "; ".join(
                f"{v['chrom']}:{v['pos']}{v['ref']}>{v['alt']} {v['gt']} "
                f"{v['klass']} {v['hgvsp'] or v['hgvsc']} AF={v['af']:g}"
                for v in sorted(vs, key=lambda x: x["pos"])
            )
            fh.write(f"{score}\t{gene}\t{vs[0]['tier']}\t{model}\t{len(vs)}\t{detail}\n")

    print("wrote candidates.tsv\n")
    print("=== triage ===")
    if not candidates:
        print("NO candidate gene met the recessive pattern.")
        print("Do not stop here. Next moves, in order:")
        print("  1. rerun with --keep-nonpass (the causal variant may be FILTER-flagged)")
        print("  2. relax --af-max to 0.02")
        print("  3. inspect panel_annotated.tsv for intronic/UTR variants by hand")
        print("  4. escalate to the genome-wide backstop: the panel may be wrong")
    else:
        for score, gene, model, vs in candidates[:8]:
            print(f"  {score:>4}  {gene:<9} {vs[0]['tier']}  {model:<22} "
                  f"{len(vs)} rare damaging allele(s)")
            for v in sorted(vs, key=lambda x: x["pos"]):
                phase = f" PID={v['pid']}" if v["pid"] not in ("", ".") else ""
                print(f"          {v['chrom']}:{v['pos']} {v['ref']}>{v['alt']} "
                      f"{v['gt']} {v['klass']:<13} {v['hgvsp'] or v['hgvsc']} "
                      f"AF={v['af']:g} ({v['af_pop']}){phase}")

    print("\nReminder: a shared non-empty PID on two variants means GATK phased "
          "them physically. Different or empty PID means phase is UNKNOWN, not trans.")


if __name__ == "__main__":
    main()
