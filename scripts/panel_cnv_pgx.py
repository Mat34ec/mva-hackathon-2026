#!/usr/bin/env python3
"""
Two checks the project has never run.

1. CRYPTIC DELETION SCAN across the 12 MVA panel genes plus LZTR1.
   The VCF is GATK HaplotypeCaller output, which calls SNVs and small
   indels and does NOT call structural variants. For a chromosomal
   instability syndrome that is a conspicuous gap: a heterozygous deletion
   removing part of a panel gene would be invisible to every analysis run
   so far, and would change the genetic model (a deletion in trans with a
   point mutation is still biallelic loss).
   Without a BAM we cannot call SVs directly, but a deletion leaves three
   fingerprints in a variants-only VCF:
     heterozygous deletion -> depth near half, and NO heterozygous calls
                              across the deleted span
     homozygous deletion   -> no calls of any kind across the span
   So we report depth ratio, heterozygous fraction and call density per
   gene against the genome-wide baseline.

2. PHARMACOGENOMICS by rsID.
   The mechanism map states plainly that metaboliser status is NOT
   established, because the earlier search covered rare coding variants
   and star alleles are mostly common. The VCF turns out to carry
   populated rsIDs, so the actionable loci can be read directly.
   This matters for Track 2: proposing a drug without knowing whether the
   child can clear it is half an answer. Vincristine, cyclophosphamide and
   actinomycin D are the agents actually used in the one published MVA
   case, so their pharmacogenes are the relevant ones.

   IMPORTANT on interpretation. This is a variants-only VCF, so a locus
   with no record means the child is homozygous reference there, i.e.
   normal-function. That is informative. The caveat is that a no-call from
   low coverage is indistinguishable from a reference call without a BAM;
   at 44x mean depth most positions are well covered, so absence is very
   likely reference, but it is not proof. Genotypes are reported as data.
   Star-allele diplotypes and CPIC phenotypes are NOT asserted here and
   must be confirmed against the current CPIC tables before use.
"""
import json
import glob
import os

import numpy as np
import pysam

# The dataset is controlled access and the proband accession is not
# published here. Set MVA_VCF to the file, or MVA_DATA to its folder.
D = os.environ.get("MVA_DATA",
                   os.environ.get("HOME", "") + "/mnt/MVA Hackathon 2026 - Dataset")
_hits = sorted(glob.glob(os.path.join(D, "*.vcf.gz")))
VCF = os.environ.get("MVA_VCF") or (_hits[0] if _hits else "")

GENES = [
    ("BUB1",     "2",  110630467, 110683098, "T2"),
    ("TRIP13",   "5",  887848,    924357,    "T1"),
    ("MAD2L1BP", "6",  43624493,  43645960,  "T2"),
    ("MAD1L1",   "7",  1810786,   2238243,   "T1"),
    ("SMC5",     "9",  70253269,  70359874,  "T1"),
    ("SLF2",     "10", 100907962, 100970134, "T1"),
    ("BUB3",     "10", 123149364, 123318144, "T3"),
    ("CEP57",    "11", 95784964,  95842070,  "T1"),
    ("CENATAC",  "11", 118993050, 119020811, "T1"),
    ("BUB1B",    "15", 40155983,  40226137,  "T1"),
    ("TRIM37",   "17", 58977632,  59111921,  "T3"),
    ("CEP192",   "18", 12986282,  13130053,  "T2"),
    ("LZTR1",    "22", 20979980,  21001037,  "--"),
]

# Clinically actionable pharmacogenes. Chosen around the agents actually
# used in paediatric soft tissue sarcoma and in the one published MVA
# chemotherapy case: vincristine, cyclophosphamide, actinomycin D, plus
# the anthracycline and platinum loci that a protocol change would raise.
PGX = {
    # thiopurines, relevant because MVA predisposes to leukaemia
    "rs1800462": "TPMT *2",
    "rs1800460": "TPMT *3B",
    "rs1142345": "TPMT *3C",
    "rs116855232": "NUDT15 *3",
    # fluoropyrimidines
    "rs3918290": "DPYD *2A",
    "rs55886062": "DPYD *13",
    "rs67376798": "DPYD p.Asp949Val",
    "rs56038477": "DPYD HapB3",
    # irinotecan
    "rs4148323": "UGT1A1 *6",
    "rs887829": "UGT1A1 *28 proxy",
    # methotrexate, used in sarcoma protocols
    "rs4149056": "SLCO1B1 *5",
    "rs1801133": "MTHFR C677T",
    "rs1801131": "MTHFR A1298C",
    # general oxidative metabolism
    "rs4244285": "CYP2C19 *2",
    "rs4986893": "CYP2C19 *3",
    "rs12248560": "CYP2C19 *17",
    "rs3892097": "CYP2D6 *4",
    "rs35742686": "CYP2D6 *3",
    "rs5030655": "CYP2D6 *6",
    "rs776746": "CYP3A5 *3",
    "rs1799853": "CYP2C9 *2",
    "rs1057910": "CYP2C9 *3",
    # cyclophosphamide activation
    "rs2279343": "CYP2B6 *6 (516/785)",
    "rs3745274": "CYP2B6 *9 G516T",
    # anthracycline cardiotoxicity
    "rs2229774": "RARG, anthracycline cardiotox",
    "rs7853758": "SLC28A3, anthracycline cardiotox",
    "rs1056892": "CBR3, anthracycline cardiotox",
    # platinum ototoxicity, relevant alongside the nephrocalcinosis
    "rs1872328": "ACYP2, cisplatin ototox",
    # vincristine neuropathy; CEP72 is itself a centrosomal protein
    "rs924607": "CEP72, vincristine neuropathy",
    # G6PD, hemizygous in a male
    "rs1050828": "G6PD A-",
    "rs1050829": "G6PD A",
}


def gene_stats(v, s, chrom, a, b):
    dps, n_het, n_hom, n = [], 0, 0, 0
    for rec in v.fetch(chrom, a, b):
        if rec.filter.keys() and "PASS" not in rec.filter.keys():
            continue
        if not rec.alts or len(rec.alts) != 1:
            continue
        if len(rec.ref) != 1 or len(rec.alts[0]) != 1:
            continue
        g = rec.samples[s]
        gt = g.get("GT")
        if gt is None or None in gt or all(x == 0 for x in gt):
            continue
        ad = g.get("AD")
        if not ad or len(ad) < 2:
            continue
        dp = int(ad[0]) + int(ad[1])
        if dp < 5:
            continue
        n += 1
        dps.append(dp)
        if len(set(gt)) > 1:
            n_het += 1
        else:
            n_hom += 1
    return dps, n_het, n_hom, n


if __name__ == "__main__":
    v = pysam.VariantFile(VCF)
    s = list(v.header.samples)[0]

    print("=" * 74)
    print("1. CRYPTIC DELETION SCAN, 12 MVA panel genes plus LZTR1")
    print("=" * 74)
    print("Baseline from the aneuploidy scan: autosomal median depth 44x,")
    print("heterozygous fraction about 0.62 of carried sites.\n")
    print(f"{'gene':<9} {'tier':<5} {'kb':>6} {'sites':>6} {'per kb':>7} "
          f"{'medDP':>6} {'ratio':>6} {'het':>5} {'het_frac':>8}  flag")
    out = {}
    for name, c, a, b in [(g[0], g[1], g[2], g[3]) for g in GENES]:
        tier = dict((g[0], g[4]) for g in GENES)[name]
        dps, nh, nhom, n = gene_stats(v, s, c, a, b)
        kb = (b - a) / 1000
        if n == 0:
            print(f"{name:<9} {tier:<5} {kb:>6.1f} {0:>6} {'--':>7} "
                  f"{'--':>6} {'--':>6} {'--':>5} {'--':>8}  "
                  f"NO CALLS -> possible homozygous deletion")
            continue
        med = float(np.median(dps))
        ratio = med / 44.0
        hf = nh / n
        flag = ""
        if ratio < 0.75:
            flag = "DEPTH LOW -> possible het deletion"
        if hf < 0.25 and n >= 20:
            flag = (flag + " ").strip() + " NO HETS -> possible het deletion"
        if ratio > 1.3:
            flag = "depth high -> duplication or paralogy"
        print(f"{name:<9} {tier:<5} {kb:>6.1f} {n:>6} {n/kb:>7.2f} "
              f"{med:>6.1f} {ratio:>6.3f} {nh:>5} {hf:>8.3f}  {flag}")
        out[name] = {"n": n, "per_kb": n / kb, "med_dp": med,
                     "ratio": ratio, "n_het": nh, "het_frac": hf}

    print("\n" + "=" * 74)
    print("2. PHARMACOGENOMICS by rsID")
    print("=" * 74)
    found = {}
    want = set(PGX)
    for rec in v.fetch():
        rid = rec.id
        if not rid or rid == ".":
            continue
        for part in rid.split(";"):
            if part in want:
                g = rec.samples[s]
                gt = g.get("GT")
                ad = g.get("AD") or []
                zyg = ("het" if gt and len(set(gt)) > 1
                       else "hom_alt" if gt and gt[0] not in (None, 0)
                       else "ref")
                found[part] = {
                    "label": PGX[part], "chrom": rec.chrom, "pos": rec.pos,
                    "ref": rec.ref, "alt": ",".join(rec.alts or []),
                    "zyg": zyg, "gt": str(gt),
                    "ad": ",".join(str(x) for x in ad),
                    "filter": ";".join(rec.filter.keys()) or "PASS",
                    "gq": g.get("GQ"),
                }
    print(f"\nCARRIED ({len(found)} of {len(PGX)} loci queried):")
    if not found:
        print("  none")
    for r, d_ in sorted(found.items(), key=lambda kv: kv[1]["label"]):
        print(f"  {d_['label']:<34} {r:<12} {d_['chrom']}:{d_['pos']} "
              f"{d_['ref']}>{d_['alt']:<4} {d_['zyg']:<8} "
              f"AD {d_['ad']:<9} GQ {d_['gq']} {d_['filter']}")
    missing = sorted(set(PGX) - set(found), key=lambda r: PGX[r])
    print(f"\nNO RECORD, therefore homozygous reference / normal function "
          f"({len(missing)}):")
    for r in missing:
        print(f"  {PGX[r]:<34} {r}")

    json.dump({"genes": out, "pgx_found": found,
               "pgx_absent": [PGX[r] for r in missing]},
              open(os.environ["HOME"] + "/an/panel_pgx.json", "w"), indent=1)
