#!/usr/bin/env python3
"""
Build a synthetic VCF that exercises the copy-number discriminator, so the
method can be validated by anyone, with no access to the patient data and no
API key.

Why this exists. The central methodological claim of this submission is that
depth alone cannot distinguish a real mosaic chromosome gain from paralogous
mismapping, and that the SHAPE of the heterozygous allele-balance distribution
can. A reviewer has no quick way to check that claim against controlled-access
patient data, so this fixture reproduces both situations with known ground
truth and gives them the identical depth lift. If the discriminator works, it
must separate them anyway.

Four synthetic contigs, all with the same number of sites:

  chrNORM   clean diploid.          depth 1.00x, het AB ~ Binomial(DP, 0.5)
  chrGAIN   mosaic trisomy, f=0.30. depth 1.15x, het AB drawn from the two
                                    modes a real gain produces,
                                    (1+f)/(2+f)=0.565 and 1/(2+f)=0.435
  chrMAP    paralogous mismapping.  depth 1.15x  <-- deliberately identical
                                    to chrGAIN, but the extra reads all match
                                    the reference, so AB is pushed LOW only
  chrHEMI   hemizygous, male chrX.  depth 0.50x, almost no heterozygous calls

Expected discriminator behaviour, asserted in test_discriminator.py:

  chrGAIN   depth up AND both tails up        -> copy-number change
  chrMAP    depth up AND only the low tail up -> mismapping
  chrHEMI   depth ratio ~ 0.50                -> one-copy change detected

chrGAIN and chrMAP differ in no summary statistic except tail symmetry. That
is the whole point: a method that only looks at depth calls them the same.

Usage:
    python tests/make_synthetic_vcf.py --out tests/synthetic.vcf.gz
"""

import argparse
import os
import random

CONTIGS = {
    "chrNORM": 5_000_000,
    "chrGAIN": 5_000_000,
    "chrMAP": 5_000_000,
    "chrHEMI": 5_000_000,
}

HEADER = """##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##FILTER=<ID=LowQual,Description="Low quality">
##INFO=<ID=AF,Number=A,Type=Float,Description="Per-sample allele fraction (NOT population frequency)">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
"""

BASES = "ACGT"


def draw_depth(base_depth, ratio, rng):
    """Poisson-ish depth around base_depth * ratio."""
    mean = base_depth * ratio
    # sum of a few uniforms approximates the spread of real coverage well
    # enough for this purpose and avoids a scipy dependency
    d = int(round(mean + rng.gauss(0, mean ** 0.5)))
    return max(10, min(d, 300))


def binom(n, p, rng):
    return sum(1 for _ in range(n) if rng.random() < p)


def emit(fh, chrom, pos, ref, alt, gt, ref_n, alt_n, flt="PASS", gq=99):
    dp = ref_n + alt_n
    af = alt_n / dp if dp else 0.0
    fh.write(
        f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t500\t{flt}\tAF={af:.4f}\t"
        f"GT:AD:DP:GQ\t{gt}:{ref_n},{alt_n}:{dp}:{gq}\n"
    )


def build(path, n_sites=4000, base_depth=44, seed=20261024):
    """Write a BGZF-compressed, tabix-indexed VCF.

    Plain gzip will not do: tabix requires BGZF block compression, and pysam
    raises "seek not implemented in files compressed by method 1" on a plain
    gzip file. Written uncompressed first, then compressed and indexed with
    pysam's own helpers so the fixture needs no external bcftools/tabix.
    """
    import pysam
    rng = random.Random(seed)
    f = 0.30                       # mosaic fraction for chrGAIN
    hi = (1 + f) / (2 + f)         # 0.5652
    lo = 1 / (2 + f)               # 0.4348

    plain = path[:-3] if path.endswith(".gz") else path + ".vcf"
    with open(plain, "w") as fh:
        fh.write(HEADER)
        for c, ln in CONTIGS.items():
            fh.write(f"##contig=<ID={c},length={ln}>\n")
        fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSYNTH01\n")

        for chrom in CONTIGS:
            pos = 1000
            for _ in range(n_sites):
                pos += rng.randint(400, 1600)
                ref = rng.choice(BASES)
                alt = rng.choice([b for b in BASES if b != ref])

                if chrom == "chrNORM":
                    dp = draw_depth(base_depth, 1.00, rng)
                    if rng.random() < 0.62:
                        a = binom(dp, 0.5, rng)
                        emit(fh, chrom, pos, ref, alt, "0/1", dp - a, a)
                    else:
                        emit(fh, chrom, pos, ref, alt, "1/1", 0, dp)

                elif chrom == "chrGAIN":
                    # real mosaic trisomy: depth up, het AB splits into TWO
                    # modes, symmetric about 0.5
                    dp = draw_depth(base_depth, 1.15, rng)
                    if rng.random() < 0.62:
                        p = hi if rng.random() < 0.5 else lo
                        a = binom(dp, p, rng)
                        emit(fh, chrom, pos, ref, alt, "0/1", dp - a, a)
                    else:
                        emit(fh, chrom, pos, ref, alt, "1/1", 0, dp)

                elif chrom == "chrMAP":
                    # paralogous mismapping: SAME depth lift, but every extra
                    # read matches the reference, so AB only ever drops
                    dp_true = draw_depth(base_depth, 1.00, rng)
                    extra_ref = int(round(dp_true * 0.15))
                    if rng.random() < 0.62:
                        a = binom(dp_true, 0.5, rng)
                        emit(fh, chrom, pos, ref, alt, "0/1",
                             dp_true - a + extra_ref, a)
                    else:
                        emit(fh, chrom, pos, ref, alt, "1/1", extra_ref, dp_true)

                elif chrom == "chrHEMI":
                    # hemizygous: half depth, essentially no true hets
                    dp = draw_depth(base_depth, 0.50, rng)
                    if rng.random() < 0.04:
                        a = binom(dp, 0.5, rng)
                        emit(fh, chrom, pos, ref, alt, "0/1", dp - a, a)
                    else:
                        emit(fh, chrom, pos, ref, alt, "1/1", 0, dp)
    pysam.tabix_compress(plain, path, force=True)
    pysam.tabix_index(path, preset="vcf", force=True)
    os.remove(plain)
    return path


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="tests/synthetic.vcf.gz")
    ap.add_argument("--sites", type=int, default=4000)
    ap.add_argument("--depth", type=int, default=44)
    ap.add_argument("--seed", type=int, default=20261024)
    a = ap.parse_args()
    p = build(a.out, a.sites, a.depth, a.seed)
    print(f"wrote {p}: {len(CONTIGS)} contigs x {a.sites} sites at ~{a.depth}x")
    print("  chrNORM  diploid            expect ratio 1.00, symmetric tails")
    print("  chrGAIN  mosaic trisomy 30% expect ratio >1,   BOTH tails up")
    print("  chrMAP   mismapping         expect ratio >1,   LOW tail only")
    print("  chrHEMI  hemizygous         expect ratio ~0.50")
