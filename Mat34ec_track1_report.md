# Track 1 — Variant Prediction

**Rare Disease, Real Kid: The MVA Hackathon 2026**
Participant: `Mat34ec` · Proband `WGS_EX2312012` · GRCh38 · Submission `Mat34ec_track1_panel-vep.csv`

---

## 1. The prediction

**Primary — *BUB1B* compound heterozygote, `NM_001211.6`, EPCR 0.95**

| | Variant | Consequence |
|---|---|---|
| Allele 1 | `chr15:40209701 T>G` · `c.2210T>G` · `p.Leu737Ter` | Nonsense, NMD-competent |
| Allele 2 | `chr15:40220612 T>G` · `c.3006T>G` · `p.Asn1002Lys` | Missense, pseudokinase domain |

Mosaic variegated aneuploidy type 1 (MIM 257300), autosomal recessive, biallelic *BUB1B*.

**Secondary** — *LZTR1* `chr22:20996720 C>G` `p.Tyr748Ter`, monoallelic (EPCR 0.25); *SLC34A1*
`chr5:177386432 C>T` `p.Ala133Val`, heterozygous (EPCR 0.04, submitted as a **documented exclusion**,
§6.2); and a pharmacogenomic profile reported in §6.3.

---

## 2. Approach

A deliberately small pipeline, built so that someone else can re-run every step from this repository,
and so that no coordinate, threshold or gene list is a number typed in by hand.

```
  ①  panel construction      make_panel_bed.py       Ensembl REST, fetched live
  ②  region extraction       bcftools view -R
  ③  annotation              annotate_and_triage.py  Ensembl VEP REST
  ④  recessive triage        annotate_and_triage.py  both inheritance models
  ⑤  independent validation  §4, methods not sharing assumptions with ①–④
```

**① The panel is built, not pasted.** Twelve genes, tiered by how firmly each is established in MVA:
Tier 1, the numbered MVA types (*BUB1B, CEP57, TRIP13, CENATAC, SLF2, SMC5, MAD1L1*); Tier 2, reported
candidates (*MAD2L1BP, BUB1, CEP192*); Tier 3, exploratory (*TRIM37, BUB3*). ***MVK* is deliberately
excluded**: mevalonic aciduria shares the "MVA" abbreviation and is an unrelated disorder, and a
pipeline that includes it on string similarity has already made a category error.

Coordinates are fetched from Ensembl REST at run time, with alias fallback (*CENATAC*→*CCDC84*,
*SLF2*→*FAM178A*, *MAD2L1BP*→*CMT2*). The script emits GRCh37 and GRCh38, `chr`-prefixed and bare-contig
BEDs, and a provenance TSV recording Ensembl gene ID, assembly, padding and retrieval date per gene.
5 kb of padding captures promoter-proximal and UTR variants. Total footprint 1.3 Mb, 0.04% of the genome.

---

## 3. Three decisions that determined the answer

Stated explicitly because each could silently have produced the wrong result.

**`INFO/AF` in this VCF is not a population frequency.** It is GATK HaplotypeCaller output, where `AF` is
the **per-sample allele fraction**. Filtering "rare" variants on `INFO/AF < 0.01` discards essentially
every true heterozygous call, because a clean germline het sits near 0.5. Population frequency was taken
only from gnomAD via VEP, never from the VCF. **This is the single most consequential trap in the
dataset and it is invisible if the field name is trusted.**

**A variant absent from gnomAD is treated as AF 0.0, not as unknown.** Dropping unseen alleles discards
precisely the class of ultra-rare variant that causes a severe recessive disorder. `p.Leu737Ter` is
absent from gnomAD; a "require a frequency record" filter loses the primary finding outright.

**Both compound-heterozygous and homozygous models are scored.** A triage that only ranks homozygous
candidates cannot return a compound heterozygote, which is the answer here. The pipeline pairs same-gene
rare damaging alleles explicitly, as a design requirement rather than a post-hoc capability.

---

## 4. Independent validation

The triage produces a candidate. Everything below was run to break it, using methods that do not share
assumptions with the triage.

**Impact is attributed to the protein, not to regulation.** AlphaGenome Atlas AVI feature decomposition
for `p.Asn1002Lys`: **63% AlphaMissense, 19% PhastCons-470, 18% Cactus-241, 1.2% merged splicing, ~2.8%
across all nine regulatory tracks combined.** The model itself attributes the impact to protein coding.
This matters because AVI blends coding and regulatory signal, so a high AVI score alone would not have
separated the two.

**Splicing is excluded by two model families, with an internal positive control.** AlphaGenome ranked
`N1002K` **21st of 22** synthetic local controls (other substitutions at the same position, plus
transitions at ±5/10/20/30 bp), empirical *p* ≈ 0.95. SpliceAI v1.3.1 returned 0.02/0.00/0.00/0.00
against a 0.2 high-recall threshold — **in the same run** where a local control scored 1.00 donor gain
with 0.81 donor loss. A negative is only as good as the demonstration that the instrument registers a
positive.

**The locus is saturated, not spot-checked.** Atlas scored every possible SNV across *BUB1B* (192,459)
and *LZTR1* (63,174), plus 2,626 carried SNVs across both 1 Mb context windows. The two *BUB1B* variants
rank **1 and 2 of 2,626**. Thirteen indels at the loci are all intronic or 3'UTR.

**The parents are not consanguineous, which independently supports compound heterozygosity.** Excluding
chrX — which scores homozygous throughout in a male and is an artefact here — the longest autosomal run
of homozygosity is **under 2.6 Mb and no segment exceeds 5 Mb**. Recent consanguinity produces long
tracts. Unrelated parents are exactly what two different alleles arriving from two sides requires. This
argument needs no parental samples, and it also removes consanguinity as a rival explanation for the
family's recurrent pregnancy loss.

**No cryptic deletion masquerades as a second hit.** HaplotypeCaller does not call structural variants,
so a deletion in trans with a point mutation would be invisible to everything above. Twelve of thirteen
genes sit at 0.93–1.21× baseline depth with normal heterozygous fractions. The low variant density in
*BUB1B*, *BUB1*, *CEP57* and *CENATAC* is purifying selection on essential mitotic genes, not deletion —
their depth is normal. *TRIP13* is the exception and is reported as **unassessed** (§7.2).

**Both calls are clean germline heterozygotes.** `L737Ter` reads 21 reference / 25 alternate at depth 46;
`N1002K` reads 15 / 13 at depth 28. Neither is a low-fraction artefact, neither looks mosaic.

**The stop codon was confirmed against the reference, not inferred from the HGVS.** Reference codon
`TTA` (Leu) on the plus strand; `c.2210T>G` hits base 2 giving `TGA`, i.e. **UGA**, +4 context **A**.
Reconstructed from the Atlas per-position reference bases, independently of the annotation. This matters
downstream, because stop-codon identity determines whether readthrough is even plausible, and UGA is the
permissive one.

---

## 5. Tested and negative, with detection floors

A negative without a floor is not a result.

**Mosaic aneuploidy is below detection in this sample.** Scored per chromosome across 3.6M filtered
biallelic sites using depth and allele balance jointly. Chromosomes 20, 21 and 22 initially looked
positive: raised depth, raised heterozygous fraction, wider allele balance.

The discriminating test is the **shape** of the widening. A mosaic gain must widen the heterozygous
distribution **symmetrically**, because the extra copy duplicates either allele with equal probability,
giving modes at (1+f)/(2+f) and 1/(2+f). Mismapping is one-sided low, because mismapped reads only ever
add reference-matching reads.

The fraction of sites above 0.60 is flat at **0.095–0.126 on every chromosome**, including the suspicious
ones, while only the fraction below 0.40 moves. The excess is one-sided. Binning at 10 Mb localises it:
chr21's first bin runs at **2.34×** baseline, chr22's at 1.30×, chr20's mid-arm bin at 1.14×, all
acrocentric or pericentromeric, every other bin at 1.00–1.02×. Depth and low-tail correlate at **+0.63
across 286 bins** — one mechanism, not two. chrX returns a ratio of exactly **0.500**, a free positive
control proving the method resolves a real one-copy change.

**Floor ≈ 15% mosaic.** MVA aneuploidy is scored cytogenetically in cultured lymphocytes or fibroblasts;
this is uncultured blood, and aneuploid cells are plausibly selected against in rapidly dividing blood
lineages. **The negative is a statement about this sample, not about the child, and does not substitute
for a karyotype.**

**No subclonal LOH anywhere in the genome.** Tested as E[z²] by sign across 292 autosomal bins, where z
is the binomial deviate at heterozygous sites. Copy-number imbalance lifts both signs; mismapping lifts
only the negative. The positive side sits at the null (median 1.19, IQR 1.16–1.24); the negative side is
inflated genome-wide and one-sided. Every top-ranked bin is acrocentric or pericentromeric with strongly
negative mean z. **Floor 12–15% tumour fraction**, far above the circulating tumour DNA fraction expected
from a localized paediatric sarcoma, so the question of a somatic second hit is *untestable from blood*
rather than excluded.

---

## 6. Secondary and incidental findings

### 6.1 *LZTR1* `p.Tyr748Ter`, monoallelic, EPCR 0.25

A genuine inactivating allele in a gene with a clear genotype–phenotype split, submitted because the
proband's presentation included a malignancy. The codon is instructive: Tyr748 is `TAC`, and C>G makes
`TAG`, which is simultaneously a stop codon **and** an `AG` splice-acceptor consensus — which is why
SpliceAI returns DS_AG 0.92 here, and why an apparent disagreement between tools resolved with both
being right.

There is **no germline second hit** anywhere in the locus or the surrounding 1 Mb. Under two-hit logic a
monoallelic inactivating *LZTR1* allele requires a somatic second hit, which cannot be assessed without
tumour tissue (§5).

### 6.2 *SLC34A1* `p.Ala133Val`, EPCR 0.04 — a documented exclusion

Reported not as a candidate but as a closed negative, because the reasoning may be reusable. ClinVar is
**conflicting** (likely benign / likely pathogenic / uncertain significance); gnomAD max 0.0101.

The locus was swept on both sides. No second damaging coding allele exists. An AlphaGenome saturation
scan of all **94,839** possible SNVs across the gene placed the only rare upstream allele at the
**67th percentile**, rank 30,966 — below the gene's own 90th percentile and 35-fold below its maximum.
No phase tag exists anywhere in the gene, so trans configuration would be unprovable even for a positive.
A heterozygous, conflicting-classification allele in a recessive gene is a carrier finding.

### 6.3 Pharmacogenomic profile

Reported because the proband's syndrome carries a high tumour risk, which makes future systemic therapy
foreseeable, and because the VCF carries populated rsIDs so these loci can be read directly.

**Reference, therefore normal function:** all four *DPYD* risk alleles (`rs3918290`, `rs55886062`,
`rs67376798`, `rs56038477`); all three *TPMT* alleles and *NUDT15* `rs116855232`; *G6PD* (hemizygous in a
male, so a single variant would be fully penetrant); *SLCO1B1\*5*; *UGT1A1*; *RARG* and *SLC28A3*
(anthracycline cardiotoxicity); *ACYP2* (cisplatin ototoxicity). **Clean thiopurine and fluoropyrimidine
status in a child whose syndrome predisposes to leukaemia is usable information now.**

**Carried:** *CEP72* `rs924607` homozygous; *MTHFR* C677T homozygous; *CYP2B6* `*9` and the `*6` tag, both
heterozygous; *CYP2C19* `*2` and `*17`, both heterozygous; *GSTP1* I105V homozygous; *XPC* `rs2228001`
homozygous; the *ABCB1* 1236/2677/3435 haplotype heterozygous throughout; *SOD2* A16V, *CBR3*, *SLCO1B3*
`rs4149117` and *ABCC4* `rs11045585` heterozygous.

**Mitochondrial genome:** the aminoglycoside-ototoxicity alleles `m.1555A>G` and `m.1494C>T` are
**absent**; the only variant in that region is the common benign `m.1438A>G`, homoplasmic at depth 5,212.

**One result that reversed on checking, reported because of it.** *CEP72* `rs924607` came back homozygous
for the reported vincristine-neuropathy risk genotype, and CEP72 is a centrosomal protein, so the
convergence with a mitotic lesion looked striking. Checking the literature rather than asserting it
changed the reading twice. The association is contested: a Spanish B-ALL cohort found none during
induction, and an 817-patient meta-analysis could not confirm it globally, finding it only in the
continuation phase, while the largest analysis at 1,095 patients supports it at an odds ratio of 2.15.
More decisively, **the same literature names loci with larger effects and the proband is reference at all
of them** — *ABCC1* `rs3784867` (reported OR 4.9), *SLC5A7* `rs1013940` (8.6), *TTPA* `rs10504361` (2.0),
and the *CEP72* missense `rs12522955` (2.3). Of five reported vincristine-neuropathy loci he carries one,
and it is the weakest and most disputed. **A caution flag, not a contraindication**, and considerably
softer than the single-locus result looked on its own. Caveats in §7.4.

---

## 7. Limitations

1. **Phase is inferred, not observed.** No parental samples, and the variants are 10,911 bp apart while
   reads are 2×149 bp, so no read or read pair can span them. Read-backed phasing is impossible **in
   principle** here, not merely unavailable. Phase is inferred from the phenotype matching a recessive
   disorder. **This is the largest unverified assumption in the submission.**
2. ***TRIP13* is effectively unassessed.** Median depth 7 against a baseline of 44, three variant sites
   across 36.5 kb. A subtelomeric coverage failure rather than a deletion — heterozygous calls cannot
   exist in a homozygously deleted span — but a Tier 1 MVA gene has no usable data in this dataset.
3. **No structural variant calling.** Only SNV and indel calls are available. The depth-and-heterozygosity
   scan in §4 is a partial substitute, not a replacement for an SV caller on aligned reads.
4. **Pharmacogenomics is a single-variant screen, not star-allele diplotyping.** *CYP2D6* in particular
   cannot be genotyped reliably from short reads because of *CYP2D7* paralogy and copy-number variation;
   the `rs3892097` call shows an allele balance of 0.28, consistent with that paralogy, and no *CYP2D6*
   diplotype is asserted. Absence of a record means homozygous reference **only if** the position was
   covered, unconfirmable without a BAM. §6.3 is a screen, not a clinical pharmacogenomic report, and
   should be confirmed against current CPIC tables before use.
5. **The §5 negatives carry the floors stated there** and are not claims that the phenomena are absent in
   the child.

---

## 8. Methods description

*Responses to the provided Track 1 methods template.*

**Team name.** `Mat34ec` (individual participant).

**Model number.** 1 of 1.

**Describe your model/approach in detail.** §2–§5 above.

**Generative AI disclosure (required).** Anthropic Claude Pro and Claude Science, Claude Opus 5,
non-commercial terms, no training on customer content. Google Gemini 3.1 Pro, non-commercial terms, no
training on customer content. Perplexity (free tier) as a search engine, non-commercial, no training on
customer content. AI assistance was used for code authoring, literature retrieval, analytical design and
adversarial review of the author's own hypotheses. **Every quantitative claim in this report derives from
a script in this repository and can be re-run; no number is reported from model recall.** Claims drawn
from literature were checked against the primary source, and §9 records where that check overturned an
initial reading.

**Is the submission file automated output or manually curated?** The variant identities are the
**automated output** of the triage — the pipeline returns the *BUB1B* pair without human intervention.
Manual review was applied to the EPCR values and to the decision to include the two secondary rows.

**Describe any downstream manual review.** EPCR values were assigned by the author rather than computed:
the ordering is the pipeline's, the calibration is a judgement. The *SLC34A1* row was added by hand,
specifically to document a closed negative rather than to propose a candidate.

**Public or proprietary data?** Public only, exclusively. No proprietary sources.

**Describe the public data used.** Ensembl REST (gene coordinates, sequence, VEP consequences); gnomAD
allele frequencies via VEP; ClinVar via VEP; GENCODE (SpliceAI annotation); AlphaGenome and AlphaGenome
Atlas (AVI scores, feature decomposition, saturation scans); SpliceAI v1.3.1; dbSNP rsIDs as carried in
the VCF; PubMed for the literature in §6.3 and §9.

**Describe any proprietary data used.** None.

**Can the approach output compound heterozygous pairs?** Yes. Both compound-heterozygous and homozygous
models are scored and same-gene pairs are emitted explicitly. This is a design requirement (§3), not an
incidental capability: a pipeline without it cannot return the correct answer in this case.

**How were secondary or incidental findings handled?** Flagged via `finding_type=secondary` with rationale
in `notes`. Three classes are represented deliberately: a genuine candidate contributor (*LZTR1*), a
documented exclusion (*SLC34A1*), and a pharmacogenomic profile (§6.3) reported in the text rather than
as CSV rows, because it is not a causal-variant claim.

**Run time and cost.** Panel construction ~30 s, network-bound. Region extraction and triage ~5 min on a
laptop. The validation layers are the expensive part: Atlas saturation scans ~45 min across three genes,
the per-chromosome aneuploidy scan ~50 s, the genome-wide LOH scan ~60 s, SpliceAI ~3 min on CPU, the
panel CNV and pharmacogenomic scan ~17 s. No specialised hardware; a consumer laptop plus API access.

**Method abstract (≤500 words).** A twelve-gene MVA panel is constructed live from Ensembl REST with full
provenance rather than from hardcoded coordinates, and *MVK* is excluded on disease grounds despite the
shared abbreviation. Panel regions are extracted from the proband's VCF and annotated through Ensembl
VEP. Rare damaging alleles are triaged under both homozygous and compound-heterozygous models, with
population frequency drawn only from gnomAD — never from the VCF's `INFO/AF`, which in GATK output is
per-sample allele fraction and whose misuse would discard nearly every true heterozygous call — and with
gnomAD-absent variants treated as frequency zero rather than dropped. The pipeline returns a *BUB1B*
compound heterozygote: `p.Leu737Ter`, NMD-competent, with `p.Asn1002Lys` in the pseudokinase domain.
Five independent validation layers follow, each chosen not to share assumptions with the triage: AVI
feature decomposition attributing the missense impact to protein coding rather than regulation;
exclusion of splicing by AlphaGenome against 22 synthetic local controls and by SpliceAI with a positive
control in the same run; saturation scanning of every possible SNV across both loci plus 2,626 carried
SNVs in the surrounding 1 Mb, where the two *BUB1B* variants rank first and second; a runs-of-homozygosity
analysis showing no consanguinity, which supports compound heterozygosity without parental samples; and a
depth-and-heterozygosity scan excluding cryptic deletions across the panel. Two negatives are reported
with explicit detection floors. Mosaic aneuploidy is below roughly 15% in this uncultured blood sample,
established by a symmetry test that distinguishes true copy-number change from paralogous mismapping
using allele balance alone — a mosaic gain widens the heterozygous distribution symmetrically while
mismapping is one-sided low — with chrX as a built-in positive control at a depth ratio of exactly 0.500.
No subclonal loss of heterozygosity is detectable above 12–15% tumour fraction genome-wide. Incidental
findings include a pharmacogenomic profile showing normal-function genotypes at *DPYD*, *TPMT*, *NUDT15*,
*G6PD* and *SLCO1B1*. **Strengths:** full reproducibility from a live gene-coordinate source; validation
by orthogonal methods rather than repetition; negatives reported with detection floors; explicit
documentation of five instances where the pipeline's own first reading was wrong. **Limitations:** phase
is inferred rather than observed and cannot be observed with 2×149 bp reads across an 10,911 bp
separation; *TRIP13* is unassessed owing to subtelomeric coverage failure; no structural variant calling
is possible without aligned reads; and the pharmacogenomic screen is single-variant rather than
star-allele diplotyping.

---

## 9. Errors this pipeline caught in its own reasoning

Included because a method that never reports its own failures gives a reviewer no way to calibrate it.

- An apparent **11-fold splice-junction signal** at `L737Ter` looked real until synthetic local controls
  were run; a control 30 bp away scored *higher*, empirical *p* ≈ 0.24. The signal was an artefact of
  comparing two loci instead of comparing a variant to its own neighbourhood.
- AVI was initially judged to have a **gene-model error at *LZTR1***. SpliceAI later returned DS_AG 0.92
  at the same position, vindicating AVI. The codon explains both (§6.1).
- An apparent **mosaic aneuploidy signal** on chromosomes 20, 21 and 22 survived until the symmetry test
  was applied (§5).
- A locus-sweep script reported three "second candidate alleles" at *SLC34A1* that lie **5–6 kb beyond the
  gene's 3' end** and belong to a neighbouring transcript; the script had fallen back to
  `most_severe_consequence` when no *SLC34A1* transcript matched. The conclusion survived, but for the
  wrong reason until corrected.
- A single-locus pharmacogenomic lookup suggested a **vincristine-neuropathy risk genotype**. Running the
  full panel reversed the reading (§6.3).

---

## 10. Reproducing this

```bash
git clone https://github.com/Mat34ec/mva-hackathon-2026 && cd mva-hackathon-2026
pip install -r requirements.txt

# 1. build the panel from Ensembl — no hardcoded coordinates
python scripts/make_panel_bed.py --build GRCh38 --pad 5000 --out mva_panel

# 2. extract panel regions (data not included — see DATA.md)
bcftools view -R mva_panel.GRCh38.nochr.bed -Oz -o panel_raw.vcf.gz "$VCF"
tabix -f -p vcf panel_raw.vcf.gz

# 3. annotate and triage, both inheritance models
python scripts/annotate_and_triage.py --vcf panel_raw.vcf.gz --out candidates.tsv

# 4. validation layers (ALPHAGENOME_API_KEY lives in the environment, never in a file)
python scripts/atlas_percentile.py --gene BUB1B --mark "chr15:40220612:T>G=N1002K"
python scripts/atlas_decompose.py  --gene BUB1B
python scripts/aneu.py      1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y
python scripts/loh_scan.py  1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
python scripts/panel_cnv_pgx.py
```

**The patient data is not in this repository and must not be.** See [`DATA.md`](DATA.md) for the
controlled-access terms and the deletion obligation that comes with them.

---

*Released under CC BY 4.0 in accordance with the hackathon rules.*
