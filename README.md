# MVA Hackathon 2026, Tracks 1 and 2

Submission for **Rare Disease, Real Kid: The MVA Hackathon 2026**
(Sage Bionetworks · MVA Society · Hugging Face).

Participant HF **`Mat34ec`** · GRCh38

> **No patient data is in this repository and none may be added.** The genome and clinical record
> belong to a real, living child, released under a controlled-access protocol (WCG IRB #20252010).
> > The `.gitignore` denies by file extension and
> directory rather than by filename, deliberately.

---

## Verify the method in thirty seconds, without the data

The central methodological claim of this submission is that **depth alone cannot separate a real mosaic
chromosome gain from paralogous mismapping, and the shape of the heterozygous allele-balance
distribution can.** A reviewer has no quick way to check that against controlled-access data, so the
claim is made testable here with no data and no API key:

```bash
pip install pysam numpy
make validate            # or: python3 tests/test_discriminator.py
```

This builds a synthetic VCF with four contigs of known composition and runs the **unmodified analysis
code** from `scripts/aneu.py` against it. `chrGAIN` (a simulated 30% mosaic trisomy) and `chrMAP` (pure
paralogous mismapping) are given the **identical depth lift on purpose**:

| contig | ground truth | depth ratio | `AB<0.40` | `AB>0.60` | both tails | **hi/lo** |
|---|---|---|---|---|---|---|
| `chrNORM` | clean diploid | 1.000 | 0.084 | 0.092 | 0.176 | 1.10 |
| `chrGAIN` | mosaic trisomy, f=0.30 | **1.159** | 0.161 | 0.158 | **0.320** | **0.98** |
| `chrMAP` | paralogous mismapping | **1.159** | 0.306 | 0.005 | **0.311** | **0.02** |
| `chrHEMI` | hemizygous | 0.500 | 0.160 | 0.154 | 0.314 | 0.96 |

**Identical depth. Near-identical total spread. Forty-nine-fold apart on the discriminating statistic.**
A gain duplicates either allele with equal probability, so it widens the distribution symmetrically;
mismapping only ever adds reference-matching reads, so it widens the low side alone. `chrHEMI` returning
exactly 0.500 shows the method does resolve a real one-copy change, the role chrX plays as a built-in
positive control in the real analysis.

Seven assertions, two of which would fail if the discriminator were secretly just reading depth.
Committed output: [`results/discriminator_validation.txt`](results/discriminator_validation.txt).

---

## What this is

**Track 1, variant prediction.** *BUB1B* compound heterozygote, `NM_001211.6`:
`c.2210T>G p.Leu737Ter` (nonsense, NMD-competent) in trans with `c.3006T>G p.Asn1002Lys` (missense,
pseudokinase domain). Mosaic variegated aneuploidy type 1.
→ [`Mat34ec_track1_report.md`](Mat34ec_track1_report.md)

**Track 2, drug repurposing.** Mechanism characterization first, then five derisked candidate axes,
then the drugs this genome argues *against*. Lead candidate is nicotinamide, reaching BubR1 through
NAD⁺ and SIRT2, with Phase 3 human chemoprevention evidence at the same dose in a cancer-prone
population.
→ [`Mat34ec_track2_report.md`](Mat34ec_track2_report.md)

---

## Three findings that transfer beyond this case

**A trap that silently breaks rare-variant filtering on this dataset.** `INFO/AF` in this VCF is GATK's
per-sample *allele fraction*, not a population frequency. Filtering `AF < 0.01` against it discards
nearly every true heterozygous call, because a clean germline het sits near 0.5. The field name invites
the error. Population frequency must come from gnomAD.

**The symmetry discriminator above**, which needs no reference panel and carries its own positive
control in any male sample.

**AVI is frequency-naive by design.** In this genome a variant 8 kb from a candidate gene scored 1.149
as a genuine splice-acceptor variant and is carried by 52.7% of the population. Pairing AlphaGenome's
AVI with population frequency is not optional; used alone it generates confident false leads next to
the right gene.

---

## Layout

```
├── validate.py                    make validate | make panel | make clean
├── Mat34ec_track1_report.md       Track 1 report + methods description
├── Mat34ec_track2_report.md       Track 2 report + methods description
├── Mat34ec_track1_panel-vep.csv   Track 1 submission file
├── requirements.txt
├── tests/
│   ├── make_synthetic_vcf.py      builds the ground-truth fixture
│   └── test_discriminator.py      7 assertions, no data or API key needed
├── results/
│   └── discriminator_validation.txt
└── scripts/
    ├── make_panel_bed.py          builds the 12-gene panel live from Ensembl
    ├── annotate_and_triage.py     VEP annotation + recessive triage, both models
    ├── aneu.py                    per-chromosome and 10 Mb-bin aneuploidy scan
    ├── loh_scan.py                genome-wide subclonal LOH scan
    ├── panel_cnv_pgx.py           cryptic-deletion scan + pharmacogenomics by rsID
    ├── atlas_percentile.py        AlphaGenome Atlas saturation scan + percentile
    ├── atlas_decompose.py         AVI feature attribution
    ├── atlas_context_scan.py      1 Mb context scan around each locus
    ├── score_bub1b_alphagenome.py live-model scoring of the two candidates
    ├── score_controls_alphagenome.py synthetic local controls + empirical p
    └── check_*.py                 targeted VEP verification scripts
```

---

## Running the real analyses

These need the dataset and, for the AlphaGenome layers, an API key in the
environment. `aneu.py` and `loh_scan.py` take `--vcf PATH` or read `$MVA_VCF`, so nothing hardcodes a
directory layout.

```bash
export MVA_VCF=/path/to/proband.vcf.gz
export ALPHAGENOME_API_KEY=...          # never commit this

make panel                               # rebuild the panel from Ensembl
python3 scripts/annotate_and_triage.py --vcf panel_raw.vcf.gz --out candidates.tsv
python3 scripts/aneu.py 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y
python3 scripts/loh_scan.py 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
python3 scripts/panel_cnv_pgx.py
python3 scripts/atlas_percentile.py --gene BUB1B --mark "chr15:40220612:T>G=N1002K"
python3 scripts/atlas_decompose.py --gene BUB1B
```

Full sequence with timings in [`Mat34ec_track1_report.md` section 10](Mat34ec_track1_report.md). Everything
runs on a consumer laptop; the expensive steps are API-bound, not compute-bound.

SpliceAI needs its own environment, because it uses the Keras 2 API and the import fails outright under
Keras 3:

```bash
pip install spliceai "tensorflow-cpu>=2.13,<2.16" "keras>=2.13,<3"
```

---

## Tools and data sources

Ensembl REST and VEP · gnomAD · ClinVar · GENCODE · AlphaGenome and AlphaGenome Atlas · SpliceAI v1.3.1 ·
bcftools / htslib · pysam · ChEMBL · ClinicalTrials.gov · PubMed · Lens.org · EMA public assessment
reports.

**Generative AI disclosure.** Anthropic Claude Pro and Claude Science, Claude Opus 5, non-commercial
terms, no training on customer content. Google Gemini 3.1 Pro, non-commercial terms, no training on
customer content. Perplexity free tier as a search engine, non-commercial, no training on customer
content. Used for code authoring, literature retrieval, analytical design and adversarial review of the
author's own hypotheses. **Every quantitative claim in both reports derives from a script in this
repository or from a named primary source; no figure is reported from model recall.** Three cases where
checking a model-supplied claim against the primary source overturned it are recorded in
[`Mat34ec_track2_report.md`](Mat34ec_track2_report.md) section 10, and five cases where the pipeline's own
first reading was wrong are in [`Mat34ec_track1_report.md`](Mat34ec_track1_report.md) section 9.

---

## License

Code and reports released under **CC BY 4.0**, per the hackathon rules, so that methods developed for
this case can be reused for other undiagnosed individuals.

---

## Acknowledgement

To the child and family who chose to share a genome with strangers in the hope that someone could help.
This work exists only because of that decision.

This work was made possible through the Hackathon, organized by Sage Bionetworks in partnership with the MVA Society, Hugging Face, and BEACON (The Benchmarking, Evaluation, and Assessment Consortium for Science), with prize sponsorship from AWS and Anthropic.
