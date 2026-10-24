# Track 1: Variant Prediction

**Rare Disease, Real Kid: The MVA Hackathon 2026**
Participant `Mat34ec` · GRCh38 · Submission file `Mat34ec_track1_panel-vep.csv`

## The answer, up front

**A *BUB1B* compound heterozygote.** `p.Leu737Ter`, a nonsense allele that clears nonsense mediated decay, in trans with `p.Asn1002Lys` in the pseudokinase domain. Mosaic variegated aneuploidy type 1. EPCR 0.95.

It was found by examining 0.04% of the genome rather than filtering all of it.

That choice is the submission. The VCF holds 5,012,204 variants and its `INFO/AF` field is a per sample allele fraction, not a population frequency. The standard rare variant filter against that field deletes both causal alleles and reports no error. Subtractive pipelines fail silently. A twelve gene panel cannot, and being small is what made it affordable to score every possible variant across the candidate genes instead of a filtered subset.

Three things follow from that, and they are what this report is for:

- **The panel is built, not pasted.** Coordinates fetched live from Ensembl with provenance per gene, tiered by strength of evidence, with *MVK* deliberately excluded because it shares the MVA abbreviation and nothing else.
- **Five validation layers, none sharing assumptions with the triage.** Saturation scoring of 192,459 possible variants in *BUB1B*, splicing excluded by two model families with a positive control inside the same run, and consanguinity ruled out without parental samples.
- **Two negatives carrying detection floors, and five documented errors.** Including a mosaic aneuploidy signal on three chromosomes that looked real until a symmetry test killed it, and a pharmacogenomic result that reversed when the full panel was run.

The method can be checked in thirty seconds without the patient data: `python validate.py` in the repository.

## Contents

| | |
|---|---|
| [1. The prediction](#1-the-prediction) | What is being submitted, and the two secondary rows |
| [2. Why twelve genes](#2-why-twelve-genes) | The core argument: why a narrow panel beats a wide filter |
| [3. Three decisions that decided the answer](#3-three-decisions-that-decided-the-answer) | The traps, stated rather than buried in code |
| [4. Trying to break it](#4-trying-to-break-it) | Seven orthogonal attacks on the candidate |
| [5. Tested and negative, with floors](#5-tested-and-negative-with-detection-floors) | Mosaic aneuploidy and subclonal LOH, and how far down the test sees |
| [6. Secondary and incidental findings](#6-secondary-and-incidental-findings) | *LZTR1*, a documented exclusion, pharmacogenomics |
| [7. Limitations](#7-limitations) | What this cannot establish, starting with phase |
| [8. Methods description](#8-methods-description) | Template responses and the method abstract |
| [9. Where I was wrong](#9-where-i-was-wrong) | Five errors the validation layers caught |
| [10. Reproducing this](#10-reproducing-this) | Commands, and how to verify without the data |

## 1. The prediction

**Primary: a *BUB1B* compound heterozygote.** Transcript `NM_001211.6`, EPCR 0.95.

| | Variant | Consequence |
|---|---|---|
| Allele 1 | `chr15:40209701 T>G` · `c.2210T>G` · `p.Leu737Ter` | Nonsense, NMD competent |
| Allele 2 | `chr15:40220612 T>G` · `c.3006T>G` · `p.Asn1002Lys` | Missense, pseudokinase domain |

Mosaic variegated aneuploidy type 1 (MIM 257300): autosomal recessive, biallelic *BUB1B*.

Two secondary rows are also submitted. *LZTR1* `chr22:20996720 C>G` `p.Tyr748Ter`, monoallelic, EPCR 0.25. And *SLC34A1* `chr5:177386432 C>T` `p.Ala133Val`, heterozygous, EPCR 0.04, submitted as a documented exclusion rather than a candidate for the reasons in section 6.2. A pharmacogenomic profile is reported in section 6.3, in the text rather than as CSV rows, because it is not a causal variant claim.

## 2. Why twelve genes

The VCF holds 5,012,204 variants. The obvious move is to filter: drop the common ones, drop the benign consequences, drop whatever is left that looks uninteresting, and see what survives.

Filtering has a property worth naming. Every filter is an opportunity to delete the answer, and when the answer is deleted there is no error message. There is a clean empty result that looks like diligence.

In this dataset that is not hypothetical. The VCF's `INFO/AF` field is GATK's per sample allele fraction, not a population frequency. A standard rare variant filter of `AF < 0.01` against that field removes essentially every true heterozygous call in the genome, because a clean germline het sits near 0.5. Run it and both causal alleles are gone, and the pipeline still finishes without complaint. Section 3 covers the mechanics. The point here is the design consequence: subtractive pipelines fail silently.

So the problem is inverted. Rather than removing five million things and hoping the answer survived, select the small region where the answer has to be if the diagnosis is MVA, then examine everything inside it without filtering anything out.

The panel covers twelve genes and 1.3 Mb, 0.04% of the genome. That smallness is not a shortcut. It is what buys the depth. Inside 1.3 Mb it was affordable to be exhaustive in ways unavailable genome wide: annotate every variant through VEP rather than a filtered subset, score every possible single nucleotide variant across the candidate genes with AlphaGenome (192,459 of them in *BUB1B* alone), check every indel individually, test depth and heterozygosity gene by gene for cryptic deletions, and build synthetic local controls at the exact positions of interest. None of that is affordable across three gigabases. All of it was affordable here, and most of section 4 exists only because the search space was small enough to interrogate completely.

A narrow panel analysed exhaustively tells you more than a wide net filtered aggressively, because the narrow panel has no silent failure mode. If the answer is not in the panel, that is a visible, arguable, checkable claim about gene selection. If the answer was filtered out of a genome wide scan, nobody ever finds out.

**Which twelve, and why those.** Tiered by how firmly each gene is established in MVA rather than treated as a flat list.

Tier 1 is the numbered MVA types and the genes with established biallelic disease: *BUB1B*, *CEP57*, *TRIP13*, *CENATAC*, *SLF2*, *SMC5*, *MAD1L1*. This is where the diagnosis lives.

Tier 2 is reported candidates with weaker or single family support: *MAD2L1BP*, *BUB1*, *CEP192*.

Tier 3 is exploratory, included because the mechanism makes them plausible rather than because the literature is settled: *TRIM37*, *BUB3*.

The tiering decides what a hit means. A biallelic finding in Tier 1 is a diagnosis. The same finding in Tier 3 is a hypothesis needing a second family. A flat list hides that distinction from whoever reads the output.

**One gene deliberately left out.** *MVK* is not in the panel. Mevalonic aciduria abbreviates to MVA and is an unrelated metabolic disorder. A pipeline that includes *MVK* here selected a gene on string similarity rather than biology, and has made a category error before running a line of analysis. The error would be invisible in a result table while quietly contaminating the frequency statistics of everything else in the panel.

**The panel is built, not pasted.** Coordinates are fetched from Ensembl REST at run time, with alias fallback where the current symbol differs from older literature (*CENATAC* to *CCDC84*, *SLF2* to *FAM178A*, *MAD2L1BP* to *CMT2*). The script writes GRCh37 and GRCh38, both `chr` prefixed and bare contig BED files, plus a provenance TSV recording Ensembl gene ID, assembly, padding and retrieval date for every gene. Padding is 5 kb to catch promoter proximal and UTR variants.

A hardcoded coordinate list is faster to write and is the normal thing to do. It is also another silent failure: a pasted coordinate can be the wrong assembly, can be stale, can be off by a transcript, and nothing downstream reveals it. Fetching live makes the gene list auditable and regenerates the whole panel in about thirty seconds on any machine, which is what makes the rest of this report reproducible by someone else.

The pipeline is five steps:

```
1  panel construction      make_panel_bed.py       Ensembl REST, fetched live
2  region extraction       bcftools view -R
3  annotation              annotate_and_triage.py  Ensembl VEP REST
4  recessive triage        annotate_and_triage.py  both inheritance models
5  independent validation  section 4, methods not sharing assumptions with 1 to 4
```

## 3. Three decisions that decided the answer

Each could have silently produced the wrong result, so each is stated rather than buried in code.

**`INFO/AF` in this VCF is not a population frequency.** This is HaplotypeCaller output, where `AF` is the per sample allele fraction. Population frequency came only from gnomAD through VEP, never from the VCF. This is the single most consequential trap in the dataset and it is invisible if the field name is trusted.

**A variant absent from gnomAD is frequency zero, not unknown.** Dropping unseen alleles throws away exactly the class of ultra rare variant that causes severe recessive disease. `p.Leu737Ter` is absent from gnomAD. A filter requiring a frequency record loses the primary finding outright.

**Both inheritance models are scored, compound heterozygous and homozygous.** A triage that only ranks homozygous candidates cannot return a compound heterozygote, and a compound heterozygote is the answer. The pipeline pairs rare damaging alleles in the same gene explicitly, as a design requirement rather than an incidental capability.

## 4. Trying to break it

The triage produced a candidate in minutes. Everything below was run to attack it, using methods chosen not to share assumptions with the triage. Most of it was only possible because the panel was small.

**The impact is attributed to the protein, not to regulation.** AlphaGenome Atlas AVI feature decomposition for `p.Asn1002Lys`: 63% AlphaMissense, 19% PhastCons 470, 18% Cactus 241, 1.2% merged splicing, and about 2.8% across all nine regulatory tracks combined. The model attributes its own score to protein coding. This matters because AVI blends coding and regulatory signal, so a high AVI number alone would not have separated those explanations.

**Splicing is excluded by two model families, with a positive control inside the same run.** AlphaGenome ranked `N1002K` 21st of 22 synthetic local controls, those being other substitutions at the same position plus transitions at 5, 10, 20 and 30 bp either side, giving an empirical *p* of about 0.95. SpliceAI v1.3.1 returned 0.02, 0.00, 0.00, 0.00 against its 0.2 high recall threshold, in the same run where one local control scored 1.00 donor gain with 0.81 donor loss. A negative is worth only as much as the demonstration that the instrument registers a positive.

**The locus is saturated, not spot checked.** Atlas scored every possible SNV across *BUB1B* (192,459) and *LZTR1* (63,174), plus the 2,626 variants actually carried across both 1 Mb context windows. The two *BUB1B* variants rank 1 and 2 of those 2,626. The thirteen indels at the loci are all intronic or 3'UTR. This is the clearest return on the narrow panel: exhaustive scoring of a locus is available for twelve genes and not for twenty thousand.

**The parents are not consanguineous, which independently supports compound heterozygosity.** Excluding chrX, which scores homozygous throughout in a male and is an artefact in this context, the longest autosomal run of homozygosity is under 2.6 Mb and no segment exceeds 5 Mb. Recent consanguinity produces long tracts and there are none. Unrelated parents are exactly what two different alleles arriving from two sides requires. The argument needs no parental samples, which matters because there are none, and it also removes consanguinity as a rival explanation for the family's recurrent pregnancy loss.

**No cryptic deletion is masquerading as a second hit.** HaplotypeCaller does not call structural variants, so a deletion in trans with a point mutation would be invisible to everything above. Twelve of thirteen genes sit at 0.93 to 1.21 times baseline depth with normal heterozygous fractions. The low variant density in *BUB1B*, *BUB1*, *CEP57* and *CENATAC* is purifying selection on essential mitotic genes rather than deletion, and the giveaway is that their depth is normal. *TRIP13* is the exception, reported as unassessed in section 7.

**Both calls are clean germline heterozygotes.** `L737Ter` reads 21 reference and 25 alternate at depth 46. `N1002K` reads 15 and 13 at depth 28. Neither is a low fraction artefact and neither looks mosaic.

**The stop codon was confirmed against the reference, not inferred from the HGVS string.** The reference codon is `TTA` (Leu) on the plus strand, and `c.2210T>G` hits base 2, giving `TGA`, which is UGA, with an A in the +4 position. Reconstructed from the Atlas per position reference bases, independently of the annotation. It matters downstream because stop codon identity determines whether readthrough is worth considering at all, and UGA is the permissive one.

## 5. Tested and negative, with detection floors

A negative without a floor is not a result, so both of these carry one.

**Mosaic aneuploidy is below detection in this sample.** Every chromosome was scored across 3.6 million filtered biallelic sites, using depth and allele balance together. Chromosomes 20, 21 and 22 initially looked positive: raised depth, raised heterozygous fraction, wider allele balance. That is what a mosaic gain should look like.

The test that settles it is the shape of the widening, not its size. A real mosaic gain must widen the heterozygous distribution symmetrically, because the extra copy duplicates either allele with equal probability, putting modes at (1+f)/(2+f) and 1/(2+f). Mismapping cannot do that. Mismapped reads only ever add reference matching reads, so mismapping pushes allele balance down and never up.

The data is one sided. The fraction of sites above allele balance 0.60 is flat at 0.095 to 0.126 on every chromosome including the suspicious ones, while only the fraction below 0.40 moves. Binning at 10 Mb localises the excess: chr21's first bin runs at 2.34 times baseline, chr22's at 1.30, chr20's mid arm bin at 1.14, all acrocentric or pericentromeric, every other bin in the genome at 1.00 to 1.02. Depth and the low tail correlate at +0.63 across 286 bins, which is one mechanism rather than two. chrX returns a depth ratio of exactly 0.500, a free positive control proving the method resolves a real one copy change.

The floor is roughly 15% mosaic. MVA aneuploidy is scored cytogenetically in cultured lymphocytes or fibroblasts, this is uncultured blood, and aneuploid cells are plausibly selected against in rapidly dividing blood lineages. The negative is a statement about this sample and not about the child, and it does not substitute for a karyotype.

**No subclonal loss of heterozygosity anywhere in the genome.** Tested as E[z²] by sign across 292 autosomal bins, where z is the binomial deviate at heterozygous sites. The logic is the same asymmetry: copy number imbalance lifts both signs, mismapping lifts only the negative. The positive side sits at the null (median 1.19, IQR 1.16 to 1.24) while the negative side is inflated genome wide and one sided. Every top ranked bin is acrocentric or pericentromeric with a strongly negative mean z.

The floor here is 12 to 15% tumour fraction, far above the circulating tumour DNA fraction expected from a localized paediatric sarcoma. The question of a somatic second hit at *LZTR1* is therefore untestable from blood rather than excluded.

## 6. Secondary and incidental findings

### 6.1 *LZTR1* `p.Tyr748Ter`, monoallelic, EPCR 0.25

A genuine inactivating allele in a gene with a clean genotype to phenotype split, submitted because the presentation included a malignancy.

The codon is instructive. Tyr748 is `TAC`, and C>G makes it `TAG`, which is simultaneously a stop codon and an `AG` splice acceptor consensus. That is why SpliceAI returns DS_AG 0.92 at this position, and why an apparent disagreement between two tools resolved with both being right about different things.

There is no germline second hit anywhere in the locus or the surrounding 1 Mb. Under two hit logic a monoallelic inactivating *LZTR1* allele needs a somatic second hit, which cannot be assessed without tumour tissue. Section 5 gives the floor on why blood cannot answer it.

### 6.2 *SLC34A1* `p.Ala133Val`, EPCR 0.04, a documented exclusion

Submitted as a closed negative rather than a candidate, because the reasoning may be reusable. ClinVar is conflicting, carrying likely benign, likely pathogenic and uncertain significance submissions simultaneously, and gnomAD maximum frequency is 0.0101, too common for a severe recessive presentation.

The locus was swept on both sides. No second damaging coding allele exists. An AlphaGenome saturation scan of all 94,839 possible SNVs across the gene placed the only rare upstream allele at the 67th percentile, rank 30,966, below the gene's own 90th percentile and 35 times below its maximum. No phase tag exists anywhere in the gene, so a trans configuration would be unprovable even if there had been a second allele to phase.

A heterozygous allele with a conflicting classification in a recessive gene is a carrier finding. The nephrocalcinosis is better explained by prematurity and neonatal diuretic exposure.

### 6.3 Pharmacogenomic profile

Reported because the syndrome carries a high tumour risk, which makes future systemic therapy foreseeable, and because the VCF carries populated rsIDs so these loci can be read directly.

Reference and therefore normal function: all four *DPYD* risk alleles (`rs3918290`, `rs55886062`, `rs67376798`, `rs56038477`), all three *TPMT* alleles and *NUDT15* `rs116855232`, *G6PD* (hemizygous in a male, so a single variant would be fully penetrant), *SLCO1B1\*5*, *UGT1A1*, *RARG* and *SLC28A3* for anthracycline cardiotoxicity, and *ACYP2* for cisplatin ototoxicity. Clean thiopurine and fluoropyrimidine status in a child whose syndrome predisposes to leukaemia is usable information now rather than later.

Carried: *CEP72* `rs924607` homozygous, *MTHFR* C677T homozygous, *CYP2B6* `*9` and the `*6` tag both heterozygous, *CYP2C19* `*2` and `*17` both heterozygous, *GSTP1* I105V homozygous, *XPC* `rs2228001` homozygous, the *ABCB1* 1236, 2677 and 3435 haplotype heterozygous throughout, and *SOD2* A16V, *CBR3*, *SLCO1B3* `rs4149117` and *ABCC4* `rs11045585` heterozygous.

Mitochondrial genome: the aminoglycoside ototoxicity alleles `m.1555A>G` and `m.1494C>T` are both absent. The only variant in that region is the common benign `m.1438A>G`, homoplasmic at depth 5,212.

**One result that reversed on checking, reported because it reversed.** *CEP72* `rs924607` came back homozygous for the reported vincristine neuropathy risk genotype. CEP72 is a centrosomal protein, so the convergence with a mitotic lesion looked striking, and I wanted it to be real.

Checking the literature rather than asserting it changed the reading twice. The association is contested: a Spanish B-ALL cohort found none during induction, and an 817 patient meta analysis could not confirm it globally, finding it only in the continuation phase, while the largest single analysis at 1,095 patients does support it at an odds ratio of 2.15.

More decisively, the same literature names loci with larger effects and the proband is reference at all of them: *ABCC1* `rs3784867` at a reported odds ratio of 4.9, *SLC5A7* `rs1013940` at 8.6, *TTPA* `rs10504361` at 2.0, and the *CEP72* missense `rs12522955` at 2.3. Of five reported vincristine neuropathy loci he carries one, and it is the weakest and most disputed of the five.

A caution flag rather than a contraindication, and considerably softer than the single locus result looked on its own. Caveats in section 7.

## 7. Limitations

1. **Phase is inferred, not observed.** There are no parental samples, and the two variants are 10,911 bp apart while the reads are 2×149 bp, so no read or read pair can span them. Read backed phasing is impossible in principle here rather than merely unavailable. Phase is inferred from the phenotype matching a recessive disorder. This is the largest unverified assumption in the submission.

2. ***TRIP13* is effectively unassessed.** Median depth 7 against a baseline of 44, with three variant sites across 36.5 kb. This is a subtelomeric coverage failure rather than a deletion, and the reason that can be stated is that heterozygous calls cannot exist inside a homozygously deleted span. But a Tier 1 MVA gene has no usable data in this dataset.

3. **No structural variant calling.** Only SNV and indel calls are available. The depth and heterozygosity scan in section 4 is a partial substitute, not a replacement for running an SV caller on aligned reads.

4. **Pharmacogenomics is a single variant screen, not star allele diplotyping.** *CYP2D6* in particular cannot be genotyped reliably from short reads because of *CYP2D7* paralogy and copy number variation. The `rs3892097` call shows an allele balance of 0.28, consistent with exactly that paralogy, and no *CYP2D6* diplotype is asserted. Absence of a record means homozygous reference only if the position was covered, unconfirmable without a BAM. Section 6.3 is a screen rather than a clinical pharmacogenomic report and should be confirmed against current CPIC tables before use.

5. **The negatives in section 5 carry the floors stated there** and are not claims that the phenomena are absent in the child.

## 8. Methods description

*Responses to the provided Track 1 methods template.*

**Team name.** `Mat34ec`, individual participant.

**Model number.** 1 of 1.

**Describe your model or approach in detail.** Sections 2 through 5 above.

**Generative AI disclosure (required).** Anthropic Claude Pro and Claude Science, Claude Opus 5, non commercial terms, no training on customer content. Google Gemini 3.1 Pro, non commercial terms, no training on customer content. Perplexity free tier as a search engine, non commercial, no training on customer content. AI assistance was used for code authoring, literature retrieval, analytical design and adversarial review of my own hypotheses. Every quantitative claim in this report comes from a script in the repository and can be re-run; no number here is reported from model recall. Claims drawn from literature were checked against the primary source, and section 9 records where that check overturned an initial reading.

**Is the submission file automated output or manually curated?** The variant identities are automated output. The pipeline returns the *BUB1B* pair without human intervention. Manual review was applied to the EPCR values and to the decision to include the two secondary rows.

**Describe any downstream manual review.** EPCR values were assigned by hand rather than computed. The ordering is the pipeline's, the calibration is a judgement. The *SLC34A1* row was added manually, specifically to document a closed negative rather than to propose a candidate.

**Public or proprietary data?** Public only, exclusively. No proprietary sources.

**Describe the public data used.** Ensembl REST for gene coordinates, sequence and VEP consequences; gnomAD allele frequencies through VEP; ClinVar through VEP; GENCODE for SpliceAI annotation; AlphaGenome and AlphaGenome Atlas for AVI scores, feature decomposition and saturation scans; SpliceAI v1.3.1; dbSNP rsIDs as carried in the VCF; PubMed for the literature in sections 6.3 and 9.

**Describe any proprietary data used.** None.

**Can the approach output compound heterozygous pairs?** Yes. Both models are scored and same gene pairs are emitted explicitly. This was a design requirement rather than an incidental capability, because a pipeline without it cannot return the correct answer in this case.

**How were secondary or incidental findings handled?** Flagged with `finding_type=secondary` and a rationale in the `notes` column. Three classes are represented deliberately: a genuine candidate contributor (*LZTR1*), a documented exclusion (*SLC34A1*), and a pharmacogenomic profile reported in the text rather than as CSV rows because it is not a causal variant claim.

**Run time and cost.** Panel construction about 30 seconds, network bound. Region extraction and triage about 5 minutes on a laptop. The validation layers are the expensive part: Atlas saturation scans about 45 minutes across three genes, the per chromosome aneuploidy scan about 50 seconds, the genome wide LOH scan about 60 seconds, SpliceAI about 3 minutes on CPU, and the panel CNV and pharmacogenomic scan about 17 seconds. No specialised hardware. A consumer laptop plus API access.

**Method abstract.** The VCF holds five million variants, and the usual approach is to filter them down. This submission inverts that, because every filter can delete the answer without reporting an error, and this dataset demonstrates the risk directly: its `INFO/AF` field is a per sample allele fraction rather than a population frequency, so a routine rare variant filter against it discards nearly every true heterozygous call, including both causal alleles, and terminates cleanly.

Instead of filtering a genome, a twelve gene MVA panel covering 1.3 Mb, or 0.04% of the genome, was built and examined exhaustively. The panel is fetched live from Ensembl REST with full provenance rather than pasted from hardcoded coordinates, tiered by how firmly each gene is established in the disease, and deliberately excludes *MVK*, which shares the MVA abbreviation but is an unrelated disorder.

The smallness is what made the rest possible. Within 1.3 Mb everything could be annotated through VEP, triaged under both homozygous and compound heterozygous models, and then attacked with methods that do not share the triage's assumptions. The pipeline returns a *BUB1B* compound heterozygote: `p.Leu737Ter`, NMD competent, together with `p.Asn1002Lys` in the pseudokinase domain.

Five validation layers follow. AVI feature decomposition attributes the missense impact to protein coding rather than regulation. Splicing is excluded by AlphaGenome against 22 synthetic local controls and by SpliceAI with a positive control in the same run. Every possible SNV across both loci was scored, plus the 2,626 variants actually carried in the surrounding 1 Mb, and the two *BUB1B* variants rank first and second. A runs of homozygosity analysis shows no consanguinity, which supports compound heterozygosity without parental samples. A depth and heterozygosity scan rules out cryptic deletions across the panel.

Two negatives are reported with explicit floors. Mosaic aneuploidy is below roughly 15% in this uncultured blood sample, established by a symmetry test: a true copy number gain widens the heterozygous allele balance distribution symmetrically, while paralogous mismapping only ever pushes it low, and chrX provides a built in positive control at a depth ratio of exactly 0.500. No subclonal loss of heterozygosity is detectable above 12 to 15% tumour fraction genome wide.

Strengths: the gene list is auditable and regenerates in thirty seconds; validation is by orthogonal methods rather than repetition; negatives carry detection floors; and five instances where the pipeline's own first reading was wrong are documented. Limitations: phase is inferred rather than observed and cannot be observed with 2×149 bp reads across an 10,911 bp gap; *TRIP13* is unassessed because of subtelomeric coverage failure; no structural variant calling is possible without aligned reads; and the pharmacogenomic screen is single variant rather than star allele diplotyping.

## 9. Where I was wrong

A method that never reports its own failures gives a reviewer no way to calibrate it. All five of these were caught by the validation layers rather than by luck, which is the argument for building them.

An apparent 11 fold splice junction signal at `L737Ter` looked real until synthetic local controls were run, at which point a control 30 bp away scored higher and the empirical *p* came out at 0.24. The signal was an artefact of comparing two different loci instead of comparing a variant to its own neighbourhood.

AVI was initially judged to have a gene model error at *LZTR1*. SpliceAI later returned DS_AG 0.92 at the same position, vindicating AVI. The codon explains both tools at once, as described in section 6.1.

An apparent mosaic aneuploidy signal on chromosomes 20, 21 and 22 survived until the symmetry test in section 5 was applied.

A locus sweep script reported three "second candidate alleles" at *SLC34A1* which turned out to lie 5 to 6 kb beyond the gene's 3' end and belong to a neighbouring transcript. The script had fallen back to `most_severe_consequence` when no *SLC34A1* transcript matched. The conclusion survived, but it had been right for the wrong reason until corrected.

A single locus pharmacogenomic lookup suggested a vincristine neuropathy risk genotype. Running the full panel reversed the reading, as described in section 6.3.

## 10. Reproducing this

```bash
git clone https://github.com/Mat34ec/mva-hackathon-2026 && cd mva-hackathon-2026
pip install -r requirements.txt

# 1. build the panel from Ensembl, no hardcoded coordinates
python scripts/make_panel_bed.py --build GRCh38 --pad 5000 --out mva_panel

# 2. extract panel regions
bcftools view -R mva_panel.GRCh38.nochr.bed -Oz -o panel_raw.vcf.gz "$MVA_VCF"
tabix -f -p vcf panel_raw.vcf.gz

# 3. annotate and triage, both inheritance models
python scripts/annotate_and_triage.py --vcf panel_raw.vcf.gz --out candidates.tsv

# 4. validation layers (the API key lives in the environment, never in a file)
python scripts/atlas_percentile.py --gene BUB1B --mark "chr15:40220612:T>G=N1002K"
python scripts/atlas_decompose.py  --gene BUB1B
python scripts/aneu.py      1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y
python scripts/loh_scan.py  1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
python scripts/panel_cnv_pgx.py
```

The method itself can be checked without the data. `python validate.py` builds a synthetic VCF with four contigs of known composition and runs the unmodified aneuploidy code against it, confirming that the symmetry test in section 5 separates a real mosaic gain from paralogous mismapping even when both are given an identical depth lift. Seven assertions, no patient data and no API key needed. The repository holds the full instructions and the access terms.

---

*Released under CC BY 4.0 in accordance with the hackathon rules.*
