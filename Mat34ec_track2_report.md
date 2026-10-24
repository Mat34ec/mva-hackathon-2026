# Track 2 — Drug Repurposing

**Rare Disease, Real Kid: The MVA Hackathon 2026**
Participant: `Mat34ec` · Proband `WGS_EX2312012` · GRCh38

> These are **hypotheses for follow-up**, not evidence that any medicine works, and not medical advice.
> Every candidate is stated with what must be true for it to work, and with the experiment that would
> show it does not.

---

## Summary

The child is a *BUB1B* compound heterozygote: one null allele and one destabilizing missense allele in
the pseudokinase domain. The lesion is therefore **partial loss of BubR1, not absence** — a hypomorph.
Everything downstream follows from *how much* functional BubR1 a cell has, which is what makes the
disease addressable at all: there is a quantity to move, not an absence to replace.

That quantity has a named biochemistry. BubR1 recruits PP2A-B56 to the kinetochore, PP2A-B56 opposes
Aurora B, and losing that counterweight is what produces premature chromatid separation and mosaic
aneuploidy. **The defect is not a missing function. It is a lost brake on a function that is supposed to
be there** — which is why restoring a balance, rather than replacing a protein, is the tractable move.

Working from that mechanism rather than from the disease name produces three things:

1. **A lead candidate that is approved, cheap, and has Phase 3 human cancer-prevention evidence** —
   nicotinamide, reaching BubR1 through NAD⁺ and SIRT2.
2. **A list of drugs that must not be given**, including an entire class that the obvious analysis would
   have recommended.
3. **An honest map of where the field actually is**: in this disease, evidence quality and clinical
   availability run in opposite directions, and saying so is more useful to a clinician than picking a
   point on that curve and calling it an answer.

---

## 1. Mechanism characterization

*Required by the track: loss- or gain-of-function, pathway disrupted, downstream consequence.*

### 1.1 The two alleles

| | Variant | Class | Consequence |
|---|---|---|---|
| **Allele 1** | `chr15:40209701 T>G` · `NM_001211.6:c.2210T>G` · `p.Leu737Ter` | Loss of function, complete | Reference codon `TTA`; `T>G` at base 2 gives `TGA`, i.e. **UGA**, +4 context **A**. More than 50 nt upstream of the last exon–exon junction, therefore NMD-competent. Functionally a null. |
| **Allele 2** | `chr15:40220612 T>G` · `c.3006T>G` · `p.Asn1002Lys` | Loss of function, partial | Missense in the pseudokinase domain. Asparagine to lysine: neutral polar replaced by longer and positively charged. |

The stop codon was read off the reference sequence rather than inferred from the HGVS string, because
the codon identity determines whether one of the candidate strategies is even plausible (§3.5).

**The lesion is a hypomorph.** Complete BubR1 loss is embryonic lethal; this child is alive, which
constrains the model independently of any assay. Allele 1 contributes no protein. Allele 2 contributes
protein of reduced competence. **Allele 2 is therefore the only rescuable allele in the genome**, and
total functional BubR1 is the variable every candidate below is trying to move.

### 1.2 The chain, link by link

1. **BubR1 is a pseudokinase.** Its C-terminal domain adopts a kinase fold, but catalysis is
   dispensable and the nucleotide-binding residues serve conformational stability. **The protein is a
   scaffold, not an enzyme.** This single fact decides the shape of the intervention: on a scaffold, the
   mechanistically correct ligand is a *stabilizer*, not an inhibitor, because there is no catalysis to
   inhibit.
2. **The KARD motif is phosphorylated**, licensing recruitment of **PP2A-B56** to the kinetochore.
3. **PP2A-B56 exists to oppose Aurora B.** Aurora B phosphorylates kinetochore substrates to loosen
   microtubule attachments, which is how a cell discards incorrect ones. PP2A-B56 dephosphorylates the
   same substrates to lock in correct, tension-bearing attachments. The kinetochore needs both.
4. **Less BubR1 means less PP2A-B56, so Aurora B runs unopposed.** Ndc80 phosphorylation rises and
   attachments that should have been kept are torn down along with the bad ones.
5. **Chromosomes fail to congress, sister chromatids separate prematurely**, and daughter cells inherit
   wrong chromosome counts. Because the failure is probabilistic and recurs at every division, different
   cells acquire different aneuploidies — hence *mosaic variegated*, a patchwork rather than one uniform
   trisomy.

**Link 3–4 is supported directly and the support is specific.** The BUBR1:B56-PP2A complex opposes
Aurora B kinase activity, **loss of the complex can be reverted by inhibiting Aurora B**, and failure to
recruit B56-PP2A **contributes to the chromosome congression defects found in cells derived from
patients with MVA syndrome** (Xu et al. 2013, *Biology Open*, PMID 23789096,
[DOI](https://doi.org/10.1242/bio.20134051)).

### 1.3 An independent arrival at the same lesion

The chain above was built from the checkpoint literature. A completely separate literature — on how
BubR1 abundance is controlled — arrives at the same place, which is the strongest internal check this
model has.

BubR1 is acetylated at **Lys250** by PCAF, and that acetylation **blocks ubiquitin-dependent BubR1
degradation** by APC/C-Cdc20. The acetyl-mimetic K250Q is stabilized; the acetyl-deficient K250R is
unstable. Acetylation is "a molecular switch that converts BubR1 from an inhibitor to a substrate of the
APC/C complex" (Choi et al. 2009, *EMBO J*, [DOI](https://doi.org/10.1038/emboj.2009.123)).

Mice heterozygous for acetylation-deficient BubR1 **spontaneously develop tumours with massive
chromosome missegregation** and a weakened checkpoint, caused by premature degradation of the protein.
Decisively: "the association of acetylation-deficient BubR1 with **PP2A-B56α** phosphatase was reduced"
with elevated phosphorylated Ndc80 at the kinetochore, i.e. **excessive Aurora B activity**
(Park et al. 2013, *J Cell Biol*, [DOI](https://doi.org/10.1083/jcb.201210099)).

**Two independent routes to the same lesion.** Reduce BubR1 by mutation, or reduce it by blocking its
acetylation, and you get the same biochemistry and the same disease phenotype. That is what a correct
mechanism looks like. It also establishes the therapeutic principle that organizes §3: **BubR1 abundance
is the pharmacologically meaningful variable.**

### 1.4 Caretaker, not gatekeeper

MVA is a **caretaker** defect: it raises the rate of missegregation and mutation rather than directly
driving proliferation. A gatekeeper lesion removes a brake on growth; a caretaker lesion corrupts the
machinery that keeps the genome accurate.

The practical consequence governs this whole report. **Correcting the caretaker defect does not treat an
existing tumour, and treating a tumour does not correct the caretaker defect.** They are two problems
and this document keeps them apart. It also points at the therapeutic category that actually fits a
caretaker disease, which is **chemoprevention** rather than cure — and that is where the lead candidate
turns out to live.

Checkpoint severity tracks tumour risk across the MVA genes: biallelic *TRIP13* or *BUB1B* patients show
severe checkpoint impairment and high embryonal tumour risk, while *CEP57* and unexplained MVA show
minimal deficiency and no embryonal tumours (Yost et al. 2017, *Nature Genetics*,
[DOI](https://doi.org/10.1038/ng.3883)). This child is in the high-risk group, which is why the
oncological axis cannot be dropped.

### 1.5 The second locus, and what it is actually for

*LZTR1* `p.Tyr748Ter` is monoallelic. LZTR1 is a CUL3 substrate adaptor that ubiquitinates RAS for
degradation; monoallelic inactivating alleles are the schwannomatosis genotype and require a somatic
second hit under two-hit logic.

**It is not downstream of BUB1B.** There is no signalling link between the spindle checkpoint and the
RAS–MAPK cascade and this report claims none. The only connection is mutational: chromosomal instability
plausibly supplies the somatic second hit that a monoallelic *LZTR1* allele needs.

Its value is easy to misstate, so stated precisely: **not** that it explains the tumour better than
*BUB1B* does, but that if the somatic hit is real, RAS–MAPK activation exists **only in tumour cells**.
In a child whose every other cell is uniformly abnormal, that is the single place where tumour and
normal tissue genuinely differ. **It is the only surviving selectivity axis in this genome.**

---

## 2. How candidates were chosen

Not by searching for drugs associated with MVA — there are none — but by asking, for each link in §1.2,
whether an approved medicine moves it, and then asking what would have to be true for that medicine to
help **this** child rather than a generic patient.

Three filters, applied in order:

1. **Does it act on a named link in the chain?** A drug that treats a symptom is not a repurposing
   hypothesis.
2. **Does an approved medicine exist?** The track asks for existing approved medications. Where none
   exists, the axis is reported as a mechanistic finding rather than dressed up as a candidate.
3. **What breaks it in this child specifically?** His kidneys, his pharmacogenomics, his cancer
   predisposition, the one published treatment record in his disease.

Filter 3 is the one that did the most work, and it is also what produced §4, the list of drugs that must
not be given.

---

## 3. The candidates, derisked

### 3.1 Nicotinamide — raise BubR1 abundance · **lead**

**The mechanism, and why the direction is right.** BubR1 abundance falls as NAD⁺ declines, through
SIRT2 acting at **Lys668**, counteracted by the acetyltransferase CBP. **Overexpression of SIRT2, or
treatment with the NAD⁺ precursor NMN, increases BubR1 abundance in vivo**, and SIRT2 overexpression in
**BubR1^H/H hypomorphic mice** increases median lifespan (North et al. 2014, *EMBO J*,
[DOI](https://doi.org/10.15252/embj.201386907)).

That mouse is the right mouse. It is hypomorphic, not knockout — reduced BubR1 rather than absent —
which is this child's genotype expressed in an animal.

**Four independent lines converge, and none was constructed to meet the others.**

**(i) Mechanism is BubR1-specific and in the correct direction.** Not a generic proteostasis argument.
The target is the protein itself and the effect is an increase.

**(ii) There is Phase 3 human evidence, at this dose, in a cancer-prone population.** ONTRAC: 386
participants with at least two prior non-melanoma skin cancers, randomised double-blind to
**nicotinamide 500 mg twice daily** or placebo for 12 months. New cancers fell by **23%** (95% CI 4–38,
p=0.02), actinic keratoses by 13–20%, and there were **"no noteworthy between-group differences... with
respect to the number or types of adverse events"** (Chen et al. 2015, *NEJM*,
[DOI](https://doi.org/10.1056/NEJMoa1506197)).

This matters because of §1.4. A caretaker disease's therapeutic category is chemoprevention, and
chemoprevention with this specific agent has a Phase 3 precedent in people already prone to cancer.

**(iii) It may act on the one phenotype *BUB1B* does not explain.** Nephrocalcinosis is
calcium-phosphate deposition. Nicotinamide lowers serum phosphate, the calcium-phosphorus product and
iPTH, and improves calcification propensity, in randomised trials
([DOI](https://doi.org/10.1016/j.ekir.2020.12.012);
[DOI](https://doi.org/10.3265/Nefrologia.pre2010.Nov.10734)). The COMBINE trial found no significant
effect over 12 months in CKD and describes it as "an inhibitor of active intestinal phosphate transport"
([DOI](https://doi.org/10.1681/ASN.2018101058)), so the effect is real but not reliable.

This is **not** proposed as a treatment for the nephrocalcinosis. It is flagged because a drug selected
on mitotic grounds landing on the unexplained renal phenotype is a coincidence that should be stated and
tested, not quietly used as support.

**(iv) Availability is not a constraint.** Approved 1982, oral and parenteral, ATC A11HA01, no black box
warning, decades of paediatric use, and cheap enough to reach every family the MVA Society serves.

**Now the three things that must be true and are not yet known.**

**(a) The sirtuin tension — the sharpest objection to the whole idea.** The published mechanism requires
SIRT2 **activity**. Nicotinamide is an NAD⁺ precursor *and* a sirtuin inhibitor by product inhibition,
so above some concentration it opposes the enzyme the mechanism depends on. This is a dose-window
question, not a yes/no one, and it is unresolved. **Nicotinic acid reaches NAD⁺ by the Preiss-Handler
pathway without that liability** and must be evaluated head to head. Any proposal naming nicotinamide
without naming this has not read the mechanism.

**(b) Thrombocytopenia — patient-specific, not generic.** A randomised trial in haemodialysis patients
found falling platelet counts at **500–1000 mg/day, lower doses than other studies**, with the authors
explicitly asking for safety work in other populations
([DOI](https://doi.org/10.3265/Nefrologia.pre2010.Nov.10734)). ONTRAC saw no adverse-event signal at the
same dose in a healthier group, so population matters — and it matters here in a specific way. **The
single published account of chemotherapy in an MVA patient was dose-limited by thrombocytopenia after
cycle 2.** A child whose disease's only treatment record ran into platelets does not receive a chronic
agent with a platelet signal without monitoring from day one.

**(c) ONTRAC's protection is probably not the mitotic mechanism.** Nicotinamide's effect against
UV-driven skin cancer is attributed to NAD⁺ supporting DNA repair and limiting photoimmunosuppression.
That is a different route to cancer than missegregation. **ONTRAC establishes safety and the
plausibility of chemoprevention with this agent; it does not establish prevention of aneuploidy-driven
tumours.** Treating it as if it did would be the same category error this report rejects elsewhere.

**The gate.** Before: platelet count, serum phosphate, calcium, PTH, 1,25-vitamin D — because the
phosphate effect could help or harm depending on which mechanism drives his nephrocalcinosis, and that
is unknown. During: serial platelets. The decisive preclinical experiment: does the agent raise BubR1
protein **and** restore PP2A-B56 kinetochore recruitment in patient-derived cells, with nicotinamide and
nicotinic acid compared across a dose range to find the sirtuin window.

---

### 3.2 MEK inhibition — the tumour-restricted axis · **gated**

**Rationale.** If the somatic second hit at *LZTR1* is real, RAS protein accumulates and MAPK signalling
is elevated **in tumour cells only** (§1.5). MEK inhibition then acts on the one lesion in this child
that is not constitutional. Separately, rhabdomyosarcoma is a differentiation-blocked muscle tumour and
ERK signalling restrains MYOG-driven differentiation, so there is a differentiation rationale as well as
an anti-proliferative one.

**A second leg for the differentiation argument, native to BUB1B.** BubR1 insufficiency *by itself*
suppresses myogenic differentiation: BubR1-hypomorphic C2C12 myoblasts show markedly reduced myotube
formation with lower **myogenin and MyoD**, and raised Atrogin-1 and MuRF-1 driving proteasomal muscle
protein loss (Nozaki et al. 2015, *Cytotechnology*,
[DOI](https://doi.org/10.1007/s10616-015-9920-7)). The differentiation argument does not depend on RAS
alone.

**Availability is the strongest part of this candidate.** Selumetinib is approved for paediatric NF1
plexiform neurofibroma with a mature paediatric safety dataset, including a post-authorisation safety
study and a 409-patient observational cohort. MEK inhibition in germline RASopathy is itself in trial
(trametinib for RASopathy-associated infantile hypertrophic cardiomyopathy, NCT07817186).

**Why it is gated rather than proposed.** The somatic hit is a hypothesis, not an observation. A
genome-wide subclonal LOH scan of blood was negative with a **floor of 12–15% tumour fraction**, far
above the circulating tumour DNA fraction expected from a localized paediatric sarcoma. At the
`p.Tyr748Ter` position the read split is **24 reference to 24 mutant at depth 48** — dead even. This is
**untestable from blood, not excluded.**

> **Gate.** If tumour tissue shows *LZTR1* loss of heterozygosity or a second somatic mutation, this arm
> opens. If not, it stays closed. Proposing MEK inhibition without that test is proposing to treat a
> state that has not been shown to exist in this child.

Stating the gate is the proposal's most rigorous feature, not its weakness: it converts the submission's
weakest assumption into a decision with a defined answer.

**Derisking.** MEK inhibition in children is not quality-of-life-neutral — rash, paronychia, GI
toxicity, CK elevation, reduced ejection fraction, ocular toxicity, with dose interruptions common. In a
child with nephrocalcinosis and a prematurity history, chronic exposure needs its own justification and
a specified duration.

---

### 3.3 Chemical chaperones — stabilize the allele 2 fold

**Rationale from first principles.** Allele 1 is a null and cannot be stabilized. Allele 2 produces
protein whose pseudokinase domain is destabilized, and since that domain's role is **conformational
rather than catalytic** (§1.2, link 1), a ligand occupying the nucleotide pocket would act as a
*pharmacological chaperone* rather than an inhibitor. There is no catalysis to inhibit.

**The fold is demonstrably ligand-stabilizable.** A 2025 cellular thermal shift dataset against BUB1B in
CAKI-1 cells reports four compounds: two reduce BubR1 surviving a 45 °C challenge (fold change 0.29 and
0.65) and **two increase it (1.33 and 2.19, annotated "Increase In Heat Resistance Observed")** (ChEMBL
assay `CHEMBL6125447`, document `CHEMBL6121040`). Both stabilizers are preclinical compounds. The only
approved molecule with measured BubR1 engagement is **ceritinib**, at IC50 13 µM — weak, and an
inhibition assay, which is the wrong direction. Dasatinib was tested and does not engage.

| Agent | Status | Objection |
|---|---|---|
| **Arimoclomol** (Miplyffa) | Approved 2024, oral, paediatric NPC, ATC N07XX17; USAN stem `-clomol` = "heat-shock protein inducers" | No evidence BubR1 is an HSP90/HSP70 client. Searched specifically; the literature at that intersection is 22 papers and none establish it. Rationale is class-level |
| **Ursodeoxycholic acid** | Approved, decades of paediatric cholestasis use | Non-specific; no BubR1 data |
| **Betaine** (Cystadane) | Approved, infancy onward | Non-specific; no BubR1 data |

**The trap, disarmed before a reviewer finds it.** Thermal stabilization is **not** functional rescue.
Two of the four compounds in that CETSA dataset are annotated as BUB1B inhibitors. A ligand can rigidify
the fold and simultaneously lock out the conformation KARD presentation and PP2A-B56 recruitment depend
on. One would measure a textbook thermal shift and have made the child worse.

**So the success criterion is defined in advance: CETSA *plus* a functional readout.** A candidate passes
only with a positive thermal shift on BubR1 **and** preserved or improved KARD phosphorylation or
PP2A-B56 recruitment in the same system. A thermal shift alone is a disqualifying result, not a
promising one.

---

### 3.4 PP2A activation — the approved way onto the Aurora B axis

**The axis with the best evidence has no approved drug.** ChEMBL lists **nine** Aurora B inhibitors with
annotated mechanism; the maximum clinical phase is **3** and none is approved. Aurora B inhibition is
therefore reported here as a **mechanistic keystone, not a candidate** — it is what establishes the
PP2A-B56/Aurora B axis as a *correctable* lesion rather than merely a described one (§1.2).

**But the lesion is lost PP2A-B56 recruitment, and the damaged side has approved agents.** The
sphingolipid analogue **fingolimod (FTY720)** binds and targets the endogenous PP2A inhibitor
**I2PP2A/SET**, producing **PP2A reactivation** ([DOI](https://doi.org/10.1002/emmm.201201283)),
reproduced independently ([DOI](https://doi.org/10.1292/jvms.15-0193)) and reviewed as a strategy
([DOI](https://doi.org/10.3389/fonc.2014.00388)). Fingolimod is approved, oral, licensed from age 10.

**Three objections, all more informative than the proposal.**

1. **Targeting failure, not activity failure.** BubR1 cannot bring B56-PP2A *to the kinetochore*.
   Relieving SET inhibition raises global PP2A activity; it does not obviously deliver B56-PP2A to the
   place the lesion emptied. A localization problem, not a dose problem, and unresolved.
2. **Immunosuppression in a cancer-predisposition syndrome.** Fingolimod's approved mechanism is S1P
   receptor modulation and lymphocyte sequestration. This may be disqualifying by itself.
3. **The PP2A evidence comes from tumour suppression**, not kinetochore biology. Right phosphatase,
   studied for a different reason.

**What this axis yields is a specification rather than a prescription**, and that is a legitimate output
of a repurposing analysis: the agent this child needs restores B56-PP2A **at the kinetochore
specifically**, without global phosphatase activation and without touching lymphocytes. It does not
exist. Writing the specification down is how one gets made.

**And one caution that applies to the whole axis.** Aneuploid cells are more sensitive to **genetic**
perturbation of checkpoint components and *less* sensitive to short-term **pharmacological** inhibition
of the same components (Cohen-Sharir et al. 2021, *Nature*,
[DOI](https://doi.org/10.1038/s41586-020-03114-6)). The rescue evidence sits on the genetic side. Any
pharmacological attempt inherits that gap.

---

### 3.5 Readthrough of allele 1 — an enabling experiment without an agent

**Rationale.** Allele 1 contributes nothing. Recovering a fraction converts a null into a partial
contributor, and in a dose-dependent hypomorphic disease that is the difference between two
half-measures and one.

**What the sequence says, measured rather than assumed.** The premature stop is **UGA**, +4 context
**A**. UGA is the most readthrough-permissive of the three stop codons, reported as UGA > UAG > UAA
([DOI](https://doi.org/10.1016/j.bbrc.2022.10.081)), which favours the approach. The +4 nucleotide is an
established strong modifier but its hierarchy is **system-dependent rather than universal**: in *HERG*
the best-responding context was TGA-T ([DOI](https://doi.org/10.3892/ijmm.2019.4360)), and stop-codon
context determines not only whether readthrough works but which drug and what concentration
([DOI](https://doi.org/10.3390/biom11071006); review,
[DOI](https://doi.org/10.1080/15476286.2015.1068497)).

**There is currently no usable approved agent, and that is the finding.**

| Agent | Status |
|---|---|
| **Ataluren** (Translarna) | **EU authorisation not renewed.** The EMA recommended non-renewal, confirmed it, and re-confirmed it on re-examination; never FDA-approved. ChEMBL still records max phase 4 with `withdrawn_flag: false`, which is stale — a caution about trusting a database field over the regulator |
| Aminoglycosides (gentamicin, amikacin) | Approved, but **vetoed in this child**: nephrotoxic against existing nephrocalcinosis. The genome refines rather than compounds this — MT-RNR1 ototoxicity alleles `m.1555A>G` and `m.1494C>T` are **absent**, so the veto rests on renal risk alone |

So this is honestly **an enabling experiment in search of an agent**, not a treatment proposal. It stays
in the report because the codon is the favourable one and the resolving experiment is cheap and
decisive: a **dual-reporter readthrough assay carrying this exact BUB1B UGA-A context**, screened across
a panel of readthrough agents, determines feasibility and agent choice at low cost, and a negative
retires the hypothesis permanently.

**Further derisking.** Readthrough is not selective for premature stops; it also reads through normal
termination codons, and that burden must be measured. NMD interacts directly: an NMD-competent
transcript may be degraded before it can be read through, so the achievable ceiling depends on
transcript abundance, not only on readthrough efficiency.

---

## 4. What this genome argues *against*

Repurposing analyses propose. In a child with a caretaker disease, **the drugs not given matter as much
as the drugs given**, and this section is the half most submissions will omit.

### 4.1 The aneuploidy-suppression playbook targets the opposite lesion

This is the least obvious de-selection here and it follows directly from §1.2.

In CIN cancer cells, kinetochore-microtubule attachments are **inherently more stable** than in normal
diploid cells, errors therefore escape correction and persist into anaphase, and artificially increasing
attachment stability in normal cells is by itself sufficient to cause missegregation (Bakhoum, Genovese
& Compton 2009, *Current Biology*, [DOI](https://doi.org/10.1016/j.cub.2009.09.055)). Aged human
fibroblasts show the same direction, and missegregation there **fell** when MT-depolymerizing kinesin-13
activity was enhanced genetically and with a small molecule, which also delayed senescence
(Barroso-Vilares, Compton & Logarinho 2020, *EMBO Reports*,
[DOI](https://doi.org/10.15252/embr.201949248)).

**This child sits at the other end of the same axis.** Reduced BubR1 means reduced PP2A-B56 and
excessive Aurora B, which elevates Ndc80 phosphorylation and makes attachments *less* stable: BubR1
acetylation is required "to counteract excessive Aurora B activity at the KT **for stable**
chromosome-spindle attachments" (Park et al. 2013). Both states missegregate, by opposite routes — one
fails to correct errors, the other fails to hold correct attachments.

So **kinesin-13 / MCAK enhancement, and the small-molecule strategy validated in ageing cells, is
contraindicated here**, because it pushes the dial the lesion has already pushed too far. An analysis
matching on the phrase "reduces chromosomal instability" rather than on **mechanism direction** would
have proposed exactly this, and been wrong.

### 4.2 Microtubule-targeting agents, with the counterexample stated

Paclitaxel and vincristine act by prolonging spindle checkpoint activation and delaying mitotic exit.
The checkpoint they prolong is the one BubR1 operates. **Inhibition of BubR1 promotes resistance to
microtubule inhibitors in vitro** (Swanton, Tomlinson & Downward 2006, *Cell Cycle* — a **review**, so
the primary data sit in its references, [DOI](https://doi.org/10.4161/cc.5.8.2682)), and independently,
suppressing BubR1 function confers resistance to anti-microtubule drug-induced apoptosis in a manner
**reversed by BubR1 overexpression** ([DOI](https://doi.org/10.1007/s10549-005-6938-0)).

A patent landscape search adds a commercial line of evidence. Two unrelated owners a decade apart on
unrelated chemistry — Aventis on taxanes (`EP 1831398 B1`, expired) and Basilea on
**BAL27862 / BAL101553 (lisavanbulin)** (`WO 2012/098207`) — **independently claimed that BubR1 level
determines response to their microtubule agent**, the Basilea claim specifying that low BUBR1 predicts
*resistance*. Convergent commercial investment on one biomarker, from different chemistry and different
decades, is market-grade corroboration that **BubR1 abundance is the pharmacologically meaningful
variable** — which is also the premise §3.1 depends on.

**The counterexample, stated because it is the only human evidence that exists.** The one published
account of chemotherapy in an MVA patient used VAC. That patient started at the usual infant dose,
developed thrombocytopenia after cycle 2 requiring platelet transfusions, had vincristine and
cyclophosphamide reduced by 50%, **completed eight cycles, and the tumour was still shrinking six months
after stopping.** A vincristine-containing regimen worked in the only relevant human case on record.

So the defensible statement is narrow: **not "avoid vincristine", but that microtubule agents are the
one drug class whose efficacy in this child is predicted to be altered by his own causal gene**, which
makes response monitoring in any future regimen a genome-informed decision rather than a routine one.

### 4.3 Three further de-selections

- **Aminoglycosides**, vetoed on nephrotoxicity against existing nephrocalcinosis (§3.5).
- **Aneuploidy-selective agents face a window problem specific to this child.** Five KIF18A inhibitors
  are recruiting, two listing chromosomal instability as an enrolment condition. Their selectivity is a
  contrast between unstable tumour and stable normal tissue. **This child has no chromosomally stable
  tissue**, so the window narrows from tumour-versus-normal-body to tumour-versus-already-unstable-body.
  Nobody has measured that difference, and until someone does, an otherwise perfectly-indicated class
  cannot be proposed.
- **MPS1 inhibition is the wrong direction.** `EP 3645535 B1` covers an MPS1 inhibitor for "highly
  aneuploid tumors and tumors which overexpress mitotic checkpoint components" — further checkpoint
  weakening in a child whose checkpoint is already hypomorphic.

### 4.4 One de-selection with an open counter-question

**4-phenylbutyrate is deprioritized within its own chaperone class** because it is also an HDAC
inhibitor and HDAC inhibition perturbs mitotic progression — a specific hazard, not a generic side
effect, in a child whose pathology *is* chromosome segregation.

**But BubR1 is deacetylated by HDAC11** (Watanabe et al. 2014, *Cell Reports*,
[DOI](https://doi.org/10.1016/j.celrep.2014.03.022)), and §1.3 establishes that acetylation *stabilizes*
BubR1. If so, HDAC inhibition could raise BubR1 — making 4-PBA's HDAC activity a **feature** rather than
the hazard treated here. The counter is that 4-PBA is a weak class I/II inhibitor while HDAC11 is class
IV, so the relevant selectivity is unknown. **Left open rather than resolved**, because the experiment
is small and the answer either rehabilitates a cheap approved drug or confirms its exclusion.

### 4.5 What is cleared

From the pharmacogenomic screen of this genome: reference, and therefore normal function, at all four
*DPYD* risk alleles, all three *TPMT* alleles and *NUDT15*, *G6PD* (hemizygous in a male, so a single
variant would be fully penetrant), and *SLCO1B1\*5*. **Clean thiopurine and fluoropyrimidine status in a
child whose syndrome predisposes to leukaemia is usable information today**, before any hypothesis above
is tested. Caveats in §8.

---

## 5. Where the field actually is

**Does each hypothesis meet this track's goal line?** The goal asks for existing, approved medications
that are plausible candidates for further investigation. Answered before the argument:

| | Approved medicine? | Clears this child's own picture? |
|---|---|---|
| **Nicotinamide**, abundance | **Yes**, 1982, oral + IV, paediatric | **Yes**, with a defined gate on platelets and phosphate |
| **MEK**, tumour-restricted | **Yes**, paediatric, mature safety | Yes on the drug; gated on evidence |
| **Chaperones** (arimoclomol, UDCA, betaine) | **Yes**, all dosed in children | Partly — no BubR1-specific evidence |
| **PP2A** (fingolimod) | **Yes**, from age 10 | **Probably not** — immunosuppression |
| **Aurora B** | **No** — nine agents, max Phase 3 | Best mechanistic evidence, no drug |
| **Readthrough** | **No longer** — ataluren not renewed | **No** — aminoglycosides vetoed |

**Read the two columns against each other. They run in opposite directions.** The axis with the best
biology has no medicine; the medicines with the best availability have the weakest link to this protein.

That is not a flaw in the analysis. It is the condition of drug repurposing in an ultra-rare disease, and
a submission that conceals it is hiding the thing a clinician most needs to know. **The one candidate
that breaks the pattern is nicotinamide** — approved, cheap, paediatric, with a BubR1-specific mechanism
and Phase 3 human chemoprevention data — which is why it leads.

The practical consequence is a staged protocol rather than a single pick:

- **Now, zero toxicity.** The pharmacogenomic clearances in §4.5, and a surveillance schedule.
- **Cheap and decisive next.** The §3.1 dose-window experiment and the §3.3 dual readout. Weeks, not
  years, and either can retire its hypothesis.
- **Gated.** MEK, on a tumour tissue result.
- **Blocked on drug availability, not on biology.** Aurora B.

---

## 6. What would falsify this

- BubR1 protein is **not** reduced in patient cells → the hypomorph model is wrong and §3.1 and §3.3
  lose their rationale.
- PP2A-B56 kinetochore recruitment is normal in patient cells → links 2–4 fail.
- `p.Asn1002Lys` protein is **not** destabilized → §3.3 is treating a lesion that is not there.
- Nicotinamide raises BubR1 in no dose window, or only at doses where sirtuin inhibition dominates →
  §3.1 collapses to nicotinic acid or to nothing.
- The two variants are **in cis** → the entire genetic model collapses. Untestable with this data
  (10,911 bp apart, 2×149 bp reads, no parental samples) and **the single largest unverified assumption
  in this submission**.
- Tumour tissue shows no *LZTR1* second hit → §3.2 closes.

---

## 7. What transfers beyond this case

*Scalability: the methods, not the drug.*

- **The `INFO/AF` trap.** In a GATK VCF this field is per-sample allele fraction, not population
  frequency. Filtering rare variants against it discards nearly every true heterozygous call. Applies to
  any rare-disease analysis handed a HaplotypeCaller VCF.
- **A symmetry discriminator for copy number versus mismapping**, from allele balance alone, no
  reference panel. A mosaic gain widens the heterozygous distribution symmetrically; mismapping is
  one-sided low. Testing the two tails separately turns an ambiguous depth signal into a decision, with
  chrX as a built-in positive control in any male sample.
- **AVI requires a frequency pairing.** AlphaGenome's AVI is frequency-naive by design. Here a variant
  8 kb from a candidate gene scored 1.149 as a genuine splice-acceptor variant and is carried by 52.7%
  of the population. Used alone it generates confident false leads next to the right gene.
- **Mechanism direction, not disease label, decides transferability.** §4.1 is the general lesson: two
  diseases sharing the phrase "chromosomal instability" can sit at opposite ends of the same axis, and a
  therapy validated for one can be contraindicated in the other. Any repurposing pipeline that matches
  on phenotype vocabulary rather than on the sign of the perturbation will make this error.
- **The pseudokinase-as-chaperone-target pattern.** Pseudokinase pockets are structural rather than
  catalytic, so any destabilizing missense in any pseudokinase is a candidate for the CETSA-plus-function
  screen in §3.3. A template, not a one-off.
- **Stop-codon-context-first readthrough triage.** Read the PTC identity and +4 context off the
  reference, then choose the assay and agent panel, rather than assuming a universal hierarchy.

---

## 8. Limitations

1. **Phase is inferred, not observed**, and cannot be observed with this data. The largest unverified
   assumption here.
2. **No tumour tissue**, so the *LZTR1* somatic hit is untestable, which is why §3.2 is gated.
3. **No approved Aurora B inhibitor at any phase**, which is why the axis with the best evidence is a
   keystone rather than a candidate.
4. **The BubR1-specific evidence is about abundance, not folding.** §1.3 and §3.1 show BubR1 *levels*
   can be raised in vivo; they do not show a destabilized `p.Asn1002Lys` fold can be rescued. Adjacent
   hypotheses, not the same one.
5. **The in vivo abundance work is mouse, the readout was lifespan in a progeroid model, and the SIRT2
   arm was genetic.** NMN was the only pharmacological arm.
6. **Aneuploidy was not directly observed** in this blood sample; the scan was negative with a ~15%
   mosaic floor and does not substitute for a karyotype in cultured cells.
7. **Pharmacogenomics is a single-variant screen, not star-allele diplotyping.** *CYP2D6* cannot be
   called reliably from short reads. Absence of a record means homozygous reference *only if* the
   position was covered, unconfirmable without a BAM.
8. ***TRIP13* is unassessed** owing to subtelomeric coverage failure.
9. **This is a hypothesis document.** Nothing here has been tested in this child, and nothing should be
   acted on without the experiments named alongside each candidate.

---

## 9. Open questions

*Generated by this analysis and not resolved by it. Recorded because the questions a piece of work opens
are part of its output.*

**On the biology**

- **Why is the MVA tumour spectrum narrow when the instability is global?** Random missegregation at the
  cell level produces a recurrent, limited tumour set at the organism level. What selects it? Chaos at
  one scale is evidently order at another, and the selection rule is unknown.
- **Is BubR1 level threshold-like or graded?** The entire dose-rescue logic of §3.1 and §3.3 depends on
  this, and a threshold would define how much recovery is "enough."
- **How has the immune system coexisted with constitutional chromosomal instability without clearing
  aneuploid cells?** The likeliest answer is that instability has been present since embryogenesis, so
  central tolerance had the whole developmental window to accommodate it. If true, immune approaches
  here differ fundamentally from immune approaches to somatic CIN in cancer, and that is testable.
- **Is there a point in time at which the cells were euploid?** Constitutional aneuploidy cannot mean
  every cell at every moment. Where the first error occurs, and whether a euploid progenitor compartment
  persists, bears directly on whether any intervention has something to protect.
- **Does meiosis fail the same way?** Carrier parents show premature chromatid separation and the family
  history includes recurrent pregnancy loss, which suggests the meiotic divisions carry the same lesion.
  Whether meiotic checkpoints differ enough to be separately targetable is unexamined.

**On the therapy**

- **Does an aneuploidy-selective agent have a therapeutic window in a patient with no chromosomally
  stable tissue?** Stated in §4.3 as a de-selection; stated here as the quantitative question that would
  resolve it, which is answerable by measuring the difference in instability between tumour and normal
  tissue in the same patient.
- **Does HDAC11 inhibition raise BubR1?** §4.4. The answer either rehabilitates a cheap approved drug or
  confirms its exclusion.
- **Is the nephrocalcinosis phosphate-driven or calcium-driven?** It determines whether §3.1's phosphate
  effect is a benefit or a hazard, and it is answerable from routine bloods.
- **Could a ligand restore B56-PP2A at the kinetochore without global phosphatase activation?** §3.4
  ends in a specification for a molecule that does not exist. Writing it down is the first step.

**On the patent landscape**

- **BubR1 is essentially unclaimed as a therapeutic target.** A Lens.org search across 306 patent records
  found **no patent claiming a chaperone or proteostasis rescue of BubR1** — all eleven hits at that
  intersection are incidental gene-panel or PROTAC-target lists. Seoul National University's
  acetylated-BubR1 screening method patent (`US 2008/0102477`) is **discontinued**, so it is prior art
  rather than an obstacle; only a mouse-specific reagent patent survives.
- **One freedom-to-operate item remains open.** Basilea's `WO 2012/098207` contains a claim to
  "increasing the level of BUBR1 in a subject that has a sample with a lower level of BUBR1" (US member
  `US 10222377`, adjusted expiry 2032-01-25). The granted US claims that could be read are
  measure/determine/kit claims, so the broad rescue claim most likely died in prosecution, **but the
  complete granted claim set should be read before anyone asserts freedom to operate on raising BubR1.**

**Validated target, open mechanism.** That combination is unusual, and it is the structural reason this
case is worth more work than one child's genome would normally justify.

---

## 10. Methods description

*Responses to the provided Track 2 methods template.*

**Team name.** `Mat34ec` (individual participant).

**Approach, from variant to candidate.** Variants were recovered by the Track 1 pipeline. Mechanism was
characterized by combining AVI feature attribution (to establish the impact is protein-coding rather
than regulatory), two independent model families (to exclude splicing), and primary literature for the
pathway. Candidates were then derived **from the characterized mechanism**: each one names the specific
link in §1.2 it acts on, and was then tested against this child's own clinical and genomic picture
(§2, filter 3). Agents were identified through ChEMBL (target records, mechanism of action, clinical
phase, CETSA bioactivity), ClinicalTrials.gov (clinical stage, paediatric precedent), PubMed, and
Lens.org (patent landscape and freedom-to-operate).

**Generative AI disclosure (required).** Anthropic Claude Pro and Claude Science, Claude Opus 5,
non-commercial terms, no training on customer content. Google Gemini 3.1 Pro, non-commercial terms, no
training on customer content. Perplexity (free tier) as a search engine, non-commercial, no training on
customer content.

**The division of labour between the human and the models, since the track asks.** The models did
retrieval at a scale a single person cannot, wrote the analysis code, and were used deliberately as
adversaries against the author's own hypotheses. The author set the direction, chose the decision
criteria, supplied the biological framing, and made every call about what survived. **Every quantitative
claim in this report derives from a script in the repository or from a named primary source; no figure
is reported from model recall.** Where a model-supplied citation could not be verified against the
primary source, it was removed — three such cases are recorded below.

**Automated or manual candidate identification?** Both, in sequence. Database queries were automated;
mechanism characterization, hypothesis generation, the ordering in §3 and every objection raised against
the author's own candidates are manual expert reasoning.

**Manual literature review and curation.** Extensive and deliberate, and it changed conclusions three
times. A proposed LZTR1→KARD→PP2A→Aurora B connection was **discarded**: it does not exist, and it would
have been the most damaging error in the submission. A single-locus pharmacogenomic result suggesting a
vincristine-neuropathy risk genotype **reversed** when the full panel was run — of five reported loci the
proband carries only the weakest and most contested, and is reference at the two with the largest
reported effects. A cited PDB structure proved to be **an unrelated bacterial enzyme**, and the real
BUB1B entries were substituted (`2WVI`, `3SI5`, `4GGD`, `5JJA`, `5K6S`, `5KHU`, `5LCW`, `5SWF`, `6TLJ`,
UniProt O60566). These are recorded rather than quietly removed, because a method that never reports its
own failures gives a reviewer nothing to calibrate against.

**Public data only?** Yes, exclusively. No proprietary sources.

**Public sources used.** ChEMBL, ClinicalTrials.gov, PubMed, Lens.org, Ensembl REST and VEP, gnomAD,
ClinVar, GENCODE, AlphaGenome and AlphaGenome Atlas, SpliceAI v1.3.1, EMA public assessment reports.

**How was the mechanism characterized?** See §1. Loss of function, partial, hypomorphic. Pathway: BubR1
→ KARD phosphorylation → PP2A-B56 recruitment → opposition of Aurora B at the kinetochore. Downstream
consequence: premature chromatid separation and mosaic aneuploidy. Disease class: caretaker rather than
gatekeeper, which determines that the fitting therapeutic category is chemoprevention.

**Effort.** Roughly two months part-time, of which the majority was validation and self-refutation rather
than hypothesis generation.

---

## 11. Sources

Xu P, Raetz EA, Kitagawa M, Virshup DM, Lee SH (2013). BUBR1 recruits PP2A via the B56 family of
targeting subunits to promote chromosome congression. *Biology Open* 2(5):479-86. PMID 23789096.
[DOI](https://doi.org/10.1242/bio.20134051)

Choi E, Choe H, Min J, Choi JY, Kim J, Lee H (2009). BubR1 acetylation at prometaphase is required for
modulating APC/C activity and timing of mitosis. *EMBO J* 28(14):2077-89. PMID 19407811.
[DOI](https://doi.org/10.1038/emboj.2009.123)

Park I, et al. (2013). Loss of BubR1 acetylation causes defects in spindle assembly checkpoint signaling
and promotes tumor formation. *J Cell Biol* 202(2):295-309. PMID 23878276.
[DOI](https://doi.org/10.1083/jcb.201210099)

North BJ, et al. (2014). SIRT2 induces the checkpoint kinase BubR1 to increase lifespan. *EMBO J*
33(13):1438-53. PMID 24825348. [DOI](https://doi.org/10.15252/embj.201386907)

Chen AC, et al. (2015). A phase 3 randomized trial of nicotinamide for skin-cancer chemoprevention.
*NEJM* 373(17):1618-26. PMID 26488693. [DOI](https://doi.org/10.1056/NEJMoa1506197)

Ketteler M, et al. (2021). Efficacy and safety of a novel nicotinamide modified-release formulation in
refractory hyperphosphatemia in haemodialysis. *Kidney Int Rep* 6(3):594-604. PMID 33732974.
[DOI](https://doi.org/10.1016/j.ekir.2020.12.012)

Shahbazian H, et al. (2011). Oral nicotinamide reduces serum phosphorus, increases HDL, and induces
thrombocytopenia in hemodialysis patients. *Nefrologia* 31(1):58-65. PMID 21270914.
[DOI](https://doi.org/10.3265/Nefrologia.pre2010.Nov.10734)

Ix JH, et al. (2019). Effects of nicotinamide and lanthanum carbonate on serum phosphate and FGF-23 in
CKD: the COMBINE trial. *JASN* 30(6):1096-1108. PMID 31085679.
[DOI](https://doi.org/10.1681/ASN.2018101058)

Bakhoum SF, Genovese G, Compton DA (2009). Deviant kinetochore microtubule dynamics underlie chromosomal
instability. *Current Biology* 19(22):1937-42. PMID 19879145.
[DOI](https://doi.org/10.1016/j.cub.2009.09.055)

Barroso-Vilares M, Macedo JC, Reis M, Warren JD, Compton D, Logarinho E (2020). Small-molecule
inhibition of aging-associated chromosomal instability delays cellular senescence. *EMBO Reports*
21(5):e49248. PMID 32134180. [DOI](https://doi.org/10.15252/embr.201949248)

Cohen-Sharir Y, et al. (2021). Aneuploidy renders cancer cells vulnerable to mitotic checkpoint
inhibition. *Nature*. [DOI](https://doi.org/10.1038/s41586-020-03114-6)

Yost S, et al. (2017). Biallelic TRIP13 mutations predispose to Wilms tumor and chromosome
missegregation. *Nature Genetics*. [DOI](https://doi.org/10.1038/ng.3883)

Matsuura S, et al. (2006). Monoallelic BUB1B mutations and defective mitotic-spindle checkpoint in seven
families with premature chromatid separation syndrome. *Am J Med Genet A*.
[DOI](https://doi.org/10.1002/ajmg.a.31069)

García-Castillo H, et al. (2008). Clinical and genetic heterogeneity in patients with mosaic variegated
aneuploidy. *Am J Med Genet A*. [DOI](https://doi.org/10.1002/ajmg.a.32315)

Nozaki T, et al. (2015). A novel in vitro model of sarcopenia using BubR1 hypomorphic C2C12 myoblasts.
*Cytotechnology* 68(5):1705-15. [DOI](https://doi.org/10.1007/s10616-015-9920-7)

Watanabe Y, Khodosevich K, Monyer H (2014). Dendrite development regulated by the schizophrenia-associated
gene FEZ1 involves the ubiquitin proteasome system. *Cell Reports* 7(2):552-64. PMID 24726361.
[DOI](https://doi.org/10.1016/j.celrep.2014.03.022)

Saddoughi SA, et al. (2013). Sphingosine analogue drug FTY720 targets I2PP2A/SET and mediates lung
tumour suppression via activation of PP2A-RIPK1-dependent necroptosis. *EMBO Mol Med* 5(1):105-21.
PMID 23180565. [DOI](https://doi.org/10.1002/emmm.201201283)

Enjoji S, et al. (2015). The therapeutic effects of SET/I2PP2A inhibitors on canine melanoma.
*J Vet Med Sci* 77(11):1451-6. [DOI](https://doi.org/10.1292/jvms.15-0193)

Oaks J, Ogretmen B (2015). Regulation of PP2A by sphingolipid metabolism and signaling.
*Frontiers in Oncology* 4:388. [DOI](https://doi.org/10.3389/fonc.2014.00388)

Swanton C, Tomlinson I, Downward J (2006). Chromosomal instability, colorectal cancer and taxane
resistance. *Cell Cycle* 5(8):818-23. **Review.** PMID 16628000.
[DOI](https://doi.org/10.4161/cc.5.8.2682)

Inaba S, et al. (2005). Synuclein gamma inhibits the mitotic checkpoint function and promotes
chromosomal instability of breast cancer cells. *Breast Cancer Res Treat* 94(1):25-35. PMID 16142440.
[DOI](https://doi.org/10.1007/s10549-005-6938-0)

Ngiwsara L, et al. (2022). Effects of gentamicin inducing readthrough premature stop codons.
*Biochem Biophys Res Commun*. [DOI](https://doi.org/10.1016/j.bbrc.2022.10.081)

Yu H, et al. (2019). Stop codons and the +4 nucleotide may influence the efficiency of G418 in rescuing
nonsense mutations of the HERG gene. *Int J Mol Med*. [DOI](https://doi.org/10.3892/ijmm.2019.4360)

Schilff M, Sargsyan Y, Hofhuis J, Thoms S (2021). Stop codon context-specific induction of translational
readthrough. *Biomolecules* 11(7):1006. [DOI](https://doi.org/10.3390/biom11071006)

Dabrowski M, Bukowy-Bieryllo Z, Zietkiewicz E (2015). Translational readthrough potential of natural
termination codons in eucaryotes. *RNA Biology*. [DOI](https://doi.org/10.1080/15476286.2015.1068497)

Tabe S, et al. (2026). Spindle assembly checkpoint competency determines sensitivity to KIF18A
inhibition in small cell lung cancer. *Cancer Research Communications*.
[DOI](https://doi.org/10.1158/2767-9764.CRC-26-0560)

**Database and regulatory records.** ChEMBL targets `CHEMBL4295998` (BUB1B), `CHEMBL2185` (Aurora kinase
B), `CHEMBL4523403` (KIF18A); CETSA assay `CHEMBL6125447`, document `CHEMBL6121040`; compounds
`CHEMBL1140` (nicotinamide), `CHEMBL2107726` (arimoclomol), `CHEMBL1746` / `CHEMBL2105745`
(phenylbutyrate forms), `CHEMBL256997` (ataluren). EMA public assessment reports and news releases on
the non-renewal of Translarna's marketing authorisation. Lens.org patent records `WO 2012/098207`,
`US 10222377`, `EP 1831398 B1`, `US 2008/0102477`, `US 11186633`, `EP 3645535 B1`. ClinicalTrials.gov
NCT06084416, NCT05902988, NCT06799065, NCT06772415, NCT07226427 (KIF18A); NCT03155620, NCT04924608,
NCT05849662, NCT07817186 (MEK and RASopathy).

---

*Released under CC BY 4.0 in accordance with the hackathon rules.*
