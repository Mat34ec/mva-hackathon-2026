# Track 2: Drug Repurposing

**Rare Disease, Real Kid: The MVA Hackathon 2026**
Participant `Mat34ec` · GRCh38

> These are hypotheses for follow up, not evidence that any medicine works, and not medical advice. Every candidate is stated with what must be true for it to work, and with the experiment that would show it does not.

## The answer, up front

**The two *BUB1B* alleles found in this genome make a hypomorph, and that one word decides everything else.** One null allele, one destabilizing missense allele in the pseudokinase domain. Complete BubR1 loss is embryonic lethal, so a living patient means reduced protein rather than absent protein. There is a quantity to move, not an absence to replace. That is what makes the disease addressable at all.

The quantity has a named biochemistry running in one direction. BubR1 recruits PP2A-B56 to the kinetochore. PP2A-B56 opposes Aurora B. Less BubR1 means less PP2A-B56, so Aurora B runs unopposed and attachments that should have been kept are torn down along with the bad ones. **The defect is not a missing function. It is a lost brake on a function that is supposed to be there.** Restoring a balance is tractable in a way that replacing a protein is not.

Working from that mechanism rather than from the disease name produces three things:

- **A lead candidate that is approved, cheap, paediatric, and carries Phase 3 human cancer prevention evidence at the dose in question.** Nicotinamide, reaching BubR1 through NAD⁺ and SIRT2, in the right direction, demonstrated in the mouse that carries this genotype.
- **A list of drugs that must not be given, including a class the obvious analysis would have recommended.** The aneuploidy suppression playbook targets the opposite end of the same axis.
- **One selectivity axis in a genome that has almost none**, since every cell here is abnormal and there is no tumour versus normal contrast left to exploit. The monoallelic *LZTR1* allele is the exception.

**The mechanism map is the deliverable. The candidates are what it predicts.** If the map is right it keeps producing hypotheses after this submission closes, which is why section 1 is the longest part of this report and section 2 is spent attacking it.

## Contents

| | |
|---|---|
| [1. The mechanism, link by link](#1-the-mechanism-link-by-link) | Genotype to phenotype, and the targets it exposes |
| [2. Why this mechanism holds up](#2-why-this-mechanism-holds-up) | Three kinds of evidence converging on one lesion |
| [3. The targets, one by one](#3-the-targets-one-by-one) | Nicotinamide leading, then four more with their gates |
| [4. What this genome argues against](#4-what-this-genome-argues-against) | Including the therapy the obvious analysis would pick |
| [5. Where the field actually is](#5-where-the-field-actually-is) | Evidence quality and drug availability run opposite |
| [6. What would falsify this](#6-what-would-falsify-this) | Six results that would end it |
| [7. What transfers beyond this case](#7-what-transfers-beyond-this-case) | Six methods usable on other genomes |
| [8. Limitations](#8-limitations) | What this cannot establish |
| [9. Open questions](#9-open-questions) | The ones worth someone else's career |
| [10. Methods description](#10-methods-description) | Template responses |
| [11. Sources](#11-sources) | Primary literature, databases, regulatory records |

## 1. The mechanism, link by link

*Required by the track: loss or gain of function, pathway disrupted, downstream consequence.*

### 1.1 The two alleles

| | Variant | Class | Consequence |
|---|---|---|---|
| **Allele 1** | `chr15:40209701 T>G` · `NM_001211.6:c.2210T>G` · `p.Leu737Ter` | Loss of function, complete | Reference codon `TTA`; `T>G` at base 2 gives `TGA`, which is UGA, +4 context A. More than 50 nt upstream of the last exon junction, therefore NMD competent. Functionally a null. |
| **Allele 2** | `chr15:40220612 T>G` · `c.3006T>G` · `p.Asn1002Lys` | Loss of function, partial | Missense in the pseudokinase domain. Asparagine to lysine: neutral polar replaced by longer and positively charged. |

The stop codon was read off the reference sequence rather than inferred from the HGVS string, because codon identity determines whether one of the candidate strategies is plausible at all (section 3.5).

**The lesion is a hypomorph, and the constraint is clinical rather than biochemical.** Complete BubR1 loss is embryonic lethal, and he is alive, which fixes the model independently of any assay. Allele 1 contributes no protein. Allele 2 contributes protein of reduced competence. **Allele 2 is therefore the only rescuable allele in this genome**, and total functional BubR1 is the variable every candidate below is trying to move.

### 1.2 The chain

**Link 1. BubR1 is a pseudokinase.** Its C terminal domain adopts a kinase fold, but catalysis is dispensable and the nucleotide binding residues serve conformational stability instead. The protein is a scaffold, not an enzyme.

This decides the shape of every intervention downstream: on a scaffold the mechanistically correct ligand is a **stabilizer, not an inhibitor**, because there is no catalysis to inhibit.

**Link 2. The KARD motif is phosphorylated**, licensing recruitment of PP2A-B56 to the kinetochore. This is the step that converts BubR1 abundance into phosphatase presence at a specific location.

**Link 3. PP2A-B56 exists to oppose Aurora B.** Aurora B phosphorylates kinetochore substrates to loosen microtubule attachments, which is how a cell discards incorrect ones. PP2A-B56 dephosphorylates the same substrates to lock in correct, tension bearing attachments. The kinetochore needs both, in balance.

**Link 4. Less BubR1 means less PP2A-B56, so Aurora B runs unopposed.** Ndc80 phosphorylation rises, and attachments that should have been kept are torn down along with the bad ones.

**Link 5. Chromosomes fail to congress and sister chromatids separate prematurely.** Daughter cells inherit wrong chromosome counts. Because the failure is probabilistic and recurs at every division, different cells acquire different aneuploidies. Hence mosaic variegated: a patchwork rather than one uniform trisomy.

### 1.3 What the chain exposes

Each link is a place a drug could act. Reading the chain as a target list gives five, and the direction required at each is fixed by the link itself:

| Link | Target | Required direction | Section |
|---|---|---|---|
| Allele 2 protein level | BubR1 abundance | **Raise** | 3.1 |
| Allele 2 protein fold | Pseudokinase domain stability | **Stabilize**, never inhibit | 3.3 |
| Link 2 to 3 | PP2A-B56 at the kinetochore | **Restore recruitment** | 3.4 |
| Link 3 | Aurora B | **Oppose** | 3.4 |
| Allele 1 | Premature stop codon | **Read through** | 3.5 |

The directions matter more than the targets. Four of the five require increasing or restoring something, and mitotic drug development has overwhelmingly built inhibitors. Section 4 is what happens when the direction is ignored.

### 1.4 Caretaker, not gatekeeper

MVA is a **caretaker** defect: it raises the rate of missegregation and mutation rather than directly driving proliferation. A gatekeeper lesion removes a brake on growth; a caretaker lesion corrupts the machinery that keeps the genome accurate.

The consequence governs this whole report. **Correcting the caretaker defect does not treat an existing tumour, and treating a tumour does not correct the caretaker defect.** Two problems, kept apart throughout.

It also points at the therapeutic category that fits a caretaker disease, which is **chemoprevention** rather than cure. That is where the lead candidate turns out to live, and the category was derived from the disease class before any drug was considered.

Checkpoint severity tracks tumour risk across the MVA genes. Biallelic *TRIP13* or *BUB1B* patients show severe checkpoint impairment and high embryonal tumour risk, while *CEP57* and unexplained MVA show minimal deficiency and no embryonal tumours (Yost et al. 2017, *Nature Genetics*, [DOI](https://doi.org/10.1038/ng.3883)). He falls in the high risk group, which is why the oncological axis cannot be dropped even though the caretaker axis is what explains the disease.

### 1.5 The second locus, and what it is actually for

*LZTR1* `p.Tyr748Ter` is monoallelic. LZTR1 is a CUL3 substrate adaptor that ubiquitinates RAS for degradation, and monoallelic inactivating alleles are the schwannomatosis genotype, requiring a somatic second hit under two hit logic.

**It is not downstream of *BUB1B*.** There is no signalling link between the spindle checkpoint and the RAS MAPK cascade and this report claims none. The only connection is mutational: chromosomal instability plausibly supplies the somatic second hit a monoallelic *LZTR1* allele needs. Instability raises the probability of the hit; it does not transmit a signal.

Its value stated precisely. Not that it explains the tumour better than *BUB1B* does, but that **if the somatic hit is real, RAS MAPK activation exists only in tumour cells.**

Ordinary oncology leans on a contrast between abnormal tumour and normal body, and a constitutional instability syndrome has deleted that contrast everywhere else. **This locus is the only surviving selectivity axis in this genome**, which is why section 3.2 exists and why the gate on it is worth a tissue test.

## 2. Why this mechanism holds up

A mechanism map is only worth building if it can be attacked. Three independent kinds of evidence converge on the same lesion, none assembled to agree with the others.

**Biological, and in patient cells.** Links 3 and 4 carry the weight of the model, and the support is specific rather than inferential. The BUBR1:B56-PP2A complex opposes Aurora B kinase activity, **loss of the complex can be reverted by inhibiting Aurora B**, and failure to recruit B56-PP2A **contributes to the chromosome congression defects found in cells derived from patients with MVA syndrome** (Xu et al. 2013, *Biology Open*, PMID 23789096, [DOI](https://doi.org/10.1242/bio.20134051)). Three things in one result: the complex opposes Aurora B, the defect is reversible, and it was seen in MVA patient derived cells rather than an engineered line. Reversibility is the part that matters therapeutically, because it makes the lesion correctable rather than merely described.

**Genetic, from a literature that was not looking for this.** The chain in section 1.2 came from the checkpoint literature. A separate literature on how BubR1 abundance is controlled arrives at the same place. BubR1 is acetylated at **Lys250** by PCAF, and that acetylation **blocks ubiquitin dependent BubR1 degradation** by APC/C-Cdc20: the acetyl mimetic K250Q is stabilized, the acetyl deficient K250R is unstable (Choi et al. 2009, *EMBO J*, [DOI](https://doi.org/10.1038/emboj.2009.123)). Mice heterozygous for acetylation deficient BubR1 **spontaneously develop tumours with massive chromosome missegregation**, and decisively, "the association of acetylation deficient BubR1 with **PP2A-B56α** phosphatase was reduced" with elevated phosphorylated Ndc80 at the kinetochore (Park et al. 2013, *J Cell Biol*, [DOI](https://doi.org/10.1083/jcb.201210099)).

Two independent routes to one lesion. Reduce BubR1 by mutation, as here, or reduce it by blocking its acetylation, as in those mice, and the result is the same biochemistry and the same phenotype. That convergence raises confidence the map describes reality rather than a plausible story, and it establishes the principle that organizes section 3: **BubR1 abundance is the pharmacologically meaningful variable.** Not its activity, which it does not have.

**Commercial, which is the unusual one.** Two unrelated owners a decade apart on unrelated chemistry both staked patent claims on BubR1 level predicting response to their drug: Aventis on taxanes (`EP 1831398 B1`, expired) and Basilea on BAL27862 / BAL101553 (lisavanbulin) (`WO 2012/098207`), the Basilea claim specifying that low BUBR1 predicts resistance. Companies do not spend prosecution money on a biomarker they think is noise. Convergent investment on one variable, from different chemistry in different decades, corroborates from the market what the biology already says.

## 3. The targets, one by one

Candidates were not found by searching for drugs associated with MVA, because there are none. Each was found by asking whether an approved medicine moves a named link in section 1.2 **in the required direction**, then asking what would break it in this patient specifically: his kidneys, his pharmacogenomics, his cancer predisposition, and the single published treatment record that exists in this disease. That last filter did the most work, moving the original favourite out of the lead position and producing section 4. Where no approved medicine exists for an axis, it is reported as a mechanistic finding rather than dressed up as a candidate, and two axes end that way.

### 3.1 Nicotinamide, to raise BubR1 abundance. **Lead candidate**

Strong from every direction at once: mechanism, human evidence, availability, cost, and monitoring that costs nothing extra.

**The mechanism, and why the direction is right.** BubR1 abundance falls as NAD⁺ declines, through SIRT2 acting at **Lys668**, counteracted by the acetyltransferase CBP. **Overexpression of SIRT2, or treatment with the NAD⁺ precursor NMN, increases BubR1 abundance in vivo**, and SIRT2 overexpression in **BubR1^H/H hypomorphic mice** increases median lifespan (North et al. 2014, *EMBO J*, [DOI](https://doi.org/10.15252/embj.201386907)).

That mouse is the right mouse: hypomorphic rather than knockout, which is this genotype expressed in an animal. The intervention has already been run against the correct model organism rather than an approximation.

**Phase 3 human evidence, at this dose, in a cancer prone population.** ONTRAC randomised 386 participants with at least two prior nonmelanoma skin cancers, double blind, to **nicotinamide 500 mg twice daily** or placebo for 12 months. New cancers fell by **23%** (95% CI 4 to 38, p=0.02), actinic keratoses by 13 to 20%, and there were **"no noteworthy between-group differences... with respect to the number or types of adverse events"** (Chen et al. 2015, *NEJM*, [DOI](https://doi.org/10.1056/NEJMoa1506197)).

Section 1.4 derived chemoprevention as the category fitting a caretaker disease before any drug was on the table. Chemoprevention with this agent then turns out to have a Phase 3 precedent in people already prone to cancer. Category and evidence were reached separately and met.

**It may also act on the one phenotype *BUB1B* does not explain.** Nephrocalcinosis is calcium phosphate deposition, and nicotinamide lowers serum phosphate, the calcium phosphorus product and iPTH in randomised trials ([DOI](https://doi.org/10.1016/j.ekir.2020.12.012); [DOI](https://doi.org/10.3265/Nefrologia.pre2010.Nov.10734)), though COMBINE found no significant effect in CKD ([DOI](https://doi.org/10.1681/ASN.2018101058)), so the effect is real but not reliable. Not proposed as a treatment for it. Flagged because a drug selected on mitotic grounds landing on the one unexplained phenotype is a coincidence worth testing rather than quietly using as support.

**Scale and availability are the strongest practical argument in this report.** Approved 1982, oral and parenteral, ATC A11HA01, no black box warning, decades of paediatric use, manufactured globally as a vitamin at commodity prices with no cold chain and no specialist administration.

That matters beyond one patient. MVA has fewer than a hundred reported cases worldwide, scattered across very different health systems, and not every family reaches a tertiary centre. **An intervention costing pennies and already on formulary almost everywhere is the only kind that reaches all of them.** A brilliant candidate reaching one family in ten is worth less here than a good one reaching ten in ten.

**Three things must be true and are not yet known.** Each is followed by what resolves it.

**(a) The sirtuin tension, the sharpest objection to the whole idea.** The published mechanism requires SIRT2 **activity**, and nicotinamide is an NAD⁺ precursor *and* a sirtuin inhibitor by product inhibition, so above some concentration it opposes the enzyme the mechanism depends on. A dose window question, unresolved.

*Resolved by:* a dose response measurement of BubR1 protein in patient derived cells, head to head against **nicotinic acid**, which reaches NAD⁺ by the Preiss-Handler pathway without that liability. One experiment, two agents, a few weeks. It finds the window, moves the proposal to nicotinic acid, or closes the axis, and all three outcomes are useful.

**(b) Thrombocytopenia, which is patient specific rather than generic.** A randomised trial in haemodialysis patients found falling platelet counts at **500 to 1000 mg/day, lower doses than other studies**, with the authors asking for safety work in other populations ([DOI](https://doi.org/10.3265/Nefrologia.pre2010.Nov.10734)). ONTRAC saw no signal at the same dose in a healthier group, so population matters, and it matters here specifically: **the single published account of chemotherapy in an MVA patient was dose limited by thrombocytopenia after cycle 2.**

*Resolved by:* the schedule below, which is the part that costs nothing. A patient with this syndrome is already under cancer surveillance, so bloods are already being drawn. This monitoring rides along on venipuncture that is happening anyway.

| When | What | Why |
|---|---|---|
| Baseline | Platelets, serum phosphate, calcium, PTH, 1,25-vitamin D | The phosphate effect could help or harm depending on which mechanism drives the nephrocalcinosis, and that is unknown |
| Weeks 2 and 4 | Platelets | The haemodialysis signal appeared within weeks, so the early window is the informative one |
| Monthly to month 3 | Platelets, phosphate, calcium | Establishes whether a trend exists before it matters |
| Then quarterly | Platelets, phosphate, calcium, PTH | Aligns with existing tumour surveillance intervals |
| Stop rule, fixed in advance | A confirmed fall below the age specific threshold, or any bleeding symptom | Decided before starting rather than during |

**(c) ONTRAC's protection is probably not the mitotic mechanism.** Nicotinamide's effect against UV driven skin cancer is attributed to NAD⁺ supporting DNA repair and limiting photoimmunosuppression, a different route to cancer than missegregation. ONTRAC establishes safety and the plausibility of chemoprevention with this agent, not prevention of aneuploidy driven tumours.

*Resolved by:* the same experiment as (a), with a second readout. Does the agent raise BubR1 protein **and** restore PP2A-B56 kinetochore recruitment in patient derived cells? Abundance alone is not enough, because link 2 is where abundance becomes function. A positive on both connects ONTRAC's safety record to the actual lesion; a positive on abundance alone means the agent raises a protein that still cannot do its job.

Every one of those three unknowns can kill this candidate, and each has an experiment that settles it in weeks. That is the point: built to be tested quickly and cheaply rather than to survive scrutiny.

### 3.2 MEK inhibition, the tumour restricted axis. **Gated on one test**

**Rationale.** If the somatic second hit at *LZTR1* is real, RAS protein accumulates and MAPK signalling is elevated **in tumour cells only** (section 1.5). MEK inhibition then acts on the one lesion that is not constitutional, the only place in this genome where a selective therapy has somewhere to be selective. Separately, rhabdomyosarcoma is a differentiation blocked muscle tumour and ERK signalling restrains MYOG driven differentiation.

**The differentiation argument has a second leg native to *BUB1B*.** BubR1 insufficiency by itself suppresses myogenic differentiation: BubR1 hypomorphic C2C12 myoblasts show markedly reduced myotube formation with lower **myogenin and MyoD**, and raised Atrogin-1 and MuRF-1 driving proteasomal muscle protein loss (Nozaki et al. 2015, *Cytotechnology*, [DOI](https://doi.org/10.1007/s10616-015-9920-7)). So it does not rest on RAS alone.

**Availability is strong.** Selumetinib is approved for paediatric NF1 plexiform neurofibroma with a mature safety dataset behind it, including a post authorisation safety study and a 409 patient observational cohort, and MEK inhibition in germline RASopathy is itself in trial (NCT07817186). This is a drug paediatric oncology already knows how to give.

**Why it is gated.** The somatic hit is a hypothesis, not an observation. A genome wide subclonal LOH scan of blood was negative with a floor of 12 to 15% tumour fraction, far above the circulating tumour DNA fraction expected from a localized paediatric sarcoma. At the `p.Tyr748Ter` position the read split is 24 reference to 24 mutant at depth 48, dead even. **Untestable from blood rather than excluded.**

> **The gate.** If tumour tissue shows *LZTR1* loss of heterozygosity or a second somatic mutation, this arm opens. If not, it stays closed.

That converts the weakest assumption in the submission into a decision with a defined answer, and the test is a single assay on tissue obtained anyway in the course of managing a tumour.

**The known risks, and how each is managed.** MEK inhibition in children is not free: rash, paronychia, gastrointestinal toxicity, CK elevation, reduced ejection fraction, ocular toxicity, with dose interruptions common in the NF1 experience. None of it is unfamiliar territory and none of it is unmanaged in current practice:

| Risk | How paediatric MEK use already handles it |
|---|---|
| Rash and paronychia | Commonest events; topical therapy and dose interruption; reversible |
| Reduced ejection fraction | Baseline echocardiogram and scheduled surveillance, with defined interruption thresholds |
| Ocular toxicity | Baseline and periodic ophthalmology review |
| CK elevation | Routine monitoring, usually asymptomatic |
| Renal exposure concern | A fixed duration course with renal function monitoring rather than indefinite dosing, decided before starting |

**So this is a known risk picture with established management, on a drug paediatric oncologists use routinely, given only if a tissue test first shows the target exists.** The gate comes first, the duration is fixed in advance, and every listed toxicity already has a surveillance protocol attached. That is the best characterised risk profile of any candidate here, which is an argument for it.

### 3.3 Chemical chaperones, to stabilize the allele 2 fold

**Rationale from first principles.** Allele 1 is a null and cannot be stabilized. Allele 2 produces protein whose pseudokinase domain is destabilized, and since that domain's role is conformational rather than catalytic (link 1), a ligand occupying the nucleotide pocket would act as a **pharmacological chaperone** rather than an inhibitor.

**The fold is demonstrably ligand stabilizable.** A 2025 cellular thermal shift dataset against BUB1B in CAKI-1 cells reports four compounds: two reduce the BubR1 surviving a 45 °C challenge (fold change 0.29 and 0.65) and **two increase it (1.33 and 2.19, annotated "Increase In Heat Resistance Observed")** (ChEMBL assay `CHEMBL6125447`, document `CHEMBL6121040`). Both stabilizers are preclinical. The only approved molecule with measured BubR1 engagement is **ceritinib** at IC50 13 µM, weak and from an inhibition assay, the wrong direction.

So the fold responds to ligands. The approved cabinet is where it gets thin:

| Agent | Status | Objection |
|---|---|---|
| **Arimoclomol** (Miplyffa) | Approved 2024, oral, paediatric NPC, ATC N07XX17; USAN stem `-clomol` denotes heat shock protein inducers | No evidence BubR1 is an HSP90 or HSP70 client. Searched specifically; 22 papers at that intersection and none establish it. Class level rationale only |
| **Ursodeoxycholic acid** | Approved, decades of paediatric cholestasis use | Nonspecific; no BubR1 data |
| **Betaine** (Cystadane) | Approved, infancy onward | Nonspecific; no BubR1 data |

**There is a trap in this axis.** Thermal stabilization is not functional rescue. Two of the four compounds in that CETSA dataset are annotated as BUB1B inhibitors, and a ligand can rigidify the fold while locking out the very conformation KARD presentation and PP2A-B56 recruitment depend on. Measured by thermal shift alone that reads as success while having made things worse.

**So the success criterion is defined in advance: CETSA plus a functional readout.** A candidate passes only with a positive thermal shift on BubR1 **and** preserved or improved KARD phosphorylation or PP2A-B56 recruitment in the same system. A thermal shift on its own is a disqualifying result here, not a promising one. That requirement is cheap to implement, and it generalizes to any destabilizing missense in any pseudokinase (section 7).

### 3.4 PP2A activation, the approved way onto the Aurora B axis

**The axis with the best evidence has no approved drug.** ChEMBL lists nine Aurora B inhibitors with annotated mechanism; maximum clinical phase 3, none approved. Aurora B inhibition is therefore a **mechanistic keystone rather than a candidate**: it is what establishes the PP2A-B56 and Aurora B axis as correctable rather than merely described.

**But the lesion is lost PP2A-B56 recruitment, and the damaged side has approved agents.** The sphingolipid analogue **fingolimod (FTY720)** binds and targets the endogenous PP2A inhibitor I2PP2A/SET, producing PP2A reactivation ([DOI](https://doi.org/10.1002/emmm.201201283)), reproduced independently ([DOI](https://doi.org/10.1292/jvms.15-0193)) and reviewed as a strategy ([DOI](https://doi.org/10.3389/fonc.2014.00388)). Fingolimod is approved, oral, licensed from age 10.

Three objections, each more informative than the proposal. This is a **targeting failure, not an activity failure**: BubR1 cannot bring B56-PP2A *to the kinetochore*, and relieving SET inhibition raises global PP2A activity without obviously delivering B56-PP2A to the place the lesion emptied. **Immunosuppression** in a cancer predisposition syndrome may be disqualifying on its own, since fingolimod's approved mechanism is S1P receptor modulation and lymphocyte sequestration. And **the PP2A evidence comes from tumour suppression**, not kinetochore biology: right phosphatase, studied for a different reason.

**What this axis yields is a specification rather than a prescription**, which is a legitimate output of a repurposing analysis. The agent needed here restores B56-PP2A at the kinetochore specifically, without global phosphatase activation and without touching lymphocytes. It does not exist. Writing the specification down, with the three failure modes that would disqualify a candidate, is how one gets made.

**One caution applies to the whole axis.** Aneuploid cells are more sensitive to **genetic** perturbation of checkpoint components and less sensitive to short term **pharmacological** inhibition of the same components (Cohen-Sharir et al. 2021, *Nature*, [DOI](https://doi.org/10.1038/s41586-020-03114-6)). The rescue evidence sits on the genetic side, so any pharmacological attempt inherits that gap. The way across is to measure the pharmacological effect directly in patient cells rather than assume the genetic result transfers.

### 3.5 Readthrough of allele 1, an enabling experiment in search of an agent

**Rationale.** Allele 1 contributes nothing. Recovering even a fraction converts a null into a partial contributor, and in a dose dependent hypomorphic disease that is the difference between two half measures and one.

**What the sequence says, measured rather than assumed.** The premature stop is **UGA**, +4 context **A**. UGA is the most readthrough permissive of the three stop codons, reported as UGA then UAG then UAA ([DOI](https://doi.org/10.1016/j.bbrc.2022.10.081)), which favours the approach. The +4 nucleotide is a strong modifier but its hierarchy is **system dependent rather than universal**: in *HERG* the best responding context was TGA-T ([DOI](https://doi.org/10.3892/ijmm.2019.4360)), and stop codon context determines not only whether readthrough works but which drug and at what concentration ([DOI](https://doi.org/10.3390/biom11071006); review, [DOI](https://doi.org/10.1080/15476286.2015.1068497)).

**There is currently no usable approved agent, and that is the finding.**

| Agent | Status |
|---|---|
| **Ataluren** (Translarna) | **EU authorisation not renewed.** The EMA recommended non renewal, confirmed it, and re confirmed it on re examination; never FDA approved. ChEMBL still records max phase 4 with `withdrawn_flag: false`, which is stale. A caution about trusting a database field over the regulator |
| Aminoglycosides (gentamicin, amikacin) | Approved, but **not usable here**: nephrotoxic against existing nephrocalcinosis. The genome refines rather than compounds this, since the MT-RNR1 ototoxicity alleles `m.1555A>G` and `m.1494C>T` are both absent, so the exclusion rests on renal risk alone |

So this is an enabling experiment rather than a treatment proposal. It stays in because the codon is the favourable one and the resolving experiment is cheap and decisive: a **dual reporter readthrough assay carrying this exact BUB1B UGA-A context**, screened across a panel of readthrough agents. A positive makes the axis worth pursuing when a usable agent appears. A negative retires the hypothesis permanently, which is worth the cost of the assay on its own.

Two further issues belong in that same assay rather than in an assumption: readthrough is not selective for premature stops, so the burden on normal termination codons must be quantified, and an NMD competent transcript may be degraded before it can be read through, so the ceiling depends on transcript abundance rather than readthrough efficiency alone.

## 4. What this genome argues against

Repurposing analyses propose. In a caretaker disease the drugs not given matter as much as the drugs given.

### 4.1 The aneuploidy suppression playbook targets the opposite lesion

The least obvious de-selection here, and it follows directly from section 1.2.

In CIN cancer cells, kinetochore microtubule attachments are **inherently more stable** than in normal diploid cells, so errors escape correction and persist into anaphase, and artificially increasing attachment stability in normal cells is by itself sufficient to cause missegregation (Bakhoum, Genovese & Compton 2009, *Current Biology*, [DOI](https://doi.org/10.1016/j.cub.2009.09.055)). Aged human fibroblasts show the same direction, and missegregation there **fell** when MT depolymerizing kinesin-13 activity was enhanced, genetically and with a small molecule (Barroso-Vilares, Compton & Logarinho 2020, *EMBO Reports*, [DOI](https://doi.org/10.15252/embr.201949248)).

**This genome sits at the other end of the same axis.** Reduced BubR1 means reduced PP2A-B56 and excessive Aurora B, which elevates Ndc80 phosphorylation and makes attachments *less* stable: BubR1 acetylation is required "to counteract excessive Aurora B activity at the KT **for stable** chromosome-spindle attachments" (Park et al. 2013). Both states missegregate, by opposite routes. One fails to correct errors, the other fails to hold correct attachments.

So kinesin-13 and MCAK enhancement, and the small molecule strategy validated in ageing cells, is **contraindicated here**, because it pushes the dial the lesion has already pushed too far. An analysis matching on the phrase "reduces chromosomal instability" rather than on mechanism direction would have proposed exactly this, with good literature support, and been wrong.

### 4.2 Microtubule targeting agents, with the counterexample stated

Paclitaxel and vincristine act by prolonging spindle checkpoint activation and delaying mitotic exit, and the checkpoint they prolong is the one BubR1 operates. **Inhibition of BubR1 promotes resistance to microtubule inhibitors in vitro** (Swanton, Tomlinson & Downward 2006, *Cell Cycle*, a review, [DOI](https://doi.org/10.4161/cc.5.8.2682)), and independently, suppressing BubR1 function confers resistance to antimicrotubule drug induced apoptosis in a manner **reversed by BubR1 overexpression** ([DOI](https://doi.org/10.1007/s10549-005-6938-0)).

**The counterexample, stated because it is the only human evidence that exists.** The one published account of chemotherapy in an MVA patient used VAC. That patient started at the usual infant dose, developed thrombocytopenia after cycle 2 requiring platelet transfusions, had vincristine and cyclophosphamide reduced by 50%, **completed eight cycles, and the tumour was still shrinking six months after stopping.** A vincristine containing regimen worked in the only relevant human case on record.

So the defensible statement is narrow. Not "avoid vincristine", but that microtubule agents are **the one drug class whose efficacy here is predicted to be altered by the causal gene itself.** That makes response monitoring in any future regimen a genome informed decision rather than a routine one, which is actionable today and costs nothing.

### 4.3 Three further de-selections

**Aminoglycosides**, excluded on nephrotoxicity against existing nephrocalcinosis (section 3.5).

**Aneuploidy selective agents face a window problem specific to this case.** Five KIF18A inhibitors are recruiting, two listing chromosomal instability as an enrolment condition. Their selectivity rests on a contrast between unstable tumour and stable normal tissue. **There is no chromosomally stable tissue here**, so the window narrows from tumour versus normal body to tumour versus already unstable body. Nobody has measured that difference, and until someone does an otherwise well indicated class cannot be proposed.

**MPS1 inhibition is the wrong direction.** `EP 3645535 B1` covers an MPS1 inhibitor for "highly aneuploid tumors and tumors which overexpress mitotic checkpoint components", which means further checkpoint weakening in a checkpoint that is already hypomorphic.

### 4.4 One de-selection that may reverse

**4-phenylbutyrate is deprioritized within its own chaperone class** because it is also an HDAC inhibitor, and HDAC inhibition perturbs mitotic progression. A specific hazard rather than a generic side effect where the pathology is chromosome segregation.

**But BubR1 is deacetylated by HDAC11** (Watanabe et al. 2014, *Cell Reports*, [DOI](https://doi.org/10.1016/j.celrep.2014.03.022)), and section 2 establishes that acetylation stabilizes BubR1. If that holds, HDAC inhibition could raise BubR1, making 4-PBA's HDAC activity a feature rather than the hazard treated here. The counter is that 4-PBA is a weak class I and II inhibitor while HDAC11 is class IV, so the relevant selectivity is unknown. Left open rather than resolved, because the experiment is small and either outcome is valuable.

### 4.5 What is cleared

From the pharmacogenomic screen: reference, and therefore normal function, at all four *DPYD* risk alleles, all three *TPMT* alleles and *NUDT15*, *G6PD* (hemizygous in a male, so a single variant would be fully penetrant), and *SLCO1B1\*5*.

**Clean thiopurine and fluoropyrimidine status in a syndrome that predisposes to leukaemia is usable information today**, before any hypothesis here is tested. Caveats in section 8.

## 5. Where the field actually is

**Does each hypothesis meet this track's goal line?** The goal asks for existing, approved medications that are plausible candidates for further investigation. Answered before the argument:

| | Approved medicine? | Clears this patient's own picture? |
|---|---|---|
| **Nicotinamide**, abundance | **Yes**, 1982, oral and IV, paediatric | **Yes**, with a defined schedule on platelets and phosphate |
| **MEK**, tumour restricted | **Yes**, paediatric, mature safety | Yes on the drug, gated on one tissue test |
| **Chaperones** (arimoclomol, UDCA, betaine) | **Yes**, all dosed in children | Partly, no BubR1 specific evidence yet |
| **PP2A** (fingolimod) | **Yes**, from age 10 | Probably not, on immunosuppression |
| **Aurora B** | **No**, nine agents, max Phase 3 | Best mechanistic evidence, no drug |
| **Readthrough** | **No longer**, ataluren not renewed | No, aminoglycosides excluded on renal risk |

**Read the two columns against each other. They run in opposite directions.** The axis with the best biology has no medicine. The medicines with the best availability have the weakest link to this protein.

That is not a flaw in the analysis. It is the condition of drug repurposing in an ultra rare disease, and concealing it would hide the thing a clinician most needs to know. **Nicotinamide is the one candidate that breaks the pattern**, which is why it leads.

The practical consequence is a staged protocol rather than a single pick:

- **Now, at zero toxicity.** The pharmacogenomic clearances in section 4.5 and the microtubule monitoring note in section 4.2.
- **Cheap and decisive next.** The section 3.1 dose window experiment and the section 3.3 dual readout. Weeks rather than years, and either can retire its own hypothesis.
- **Gated.** MEK, on a tumour tissue result.
- **Blocked on drug availability rather than on biology.** Aurora B.

## 6. What would falsify this

- BubR1 protein is **not** reduced in patient cells, which makes the hypomorph model wrong and removes the rationale from 3.1 and 3.3.
- PP2A-B56 kinetochore recruitment is normal in patient cells, which breaks links 2 through 4.
- `p.Asn1002Lys` protein is **not** destabilized, which makes 3.3 a treatment for a lesion that is not there.
- Nicotinamide raises BubR1 in no dose window, or only where sirtuin inhibition dominates, which collapses 3.1 to nicotinic acid or to nothing.
- The two variants are **in cis**, which collapses the entire genetic model. Untestable with this data and the largest unverified assumption in this submission.
- Tumour tissue shows no *LZTR1* second hit, which closes 3.2.

## 7. What transfers beyond this case

*Scalability: the methods, not the drug.*

**The `INFO/AF` trap.** In a GATK VCF this field is per sample allele fraction, not population frequency. Filtering rare variants against it discards nearly every true heterozygous call. Applies to any rare disease analysis handed a HaplotypeCaller VCF.

**A symmetry discriminator for copy number versus mismapping**, from allele balance alone, no reference panel. A mosaic gain widens the heterozygous distribution symmetrically; mismapping is one sided low. Testing the two tails separately turns an ambiguous depth signal into a decision, with chrX as a built in positive control in any male sample.

**AVI requires a frequency pairing.** AlphaGenome's AVI is frequency naive by design. Here a variant 8 kb from a candidate gene scored 1.149 as a genuine splice acceptor variant and is carried by 52.7% of the population. Used alone it generates confident false leads next to the right gene.

**Mechanism direction, not disease label, decides transferability.** Section 4.1 is the general lesson: two diseases sharing the phrase "chromosomal instability" can sit at opposite ends of one axis, and a therapy validated for one can be contraindicated in the other.

**The pseudokinase as chaperone target pattern.** Pseudokinase pockets are structural rather than catalytic, so any destabilizing missense in any pseudokinase is a candidate for the CETSA plus function screen in section 3.3. A template rather than a one off.

**Stop codon context first readthrough triage.** Read the premature stop identity and +4 context off the reference, then choose the assay and agent panel, rather than assuming a universal hierarchy.

## 8. Limitations

1. **Phase is inferred, not observed**, and cannot be observed with this data: the variants are 10,911 bp apart, the reads are 2×149 bp, and there are no parental samples. The largest unverified assumption here.
2. **No tumour tissue**, so the *LZTR1* somatic hit is untestable, which is why 3.2 is gated.
3. **No approved Aurora B inhibitor at any phase**, which is why the axis with the best evidence is a keystone rather than a candidate.
4. **The BubR1 specific evidence is about abundance, not folding.** Sections 2 and 3.1 show BubR1 levels can be raised in vivo. They do not show a destabilized `p.Asn1002Lys` fold can be rescued. Adjacent hypotheses, not the same one.
5. **The in vivo abundance work is mouse, the readout was lifespan in a progeroid model, and the SIRT2 arm was genetic.** NMN was the only pharmacological arm.
6. **Aneuploidy was not directly observed** in this blood sample. The scan was negative with a roughly 15% mosaic floor and does not substitute for a karyotype in cultured cells.
7. **Pharmacogenomics is a single variant screen, not star allele diplotyping.** *CYP2D6* cannot be called reliably from short reads, and absence of a record means homozygous reference only if the position was covered.
8. ***TRIP13* is unassessed** because of subtelomeric coverage failure.
9. **One freedom to operate item is open.** Basilea's `WO 2012/098207` claims "increasing the level of BUBR1 in a subject that has a sample with a lower level of BUBR1" (`US 10222377`, adjusted expiry 2032). The granted US claims that can be read are measure, determine and kit claims, so the broad rescue claim most likely died in prosecution, but the full granted set should be checked before anyone asserts a clear path.
10. **This is a hypothesis document.** Nothing here has been tested in this patient, and nothing should be acted on without the experiments named alongside each candidate.

## 9. Open questions

*Came out of the work, not resolved by it.*

**Why is the MVA tumour spectrum narrow when the instability is global?** Missegregation is random at the cell level yet produces a recurrent, limited tumour set at the level of a person. Something selects, and nobody knows what.

**Was there ever a moment when these cells were euploid?** Constitutional aneuploidy cannot mean every cell at every instant. Whether a euploid progenitor compartment survives anywhere decides whether an intervention has something left to protect.

**Is BubR1 level a threshold or a gradient?** A threshold means partial recovery is worthless. A gradient means every increment counts.

**How has the immune system coexisted with this for years without clearing the aneuploid cells?** Somatic instability in cancer provokes recognition. Constitutional instability apparently does not.

**Does meiosis fail the same way mitosis does?** Carrier parents show premature chromatid separation, and the family history includes recurrent pregnancy loss.

**Can a molecule restore a phosphatase to one location without activating it everywhere?** Section 3.4 ends in a specification for an agent that does not exist.

**Does an aneuploidy selective agent have any window in a patient with no stable tissue?** The whole class assumes a contrast this genome does not have.

**Does HDAC11 inhibition raise BubR1?** A small experiment on a cheap approved drug, with both answers useful.

**Why is BubR1 essentially unclaimed as a therapeutic target?** A search across 306 patent records found no patent claiming a chaperone or proteostasis rescue of it. A validated target with an open mechanism and an empty landscape is an unusual combination.

## 10. Methods description

*Responses to the provided Track 2 methods template.*

**Team name.** `Mat34ec`, individual participant.

**Approach, from variant to candidate.** Variants were recovered by the Track 1 pipeline. Mechanism was characterized by combining AVI feature attribution, to establish that the impact is protein coding rather than regulatory, two independent model families to exclude splicing, and primary literature for the pathway. Candidates were then derived **from the characterized mechanism**: each names the specific link in section 1.2 it acts on, and was then tested against this patient's own clinical and genomic picture. Agents were identified through ChEMBL for target records, mechanism of action, clinical phase and CETSA bioactivity; ClinicalTrials.gov for clinical stage and paediatric precedent; PubMed; and Lens.org for the patent landscape and freedom to operate.

**Generative AI disclosure (required).** Anthropic Claude Pro and Claude Science, Claude Opus 5, non commercial terms, no training on customer content. Google Gemini 3.1 Pro, non commercial terms, no training on customer content. Perplexity free tier as a search engine, non commercial, no training on customer content.

**The division of labour between the human and the models, since the track asks.** The models did retrieval at a scale a single person cannot, wrote the analysis code, and were used deliberately as adversaries against my own hypotheses. I set the direction, chose the decision criteria, supplied the biological framing, and made every call about what survived. **Every quantitative claim in this report derives from a script in the repository or from a named primary source; no figure is reported from model recall.** Where a model supplied citation could not be verified against the primary source it was removed, and three such cases are recorded below.

**Automated or manual candidate identification?** Both, in sequence. Database queries were automated. Mechanism characterization, hypothesis generation, the ordering in section 3 and every objection raised against my own candidates are manual expert reasoning.

**Manual literature review and curation.** Extensive, and it changed conclusions three times. A proposed LZTR1 to KARD to PP2A to Aurora B connection was **discarded**: it does not exist, and it would have been the most damaging error in the submission. A single locus pharmacogenomic result suggesting a vincristine neuropathy risk genotype **reversed** when the full panel was run, because of five reported loci only the weakest and most contested is carried, with reference genotypes at the two with the largest reported effects. A cited PDB structure proved to be **an unrelated bacterial enzyme**, and the real BUB1B entries were substituted (`2WVI`, `3SI5`, `4GGD`, `5JJA`, `5K6S`, `5KHU`, `5LCW`, `5SWF`, `6TLJ`, UniProt O60566). These are recorded rather than quietly removed, because a method that never reports its own failures gives a reviewer nothing to calibrate against.

**Public data only?** Yes, exclusively. No proprietary sources.

**Public sources used.** ChEMBL, ClinicalTrials.gov, PubMed, Lens.org, Ensembl REST and VEP, gnomAD, ClinVar, GENCODE, AlphaGenome and AlphaGenome Atlas, SpliceAI v1.3.1, EMA public assessment reports.

**How was the mechanism characterized?** See section 1. Loss of function, partial, hypomorphic. Pathway: BubR1, then KARD phosphorylation, then PP2A-B56 recruitment, then opposition of Aurora B at the kinetochore. Downstream consequence: premature chromatid separation and mosaic aneuploidy. Disease class: caretaker rather than gatekeeper, which determines that the fitting therapeutic category is chemoprevention.

**Effort.** Roughly 20 days part time, of which the majority went on validation and self refutation rather than on hypothesis generation.

## 11. Sources

Xu P, Raetz EA, Kitagawa M, Virshup DM, Lee SH (2013). BUBR1 recruits PP2A via the B56 family of targeting subunits to promote chromosome congression. *Biology Open* 2(5):479-86. PMID 23789096. [DOI](https://doi.org/10.1242/bio.20134051)

Choi E, Choe H, Min J, Choi JY, Kim J, Lee H (2009). BubR1 acetylation at prometaphase is required for modulating APC/C activity and timing of mitosis. *EMBO J* 28(14):2077-89. PMID 19407811. [DOI](https://doi.org/10.1038/emboj.2009.123)

Park I, et al. (2013). Loss of BubR1 acetylation causes defects in spindle assembly checkpoint signaling and promotes tumor formation. *J Cell Biol* 202(2):295-309. PMID 23878276. [DOI](https://doi.org/10.1083/jcb.201210099)

North BJ, et al. (2014). SIRT2 induces the checkpoint kinase BubR1 to increase lifespan. *EMBO J* 33(13):1438-53. PMID 24825348. [DOI](https://doi.org/10.15252/embj.201386907)

Chen AC, et al. (2015). A phase 3 randomized trial of nicotinamide for skin-cancer chemoprevention. *NEJM* 373(17):1618-26. PMID 26488693. [DOI](https://doi.org/10.1056/NEJMoa1506197)

Ketteler M, et al. (2021). Efficacy and safety of a novel nicotinamide modified-release formulation in refractory hyperphosphatemia in haemodialysis. *Kidney Int Rep* 6(3):594-604. PMID 33732974. [DOI](https://doi.org/10.1016/j.ekir.2020.12.012)

Shahbazian H, et al. (2011). Oral nicotinamide reduces serum phosphorus, increases HDL, and induces thrombocytopenia in hemodialysis patients. *Nefrologia* 31(1):58-65. PMID 21270914. [DOI](https://doi.org/10.3265/Nefrologia.pre2010.Nov.10734)

Ix JH, et al. (2019). Effects of nicotinamide and lanthanum carbonate on serum phosphate and FGF-23 in CKD: the COMBINE trial. *JASN* 30(6):1096-1108. PMID 31085679. [DOI](https://doi.org/10.1681/ASN.2018101058)

Bakhoum SF, Genovese G, Compton DA (2009). Deviant kinetochore microtubule dynamics underlie chromosomal instability. *Current Biology* 19(22):1937-42. PMID 19879145. [DOI](https://doi.org/10.1016/j.cub.2009.09.055)

Barroso-Vilares M, Macedo JC, Reis M, Warren JD, Compton D, Logarinho E (2020). Small-molecule inhibition of aging-associated chromosomal instability delays cellular senescence. *EMBO Reports* 21(5):e49248. PMID 32134180. [DOI](https://doi.org/10.15252/embr.201949248)

Cohen-Sharir Y, et al. (2021). Aneuploidy renders cancer cells vulnerable to mitotic checkpoint inhibition. *Nature*. [DOI](https://doi.org/10.1038/s41586-020-03114-6)

Yost S, et al. (2017). Biallelic TRIP13 mutations predispose to Wilms tumor and chromosome missegregation. *Nature Genetics*. [DOI](https://doi.org/10.1038/ng.3883)

Matsuura S, et al. (2006). Monoallelic BUB1B mutations and defective mitotic-spindle checkpoint in seven families with premature chromatid separation syndrome. *Am J Med Genet A*. [DOI](https://doi.org/10.1002/ajmg.a.31069)

García-Castillo H, et al. (2008). Clinical and genetic heterogeneity in patients with mosaic variegated aneuploidy. *Am J Med Genet A*. [DOI](https://doi.org/10.1002/ajmg.a.32315)

Nozaki T, et al. (2015). A novel in vitro model of sarcopenia using BubR1 hypomorphic C2C12 myoblasts. *Cytotechnology* 68(5):1705-15. [DOI](https://doi.org/10.1007/s10616-015-9920-7)

Watanabe Y, Khodosevich K, Monyer H (2014). Dendrite development regulated by the schizophrenia-associated gene FEZ1 involves the ubiquitin proteasome system. *Cell Reports* 7(2):552-64. PMID 24726361. [DOI](https://doi.org/10.1016/j.celrep.2014.03.022)

Saddoughi SA, et al. (2013). Sphingosine analogue drug FTY720 targets I2PP2A/SET and mediates lung tumour suppression via activation of PP2A-RIPK1-dependent necroptosis. *EMBO Mol Med* 5(1):105-21. PMID 23180565. [DOI](https://doi.org/10.1002/emmm.201201283)

Enjoji S, et al. (2015). The therapeutic effects of SET/I2PP2A inhibitors on canine melanoma. *J Vet Med Sci* 77(11):1451-6. [DOI](https://doi.org/10.1292/jvms.15-0193)

Oaks J, Ogretmen B (2015). Regulation of PP2A by sphingolipid metabolism and signaling. *Frontiers in Oncology* 4:388. [DOI](https://doi.org/10.3389/fonc.2014.00388)

Swanton C, Tomlinson I, Downward J (2006). Chromosomal instability, colorectal cancer and taxane resistance. *Cell Cycle* 5(8):818-23. **Review.** PMID 16628000. [DOI](https://doi.org/10.4161/cc.5.8.2682)

Inaba S, et al. (2005). Synuclein gamma inhibits the mitotic checkpoint function and promotes chromosomal instability of breast cancer cells. *Breast Cancer Res Treat* 94(1):25-35. PMID 16142440. [DOI](https://doi.org/10.1007/s10549-005-6938-0)

Ngiwsara L, et al. (2022). Effects of gentamicin inducing readthrough premature stop codons. *Biochem Biophys Res Commun*. [DOI](https://doi.org/10.1016/j.bbrc.2022.10.081)

Yu H, et al. (2019). Stop codons and the +4 nucleotide may influence the efficiency of G418 in rescuing nonsense mutations of the HERG gene. *Int J Mol Med*. [DOI](https://doi.org/10.3892/ijmm.2019.4360)

Schilff M, Sargsyan Y, Hofhuis J, Thoms S (2021). Stop codon context-specific induction of translational readthrough. *Biomolecules* 11(7):1006. [DOI](https://doi.org/10.3390/biom11071006)

Dabrowski M, Bukowy-Bieryllo Z, Zietkiewicz E (2015). Translational readthrough potential of natural termination codons in eucaryotes. *RNA Biology*. [DOI](https://doi.org/10.1080/15476286.2015.1068497)

Tabe S, et al. (2026). Spindle assembly checkpoint competency determines sensitivity to KIF18A inhibition in small cell lung cancer. *Cancer Research Communications*. [DOI](https://doi.org/10.1158/2767-9764.CRC-26-0560)

**Database and regulatory records.** ChEMBL targets `CHEMBL4295998` (BUB1B), `CHEMBL2185` (Aurora kinase B), `CHEMBL4523403` (KIF18A); CETSA assay `CHEMBL6125447`, document `CHEMBL6121040`; compounds `CHEMBL1140` (nicotinamide), `CHEMBL2107726` (arimoclomol), `CHEMBL1746` and `CHEMBL2105745` (phenylbutyrate forms), `CHEMBL256997` (ataluren). EMA public assessment reports and news releases on the non renewal of Translarna's marketing authorisation. Lens.org patent records `WO 2012/098207`, `US 10222377`, `EP 1831398 B1`, `US 2008/0102477`, `US 11186633`, `EP 3645535 B1`. ClinicalTrials.gov NCT06084416, NCT05902988, NCT06799065, NCT06772415, NCT07226427 (KIF18A); NCT03155620, NCT04924608, NCT05849662, NCT07817186 (MEK and RASopathy).

---

Thanks for the opportunity to work on a meaningful challenge. I am really grateful for it, and if I can help further in any way, let me know.

**Mateo Cueva** · [mateo@mateocueva.com](mailto:mateo@mateocueva.com)

*Released under CC BY 4.0 in accordance with the hackathon rules.*
