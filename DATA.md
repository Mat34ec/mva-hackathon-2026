# Data access, handling and deletion

**No patient data is included in this repository, and none may be added to it.**

## What the data is

The genome and clinical record used in this work belong to a real, living child with Mosaic Variegated
Aneuploidy, and to their family. They were released by Sage Bionetworks under a controlled-access
protocol approved by WCG IRB (**protocol #20252010**), with the family's explicit permission, in the
hope that the research community could help.

That provenance sets the terms. The rules below are not repository conventions; they are the conditions
under which access was granted.

## Obtaining the data

Request access on the dataset page and agree to the hackathon data-use terms:

```
https://huggingface.co/datasets/SageBio/mva-hackathon-2026-data
```

```bash
hf download SageBio/mva-hackathon-2026-data \
  --repo-type dataset --local-dir ./data \
  --exclude "README.md" ".gitattributes"
```

The full dataset is ~85 GB across 11 files (8 FASTQ lanes, a VCF with index, a clinical phenotype
document). Allow 100–150 GB of working space for caching, indexing and intermediates.

## Obligations that come with it

Agreed to as a condition of access:

- **Delete all data after the hackathon**, from every personal compute environment.
- **Notify Sage Bionetworks by email once deleted**, at `MVAHackathon2026@synapse.org`.
- **Do not attempt to recontact** the data subject, their family, or any MVA Society point of contact.
- **Do not redistribute** or grant access to anyone; every team member registers individually.
- **Report any unauthorized disclosure** to the Sage Bionetworks Privacy and Compliance Office.
- Nothing published may contain information that could re-identify the child or the family.

## How this repository is kept clean

`.gitignore` denies by **file extension and directory**, not by filename, so that a differently-named
copy of the same data cannot slip through. It excludes all FASTQ, BAM, CRAM, VCF, FASTA and index
formats; any `data/`, `Dataset/`, `raw/` or `RESTRICTED*/` directory; all `.docx` files; and the
per-variant genotype tables derived from the genome, which are individual-level data even though they
are plain text.

Results committed here are **derived summaries only** — distributions, percentiles, per-chromosome
statistics, counts — never per-variant genotypes, and never anything from the clinical document beyond
mechanism-level statements.

The two causal variants and the secondary findings appear in the Track 1 submission file and the
reports. That is unavoidable: identifying them is the task, and the hackathon requires submissions to be
released under CC BY 4.0.

## API credentials

`ALPHAGENOME_API_KEY` and any other credential live in the shell environment and never in a script, a
notebook, a config file or a commit:

```bash
export ALPHAGENOME_API_KEY=...                    # bash/zsh
$env:ALPHAGENOME_API_KEY = '...'                  # PowerShell, current session
```
