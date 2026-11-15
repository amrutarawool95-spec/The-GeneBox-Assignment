# Bioinformatics / Research Assignment: Submission

**Candidate:** Amruta Rawool

## What is in this folder

```
submission/
├── README.md              this file
├── part1.md               Part 1: 5 technical short answers (1.1, 1.2, 1.5, 1.6, 1.7)
├── part3.md               Part 3: all 10 mathematical problems with derivations
├── FINDINGS.md            Part 2 (Problem C): data defects, filter evaluation, normal_type
├── src/
│   └── problemC.py        Part 2 pipeline (Python standard library only)
├── out/
│   ├── problemC_arguments.tsv   C1: workflow argument specification
│   ├── problemC_calls.tsv       C2: all 48 sites with one filter each
│   ├── problemC_summary.json    C2: counts, TMB, median PASS VAF
│   └── problemC_run.md          C2: command, metrics, summary, TMB verdict
├── problemC_somatic/      input files exactly as provided (never edited)
└── tests/
    └── test_pipeline.py   9 checks on the pipeline (optional)
```

## How to run

From inside `submission/`:

```
python3 src/problemC.py problemC_somatic out
```

- Needs Python 3 and nothing else (standard library only, no install step).
- Reads the 5 TSV files, `target_intervals.bed` and `sample_meta.json` from `problemC_somatic/`.
- Writes `problemC_calls.tsv` and `problemC_summary.json` to `out/`.
- Deterministic: no randomness and no network. Repeated runs give byte-identical files.
- Data problems are printed as warnings on screen; ambiguous input (missing columns, conflicting duplicate rows) stops the run with an error instead of guessing.

To run the checks:

```
python3 -m unittest tests.test_pipeline -v
```

## Why I chose Problem C

- **It is the most complete test of the whole workflow.** Problem C covers read-level evidence, reference integration (germline, panel of normals, repeats) and final filtering, which is closest to real clinical bioinformatics work.
- **It fits my background.** My focus areas are computational biology and variant-level analysis, so I could check the output against biological reasoning and not only against the rules.
- **The data had several traps worth finding,** for example multiallelic loci where the pileup allele differs from the resource allele, BED coordinate conventions, and two alt alleles at one position. A careful implementation shows its value there.
- The assignment prefers one thorough problem over partial attempts, so I did not submit Problems A or B.

## Time allocation

| Part | Hours |
|---|---|
| Part 1: 5 short answers and source checking | 5 |
| Part 2: Problem C pipeline, outputs, FINDINGS | 1 |
| Part 3: 10 mathematical problems | 10 |
| Checking, testing and write-up | 1 |
| **Total** | 17 |

## Omissions and limits (stated honestly)

- **Problems A and B were not attempted,** by choice (one problem is required).
- **Part 1:** I answered 5 of 10 questions: 1.1, 1.2, 1.5, 1.6 and 1.7. The other five were not written up.
- **C1 is a specification, not an executed workflow.** No raw sequencing reads were provided, so the GATK, bwa-mem2, fastp and samtools commands in `problemC_arguments.tsv` were not run. The flags come from tool documentation and my own knowledge and should be checked against the tool versions in use.
- **TMB is reported but judged unsuitable for clinical use** (see `problemC_run.md`). I made no attempt to calibrate it against whole-exome sequencing, since no reference data was given.
- **The filter thresholds are the ones given in the assignment.** I did not tune them. Where I think a threshold misbehaves, it is described in `FINDINGS.md` and not changed in the code.
- **CYP2D6 values (Part 1, Q1.5):** I give the 2020 CPIC/DPWG values and the March 2023 update for \*41. The assignment's own allele table uses the 2020 value (\*41 = 0.5).

## Tools and verification

External tools and LLMs were used, as the guidelines allow. Every quantitative result was checked before submission:
- Part 3 numbers were recomputed with a short Python script (scipy).
- The Problem C pipeline was run repeatedly to confirm identical output, tested on shuffled and reformatted copies of the inputs, and its results were compared against an independent implementation.
- Part 1 citations are published papers, tool manuals or the VCF specification. I did not cite source-code line numbers I could not verify.
