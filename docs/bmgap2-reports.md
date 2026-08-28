# BMGAP2 reports

[BMGAP2](https://github.com/CDCgov/BMGAP2) characterises *Neisseria meningitidis* and
*Haemophilus influenzae*: AMR variants, vaccine antigens and a species confirmation.
It produces no summary of its own, so Sanibel derives one. `nm_sum_report.txt` and
`hi_sum_report.txt` are that derived layer, built by
[`parse_bmgap2()`](../bin/summary_report.py) from the outputs of the three
[`bmgap2_*` modules](../modules).

Because the reports are derived rather than produced by BMGAP2, the distinction between
what BMGAP2 said and what Sanibel concluded matters when reading them. Every column below
notes which it is.

## What runs, and for which samples

PMGA and the BMGAP2 chain execute for any sample whose genus is *Neisseria* or
*Haemophilus* ([`sanibel.nf`](../sanibel.nf)), and the three BMGAP2 wrapper scripts skip
anything whose MLST scheme is not `neisseria` or `hinfluenzae`
([`bmgap2_helpers.py`](../bin/bmgap2_helpers.py)).

**Neither of those gates decides what reaches the reports.** `mlst` maps the whole
*Neisseria* genus to the `neisseria` scheme in
[`mlst_schemes.tsv`](../assets/mlst_schemes.tsv), so a *N. gonorrhoeae* or *N. lactamica*
isolate passes both. Report routing instead uses the skani species call, with BMScan as a
fallback, via `meningitis_organism()` in [`summary_report.py`](../bin/summary_report.py).
Only a confirmed *N. meningitidis* or *H. influenzae* gets a row.

A non-meningitidis *Neisseria* therefore appears in `sum_report.txt` with `No data` in the
serotype column and is absent from `nm_sum_report.txt`. Its BMGAP2 output still exists
under the sample's output directory; it is simply not reported, because a meningococcal
serogroup and a meningococcal AMR interpretation do not apply to it.

`bmgap2_locusextractor` and its siblings run
[**two at a time**](../nextflow.config), matching the limit in BMGAP2's own SGE runner.

## Reading the AMR columns

Each gene contributes an allele, its curated mutations and a phenotype, read through
BMGAP2's own per-gene `status` field. The gene sets differ by organism: `runAST` screens
`ponA` for *N. meningitidis* but not *H. influenzae*, and `folA` and `acrR` the other way
round. A gene missing from BMGAP2's output was never screened for that species, which is
why those cells read `No data` rather than `Not detected`.

Two columns describe resistance and they can disagree:

| Column | Source |
| ------ | ------ |
| `predicted_resistance` | BMGAP2's own `summary.predicted_resistance`, passed through verbatim |
| `resistance_markers` | The drug names BMGAP2 attaches to each curated mutation it found, deduplicated |

They disagree because BMGAP2's summary does not always reflect its own annotations. On one
*H. influenzae* isolate `predicted_resistance` read `None` while gyrA carried S84L annotated
`Ciprofloxacin; Levofloxacin`, which `resistance_markers` reported. Both columns are BMGAP2
data. Sanibel aggregates the second rather than adjudicating between them, so a
disagreement is visible instead of hidden.

`resistance_markers` inherits BMGAP2's vocabulary, which mixes drug names and class names,
so a cell can read `Ampicillin;Cephalosporins;Penicillin;Penicillins`.

## Reading the capsule and antigen columns

`nm_genogroup` and `hi_capsule_genotype` come from PMGA and are **genotypic capsule
predictions, not phenotypic serogrouping**. `serotype_notes` carries PMGA's assessment of
the capsule locus, which is where a non-groupable isolate declares itself, for example
`E backbone: missing cseA;missing cseB;...` alongside a genogroup of `NG`.

`vaccine_antigens_present` lists which of fHbp, NHBA, NadA and PorA were found, annotating
a novel allele as `(new)` and a disrupted one as `(disrupted)`. **It is not a vaccine
coverage prediction.** Coverage depends on peptide identity and cross-reactivity, which
this pipeline does not assess. An antigen counts as present if any of its lookup columns
resolves, because the Novartis, Oxford and Pfizer schemes fail independently and a novel
allele often resolves in only one of them.

## The `bmgap2_status` column

The last column of both reports. It reads `ok` when all three BMGAP2 components produced
usable output, otherwise a `;`-joined list of the components that did not, in the form
`amr:<state>`, `le:<state>`, `bmscan:<state>`.

| State | Meaning |
| ----- | ------- |
| `ok` | the component produced usable output |
| `no_output` | the component produced no output file |
| `no_calls` | LocusExtractor wrote a results file with every field empty |
| `ran_no_csv` | the `LE_*` directory exists but holds no `molecular_data` CSV |
| `failed` | the file exists but could not be parsed |
| `runAST_exit_<N>` | `runAST.py` exited non-zero and produced nothing |

**Anything other than `ok` means that row's BMGAP2 columns are incomplete and should not
be used.** The column exists because none of these states raises an error on its own.

## Known limitations

**LocusExtractor can fail without failing.** It catches its own errors, writes a results
file with every field blank and exits 0. A run has been observed where `makeblastdb` failed
for all three staged inputs and every layer still reported success. `bmgap2_status` reads
`no_calls` in that case, which is the only signal, and the antigen and MLST columns will
read `No data`. The root cause is not established; LocusExtractor discards `makeblastdb`'s
error output, so it cannot be recovered from the logs. It has been observed on roughly one
sample in ten, affecting a different sample each run.

## Checking a run

`params.output/pipeline_info/trace.txt` records every task. A healthy
`bmgap2_locusextractor` task runs for minutes at 85 to 97 percent CPU; the observed failure
ran 23 seconds at 29 percent.

```bash
# any BMGAP2 row whose columns are incomplete
cut -f1,32 <output>/nm_sum_report.txt
cut -f1,27 <output>/hi_sum_report.txt

# tasks that failed and were ignored
awk -F'\t' 'NR==1 || $5!="COMPLETED"' <output>/pipeline_info/trace.txt
```
