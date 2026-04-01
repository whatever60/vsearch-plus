# Paired `fastq_filter` Parity (Stock vs Extension Overlay)

This note documents parity and intentional divergence for `vsearch --fastq_filter` after adding the paired extension overlay.

## Scope and intent

- Keep stock paired behavior untouched when users explicitly use `--reverse`.
- Add extension behavior on top of stock behavior when users provide paired extension input (`R1 R2` positional input or `--interleaved`).
- Keep all non-EE filter logic stock-like.
- Change `--fastq_maxee` and `--fastq_maxee_rate` to true pair-level criteria in extension mode.

## Mode split contract

### Stock mode (unchanged)

Trigger:
- `--fastq_filter <R1>` with explicit `--reverse <R2>`

Output conventions:
- Uses stock paired output options:
  - `--fastqout_rev`
  - `--fastqout_discarded_rev`
  - `--fastaout_rev`
  - `--fastaout_discarded_rev`

Execution path:
- `cmd_fastq_filter()` routes to stock `fastq_filter()`.
- Stock `filter(true, ...)` path is used.

### Extension mode (new)

Trigger:
- split paired extension input: `--fastq_filter <R1> <R2>`
- interleaved paired extension input: `--fastq_filter <interleaved> --interleaved`

Output conventions:
- Uses paired split output options:
  - `--fastqout` + `--fastqout2`
  - `--fastqout_discarded` + `--fastqout_discarded2`
  - `--fastaout` + `--fastaout2`
  - `--fastaout_discarded` + `--fastaout_discarded2`
- `*_rev` options are rejected in extension mode.

Execution path:
- `cmd_fastq_filter()` routes to `fastq_filter_paired_ext()`.
- Extension implementation is `filter_paired_ext_fastq(...)`.

## Filtering semantics

### Shared behavior (stock + extension)

- Per-read trimming:
  - `--fastq_stripleft`
  - `--fastq_stripright`
  - `--fastq_trunclen`
  - `--fastq_trunclen_keep`
  - quality/EE truncation (`--fastq_truncqual`, `--fastq_truncee`, `--fastq_truncee_rate`)
- Per-read quality/length/N/abundance guards:
  - `--fastq_minqual`, `--fastq_minlen`, `--fastq_maxlen`, `--fastq_maxns`
  - `--minsize`, `--maxsize`

## Intentional extension divergence: EE thresholds are pair-level

In extension mode only:

- `pair_ee = ee_r1 + ee_r2`
- `pair_len = len_r1 + len_r2`
- Pair is discarded if either:
  - `pair_ee > --fastq_maxee`
  - `pair_len > 0` and `(pair_ee / pair_len) > --fastq_maxee_rate`

This replaces stock paired behavior for these two criteria (which effectively applies per-read EE checks and then combines with pair keep/discard logic).

## Guardrails and validation rules

- Positional R2 for `fastq_filter` is extension mode only.
- `--reverse` + positional R2 together are rejected.
- In stock mode, `*2` options are rejected.
- In extension mode, `*_rev` options are rejected.
- In extension mode, paired split outputs must be provided as complete pairs (for example `--fastqout` with `--fastqout2`).

## Implementation map

- CLI routing:
  ```text
  cmd_fastq_filter(...)
    -> stock mode: fastq_filter(...)
    -> extension mode: fastq_filter_paired_ext(...)
  ```
- Stock EE decision call stack:
  ```text
  fastq_filter(...)
    -> filter(true, ...)
    -> analyse(..., apply_ee_filters=true)
    -> ee_thresholds_pass(...)
  ```
- Extension EE decision call stack:
  ```text
  fastq_filter_paired_ext(...)
    -> filter_paired_ext_fastq(...)
    -> analyse(..., apply_ee_filters=false) [R1 + R2]
    -> pair aggregation
    -> ee_thresholds_pass(pair_ee, pair_len)
  ```
- Shared low-level kernel location:
  `cpp/src/filter.cc` -> `ee_thresholds_pass(...)`

## Orientation note

- `fastq_filter` does not reverse-complement R2.
- It follows the same external orientation policy as the unified paired pipeline: keep R2 in native orientation.
