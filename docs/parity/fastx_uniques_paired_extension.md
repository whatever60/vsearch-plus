# Native paired-end support for `fastx_uniques`

This document describes the native paired-end path for `vsearch --fastx_uniques` implemented in this workspace.

## 1) What is the extension

### Scope

Stock `fastx_uniques` dereplicates one sequence stream. The native paired path adds a paired-end mode when paired input is provided (second positional input for `R2`, or `--interleaved`), so each unique unit is a pair of anchors:

- left anchor: from forward read (R1)
- right anchor: from reverse read (R2), kept in native read orientation

A pair is considered identical only when both anchors are exactly identical.

### CLI trigger

Paired extension is selected by running `--fastx_uniques` with paired input:

- split paired input: provide `R2` as the second positional input
- interleaved paired input: set `--interleaved` and provide one input stream

Example:

```bash
./scripts/vsearch-plus \
  --fastx_uniques reads_R1.fastq.gz reads_R2.fastq.gz \
  --fastaout tav_left.fasta \
  --fastaout_rev tav_right.fasta \
  --tabbedout tav_catalog.tsv
```

Without paired input (no second positional input and no `--interleaved`), command dispatch remains on stock derep behavior in `cpp/src/derep.cc`.

### Input expectations

- Split mode: R1 and R2 are read in lockstep from the two input streams.
- Interleaved mode: input is consumed as `(R1,R2,R1,R2,...)`; odd record count is an error.
- If one file has extra records, execution fails.
- Anchor length is `min(len(R1), len(R2))`, additionally bounded by `--fastq_trunclen` when provided.

Low-level relationship in current code:

- `fastx_uniques` input call stack:
  ```text
  main command switch in vsearch.cc
    -> derep_paired(...)
       -> paired FASTX read loop
       -> anchor truncation + pair-key derep
  ```
- `cluster_unoise` / `uchime3_denovo` input call stack:
  ```text
  cluster_unoise_paired(...)
    -> native paired FASTX load in cluster_paired.cc
    -> per-end qmask/hardmask during record ingestion

  uchime3_denovo_paired(...)
    -> native paired FASTX load in chimera_paired.cc
    -> per-end qmask/hardmask during record ingestion
  ```

### Deduplication rule

For each synchronized pair:

1. extract left anchor from R1
2. extract right anchor from R2 in native orientation
3. define key as `(left_anchor, right_anchor)`
4. aggregate abundance into the same key only when both sides match exactly

### Orientation policy

R2 is kept in native read orientation for external paired inputs/outputs. This removes implicit orientation flips across command boundaries. If a downstream algorithm needs orientation normalization, it must do that internally during computation.

### Outputs

#### Catalog TSV (`--tabbedout`)

Custom paired catalog with columns:

- `tav_id`
- `abundance`
- `header`
- `left_anchor`
- `right_anchor`

Notes:

- `tav_id` currently stores the representative header (stock-like label behavior), not an auto-generated `TAVxxxxxx` token.
- `right_anchor` is written in native R2 orientation.

#### Left FASTA (`--fastaout`)

Contains unique left anchors, one per paired unique, with abundance annotations handled by stock FASTA printer options.

#### Right FASTA (`--fastaout_rev`)

Contains unique right anchors in native R2 orientation.

Header labels for left and right FASTA are intentionally the same representative ID for each paired unit.

## 2) Stock behavior parity for non-core behavior

This section documents what is intentionally kept aligned with stock `fastx_uniques` behavior.

### Tie-breaking and deterministic order

Sort order is:

1. abundance descending
2. representative header lexicographic ascending
3. first-seen order in input

This mirrors stock derep intent (abundance first, then header, then stable first-seen ordering).

### Representative header convention

Representative header is the header of the first observed read pair for a unique paired key. This mirrors stock representative-header semantics.

### Header truncation behavior (`--notrunclabels`)

Native paired mode now follows stock parser behavior:

- default: truncate at first whitespace
- with `--notrunclabels`: keep full header line (excluding newline)

### Relabel and header-format options

FASTA header rendering goes through stock print paths, so stock options continue to apply in paired FASTA outputs, including:

- `--relabel`, `--relabel_self`, `--relabel_md5`, `--relabel_sha1`, `--relabel_keep`
- `--label_suffix`
- `--sizeout`, `--xee`, `--xlength`, `--xsize`

### Abundance semantics

Abundance aggregation uses the same size-input concept as stock (`--sizein`-aware read abundance from FASTX parser), then emits abundance through stock output formatting where applicable.

## 3) Differences from stock to be aware of

- The paired catalog (`--tabbedout`) is a custom paired schema, not stock derep tabbed format.
- Paired mode writes two coordinated FASTA streams (`--fastaout` and `--fastaout_rev`) for left/right anchors.
- In paired mode, right anchor sequences keep native R2 orientation.

## 4) Ownership status

- Stock single-end owner: `cpp/src/derep.cc`
- Native paired owner: `cpp/src/derep_paired.cc`
- Legacy `tav_fastx_uniques(...)` is retired from the CLI path and is no longer part of the build.

## 5) Practical guidance

- Use paired mode whenever downstream steps are paired-aware.
- Keep `--fastq_trunclen` consistent between runs if reproducibility of anchors is important.
- Use `--notrunclabels` if full Illumina header suffix fields must be preserved in representative labels.
