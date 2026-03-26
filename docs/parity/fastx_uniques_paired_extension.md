# Paired-end extension for fastx_uniques

This document describes the paired-end extension path for `vsearch --fastx_uniques` implemented in this workspace.

## 1) What is the extension

### Scope

Stock `fastx_uniques` dereplicates one sequence stream. The extension adds a paired-end mode when `--reverse` is provided, so each unique unit is a pair of anchors:

- left anchor: from forward read (R1)
- right anchor: from reverse read (R2), reverse-complemented before pairing

A pair is considered identical only when both anchors are exactly identical.

### CLI trigger

Paired extension is selected by running `--fastx_uniques` with `--reverse`.

Example:

```bash
./bin/vsearch \
  --fastx_uniques reads_R1.fastq.gz \
  --reverse reads_R2.fastq.gz \
  --fastaout tav_left.fasta \
  --fastaout_rev tav_right.fasta \
  --tabbedout tav_catalog.tsv
```

Without `--reverse`, command dispatch remains on stock derep behavior.

### Input expectations

- R1 and R2 are read in lockstep.
- If one file has extra records, execution fails.
- Anchor length is `min(len(R1), len(R2))`, additionally bounded by `--fastq_trunclen` when provided.

### Deduplication rule

For each synchronized pair:

1. extract left anchor from R1
2. reverse-complement R2, then extract right anchor
3. define key as `(left_anchor, right_anchor)`
4. aggregate abundance into the same key only when both sides match exactly

### Why R2 is reverse-complemented

R2 is reverse-complemented so both anchors are normalized to the same biological orientation before downstream distance/chimera/search logic. This keeps pairwise scoring simpler and consistent across pipeline stages.

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
- `right_anchor` is already reverse-complemented sequence (normalized orientation).

#### Left FASTA (`--fastaout`)

Contains unique left anchors, one per paired unique, with abundance annotations handled by stock FASTA printer options.

#### Right FASTA (`--fastaout_rev`)

Contains unique right anchors, where sequence is reverse-complemented R2 anchor.

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

Paired extension now follows stock parser behavior:

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
- In paired mode, right anchor sequences are reverse-complemented by design.

## 4) Practical guidance

- Use paired mode whenever downstream steps are paired-aware.
- Keep `--fastq_trunclen` consistent between runs if reproducibility of anchors is important.
- Use `--notrunclabels` if full Illumina header suffix fields must be preserved in representative labels.
