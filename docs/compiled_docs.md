# VSEARCH Plus Documentation

_Generated on 2026-03-27 08:26 UTC from Markdown files under `docs/parity/`._

## Source: `docs/parity/fastq_filter_paired_parity.md`

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
  `src/filter.cc` -> `ee_thresholds_pass(...)`

## Orientation note

- `fastq_filter` does not reverse-complement R2.
- It follows the same external orientation policy as the unified paired pipeline: keep R2 in native orientation.

\newpage

## Source: `docs/parity/fastx_uniques_paired_extension.md`

# Paired-end extension for fastx_uniques

This document describes the paired-end extension path for `vsearch --fastx_uniques` implemented in this workspace.

## 1) What is the extension

### Scope

Stock `fastx_uniques` dereplicates one sequence stream. The extension adds a paired-end mode when paired input is provided (second positional input for `R2`, or `--interleaved`), so each unique unit is a pair of anchors:

- left anchor: from forward read (R1)
- right anchor: from reverse read (R2), kept in native read orientation

A pair is considered identical only when both anchors are exactly identical.

### CLI trigger

Paired extension is selected by running `--fastx_uniques` with paired input:

- split paired input: provide `R2` as the second positional input
- interleaved paired input: set `--interleaved` and provide one input stream

Example:

```bash
./bin/vsearch \
  --fastx_uniques reads_R1.fastq.gz reads_R2.fastq.gz \
  --fastaout tav_left.fasta \
  --fastaout_rev tav_right.fasta \
  --tabbedout tav_catalog.tsv
```

Without paired input (no second positional input and no `--interleaved`), command dispatch remains on stock derep behavior.

### Input expectations

- Split mode: R1 and R2 are read in lockstep from the two input streams.
- Interleaved mode: input is consumed as `(R1,R2,R1,R2,...)`; odd record count is an error.
- If one file has extra records, execution fails.
- Anchor length is `min(len(R1), len(R2))`, additionally bounded by `--fastq_trunclen` when provided.

Low-level relationship in current code:

- `fastx_uniques` input call stack:
  ```text
  tav_fastx_uniques(...)
    -> load_paired_records_from_fastx(..., apply_qmask=false)
    -> anchor truncation + pair-key derep
  ```
- `cluster_unoise` / `uchime3_denovo` input call stack (same loader, different mode):
  ```text
  tav_cluster_unoise(...) or tav_uchime3_denovo(...)
    -> load_paired_records_from_fastx(..., apply_qmask=true)
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
- In paired mode, right anchor sequences keep native R2 orientation.

## 4) Practical guidance

- Use paired mode whenever downstream steps are paired-aware.
- Keep `--fastq_trunclen` consistent between runs if reproducibility of anchors is important.
- Use `--notrunclabels` if full Illumina header suffix fields must be preserved in representative labels.

\newpage

## Source: `docs/parity/cluster_unoise_paired_parity.md`

# Cluster UNOISE Parity Notes: Stock vs Paired Extension

This document captures parity between stock `vsearch --cluster_unoise` and the paired-end extension in this workspace.

## Scope and intent

- Stock `cluster_unoise` denoises one sequence stream.
- Paired extension denoises synchronized sequence pairs `(R1, R2)`.
- Design target: preserve stock UNOISE logic where possible, and extend each operation to paired data with explicit, deterministic rules.

## Input and output semantics

### Stock mode: single-end `cluster_unoise`

- Input: FASTA file to `--cluster_unoise`.
- Centroid output semantics: `--centroids` writes denoised centroid FASTA sequences.

### Paired mode: `cluster_unoise` with paired input

- Input: two synchronized sequence files:
  - `--cluster_unoise` for left reads
  - second positional input for right reads (`R2`)
- Alternative input mode: `--interleaved` with one interleaved paired FASTA/FASTQ stream
- Centroid output semantics in paired mode:
  - `--centroids` writes denoised left-centroid FASTA (R1)
  - `--fastaout_rev` writes denoised right-centroid FASTA (R2)
- Optional `--tabbedout` writes paired tabular summary.

## Step-by-step parity map

### Step 1: Input ingestion and ordering

Stock: UNOISE runs on one sequence stream, sorted by decreasing abundance before clustering.
Paired extension: runs on synchronized R1 and R2 streams, then processes paired records by decreasing abundance.
Extension mechanism: each item is a paired record with left and right sequences plus one abundance field.

### Step 2: Candidate generation by k-mer evidence

Stock: for each query sequence, extract unique query k-mers, count k-mer hits against DB/indexed targets, keep high-scoring candidates for downstream filtering/alignment.
Paired extension: build and incrementally maintain a paired centroid k-mer index (left/right postings) and score candidates from that index using stock `unique_count`-style k-mer extraction on both ends, then combine:
score_left = shared k-mers between query.left and centroid.left
score_right = shared k-mers between query.right and centroid.right
score_total = score_left + score_right
Candidate ordering follows stock minheap comparator behavior: score_total higher first, then shorter combined length first, then earlier centroid index first.
Extension mechanism: additive two-end evidence while preserving stock idea of index-driven k-mer pre-screening.

### Step 3: Unaligned pre-filters (explicit parity list)

Stock pre-alignment filters in `search_acceptable_unaligned` and their defaults are:

| Filter / gate | Condition / meaning | Default |
| :-- | :-- | :-- |
| `maxqsize` | `qsize <= maxqsize` | `int_max` (effectively unbounded; `std::numeric_limits<int>::max()`) |
| `mintsize` | `tsize >= mintsize` | `0` |
| `minsizeratio` | `qsize >= minsizeratio * tsize` | `0.0` |
| `maxsizeratio` | `qsize <= maxsizeratio * tsize` | `dbl_max` (effectively unbounded) |
| `minqt` | `qlen >= minqt * tlen` | `0.0` |
| `maxqt` | `qlen <= maxqt * tlen` | `dbl_max` (effectively unbounded) |
| `minsl` | `shorter_len >= minsl * longer_len` | `0.0` |
| `maxsl` | `shorter_len <= maxsl * longer_len` | `dbl_max` (effectively unbounded) |
| `idprefix` | first `N` bases must match exactly | `0` (disabled) |
| `idsuffix` | last `N` bases must match exactly | `0` (disabled) |
| `self` | reject same header | `0` (disabled) |
| `selfid` | reject perfect self sequence | `0` (disabled) |

Paired unaligned extension preserves the stock pre-alignment filter semantics, but replaces single-end abundance and length variables with paired aggregated values in `paired_unaligned_filters_pass`:

Paired aggregated variables used by filtering

| Paired metric | Definition |
| :-- | :-- |
| `qsize_pair` | `query.abundance` |
| `tsize_pair` | `target.abundance` |
| `qlen_pair` | `len(query.left) + len(query.right)` |
| `tlen_pair` | `len(target.left) + len(target.right)` |
| `shorter_pair` | `min(qlen_pair, tlen_pair)` |
| `longer_pair` | `max(qlen_pair, tlen_pair)` |

Filtering behavior

The following pre-alignment filters are the same as in the stock pre-alignment filter table, but evaluated on the paired aggregated variables above:

- `maxqsize`
- `mintsize`
- `minsizeratio`
- `maxsizeratio`
- `minqt`
- `maxqt`
- `minsl`
- `maxsl`

Low-level relationship in current code:

- Stock Step 3 call stack:
  ```text
  search_acceptable_unaligned(...)
    -> search_unaligned_numeric_filters_pass(...)
    -> single-end idprefix/idsuffix/self/selfid checks
  ```
- Paired Step 3 call stack:
  ```text
  paired_unaligned_filters_pass(...)
    -> search_unaligned_numeric_filters_pass(...)
    -> paired idprefix(R1 prefix)/idsuffix(R2 prefix)/self/selfid checks
  ```

### Paired-specific differences

| Filter | Paired-specific rule |
| :-- | :-- |
| `idprefix` | checked only on the R1 prefix (first `N` bases of left read) |
| `idsuffix` | checked only on the R2 prefix (first `N` bases of right read; treated as merged-sequence suffix analog) |
| `self`     | reject by checking if paired headers are identical        |
| `selfid`   | reject by checking if both paired sequences are identical |


Relevant code pointers for Step 3:

- stock: `src/searchcore.cc` -> `search_acceptable_unaligned(...)`
- shared low-level numeric kernel: `src/searchcore.cc` -> `search_unaligned_numeric_filters_pass(...)`
- paired ext: `src/tav_extension.cc` -> `paired_unaligned_filters_pass(...)`

### Step 4: Delayed alignment batching behavior

Stock: candidates that pass unaligned filters are queued and aligned in delayed batches (MAXDELAYED), not aligned one-by-one immediately.
Paired extension: same staging pattern; candidate hits are delayed then aligned in batch-processing loop.
Extension mechanism: batching kept; each delayed hit triggers two end-alignments.

### Step 5: Alignment statistics used in UNOISE acceptance

Stock: alignment-derived stats are used; UNOISE distance uses mismatch count from alignment, ignoring gaps.
Paired extension: each end is globally aligned separately; mismatch counts are extracted per end; paired UNOISE distance is:
d_pair = mismatches_left + mismatches_right
gaps are not added into d_pair for UNOISE acceptance
Extension mechanism: direct sum of two stock-style mismatch counts, one per end.

### Step 6: Aligned filters (explicit parity list)

Stock aligned filters in search_acceptable_aligned are:

| Aligned filter | Condition / meaning | Default |
| :-- | :-- | :-- |
| `weak_id` | `id >= 100 * weak_id` | `10.0` (percent) |
| `maxsubs` | `mismatches <= maxsubs` | `int_max` (effectively unbounded) |
| `maxgaps` | `internal_gaps <= maxgaps` | `int_max` (effectively unbounded) |
| `mincols` | `internal_alignmentlength >= mincols` | `0` |
| `leftjust` | no left terminal gaps | `0` (disabled) |
| `rightjust` | no right terminal gaps | `0` (disabled) |
| `query_cov` | `(matches + mismatches) >= query_cov * qlen` | `0.0` |
| `target_cov` | `(matches + mismatches) >= target_cov * tlen` | `0.0` |
| `maxid` | `id <= 100 * maxid` | `1.0` (used as `100 * maxid`, so default upper bound is `100%`) |
| `mid` | `100 * matches / (matches + mismatches) >= mid` | `0.0` |
| `maxdiffs` | `mismatches + internal_indels <= maxdiffs` | `int_max` (effectively unbounded) |

Paired aligned extension preserves the stock aligned-filter semantics, but replaces single-end counts with paired aggregated totals from the left and right end alignments.

Paired aggregated variables used by filtering

| Paired metric | Definition |
| :-- | :-- |
| `mismatches_total` | `left.mismatches + right.mismatches` |
| `nwgaps_total` | `left.nwgaps + right.nwgaps` |
| `nwalignment_cols_total` | `left.nwalignmentlength + right.nwalignmentlength` |
| `internal_alignment_cols_total` | `left.internal_alignmentlength + right.internal_alignmentlength` |
| `internal_gaps_total` | `left.internal_gaps + right.internal_gaps` |
| `internal_indels_total` | `left.internal_indels + right.internal_indels` |
| `matches_total` | `left.matches + right.matches` |
| `ungapped_cols_total` | `matches_total + mismatches_total` |
| `query_len_pair` | `len(query.left) + len(query.right)` |
| `target_len_pair` | `len(target.left) + len(target.right)` |
| `shortest_pair` | `min(query_len_pair, target_len_pair)` |
| `longest_pair` | `max(query_len_pair, target_len_pair)` |

Paired combined identity metrics

| Identity / metric | Definition |
| :-- | :-- |
| `id0` | `100 * matches_total / shortest_pair` |
| `id1` | `100 * matches_total / nwalignment_cols_total` |
| `id2` | `100 * matches_total / internal_alignment_cols_total` |
| `id3` | `max(0, 100 * (1 - (mismatches_total + nwgaps_total) / longest_pair))` |
| `id4` | `id1` |
| `hit.id` | `id{opt_iddef}`; default `opt_iddef = 2`, so default identity is `id2` |
| `hit.mid` | `100 * matches_total / ungapped_cols_total` when denominator > 0 |

Filtering behavior

The following aligned filters are the same as in the stock aligned-filter table, but evaluated on the paired aggregated variables above:

- `weak_id`
- `maxsubs`
- `maxgaps`
- `mincols`
- `query_cov`
- `target_cov`
- `maxid`
- `mid`
- `maxdiffs`

Low-level relationship in current code:

- Stock Step 6 call stack:
  ```text
  search_acceptable_aligned(...)
    -> search_aligned_compute_identity_metrics(...)
    -> search_aligned_threshold_filters_pass(...)
    -> stock accept/weak decision (UNOISE skew-beta branch for cluster_unoise)
  ```
- Paired Step 6 call stack:
  ```text
  paired_aligned_filters_pass(...)
    -> aggregate R1/R2 alignment totals
    -> search_aligned_compute_identity_metrics(...)
    -> search_aligned_threshold_filters_pass(...)
    -> paired accept/weak labeling
  ```

Paired-specific differences

| Filter | Paired-specific rule |
| :-- | :-- |
| `leftjust` | `(left.trim_q_left + left.trim_t_left + right.trim_q_left + right.trim_t_left) == 0` |
| `rightjust` | `(left.trim_q_right + left.trim_t_right + right.trim_q_right + right.trim_t_right) == 0` |

Terminal gap trims are the leading and trailing gap runs from each end-alignment, produced by the stock `align_trim` logic: `trim_q_left`, `trim_t_left`, `trim_q_right`, and `trim_t_right`.

Relevant code pointers for Step 6:

- stock: `src/searchcore.cc` -> `search_acceptable_aligned(...)`
- shared low-level metric kernel: `src/searchcore.cc` -> `search_aligned_compute_identity_metrics(...)`
- shared low-level threshold kernel: `src/searchcore.cc` -> `search_aligned_threshold_filters_pass(...)`
- stock terminal-gap trimming: `src/searchcore.cc` -> `align_trim(...)`
- paired ext: `src/tav_extension.cc` -> `paired_aligned_filters_pass(...)`
- paired per-end alignment backend feeding these metrics: `src/tav_extension.cc` -> `align_one_end_stock_style(...)`

### Step 7: UNOISE skew-beta acceptance rule

Stock: for cluster_unoise, after aligned filters:
skew = q_abundance / target_abundance
beta = 2^(-1 - alpha * mismatches)
Accept if skew <= beta, or mismatches == 0.
Paired extension: same formula and same acceptance logic, with:
mismatches replaced by d_pair = mismatches_left + mismatches_right
Extension mechanism: formula unchanged, distance generalized to paired sum.

### Step 8: Weak-hit bookkeeping and early-stop behavior

Stock:
accepted hits are marked accepted
hits passing aligned filters but failing UNOISE skew-beta are marked weak
reject/accept counters drive stopping via maxaccepts/maxrejects

Paired extension:
same accepted/weak/rejected bookkeeping
same maxaccepts/maxrejects-driven early stop behavior
Extension mechanism: identical control-flow pattern, but per-candidate alignment work is paired.

Low-level ext control flow details (`tav_cluster_unoise` main loop):

- Candidate generation yields ordered centroid candidates and starts local counters:
  - `accepts = 0`
  - `rejects = 0`
  - delayed queue `delayed` (size-bounded by `MAXDELAYED`)
- Unaligned stage:
  - if `paired_unaligned_filters_pass(...)` fails, candidate is immediately marked:
    - `hit.rejected = true`
    - `hit.weak = false`
    - `rejects++`
  - if unaligned checks pass, candidate enters `delayed`.
- Delayed aligned stage (`process_delayed`):
  - before each hit, short-circuit if:
    - `accepts >= maxaccepts_effective`, or
    - `rejects >= maxrejects_effective`
  - run `paired_aligned_filters_pass(...)`:
    - if false: mark hard reject (`rejected=true`, `weak=false`), `rejects++`
  - if aligned filters pass, evaluate UNOISE skew-beta:
    - `skew = q.abundance / centroid.abundance`
    - `beta = 2^(-1 - alpha * mismatches_total)`
    - accept when `(skew <= beta) or (mismatches_total == 0)`
      - mark `accepted=true`, `weak=false`, `accepts++`
    - otherwise weak reject
      - mark `rejected=true`, `weak=true`, `rejects++`
- Early-stop semantics:
  - scanning candidates stops once accepts/rejects limits are reached.
  - delayed hits that were never finalized due early stop are ignored.
- This is why Step 8 parity is about state transitions and stopping rules, not just final labels.

Relevant code pointers for Step 8:

- paired ext orchestration: `src/tav_extension.cc` -> `tav_cluster_unoise(...)`
  - `process_delayed` lambda (weak-hit bookkeeping + early stop)
  - candidate loop guards using `maxaccepts_effective`/`maxrejects_effective`
- stock counterpart behavior:
  - `src/searchcore.cc` -> `align_delayed(...)`
  - `src/searchcore.cc` -> `search_acceptable_aligned(...)`

### Step 9: Best-hit selection and centroid update

Stock: best accepted candidate for UNOISE path is selected by bysize comparator family (abundance first, then identity, then target index), then query is assigned to existing centroid or starts a new centroid.
Paired extension: among accepted paired hits, best is chosen with the same ordering (abundance first, then identity, then centroid index), then query abundance is merged or new centroid is created.
Extension mechanism: same centroid-update semantics, paired metrics for ranking.

### Step 10: Output semantics for denoised centroids

Stock: centroids output means denoised FASTA centroid sequences.
Paired extension (now aligned to your requested semantics):
centroids means denoised paired centroid FASTA output:
centroids path writes R1 centroid FASTA
fastaout_rev writes corresponding R2 centroid FASTA
Extension mechanism: centroid concept lifted from single sequence to paired sequence tuple.


## Operational constraints in paired mode

- If `--fastaout` is provided, `--fastaout_rev` must also be provided.
- If `--centroids` is provided, `--fastaout_rev` must also be provided.
- At least one output path must be requested among:
  - paired FASTA (`--fastaout` + `--fastaout_rev`)
  - paired centroids (`--centroids` + `--fastaout_rev`)
  - `--tabbedout`


## Parity status summary

- Current status: cluster_unoise paired extension is now on the stock backbone for candidate generation, delayed alignment, aligned/unaligned filters, UNOISE skew-beta acceptance, and bysize best-hit tie-breaking.
- Intentional paired-only extensions are:
  - two-end aggregation semantics
  - no extra loader-level pair-drop gate from `--minseqlength`/`--maxseqlength`/`--filter`
  - paired output requirements (`--fastaout_rev` companion output)


## Example paired command

```bash
./bin/vsearch \
  --cluster_unoise uniques_r1.fasta uniques_r2.fasta \
  --centroids denoised_r1.fasta \
  --fastaout_rev denoised_r2.fasta \
  --tabbedout denoised_pairs.tsv
```

\newpage

## Source: `docs/parity/uchime3_denovo_paired_parity.md`

# UCHIME3-Denovo Parity Notes: Stock vs Paired Extension

This document captures parity status between stock `vsearch --uchime3_denovo` and the paired-end extension in this workspace.

## Scope and intent

- Stock `uchime3_denovo` evaluates one denoised sequence stream.
- Paired extension evaluates synchronized denoised sequence pairs `(R1, R2)`.
- Design target: preserve stock command intent (abundance-aware de novo chimera detection), while extending to paired data with explicit breakpoint logic.

## Input and output semantics

### Stock mode: single-end `uchime3_denovo`

- Input: denoised FASTA file to `--uchime3_denovo`.
- Outputs: single-end chimera/nonchimera/report outputs.

### Paired mode: `uchime3_denovo` with paired input

- Input: two synchronized sequence files:
  - `--uchime3_denovo` for left reads
  - second positional input for right reads (`R2`)
- Alternative input mode: `--interleaved` with one interleaved paired FASTA/FASTQ stream
- Input preprocessing reuses paired loader behavior:
  - synchronized paired FASTX loading and record-count guards
  - qmask/hardmask preprocessing on both ends
- Current paired outputs:
  - `--tabbedout`: paired-extension report (`parent_a`, `parent_b`, `breakpoint_class`, scores, delta)
  - `--uchimeout`: stock-shaped UCHIME tabular output adapted to paired scoring
  - `--uchimealns`: paired chimera alignment summaries (left/right context)
  - `--borderline`: accepted for CLI parity (typically empty with uchime3 logic)
  - `--nonchimeras` + `--nonchimeras2`: paired FASTA outputs (split `R1` and `R2`)
  - `--chimeras` + `--chimeras2`: paired FASTA outputs (split `R1` and `R2`)
  - `--nonchimeras_tsv`: paired TSV catalog rows for non-chimeras
  - `--chimeras_tsv`: paired TSV catalog rows for chimeras

### Clarifications: catalog input and masking

- Catalog input note:
  - Stock `uchime3_denovo` does not accept TSV catalog input.
  - Paired `uchime3_denovo` also only accepts FASTA/FASTQ pairs (split positional `R1/R2` input or `--interleaved`).
  - `--chimeras_tsv` / `--nonchimeras_tsv` are paired output formats, not input formats.

- `qmask` / `hardmask` note:
  - `--qmask dust` runs DUST low-complexity masking on each end before candidate scoring/classification.
  - In this implementation, DUST scans 64 nt windows (stepping by 32 nt), computes a low-complexity score from repeated 3-mer composition, and masks windows whose score exceeds an internal threshold (`dust_level = 20`).
  - With soft masking (`--qmask dust` and no `--hardmask`), masked positions are lowercased.
  - With hard masking (`--hardmask` active), masked/lowercase positions are converted to `N`.
  - Masked positions are then excluded from downstream k-mer matching behavior (`MASK_SOFT` path treats lowercase as masked; hardmasked `N` is also excluded as ambiguous).

## Step-by-step parity map

### Step 1: Command dispatch and hard guards

Stock implementation: `cmd_chimera()` routes plain `--uchime3_denovo` (without paired input) to `chimera()`. Parameter guards enforce `abskew >= 1`, `xn > 1`, and `dn > 0`.

Paired extension implementation: `cmd_chimera()` detects paired input (`R2` as second positional input or `--interleaved`) and routes to `tav_uchime3_denovo()`. It enforces the same `abskew/xn/dn` guards, enforces split FASTA pairing (`--chimeras` with `--chimeras2`, `--nonchimeras` with `--nonchimeras2`), and requires at least one output among paired FASTA, TSV catalog, `--tabbedout`, `--uchimeout`, `--uchimealns`, or `--borderline`.

Extension mechanism: paired mode is an explicit alternate execution path, not a flag inside the stock `chimera()` loop.

### Step 2: Input loading, masking, and deterministic processing order

Stock implementation: `chimera()` reads the denoised FASTA from `--uchime3_denovo`, optionally applies dust/hardmask (`opt_qmask` path), and then calls `db_sortbyabundance()`. Query processing order is abundance descending.

Paired extension implementation: `tav_uchime3_denovo()` calls `load_paired_records_from_fastx()` on paired FASTX input (split positional `R1/R2` or interleaved via `--interleaved`). The loader enforces synchronized record counts, applies the same qmask/hardmask operations to both ends. Records are then sorted by:
1. abundance descending
2. header ascending
3. left sequence ascending
4. right sequence ascending

Extension mechanism: each processing unit is a `TavRecord` pair `(left,right,abundance,header)` rather than a single sequence.

### Step 3: Parent search space construction and growth rule

Stock implementation: in denovo mode, `dbindex_prepare()` initializes an initially empty parent index. `chimera()` sets `opt_self=1`, `opt_selfid=1`, and `opt_maxsizeratio = 1/abskew` so candidate parents must be sufficiently more abundant. After each query is classified as non-chimeric (`status < suspicious`), `dbindex_addsequence(seqno, opt_qmask)` inserts that sequence for future queries.

Paired extension implementation: `tav_uchime3_denovo()` now reuses stock dbindex primitives directly:
- `dbindex_prepare(0, opt_qmask)` initializes one stock k-mer index over the in-memory paired sequence store (left and right ends are stored as alternating sequence slots).
- `parent_pool_indices` tracks non-chimeric parent record IDs in abundance order.
- After each non-chimeric query, `dbindex_addsequence(left_seqno, opt_qmask)` and `dbindex_addsequence(right_seqno, opt_qmask)` are called to incrementally expose that pair as future parent material.

Extension mechanism: both modes use a monotonic non-chimeric parent pool and incremental `dbindex_addsequence` growth; paired mode keeps end-specific semantics by mapping left/right ends onto deterministic sequence-slot parity (even = R1, odd = R2).

### Step 4: Candidate enumeration for the current query

Stock implementation: each query is split into `parts=4` (`partition_query()` for uchime/uchime2/uchime3). `search_onequery()` plus `search_joinhits()` collects accepted hits from each part, deduplicates targets, and builds `cand_list`. Full-query global alignments are then computed for every candidate (`search16`, with linear-memory fallback when SIMD overflows).

Paired extension implementation:
- Query kmers are looked up through stock dbindex match lists (`dbindex_getmatchlist`/`dbindex_getmatchcount`/`dbindex_getmapping`).
- R1 k-mers only contribute from even mapped sequence slots (left-end parents), and R2 k-mers only contribute from odd mapped sequence slots (right-end parents), then both contributions are summed to a pair-level parent score.
- Candidate parents are scored by summed left/right k-mer overlap.
- The abundance gate (`p.abundance >= abskew * q.abundance`) is enforced.
- Top candidates are retained with stock-like bounded heap behavior (`maxaccepts + maxrejects` style cap).
- Full-query alignments are then computed for each retained candidate to build match vectors used by parent selection.

Extension mechanism: this is now a paired analog of stock candidate discovery + full-query alignment materialization, rather than direct iteration over all eligible parents.

### Step 5: One-parent baseline (`best_one`)

Stock implementation: there is no explicit pre-parent-search `best_one` variable in `chimera.cc`. The closest stock counterpart is `QT` in `eval_parents()`, where:
- `QA` = identity(query, parent A)
- `QB` = identity(query, parent B)
- `QT = max(QA, QB)`
This one-parent baseline exists only after two parents have already been selected/evaluated. If two parents are not found (`find_best_parents() == 0`), stock returns `Status::no_parents` and emits no-parent outputs (placeholder parent fields in `uchimeout`).

Paired extension implementation: `best_one_score` is now computed from the selected two-parent evaluation itself as:
`best_one = max(match_QA, match_QB)` (count scale), mirroring stock's `QT = max(QA,QB)` baseline semantics within the same evaluation stage.
For no-parent cases, `best_one` falls back to full query-length count for reporting continuity.

Extension mechanism: the paired path now aligns with stock by deriving one-parent baseline from the evaluated parent pair, rather than running a separate pre-model one-parent optimization sweep.

### Step 6: Two-parent preselection (stock-like 32 bp voting)

Stock implementation detail (`find_matches()` + `find_best_parents()`):
- Build `match[candidate][qpos]` from the full-query alignment CIGAR:
  - `1` if aligned query/parent symbols overlap at `qpos`
  - `0` otherwise
- For each candidate, compute `smooth[candidate][qpos]` = number of matches in the trailing 32 bp window ending at `qpos`.
- For each `qpos`, compute `maxsmooth[qpos]` = best window score among candidates not yet selected.
- "Win" definition: a candidate wins position `qpos` if `smooth[candidate][qpos] == maxsmooth[qpos]` and `maxsmooth[qpos] > 0`.
- Parent 1 = candidate with maximum total wins.
- Before picking parent 2, for every `qpos` where parent 1 was tied for best, zero all candidates' `match` entries inside that winning 32 bp window (stock's way to remove already-explained windows).
- Recompute smoothing/wins on the masked matrix and pick parent 2.

Paired extension implementation detail:
- Use concatenated query axis:
  `[left query positions in forward order][right query positions in forward order]`.
- Build candidate match vectors from full-query alignment CIGARs (paired analog of stock `find_matches()`).
- Run the same 32 bp smoothing + per-position winner counting + winner-window wipeout to choose parent 1 then parent 2.
- Parent 1/2 winner selection now calls the same low-level helper (`select_best_two_parents_from_match_matrix`) from both stock `find_best_parents()` and paired Step 6.
- No fallback all-pairs sweep is performed when no two parents are found; behavior is stock-like no-parent handling.

Low-level relationship in current code:

- Stock Step 6 call stack:
  ```text
  find_matches(...)
    -> find_best_parents(...)
    -> select_best_two_parents_from_match_matrix(...)
  ```
- Paired Step 6 call stack:
  ```text
  compute_match_vector(...)
    -> parse_cigar_operations_from_string(...)
    -> parse_cigar_string(...)
    -> select_best_parent_pair(...)
    -> select_best_two_parents_from_match_matrix(...)
  ```

Extension mechanism: paired parent preselection now uses alignment-derived match vectors and shared stock winner-window selection logic rather than direct character-equality heuristics or all-pairs fallback.

### Step 7: Diff coding and `h` optimization for two-parent models

Stock implementation: `eval_parents()` builds aligned strings for query and two parents, marks ignored columns (gap columns, neighbors of gaps, and ambiguous symbols), then emits `diffs` symbols (`A`, `B`, `N`, `?`, or space). It scans breakpoint position `i` across the alignment and tracks left/right counts:
- `Y` (yes) votes
- `N` (no) votes
- `A` (abstain/other)
Here, `Y/N/A` are vote roles (support, contradiction, abstain), not nucleotide letters.
It computes:
`left_h = left_y / (xn*(left_n + dn) + left_a)`
`right_h = right_y / (xn*(right_n + dn) + right_a)`
`h = left_h * right_h`
and also checks the reverse orientation (swap yes/no) when anti-chimera orientation is stronger. Best `h` determines the model.

Paired extension implementation:
- For each selected parent pair, both ends are globally aligned and converted into unified gapped strings (`Q`, `A`, `B`) using stock-like insertion normalization.
- Ambiguous and gap-adjacent columns are marked ignored exactly in stock style.
- `diffs` are emitted with stock symbols (`A/B/N/?/space`), with reverse-orientation swap handling (`A<->B`) when anti-chimera orientation is stronger.
- `h` optimization scans one concatenated paired axis with right-end scan direction reversed (left: start->end, right: end->start), yielding one `best_i`.
- Paired breakpoint class is then derived from `best_i` location:
  - in left end -> `LEFT_BREAK`
  - at end boundary -> `MIDDLE_BREAK`
  - in right end -> `RIGHT_BREAK`

Low-level relationship in current code:

- Paired Step 7 call stack:
  ```text
  evaluate_parent_pair(...)
    -> build_end_alignments(...)
    -> parse_cigar_operations_from_string(...)
    -> parse_cigar_string(...)
    -> stock-style diff/vote/model scoring
  ```

Extension mechanism: the paired path now uses stock-style diff coding and one-pass breakpoint scan; paired class labels are metadata derived from the chosen breakpoint position.

### Step 8: Final chimeric decision rule

Stock implementation (`uchime2/3` path in `eval_parents()`): after best model selection, define:
- `cols`: number of non-ignored alignment columns
- `match_QM`: count of columns where query matches the selected two-parent model
- `QA`: `100 * match_QA / cols` (query vs parent A identity)
- `QB`: `100 * match_QB / cols` (query vs parent B identity)
- `QT = max(QA, QB)` (best one-parent identity)
- `QM = 100 * match_QM / cols` (two-parent model identity)
`uchime3_denovo` marks chimeric when the model is perfect and best single parent is imperfect:
`match_QM == cols` and `QT < 100.0`.

Paired extension implementation: after selecting/evaluating two parents on the concatenated paired axis, it uses the same uchime3 condition in paired form:
- `match_QM`: non-ignored columns explained by the selected two-parent model
- `cols`: total non-ignored columns
- `QT = max(QA, QB)`: best one-parent identity
Classify as chimeric iff:
- `match_QM == cols`
- `QT < 100.0`

Extension mechanism: parity is kept in the decision rule itself; paired breakpoint class (`LEFT_BREAK`/`MIDDLE_BREAK`/`RIGHT_BREAK`) is retained as annotation in paired reports, not as an extra classification gate.

### Step 9: Parent-pool update after classification

Stock implementation: all non-chimeric outcomes (`status < suspicious`) are added to the denovo parent dbindex, enabling incremental parent search for later, lower-abundance queries.

Paired extension implementation: all non-chimeric pairs are appended to `parent_pool_indices`; chimeric pairs are never added.

Extension mechanism: identical growth policy, but parent units are paired records.

### Step 10: Output materialization

Stock implementation file/schema summary:
- `--chimeras`, `--nonchimeras`: single-end FASTA (one record per query sequence).
- `--uchimeout`: stock UCHIME tabular schema.
  - default form includes columns equivalent to:
    `h, query, parentA, parentB, top, QM, QA, QB, AB, QT, LY, LN, LA, RY, RN, RA, divdiff, YN`
  - with `--uchimeout5`, the `top` column is omitted (stock compatibility mode).
- `--uchimealns`: stock textual alignment report with aligned query/parents plus `Diffs`, `Votes`, and `Model` lines.
- `--borderline`: FASTA output for borderline calls (normally unused in uchime2/3 decision path).

Paired extension comparison against stock:
- Same option name, mostly same schema family:
  - `--uchimeout` keeps stock-like column semantics (including `--uchimeout5` behavior), but values come from paired scoring/counts.
- Same option name, paired-adapted body format:
  - `--uchimealns` now includes paired alignment blocks with stock-style `Diffs/Votes/Model` lines split by side (`*_LEFT`, `*_RIGHT`), while keeping paired context labels.
- Same option name, different FASTA layout:
  - `--chimeras`/`--chimeras2` and `--nonchimeras`/`--nonchimeras2` are split paired FASTA outputs (`R1` file and `R2` file).
- Same option name, effectively empty in paired uchime3:
  - `--borderline` accepted for CLI parity but not populated by current paired decision path.
- Extra paired-only outputs not present in stock:
  - `--tabbedout` paired report with 7 columns:
    `query_tav, parent_a, parent_b, breakpoint_class, best_one_score, best_two_score, delta`
  - `--chimeras_tsv` / `--nonchimeras_tsv` paired catalog TSV with 5 columns:
    `tav_id, abundance, header, left_anchor, right_anchor`

Extension mechanism: paired outputs carry more metrics than stock because the paired algorithm has extra latent variables (two ends + breakpoint class) that do not exist in stock's single-axis model.

## Operational constraints in paired mode

- Paired extension is entered when paired input is provided (`R2` as second positional input or `--interleaved`).
- At least one output must be specified among:
  - `--nonchimeras` + `--nonchimeras2`
  - `--chimeras` + `--chimeras2`
  - `--tabbedout`
  - `--uchimeout`
  - `--uchimealns`
  - `--borderline`
  - `--chimeras_tsv`
  - `--nonchimeras_tsv`

## Parity status summary

- Current status: paired `uchime3_denovo` follows paired FASTX ingestion routing (split positional input or `--interleaved`) and removes catalog-mode input.
- Current status: paired mode now mirrors key stock denovo behavior for parent-pool growth (nonchimera-only), index-backed candidate discovery, 32bp winner-vote parent preselection, and stock-style one-pass `h` breakpoint optimization.
- Current algorithm status: paired chimera core now follows stock diff/vote/model mechanics on a concatenated paired axis; paired breakpoint classes are retained as derived annotations.
- Main open parity work items:
  - calibrate/validate paired classification behavior against broader real datasets
  - decide whether to further tighten candidate discovery to an even closer `search_onequery()`/`search_joinhits()` analogue
  - continue harmonizing paired `--uchimealns` body formatting with stock text layout where practical

## Example paired command

```bash
./bin/vsearch \
  --uchime3_denovo denoised_r1.fasta denoised_r2.fasta \
  --tabbedout chim_report.tsv \
  --uchimeout chim_uchime.tsv \
  --uchimealns chim_alns.txt \
  --nonchimeras_tsv nonchim_pairs.tsv \
  --chimeras_tsv chim_pairs.tsv \
  --nonchimeras nonchim_pairs.fa \
  --nonchimeras2 nonchim_pairs_r2.fa \
  --chimeras chim_pairs.fa \
  --chimeras2 chim_pairs_r2.fa
```

\newpage

## Source: `docs/parity/usearch_global_paired_parity.md`

# Paired `usearch_global` Parity (Stock VSEARCH vs Extension)

This document captures parity status between stock `vsearch --usearch_global` and the paired-end extension path in this workspace.

## Scope and intent

- Command: `vsearch --usearch_global` with paired input (second positional input for `R2`, or `--interleaved`).
- Design target: preserve stock search backbone semantics, then extend only what is necessary to treat one read pair as one query unit.

## Current paired-mode contract

### Required paired inputs

- Query input:
  - split mode: `--usearch_global` (R1) + second positional input (R2)
  - interleaved mode: `--usearch_global` with `--interleaved`
- Database input:
  - one interleaved paired FASTA/FASTQ file via `--db` (left record then right record), or
  - split paired FASTA/FASTQ files via `--db` (left) + `--db2` (right).
- Catalog TSV input is intentionally removed for paired `usearch_global` parity.

### Supported paired outputs

- `--alnout`
- `--userout`
- `--uc`
- `--blast6out`
- `--matched` + `--match2` (split paired outputs; required together)
- `--notmatched` + `--notmatched2` (split paired outputs; required together)
- `--dbmatched` + `--dbmatched2` (split paired outputs; required together)
- `--dbnotmatched` + `--dbnotmatched2` (split paired outputs; required together)
- `--otutabout`, `--biomout`, `--mothur_shared_out`

### Not supported in paired mode (current)

- `--fastapairs`
- `--qsegout`
- `--tsegout`
- `--samout`
- `--lcaout`
- Custom `--userfields`

## Step-by-step parity map

### Step 1: Command dispatch and mode guards

Stock implementation: `cmd_usearch_global()` routes non-paired runs to `usearch_global()`. It requires `--db`, requires `--id` in `[0,1]`, and accepts the full stock output set (including SAM/segment/alignment-rich outputs).

Paired extension implementation: `cmd_usearch_global()` routes paired input (`R2` positional input or `--interleaved`) to `tav_usearch_global()`. It enforces:
- `--strand plus` only
- no `--fastapairs`, `--qsegout`, `--tsegout`, `--samout`, `--lcaout`
- no custom `--userfields`
- split-output pairing constraints (paired FASTA outputs require both `R1` and `R2` file arguments)
- at least one paired-supported output
- `--db` and valid `--id`

Extension mechanism: paired mode is a dedicated execution path with explicit CLI constraints.

### Step 2: Database ingestion and target construction

Stock implementation: `search_prep()` loads DB from UDB or FASTA/FASTQ, applies optional DB masking (`dust`/`hardmask`), then builds `dbindex` for k-mer lookup (`dbindex_prepare`, `dbindex_addallsequences`).

Paired extension implementation: `tav_usearch_global()` loads paired targets from:
- interleaved `--db` (left/right alternating), or
- split `--db` (left) + `--db2` (right)
Catalog TSV input is rejected. Each DB pair is converted to one `TavRecord`:
- anchor length `anchor_len = min(left_len, right_len)` (or bounded by `--fastq_trunclen` via `get_anchor_len`)
- min/max length filtering with `--filter any|both`
- DB masking on both ends (`opt_dbmask`)
- header base normalized by removing trailing `/1` or `/2`
- abundance taken from left record (warning on mismatch)

Extension mechanism: one paired target record replaces one single-sequence target.

### Step 3: Query ingestion and paired normalization

Stock implementation: each query sequence is read once (or twice with reverse strand search when enabled), optionally query-masked, then searched.

Paired extension implementation: split queries (`R1`/`R2`) are read synchronously, or interleaved queries are consumed as `(R1,R2,R1,R2,...)`. For each pair:
- right read is kept in native R2 orientation
- both ends are truncated to `anchor_len`
- query masking (`opt_qmask`) is applied on both ends
- abundance is `max(fwd_abundance, 1)`
- query unit becomes one paired `TavRecord`

Extension mechanism: paired search operates directly on native-orientation left/right anchors.

### Step 4: K-mer candidate scoring and top-hit preselection

Stock implementation: `search_onequery()` extracts unique query k-mers (`unique_count`), accumulates DB k-mer hit counts through `dbindex`, and keeps top candidates via min-heap (`search_topscores`).

Paired extension implementation:
- builds paired postings once (`left_postings`, `right_postings`) from DB k-mers
- for each query, computes `db_kmer_scores[i] = left_overlap_i + right_overlap_i`
- keeps candidates with `score >= min(opt_minwordmatches, qk_left + qk_right)`
- uses the same heap comparator behavior (higher k-mer count first, then shorter target length, then earlier index)

Extension mechanism: stock k-mer backbone is preserved with additive left+right evidence.

### Step 5: Unaligned filtering and delayed-alignment staging

Stock implementation: candidate loop applies `search_acceptable_unaligned()`, pushes survivors into delayed alignment batches (`MAXDELAYED`), and stops by `maxaccepts/maxrejects` limits.

Paired extension implementation: mirrors this control flow with:
- `paired_unaligned_filters_pass()` (paired analog of stock unaligned filters:
  `maxqsize`, `mintsize`, `minsizeratio`, `maxsizeratio`, `minqt`, `maxqt`, `minsl`, `maxsl`, `idprefix`, `idsuffix`, `self`, `selfid`)
  - `idprefix` is checked on R1 prefix only
  - `idsuffix` is checked on R2 prefix only (merged-suffix analog under paired orientation contract)
- pending batch processed at `MAXDELAYED` or loop end
- effective limits:
  - `maxaccepts_effective`
  - `maxrejects_effective`
  - `maxhits_considered = maxaccepts + maxrejects - 1`

Low-level relationship in current code:

- Stock Step 5 call stack:
  ```text
  search_acceptable_unaligned(...)
    -> search_unaligned_numeric_filters_pass(...)
    -> single-end idprefix/idsuffix/self/selfid checks
  ```
- Paired Step 5 call stack:
  ```text
  paired_unaligned_filters_pass(...)
    -> search_unaligned_numeric_filters_pass(...)
    -> paired idprefix(R1 prefix)/idsuffix(R2 prefix)/self/selfid checks
  ```

Extension mechanism: delayed-batch semantics are kept, but filters evaluate paired totals/per-end equivalents.

### Step 6: Alignment, aligned filters, and accept/weak/reject assignment

Stock implementation: delayed hits are globally aligned (`search16` SIMD batch, with linear-memory fallback), then `align_trim()` computes terminal/internal stats, `search_acceptable_aligned()` applies aligned thresholds, and hits are classified:
- accepted if `id >= 100 * opt_id`
- weak if aligned-filters pass but `id < 100 * opt_id`
- rejected otherwise

Paired extension implementation: each pending hit is aligned per end with `align_one_end_stock_style()` (linear-memory aligner), then aggregated:
- mismatches/gaps/columns/indels summed across ends
- identity metrics `id0..id4` computed and selected by `opt_iddef`
- `mid` and coverage metrics computed on combined totals
- aligned filter family applied in paired form:
  `weak_id`, `maxsubs`, `maxgaps`, `mincols`, `leftjust`, `rightjust`, `query_cov`, `target_cov`, `maxid`, `mid`, `maxdiffs`
- final accepted/weak/rejected decision uses the same `opt_id` cutoff semantics as stock

Extension mechanism: accept/weak/reject logic is stock-equivalent after two-end stat aggregation.

Low-level relationship in current code:

- Stock Step 6 call stack:
  ```text
  search_acceptable_aligned(...)
    -> search_aligned_compute_identity_metrics(...)
    -> search_aligned_threshold_filters_pass(...)
    -> stock accept/weak labeling
  ```
- Paired Step 6 call stack:
  ```text
  paired_aligned_filters_pass(...)
    -> aggregate R1/R2 alignment totals
    -> search_aligned_compute_identity_metrics(...)
    -> search_aligned_threshold_filters_pass(...)
    -> paired accept/weak labeling
  ```

### Step 7: Hit retention, ordering, and report window

Stock implementation: `search_joinhits()` keeps only accepted/weak hits and sorts by stock comparator:
1. accepted before weak
2. aligned before non-aligned (if present)
3. higher identity
4. lower target index
Reporting applies `--maxhits`, `--top_hits_only`, and `--uc_allhits`.

Paired extension implementation: same retention and ordering semantics over `TavHit`, then:
- `toreport = min(opt_maxhits, report_hits.size())`
- `--top_hits_only` truncates below top identity
- UC output respects `--uc_allhits`

Extension mechanism: stock ranking semantics preserved on paired hit objects.

### Step 8: Per-query output emission and match accounting

Stock implementation: writes configured outputs (`alnout`, `userout`, `blast6`, `uc`, etc.); query is matched if any kept hit exists; `dbmatched[target]` increments for every accepted/weak hit by query abundance (`--sizein`) or by 1.

Paired extension implementation:
- matched query if `report_hits` non-empty; unmatched otherwise
- `dbmatched_abundance[target]` increments for every accepted/weak hit with stock-like `--sizein` weighting
- top hit only drives OTU assignment
- paired output renderers:
  - `--userout`: fixed schema (`query`, `target`, `id_left`, `id_right`, `id_total`, `d_left`, `d_right`, `d_total`)
  - `--blast6out`, `--uc`: stock-shaped paired adaptations
  - `--alnout`: paired summary plus separate R1/R2 alignments
  - `--matched`/`--notmatched`: split paired FASTA only (`--matched` + `--match2`, `--notmatched` + `--notmatched2`)
- if unmatched and table output enabled:
  - default: unassigned row behavior (`nullptr` target)
  - with `--unknown_name`: assign unmatched pairs to that explicit target name

Extension mechanism: stock accounting rules are preserved; formatting is paired-specific.

### Step 9: End-of-run table/database outputs

Stock implementation: adds zero-count rows for DB targets not observed in matches, writes OTU/BIOM/mothur tables, and writes `dbmatched`/`dbnotmatched` FASTA outputs from accumulated counts.

Paired extension implementation:
- adds all paired DB targets to OTU table with zero-count completion
- optionally adds `--unknown_name` row
- writes OTU/BIOM/mothur tables
- writes paired DB outputs as split FASTA only:
  - matched targets use accumulated `dbmatched_abundance`
  - non-matched targets use original DB abundance

Extension mechanism: terminal aggregation semantics are stock-aligned with paired output encoding.

## Remaining paired-only deltas

- Paired `--userout` currently emits a fixed paired summary schema:
  - `query`, `target`, `id_left`, `id_right`, `id_total`, `d_left`, `d_right`, `d_total`
- Alignment-rich single-sequence outputs (`fastapairs`, `qsegout`, `tsegout`, `samout`, `lcaout`) are intentionally unsupported in paired mode.

## Example paired command

```bash
vsearch --usearch_global query_r1.fa query_r2.fa \
  --db tav_db_left.fa \
  --db2 tav_db_right.fa \
  --id 0.97 \
  --alnout paired.aln \
  --userout paired_userout.tsv \
  --uc paired.uc \
  --blast6out paired.b6 \
  --otutabout table.tsv
```

\newpage

## Source: `docs/parity/five_command_lowlevel_and_header_parity.md`

# Low-Level Runtime and Header Naming Parity (Five Extended Commands)

This note summarizes parity status for:

- `vsearch --fastq_filter`
- `vsearch --fastx_uniques`
- `vsearch --cluster_unoise`
- `vsearch --uchime3_denovo`
- `vsearch --usearch_global`

It focuses on:

1. low-level runtime behavior (`--threads`, SIMD Needleman-Wunsch, linear-memory fallback)
2. output FASTA/FASTQ header naming conventions
3. paired-orientation policy for R2

## Scope of concrete examples

Examples below use real IDs from `real_sub/r1_10k.fq.gz` and `real_sub/r2_10k.fq.gz` (first 1,000 pairs), with `--sample realsub --sizeout` to make header conventions visible.

Representative input IDs:

```text
@sample=CRLM1A_A1 1 494035 LH00659:244:232WKKLT3:4:1123:33117:29695
@sample=CRLM1A_A1 2 494035 LH00659:244:232WKKLT3:4:1123:33117:29695
```

## 1) Low-level runtime parity

| Command | Stock multicore | Ext multicore | Stock SIMD NW + LMA fallback | Ext SIMD NW + LMA fallback | Why in ext |
|---|---|---|---|---|---|
| `fastq_filter` | No (single stream) | No | N/A (no NW alignment path) | N/A | Stock mode (`--reverse`) remains on stock path; extension mode adds paired input/output CLI conventions and pair-level EE threshold evaluation. |
| `fastx_uniques` | No (single-thread derep path) | No | N/A (no NW alignment path) | N/A | Paired mode is custom (`tav_fastx_uniques`), but operation is counting/aggregation only. |
| `cluster_unoise` | Yes (`cluster.cc` thread workers) | Yes (`--threads` parallel delayed candidate alignments in `tav_cluster_unoise`) | Yes (`search16` SIMD first, LMA fallback on overflow) | Yes (`search16` SIMD first, LMA fallback on overflow) | Paired mode remains custom, but now has threaded candidate-alignment execution plus stock-like SIMD/LMA low-level alignment behavior. |
| `uchime3_denovo` | Yes (`chimera.cc` worker threads) | Yes (`--threads` parallel parent-candidate alignments in `tav_uchime3_denovo`) | Yes (candidate/full-query alignments use SIMD path then LMA fallback) | Yes (candidate/full-query paired-end alignments now use SIMD-first with stock overflow fallback) | Paired mode remains custom, but now has threaded parent-candidate alignment execution while preserving stock-style scoring/classification. |
| `usearch_global` | Yes (`search.cc` worker threads) | Yes (`--threads` parallel delayed candidate alignments in `tav_usearch_global`) | Yes (`searchcore.cc` uses `search16` with LMA fallback) | Yes (`search16` SIMD first, LMA fallback on overflow, per end) | Paired mode remains custom, but now has threaded candidate-alignment execution plus stock-like SIMD/LMA backend behavior. |

### Key detail: what is reused vs engineered

- Reused low-level naming/formatting:
  - `fasta_print_general`
  - `fastq_print_general`
  - `header_fprint_strip`
- Reused shared filter kernels across stock/ext:
  - `search_unaligned_numeric_filters_pass`
  - `search_aligned_compute_identity_metrics`
  - `search_aligned_threshold_filters_pass`
- Reused low-level alignment post-processing in ext:
  - `align_trim` is called by `align_one_end_stock_style`.
- Reused shared CIGAR traversal kernel:
  - `parse_cigar_string` (via `parse_cigar_operations_from_string` wrapper in paired `uchime3_denovo`)
- Reused shared fastq EE threshold kernel:
  - `ee_thresholds_pass` (stock per-read and extension pair-level aggregation both call it)
- Engineered in ext:
  - paired candidate enumeration, paired filtering, paired output writing for `cluster_unoise`, `uchime3_denovo`, `usearch_global`.
- Not currently reused by ext paired paths:
  - stock threaded execution engines in `cluster.cc`, `chimera.cc`, `search.cc`

## 2) Output header naming parity (with real IDs)

Important context:

- Most naming consistency comes from shared low-level printers (`fasta_print_general` / `fastq_print_general`), not from duplicated custom string logic.
- Because examples were run as a pipeline with `--sample realsub`, `;sample=realsub` can appear multiple times when a previous output header already contains `;sample=...` and the next command appends it again.

### A) `fastq_filter`

Status: stock paired behavior is preserved, and extension mode is intentionally additive.

```text
Stock output R1:
@sample=CRLM1A_A1 1 494035 LH00659:244:232WKKLT3:4:1123:33117:29695;sample=realsub;size=1

Stock output R2:
@sample=CRLM1A_A1 2 494035 LH00659:244:232WKKLT3:4:1123:33117:29695;sample=realsub;size=1
```

Why:

- Stock mode (`--reverse`) still reuses stock `filter.cc` output calls to `fastq_print_general`.
- Extension mode (`R1 R2` positional input or `--interleaved`) adds paired split outputs (`*2`) and uses pair-level `--fastq_maxee` / `--fastq_maxee_rate` checks.
- See `docs/parity/fastq_filter_paired_parity.md` for full mode/semantics details.

### B) `fastx_uniques`

Status: nearly same naming convention; paired extension writes two coordinated FASTA files.

```text
Stock single-end (`--fastaout`):
>sample=CRLM1A_A1;sample=realsub;size=803

Extension paired R1 (`--fastaout`):
>sample=CRLM1A_A1;sample=realsub;size=749

Extension paired R2 (`--fastaout_rev`):
>sample=CRLM1A_A1;sample=realsub;size=749
```

Why:

- Base IDs come from extension logic (paired keying and record tracking).
- Final header formatting (`;sample=...`, `;size=...`, relabel/strip rules) comes from shared `fasta_print_general`.

### C) `cluster_unoise`

Status: naming convention is mostly aligned; paired extension emits synchronized left/right centroid files.

```text
Stock centroid (`--centroids`):
>sample=CRLM1A_A1;sample=realsub;sample=realsub;size=33

Extension centroid R1 (`--fastaout`):
>sample=CRLM1A_A1;sample=realsub;sample=realsub;size=839

Extension centroid R2 (`--fastaout_rev`):
>sample=CRLM1A_A1;sample=realsub;sample=realsub;size=839
```

Why:

- Paired centroid identity bookkeeping is extension logic.
- Header suffix conventions (`;sample`, `;size`, relabel/strip) are still from shared FASTA printer.

### D) `uchime3_denovo`

Status: shared convention; paired split FASTA outputs now keep the same base header on both R1 and R2 files (no extra `/1` or `/2`).

```text
Stock nonchimera (`--nonchimeras`):
>sample=CRLM1A_A1;sample=realsub;sample=realsub;sample=realsub;size=33

Extension nonchimera pair split output (`--nonchimeras` + `--nonchimeras2`):
R1 file (`--nonchimeras`):
>sample=CRLM1A_A1;sample=realsub;sample=realsub;size=839;sample=realsub;size=839
R2 file (`--nonchimeras2`):
>sample=CRLM1A_A1;sample=realsub;sample=realsub;size=839;sample=realsub;size=839

Extension chimera pair split output (`--chimeras` + `--chimeras2`):
R1 file (`--chimeras`):
>sample=CRLR1A_A10;sample=realsub;sample=realsub;size=2;sample=realsub;size=2
R2 file (`--chimeras2`):
>sample=CRLR1A_A10;sample=realsub;sample=realsub;size=2;sample=realsub;size=2
```

Why:

- Paired split FASTA writers now emit the same base header in both files.
- Final header formatting still comes from shared `fasta_print_general`.

### E) `usearch_global`

Status: shared convention; paired split FASTA outputs now keep the same base header on both R1 and R2 files.

```text
Stock matched (`--matched`):
>sample=CRLM1A_A1;sample=realsub;sample=realsub;size=803

Stock notmatched (`--notmatched`):
>sample=CRLM1A_A1;sample=realsub;sample=realsub;size=1

Extension notmatched pair split output (`--notmatched` + `--notmatched2`):
R1 file (`--notmatched`):
>sample=CRLM1A_A1;sample=realsub;size=749;sample=realsub;size=749
R2 file (`--notmatched2`):
>sample=CRLM1A_A1;sample=realsub;size=749;sample=realsub;size=749

Extension matched pair split output example (`--matched` + `--match2`; from `--id 0.80` run):
R1 file (`--matched`):
>sample=CRLM1A_A1;sample=realsub;size=749;sample=realsub;size=749
R2 file (`--match2`):
>sample=CRLM1A_A1;sample=realsub;size=749;sample=realsub;size=749
```

Extension matched naming follows the same pattern as extension notmatched (same base header in R1/R2 split files).

Why:

- Base strip/relabel/sample/size behavior is inherited from shared FASTA output code.
- No extra `/1` or `/2` suffix is appended by paired output helpers.

## Bottom line

- Naming parity is generally strong because both stock and extension rely on the same low-level header-printing functions.
- Runtime backend parity is mixed:
  - `fastq_filter` has mixed parity: stock paired mode is unchanged, while extension mode intentionally changes CLI/output conventions and EE-threshold semantics,
  - mostly N/A for `fastx_uniques` (no alignment backend),
  - `cluster_unoise`, `uchime3_denovo`, and `usearch_global` now have SIMD-first + overflow-fallback parity in paired mode,
  - these three commands now also honor `--threads` for candidate-alignment-heavy stages, while still using custom paired orchestration code paths.

## Orientation policy status

- Extension pipeline policy is now unified: R2 remains in original read orientation at external I/O for paired
  `fastx_uniques`, `cluster_unoise`, `uchime3_denovo`, and `usearch_global`.
- Previous behavior that persisted R2 as reverse-complemented is removed.
- Internal orientation handling may still occur where algorithmically required (for example RDP classifier internals),
  but this does not change external orientation contracts.
- See `docs/parity/paired_orientation_unification.md` for developer details.

\newpage

## Source: `docs/parity/paired_orientation_unification.md`

# Paired Orientation Unification (R2 Kept Native)

This note documents the paired-orientation policy for the extension pipeline.

## Policy

- External input/output orientation is now unified: R1 and R2 stay in their original read orientation.
- We no longer reverse-complement R2 as a persistent data transformation in the paired VSEARCH extension path.
- Internal reverse handling is allowed only when required by a specific algorithm.

## Command-level behavior

For the paired extension path (`fastx_uniques -> cluster_unoise -> uchime3_denovo -> usearch_global`) and paired RDP classification:

1. `fastx_uniques` (paired extension)
- R2 is not reverse-complemented at ingestion.
- Right anchors are written in native R2 orientation.

2. `cluster_unoise` (paired extension)
- Uses paired records from `fastx_uniques` in native orientation.
- No persistent R2 RC transform.

3. `uchime3_denovo` (paired extension)
- Uses paired records in native orientation.
- Breakpoint scan logic is orientation-aware for the right segment:
  - left segment scanned from start to end
  - right segment scanned from end to start
  This preserves the intended combined-axis chimera model while keeping stored R2 unchanged.

4. `usearch_global` (paired extension)
- Query and database right-end anchors are compared in native R2 orientation.
- No query-side forced R2 RC at ingestion.

5. RDP paired classifier (Java extension)
- Query-side orientation normalization remains internal to stock-RDP logic:
  - `trainingInfo.isSeqReversed(...)`
  - `ClassifierSequence.getReversedSeq()`
- This is an internal scoring detail; it does not change external FASTA/FASTQ orientation contracts.

## Developer guidance

- If you add new paired steps, keep persisted/output R2 in native orientation unless there is a hard compatibility reason not to.
- If a method needs orientation normalization for scoring, do it inside that method and keep it local to computation.
- Document any internal RC usage explicitly in parity notes.

\newpage

## Source: `docs/parity/rdp_tav_taxonomy_paired_parity.md`

# RDP Taxonomy Paired-End Parity (TAV)

## Scope

This note compares stock single-sequence RDP classifier behavior with the native paired-end TAV extension.

Goal: preserve the stock Naive Bayes backbone while changing the query object from ASV to TAV `(R1,R2)`.

## Key decisions

1. No catalog mode in extension command surface.
2. No posthoc per-anchor merge policy.
3. Paired support starts inside NB likelihood and bootstrap computation.
4. One primary TAV output stream, stock formatter semantics.

## Why Python launcher + Java core

- The classifier core is Java (`PairedClassifierMain`, `PairedNaiveBayesClassifier`) because it directly reuses stock RDP Java internals.
- The Python entrypoint is only a thin launcher for:
  - resolving RDP jar/model paths from `manifest.json`
  - compiling local extension Java sources
  - constructing Java classpath/JVM options consistently
- Classification math and formatting still run in Java.

## Command surface

```bash
python3 rdp_tav_taxonomy.py \
  --input tav_left.fa \
  --input2 tav_right.fa \
  --output tav_taxonomy.tsv
```

Interleaved mode:

```bash
python3 rdp_tav_taxonomy.py \
  --input tav_interleaved.fa \
  --interleaved \
  --output tav_taxonomy.tsv
```

Stock-classify option aliases are supported for core settings:

- `-o/--output/--outputFile`
- `-f/--format`
- `-c/--conf`
- `-w/--min-words/--minWords`
- `-s/--shortseq-outfile/--shortseq_outfile`
- `-t/--train-prop/--train_propfile`
- `-g/--gene`
- `-q/--queryFile` (accepted and ignored, matching stock legacy behavior)
- `-b/--bootstrap_outfile` (supported; stock-style bootstrap summary output)
- `-h/--hier_outfile` (supported; stock-style hierarchy count output)

Stock classify metadata/biom options currently parsed for parity but intentionally unsupported in paired TAV mode:

- `-d/--metadata`
- `-m/--biomFile`
- `--format biom`

## Backbone parity mapping

### 1. Training/model loading parity: preserved

- Uses stock `ClassifierFactory` and stock training property files.
- Reuses pretrained `rRNAClassifier.properties` assets downloaded from SourceForge.

### 2. Feature model parity: preserved with paired generalization

- Stock: overlapping k-mers from one sequence.
- Paired extension: overlapping k-mers from both anchors, treated as one combined evidence pool.
- Orientation normalization remains internal to stock-RDP mechanics (`isSeqReversed`/`getReversedSeq`) and does not change external paired FASTA orientation contracts.

### 3. Likelihood/posterior parity: preserved with paired generalization

- Stock posterior uses sum of per-word log-probabilities for one sequence.
- Paired posterior uses the same equation over concatenated word rows from R1 and R2.
- Operationally: `log P(TAV|g) = log P(R1|g) + log P(R2|g)`.

### 4. Bootstrap parity: preserved with paired generalization

- Stock: bootstrap draws words from one query word pool.
- Paired extension: bootstrap draws from combined R1+R2 word pool.
- This is native evidence pooling, not posthoc confidence merging.

### 5. Output formatter parity: preserved

- Uses stock `ClassificationResultFormatter` (`allrank`, `fixrank`, `filterbyconf`, `db`).
- Produces one line per TAV ID.

## Real data example

From `data/real_out2/tav_taxonomy_native.tsv`:

- `TAV000001;size=8242` -> `... Micrococcaceae family 1.0 ... Pseudarthrobacter genus 0.98`
- `TAV000002;size=310` -> `... Pseudomonadaceae family 1.0 ... Pseudomonas genus 0.99`
- `TAV000003;size=25` -> `... Oxalobacteraceae family 0.98 ... Collimonas genus 0.89`

These are single TAV-level calls generated by one paired NB path.

## Intentional disparities

- New paired IO mode (`--input/--input2` or `--interleaved`).
- New paired classifier class that extends stock NB logic to paired query units.
- Stock classify metadata/biom flows (`--metadata`, `--biomFile`, `--format biom`) are not implemented in paired TAV mode yet.

## Conclusion

Parity is preserved for the stock RDP model, scoring, and formatter backbone. The extension is a native paired-query NB generalization, not a posthoc two-run merge.
