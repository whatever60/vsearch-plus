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

Low-level relationship in current code:

- Paired Step 2 call stack:
  ```text
  search_onequery_paired(...)
    -> unique_count(...) + paired postings accumulation
  search_joinhits_paired(...)
    -> minheap_add/sort/pop
  ```

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
