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

### Paired mode: `cluster_unoise` + `--reverse`

- Input: two synchronized sequence files:
  - `--cluster_unoise` for left reads
  - `--reverse` for right reads
- Catalog mode input is intentionally removed for `cluster_unoise` parity.
- Length filtering honors `--minseqlength`/`--maxseqlength` in paired mode with:
  - `--filter any` (default): discard pair if either end fails threshold
  - `--filter both`: discard pair only if both ends fail threshold
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

Stock pre-alignment filters in search_acceptable_unaligned are:

maxqsize: qsize <= maxqsize
mintsize: tsize >= mintsize
minsizeratio: qsize >= minsizeratio * tsize
maxsizeratio: qsize <= maxsizeratio * tsize
minqt: qlen >= minqt * tlen
maxqt: qlen <= maxqt * tlen
minsl: shorter_len >= minsl * longer_len
maxsl: shorter_len <= maxsl * longer_len
idprefix: first N bases must match exactly
idsuffix: last N bases must match exactly
self: reject same header when self enabled
selfid: reject perfect self sequence when selfid enabled
Paired extension applies paired analogs of the same family:

maxqsize/mintsize/size ratios on paired abundances
minqt/maxqt/minsl/maxsl on combined paired lengths (left+right)
idprefix/idsuffix required on both ends
self/selfid checks on paired header/paired sequence identity
minseqlength/maxseqlength paired drop policy controlled by `--filter any|both`
Extension mechanism: same filter categories, generalized to two ends.

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

weak_id threshold: id >= 100 * weak_id
maxsubs: mismatches <= maxsubs
maxgaps: internal_gaps <= maxgaps
mincols: internal_alignmentlength >= mincols
leftjust (if enabled): no left terminal gaps
rightjust (if enabled): no right terminal gaps
query_cov: (matches + mismatches) >= query_cov * qlen
target_cov: (matches + mismatches) >= target_cov * tlen
maxid: id <= 100 * maxid
mid: 100 * matches / (matches + mismatches) >= mid
maxdiffs: mismatches + internal_indels <= maxdiffs
Paired extension uses paired analog checks over aggregated two-end stats:

weak_id on combined identity
maxsubs on total mismatches (left+right)
maxgaps on total internal gaps (left+right)
mincols on aggregated internal alignment columns
query_cov/target_cov on paired totals
maxid and mid on paired-combined metrics
maxdiffs on paired mismatch+internal-indel totals
Extension mechanism: preserve stock filter intent by replacing single-end counts with summed paired counts.
leftjust/rightjust are mirrored with paired terminal-gap aggregation.

Terminal gap trims are the leading/trailing gap runs in each end-alignment (`trim_q_left`, `trim_t_left`, `trim_q_right`, `trim_t_right`) produced by the stock `align_trim` logic; leftjust/rightjust require these summed terminal trims to be zero on the respective side.

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
  - paired min/max length policy with `--filter any|both`
  - paired output requirements (`--fastaout_rev` companion output)


## Example paired command

```bash
./bin/vsearch \
  --cluster_unoise uniques_r1.fasta \
  --reverse uniques_r2.fasta \
  --centroids denoised_r1.fasta \
  --fastaout_rev denoised_r2.fasta \
  --tabbedout denoised_pairs.tsv
```
