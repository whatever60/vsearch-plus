# Native Paired `cluster_unoise`: SOP and Implementation Checklist

This document captures the working SOP and the concrete implementation checklist for extending `vsearch --cluster_unoise` to support paired-end reads natively.

It is based on:

- `docs/parity/cluster_unoise_paired_parity.md`
- `docs/callstacks/cluster_unoise.txt`
- `src/cluster.cc`
- `src/searchcore.cc`

## Current Status

- The native paired implementation lives in:
  - `src/cluster_paired.cc`
  - `src/searchcore_paired.cc`
  - `src/dbindex_paired.cc`
- The full project currently builds with `make -C src -j8`.
- The paired native path is past the bootstrap stage:
  - stock-shaped paired search and cluster modules are present
  - paired dbindex support is present
  - the paired native path is linked into the normal `vsearch` binary
- The paired core engine parity checklist is now closed for the current supported paired command surface:
  - `searchinfo_s_paired` uses stock-style raw query/k-mer/hit storage
  - `hit_paired_s` owns raw per-end stock `struct hit` alignment state
  - `clusterinfo_s_paired` retains per-end CIGAR ownership in stock style
  - both `--threads 2` and `--threads 1` paired smoke outputs now agree on the hit-case and OTU-table smoke inputs
- A native paired smoke test has succeeded with:
  - `./bin/vsearch --cluster_unoise .tmp_refactor_smoke/r1.fa .tmp_refactor_smoke/r2.fa --minsize 1 --fastaout .tmp_native_cluster_smoke/out_r1.fa --fastaout_rev .tmp_native_cluster_smoke/out_r2.fa --threads 2`
- A paired core-result smoke test has also succeeded with:
  - `--uc`
  - `--otutabout`
  - a tiny hit-case input that produced real `H` lines in UC output
- A serial-only paired dbindex insertion bug was fixed during the parity closure pass:
  - `dbindex_addsequence_paired()` now consumes the R1 unique-kmer list before the shared unique-count handle is reused for R2
  - this removed a real serial/parallel behavior mismatch where same-round parallel reconciliation had been masking serial misses

## Purpose

The goal is not to design a brand new paired denoiser. The goal is to reproduce the stock `cluster_unoise` structure and behavior as faithfully as possible, while replacing one read with one synchronized pair of reads.

The guiding rule is:

```text
copy stock first
rename with _paired
replace one read with two reads
reuse stock primitives
justify every divergence
```

## SOP

### 1. Define paired semantics before touching the code

Write down how each stock single-read step becomes paired-read logic in the parity docs before changing implementation details.

For `cluster_unoise`, the main reference is:

- `docs/parity/cluster_unoise_paired_parity.md`

Cross-check against other paired command docs for consistency:

- `docs/parity/fastq_filter_paired_parity.md`
- `docs/parity/usearch_global_paired_parity.md`
- `docs/parity/uchime3_denovo_paired_parity.md`

This prevents code from drifting into undocumented behavior.

### 2. Freeze the complete stock callstack first

Treat stock `cluster_unoise` as the ground truth and document the complete callstack before implementing the paired port.

The stock callstack reference is:

- `docs/callstacks/cluster_unoise.txt`

We care about both execution modes:

- `--threads > 1`
- `--threads = 1`

The target is not vague behavioral similarity. The target is one-to-one structural parity wherever paired-end semantics do not force a difference.

### 3. Port by stock module boundary

Duplicate stock source files that are on the `cluster_unoise` callstack and need paired-aware state.

Current module mapping:

- `src/cluster.cc` -> `src/cluster_paired.cc`
- `src/searchcore.cc` -> `src/searchcore_paired.cc`
- stock dbindex layer -> `src/dbindex_paired.cc`

If more stock modules become necessary for parity, duplicate those too as `*_paired.cc`.

Do not keep expanding `src/tav_extension.cc` as the final native implementation path. That file is useful as a semantic reference, but the native paired command should live in stock-shaped paired modules.

### 4. Keep paired names and function shapes as close to stock as possible

The paired port should look like the stock code with `_paired` appended, not like a separate redesign.

Examples:

- `searchinfo_s` -> `searchinfo_s_paired`
- `search_onequery` -> `search_onequery_paired`
- `search_topscores` -> `search_topscores_paired`
- `cluster_core_parallel` -> `cluster_core_parallel_paired`
- `cluster_core_serial` -> `cluster_core_serial_paired`

The signature and implementation body should stay almost exactly the same unless the difference is explicitly about one end versus two ends.

### 5. Replace one read with one pair of reads systematically

This is the central transformation rule.

Examples:

- `qsequence` -> `qsequence_r1` and `qsequence_r2`
- `qseqlen` -> `qseqlen_r1` and `qseqlen_r2`
- one query k-mer sample -> one R1 sample plus one R2 sample
- one target sequence identifier -> one paired target identifier plus mapping to R1 and R2 sequence numbers
- one alignment stat block -> per-end stats plus pair-level aggregates only where stock logic expects a single value

The stock control flow should stay intact. Only the data flowing through it becomes paired.

### 6. Reuse stock low-level primitives whenever possible

The preferred approach is to reuse stock building blocks instead of re-implementing them.

Important shared primitives include:

- `unique_init`
- `unique_exit`
- `unique_count`
- `unique_count_shared`
- `search16_init`
- `search16_qprep`
- `search16`
- `search16_exit`
- `LinearMemoryAligner`
- `align_trim`
- `search_unaligned_numeric_filters_pass`
- `search_aligned_compute_identity_metrics`
- `search_aligned_threshold_filters_pass`
- `db_getsequence`
- `db_getsequencelen`
- `db_getabundance`
- `db_getheader`

If stock `db_*` functions work with stock sequence numbers, the paired port should translate from logical paired targets to the underlying R1 and R2 sequence numbers and then reuse those stock functions.

The paired dbindex layer should follow the same rule:

- keep a stock-shaped API
- keep paired-only logic at the mapping/keying layer
- verify both serial and parallel paths against the same smoke inputs, because same-round parallel reconciliation can hide dbindex insertion bugs that the serial path exposes immediately

### 7. Avoid non-stock helpers unless they are truly necessary

If stock `searchcore.cc` does not define a helper, the default should be that `searchcore_paired.cc` also does not define one.

When possible:

- inline the logic back into the stock-shaped function
- or move the implementation closer to the stock structure

The paired port should be explainable as a stock port, not as a new helper framework.

### 8. Reproduce both stock execution modes

The native paired implementation must mirror the stock structure for:

- `--threads > 1`
- `--threads = 1`

This means the paired port should preserve the stock distinction between:

- parallel worker execution
- main-thread reconciliation
- best-hit selection
- result handling
- tail/output materialization

This matters because stock `cluster_unoise` is not just serial logic wrapped in threads. It has a specific scheduling and reconciliation pattern that needs to be mirrored.

### 9. Audit every paired function line by line

For each paired function and struct:

1. identify the exact stock counterpart
2. compare signature shape
3. compare body structure
4. justify every extra line by paired semantics
5. verify stock primitives are reused where possible
6. remove paired-only helpers that exist only for convenience

If a difference cannot be explained by replacing one read with one pair of reads, it is probably a parity bug.

### 10. Only optimize after parity is structurally correct

Performance work should come after semantic and structural parity.

An optimization is only a parity fix when it corrects a real stock/ext structural mismatch, such as:

- failing to use the stock threading model
- failing to reuse the stock SIMD path
- introducing a serial bottleneck where stock is parallel

The preferred order is:

1. define semantics
2. lock stock callstack
3. implement one-to-one paired port
4. audit divergences
5. benchmark on full data

## Target Callstack

### `--threads > 1`

```text
cmd_cluster
-> cluster_unoise_paired
-> cluster_paired
-> db setup / sorting / dbindex_prepare_paired
-> cluster_core_parallel_paired
-> threads_init_paired
-> cluster_query_init_paired
-> threads_wakeup_paired
-> cluster_worker_paired
-> cluster_query_core_paired
-> search_onequery_paired
-> search_topscores_paired
-> search_acceptable_unaligned_paired
-> align_delayed_paired
-> search_acceptable_aligned_paired
-> main-thread extra_list reconciliation
-> search_findbest2_byid_paired or search_findbest2_bysize_paired
-> cluster_core_results_hit/nohit_paired
-> cluster tail/output machinery
-> dbindex_free_paired
-> db_free
```

### `--threads = 1`

```text
cmd_cluster
-> cluster_unoise_paired
-> cluster_paired
-> db setup / sorting / dbindex_prepare_paired
-> cluster_core_serial_paired
-> cluster_query_init_paired
-> cluster_query_core_paired
-> search_onequery_paired
-> search_topscores_paired
-> search_acceptable_unaligned_paired
-> align_delayed_paired
-> search_acceptable_aligned_paired
-> search_findbest2_byid_paired or search_findbest2_bysize_paired
-> cluster_core_results_hit/nohit_paired
-> cluster tail/output machinery
-> dbindex_free_paired
-> db_free
```

## Implementation Checklist

### 1. Freeze the target behavior

- Treat `src/cluster.cc` and `src/searchcore.cc` as the only behavioral ground truth.
- Treat `docs/parity/cluster_unoise_paired_parity.md` as the paired semantic contract.
- Require every paired-only difference to be justified as a one-read to one-pair transformation.

### 2. Keep the file split aligned to stock modules

- `src/cluster.cc` -> `src/cluster_paired.cc`
- `src/searchcore.cc` -> `src/searchcore_paired.cc`
- stock dbindex layer -> `src/dbindex_paired.cc`

Do not continue growing `src/tav_extension.cc` for the final native paired `cluster_unoise` path.

### 3. Struct mapping checklist

Target one-to-one struct mapping:

- `searchinfo_s` -> `searchinfo_s_paired`
- `hit_s` -> `hit_paired_s`
- `clusterinfo_s` -> `clusterinfo_s_paired`
- stock db/index state -> paired dbindex state
- single sequence record -> `record_paired_s`

Required paired substitutions:

- `qsequence` -> `qsequence_r1` and `qsequence_r2`
- `qseqlen` -> `qseqlen_r1` and `qseqlen_r2`
- `kmersample` -> `kmersample_r1` and `kmersample_r2`
- `uh` -> `uh_r1` and `uh_r2`
- `s16info` -> `s_r1` and `s_r2`
- `lma` -> `lma_r1` and `lma_r2`
- target seqno -> logical paired target plus `target_seqnos_r1` and `target_seqnos_r2`

### 4. Function mapping checklist

Cluster layer:

- `cluster_unoise` -> `cluster_unoise_paired`
- `cluster` -> `cluster_paired`
- `cluster_core_parallel` -> `cluster_core_parallel_paired`
- `cluster_core_serial` -> `cluster_core_serial_paired`
- `cluster_query_init` -> `cluster_query_init_paired`
- `cluster_query_exit` -> `cluster_query_exit_paired`
- `cluster_query_core` -> `cluster_query_core_paired`
- `cluster_worker` -> `cluster_worker_paired`
- `threads_init` -> `threads_init_paired`
- `threads_wakeup` -> `threads_wakeup_paired`
- `threads_exit` -> `threads_exit_paired`

Search core:

- `search_topscores` -> `search_topscores_paired`
- `search_onequery` -> `search_onequery_paired`
- `search_acceptable_unaligned` -> `search_acceptable_unaligned_paired`
- `search_acceptable_aligned` -> `search_acceptable_aligned_paired`
- `search_findbest2_byid` -> `search_findbest2_byid_paired`
- `search_findbest2_bysize` -> `search_findbest2_bysize_paired`
- `search_joinhits` -> `search_joinhits_paired`
- `search_enough_kmers` -> `search_enough_kmers_paired`
- `align_delayed` -> `align_delayed_paired`

Paired dbindex layer:

- `dbindex_prepare` -> `dbindex_prepare_paired`
- `dbindex_addsequence` -> `dbindex_addsequence_paired`
- `dbindex_addallsequences` -> `dbindex_addallsequences_paired`
- `dbindex_getcount` -> `dbindex_getcount_paired`
- `dbindex_getmapping` -> `dbindex_getmapping_paired`
- `dbindex_getmatchcount` -> `dbindex_getmatchcount_paired`
- `dbindex_getmatchlist` -> `dbindex_getmatchlist_paired`
- `dbindex_getbitmap` -> `dbindex_getbitmap_paired`
- `dbindex_free` -> `dbindex_free_paired`

### 5. Shared low-level primitive reuse checklist

These should remain shared rather than re-implemented:

- `unique_init`
- `unique_exit`
- `unique_count`
- `unique_count_shared`
- `search16_init`
- `search16_qprep`
- `search16`
- `search16_exit`
- `LinearMemoryAligner`
- `align_trim`
- `search_unaligned_numeric_filters_pass`
- `search_aligned_compute_identity_metrics`
- `search_aligned_threshold_filters_pass`
- `db_getsequence`
- `db_getsequencelen`
- `db_getabundance`
- `db_getheader`
- `reverse_complement`
- `otutable_*`
- `fasta_print_general`
- `progress_*`

### 6. Paired-only logic points checklist

These are the places where divergence from stock is expected:

- `cluster_query_core_paired`
  - load two ends instead of one
  - reverse-complement both ends appropriately for reverse-strand search

- `search_topscores_paired`
  - accumulate k-mer evidence from both ends
  - keep stock heap and comparator structure
  - use paired dbindex with stock-like control flow

- `search_acceptable_unaligned_paired`
  - evaluate stock numeric filters on paired aggregate length and abundance metrics
  - apply paired prefix, suffix, self, and selfid rules

- `align_delayed_paired`
  - align R1 and R2 separately using stock SIMD and LMA paths
  - aggregate the per-end stats afterward

- `search_acceptable_aligned_paired`
  - compute stock aligned thresholds on aggregated totals

- `search_findbest2_*_paired`
  - preserve stock ranking behavior
  - resolve abundance and length through stock `db_*` using the R1 and R2 mapping

- cluster tail and output code
  - paired centroid outputs
  - paired UC, clusters, tabbed, and OTU-table behavior

### 7. Line-by-line review checklist

For each `_paired` function:

1. identify the exact stock counterpart
2. compare the signature
3. compare the body structure
4. justify every extra line by paired semantics
5. verify stock primitive reuse
6. remove or refold helpers that stock does not have

If a difference cannot be explained by replacing one read with one synchronized read pair, treat it as a likely parity issue.

### 8. Recommended implementation order

1. lock semantics in `docs/parity/cluster_unoise_paired_parity.md`
2. lock the stock callstack in `docs/callstacks/cluster_unoise.txt`
3. build or clean `searchinfo_s_paired`, `hit_paired_s`, and `clusterinfo_s_paired`
4. build or clean `dbindex_paired` so `search_topscores_paired` can stay stock-shaped
5. port `src/searchcore.cc` to `src/searchcore_paired.cc` line by line
6. port `src/cluster.cc` to `src/cluster_paired.cc` line by line
7. port tail and output machinery last
8. benchmark only after structural parity is in place

### 9. Definition of done

The native paired `cluster_unoise` port is only done when all of the following are true:

- the `--threads > 1` callstack matches stock structure one to one
- the `--threads = 1` callstack matches stock structure one to one
- `src/searchcore_paired.cc` is stock-shaped rather than extension-shaped
- `src/cluster_paired.cc` is stock-shaped rather than a shim into `tav_extension`
- stock `db_*` accessors are reused whenever possible
- paired `dbindex_*` closely mirrors stock `dbindex_*`
- all paired-only differences are documented in the parity markdown
- full-data behavior and timing can be compared honestly against stock

## Short Version

The workflow is:

1. define paired semantics in markdown first
2. document the complete stock callstack
3. duplicate stock modules into `*_paired.cc`
4. rename stock structs and functions with `_paired`
5. replace one read with two reads, systematically
6. reuse stock primitives and stock `db_*` behavior
7. keep both serial and parallel stock execution structure
8. audit every divergence line by line
9. benchmark only after parity is structurally correct
