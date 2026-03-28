# Native Paired `usearch_global`: SOP and Implementation Checklist

This document captures the working SOP and the concrete implementation checklist for extending `vsearch --usearch_global` to support paired-end reads natively in the same stock-shaped style we used for paired `cluster_unoise`.

It is based on:

- `docs/parity/usearch_global_paired_parity.md`
- `docs/callstacks/usearch_global.txt`
- `src/search.cc`
- `src/searchcore.cc`

## Current Status

- Paired CLI routing in `cmd_usearch_global()` now splits cleanly:
  - paired input is routed to `usearch_global_paired()`
  - single-end input is routed to stock `usearch_global()`
- A stock-shaped native paired wrapper now exists in `src/search_paired.cc`.
  - the long-term paired engine is no longer `tav_usearch_global()`
  - the wrapper now follows the stock `src/search.cc` layout:
    - `search_output_results_paired`
    - `search_query_paired`
    - `search_thread_run_paired`
    - `search_thread_init_paired`
    - `search_thread_exit_paired`
    - `search_thread_worker_paired`
    - `search_thread_worker_run_paired`
    - `search_prep_paired`
    - `search_done_paired`
    - `usearch_global_paired`
- The shared paired search machinery is now reused directly:
  - `src/searchcore_paired.cc`
    - `search_topscores_paired`
    - `search_acceptable_unaligned_paired`
    - `search_acceptable_aligned_paired`
    - `align_delayed_paired`
    - `search_onequery_paired`
    - `search_joinhits_paired`
  - `src/dbindex_paired.cc`
    - `dbindex_prepare_paired`
    - `dbindex_addallsequences_paired`
    - `dbindex_addsequence_paired`
    - paired bitmap/matchlist/mapping accessors
- The paired DB/query loaders now live under `search_prep_paired()` / `search_thread_run_paired()`:
  - split query input
  - interleaved query input
  - split paired DB input
  - interleaved paired DB input
  - split paired FASTA outputs
  - paired UC / userout / blast6 / alnout / OTU-table outputs
- The serial-only paired dbindex insertion bug we fixed during paired `cluster_unoise` parity closure also matters here:
  - paired dbindex behavior must always be checked in both `--threads 1` and `--threads > 1`
  - same-round parallel behavior can hide insertion bugs that the one-worker path exposes immediately
- Verification completed for the native paired `usearch_global` path:
  - full `make -C src -j8`
  - split-input smoke test with `--threads 1` and `--threads 2`
  - exact output agreement for `uc`, `otutabout`, `userout`, `blast6out`, `matched/match2`, `notmatched/notmatched2`, `dbmatched/dbmatched2`, `dbnotmatched/dbnotmatched2`
  - `alnout` body agreement between `--threads 1` and `--threads 2` (command header line differs because output filenames and thread count differ)
  - interleaved-input smoke test matching the split-input OTU table
- One convenience helper was refolded to keep the wrapper closer to stock shape:
  - `mask_sequence_paired` was folded into the DB-load block inside `search_prep_paired()`
- Paired loaders now preserve separate R1 and R2 headers while enforcing that the first whitespace-delimited token matches between ends.
  - this check is independent of `--notrunclabels`
  - there is no extra `/1` or `/2` normalization layer
- The remaining non-stock helper functions are narrow paired-only helpers:
  - `paired_header_key_paired`
  - `get_anchor_len_paired`
  - `is_catalog_file_paired`
  - `write_fasta_record_paired`
  - `write_query_pair_split_paired`
  - `write_db_pair_split_paired`

## Purpose

The goal is not to keep growing the extension path in the now-removed legacy `tav_extension.cc` path.

The goal is to reproduce stock `usearch_global` structure and behavior as faithfully as possible, while replacing one query read and one target read with one synchronized query pair and one synchronized target pair.

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

Write down how each stock single-read `usearch_global` step becomes paired-read logic before moving code.

Main semantic reference:

- `docs/parity/usearch_global_paired_parity.md`

Cross-check with the other paired docs so the command family stays internally consistent:

- `docs/parity/fastq_filter_paired_parity.md`
- `docs/parity/cluster_unoise_paired_parity.md`
- `docs/parity/uchime3_denovo_paired_parity.md`

### 2. Freeze the complete stock callstack first

Treat stock `usearch_global` as the only structural ground truth and document the complete callstack before porting.

Stock callstack reference:

- `docs/callstacks/usearch_global.txt`

Important nuance relative to `cluster_unoise`:

- stock `search.cc` has no separate serial core function
- `--threads 1` still uses the same threaded driver, just with one worker thread

So the paired native port should mirror that fact instead of inventing a separate serial architecture.

### 3. Port by stock module boundary

Duplicate stock modules that sit on the `usearch_global` callstack and need paired-aware state.

Target module mapping:

- `src/search.cc` -> `src/search_paired.cc`
- `src/searchcore.cc` -> `src/searchcore_paired.cc`
- `src/dbindex.cc`/`src/dbindex.h` surface -> `src/dbindex_paired.cc`/`src/dbindex_paired.h`

Do not treat the now-removed legacy `tav_extension.cc` path as the final native location for paired `usearch_global`.
It was a semantic reference and temporary implementation source, not the end-state architecture.

### 4. Keep paired names and function shapes close to stock

The paired port should look like the stock code with `_paired` appended.

Examples for `usearch_global`:

- `search_output_results` -> `search_output_results_paired`
- `search_query` -> `search_query_paired`
- `search_thread_run` -> `search_thread_run_paired`
- `search_thread_init` -> `search_thread_init_paired`
- `search_thread_exit` -> `search_thread_exit_paired`
- `search_thread_worker` -> `search_thread_worker_paired`
- `search_thread_worker_run` -> `search_thread_worker_run_paired`
- `search_prep` -> `search_prep_paired`
- `search_done` -> `search_done_paired`
- `usearch_global` -> `usearch_global_paired`

The body should stay almost exactly the same unless the difference is explicitly caused by replacing one read with one pair.

### 5. Replace one query read with one query pair systematically

This is the main transformation rule.

Examples:

- one query FASTX record -> one synchronized `(R1, R2)` query unit
- one DB sequence -> one paired DB target unit
- one target seqno -> one paired target index plus mapping to R1/R2 seqnos
- one aligned hit -> one paired hit with two embedded stock per-end alignments plus aggregate ranking/filter state

The stock control flow should stay intact. Only the data model becomes paired.

### 6. Reuse the shared paired search wheels that already exist

For `usearch_global`, the most important lesson from paired `cluster_unoise` is that we should not rebuild the search core again.

Shared paired search pieces already completed:

- `search_topscores_paired`
- `search_acceptable_unaligned_paired`
- `search_acceptable_aligned_paired`
- `align_delayed_paired`
- `search_onequery_paired`
- `search_joinhits_paired`
- paired `hit_paired_s`
- paired `searchinfo_s_paired`
- paired `dbindex_*` API

So the native `usearch_global` port should concentrate on the stock `search.cc` wrapper/orchestration layer:

- paired query reading
- stock-shaped worker orchestration
- stock-shaped per-query output emission
- stock-shaped summary / dbmatched / OTU-table tail

### 7. Reuse stock low-level primitives aggressively

The paired port should continue reusing stock building blocks:

- `unique_init`, `unique_exit`, `unique_count`
- `search16_init`, `search16_qprep`, `search16`, `search16_exit`
- `LinearMemoryAligner`
- `align_trim`
- `search_unaligned_numeric_filters_pass`
- `search_aligned_compute_identity_metrics`
- `search_aligned_threshold_filters_pass`
- `db_getsequence`, `db_getsequencelen`, `db_getabundance`, `db_getheader`
- OTU-table and results helpers where their output surface is still supported

If a stock helper already does the right thing once given stock-like per-end state, reuse it.

### 8. Avoid non-stock helpers unless they are truly necessary

The default should be to delete/refold usearch-global-specific helpers from the now-removed legacy `tav_extension.cc` path when the native port lands.

Examples that should usually disappear into stock-shaped paired code rather than survive as long-term helpers:

- custom paired candidate sorting wrappers
- custom paired aligned-filter wrappers when `searchcore_paired.cc` already covers that behavior
- one-off load/write helpers that only exist to support the extension path

### 9. Mirror stock thread behavior exactly

For `usearch_global`, stock behavior matters here:

- there is one shared query reader guarded by `mutex_input`
- output/stat/progress updates are guarded by `mutex_output`
- each worker owns reusable per-thread search state
- `--threads 1` uses the same mechanism, just with one worker

The paired native port should preserve that architecture exactly.

### 10. Audit every paired function line by line

For each paired function and struct:

1. identify the exact stock counterpart
2. compare signature shape
3. compare body structure
4. justify every extra line by paired semantics
5. verify stock/shared paired primitives are reused where possible
6. delete helper layers that only exist because the current code came from the now-removed legacy `tav_extension.cc` path

### 11. Only optimize after structural parity is correct

Performance work comes after semantic and structural parity.

In particular, if we see a serial/parallel mismatch or a speed anomaly, we should first assume a structural parity bug in:

- paired dbindex insertion
- worker orchestration
- repeated per-query setup/teardown
- unnecessary divergence from stock hit/output flow

## Target Callstack

### `--threads > 1`

```text
cmd_usearch_global
-> usearch_global_paired
-> search_prep_paired
-> open dbmatched/dbnotmatched outputs
-> otutable_init
-> open paired query input
-> allocate si_plus_paired / si_minus_paired / pthread
-> init mutex_input / mutex_output
-> progress_init
-> search_thread_worker_run_paired
   -> search_thread_init_paired
   -> xpthread_create(..., search_thread_worker_paired)
      -> search_thread_run_paired
         -> read one paired query unit under mutex_input
         -> populate plus/minus searchinfo_s_paired state
         -> search_query_paired
            -> search_onequery_paired (plus[/minus])
            -> search_joinhits_paired
            -> search_output_results_paired
         -> free kept-hit alignments not retained by outputs/state
         -> update stats/progress under mutex_output
   -> xpthread_join
   -> search_thread_exit_paired
-> progress_done
-> print/log match summaries
-> zero-complete OTU table rows
-> write dbmatched/dbnotmatched paired outputs
-> search_done_paired
   -> dbindex_free_paired
   -> db_free
   -> close remaining outputs
```

### `--threads = 1`

```text
Same callstack as above, with one worker thread.
```

There should be no separate paired serial search architecture if we want true stock parity.

## Implementation Checklist

### 1. Freeze the target behavior

- Treat `src/search.cc` and `src/searchcore.cc` as the only behavioral ground truth for stock `usearch_global`.
- Treat `docs/parity/usearch_global_paired_parity.md` as the paired semantic contract.
- Treat `docs/callstacks/usearch_global.txt` as the structural target.

### 2. Build the stock-shaped paired wrapper layer in `src/search_paired.cc`

Implement paired analogs of the stock `search.cc` driver functions:

- `search_output_results_paired`
- `search_query_paired`
- `search_thread_run_paired`
- `search_thread_init_paired`
- `search_thread_exit_paired`
- `search_thread_worker_paired`
- `search_thread_worker_run_paired`
- `search_prep_paired`
- `search_done_paired`
- `usearch_global_paired`

### 3. Reuse `searchcore_paired.cc` instead of rebuilding search logic

Native paired `usearch_global` should call:

- `search_onequery_paired`
- `search_joinhits_paired`

It should not keep using the older legacy `tav_extension.cc` candidate/hit pipeline as its long-term engine.

### 4. Reuse `dbindex_paired.cc` as the paired DBAccel layer

Native paired `usearch_global` should load paired DB targets into `record_paired_s`-like state and then use:

- `dbindex_prepare_paired`
- `dbindex_addallsequences_paired`
- paired bitmap/matchlist accessors through `search_topscores_paired`

This is one of the biggest shared wins from the paired `cluster_unoise` work.

### 5. Move paired query reading into stock-shaped thread code

The stock-shaped paired worker path should handle:

- split paired queries (`R1` + positional `R2`)
- interleaved paired queries (`--interleaved`)
- per-thread reusable query buffers/state
- optional minus-strand preparation only if paired mode ever widens past `--strand plus`

For the current paired contract, `cmd_usearch_global()` still rejects non-plus strand, so the paired thread path can stay simpler than stock while preserving stock structure.

### 6. Move paired DB loading into `search_prep_paired`

`search_prep_paired` should become the stock-shaped home for:

- paired DB input opening/loading
- paired DB masking
- paired DB header consistency checks using the first whitespace-delimited token
- paired DBAccel/index preparation
- output file opening that matches stock `search_prep`

This is the right place to retire the current monolithic `tav_usearch_global()` setup section.

### 7. Port stock output flow instead of extension-specific output flow

`search_output_results_paired` should mirror stock responsibilities as closely as possible for the supported paired surface:

- decide `toreport`
- emit paired `alnout` / `userout` / `blast6out` / `uc`
- emit paired matched / notmatched FASTA outputs
- perform OTU-table top-hit or no-hit accounting
- update `dbmatched` abundance exactly once per accepted/weak hit

### 8. Keep unsupported outputs rejected in command dispatch

The native paired port does not need to force support for outputs that are still intentionally out of scope.

If paired mode continues to reject:

- `--fastapairs`
- `--qsegout`
- `--tsegout`
- `--samout`
- `--lcaout`
- custom `--userfields`

that is fine, as long as the rejection happens cleanly in `cmd_usearch_global()` and the native paired engine beneath it is stock-shaped for the supported surface.

### 9. Audit the extension code for reuse vs retirement

When porting, each `tav_usearch_global` helper should be classified as one of:

- reusable as-is in a stock-shaped module
- reusable after refolding into `search_paired.cc`
- should be deleted because `searchcore_paired.cc` or `dbindex_paired.cc` already covers it

### 10. Verify both thread modes early

Required smoke checks for the native paired port:

- `--threads 1`
- `--threads > 1`
- at least one hit-case where serial and parallel UC output agree
- at least one OTU-table smoke case where serial and parallel tables agree

This is not optional; the paired dbindex bug we already fixed is exactly the kind of issue that this catches.
