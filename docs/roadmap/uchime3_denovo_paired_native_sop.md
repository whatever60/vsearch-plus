# Native Paired `uchime3_denovo`: SOP and Implementation Checklist

This document captures the working SOP and the concrete implementation checklist for extending `vsearch --uchime3_denovo` to support paired-end reads natively in the same stock-shaped style we used for paired `cluster_unoise` and paired `usearch_global`.

It is based on:

- `docs/parity/uchime3_denovo_paired_parity.md`
- `docs/callstacks/uchime3_denovo.txt`
- `cpp/src/chimera.cc`
- `cpp/src/searchcore.cc`
- `cpp/src/searchcore_paired.cc`
- `cpp/src/dbindex_paired.cc`

## Current Status

- Paired CLI routing in `cmd_chimera()` now splits cleanly:
  - paired `--uchime3_denovo` input is routed to `uchime3_denovo_paired()`
  - single-end input is routed to stock `chimera()`
- The native paired implementation now lives in `cpp/src/chimera_paired.cc`.
  - the long-term paired engine is no longer `tav_uchime3_denovo()`
  - the paired command path now follows the stock `cpp/src/chimera.cc` structure:
    - `query_init_paired`
    - `query_exit_paired`
    - `partition_query_paired`
    - `chimera_thread_init_paired`
    - `chimera_thread_exit_paired`
    - `chimera_thread_core_paired`
    - `chimera_thread_worker_paired`
    - `chimera_threads_run_paired`
    - `uchime3_denovo_paired`
- The paired engine reuses the shared paired search/index wheels directly:
  - `cpp/src/searchcore_paired.cc`
  - `cpp/src/dbindex_paired.cc`
- The stock parent-selection primitive is now shared directly:
  - `select_best_two_parents_from_match_matrix(...)` is declared in `cpp/src/chimera.h`
  - both stock `find_best_parents()` and paired `find_best_parents_paired()` call it
- Full build passes with `make -C src -j8`.
- Native paired smoke tests now pass for:
  - split paired input
  - interleaved paired input
  - `--tabbedout`
  - `--uchimeout`
  - `--uchimealns`
  - `--chimeras` / `--chimeras2`
  - `--nonchimeras` / `--nonchimeras2`
  - `--chimeras_tsv`
  - `--nonchimeras_tsv`
- Split and interleaved paired input now preserve separate R1 and R2 headers while requiring the same first whitespace-delimited token.
- A real heap-corruption bug in the paired native engine was fixed:
  - query-part buffers were initially sized from per-end ceilings instead of the concatenated partition length
  - that was unsafe for concatenated paired partitioning and caused teardown-time allocator crashes
  - `realloc_arrays_paired()` now sizes part buffers from the full concatenated query partition length

## Purpose

The goal is not to keep growing the now-removed legacy `tav_extension.cc` path.

The goal is to reproduce stock `uchime3_denovo` structure and behavior as faithfully as possible, while replacing one denoised sequence with one synchronized denoised pair.

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

Write down how each stock single-sequence `uchime3_denovo` step becomes paired-read logic before moving code.

Main semantic reference:

- `docs/parity/uchime3_denovo_paired_parity.md`

Cross-check against the other paired docs so the command family stays internally consistent:

- `docs/parity/fastq_filter_paired_parity.md`
- `docs/parity/cluster_unoise_paired_parity.md`
- `docs/parity/usearch_global_paired_parity.md`

### 2. Freeze the complete stock callstack first

Treat stock `uchime3_denovo` as the only structural ground truth and document the complete callstack before porting.

Stock callstack reference:

- `docs/callstacks/uchime3_denovo.txt`

Important nuance:

- stock denovo `chimera()` forces `opt_threads = 1`
- so the paired native port should preserve the stock threaded-driver shape, but it does not need a separate multi-worker denovo architecture

### 3. Port by stock module boundary

Duplicate stock modules that sit on the `uchime3_denovo` callstack and need paired-aware state.

Target module mapping:

- `cpp/src/chimera.cc` -> `cpp/src/chimera_paired.cc`
- shared `cpp/src/searchcore.cc` paired search pieces -> `cpp/src/searchcore_paired.cc`
- shared `cpp/src/dbindex.cc` paired index surface -> `cpp/src/dbindex_paired.cc`

Do not treat the now-removed legacy `tav_extension.cc` path as the final native location for paired `uchime3_denovo`.
It was a semantic reference and temporary implementation source, not the end-state architecture.

### 4. Keep paired names and function shapes close to stock

The paired port should look like the stock code with `_paired` appended.

Examples:

- `query_init` -> `query_init_paired`
- `query_exit` -> `query_exit_paired`
- `partition_query` -> `partition_query_paired`
- `find_matches` -> `find_matches_paired`
- `find_best_parents` -> `find_best_parents_paired`
- `eval_parents` -> `eval_parents_paired`
- `chimera_thread_init` -> `chimera_thread_init_paired`
- `chimera_thread_exit` -> `chimera_thread_exit_paired`
- `chimera_thread_core` -> `chimera_thread_core_paired`
- `chimera_threads_run` -> `chimera_threads_run_paired`
- `chimera` denovo route -> `uchime3_denovo_paired`

The body should stay almost exactly the same unless the difference is explicitly caused by replacing one sequence with one synchronized pair.

### 5. Replace one denoised sequence with one denoised pair systematically

This is the main transformation rule.

Examples:

- one query sequence -> `query_seq_r1` + `query_seq_r2`
- one query length -> `query_len_r1` + `query_len_r2`
- one candidate seqno -> one paired target index plus mapping to `target_seqno_r1` / `target_seqno_r2`
- one alignment CIGAR array -> one per-end CIGAR array for R1 and R2
- one query axis -> one concatenated paired axis for parent selection and breakpoint scoring

The stock control flow should stay intact. Only the data model becomes paired.

### 6. Reuse shared paired search/index wheels aggressively

For `uchime3_denovo`, the most important reuse boundary is the search/index layer we already finished for paired clustering and paired `usearch_global`.

Shared pieces already completed:

- `cpp/src/searchcore_paired.cc`
  - paired top-score screening
  - paired unaligned/aligned filters
  - paired delayed alignment path
  - stock-shaped paired hit structs and search state
- `cpp/src/dbindex_paired.cc`
  - paired k-mer index preparation
  - paired dense mapping/count/bitmap accessors
  - paired incremental insertion for denovo growth

So the native paired `uchime3_denovo` port should concentrate on the stock `cpp/src/chimera.cc` wrapper/orchestration layer:

- paired input loading
- paired partitioning
- paired candidate enumeration for each part
- full-query paired alignment to retained parents
- stock winner-window parent selection
- stock-style `h` optimization and final decision
- paired output materialization

### 7. Reuse stock low-level primitives directly where possible

The paired port should continue reusing stock building blocks:

- `unique_init`, `unique_exit`, `unique_count`
- `search16_init`, `search16_qprep`, `search16`, `search16_exit`
- `LinearMemoryAligner`
- `db_getsequence`, `db_getsequencelen`, `db_getabundance`, `db_getheader`
- `fasta_print_general`
- `select_best_two_parents_from_match_matrix`

If a stock helper already does the right thing once given stock-like paired state, reuse it instead of building a parallel helper stack.

### 8. Keep the denovo growth rule stock-shaped

For denovo chimera detection, both stock and paired native mode should grow the parent search space monotonically:

- non-chimeric queries are added to the parent index
- chimeric queries are not added

In the paired port, the unit added is a paired target record through `dbindex_addsequence_paired(...)`.

### 9. Keep paired-only differences explicit and narrow

Paired-only logic belongs in a few well-defined places:

- paired FASTX loading and synchronization
- paired first-token header consistency checks across split/interleaved input
- paired candidate scoring across R1 and R2 k-mers
- paired full-query alignment across both ends
- breakpoint classification (`LEFT_BREAK`, `MIDDLE_BREAK`, `RIGHT_BREAK`)
- split paired output files (`--chimeras2`, `--nonchimeras2`)

Everything else should remain as close to stock `chimera.cc` as possible.

### 10. Audit every paired function line by line

For each paired function and struct:

1. identify the exact stock counterpart
2. compare signature shape
3. compare body structure
4. justify every extra line by paired semantics
5. verify stock/shared paired primitives are reused where possible
6. delete helper layers that only exist because the code came from the now-removed legacy `tav_extension.cc` path

### 11. Only optimize after structural parity is correct

Performance work comes after semantic and structural parity.

For `uchime3_denovo`, correctness checks come first:

- split vs interleaved input should agree
- split output pairs should stay synchronized
- denovo parent-pool growth should match stock semantics
- paired teardown should be clean and memory-safe

## Target Callstack

```text
cmd_chimera
-> uchime3_denovo_paired
-> open outputs
-> load paired FASTX input
-> sort paired records by abundance/header/sequence
-> build paired stock db backing store
-> dbindex_prepare_paired
-> chimera_threads_run_paired
-> chimera_thread_worker_paired
-> chimera_thread_core_paired
-> chimera_thread_init_paired
-> query_init_paired
-> partition_query_paired
-> paired part screening via paired dbindex counts/bitmaps
-> full-query paired alignment to candidate parents
-> find_best_parents_paired
-> select_best_two_parents_from_match_matrix
-> eval_parents_paired
-> paired output materialization
-> dbindex_addsequence_paired for non-chimeras
-> chimera_thread_exit_paired
-> query_exit_paired
-> dbindex_free_paired
-> db_free
-> close outputs
```

## Function Mapping Checklist

```text
query_init                         -> query_init_paired
query_exit                         -> query_exit_paired
partition_query                    -> partition_query_paired
find_matches                       -> find_matches_paired
find_best_parents                  -> find_best_parents_paired
eval_parents                       -> eval_parents_paired
chimera_thread_init                -> chimera_thread_init_paired
chimera_thread_exit                -> chimera_thread_exit_paired
chimera_thread_core                -> chimera_thread_core_paired
chimera_thread_worker              -> chimera_thread_worker_paired
chimera_threads_run                -> chimera_threads_run_paired
chimera (paired denovo route)      -> uchime3_denovo_paired
```
