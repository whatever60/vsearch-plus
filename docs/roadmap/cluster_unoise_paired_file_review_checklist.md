# Native Paired `cluster_unoise`: File-by-File Review Checklist

This document reviews the current state of:

- `src/searchcore_paired.cc`
- `src/cluster_paired.cc`

The review standard is strict stock parity:

- `Good`: the paired function is already close to stock shape, with only explicit one-read to two-read differences.
- `Needs closer stock parity`: there is a clear stock counterpart, but the current paired implementation still diverges more than necessary.
- `Should be deleted/refolded`: the function should be inlined, merged into a stock-shaped function, or moved to a more appropriate shared module. This usually means it has no clean stock counterpart and adds an unnecessary abstraction layer.

## File-Level Notes

Current implementation status:

- Full project build currently passes with `make -C src -j8`.
- The paired native path now builds as part of the normal tree, rather than only through isolated syntax checks.
- Native paired smoke tests now pass in both execution modes:
  - `--threads 2`
  - `--threads 1`
- A paired core-result smoke test succeeds with `--uc` and `--otutabout`.
- A tiny paired hit-case smoke test produces real `H` lines in UC output, and the serial and parallel outputs now agree.
- A serial-only paired dbindex bug was fixed while closing this checklist:
  - `dbindex_addsequence_paired()` now consumes the R1 unique-kmer list before reusing the shared unique-count handle for R2, so serial clustering no longer misses hits that parallel same-round reconciliation used to mask.
- The checklist below is now a completion record for the core paired `cluster_unoise` engine.

### `src/searchcore_paired.cc`

Current structural state:

- `hit_paired_s` is now stock-shaped enough for the paired port.
  - It keeps the stock aggregate ranking/filter fields at pair level.
  - It embeds one stock `struct hit` for each end, including raw `char *` alignment ownership.
- `searchinfo_s_paired` is now much closer to stock `searchinfo_s`.
  - query buffers, k-mer counters, and hit storage use stock-style raw allocations
  - query header handling uses a stock-style raw pointer plus length instead of owning `std::string` state
  - paired-only differences are now the explicit two-end fields

### `src/cluster_paired.cc`

Current structural state:

- `clusterinfo_s_paired` now retains per-end CIGAR ownership in stock style.
  - this closes the earlier gap where aligned state was thrown away too early
  - the paired tail still writes paired-specific outputs, but it now does so on top of retained cluster alignment state rather than summary-only metadata
- `cluster_core_results_hit_paired()` and `cluster_core_results_nohit_paired()` perform stock-like immediate output work for the supported paired surface:
  - immediate `S/H` UC emission
  - immediate OTU-table accumulation
- query/thread working state now uses stock-like raw arrays rather than vector-backed storage.
- the remaining differences in `cluster_paired()` are now paired input/output analogs, not unresolved parity gaps in the engine core.

## `src/searchcore_paired.cc`

| Function | Stock counterpart | Status | Review note |
| :-- | :-- | :-- | :-- |
| `make_hits_span_paired` | `make_hits_span` | Good | Shape and behavior are already one-to-one. |
| `count_number_of_hits_to_keep_paired` | `count_number_of_hits_to_keep` | Good | Matches stock structure closely; only paired hit type differs. |
| `copy_over_hits_to_be_kept_paired` | `copy_over_hits_to_be_kept` | Good | Matches stock structure closely. |
| `free_rejected_alignments_paired` | `free_rejected_alignments` | Good | It now follows the stock ownership model closely by freeing rejected raw per-end CIGAR pointers before joinhits returns. |
| `hit_compare_byid_typed_paired` | `hit_compare_byid_typed` | Good | Already very close to stock shape. |
| `hit_compare_bysize_typed_paired` | `hit_compare_bysize_typed` | Good | It now uses stock `db_getabundance(...)` through the paired hit’s mapped R1 target sequence number, without the old translation-unit global workaround. |
| `hit_compare_byid_paired` | `hit_compare_byid` | Good | Thin stock-shaped wrapper. |
| `hit_compare_bysize_paired` | `hit_compare_bysize` | Good | Thin stock-shaped wrapper; the underlying typed comparator is now much closer to stock. |
| `search_enough_kmers_paired` | `search_enough_kmers` | Good | The only meaningful difference is summing `kmersamplecount_r1 + kmersamplecount_r2`. |
| `search_topscores_paired` | `search_topscores` | Good | This is now structurally close to stock and reuses a stock-shaped paired dbindex API. The paired-specific difference is explicit two-end k-mer accumulation. |
| `search_acceptable_unaligned_paired` | `search_acceptable_unaligned` | Good | This is a clean paired analog: stock numeric filters are reused, and the only differences are paired length/abundance aggregation and paired prefix/suffix/self/selfid checks. |
| `search_acceptable_aligned_paired` | `search_acceptable_aligned` | Good | This is now much closer to stock: it consumes stock-shaped aggregate hit fields and writes back stock-style identity metrics. |
| `align_delayed_paired` | `align_delayed` | Good | It now follows the stock delayed-alignment flow closely, uses raw per-end `struct hit` state directly, and trims/alignment-filters each end without temporary shim objects. The remaining duplication is the explicit paired-end analogue of stock’s single-end bookkeeping. |
| `search_onequery_paired` | `search_onequery` | Good | The function is already stock-shaped at a high level: qprep, LMA setup, unique k-mers, top scores, delayed alignments, cleanup. |
| `search_findbest2_byid_paired` | `search_findbest2_byid` | Good | Already very close to stock. |
| `search_findbest2_bysize_paired` | `search_findbest2_bysize` | Good | Now inherits the improved stock-like size comparator with direct `db_getabundance(...)` lookup through mapped paired targets. |
| `search_joinhits_paired` | `search_joinhits` | Good | The function is stock-shaped end-to-end: keep accepted/weak hits, free rejected alignments, then sort in stock order. |

## `src/cluster_paired.cc`

| Function | Stock counterpart | Status | Review note |
| :-- | :-- | :-- | :-- |
| `compare_byclusterno_paired` | `compare_byclusterno` | Good | Already close to stock. |
| `compare_byclusterabundance_paired` | `compare_byclusterabundance` | Good | Already close to stock. |
| `cluster_query_core_paired` | `cluster_query_core` | Good | This is a clean paired analog of the stock query core: load query state, handle reverse strand, call `search_onequery_paired()`. |
| `cluster_worker_paired` | `cluster_worker` | Good | Already close to stock. |
| `threads_worker_paired` | `threads_worker` | Good | Already close to stock. |
| `threads_wakeup_paired` | `threads_wakeup` | Good | Already close to stock. |
| `threads_init_paired` | `threads_init` | Good | Now uses stock-like raw thread-info allocation and the stock thread protocol. |
| `threads_exit_paired` | `threads_exit` | Good | Already close to stock. |
| `cluster_query_init_paired` | `cluster_query_init` | Good | Query header/sequence/k-mer/hit state is now initialized with stock-style raw allocation and raw pointer ownership, with only explicit two-end duplication beyond stock. |
| `cluster_query_exit_paired` | `cluster_query_exit` | Good | The cleanup structure is close to stock, with paired-specific duplication for the two ends. |
| `cluster_core_results_hit_paired` | `cluster_core_results_hit` | Good | For the currently supported paired output surface, this now behaves much more like stock: it performs immediate `H`-line UC emission and immediate OTU-table accumulation. |
| `cluster_core_results_nohit_paired` | `cluster_core_results_nohit` | Good | For the currently supported paired output surface, this now behaves much more like stock: it performs immediate `S`-line UC emission and immediate OTU-table accumulation. |
| `clear_hits_paired` | closest stock analog is the explicit post-query alignment cleanup loops in `cluster_core_parallel` and `cluster_core_serial` | Good | This helper now frees raw per-end CIGAR ownership in the same places where stock frees single-end alignment strings. |
| `cluster_core_parallel_paired` | `cluster_core_parallel` | Good | The core loop now matches stock structure closely: stock-like raw query arrays, stock-like same-round reconciliation, correct aligned-hit cleanup, retained cluster CIGAR ownership, last-length checks, explicit DB/index insertion, and immediate result handoff. |
| `cluster_core_serial_paired` | `cluster_core_serial` | Good | The serial loop now mirrors stock structure closely and matches the parallel paired result surface. A serial-only dbindex insertion bug was fixed, and the serial smoke outputs now match the parallel smoke outputs on the hit-case. |
| `cluster_paired` | `cluster` | Good | The function still has the necessary paired input/output analogs, but the engine structure now mirrors stock: paired load/filter/sort, db reset/index prep, stock-shaped core dispatch, retained cluster state, cluster sorting, paired centroid/cluster writing, OTU-table emission, and final cleanup. |
| `cluster_unoise_paired` | `cluster_unoise` | Good | Thin stock-shaped wrapper. |

## Checklist Status

Core engine parity items in this checklist are now closed.

The main paired-end differences that remain in the code are intentional command-surface differences:

- paired input loading instead of stock single-end `db_read(...)`
- paired centroid/cluster output materialization
- paired dbindex keying and pair-to-end mapping

Those are native paired extensions of the stock engine, not unresolved gaps in the `cluster_unoise` core callstack.
