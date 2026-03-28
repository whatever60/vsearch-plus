# Native Paired `uchime3_denovo`: File-by-File Review Checklist

This document reviews the current state of the code that matters for a stock-shaped native paired `uchime3_denovo` port.

Reviewed files/modules:

- `src/chimera_paired.cc`
- `src/chimera.h`
- `src/searchcore_paired.cc`
- `src/dbindex_paired.cc`
- `src/tav_extension.cc` (legacy paired `uchime3_denovo` path)

The review standard is strict stock parity:

- `Good`: the paired function/module is already close to stock shape or is a clean shared wheel we should reuse.
- `Needs closer stock parity`: there is a clear stock counterpart, but the current paired implementation still diverges more than necessary.
- `Should be deleted/refolded`: the current function/module should be folded into a stock-shaped paired module or removed because a shared paired primitive already covers its job.
- `Missing`: a stock counterpart clearly exists, but the stock-shaped paired analog has not been created yet.

## File-Level Notes

Current implementation status:

- Paired `uchime3_denovo` now runs through `uchime3_denovo_paired()` in `src/chimera_paired.cc`.
- `cmd_chimera()` no longer depends on `tav_uchime3_denovo()` for the paired CLI path.
- The native engine reuses:
  - `src/dbindex_paired.cc`
  - stock `select_best_two_parents_from_match_matrix(...)`
- Full build passes with `make -C src -j8`.
- Split paired smoke test passes.
- Interleaved paired smoke test passes.
- Split and interleaved output agreement now holds at the base paired header level because paired header normalization trims `/1` and `/2` suffixes.
- A real teardown crash was fixed by sizing paired query-part buffers from the concatenated partition length instead of per-end ceilings.
- The paired-parent-selection split is now explicit:
  - `find_matches_one_end_paired` and `find_matches_paired` project the two per-end alignments into one combined match matrix on the concatenated paired axis
  - `find_best_parents_paired` reuses the stock `select_best_two_parents_from_match_matrix(...)` winner-selection logic on that paired matrix
  - `eval_parents_paired` is the function that actually introduces the paired breakpoint semantics (`LEFT_BREAK`, `MIDDLE_BREAK`, `RIGHT_BREAK`) by scanning R1 left-to-right and R2 right-to-left on the combined paired alignment

## `src/chimera.h`

| Function / declaration | Stock counterpart | Status | Review note |
| :-- | :-- | :-- | :-- |
| `Status` enum in header | stock-local enum moved to shared header | Good | Needed so stock and paired chimera modules share the same status type cleanly. |
| `select_best_two_parents_from_match_matrix(...)` declaration | stock helper export | Good | Lets paired code reuse the exact stock parent-winner primitive. |
| `uchime3_denovo_paired(...)` declaration | paired native entrypoint | Good | Correct native routing surface. |

## `src/dbindex_paired.cc`

| Function | Stock counterpart | Status | Review note |
| :-- | :-- | :-- | :-- |
| `dbindex_prepare_paired` | `dbindex_prepare` | Good | Reused directly by paired `uchime3_denovo`. |
| `dbindex_addsequence_paired` | `dbindex_addsequence` | Good | Reused directly for denovo parent-pool growth. |
| `dbindex_addallsequences_paired` | `dbindex_addallsequences` | Good | Shared wheel, though paired denovo uses incremental growth. |
| `dbindex_getbitmap_paired` | `dbindex_getbitmap` | Good | Reused by paired per-part candidate scoring. |
| `dbindex_getmatchcount_paired` | `dbindex_getmatchcount` | Good | Reused by paired per-part candidate scoring. |
| `dbindex_getmatchlist_paired` | `dbindex_getmatchlist` | Good | Reused by paired per-part candidate scoring. |
| `dbindex_getmapping_paired` | `dbindex_getmapping` | Good | Reused to map dense slots back to logical paired targets. |
| `dbindex_getcount_paired` | `dbindex_getcount` | Good | Reused by paired per-part candidate scoring. |
| `dbindex_free_paired` | `dbindex_free` | Good | Reused in paired teardown. |

## `src/chimera_paired.cc`

| Function | Stock counterpart | Status | Review note |
| :-- | :-- | :-- | :-- |
| `realloc_arrays_paired` | `realloc_arrays` | Good | Stock-shaped allocator; fixed to size part buffers from concatenated partition length, which removed the teardown heap corruption. |
| `reset_matches_paired` | `reset_matches` | Good | Direct paired analog. |
| `query_init_paired` | `query_init` | Good | Stock-shaped per-part search-state init, adapted to R1/R2 state. |
| `query_exit_paired` | `query_exit` | Good | Stock-shaped per-part teardown for paired search state. |
| `partition_query_paired` | `partition_query` | Good | Clean paired analog using one concatenated paired axis split into four stock-style parts. |
| `find_matches_one_end_paired` | part of `find_matches` | Good | Necessary paired-only helper to project one end's CIGAR onto the concatenated match matrix. |
| `find_matches_paired` | `find_matches` | Good | Paired analog that fills the stock-style match matrix across both ends; this is the first function where the two ends are fused into one parent-selection substrate. |
| `find_best_parents_paired` | `find_best_parents` | Good | Calls the same stock parent-winner primitive after building the paired match matrix; parent choice itself is still stock-style once the paired matrix exists. |
| `eval_parents_paired` | `eval_parents` | Good | This is the main paired-breakpoint function: it concatenates R1 and R2 alignments, scans breakpoint support with R1 forward and R2 reverse scan order, and assigns `LEFT_BREAK` / `MIDDLE_BREAK` / `RIGHT_BREAK`. |
| `chimera_thread_init_paired` | `chimera_thread_init` | Good | Stock-shaped worker init. |
| `chimera_thread_exit_paired` | `chimera_thread_exit` | Good | Stock-shaped worker teardown. |
| `chimera_thread_core_paired` | `chimera_thread_core` | Good | Native paired engine core; stock-shaped orchestration with paired candidate enumeration and paired full-query alignment. |
| `chimera_thread_worker_paired` | `chimera_thread_worker` | Good | Thin stock-shaped wrapper. |
| `chimera_threads_run_paired` | `chimera_threads_run` | Good | Preserves stock thread-driver shape. |
| `open_chimera_file_paired` | `open_chimera_file` | Good | Thin paired wrapper. |
| `close_chimera_file_paired` | `close_chimera_file` | Good | Thin paired wrapper. |
| `uchime3_denovo_paired` | denovo branch of `chimera` | Good | Native paired entrypoint; handles paired input loading, stock db/dbindex setup, worker run, summaries, and teardown. |

Necessary helper logic with no direct stock counterpart:

| Helper / logic | Why it exists | Status | Review note |
| :-- | :-- | :-- | :-- |
| paired header normalization inside `append_record` | Split and interleaved input must resolve to the same base paired ID. | Good | Necessary paired-only normalization. |
| `find_matches_one_end_paired` | Stock has one sequence axis; paired mode needs to map two end-specific CIGARs into one combined match matrix. | Good | Necessary paired-only helper. |
| breakpoint class labeling (`LEFT_BREAK`, `MIDDLE_BREAK`, `RIGHT_BREAK`) | Stock single-end `uchime3_denovo` has no paired breakpoint-class output. | Good | Necessary paired-only annotation. |
| split paired output files (`--chimeras2`, `--nonchimeras2`) | Stock has no R2 companion output files. | Good | Necessary paired-only command-surface logic. |

## `src/searchcore_paired.cc`

| Function / module | Stock counterpart | Status | Review note |
| :-- | :-- | :-- | :-- |
| paired search core and hit/search state | `src/searchcore.cc` | Good | Shared wheel already completed during paired clustering and paired `usearch_global`; reused here through paired dbindex semantics and stock primitives. |

## `src/tav_extension.cc` (legacy paired `uchime3_denovo` path)

| Function / object | Stock counterpart or target home | Status | Review note |
| :-- | :-- | :-- | :-- |
| `tav_uchime3_denovo` | `uchime3_denovo_paired` in `src/chimera_paired.cc` | Should be deleted/refolded | The paired CLI path no longer routes here. It is now legacy code. |
| old paired loader/search/output helpers for uchime3 | stock-shaped paired code in `src/chimera_paired.cc` and shared paired modules | Should be deleted/refolded | Useful as migration history, but no longer the native architecture. |

## Checklist Goal

This checklist is complete for the native paired command path:

- paired `uchime3_denovo` no longer relies on `tav_uchime3_denovo()` as its engine
- stock `chimera.cc` now has a real paired analog in `src/chimera_paired.cc`
- shared paired dbindex machinery is reused directly
- stock parent-selection logic is shared directly
- split and interleaved paired smoke tests both run successfully
- the remaining paired-only differences are algorithm/data-surface differences required by two-end input, not wrapper-architecture drift
