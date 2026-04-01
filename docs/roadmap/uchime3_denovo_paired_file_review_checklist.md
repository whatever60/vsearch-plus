# Native Paired `uchime3_denovo`: File-by-File Review Checklist

Status legend:

- `Good`: already close to stock shape or a clean shared wheel we should reuse.
- `Needs closer stock parity`: there is a clear stock counterpart, but the paired implementation still diverges more than necessary.
- `Should be deleted/refolded`: should be folded into a stock-shaped paired module or removed because a shared paired primitive already covers its job.

## `cpp/src/chimera.cc` / `cpp/src/chimera.h`

| Function / declaration | Stock counterpart | Status | Review note |
| :-- | :-- | :-- | :-- |
| `Status` enum | `Status` enum | Good | Shared stock status type used cleanly by both stock and paired chimera modules. |
| `select_best_two_parents_from_match_matrix(...)` declaration | `select_best_two_parents_from_match_matrix(...)` | Good | Lets paired code reuse the exact stock parent-winner primitive. |

## `cpp/src/dbindex_paired.cc` / `cpp/src/dbindex_paired.h`

| Function / object | Stock counterpart | Status | Review note |
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

## `cpp/src/searchcore_paired.cc` / `cpp/src/searchcore_paired.h`

| Function / module | Stock counterpart | Status | Review note |
| :-- | :-- | :-- | :-- |
| paired search core and hit/search state | `cpp/src/searchcore.cc` | Good | Shared wheel already completed during paired clustering and paired `usearch_global`; reused here through paired dbindex semantics and stock primitives. |

## `cpp/src/chimera_paired.cc` / `cpp/src/chimera_paired.h`

| Function / helper | Stock counterpart | Status | Review note |
| :-- | :-- | :-- | :-- |
| `uchime3_denovo_paired` | denovo branch of `chimera` | Good | Native paired entrypoint declared in `chimera_paired.h`; handles paired input loading, stock db/dbindex setup, worker run, summaries, and teardown. |
| `realloc_arrays_paired` | `realloc_arrays` | Good | Stock-shaped allocator; fixed to size part buffers from concatenated partition length, which removed the teardown heap corruption. |
| `reset_matches_paired` | `reset_matches` | Good | Direct paired analog. |
| `query_init_paired` | `query_init` | Good | Stock-shaped per-part search-state init adapted to R1/R2 state. |
| `query_exit_paired` | `query_exit` | Good | Stock-shaped per-part teardown for paired search state. |
| `partition_query_paired` | `partition_query` | Good | Clean paired analog using one concatenated paired axis split into four stock-style parts. |
| `find_matches_one_end_paired` | part of `find_matches` | Good | Necessary paired-only helper to project one end's CIGAR onto the concatenated match matrix. |
| `find_matches_paired` | `find_matches` | Good | Fills the stock-style match matrix across both ends; this is where the two ends are fused into one parent-selection substrate. |
| `find_best_parents_paired` | `find_best_parents` | Good | Reuses the stock `select_best_two_parents_from_match_matrix(...)` winner-selection logic after building the paired match matrix. |
| `eval_parents_paired` | `eval_parents` | Good | Main paired-breakpoint function: concatenates R1 and R2 alignments, scans R1 left-to-right and R2 right-to-left, and assigns `LEFT_BREAK`, `MIDDLE_BREAK`, or `RIGHT_BREAK`. |
| `chimera_thread_init_paired` | `chimera_thread_init` | Good | Stock-shaped worker init. |
| `chimera_thread_exit_paired` | `chimera_thread_exit` | Good | Stock-shaped worker teardown. |
| `chimera_thread_core_paired` | `chimera_thread_core` | Good | Native paired engine core with stock-shaped orchestration plus paired candidate enumeration and full-query alignment. |
| `chimera_thread_worker_paired` | `chimera_thread_worker` | Good | Thin stock-shaped wrapper. |
| `chimera_threads_run_paired` | `chimera_threads_run` | Good | Preserves stock thread-driver shape. |
| `open_chimera_file_paired` | `open_chimera_file` | Good | Thin paired wrapper. |
| `close_chimera_file_paired` | `close_chimera_file` | Good | Thin paired wrapper. |
| paired first-token header consistency check inside input loading | none in stock `chimera.cc` | Good | Necessary paired-only validation so split and interleaved input preserve separate end headers while still requiring the same first whitespace-delimited token. |
| breakpoint class labeling (`LEFT_BREAK`, `MIDDLE_BREAK`, `RIGHT_BREAK`) | none in stock `chimera.cc` | Good | Necessary paired-only breakpoint annotation. |
| split paired output files (`--chimeras2`, `--nonchimeras2`) | none in stock `chimera.cc` | Good | Necessary paired-only command-surface logic. |
