# Native Paired `usearch_global`: File-by-File Review Checklist

This document reviews the current state of the code that matters for a stock-shaped native paired `usearch_global` port.

Reviewed files/modules:

- `src/searchcore_paired.cc`
- `src/dbindex_paired.cc`
- `src/search_paired.cc`
- `src/tav_extension.cc` (legacy paired `usearch_global` path)

The review standard is strict stock parity:

- `Good`: the paired function/module is already close to stock shape or is a clean shared wheel we should reuse.
- `Needs closer stock parity`: there is a clear stock counterpart, but the current paired implementation still diverges more than necessary.
- `Should be deleted/refolded`: the current function/module should be folded into a stock-shaped paired module or removed because a shared paired primitive already covers its job.
- `Missing`: a stock counterpart clearly exists, but the stock-shaped paired analog has not been created yet.

## File-Level Notes

Current implementation status:

- Paired `usearch_global` now runs through `usearch_global_paired()` in `src/search_paired.cc`.
- `cmd_usearch_global()` no longer depends on `tav_usearch_global()` for the paired command path.
- The native wrapper reuses the shared paired engine directly:
  - `src/searchcore_paired.cc`
  - `src/dbindex_paired.cc`
- Supported paired wrappers were verified with:
  - split paired input, `--threads 1`
  - split paired input, `--threads 2`
  - interleaved paired input, `--threads 2`
- Supported serial/parallel outputs matched exactly for:
  - `uc`
  - `otutabout`
  - `userout`
  - `blast6out`
  - `matched` / `match2`
  - `notmatched` / `notmatched2`
  - `dbmatched` / `dbmatched2`
  - `dbnotmatched` / `dbnotmatched2`
- `alnout` body also matched between `--threads 1` and `--threads 2`
  - only the stock command-header line differed, because the output filenames and thread count were intentionally different between the two runs

## `src/searchcore_paired.cc`

These are the shared wheels we should actively reuse for paired `usearch_global`.

| Function | Stock counterpart | Status | Review note |
| :-- | :-- | :-- | :-- |
| `make_hits_span_paired` | `make_hits_span` | Good | Shared internal utility; already stock-shaped. |
| `count_number_of_hits_to_keep_paired` | `count_number_of_hits_to_keep` | Good | Shared internal utility; already stock-shaped. |
| `copy_over_hits_to_be_kept_paired` | `copy_over_hits_to_be_kept` | Good | Shared internal utility; already stock-shaped. |
| `free_rejected_alignments_paired` | `free_rejected_alignments` | Good | Shared cleanup utility; already stock-shaped. |
| `hit_compare_byid_typed_paired` | `hit_compare_byid_typed` | Good | Reusable for stock-style paired hit ordering in `usearch_global`. |
| `hit_compare_bysize_typed_paired` | `hit_compare_bysize_typed` | Good | Mostly more relevant to clustering, but already stock-shaped and safe to keep shared. |
| `hit_compare_byid_paired` | `hit_compare_byid` | Good | Reusable wrapper. |
| `hit_compare_bysize_paired` | `hit_compare_bysize` | Good | Reusable wrapper. |
| `search_enough_kmers_paired` | `search_enough_kmers` | Good | Mainly a clustering helper, but already clean. |
| `search_topscores_paired` | `search_topscores` | Good | This is one of the most important shared wins for paired `usearch_global`; it already gives us stock-shaped paired DBAccel screening. |
| `search_acceptable_unaligned_paired` | `search_acceptable_unaligned` | Good | This should replace `paired_unaligned_filters_pass()` in the final native `usearch_global` path. |
| `search_acceptable_aligned_paired` | `search_acceptable_aligned` | Good | This should replace `paired_aligned_filters_pass()` in the final native `usearch_global` path. |
| `align_delayed_paired` | `align_delayed` | Good | This is the stock-shaped delayed batch aligner the native `usearch_global` path should use. |
| `search_onequery_paired` | `search_onequery` | Good | This is the core shared paired query search primitive we want the native port to call directly. |
| `search_findbest2_byid_paired` | `search_findbest2_byid` | Good | Not central for `usearch_global`, but already shared and stock-shaped. |
| `search_findbest2_bysize_paired` | `search_findbest2_bysize` | Good | Not central for `usearch_global`, but already shared and stock-shaped. |
| `search_joinhits_paired` | `search_joinhits` | Good | This is the stock-shaped paired hit-retention/sort path the native `usearch_global` wrapper should call directly. |

## `src/dbindex_paired.cc`

This file is already the paired DBAccel layer we want to reuse.

| Function | Stock counterpart | Status | Review note |
| :-- | :-- | :-- | :-- |
| `dbindex_prepare_paired` | `dbindex_prepare` | Good | Reusable as the paired DBAccel preparation phase in `search_prep_paired`. |
| `dbindex_addallsequences_paired` | `dbindex_addallsequences` | Good | Reusable as the paired DB population step in `search_prep_paired`. |
| `dbindex_addsequence_paired` | `dbindex_addsequence` | Good | Shared paired insertion path; the serial insertion bug has already been fixed here. |
| `dbindex_getbitmap_paired` | `dbindex_getbitmap` | Good | Reusable by `search_topscores_paired`. |
| `dbindex_getmatchcount_paired` | `dbindex_getmatchcount` | Good | Reusable by `search_topscores_paired`. |
| `dbindex_getmatchlist_paired` | `dbindex_getmatchlist` | Good | Reusable by `search_topscores_paired`. |
| `dbindex_getmapping_paired` | `dbindex_getmapping` | Good | Reusable by `search_topscores_paired` and paired result/output code. |
| `dbindex_getcount_paired` | `dbindex_getcount` | Good | Reusable by `search_topscores_paired`. |
| `dbindex_free_paired` | `dbindex_free` | Good | Reusable in `search_done_paired`. |

## `src/tav_extension.cc` (paired `usearch_global` section)

These are the main pieces that should be retired or refolded into stock-shaped paired modules.

| Function / object | Stock counterpart or target home | Status | Review note |
| :-- | :-- | :-- | :-- |
| `TavRecord` | `record_paired_s` or stock-shaped paired query/DB record state | Should be deleted/refolded | The native port should converge on one shared paired record type rather than keep a usearch-global-specific record object. |
| `TavHit` | `hit_paired_s` | Should be deleted/refolded | The native port should use the shared stock-shaped paired hit object instead of a separate extension-only hit type. |
| `mask_sequence` lambda | local stock-shaped helper inside `search_prep_paired` / `search_thread_run_paired` | Should be deleted/refolded | Useful behavior, but not as a long-term extension-only helper. |
| `trim_pair_suffix` lambda | local stock-shaped helper inside paired DB/query loader path | Should be deleted/refolded | Useful behavior, but belongs in the stock-shaped paired wrapper layer. |
| `append_db_record` lambda | `search_prep_paired` paired DB loader path | Should be deleted/refolded | This should move under stock-shaped paired DB preparation. |
| `load_paired_db_from_interleaved_fastx` | paired `search_prep_paired` DB loading | Should be deleted/refolded | Necessary behavior, wrong long-term location. |
| `load_paired_db_from_split_fastx` | paired `search_prep_paired` DB loading | Should be deleted/refolded | Necessary behavior, wrong long-term location. |
| `write_fasta_record` / `write_query_pair_split` / `write_db_pair_split` | paired `search_output_results_paired` / paired tail | Should be deleted/refolded | These are output helpers that should move under the stock-shaped paired wrapper layer. |
| old extension `search_onequery_paired(...)` | `src/searchcore_paired.cc:search_onequery_paired` | Should be deleted/refolded | The native port should stop using the older extension search core here. |
| old extension `search_joinhits_paired(...)` | `src/searchcore_paired.cc:search_joinhits_paired` | Should be deleted/refolded | Same issue as above. |
| `align_end_batch_stock_style` | `align_delayed_paired` or lower-level shared paired alignment primitive | Should be deleted/refolded | Useful during prototyping, but the native port should drive stock-shaped delayed alignment through `searchcore_paired.cc`. |
| `paired_unaligned_filters_pass` | `search_acceptable_unaligned_paired` | Should be deleted/refolded | The shared paired search core already covers this logic in stock shape. |
| `paired_aligned_filters_pass` | `search_acceptable_aligned_paired` | Should be deleted/refolded | Same issue; the native port should stop carrying a second aligned-filter stack. |
| `tav_usearch_global` | `usearch_global_paired` in `src/search_paired.cc` | Should be deleted/refolded | The paired CLI path no longer routes here. It is now legacy code that should eventually be removed after the rest of the extension file is retired. |

## `src/search_paired.cc`

This file now exists and is the native stock-shaped paired wrapper.

| Planned function | Stock counterpart | Status | Review note |
| :-- | :-- | :-- | :-- |
| `search_output_results_paired` | `search_output_results` | Good | Native paired per-query output path is in place for the supported paired surface. |
| `search_query_paired` | `search_query` | Good | Wraps `search_onequery_paired` + `search_joinhits_paired` + paired output emission in stock shape. |
| `search_thread_run_paired` | `search_thread_run` | Good | Native paired query-reader/worker loop now drives split and interleaved query input. |
| `search_thread_init_paired` | `search_thread_init` | Good | Allocates reusable per-thread paired search state in stock style. |
| `search_thread_exit_paired` | `search_thread_exit` | Good | Frees per-thread paired search state in stock style. |
| `search_thread_worker_paired` | `search_thread_worker` | Good | Thin stock-shaped wrapper. |
| `search_thread_worker_run_paired` | `search_thread_worker_run` | Good | Preserves stock thread orchestration. |
| `search_prep_paired` | `search_prep` | Good | Native home for paired DB loading, output opening, stock DB population, and paired dbindex preparation. |
| `search_done_paired` | `search_done` | Good | Native paired cleanup path. |
| `usearch_global_paired` | `usearch_global` | Good | Replaces `tav_usearch_global` on the paired CLI path. |

Necessary helper functions with no direct stock counterpart in `src/search.cc`:

| Helper | Why it exists | Status | Review note |
| :-- | :-- | :-- | :-- |
| `trim_pair_suffix_paired` | Paired DB records need one shared base header from `/1` and `/2` labels. | Good | Necessary paired-only normalization helper. `trim_header_to_id_paired` was refolded into this function. |
| `get_anchor_len_paired` | Paired reads need one shared anchor length derived from both ends and `--fastq_trunclen`. | Good | Necessary paired-only query/DB normalization helper. |
| `is_catalog_file_paired` | Paired mode explicitly rejects catalog-style DB input. | Good | Necessary paired-only input guard helper. |
| `write_fasta_record_paired` | Stock has no split paired FASTA writer for `--match2`/`--dbmatched2`-style outputs. | Good | Necessary paired-only output helper. |
| `write_query_pair_split_paired` | Paired mode writes split query FASTA outputs. | Good | Necessary paired-only output helper. |
| `write_db_pair_split_paired` | Paired mode writes split DB FASTA outputs. | Good | Necessary paired-only output helper. |

Refolded helpers removed from `src/search_paired.cc`:

- `trim_header_to_id_paired`
  - folded into `trim_pair_suffix_paired`
- `mask_sequence_paired`
  - folded into the paired DB-load path inside `search_prep_paired()`

## Port Result

- `src/search_paired.cc` is now the native paired wrapper.
- Paired DB loading now lives in `search_prep_paired()`.
- Paired query reading now lives in `search_thread_run_paired()`.
- The wrapper now calls `search_onequery_paired()` and `search_joinhits_paired()` directly.
- Supported paired output emission now lives in `search_output_results_paired()`.
- `src/tav_extension.cc` still contains the legacy paired implementation, but it is no longer on the paired CLI path.

## Checklist Goal

This checklist is complete for the native paired command path:

- paired `usearch_global` no longer relies on `tav_usearch_global()` as its engine
- stock `search.cc` now has a real paired analog in `src/search_paired.cc`
- shared paired searchcore/dbindex modules are reused directly
- serial and parallel paired smoke outputs agree on hit-case and OTU-table tests
- the remaining paired-only differences are command-surface differences, not wrapper architecture differences
