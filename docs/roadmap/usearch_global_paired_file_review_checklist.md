# Native Paired `usearch_global`: File-by-File Review Checklist

Status legend:

- `Good`: already close to stock shape or a clean shared wheel we should reuse.
- `Needs closer stock parity`: there is a clear stock counterpart, but the paired implementation still diverges more than necessary.
- `Should be deleted/refolded`: should be folded into a stock-shaped paired module or removed because a shared paired primitive already covers its job.

## `cpp/src/searchcore_paired.cc` / `cpp/src/searchcore_paired.h`

| Function / struct | Stock counterpart | Status | Review note |
| :-- | :-- | :-- | :-- |
| `hit_paired_s` | `hit_s` | Good | Shared paired hit state reused directly by the native `usearch_global` wrapper. |
| `searchinfo_s_paired` | `searchinfo_s` | Good | Shared paired query/search state reused directly by the native wrapper. |
| `make_hits_span_paired` | `make_hits_span` | Good | Shared internal utility; already stock-shaped. |
| `count_number_of_hits_to_keep_paired` | `count_number_of_hits_to_keep` | Good | Shared internal utility; already stock-shaped. |
| `copy_over_hits_to_be_kept_paired` | `copy_over_hits_to_be_kept` | Good | Shared internal utility; already stock-shaped. |
| `free_rejected_alignments_paired` | `free_rejected_alignments` | Good | Shared cleanup utility; already stock-shaped. |
| `hit_compare_byid_typed_paired` | `hit_compare_byid_typed` | Good | Reusable for stock-style paired hit ordering in `usearch_global`. |
| `hit_compare_bysize_typed_paired` | `hit_compare_bysize_typed` | Good | Mostly more relevant to clustering, but already stock-shaped and safe to keep shared. |
| `hit_compare_byid_paired` | `hit_compare_byid` | Good | Reusable wrapper. |
| `hit_compare_bysize_paired` | `hit_compare_bysize` | Good | Reusable wrapper. |
| `search_enough_kmers_paired` | `search_enough_kmers` | Good | Mainly a clustering helper, but already clean. |
| `paired_header_key_paired` | none in stock `searchcore.cc` | Good | Shared paired-only helper that compares the first whitespace-delimited header token without any `/1` or `/2` normalization. |
| `search_topscores_paired` | `search_topscores` | Good | Stock-shaped paired DBAccel screening reused directly by the native wrapper. |
| `search_acceptable_unaligned_paired` | `search_acceptable_unaligned` | Good | Native wrapper should use this instead of any extension-local unaligned filter stack. |
| `search_acceptable_aligned_paired` | `search_acceptable_aligned` | Good | Native wrapper should use this instead of any extension-local aligned filter stack. |
| `align_delayed_paired` | `align_delayed` | Good | Stock-shaped delayed batch aligner for the paired wrapper. |
| `search_onequery_paired` | `search_onequery` | Good | Core shared paired query search primitive used directly by the native port. |
| `search_findbest2_byid_paired` | `search_findbest2_byid` | Good | Not central for `usearch_global`, but already shared and stock-shaped. |
| `search_findbest2_bysize_paired` | `search_findbest2_bysize` | Good | Not central for `usearch_global`, but already shared and stock-shaped. |
| `search_joinhits_paired` | `search_joinhits` | Good | Stock-shaped paired hit-retention and sort path used directly by the native wrapper. |

## `cpp/src/dbindex_paired.cc` / `cpp/src/dbindex_paired.h`

| Function / object | Stock counterpart | Status | Review note |
| :-- | :-- | :-- | :-- |
| `dbindex_prepare_paired` | `dbindex_prepare` | Good | Reused as the paired DBAccel preparation phase in `search_prep_paired`. |
| `dbindex_addallsequences_paired` | `dbindex_addallsequences` | Good | Reused as the paired DB population step in `search_prep_paired`. |
| `dbindex_addsequence_paired` | `dbindex_addsequence` | Good | Shared paired insertion path; the serial insertion bug is already fixed here. |
| `dbindex_getbitmap_paired` | `dbindex_getbitmap` | Good | Reused by `search_topscores_paired`. |
| `dbindex_getmatchcount_paired` | `dbindex_getmatchcount` | Good | Reused by `search_topscores_paired`. |
| `dbindex_getmatchlist_paired` | `dbindex_getmatchlist` | Good | Reused by `search_topscores_paired`. |
| `dbindex_getmapping_paired` | `dbindex_getmapping` | Good | Reused to map dense paired slots back to logical paired targets. |
| `dbindex_getcount_paired` | `dbindex_getcount` | Good | Reused by `search_topscores_paired`. |
| `dbindex_free_paired` | `dbindex_free` | Good | Reused in `search_done_paired`. |

## `cpp/src/search_paired.cc` / `cpp/src/search.h`

| Function / helper | Stock counterpart | Status | Review note |
| :-- | :-- | :-- | :-- |
| `search_output_results_paired` | `search_output_results` | Good | Native paired per-query output path for the supported paired surface. |
| `search_query_paired` | `search_query` | Good | Wraps `search_onequery_paired`, `search_joinhits_paired`, and paired output emission in stock shape. |
| `search_thread_run_paired` | `search_thread_run` | Good | Native paired query-reader/worker loop for split and interleaved query input. |
| `search_thread_init_paired` | `search_thread_init` | Good | Allocates reusable per-thread paired search state in stock style. |
| `search_thread_exit_paired` | `search_thread_exit` | Good | Frees per-thread paired search state in stock style. |
| `search_thread_worker_paired` | `search_thread_worker` | Good | Thin stock-shaped wrapper. |
| `search_thread_worker_run_paired` | `search_thread_worker_run` | Good | Preserves stock thread orchestration. |
| `search_prep_paired` | `search_prep` | Good | Native home for paired DB loading, output opening, stock DB population, and paired dbindex preparation. |
| `search_done_paired` | `search_done` | Good | Native paired cleanup path. |
| `usearch_global_paired` | `usearch_global` | Good | Replaces the old extension engine on the paired CLI path. |
| `get_anchor_len_paired` | none in stock `search.cc` | Good | Necessary paired-only helper for one shared anchor length from both ends and `--fastq_trunclen`. |
| `is_catalog_file_paired` | none in stock `search.cc` | Good | Necessary paired-only input guard helper. |
| `write_fasta_record_paired` | none in stock `search.cc` | Good | Necessary paired-only FASTA output helper for split paired outputs. |
| `write_query_pair_split_paired` | none in stock `search.cc` | Good | Necessary paired-only split query FASTA output helper. |
| `write_db_pair_split_paired` | none in stock `search.cc` | Good | Necessary paired-only split DB FASTA output helper. |
