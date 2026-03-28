# Native Paired `cluster_unoise`: File-by-File Review Checklist

Status legend:

- `Good`: already close to stock shape, with only explicit one-read to two-read differences.
- `Needs closer stock parity`: there is a clear stock counterpart, but the paired implementation still diverges more than necessary.
- `Should be deleted/refolded`: should be inlined, merged into a stock-shaped function, or moved to a more appropriate shared module.

## `src/searchcore_paired.cc` / `src/searchcore_paired.h`

| Function / struct | Stock counterpart | Status | Review note |
| :-- | :-- | :-- | :-- |
| `hit_paired_s` | `hit_s` | Good | Pair-level analog of stock hit state, with one embedded stock `struct hit` per end plus aggregate ranking/filter fields. |
| `searchinfo_s_paired` | `searchinfo_s` | Good | Stock-shaped query/search state with explicit two-end duplication where needed. |
| `make_hits_span_paired` | `make_hits_span` | Good | Shape and behavior are one-to-one. |
| `count_number_of_hits_to_keep_paired` | `count_number_of_hits_to_keep` | Good | Matches stock structure closely. |
| `copy_over_hits_to_be_kept_paired` | `copy_over_hits_to_be_kept` | Good | Matches stock structure closely. |
| `free_rejected_alignments_paired` | `free_rejected_alignments` | Good | Frees rejected raw per-end CIGAR pointers in stock style. |
| `hit_compare_byid_typed_paired` | `hit_compare_byid_typed` | Good | Already very close to stock shape. |
| `hit_compare_bysize_typed_paired` | `hit_compare_bysize_typed` | Good | Uses stock `db_getabundance(...)` through mapped paired targets. |
| `hit_compare_byid_paired` | `hit_compare_byid` | Good | Thin stock-shaped wrapper. |
| `hit_compare_bysize_paired` | `hit_compare_bysize` | Good | Thin stock-shaped wrapper. |
| `search_enough_kmers_paired` | `search_enough_kmers` | Good | Paired-specific difference is explicit R1+R2 k-mer count aggregation. |
| `search_topscores_paired` | `search_topscores` | Good | Structurally close to stock, with explicit two-end k-mer accumulation on top of paired dbindex accessors. |
| `search_acceptable_unaligned_paired` | `search_acceptable_unaligned` | Good | Reuses stock numeric filters, with paired length/abundance aggregation and paired prefix/suffix/self checks. |
| `search_acceptable_aligned_paired` | `search_acceptable_aligned` | Good | Consumes stock-shaped aggregate hit fields and writes back stock-style identity metrics. |
| `align_delayed_paired` | `align_delayed` | Good | Follows the stock delayed-alignment flow closely, with explicit per-end duplication only where paired logic requires it. |
| `search_onequery_paired` | `search_onequery` | Good | Stock-shaped at a high level: qprep, LMA setup, unique k-mers, top scores, delayed alignments, cleanup. |
| `search_findbest2_byid_paired` | `search_findbest2_byid` | Good | Already very close to stock. |
| `search_findbest2_bysize_paired` | `search_findbest2_bysize` | Good | Inherits the stock-like size comparator and direct abundance lookup. |
| `search_joinhits_paired` | `search_joinhits` | Good | Stock-shaped end-to-end: keep accepted/weak hits, free rejected alignments, sort in stock order. |

## `src/cluster_paired.cc` / `src/cluster_paired.h`

| Function / struct | Stock counterpart | Status | Review note |
| :-- | :-- | :-- | :-- |
| `clusterinfo_s_paired` | `clusterinfo_s` | Good | Retains per-end alignment ownership in stock style while carrying paired-specific cluster state. |
| `compare_byclusterno_paired` | `compare_byclusterno` | Good | Already close to stock. |
| `compare_byclusterabundance_paired` | `compare_byclusterabundance` | Good | Already close to stock. |
| `cluster_query_core_paired` | `cluster_query_core` | Good | Clean paired analog: load query state, handle reverse strand, call `search_onequery_paired()`. |
| `cluster_worker_paired` | `cluster_worker` | Good | Already close to stock. |
| `threads_worker_paired` | `threads_worker` | Good | Already close to stock. |
| `threads_wakeup_paired` | `threads_wakeup` | Good | Already close to stock. |
| `threads_init_paired` | `threads_init` | Good | Uses stock-like raw thread-info allocation and thread protocol. |
| `threads_exit_paired` | `threads_exit` | Good | Already close to stock. |
| `cluster_query_init_paired` | `cluster_query_init` | Good | Initializes query/header/sequence/k-mer/hit state with stock-style raw allocation and explicit two-end duplication. |
| `cluster_query_exit_paired` | `cluster_query_exit` | Good | Cleanup structure is close to stock, with paired-specific duplication for the two ends. |
| `cluster_core_results_hit_paired` | `cluster_core_results_hit` | Good | Performs stock-like immediate `H`-line UC emission and OTU-table accumulation for the supported paired surface. |
| `cluster_core_results_nohit_paired` | `cluster_core_results_nohit` | Good | Performs stock-like immediate `S`-line UC emission and OTU-table accumulation for the supported paired surface. |
| `clear_hits_paired` | post-query alignment cleanup loops in `cluster_core_parallel` and `cluster_core_serial` | Good | Frees raw per-end CIGAR ownership in the same places where stock frees single-end alignment strings. |
| `cluster_core_parallel_paired` | `cluster_core_parallel` | Good | Core loop matches stock structure closely: raw query arrays, same-round reconciliation, retained alignment ownership, DB/index insertion, immediate result handoff. |
| `cluster_core_serial_paired` | `cluster_core_serial` | Good | Serial loop mirrors stock structure closely and matches the parallel paired result surface. |
| `cluster_paired` | `cluster` | Good | Mirrors stock engine structure while replacing stock single-end load/output with paired-native analogs. |
| `cluster_unoise_paired` | `cluster_unoise` | Good | Thin stock-shaped wrapper declared in `cluster_paired.h`. |
