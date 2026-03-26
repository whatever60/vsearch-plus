# Low-Level Runtime and Header Naming Parity (Five Extended Commands)

This note summarizes parity status for:

- `vsearch --fastq_filter`
- `vsearch --fastx_uniques`
- `vsearch --cluster_unoise`
- `vsearch --uchime3_denovo`
- `vsearch --usearch_global`

It focuses on:

1. low-level runtime behavior (`--threads`, SIMD Needleman-Wunsch, linear-memory fallback)
2. output FASTA/FASTQ header naming conventions

## Scope of concrete examples

Examples below use real IDs from `real_sub/r1_10k.fq.gz` and `real_sub/r2_10k.fq.gz` (first 1,000 pairs), with `--sample realsub --sizeout` to make header conventions visible.

Representative input IDs:

```text
@sample=CRLM1A_A1 1 494035 LH00659:244:232WKKLT3:4:1123:33117:29695
@sample=CRLM1A_A1 2 494035 LH00659:244:232WKKLT3:4:1123:33117:29695
```

## 1) Low-level runtime parity

| Command | Stock multicore | Ext multicore | Stock SIMD NW + LMA fallback | Ext SIMD NW + LMA fallback | Why in ext |
|---|---|---|---|---|---|
| `fastq_filter` | No (single stream) | Same as stock | N/A (no NW alignment path) | N/A | Directly reuses stock command path; no extension code path here. |
| `fastx_uniques` | No (single-thread derep path) | No | N/A (no NW alignment path) | N/A | Paired mode is custom (`tav_fastx_uniques`), but operation is counting/aggregation only. |
| `cluster_unoise` | Yes (`cluster.cc` thread workers) | No | Yes (`search16` SIMD first, LMA fallback on overflow) | No SIMD-first; LMA-only alignment in paired filters | Paired mode is custom (`tav_cluster_unoise`) and bypasses stock `search_onequery`/`search16` pipeline. |
| `uchime3_denovo` | Yes (`chimera.cc` worker threads) | No | Yes (candidate/full-query alignments use SIMD path then LMA fallback) | No SIMD-first; LMA-based paired alignments | Paired mode is custom (`tav_uchime3_denovo`) and does not call stock chimera threading/search pipeline. |
| `usearch_global` | Yes (`search.cc` worker threads) | No | Yes (`searchcore.cc` uses `search16` with LMA fallback) | No SIMD-first; LMA-only per-end alignments | Paired mode is custom (`tav_usearch_global`) and bypasses stock threaded searchcore. |

### Key detail: what is reused vs engineered

- Reused low-level naming/formatting:
  - `fasta_print_general`
  - `fastq_print_general`
  - `header_fprint_strip`
- Reused low-level alignment post-processing in ext:
  - `align_trim` is called by `align_one_end_stock_style`.
- Engineered in ext:
  - paired candidate enumeration, paired filtering, paired output writing for `cluster_unoise`, `uchime3_denovo`, `usearch_global`.
- Not currently reused by ext paired paths:
  - stock threaded execution engines in `cluster.cc`, `chimera.cc`, `search.cc`
  - stock `search16` SIMD-first alignment path and overflow fallback logic.

## 2) Output header naming parity (with real IDs)

Important context:

- Most naming consistency comes from shared low-level printers (`fasta_print_general` / `fastq_print_general`), not from duplicated custom string logic.
- Because examples were run as a pipeline with `--sample realsub`, `;sample=realsub` can appear multiple times when a previous output header already contains `;sample=...` and the next command appends it again.

### A) `fastq_filter`

Status: extension parity is exact because this command is stock path for both single-end and paired-end.

```text
Stock output R1:
@sample=CRLM1A_A1 1 494035 LH00659:244:232WKKLT3:4:1123:33117:29695;sample=realsub;size=1

Stock output R2:
@sample=CRLM1A_A1 2 494035 LH00659:244:232WKKLT3:4:1123:33117:29695;sample=realsub;size=1
```

Why: direct reuse of stock `filter.cc` output calls to `fastq_print_general`.

### B) `fastx_uniques`

Status: nearly same naming convention; paired extension writes two coordinated FASTA files.

```text
Stock single-end (`--fastaout`):
>sample=CRLM1A_A1;sample=realsub;size=803

Extension paired R1 (`--fastaout`):
>sample=CRLM1A_A1;sample=realsub;size=749

Extension paired R2 (`--fastaout_rev`):
>sample=CRLM1A_A1;sample=realsub;size=749
```

Why:

- Base IDs come from extension logic (paired keying and record tracking).
- Final header formatting (`;sample=...`, `;size=...`, relabel/strip rules) comes from shared `fasta_print_general`.

### C) `cluster_unoise`

Status: naming convention is mostly aligned; paired extension emits synchronized left/right centroid files.

```text
Stock centroid (`--centroids`):
>sample=CRLM1A_A1;sample=realsub;sample=realsub;size=33

Extension centroid R1 (`--fastaout`):
>sample=CRLM1A_A1;sample=realsub;sample=realsub;size=839

Extension centroid R2 (`--fastaout_rev`):
>sample=CRLM1A_A1;sample=realsub;sample=realsub;size=839
```

Why:

- Paired centroid identity bookkeeping is extension logic.
- Header suffix conventions (`;sample`, `;size`, relabel/strip) are still from shared FASTA printer.

### D) `uchime3_denovo`

Status: shared convention + paired-specific explicit suffixing for FASTA pair outputs.

```text
Stock nonchimera (`--nonchimeras`):
>sample=CRLM1A_A1;sample=realsub;sample=realsub;sample=realsub;size=33

Extension nonchimera pair (`--nonchimeras`, interleaved):
>sample=CRLM1A_A1;sample=realsub;sample=realsub;size=839/1;sample=realsub;size=839
>sample=CRLM1A_A1;sample=realsub;sample=realsub;size=839/2;sample=realsub;size=839

Extension chimera pair (`--chimeras`, interleaved):
>sample=CRLR1A_A10;sample=realsub;sample=realsub;size=2/1;sample=realsub;size=2
>sample=CRLR1A_A10;sample=realsub;sample=realsub;size=2/2;sample=realsub;size=2
```

Why:

- Extension explicitly appends `/1` and `/2` in paired interleaved writer.
- Everything else in the final header format still comes from shared `fasta_print_general`.

### E) `usearch_global`

Status: shared convention + paired-specific explicit `/1` and `/2` in paired FASTA outputs.

```text
Stock matched (`--matched`):
>sample=CRLM1A_A1;sample=realsub;sample=realsub;size=803

Stock notmatched (`--notmatched`):
>sample=CRLM1A_A1;sample=realsub;sample=realsub;size=1

Extension notmatched pair (`--notmatched`, interleaved):
>sample=CRLM1A_A1;sample=realsub;size=749/1;sample=realsub;size=749
>sample=CRLM1A_A1;sample=realsub;size=749/2;sample=realsub;size=749

Extension matched pair example (`--matched`, interleaved; from `--id 0.80` run):
>sample=CRLM1A_A1;sample=realsub;size=749/1;sample=realsub;size=749
>sample=CRLM1A_A1;sample=realsub;size=749/2;sample=realsub;size=749
```

Extension matched naming follows the same pattern as extension notmatched (paired `/1` and `/2` interleaving with shared FASTA suffix logic).

Why:

- `/1` and `/2` handling is explicit extension engineering in paired output helpers.
- Base strip/relabel/sample/size behavior is inherited from shared FASTA output code.

## Bottom line

- Naming parity is generally strong because both stock and extension rely on the same low-level header-printing functions.
- Runtime backend parity is mixed:
  - full parity for `fastq_filter` (same stock path),
  - mostly N/A for `fastx_uniques` (no alignment backend),
  - currently non-parity for `cluster_unoise`, `uchime3_denovo`, and `usearch_global` on multicore and SIMD-first alignment behavior (paired extension uses custom single-thread LMA-based paths).
