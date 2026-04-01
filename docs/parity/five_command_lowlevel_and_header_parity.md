# Low-Level Runtime and Header Naming Parity (Five Paired-Aware Commands)

This note summarizes parity status for:

- `vsearch --fastq_filter`
- `vsearch --fastx_uniques`
- `vsearch --cluster_unoise`
- `vsearch --uchime3_denovo`
- `vsearch --usearch_global`

It focuses on:

1. low-level runtime behavior (`--threads`, SIMD Needleman-Wunsch, linear-memory fallback)
2. output FASTA/FASTQ header naming conventions
3. paired-orientation policy for R2

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
| `fastq_filter` | No (single stream) | No | N/A (no NW alignment path) | N/A | Stock mode (`--reverse`) remains on stock path; extension mode adds paired input/output CLI conventions and pair-level EE threshold evaluation. |
| `fastx_uniques` | No (single-thread derep path) | No | N/A (no NW alignment path) | N/A | Paired mode now runs through native `cpp/src/derep_paired.cc`, and the operation is still counting/aggregation only. |
| `cluster_unoise` | Yes (`cluster.cc` thread workers) | Yes (`cluster_paired.cc` worker/reconciliation path) | Yes (`search16` SIMD first, LMA fallback on overflow) | Yes (`search16` SIMD first, LMA fallback on overflow) | Paired mode now runs through a native stock-shaped engine in `cluster_paired.cc`; like stock `cluster.cc`, this path aligns candidates inline with `search16(..., 1, ...)` rather than using delayed `MAXDELAYED` batches. |
| `uchime3_denovo` | Yes (`chimera.cc` worker threads; denovo forces one worker) | Yes (`chimera_paired.cc` mirrors the stock threaded driver shape) | Yes (candidate/full-query alignments use SIMD path then LMA fallback) | Yes (paired part-screening uses shared paired dbindex counts/bitmaps; full-query paired alignments use SIMD-first with stock overflow fallback) | Paired mode now runs through a native stock-shaped engine in `chimera_paired.cc`, with shared stock parent-selection logic and paired-native output handling. |
| `usearch_global` | Yes (`search.cc` worker threads) | Yes (`search_paired.cc` worker threads) | Yes (`searchcore.cc` uses `search16` with LMA fallback) | Yes (`search16` SIMD first, LMA fallback on overflow, per end) | Paired mode now runs through a native stock-shaped engine in `search_paired.cc`, reusing the shared paired search core. |

### Key detail: what is reused vs engineered

- Reused low-level naming/formatting:
  - `fasta_print_general`
  - `fastq_print_general`
  - `header_fprint_strip`
- Reused shared filter kernels across stock/ext:
  - `search_unaligned_numeric_filters_pass`
  - `search_aligned_compute_identity_metrics`
  - `search_aligned_threshold_filters_pass`
- Reused low-level alignment post-processing in ext:
  - `align_trim` is called by `align_one_end_stock_style`.
- Reused shared CIGAR traversal kernel:
  - `parse_cigar_string` (via `parse_cigar_operations_from_string` wrapper in paired `uchime3_denovo`)
- Reused shared fastq EE threshold kernel:
  - `ee_thresholds_pass` (stock per-read and extension pair-level aggregation both call it)
- Engineered in ext/native paired layers:
  - paired candidate enumeration, paired filtering, paired output writing, and paired derep key formation where single-end stock code has no direct analogue.
- Native paired wrappers now exist for the core stock owners:
  - `cpp/src/derep_paired.cc`
  - `cpp/src/cluster_paired.cc`
  - `cpp/src/chimera_paired.cc`
  - `cpp/src/search_paired.cc`

## 2) Output header naming parity (with real IDs)

Important context:

- Most naming consistency comes from shared low-level printers (`fasta_print_general` / `fastq_print_general`), not from duplicated custom string logic.
- Because examples were run as a pipeline with `--sample realsub`, `;sample=realsub` can appear multiple times when a previous output header already contains `;sample=...` and the next command appends it again.

### A) `fastq_filter`

Status: stock paired behavior is preserved, and extension mode is intentionally additive.

```text
Stock output R1:
@sample=CRLM1A_A1 1 494035 LH00659:244:232WKKLT3:4:1123:33117:29695;sample=realsub;size=1

Stock output R2:
@sample=CRLM1A_A1 2 494035 LH00659:244:232WKKLT3:4:1123:33117:29695;sample=realsub;size=1
```

Why:

- Stock mode (`--reverse`) still reuses stock `filter.cc` output calls to `fastq_print_general`.
- Extension mode (`R1 R2` positional input or `--interleaved`) adds paired split outputs (`*2`) and uses pair-level `--fastq_maxee` / `--fastq_maxee_rate` checks.
- See `docs/parity/fastq_filter_paired_parity.md` for full mode/semantics details.

### B) `fastx_uniques`

Status: nearly same naming convention; native paired mode writes two coordinated FASTA files.

```text
Stock single-end (`--fastaout`):
>sample=CRLM1A_A1;sample=realsub;size=803

Extension paired R1 (`--fastaout`):
>sample=CRLM1A_A1;sample=realsub;size=749

Extension paired R2 (`--fastaout_rev`):
>sample=CRLM1A_A1;sample=realsub;size=749
```

Why:

- Base IDs come from native paired derep logic (paired keying and record tracking).
- Final header formatting (`;sample=...`, `;size=...`, relabel/strip rules) comes from shared `fasta_print_general`.

### C) `cluster_unoise`

Status: naming convention is mostly aligned; native paired mode emits synchronized left/right centroid files.

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

Status: shared convention; paired split FASTA outputs now keep the same base header on both R1 and R2 files (no extra `/1` or `/2`).

```text
Stock nonchimera (`--nonchimeras`):
>sample=CRLM1A_A1;sample=realsub;sample=realsub;sample=realsub;size=33

Extension nonchimera pair split output (`--nonchimeras` + `--nonchimeras2`):
R1 file (`--nonchimeras`):
>sample=CRLM1A_A1;sample=realsub;sample=realsub;size=839;sample=realsub;size=839
R2 file (`--nonchimeras2`):
>sample=CRLM1A_A1;sample=realsub;sample=realsub;size=839;sample=realsub;size=839

Extension chimera pair split output (`--chimeras` + `--chimeras2`):
R1 file (`--chimeras`):
>sample=CRLR1A_A10;sample=realsub;sample=realsub;size=2;sample=realsub;size=2
R2 file (`--chimeras2`):
>sample=CRLR1A_A10;sample=realsub;sample=realsub;size=2;sample=realsub;size=2
```

Why:

- Paired split FASTA writers now emit the same base header in both files.
- Final header formatting still comes from shared `fasta_print_general`.

### E) `usearch_global`

Status: shared convention; paired split FASTA outputs now keep the same base header on both R1 and R2 files.

```text
Stock matched (`--matched`):
>sample=CRLM1A_A1;sample=realsub;sample=realsub;size=803

Stock notmatched (`--notmatched`):
>sample=CRLM1A_A1;sample=realsub;sample=realsub;size=1

Extension notmatched pair split output (`--notmatched` + `--notmatched2`):
R1 file (`--notmatched`):
>sample=CRLM1A_A1;sample=realsub;size=749;sample=realsub;size=749
R2 file (`--notmatched2`):
>sample=CRLM1A_A1;sample=realsub;size=749;sample=realsub;size=749

Extension matched pair split output example (`--matched` + `--match2`; from `--id 0.80` run):
R1 file (`--matched`):
>sample=CRLM1A_A1;sample=realsub;size=749;sample=realsub;size=749
R2 file (`--match2`):
>sample=CRLM1A_A1;sample=realsub;size=749;sample=realsub;size=749
```

Extension matched naming follows the same pattern as extension notmatched (same base header in R1/R2 split files).

Why:

- Base strip/relabel/sample/size behavior is inherited from shared FASTA output code.
- No extra `/1` or `/2` suffix is appended by paired output helpers.

## Bottom line

- Naming parity is generally strong because both stock and paired-native paths rely on the same low-level header-printing functions.
- Runtime backend parity is mixed:
  - `fastq_filter` has mixed parity: stock paired mode is unchanged, while extension mode intentionally changes CLI/output conventions and EE-threshold semantics,
  - mostly N/A for `fastx_uniques` (no alignment backend),
  - `cluster_unoise`, `uchime3_denovo`, and `usearch_global` now have SIMD-first + overflow-fallback parity in paired mode,
  - these three commands now also honor `--threads` through native paired orchestration code paths.

## Orientation policy status

- Extension pipeline policy is now unified: R2 remains in original read orientation at external I/O for paired
  `fastx_uniques`, `cluster_unoise`, `uchime3_denovo`, and `usearch_global`.
- Previous behavior that persisted R2 as reverse-complemented is removed.
- Internal orientation handling may still occur where algorithmically required (for example RDP classifier internals),
  but this does not change external orientation contracts.
- See `docs/parity/paired_orientation_unification.md` for developer details.
